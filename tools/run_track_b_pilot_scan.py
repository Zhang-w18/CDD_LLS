"""plan-023: Track B pilot scan — transparent sliding-window receiver Pareto scan.

Implements:
  - I_QAM BICM MI lookup table for unit-energy 16QAM (Gauss-Hermite quadrature)
  - Candidate families B1 (transparent small-delay CDD), B2 (PRG precoder cycling),
    B3 (CC: phase-continuous precoder cycling), B4 (CN: capped-slope N-series),
    B5 (QC folded-grid CDD benchmark), B6 (P6/P7 outage-only set, 48PRB)
  - MC outage engine (flat model g = V h, h ~ CN(0, I_8), common random numbers)
  - Closed-form mismatched-NMSE for receivers R1 (4RB sliding-window MMSE),
    R2 (wideband robust Wiener), R3 (V-aware matched LMMSE)
  - tau_span diagnostic (95% energy span of composite delay profile)
  - Phase 0a calibration (closed-form R3 vs MC with existing Alg1 estimator tool)
  - Phase 0b calibration (outage vs ideal-CSI BLER parallelism)

Stages (--stage): iqam | scan | calib0a | calib0b | analyze
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import time
from dataclasses import dataclass, field
from pathlib import Path
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# ----------------------------------------------------------------------------
# Global constants (plan-023 §2/§5)
# ----------------------------------------------------------------------------
NT = 8                       # transmit branches
SCS_HZ = 30e3                # subcarrier spacing
SF = 24                      # DMRS comb spacing (subcarriers)
R_SE = 4.0 * 553.0 / 1024.0  # spectral-efficiency threshold, bit/RE (16QAM 553/1024)
N_MC = 200_000               # outage Monte Carlo samples
SEED_MC = 20260709
SNR_GRID_DB = np.arange(0.0, 24.0 + 1e-9, 0.5)   # outage SNR grid
NMSE_SNRS_DB = (8.0, 16.0, 24.0)                  # closed-form NMSE SNR points
COND_LIMIT = 1e10
TAU_W_R1_MAIN_NS = 300.0
TAU_W_R1_SENS_NS = (100.0, 694.0)
TAU_W_R2_NS = (200.0, 694.0, 1389.0)

# I_QAM table spec (plan §3.1)
IQAM_DB_MIN = -20.0
IQAM_DB_MAX = 35.0
IQAM_DB_STEP = 0.1

# Histogram spec for the outage engine (aligned to the 0.1 dB table grid)
HIST_DB_MIN = -60.0
HIST_DB_MAX = 20.0


def db(x: float) -> float:
    return float(10.0 * np.log10(x)) if x > 0 and np.isfinite(x) else float("nan")


def write_csv(rows: List[Dict[str, object]], path: Path, fields: Optional[List[str]] = None) -> None:
    if not rows:
        return
    if fields is None:
        fields = []
        for row in rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


# ----------------------------------------------------------------------------
# I_QAM lookup table: BICM MI of unit-energy Gray 16QAM in AWGN
# ----------------------------------------------------------------------------
def iqam_16qam_bicm(rho_lin: np.ndarray, n_gh: int = 96) -> np.ndarray:
    """BICM MI (bit/RE) of unit-energy Gray-mapped 16QAM at linear receive SNR rho.

    Decomposes into two independent Gray 4-PAM dimensions.
    y_r = sqrt(rho) * a + n_r with n_r ~ N(0, 1/2), a in {+-1, +-3}/sqrt(10).
    Gauss-Hermite: y_r = sqrt(rho) a + t, weights w/sqrt(pi) (sigma*sqrt(2) = 1).
    """
    rho = np.asarray(rho_lin, dtype=np.float64).reshape(-1)
    levels = np.array([-3.0, -1.0, 1.0, 3.0]) / math.sqrt(10.0)
    # Gray mapping over levels [-3,-1,1,3]: b0 = [0,0,1,1], b1 = [0,1,1,0]
    bit_maps = [np.array([0, 0, 1, 1]), np.array([0, 1, 1, 0])]
    t_nodes, w_nodes = np.polynomial.hermite.hermgauss(n_gh)
    w_norm = w_nodes / math.sqrt(math.pi)

    out = np.zeros_like(rho)
    for gi, r in enumerate(rho):
        srt = math.sqrt(r)
        # y grid: (4 tx symbols, n_gh nodes)
        y = srt * levels[:, None] + t_nodes[None, :]
        # metric(y, a') = -(y - sqrt(r) a')^2  -> (4, n_gh, 4 hypotheses)
        d = y[:, :, None] - srt * levels[None, None, :]
        m = -d * d
        m_max = m.max(axis=2, keepdims=True)
        e = np.exp(m - m_max)
        log_num = np.log(e.sum(axis=2)) + m_max[:, :, 0]      # log sum over all a'
        mi = 0.0
        for bmap in bit_maps:
            # log sum over the hypothesis set matching the transmitted bit
            li = 1.0
            for bval in (0, 1):
                sel = bmap == bval
                es = np.exp(m[:, :, sel] - m_max)
                with np.errstate(divide="ignore"):
                    # rows whose transmitted bit != bval can underflow; they are unused below
                    log_den = np.log(es.sum(axis=2)) + m_max[:, :, 0]
                # rows (tx symbols) whose bit == bval
                tx_rows = np.where(bmap == bval)[0]
                # E over those tx symbols and noise: log2(num/den)
                contrib = (log_num[tx_rows] - log_den[tx_rows]) / math.log(2.0)
                li -= 0.25 * float((contrib * w_norm[None, :]).sum())
            mi += li
        out[gi] = 2.0 * mi
    return np.clip(out, 0.0, 4.0)


def build_iqam_table(cache_path: Path) -> np.ndarray:
    if cache_path.exists():
        data = np.load(cache_path)
        return data["table"]
    grid_db = np.arange(IQAM_DB_MIN, IQAM_DB_MAX + 1e-9, IQAM_DB_STEP)
    table = iqam_16qam_bicm(10.0 ** (grid_db / 10.0))
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(cache_path, grid_db=grid_db, table=table)
    return table


# ----------------------------------------------------------------------------
# Candidate construction (plan §4). Phase phi is (K, NT); V = exp(1j*phi).
# ----------------------------------------------------------------------------
@dataclass
class Candidate:
    cid: str
    family: str
    K: int
    phase: np.ndarray
    params: Dict[str, object] = field(default_factory=dict)
    outage_only: bool = False

    @property
    def V(self) -> np.ndarray:
        return np.exp(1j * self.phase)


def grid_step_ns(K: int) -> float:
    return 1e9 / (K * SCS_HZ)


def cdd_phase(K: int, j_by_branch: Sequence[int]) -> np.ndarray:
    k = np.arange(K, dtype=np.float64)
    j = np.asarray(j_by_branch, dtype=np.float64)
    return -2.0 * math.pi * k[:, None] * j[None, :] / float(K)


def wrap_pi(x: np.ndarray) -> np.ndarray:
    """Wrap to (-pi, pi]."""
    return -(np.mod(-x + math.pi, 2.0 * math.pi) - math.pi)


def build_b1(K: int) -> List[Candidate]:
    out = []
    for n_eff in (2, 4, 8):
        for span_ns in (50, 100, 200, 400, 800, 1389):
            delta_ns = float(span_ns) / (n_eff - 1)
            j = [int(round((n % n_eff) * delta_ns * 1e-9 * K * SCS_HZ)) for n in range(NT)]
            out.append(Candidate(
                cid=f"B1_neff{n_eff}_span{span_ns}",
                family="B1_cdd",
                K=K,
                phase=cdd_phase(K, j),
                params={"n_eff": n_eff, "span_ns": span_ns, "j_list": j},
            ))
    return out


def build_b2(K: int) -> List[Candidate]:
    out = []
    n = np.arange(NT, dtype=np.float64)
    k = np.arange(K)
    for p_rb in (2, 4, 8):
        s = (k // (12 * p_rb)) % 8
        phase = 2.0 * math.pi * n[None, :] * s[:, None].astype(np.float64) / 8.0
        out.append(Candidate(
            cid=f"B2_prg{p_rb}",
            family="B2_prg_cycling",
            K=K,
            phase=phase,
            params={"p_rb": p_rb},
        ))
    return out


BITREV = [0, 4, 2, 6, 1, 5, 3, 7]


def build_cc_phase(K: int, n_seg: int, T: int, schedule: str) -> np.ndarray:
    L = K // n_seg
    assert K % n_seg == 0
    n = np.arange(NT, dtype=np.float64)
    if schedule == "seq":
        r = np.array([s % 8 for s in range(n_seg)], dtype=np.float64)
    else:
        r = np.array([BITREV[s % 8] for s in range(n_seg)], dtype=np.float64)
    c = 2.0 * math.pi * n[None, :] * r[:, None] / 8.0   # (n_seg, NT)
    k = np.arange(K)
    phase = c[k // L, :].copy()
    half = T // 2
    for s in range(1, n_seg):
        kb = s * L
        ks = np.arange(kb - half, kb + half)
        delta = wrap_pi(c[s, :] - c[s - 1, :])
        frac = (ks - (kb - half)).astype(np.float64) / float(T)
        phase[ks, :] = c[s - 1, :][None, :] + delta[None, :] * frac[:, None]
    return phase


def build_b3(K: int) -> List[Candidate]:
    n_seg_list = {576: (8, 12, 16, 24), 288: (4, 6, 8, 12)}[K]
    out = []
    for n_seg in n_seg_list:
        L = K // n_seg
        for T in (6, 12, 24):
            if T > L / 3.0:
                continue
            for schedule in ("seq", "bitrev"):
                out.append(Candidate(
                    cid=f"B3_cc_nseg{n_seg}_T{T}_{schedule}",
                    family="B3_cc",
                    K=K,
                    phase=build_cc_phase(K, n_seg, T, schedule),
                    params={"n_seg": n_seg, "T": T, "schedule": schedule},
                ))
    return out


def build_b4(K: int) -> List[Candidate]:
    out = []
    n_seg = 8
    L = K // n_seg
    for a in (1, 2):
        for idx in range(50):
            rng = np.random.default_rng([23, a, idx])
            j = rng.integers(0, a + 1, size=(n_seg, NT))          # grid-unit slopes
            slope = j[np.arange(K) // L, :].astype(np.float64)     # (K, NT)
            phase = np.zeros((K, NT), dtype=np.float64)
            phase[1:, :] = -2.0 * math.pi / K * np.cumsum(slope[:-1, :], axis=0)
            out.append(Candidate(
                cid=f"B4_cn_a{a}_s{idx:02d}",
                family="B4_cn",
                K=K,
                phase=phase,
                params={"a": a, "seed_offset": idx, "j_list": j.tolist()},
            ))
    return out


def build_b5(K: int) -> List[Candidate]:
    j = [9 * n for n in range(NT)] if K == 576 else [n for n in range(NT)]
    return [Candidate(
        cid="B5_qc",
        family="B5_qc",
        K=K,
        phase=cdd_phase(K, j),
        params={"j_list": j},
    )]


def build_b6(K: int) -> List[Candidate]:
    if K != 576:
        return []
    out = []
    for sigma in (1, 2, 4, 9):
        j = [sigma * n for n in range(NT)]
        out.append(Candidate(
            cid=f"B6_arith_s{sigma}",
            family="B6_appendix",
            K=K,
            phase=cdd_phase(K, j),
            params={"sigma": sigma, "j_list": j},
            outage_only=True,
        ))
    sidon = [0, 1, 3, 7, 12, 20, 30, 65]
    out.append(Candidate(
        cid="B6_sidon",
        family="B6_appendix",
        K=K,
        phase=cdd_phase(K, sidon),
        params={"j_list": sidon},
        outage_only=True,
    ))
    return out


def build_all_candidates(K: int) -> List[Candidate]:
    return build_b1(K) + build_b2(K) + build_b3(K) + build_b4(K) + build_b5(K) + build_b6(K)


# ----------------------------------------------------------------------------
# tau_span: 95%-energy minimal circular window of composite delay profile
# ----------------------------------------------------------------------------
def tau_span_ns(cand: Candidate) -> float:
    V = cand.V
    prof = np.sum(np.abs(np.fft.ifft(V, axis=0)) ** 2, axis=1)   # (K,)
    total = float(prof.sum())
    K = cand.K
    target = 0.95 * total
    ext = np.concatenate([prof, prof])
    csum = np.concatenate([[0.0], np.cumsum(ext)])
    best = K
    lo = 0
    for hi in range(1, 2 * K + 1):
        while csum[hi] - csum[lo] >= target and hi - lo <= K:
            best = min(best, hi - lo)
            lo += 1
        if hi - lo > K:
            lo = hi - K
    return float((best - 1) * grid_step_ns(K))


# ----------------------------------------------------------------------------
# MC outage engine (plan §3.1)
# ----------------------------------------------------------------------------
def draw_h(n_mc: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return (rng.standard_normal((n_mc, NT)) + 1j * rng.standard_normal((n_mc, NT))) / math.sqrt(2.0)


def build_shift_table(table: np.ndarray, snr_grid_db: np.ndarray) -> np.ndarray:
    """Tmat[b, j] = I_QAM(HIST_DB_MIN + 0.1*b + snr_j), clamped to table ends."""
    n_bins = int(round((HIST_DB_MAX - HIST_DB_MIN) / IQAM_DB_STEP)) + 1
    b = np.arange(n_bins)
    base = np.round((HIST_DB_MIN - IQAM_DB_MIN) / IQAM_DB_STEP).astype(int)
    shifts = np.round(snr_grid_db / IQAM_DB_STEP).astype(int)
    idx = b[:, None] + base + shifts[None, :]
    idx = np.clip(idx, 0, len(table) - 1)
    return table[idx]   # (n_bins, n_snr)


def outage_curves(
    cands: List[Candidate],
    h_all: np.ndarray,
    table: np.ndarray,
    snr_grid_db: np.ndarray = SNR_GRID_DB,
    chunk: int = 8000,
    progress: bool = True,
) -> Dict[str, np.ndarray]:
    """Returns {cid: P_out over snr_grid}. Exact per-RE linear interpolation of the
    dB-axis lookup is preserved via fractional histogram weight splitting."""
    n_mc = h_all.shape[0]
    n_bins = int(round((HIST_DB_MAX - HIST_DB_MIN) / IQAM_DB_STEP)) + 1
    tmat = build_shift_table(table, snr_grid_db)              # (n_bins, n_snr)
    out: Dict[str, np.ndarray] = {}
    for ci, cand in enumerate(cands):
        K = cand.K
        Vt = cand.V.T.copy()                                  # (NT, K)
        counts = np.zeros(len(snr_grid_db), dtype=np.int64)
        t0 = time.time()
        for start in range(0, n_mc, chunk):
            h = h_all[start:start + chunk]
            nc = h.shape[0]
            G = h @ Vt                                        # (nc, K)
            x_db = 10.0 * np.log10(np.abs(G) ** 2 / NT + 1e-300)
            np.clip(x_db, HIST_DB_MIN, HIST_DB_MAX - IQAM_DB_STEP - 1e-6, out=x_db)
            t = (x_db - HIST_DB_MIN) / IQAM_DB_STEP
            i = np.floor(t).astype(np.int64)
            frac = t - i
            rows = np.arange(nc, dtype=np.int64)[:, None] * n_bins
            flat_lo = (rows + i).ravel()
            flat_hi = (rows + i + 1).ravel()
            w_lo = (1.0 - frac).ravel()
            w_hi = frac.ravel()
            hist = np.bincount(flat_lo, weights=w_lo, minlength=nc * n_bins)
            hist += np.bincount(flat_hi, weights=w_hi, minlength=nc * n_bins)
            hist = hist.reshape(nc, n_bins)
            mi = hist @ tmat / float(K)                       # (nc, n_snr)
            counts += (mi < R_SE).sum(axis=0)
        out[cand.cid] = counts / float(n_mc)
        if progress:
            print(f"[outage] {ci+1}/{len(cands)} {cand.cid} ({time.time()-t0:.1f}s)", flush=True)
    return out


def outage_snr_at(p_out: np.ndarray, snr_grid_db: np.ndarray, target: float) -> float:
    """Log-domain interpolation of the SNR where P_out crosses `target`."""
    p = np.asarray(p_out, dtype=np.float64)
    eps = 0.5 / N_MC
    p = np.maximum(p, eps)
    for a in range(len(p) - 1):
        if p[a] >= target > p[a + 1]:
            l0, l1 = math.log10(p[a]), math.log10(p[a + 1])
            t = (l0 - math.log10(target)) / (l0 - l1)
            return float(snr_grid_db[a] + t * (snr_grid_db[a + 1] - snr_grid_db[a]))
    return float("nan")


# ----------------------------------------------------------------------------
# Closed-form mismatched NMSE (plan §3.2)
# ----------------------------------------------------------------------------
def uniform_pdp_cov(k_a: np.ndarray, k_b: np.ndarray, tau_w_ns: float) -> np.ndarray:
    """Assumed covariance R~[k,l] = NT * exp(-j*pi*(k-l)*df*tau) * sinc((k-l)*df*tau)."""
    d = (np.asarray(k_a, dtype=np.float64)[:, None] - np.asarray(k_b, dtype=np.float64)[None, :])
    x = d * SCS_HZ * tau_w_ns * 1e-9
    return NT * np.exp(-1j * math.pi * x) * np.sinc(x)


@dataclass
class WindowSpec:
    pilots: np.ndarray   # pilot subcarrier indices in window
    data: np.ndarray     # target data subcarrier indices (block b non-pilot)


def r1_windows(K: int, prg_rb: Optional[int]) -> List[WindowSpec]:
    wins = []
    for b in range(K // 12):
        blk = np.arange(12 * b, 12 * b + 12)
        if prg_rb is None:
            ws = int(np.clip(12 * b + 6 - 24, 0, K - 48))
            win = np.arange(ws, ws + 48)
        else:
            prg_sc = 12 * prg_rb
            ps = (12 * b) // prg_sc * prg_sc
            if prg_sc <= 48:
                win = np.arange(ps, min(ps + prg_sc, K))
            else:
                ws = int(np.clip(12 * b + 6 - 24, ps, min(ps + prg_sc, K) - 48))
                win = np.arange(ws, ws + 48)
        pilots = win[win % SF == 0]
        data = blk[blk % SF != 0]
        wins.append(WindowSpec(pilots=pilots, data=data))
    return wins


@dataclass
class LinEst:
    A: np.ndarray        # (|D|, |P|) estimator
    pilots: np.ndarray
    data: np.ndarray
    cond: float
    loaded: bool


def make_estimator(
    r_tilde_pp: np.ndarray,
    r_tilde_dp: np.ndarray,
    noise_var_ls: float,
) -> Tuple[np.ndarray, float, bool]:
    n_p = r_tilde_pp.shape[0]
    M = r_tilde_pp + noise_var_ls * np.eye(n_p)
    cond = float(np.linalg.cond(M))
    loaded = False
    if cond > COND_LIMIT:
        M = M + (1e-12 * np.trace(M).real / n_p) * np.eye(n_p)
        loaded = True
        cond = float(np.linalg.cond(M))
    A = np.linalg.solve(M.conj().T, r_tilde_dp.conj().T).conj().T
    return A, cond, loaded


def mismatch_mse(A: np.ndarray, R_pp: np.ndarray, R_pd: np.ndarray, tr_R_dd: float, noise_var_ls: float) -> float:
    """tr(A (R_PP + s I) A^H) - 2 Re tr(A R_PD) + tr(R_DD)."""
    B = A @ (R_pp + noise_var_ls * np.eye(R_pp.shape[0]))
    t1 = float(np.sum(B * A.conj()).real)
    t2 = float(np.sum(A * R_pd.T).real)
    return t1 - 2.0 * t2 + tr_R_dd


def r1_nmse(
    V: np.ndarray,
    wins: List[WindowSpec],
    tau_w_ns: float,
    snr_db: float,
    est_cache: Dict[object, LinEst],
) -> Tuple[float, float, bool]:
    """Returns (NMSE_linear, worst_window_NMSE_linear, any_loading_flag)."""
    snr = 10.0 ** (snr_db / 10.0)
    noise_var_ls = NT / snr / 2.0
    num = 0.0
    den = 0.0
    worst = 0.0
    flagged = False
    for w in wins:
        key = (tuple(w.pilots.tolist()), tuple(w.data.tolist()), tau_w_ns, snr_db)
        est = est_cache.get(key)
        if est is None:
            r_pp = uniform_pdp_cov(w.pilots, w.pilots, tau_w_ns)
            r_dp = uniform_pdp_cov(w.data, w.pilots, tau_w_ns)
            A, cond, loaded = make_estimator(r_pp, r_dp, noise_var_ls)
            est = LinEst(A=A, pilots=w.pilots, data=w.data, cond=cond, loaded=loaded)
            est_cache[key] = est
        Vp = V[w.pilots, :]
        Vd = V[w.data, :]
        R_pp = Vp @ Vp.conj().T
        R_pd = Vp @ Vd.conj().T
        tr_dd = float(NT * len(w.data))
        mse = mismatch_mse(est.A, R_pp, R_pd, tr_dd, noise_var_ls)
        num += mse
        den += tr_dd
        worst = max(worst, mse / tr_dd)
        flagged = flagged or est.loaded
    return num / den, worst, flagged


def wideband_nmse_assumed(
    V: np.ndarray,
    K: int,
    tau_w_ns: float,
    snr_db: float,
    est_cache: Dict[object, LinEst],
) -> Tuple[float, bool]:
    """R2: single full-band window, robust Wiener prior."""
    snr = 10.0 ** (snr_db / 10.0)
    noise_var_ls = NT / snr / 2.0
    pilots = np.arange(0, K, SF)
    data = np.array([k for k in range(K) if k % SF != 0])
    key = ("R2", K, tau_w_ns, snr_db)
    est = est_cache.get(key)
    if est is None:
        r_pp = uniform_pdp_cov(pilots, pilots, tau_w_ns)
        r_dp = uniform_pdp_cov(data, pilots, tau_w_ns)
        A, cond, loaded = make_estimator(r_pp, r_dp, noise_var_ls)
        est = LinEst(A=A, pilots=pilots, data=data, cond=cond, loaded=loaded)
        est_cache[key] = est
    Vp = V[pilots, :]
    Vd = V[data, :]
    R_pp = Vp @ Vp.conj().T
    R_pd = Vp @ Vd.conj().T
    tr_dd = float(NT * len(data))
    mse = mismatch_mse(est.A, R_pp, R_pd, tr_dd, noise_var_ls)
    return mse / tr_dd, est.loaded


def matched_nmse(V: np.ndarray, K: int, snr_db: float, targets: str = "nonpilot") -> Tuple[float, float, bool]:
    """R3: V-aware matched LMMSE (R~ = R_g). Returns (NMSE, cond, loaded)."""
    snr = 10.0 ** (snr_db / 10.0)
    noise_var_ls = NT / snr / 2.0
    pilots = np.arange(0, K, SF)
    if targets == "all":
        data = np.arange(K)
    else:
        data = np.array([k for k in range(K) if k % SF != 0])
    Vp = V[pilots, :]
    Vd = V[data, :]
    R_pp = Vp @ Vp.conj().T
    R_dp = Vd @ Vp.conj().T
    R_pd = Vp @ Vd.conj().T
    A, cond, loaded = make_estimator(R_pp, R_dp, noise_var_ls)
    tr_dd = float(NT * len(data))
    mse = mismatch_mse(A, R_pp, R_pd, tr_dd, noise_var_ls)
    return mse / tr_dd, cond, loaded


# ----------------------------------------------------------------------------
# Stage: scan
# ----------------------------------------------------------------------------
def stage_scan(out_dir: Path, n_prbs: int, n_mc: int, table: np.ndarray, smoke: bool = False) -> None:
    K = n_prbs * 12
    cands = build_all_candidates(K)
    if smoke:
        cands = [c for c in cands if c.family in ("B1_cdd", "B2_prg_cycling", "B5_qc")][:6] + \
                [c for c in cands if c.family == "B3_cc"][:2]
    print(f"[scan] {n_prbs} PRB, {len(cands)} candidates, n_mc={n_mc}", flush=True)

    h_all = draw_h(n_mc, SEED_MC)
    curves = outage_curves(cands, h_all, table)

    est_cache: Dict[object, LinEst] = {}
    wins_free = r1_windows(K, None)
    wins_prg = {p: r1_windows(K, p) for p in (2, 4, 8)}

    rows = []
    curve_rows = []
    for ci, cand in enumerate(cands):
        V = cand.V
        p_out = curves[cand.cid]
        row: Dict[str, object] = {
            "id": cand.cid,
            "family": cand.family,
            "n_prbs": n_prbs,
            "tau_span_ns": round(tau_span_ns(cand), 2),
            "outage10_dB": round(outage_snr_at(p_out, SNR_GRID_DB, 0.10), 4),
            "outage1_dB": round(outage_snr_at(p_out, SNR_GRID_DB, 0.01), 4),
        }
        for key in ("n_eff", "span_ns", "p_rb", "n_seg", "T", "schedule", "a", "seed_offset", "sigma"):
            row[key] = cand.params.get(key, "")
        row["j_list"] = json.dumps(cand.params.get("j_list", "")) if "j_list" in cand.params else ""
        cond_flag = False
        if not cand.outage_only:
            wins = wins_prg[cand.params["p_rb"]] if cand.family == "B2_prg_cycling" else wins_free
            for snr_db in NMSE_SNRS_DB:
                nm, worst, fl = r1_nmse(V, wins, TAU_W_R1_MAIN_NS, snr_db, est_cache)
                row[f"nmse_R1_{int(snr_db)}dB"] = round(db(nm), 4)
                cond_flag |= fl
                if snr_db == 16.0:
                    row["nmse_R1_worst_window_16dB"] = round(db(worst), 4)
            for tw in TAU_W_R1_SENS_NS:
                nm, _, fl = r1_nmse(V, wins, tw, 16.0, est_cache)
                row[f"nmse_R1_tw{int(tw)}_16dB"] = round(db(nm), 4)
                cond_flag |= fl
            for tw in TAU_W_R2_NS:
                nm, fl = wideband_nmse_assumed(V, K, tw, 16.0, est_cache)
                row[f"nmse_R2_tw{int(tw)}_16dB"] = round(db(nm), 4)
                cond_flag |= fl
            if cand.family == "B5_qc":
                for snr_db in NMSE_SNRS_DB:
                    nm, cond, fl = matched_nmse(V, K, snr_db)
                    row[f"nmse_R3_{int(snr_db)}dB"] = round(db(nm), 4)
                    cond_flag |= fl
        row["cond_flag"] = int(cond_flag)
        rows.append(row)
        crow = {"id": cand.cid, "family": cand.family}
        for s, p in zip(SNR_GRID_DB, p_out):
            crow[f"p{s:g}"] = f"{p:.6g}"
        curve_rows.append(crow)
        if (ci + 1) % 20 == 0:
            print(f"[nmse] {ci+1}/{len(cands)} done", flush=True)
            write_csv(rows, out_dir / f"candidates_{n_prbs}prb_partial.csv")

    write_csv(rows, out_dir / f"candidates_{n_prbs}prb.csv")
    write_csv(curve_rows, out_dir / f"outage_curves_{n_prbs}prb.csv")
    part = out_dir / f"candidates_{n_prbs}prb_partial.csv"
    if part.exists():
        part.unlink()
    print(f"[scan] wrote {out_dir / f'candidates_{n_prbs}prb.csv'}", flush=True)


# ----------------------------------------------------------------------------
# Stage: calib0a — closed-form R3 vs MC using the existing Alg1 estimator tool
# ----------------------------------------------------------------------------
def calib_candidates() -> List[Candidate]:
    qc = build_b5(576)[0]
    cc = Candidate(
        cid="B3_cc_nseg12_T12_seq",
        family="B3_cc",
        K=576,
        phase=build_cc_phase(576, 12, 12, "seq"),
        params={"n_seg": 12, "T": 12, "schedule": "seq"},
    )
    return [qc, cc]


def stage_calib0a(out_dir: Path, trials: int = 300) -> bool:
    from cdd_lls.core.config import ResourceConfig
    from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
    from tools.search_precoder_design_alg1 import make_alg1_full_cov_estimator

    resource = ResourceConfig(
        carrier_bandwidth_mhz=100.0, scs_khz=30, n_fft=4096, n_prbs=48,
        pdsch_n_symbols=10, dmrs_symbol_indices=[2, 7], dmrs_spacing_sc=SF, prg_size_rb=4,
    )
    grid = build_resource_grid(resource)
    K = grid.n_sc
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    pdp_flat = np.array([1.0])

    rows = []
    all_pass = True
    for cand in calib_candidates():
        V = cand.V
        for snr_db in NMSE_SNRS_DB:
            snr = 10.0 ** (snr_db / 10.0)
            noise_var_ls = NT / snr / 2.0
            # closed-form R3 over all K subcarriers (matches estimator target set)
            nm_cf, cond, _ = matched_nmse(V, K, snr_db, targets="all")
            # MC with the existing tool's estimator matrix
            est = make_alg1_full_cov_estimator(grid, pdp_flat, V, noise_var_ls, 0.0)
            W = est.matrix   # (K, N_p)
            rng = np.random.default_rng([SEED_MC, 77, int(snr_db)])
            err = 0.0
            pwr = 0.0
            for _ in range(trials):
                h = (rng.standard_normal(NT) + 1j * rng.standard_normal(NT)) / math.sqrt(2.0)
                g = V @ h
                noise = (rng.standard_normal(len(pilot_local)) + 1j * rng.standard_normal(len(pilot_local))) \
                    * math.sqrt(noise_var_ls / 2.0)
                g_hat = W @ (g[pilot_local] + noise)
                err += float(np.sum(np.abs(g_hat - g) ** 2))
                pwr += float(np.sum(np.abs(g) ** 2))
            nm_mc = err / pwr
            diff = db(nm_mc) - db(nm_cf)
            ok = abs(diff) <= 0.3
            all_pass &= ok
            rows.append({
                "id": cand.cid, "snr_db": snr_db,
                "nmse_closed_form_dB": round(db(nm_cf), 4),
                "nmse_mc_dB": round(db(nm_mc), 4),
                "diff_dB": round(diff, 4),
                "trials": trials,
                "rpp_cond": f"{cond:.4g}",
                "pass": int(ok),
            })
            print(f"[0a] {cand.cid} snr={snr_db:g}: cf={db(nm_cf):.3f} dB mc={db(nm_mc):.3f} dB "
                  f"diff={diff:+.3f} dB {'PASS' if ok else 'FAIL'}", flush=True)
    calib_dir = out_dir / "calibration"
    calib_dir.mkdir(parents=True, exist_ok=True)
    write_csv(rows, calib_dir / "calib0a_nmse.csv")
    print(f"[0a] gate {'PASSED' if all_pass else 'FAILED'}", flush=True)
    return all_pass


# ----------------------------------------------------------------------------
# Stage: calib0b — outage vs ideal-CSI BLER parallelism (P3)
# ----------------------------------------------------------------------------
def stage_calib0b(out_dir: Path, table: np.ndarray, trials: int = 400, n_mc: int = N_MC) -> bool:
    from cdd_lls.core.config import ResourceConfig
    from cdd_lls.core.mcs import build_tb_layout, get_mcs
    from cdd_lls.phy.ldpc import SionnaLDPCAdapter
    from cdd_lls.phy.qam import qam_demapper_maxlog, qam_modulate
    from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
    from tools.run_v_design_piecewise_tradeoff import decode_same_tb_batch, equalize_mrc

    cands = calib_candidates()
    K = 576
    h_all = draw_h(n_mc, SEED_MC)
    curves = outage_curves(cands, h_all, table)
    calib_dir = out_dir / "calibration"
    calib_dir.mkdir(parents=True, exist_ok=True)

    targets = {c.cid: (outage_snr_at(curves[c.cid], SNR_GRID_DB, 0.10),
                       outage_snr_at(curves[c.cid], SNR_GRID_DB, 0.01)) for c in cands}
    print(f"[0b] outage anchors: {targets}", flush=True)

    lo = math.floor(min(t[0] for t in targets.values()))
    hi = math.ceil(max(t[1] for t in targets.values())) + 7
    snr_bler_grid = np.arange(lo, hi + 1e-9, 0.5)

    resource = ResourceConfig(
        carrier_bandwidth_mhz=100.0, scs_khz=30, n_fft=4096, n_prbs=48,
        pdsch_n_symbols=10, dmrs_symbol_indices=[2, 7], dmrs_spacing_sc=SF, prg_size_rb=4,
    )
    grid = build_resource_grid(resource)
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    mcs = get_mcs("nr_256qam", 8, None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(tb.cb_k_values, tb.cb_e_values, num_iter=8, llr_clip=50.0)
    print(f"[0b] TB: {tb.tb_size} info bits, {len(tb.cb_k_values)} CBs, "
          f"{grid.n_data_re} data RE, grid {len(snr_bler_grid)} SNR pts", flush=True)

    Vs = [c.V for c in cands]
    tail_trials = max(trials, 3000)

    def run_point(snr_db: float, n_trials: int, skip: int) -> List[int]:
        """Run `n_trials` link trials starting after `skip` already-consumed trials
        (deterministic continuation of the same RNG stream)."""
        snr = 10.0 ** (snr_db / 10.0)
        noise_var = NT / snr
        errors = [0 for _ in cands]
        rng = np.random.default_rng([SEED_MC, 88, int(round(snr_db * 10))])
        for t in range(skip + n_trials):
            h = (rng.standard_normal(NT) + 1j * rng.standard_normal(NT)) / math.sqrt(2.0)
            payload = [rng.integers(0, 2, size=int(kk), dtype=np.int8) for kk in tb.cb_k_values]
            if t < skip:
                continue
            coded = np.concatenate(adapter.encode(payload))
            symbols = qam_modulate(coded, int(mcs.qm))
            noise = (rng.standard_normal(len(data_local)) + 1j * rng.standard_normal(len(data_local))) \
                * math.sqrt(noise_var / 2.0)
            llrs = []
            for V in Vs:
                g = (V @ h)[data_local][None, :]
                y = g * symbols[None, :] + noise[None, :]
                z, no_eff = equalize_mrc(y, g, noise_var)
                llrs.append(qam_demapper_maxlog(z, no_eff, int(mcs.qm)))
            decs = decode_same_tb_batch(adapter, llrs, payload)
            for i, dec in enumerate(decs):
                errors[i] += int(not dec.tb_success)
        return errors

    bler_rows = []
    zero_streak = 0
    for snr_db in snr_bler_grid:
        errors = run_point(float(snr_db), trials, 0)
        n_used = trials
        # densify the tail: low error counts need more trials for a stable 1% crossing
        if any(0 <= e <= 0.08 * trials for e in errors) and max(errors) > 0:
            extra = tail_trials - trials
            errors2 = run_point(float(snr_db), extra, trials)
            errors = [a + b for a, b in zip(errors, errors2)]
            n_used = tail_trials
        for i, cand in enumerate(cands):
            bler_rows.append({
                "id": cand.cid, "snr_db": float(snr_db), "trials": n_used,
                "tb_errors": errors[i], "bler": errors[i] / n_used,
            })
        print(f"[0b] snr={snr_db:g} dB trials={n_used} bler={[f'{e/n_used:.4f}' for e in errors]}", flush=True)
        write_csv(bler_rows, calib_dir / "calib0b_bler_curves.csv")
        zero_streak = zero_streak + 1 if max(errors) == 0 else 0
        if zero_streak >= 2:
            print("[0b] both candidates at zero errors twice; stopping grid early", flush=True)
            break

    def bler_snr_at(cid: str, target: float) -> float:
        pts = sorted([(r["snr_db"], max(r["bler"], 0.5 / r["trials"])) for r in bler_rows if r["id"] == cid])
        # enforce a non-increasing envelope (MC noise can create small upticks)
        blers = np.minimum.accumulate(np.array([b for _, b in pts]))
        snrs = np.array([s for s, _ in pts])
        for a in range(len(blers) - 1):
            if blers[a] >= target > blers[a + 1]:
                l0, l1 = math.log10(blers[a]), math.log10(blers[a + 1])
                return float(snrs[a] + (l0 - math.log10(target)) / (l0 - l1) * (snrs[a + 1] - snrs[a]))
        return float("nan")

    gap_rows = []
    all_pass = True
    for cand in cands:
        o10, o1 = targets[cand.cid]
        b10 = bler_snr_at(cand.cid, 0.10)
        b1 = bler_snr_at(cand.cid, 0.01)
        gap10 = b10 - o10
        gap1 = b1 - o1
        ok = np.isfinite(gap10) and np.isfinite(gap1) and abs(gap10 - gap1) < 0.4
        all_pass &= bool(ok)
        gap_rows.append({
            "id": cand.cid,
            "outage10_dB": round(o10, 4), "outage1_dB": round(o1, 4),
            "bler10_dB": round(b10, 4), "bler1_dB": round(b1, 4),
            "gap10_dB": round(gap10, 4), "gap1_dB": round(gap1, 4),
            "gap_diff_dB": round(abs(gap10 - gap1), 4),
            "pass": int(ok),
        })
        print(f"[0b] {cand.cid}: gap10={gap10:.3f} gap1={gap1:.3f} |diff|={abs(gap10-gap1):.3f} "
              f"{'PASS' if ok else 'FAIL'}", flush=True)
    write_csv(gap_rows, calib_dir / "calib0b_gaps.csv")
    curve_rows = []
    for cand in cands:
        crow = {"id": cand.cid}
        for s, p in zip(SNR_GRID_DB, curves[cand.cid]):
            crow[f"p{s:g}"] = f"{p:.6g}"
        curve_rows.append(crow)
    write_csv(curve_rows, calib_dir / "calib0b_outage_curves.csv")
    print(f"[0b] gate {'PASSED' if all_pass else 'FAILED'}", flush=True)
    return all_pass


# ----------------------------------------------------------------------------
# Stage: calib0b_refine — top up trials near the 1% crossing with an
# independent seed stream, merge counts, re-evaluate the parallelism gate.
# ----------------------------------------------------------------------------
def stage_calib0b_refine(out_dir: Path, target_trials: int = 12000) -> bool:
    from cdd_lls.core.config import ResourceConfig
    from cdd_lls.core.mcs import build_tb_layout, get_mcs
    from cdd_lls.phy.ldpc import SionnaLDPCAdapter
    from cdd_lls.phy.qam import qam_demapper_maxlog, qam_modulate
    from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
    from tools.run_v_design_piecewise_tradeoff import decode_same_tb_batch, equalize_mrc

    calib_dir = out_dir / "calibration"
    with (calib_dir / "calib0b_bler_curves.csv").open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    with (calib_dir / "calib0b_outage_curves.csv").open("r", encoding="utf-8", newline="") as f:
        curve_rows = list(csv.DictReader(f))

    cands = calib_candidates()
    cids = [c.cid for c in cands]
    by_point: Dict[float, Dict[str, Dict[str, float]]] = {}
    for r in rows:
        by_point.setdefault(float(r["snr_db"]), {})[r["id"]] = {
            "trials": int(r["trials"]), "tb_errors": int(r["tb_errors"]),
        }
    refine_snrs = sorted(
        s for s, d in by_point.items()
        if any(0.003 <= d[c]["tb_errors"] / d[c]["trials"] <= 0.03 for c in cids)
    )
    print(f"[0b-refine] refining {refine_snrs} to {target_trials} trials", flush=True)

    resource = ResourceConfig(
        carrier_bandwidth_mhz=100.0, scs_khz=30, n_fft=4096, n_prbs=48,
        pdsch_n_symbols=10, dmrs_symbol_indices=[2, 7], dmrs_spacing_sc=SF, prg_size_rb=4,
    )
    grid = build_resource_grid(resource)
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    mcs = get_mcs("nr_256qam", 8, None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(tb.cb_k_values, tb.cb_e_values, num_iter=8, llr_clip=50.0)
    Vs = [c.V for c in cands]

    for snr_db in refine_snrs:
        extra = target_trials - min(by_point[snr_db][c]["trials"] for c in cids)
        if extra <= 0:
            continue
        snr = 10.0 ** (snr_db / 10.0)
        noise_var = NT / snr
        errors = [0 for _ in cands]
        rng = np.random.default_rng([SEED_MC, 89, int(round(snr_db * 10))])
        t0 = time.time()
        for _ in range(extra):
            h = (rng.standard_normal(NT) + 1j * rng.standard_normal(NT)) / math.sqrt(2.0)
            payload = [rng.integers(0, 2, size=int(kk), dtype=np.int8) for kk in tb.cb_k_values]
            coded = np.concatenate(adapter.encode(payload))
            symbols = qam_modulate(coded, int(mcs.qm))
            noise = (rng.standard_normal(len(data_local)) + 1j * rng.standard_normal(len(data_local))) \
                * math.sqrt(noise_var / 2.0)
            llrs = []
            for V in Vs:
                g = (V @ h)[data_local][None, :]
                y = g * symbols[None, :] + noise[None, :]
                z, no_eff = equalize_mrc(y, g, noise_var)
                llrs.append(qam_demapper_maxlog(z, no_eff, int(mcs.qm)))
            decs = decode_same_tb_batch(adapter, llrs, payload)
            for i, dec in enumerate(decs):
                errors[i] += int(not dec.tb_success)
        for i, cid in enumerate(cids):
            by_point[snr_db][cid]["trials"] += extra
            by_point[snr_db][cid]["tb_errors"] += errors[i]
        merged = ["{}/{}".format(by_point[snr_db][c]["tb_errors"], by_point[snr_db][c]["trials"])
                  for c in cids]
        print(f"[0b-refine] snr={snr_db:g} +{extra} trials ({time.time()-t0:.0f}s) -> {merged}", flush=True)

    out_rows = []
    for snr_db in sorted(by_point):
        for cid in cids:
            d = by_point[snr_db][cid]
            out_rows.append({
                "id": cid, "snr_db": snr_db, "trials": d["trials"],
                "tb_errors": d["tb_errors"], "bler": d["tb_errors"] / d["trials"],
            })
    write_csv(out_rows, calib_dir / "calib0b_bler_curves.csv")

    def bler_snr_at(cid: str, target: float) -> float:
        pts = sorted([(r["snr_db"], max(r["bler"], 0.5 / r["trials"])) for r in out_rows
                      if r["id"] == cid])
        blers = np.minimum.accumulate(np.array([b for _, b in pts]))
        snrs = np.array([s for s, _ in pts])
        for a in range(len(blers) - 1):
            if blers[a] >= target > blers[a + 1]:
                l0, l1 = math.log10(blers[a]), math.log10(blers[a + 1])
                return float(snrs[a] + (l0 - math.log10(target)) / (l0 - l1) * (snrs[a + 1] - snrs[a]))
        return float("nan")

    targets = {}
    for crow in curve_rows:
        p = np.array([float(crow[f"p{s:g}"]) for s in SNR_GRID_DB])
        targets[crow["id"]] = (outage_snr_at(p, SNR_GRID_DB, 0.10),
                               outage_snr_at(p, SNR_GRID_DB, 0.01))
    gap_rows = []
    all_pass = True
    for cid in cids:
        o10, o1 = targets[cid]
        b10 = bler_snr_at(cid, 0.10)
        b1 = bler_snr_at(cid, 0.01)
        gap10, gap1 = b10 - o10, b1 - o1
        ok = np.isfinite(gap10) and np.isfinite(gap1) and abs(gap10 - gap1) < 0.4
        all_pass &= bool(ok)
        gap_rows.append({
            "id": cid,
            "outage10_dB": round(o10, 4), "outage1_dB": round(o1, 4),
            "bler10_dB": round(b10, 4), "bler1_dB": round(b1, 4),
            "gap10_dB": round(gap10, 4), "gap1_dB": round(gap1, 4),
            "gap_diff_dB": round(abs(gap10 - gap1), 4),
            "pass": int(ok),
        })
        print(f"[0b-refine] {cid}: gap10={gap10:.3f} gap1={gap1:.3f} |diff|={abs(gap10-gap1):.3f} "
              f"{'PASS' if ok else 'FAIL'}", flush=True)
    write_csv(gap_rows, calib_dir / "calib0b_gaps.csv")
    print(f"[0b-refine] gate {'PASSED' if all_pass else 'FAILED'}", flush=True)
    return all_pass


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True,
                        choices=["iqam", "scan", "calib0a", "calib0b", "calib0b_refine", "smoke"])
    parser.add_argument("--out", default="outputs/track_b_pilot_scan/20260709_main")
    parser.add_argument("--n-prbs", type=int, default=48)
    parser.add_argument("--n-mc", type=int, default=N_MC)
    parser.add_argument("--trials", type=int, default=400)
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    table = build_iqam_table(out_dir / "iqam_table.npz")

    if args.stage == "iqam":
        grid_db = np.arange(IQAM_DB_MIN, IQAM_DB_MAX + 1e-9, IQAM_DB_STEP)
        for q in (-20, -10, 0, 5, 10, 15, 20, 30, 35):
            i = int(round((q - IQAM_DB_MIN) / IQAM_DB_STEP))
            print(f"I_QAM({q:+d} dB) = {table[i]:.4f} bit")
    elif args.stage == "smoke":
        stage_scan(out_dir, args.n_prbs, args.n_mc, table, smoke=True)
    elif args.stage == "scan":
        stage_scan(out_dir, args.n_prbs, args.n_mc, table)
    elif args.stage == "calib0a":
        ok = stage_calib0a(out_dir, trials=300)
        sys.exit(0 if ok else 2)
    elif args.stage == "calib0b":
        ok = stage_calib0b(out_dir, table, trials=args.trials, n_mc=args.n_mc)
        sys.exit(0 if ok else 2)
    elif args.stage == "calib0b_refine":
        ok = stage_calib0b_refine(out_dir, target_trials=args.trials)
        sys.exit(0 if ok else 2)


if __name__ == "__main__":
    main()
