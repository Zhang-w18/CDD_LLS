"""Experiment 024: segment-aligned CE diagnostic and Sidon/QC link comparison."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
OUT_DEFAULT = ROOT / "outputs" / "experiment024_segment_sidon_qc" / "20260716_main"
NT = 8
SF = 24
SCS_HZ = 30e3
SEED = 20260716


def write_csv(rows: List[Dict[str, object]], path: Path) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: List[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def parse_float_list(text: str) -> List[float]:
    return [float(v.strip()) for v in text.split(",") if v.strip()]


def db10(x: float) -> float:
    return float(10.0 * np.log10(x))


def segment_nmse(V: np.ndarray, n_seg: int, tau_w_ns: float, snr_db: float) -> Tuple[float, int, int]:
    from tools.run_track_b_pilot_scan import make_estimator, mismatch_mse, uniform_pdp_cov

    K = V.shape[0]
    if K % n_seg:
        raise ValueError(f"K={K} not divisible by n_seg={n_seg}")
    L = K // n_seg
    snr = 10.0 ** (snr_db / 10.0)
    noise_var_ls = NT / snr / 2.0
    total_mse = 0.0
    total_power = 0.0
    pilot_counts: List[int] = []
    for s in range(n_seg):
        idx = np.arange(s * L, (s + 1) * L)
        pilots = idx[idx % SF == 0]
        data = idx[idx % SF != 0]
        pilot_counts.append(len(pilots))
        if not len(pilots) or not len(data):
            raise ValueError(f"empty pilot/data set in segment {s}: P={len(pilots)} D={len(data)}")
        rt_pp = uniform_pdp_cov(pilots, pilots, tau_w_ns)
        rt_dp = uniform_pdp_cov(data, pilots, tau_w_ns)
        A, _, _ = make_estimator(rt_pp, rt_dp, noise_var_ls)
        Vp, Vd = V[pilots, :], V[data, :]
        R_pp = Vp @ Vp.conj().T
        R_pd = Vp @ Vd.conj().T
        tr_dd = float(NT * len(data))
        total_mse += mismatch_mse(A, R_pp, R_pd, tr_dd, noise_var_ls)
        total_power += tr_dd
    return total_mse / total_power, min(pilot_counts), max(pilot_counts)


def diagnostic_h1(rows: List[Dict[str, object]], baseline: List[Dict[str, str]]) -> Dict[str, object]:
    passed: List[str] = []
    for row in rows:
        cnmse = float(row["r4_nmse_16dB"])
        cout = float(row["outage10_dB"])
        s1 = [float(b["nmse_R1_16dB"]) for b in baseline if float(b["outage10_dB"]) <= cout + 0.05]
        s2 = [float(b["outage10_dB"]) for b in baseline if float(b["nmse_R1_16dB"]) <= cnmse + 0.2]
        margin_i = min(s1) - cnmse if s1 else float("inf")
        margin_ii = min(s2) - cout if s2 else float("inf")
        ok = margin_i >= 1.0 or margin_ii >= 0.3
        row["margin_i_nmse_dB"] = margin_i
        row["margin_ii_outage_dB"] = margin_ii
        row["diagnostic_h1_pass"] = int(ok)
        if ok:
            passed.append(str(row["id"]))
    return {"n_pass": len(passed), "passed_ids": passed}


def run_segment(out_dir: Path) -> Dict[str, object]:
    from tools.run_track_b_pilot_scan import build_b3, build_b4

    all_rows: List[Dict[str, object]] = []
    by_bw_summary: Dict[str, object] = {}
    src_dir = ROOT / "outputs" / "track_b_pilot_scan" / "20260709_main"
    for n_prbs in (48, 24):
        K = n_prbs * 12
        source = read_csv(src_dir / f"candidates_{n_prbs}prb.csv")
        source_by_id = {r["id"]: r for r in source}
        cands = build_b3(K) + build_b4(K)
        bw_rows: List[Dict[str, object]] = []
        for cand in cands:
            n_seg = int(cand.params.get("n_seg", 8))
            r4, pmin, pmax = segment_nmse(cand.V, n_seg, 300.0, 16.0)
            src = source_by_id[cand.cid]
            r1_db = float(src["nmse_R1_16dB"])
            r4_db = db10(r4)
            row: Dict[str, object] = {
                "bandwidth_prb": n_prbs,
                "id": cand.cid,
                "family": cand.family,
                "n_seg": n_seg,
                "segment_length_sc": K // n_seg,
                "pilots_per_segment_min": pmin,
                "pilots_per_segment_max": pmax,
                "r1_nmse_16dB": r1_db,
                "r4_nmse_16dB": r4_db,
                "gain_r4_vs_r1_dB": r1_db - r4_db,
                "outage10_dB": float(src["outage10_dB"]),
                "outage1_dB": float(src["outage1_dB"]),
            }
            bw_rows.append(row)
        baseline = [r for r in source if r["family"] in ("B1_cdd", "B2_prg_cycling")]
        diag = diagnostic_h1(bw_rows, baseline)
        all_rows.extend(bw_rows)
        fam_summary = {}
        for family in ("B3_cc", "B4_cn"):
            vals = [float(r["gain_r4_vs_r1_dB"]) for r in bw_rows if r["family"] == family]
            fam_summary[family] = {
                "count": len(vals),
                "improved_count": int(sum(v > 0 for v in vals)),
                "gain_median_dB": float(np.median(vals)),
                "gain_max_dB": float(np.max(vals)),
                "gain_min_dB": float(np.min(vals)),
            }
        by_bw_summary[str(n_prbs)] = {"families": fam_summary, "diagnostic_h1": diag}
    write_csv(all_rows, out_dir / "segment_nmse.csv")
    summary = {"receiver": "R4 segment-aligned mismatched MMSE, tau_w=300 ns", "by_bandwidth": by_bw_summary}
    (out_dir / "segment_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    plot_segment(all_rows, out_dir / "figures" / "segment_nmse_gain.png")
    return summary


def plot_segment(rows: List[Dict[str, object]], path: Path) -> None:
    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.2), sharex=True, sharey=True)
    colors = {"B3_cc": "#d62728", "B4_cn": "#9467bd"}
    markers = {"B3_cc": "^", "B4_cn": "."}
    for ax, n_prbs in zip(axes, (48, 24)):
        subset = [r for r in rows if int(r["bandwidth_prb"]) == n_prbs]
        for fam in ("B3_cc", "B4_cn"):
            pts = [r for r in subset if r["family"] == fam]
            ax.scatter(
                [float(r["r1_nmse_16dB"]) for r in pts],
                [float(r["r4_nmse_16dB"]) for r in pts],
                s=35 if fam == "B3_cc" else 14,
                marker=markers[fam], c=colors[fam], alpha=0.75, label=fam,
            )
        lim = [-26, 4]
        ax.plot(lim, lim, "--", color="#555555", lw=1.0, label="R4 = R1")
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_title(f"{n_prbs} PRB")
        ax.set_xlabel("R1 4RB sliding NMSE @16 dB (dB)")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("R4 segment-aligned NMSE @16 dB (dB)")
    axes[0].legend(fontsize=8)
    fig.suptitle("Experiment 024 E1: segment-aligned receiver diagnostic (below diagonal is improvement)")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def cdd_V(K: int, j_list: Sequence[int]) -> np.ndarray:
    k = np.arange(K, dtype=np.float64)
    j = np.asarray(j_list, dtype=np.float64)
    return np.exp(-2j * math.pi * k[:, None] * j[None, :] / K)


def matched_matrix(V: np.ndarray, noise_var_ls: float) -> Tuple[np.ndarray, float]:
    pilots = np.arange(0, V.shape[0], SF)
    Vp = V[pilots, :]
    R_pp = Vp @ Vp.conj().T
    R_ap = V @ Vp.conj().T
    M = R_pp + noise_var_ls * np.eye(len(pilots))
    A = np.linalg.solve(M.conj().T, R_ap.conj().T).conj().T
    return A, float(np.linalg.cond(M))


def wilson(errors: int, trials: int, z: float = 1.96) -> Tuple[float, float]:
    if trials <= 0:
        return float("nan"), float("nan")
    p = errors / trials
    den = 1.0 + z * z / trials
    center = (p + z * z / (2 * trials)) / den
    half = z * math.sqrt(p * (1 - p) / trials + z * z / (4 * trials * trials)) / den
    return max(0.0, center - half), min(1.0, center + half)


def target_snr(rows: List[Dict[str, object]], cid: str, target: float) -> float:
    pts = sorted((float(r["snr_db"]), max(float(r["bler"]), 0.5 / int(r["trials"]))) for r in rows if r["id"] == cid)
    snrs = np.array([p[0] for p in pts])
    vals = np.minimum.accumulate(np.array([p[1] for p in pts]))
    for i in range(len(vals) - 1):
        if vals[i] >= target > vals[i + 1]:
            l0, l1 = math.log10(vals[i]), math.log10(vals[i + 1])
            f = (l0 - math.log10(target)) / (l0 - l1)
            return float(snrs[i] + f * (snrs[i + 1] - snrs[i]))
    return float("nan")


def run_link(out_dir: Path, snrs: Sequence[float], trials: int, batch_size: int = 1) -> Dict[str, object]:
    from cdd_lls.core.config import ResourceConfig
    from cdd_lls.core.mcs import build_tb_layout, get_mcs
    from cdd_lls.phy.ldpc import SionnaLDPCAdapter
    from cdd_lls.phy.qam import qam_demapper_maxlog, qam_modulate
    from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
    from tools.run_v_design_piecewise_tradeoff import decode_same_tb_batch, equalize_mrc

    K = 576
    designs = [
        ("QC_arith_s9", cdd_V(K, [9 * n for n in range(NT)])),
        ("Sidon", cdd_V(K, [0, 1, 3, 7, 12, 20, 30, 65])),
    ]
    resource = ResourceConfig(
        carrier_bandwidth_mhz=100.0, scs_khz=30, n_fft=4096, n_prbs=48,
        pdsch_n_symbols=10, dmrs_symbol_indices=[2, 7], dmrs_spacing_sc=SF, prg_size_rb=4,
    )
    grid = build_resource_grid(resource)
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    mcs = get_mcs("nr_256qam", 8, None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(tb.cb_k_values, tb.cb_e_values, num_iter=8, llr_clip=50.0)
    rows: List[Dict[str, object]] = []

    for snr_db in snrs:
        snr = 10.0 ** (float(snr_db) / 10.0)
        noise_var = NT / snr
        noise_var_ls = noise_var / 2.0
        estimators = {cid: matched_matrix(V, noise_var_ls) for cid, V in designs}
        errors = {cid: 0 for cid, _ in designs}
        nmse_sum = {cid: 0.0 for cid, _ in designs}
        t0 = time.time()
        for batch_start in range(1, trials + 1, batch_size):
            batch_end = min(batch_start + batch_size - 1, trials)
            payload_rng = np.random.default_rng([SEED, int(round(snr_db * 100)), batch_start, 777])
            payload = [payload_rng.integers(0, 2, size=int(k), dtype=np.int8) for k in tb.cb_k_values]
            coded = np.concatenate(adapter.encode(payload))
            symbols = qam_modulate(coded, int(mcs.qm))
            llrs = []
            ids = []
            for trial in range(batch_start, batch_end + 1):
                rng = np.random.default_rng([SEED, int(round(snr_db * 100)), trial])
                h = (rng.standard_normal(NT) + 1j * rng.standard_normal(NT)) / math.sqrt(2.0)
                pilot_noise = (rng.standard_normal(len(pilot_local)) + 1j * rng.standard_normal(len(pilot_local))) * math.sqrt(noise_var_ls / 2.0)
                data_noise = (rng.standard_normal(len(data_local)) + 1j * rng.standard_normal(len(data_local))) * math.sqrt(noise_var / 2.0)
                for cid, V in designs:
                    g = V @ h
                    A, _ = estimators[cid]
                    g_hat = A @ (g[pilot_local] + pilot_noise)
                    nmse_sum[cid] += float(np.sum(np.abs(g_hat - g) ** 2) / np.sum(np.abs(g) ** 2))
                    y = g[data_local] * symbols + data_noise
                    z, no_eff = equalize_mrc(y[None, :], g_hat[data_local][None, :], noise_var)
                    llrs.append(qam_demapper_maxlog(z, no_eff, int(mcs.qm)))
                    ids.append(cid)
            decoded = decode_same_tb_batch(adapter, llrs, payload)
            for cid, result in zip(ids, decoded):
                errors[cid] += int(not result.tb_success)
            if batch_end % 50 < batch_size or batch_end == trials:
                print(f"[E2] snr={snr_db:g} trial={batch_end}/{trials} errors={errors} elapsed={time.time()-t0:.1f}s", flush=True)
        for cid, _ in designs:
            lo, hi = wilson(errors[cid], trials)
            rows.append({
                "id": cid, "snr_db": float(snr_db), "trials": trials,
                "tb_errors": errors[cid], "bler": errors[cid] / trials,
                "bler_wilson95_lo": lo, "bler_wilson95_hi": hi,
                "ce_nmse_mean": nmse_sum[cid] / trials,
                "ce_nmse_mean_dB": db10(nmse_sum[cid] / trials),
                "estimator_cond": estimators[cid][1],
            })
        write_csv(rows, out_dir / "sidon_qc_bler.csv")
    summary: Dict[str, object] = {"trials_per_snr": trials, "batch_size": batch_size, "snrs_dB": list(snrs), "targets": {}}
    for cid, _ in designs:
        summary["targets"][cid] = {"bler10_snr_dB": target_snr(rows, cid, 0.10), "bler1_snr_dB": target_snr(rows, cid, 0.01)}
    q10 = summary["targets"]["QC_arith_s9"]["bler10_snr_dB"]
    s10 = summary["targets"]["Sidon"]["bler10_snr_dB"]
    q1 = summary["targets"]["QC_arith_s9"]["bler1_snr_dB"]
    s1 = summary["targets"]["Sidon"]["bler1_snr_dB"]
    summary["sidon_gain_dB"] = {"bler10": q10 - s10, "bler1": q1 - s1}
    (out_dir / "sidon_qc_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    plot_link(rows, out_dir / "figures" / "sidon_qc_bler.png")
    return summary


def plot_link(rows: List[Dict[str, object]], path: Path) -> None:
    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    styles = {"QC_arith_s9": ("o-", "#f77f00"), "Sidon": ("s-", "#111111")}
    for cid, (style, color) in styles.items():
        pts = sorted((float(r["snr_db"]), float(r["bler"]), float(r["bler_wilson95_lo"]), float(r["bler_wilson95_hi"])) for r in rows if r["id"] == cid)
        x = np.array([p[0] for p in pts]); y = np.array([max(p[1], 0.5 / int(next(r["trials"] for r in rows if r["id"] == cid))) for p in pts])
        lo = np.array([max(p[2], 1e-4) for p in pts]); hi = np.array([p[3] for p in pts])
        ax.semilogy(x, y, style, color=color, label=cid, lw=1.7, ms=5)
        ax.fill_between(x, lo, hi, color=color, alpha=0.12)
    ax.axhline(0.1, color="#555555", ls=":", lw=1)
    ax.axhline(0.01, color="#555555", ls=":", lw=1)
    ax.set_xlabel("average receive SNR (dB)")
    ax.set_ylabel("estimated-CSI TB BLER")
    ax.set_title("Experiment 024 E2: V-aware matched LMMSE, 48 PRB flat channel")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("segment", "link", "all"), default="all")
    parser.add_argument("--out", type=Path, default=OUT_DEFAULT)
    parser.add_argument("--snrs", default="13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18")
    parser.add_argument("--trials", type=int, default=400)
    parser.add_argument("--batch-size", type=int, default=1)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.stage in ("segment", "all"):
        print(json.dumps(run_segment(args.out), indent=2), flush=True)
    if args.stage in ("link", "all"):
        print(json.dumps(run_link(args.out, parse_float_list(args.snrs), args.trials, args.batch_size), indent=2), flush=True)


if __name__ == "__main__":
    main()
