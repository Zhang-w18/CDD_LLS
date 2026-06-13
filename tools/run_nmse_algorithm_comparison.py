from __future__ import annotations

import argparse
import csv
import json
import math
import os
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import sys
from typing import Dict, Iterable, List, Optional

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ChannelConfig, ResourceConfig
from cdd_lls.phy.channel_tdl import generate_tdl_channel, make_exponential_pdp
from cdd_lls.phy.estimators import covariance_matrix, nmse, shifted_pdp
from cdd_lls.phy.precoding import cdd_equivalent_from_branches
from cdd_lls.phy.resource_grid import ResourceGrid, build_resource_grid, local_indices_for_subcarriers


@dataclass(frozen=True)
class CaseConfig:
    case_id: str
    delay_spread_ns: float
    cdd_delay: int
    dmrs_spacing_sc: int


@dataclass(frozen=True)
class LinearEstimator:
    label: str
    family: str
    matrix: np.ndarray
    metadata: Dict[str, object]


@dataclass(frozen=True)
class BasisEstimator:
    label: str
    support: np.ndarray
    matrix: np.ndarray
    freq_basis: np.ndarray
    cond_number: float
    effective_rank: float
    metadata: Dict[str, object]


@dataclass(frozen=True)
class BlockEstimator:
    label: str
    blocks: List[Dict[str, object]]
    block_width_sc: int
    cond_max: float
    cond_mean: float
    metadata: Dict[str, object]


def snr_points(spec: List[float]) -> List[float]:
    start, stop, step = [float(x) for x in spec]
    out = []
    x = start
    while x <= stop + 1e-9:
        out.append(round(x, 6))
        x += step
    return out


def support_indices_energy(pdp: np.ndarray, threshold: float) -> np.ndarray:
    c = np.cumsum(np.asarray(pdp, dtype=np.float64))
    cutoff = int(np.searchsorted(c, float(threshold), side="left")) + 1
    cutoff = min(max(cutoff, 1), len(pdp))
    return np.arange(cutoff, dtype=np.int64)


def freq_corr_from_pdp(pdp: np.ndarray, n_fft: int, max_delta_sc: int) -> np.ndarray:
    p = np.asarray(pdp, dtype=np.float64)
    taps = np.arange(len(p), dtype=np.float64)
    deltas = np.arange(int(max_delta_sc) + 1, dtype=np.float64)
    r = np.exp(-1j * 2.0 * np.pi * deltas[:, None] * taps[None, :] / float(n_fft)) @ p
    return np.abs(r / max(float(np.sum(p)), 1e-30))


def first_below(rho: np.ndarray, threshold: float) -> Optional[int]:
    idx = np.where(np.asarray(rho) <= float(threshold))[0]
    if len(idx) == 0:
        return None
    return int(idx[0])


def lmmse_matrix(
    target_k: np.ndarray,
    pilot_k: np.ndarray,
    pdp: np.ndarray,
    n_fft: int,
    noise_var: float,
    loading: float,
) -> np.ndarray:
    R_pp = covariance_matrix(pilot_k, pilot_k, pdp, n_fft)
    R_tp = covariance_matrix(target_k, pilot_k, pdp, n_fft)
    A = R_pp + (float(noise_var) + float(loading)) * np.eye(len(pilot_k), dtype=np.complex128)
    return R_tp @ np.linalg.solve(A, np.eye(len(pilot_k), dtype=np.complex128))


def make_direct_estimator(
    grid: ResourceGrid,
    pdp: np.ndarray,
    delays: List[int],
    noise_var_ls: float,
    loading: float,
) -> LinearEstimator:
    assumed = shifted_pdp(pdp, delays)
    W = lmmse_matrix(
        target_k=grid.subcarrier_indices,
        pilot_k=grid.pilot_subcarriers,
        pdp=assumed,
        n_fft=grid.n_fft,
        noise_var=noise_var_ls,
        loading=loading,
    )
    return LinearEstimator(
        label=f"Alg1 Direct RMMSE WB {grid.n_sc}sc",
        family="alg1_direct_rmmse",
        matrix=W,
        metadata={
            "processing_bandwidth_sc": int(grid.n_sc),
            "processing_bandwidth_rb": float(grid.n_sc / 12.0),
            "pilot_count": int(grid.pilot_count),
            "covariance": "CDD-shifted PDP",
        },
    )


def make_basis_estimator(
    grid: ResourceGrid,
    pdp: np.ndarray,
    delays: List[int],
    noise_var_ls: float,
    loading: float,
    energy_threshold: float,
) -> BasisEstimator:
    n_tx = len(delays)
    support = support_indices_energy(pdp, energy_threshold)
    p_support = np.asarray(pdp, dtype=np.float64)[support]
    pilot_k = grid.pilot_subcarriers.astype(np.float64)
    d = np.asarray(delays, dtype=np.float64)
    n_cols = n_tx * len(support)
    Phi = np.empty((len(pilot_k), n_cols), dtype=np.complex128)
    col = 0
    for m in range(n_tx):
        shifted = support.astype(np.float64) + d[m]
        Phi[:, col:col + len(support)] = np.exp(
            -1j * 2.0 * np.pi * pilot_k[:, None] * shifted[None, :] / float(grid.n_fft)
        ) / math.sqrt(float(n_tx))
        col += len(support)

    rdiag = np.tile(p_support, n_tx)
    PhiR = Phi * rdiag[None, :]
    S = PhiR @ Phi.conj().T
    S += (float(noise_var_ls) + float(loading)) * np.eye(S.shape[0], dtype=np.complex128)
    B = rdiag[:, None] * Phi.conj().T
    M = B @ np.linalg.solve(S, np.eye(S.shape[0], dtype=np.complex128))

    target_k = grid.subcarrier_indices.astype(np.float64)
    F = np.exp(-1j * 2.0 * np.pi * target_k[:, None] * support[None, :] / float(grid.n_fft))
    weighted = Phi * np.sqrt(np.maximum(rdiag, 0.0))[None, :]
    svals = np.linalg.svd(weighted, compute_uv=False)
    eps = 1e-12
    cond = float(svals[0] / max(svals[-1], eps)) if len(svals) else float("nan")
    eff_rank = float((np.sum(svals) ** 2) / max(np.sum(svals ** 2), eps)) if len(svals) else float("nan")
    return BasisEstimator(
        label=f"Alg2 Basis LMMSE E{int(round(energy_threshold * 100))}",
        support=support,
        matrix=M,
        freq_basis=F,
        cond_number=cond,
        effective_rank=eff_rank,
        metadata={
            "basis_energy_threshold": float(energy_threshold),
            "support_len": int(len(support)),
            "n_unknown_taps": int(n_cols),
            "pilot_count": int(grid.pilot_count),
        },
    )


def apply_basis(est: BasisEstimator, ls_obs: np.ndarray, grid: ResourceGrid, delays: List[int]) -> tuple[np.ndarray, np.ndarray]:
    n_rx = int(ls_obs.shape[0])
    n_tx = len(delays)
    taps = ls_obs @ est.matrix.T
    h_hat = np.empty((n_rx, n_tx, grid.n_sc), dtype=np.complex128)
    col = 0
    for m in range(n_tx):
        taps_m = taps[:, col:col + len(est.support)]
        h_hat[:, m, :] = taps_m @ est.freq_basis.T
        col += len(est.support)
    return cdd_equivalent_from_branches(h_hat, grid, delays), h_hat


def make_block_estimator(
    grid: ResourceGrid,
    pdp: np.ndarray,
    delays: List[int],
    noise_var_ls: float,
    loading: float,
    threshold: float,
) -> BlockEstimator:
    rho = freq_corr_from_pdp(pdp, grid.n_fft, grid.n_sc - 1)
    bc = first_below(rho, threshold)
    block_width = grid.n_sc if bc is None else max(1, min(grid.n_sc, int(bc)))
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    delays_arr = np.asarray(delays, dtype=np.float64)
    n_tx = len(delays)
    blocks: List[Dict[str, object]] = []
    conds = []
    for start in range(0, grid.n_sc, block_width):
        stop = min(grid.n_sc, start + block_width)
        target_local = np.arange(start, stop, dtype=np.int64)
        mask = (pilot_local >= start) & (pilot_local < stop)
        pilot_idx = np.where(mask)[0]
        pilot_k = grid.pilot_subcarriers[pilot_idx].astype(np.float64)
        if len(pilot_k) == 0:
            center = 0.5 * (grid.subcarrier_indices[start] + grid.subcarrier_indices[stop - 1])
            nearest = int(np.argmin(np.abs(grid.pilot_subcarriers.astype(np.float64) - center)))
            pilot_idx = np.asarray([nearest], dtype=np.int64)
            pilot_k = grid.pilot_subcarriers[pilot_idx].astype(np.float64)
        A = np.exp(
            -1j * 2.0 * np.pi * pilot_k[:, None] * delays_arr[None, :] / float(grid.n_fft)
        ) / math.sqrt(float(n_tx))
        S = A @ A.conj().T
        S += (float(noise_var_ls) + float(loading)) * np.eye(S.shape[0], dtype=np.complex128)
        M = A.conj().T @ np.linalg.solve(S, np.eye(S.shape[0], dtype=np.complex128))
        cond = float(np.linalg.cond(A)) if len(pilot_k) >= n_tx else float("inf")
        conds.append(cond)
        blocks.append({
            "start": int(start),
            "stop": int(stop),
            "target_local": target_local,
            "pilot_idx": pilot_idx,
            "matrix": M,
            "cond": cond,
            "pilot_count": int(len(pilot_idx)),
        })
    finite = [x for x in conds if np.isfinite(x)]
    return BlockEstimator(
        label=f"Alg2 CoherenceBlock Bc{threshold:.1f}",
        blocks=blocks,
        block_width_sc=int(block_width),
        cond_max=float(max(finite)) if finite else float("inf"),
        cond_mean=float(np.mean(finite)) if finite else float("inf"),
        metadata={
            "bc_threshold": float(threshold),
            "block_width_sc": int(block_width),
            "block_width_mhz": float(block_width * grid.scs_khz / 1000.0),
            "n_blocks": int(len(blocks)),
        },
    )


def apply_block(est: BlockEstimator, ls_obs: np.ndarray, grid: ResourceGrid, delays: List[int]) -> tuple[np.ndarray, np.ndarray]:
    n_rx = int(ls_obs.shape[0])
    n_tx = len(delays)
    h_hat = np.empty((n_rx, n_tx, grid.n_sc), dtype=np.complex128)
    for block in est.blocks:
        pilot_idx = block["pilot_idx"]
        M = block["matrix"]
        target_local = block["target_local"]
        h_block = ls_obs[:, pilot_idx] @ M.T
        h_hat[:, :, target_local] = h_block[:, :, None]
    return cdd_equivalent_from_branches(h_hat, grid, delays), h_hat


def make_port_estimators(
    grid: ResourceGrid,
    pdp: np.ndarray,
    n_tx: int,
    noise_var_ls: float,
    loading: float,
    mode: str,
) -> List[LinearEstimator]:
    out = []
    all_pilots = grid.pilot_subcarriers
    for m in range(n_tx):
        if mode == "equal_total_overhead":
            pilot_k = all_pilots[m::n_tx]
        elif mode == "equal_per_port_density":
            pilot_k = all_pilots
        else:
            raise ValueError(f"Unsupported port DMRS mode={mode}")
        W = lmmse_matrix(
            target_k=grid.subcarrier_indices,
            pilot_k=pilot_k,
            pdp=pdp,
            n_fft=grid.n_fft,
            noise_var=noise_var_ls,
            loading=loading,
        )
        out.append(LinearEstimator(
            label=f"port{m}",
            family="alg3_port_dmrs",
            matrix=W,
            metadata={
                "port": int(m),
                "pilot_count": int(len(pilot_k)),
                "pilot_k": pilot_k.astype(int).tolist(),
            },
        ))
    return out


def apply_port_dmrs(
    port_estimators: List[LinearEstimator],
    H: np.ndarray,
    grid: ResourceGrid,
    delays: List[int],
    noise_var_ls: float,
    rng: np.random.Generator,
    mode: str,
) -> tuple[np.ndarray, np.ndarray]:
    n_rx, n_tx, _ = H.shape
    h_hat = np.empty_like(H)
    all_pilots = grid.pilot_subcarriers
    for m, est in enumerate(port_estimators):
        pilot_k = np.asarray(est.metadata["pilot_k"], dtype=np.int64)
        pilot_local = local_indices_for_subcarriers(grid, pilot_k)
        noise = (
            rng.normal(size=(n_rx, len(pilot_local)))
            + 1j * rng.normal(size=(n_rx, len(pilot_local)))
        ) * math.sqrt(float(noise_var_ls) / 2.0)
        obs = H[:, m, pilot_local] + noise
        h_hat[:, m, :] = obs @ est.matrix.T
    return cdd_equivalent_from_branches(h_hat, grid, delays), h_hat


def write_csv(rows: List[Dict[str, object]], path: Path) -> None:
    if not rows:
        return
    fields = sorted(set().union(*(row.keys() for row in rows)))
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot_nmse(rows: List[Dict[str, object]], cases: List[str], metric: str, path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    algs = []
    for row in rows:
        label = str(row["algorithm"])
        if label not in algs:
            algs.append(label)
    fig, axes = plt.subplots(len(cases), 1, figsize=(10.5, max(3.1 * len(cases), 5.0)), sharex=True)
    if len(cases) == 1:
        axes = [axes]
    for ax, case_id in zip(axes, cases):
        for alg in algs:
            pts = [r for r in rows if r["case_id"] == case_id and r["algorithm"] == alg]
            pts = sorted(pts, key=lambda r: float(r["snr_db"]))
            ys = [float(r[metric]) for r in pts]
            if not pts or all(not np.isfinite(y) for y in ys):
                continue
            ax.semilogy([float(r["snr_db"]) for r in pts], ys, marker="o", markersize=3.2, linewidth=1.2, label=alg)
        ax.set_title(case_id, loc="left", fontsize=10)
        ax.set_ylabel(metric.replace("_", " "))
        ax.grid(True, which="both", alpha=0.28)
    axes[-1].set_xlabel("SNR (dB)")
    axes[0].legend(fontsize=7, ncol=2, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def mean_or_nan(xs: Iterable[float]) -> float:
    vals = [float(x) for x in xs if np.isfinite(float(x))]
    return float(np.mean(vals)) if vals else float("nan")


def median_or_nan(xs: Iterable[float]) -> float:
    vals = [float(x) for x in xs if np.isfinite(float(x))]
    return float(np.median(vals)) if vals else float("nan")


def run(args: argparse.Namespace) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out) / f"nmse_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        CaseConfig("flat_ds1_d64_sf2", 1.0, 64, 2),
        CaseConfig("ds30_d16_sf6", 30.0, 16, 6),
        CaseConfig("ds30_d32_sf6", 30.0, 32, 6),
        CaseConfig("ds100_d16_sf6", 100.0, 16, 6),
        CaseConfig("ds100_d32_sf6", 100.0, 32, 6),
        CaseConfig("ds30_d64_sf2", 30.0, 64, 2),
    ]
    snrs = snr_points([float(args.snr_start), float(args.snr_stop), float(args.snr_step)])
    rows: List[Dict[str, object]] = []
    detail_rows: List[Dict[str, object]] = []
    rng_master = np.random.default_rng(int(args.seed))

    for case in cases:
        resource = ResourceConfig(
            carrier_bandwidth_mhz=100.0,
            scs_khz=30,
            n_fft=4096,
            n_prbs=48,
            pdsch_n_symbols=10,
            dmrs_symbol_indices=[2, 7],
            dmrs_spacing_sc=int(case.dmrs_spacing_sc),
            prg_size_rb=4,
        )
        grid = build_resource_grid(resource)
        sample_period_ns = 1e9 / (float(resource.n_fft) * float(resource.scs_khz) * 1e3)
        pdp = make_exponential_pdp(case.delay_spread_ns, sample_period_ns, max_delay_factor=8.0)
        delays = [0, int(case.cdd_delay)]
        channel = ChannelConfig(delay_spread_ns=case.delay_spread_ns, max_delay_factor=8.0)
        rho = freq_corr_from_pdp(pdp, grid.n_fft, grid.n_sc - 1)
        bc09 = first_below(rho, 0.9)
        bc05 = first_below(rho, 0.5)
        case_seed = int(rng_master.integers(0, 2**31 - 1))

        for snr_db in snrs:
            noise_var = 10.0 ** (-float(snr_db) / 10.0)
            noise_var_ls = noise_var / float(len(resource.dmrs_symbol_indices))
            direct = make_direct_estimator(grid, pdp, delays, noise_var_ls, float(args.loading))
            basis_ests = [
                make_basis_estimator(grid, pdp, delays, noise_var_ls, float(args.loading), thr)
                for thr in (0.90, 0.95, 0.99)
            ]
            block = make_block_estimator(grid, pdp, delays, noise_var_ls, float(args.loading), 0.9)
            port_total = make_port_estimators(
                grid, pdp, n_tx=2, noise_var_ls=noise_var_ls, loading=float(args.loading), mode="equal_total_overhead"
            )
            port_per = make_port_estimators(
                grid, pdp, n_tx=2, noise_var_ls=noise_var_ls, loading=float(args.loading), mode="equal_per_port_density"
            )

            accum: Dict[str, Dict[str, List[float] | Dict[str, object]]] = {}

            def add_result(label: str, family: str, eff: float, branch: float, meta: Dict[str, object]) -> None:
                if label not in accum:
                    accum[label] = {"family": {"value": family}, "eff": [], "branch": [], "meta": meta}
                accum[label]["eff"].append(float(eff))
                accum[label]["branch"].append(float(branch))

            trial_rng = np.random.default_rng(case_seed + int(round((snr_db + 100.0) * 1000.0)))
            pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
            for _ in range(int(args.trials)):
                realization = generate_tdl_channel(trial_rng, grid, channel, n_tx=2, n_rx=4)
                H = realization.H
                g = cdd_equivalent_from_branches(H, grid, delays)
                noise = (
                    trial_rng.normal(size=(4, len(pilot_local)))
                    + 1j * trial_rng.normal(size=(4, len(pilot_local)))
                ) * math.sqrt(noise_var_ls / 2.0)
                ls_cdd = g[:, pilot_local] + noise

                g_hat = ls_cdd @ direct.matrix.T
                add_result(direct.label, direct.family, nmse(g, g_hat), float("nan"), direct.metadata)

                for best in basis_ests:
                    g_hat, h_hat = apply_basis(best, ls_cdd, grid, delays)
                    meta = dict(best.metadata)
                    meta.update({"cond_number": best.cond_number, "effective_rank": best.effective_rank})
                    add_result(best.label, "alg2_basis_reduced_support", nmse(g, g_hat), nmse(H, h_hat), meta)

                g_hat, h_hat = apply_block(block, ls_cdd, grid, delays)
                meta = dict(block.metadata)
                meta.update({"cond_number": block.cond_max, "cond_mean": block.cond_mean})
                add_result(block.label, "alg2_coherence_block", nmse(g, g_hat), nmse(H, h_hat), meta)

                g_hat, h_hat = apply_port_dmrs(
                    port_total, H, grid, delays, noise_var_ls, trial_rng, mode="equal_total_overhead"
                )
                add_result(
                    "Alg3 Port-DMRS RMMSE equal-total",
                    "alg3_port_dmrs",
                    nmse(g, g_hat),
                    nmse(H, h_hat),
                    {
                        "overhead_mode": "equal_total_overhead",
                        "total_pilot_re_relative": 1.0,
                        "pilot_count_per_port": int(len(grid.pilot_subcarriers[0::2])),
                    },
                )

                g_hat, h_hat = apply_port_dmrs(
                    port_per, H, grid, delays, noise_var_ls, trial_rng, mode="equal_per_port_density"
                )
                add_result(
                    "Alg3 Port-DMRS RMMSE equal-per-port",
                    "alg3_port_dmrs",
                    nmse(g, g_hat),
                    nmse(H, h_hat),
                    {
                        "overhead_mode": "equal_per_port_density",
                        "total_pilot_re_relative": 2.0,
                        "pilot_count_per_port": int(grid.pilot_count),
                    },
                )

            for label, values in accum.items():
                meta = dict(values["meta"])
                eff_vals = values["eff"]
                branch_vals = values["branch"]
                row = {
                    "case_id": case.case_id,
                    "algorithm": label,
                    "family": values["family"]["value"],
                    "snr_db": float(snr_db),
                    "trials": int(args.trials),
                    "ce_nmse_eff_mean": mean_or_nan(eff_vals),
                    "ce_nmse_eff_median": median_or_nan(eff_vals),
                    "ce_nmse_branch_mean": mean_or_nan(branch_vals),
                    "ce_nmse_branch_median": median_or_nan(branch_vals),
                    "delay_spread_ns": float(case.delay_spread_ns),
                    "cdd_delay_samples": int(case.cdd_delay),
                    "cdd_delay_ns": float(case.cdd_delay) * sample_period_ns,
                    "dmrs_spacing_sc": int(case.dmrs_spacing_sc),
                    "n_prbs": int(resource.n_prbs),
                    "n_sc": int(grid.n_sc),
                    "pilot_count_combined": int(grid.pilot_count),
                    "dmrs_symbols": int(len(resource.dmrs_symbol_indices)),
                    "noise_var_ls": float(noise_var_ls),
                    "bc09_sc_original": "" if bc09 is None else int(bc09),
                    "bc05_sc_original": "" if bc05 is None else int(bc05),
                }
                row.update(meta)
                rows.append(row)
                detail_rows.append({
                    "case_id": case.case_id,
                    "algorithm": label,
                    "snr_db": float(snr_db),
                    "metadata": json.dumps(meta, ensure_ascii=False, sort_keys=True),
                })

    write_csv(rows, out_dir / "nmse_summary.csv")
    write_csv(detail_rows, out_dir / "nmse_metadata.csv")
    with (out_dir / "run_config.json").open("w", encoding="utf-8") as f:
        json.dump({
            "snr_db": snrs,
            "trials": int(args.trials),
            "seed": int(args.seed),
            "note": "NMSE-only channel-estimation comparison; no LDPC encoding/decoding is run.",
        }, f, indent=2)

    case_ids = [c.case_id for c in cases]
    plot_nmse(rows, case_ids, "ce_nmse_eff_mean", fig_dir / "experiment9_nmse_effective.png")
    plot_nmse(rows, case_ids, "ce_nmse_branch_mean", fig_dir / "experiment9_nmse_branch.png")
    print(out_dir)
    return out_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/nmse_algorithm_comparison")
    parser.add_argument("--fig-dir", default="docs/figures")
    parser.add_argument("--snr-start", type=float, default=-6.0)
    parser.add_argument("--snr-stop", type=float, default=8.0)
    parser.add_argument("--snr-step", type=float, default=2.0)
    parser.add_argument("--trials", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260611)
    parser.add_argument("--loading", type=float, default=1e-8)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
