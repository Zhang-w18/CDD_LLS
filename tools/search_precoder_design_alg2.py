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
from typing import Dict, Iterable, List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ChannelConfig, ResourceConfig
from cdd_lls.phy.channel_tdl import generate_tdl_channel, make_exponential_pdp
from cdd_lls.phy.estimators import nmse
from cdd_lls.phy.precoding import equivalent_channel
from cdd_lls.phy.resource_grid import ResourceGrid, build_resource_grid, local_indices_for_subcarriers
from tools.run_nmse_algorithm_comparison import first_below, freq_corr_from_pdp, make_direct_estimator, mean_or_nan, median_or_nan


@dataclass(frozen=True)
class PrecoderDesign:
    label: str
    family: str
    C: np.ndarray
    metadata: Dict[str, object]


@dataclass(frozen=True)
class LocalEstimator:
    label: str
    blocks: List[Dict[str, object]]
    cond_mean: float
    cond_max: float
    metadata: Dict[str, object]


def parse_float_list(text: str) -> List[float]:
    return [float(x.strip()) for x in str(text).split(",") if x.strip()]


def parse_int_list(text: str) -> List[int]:
    return [int(float(x.strip())) for x in str(text).split(",") if x.strip()]


def stable_seed(*values: object, base: int) -> int:
    h = int(base) & 0x7FFFFFFF
    for value in values:
        for ch in str(value):
            h = (h * 131 + ord(ch)) & 0x7FFFFFFF
    return h


def write_csv(rows: List[Dict[str, object]], path: Path) -> None:
    if not rows:
        return
    fields = sorted(set().union(*(row.keys() for row in rows)))
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def gain_db(reference: float, candidate: float) -> float:
    if not np.isfinite(reference) or not np.isfinite(candidate) or reference <= 0.0 or candidate <= 0.0:
        return float("nan")
    return float(10.0 * np.log10(reference / candidate))


def two_tx_precoder_from_relative_phase(z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, dtype=np.complex128).reshape(-1)
    C = np.empty((len(z), 2), dtype=np.complex128)
    C[:, 0] = 1.0 / math.sqrt(2.0)
    C[:, 1] = z / math.sqrt(2.0)
    return C


def phase_by_pilot_bins(grid: ResourceGrid, pilot_phases: np.ndarray) -> np.ndarray:
    pilot_phases = np.asarray(pilot_phases, dtype=np.complex128)
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    local = np.arange(grid.n_sc, dtype=np.int64)
    idx = np.searchsorted(pilot_local, local, side="right") - 1
    idx = np.clip(idx, 0, len(pilot_phases) - 1)
    return pilot_phases[idx]


def cdd_design(grid: ResourceGrid, delay: int) -> PrecoderDesign:
    z = np.exp(-1j * 2.0 * np.pi * grid.subcarrier_indices.astype(np.float64) * float(delay) / float(grid.n_fft))
    return PrecoderDesign(
        label=f"CDD d={int(delay)}",
        family="cdd_linear_phase",
        C=two_tx_precoder_from_relative_phase(z),
        metadata={"cdd_delay_samples": int(delay)},
    )


def constant_design(grid: ResourceGrid) -> PrecoderDesign:
    return PrecoderDesign(
        label="CONST z=1",
        family="constant_rank1",
        C=two_tx_precoder_from_relative_phase(np.ones(grid.n_sc, dtype=np.complex128)),
        metadata={},
    )


def pilot_cycle_design(grid: ResourceGrid, phases: np.ndarray, label: str) -> PrecoderDesign:
    phases = np.asarray(phases, dtype=np.complex128)
    pilot_phases = phases[np.arange(grid.pilot_count) % len(phases)]
    z = phase_by_pilot_bins(grid, pilot_phases)
    return PrecoderDesign(
        label=label,
        family="pilot_phase_cycling",
        C=two_tx_precoder_from_relative_phase(z),
        metadata={"cycle": [complex(x) for x in phases]},
    )


def rb_qpsk_design(grid: ResourceGrid, rb_group: int, label: str) -> PrecoderDesign:
    phases = np.asarray([1.0, 1.0j, -1.0, -1.0j], dtype=np.complex128)
    z = np.empty(grid.n_sc, dtype=np.complex128)
    for local in range(grid.n_sc):
        rb = int(local // 12)
        idx = int((rb // max(int(rb_group), 1)) % len(phases))
        z[local] = phases[idx]
    return PrecoderDesign(
        label=label,
        family="rb_qpsk_cycling",
        C=two_tx_precoder_from_relative_phase(z),
        metadata={"rb_group": int(rb_group)},
    )


def local_condition_for_pilot_phases(
    grid: ResourceGrid,
    pilot_phases: np.ndarray,
    window_sc: int,
) -> Tuple[float, float, float, int]:
    C = two_tx_precoder_from_relative_phase(phase_by_pilot_bins(grid, pilot_phases))
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    conds: List[float] = []
    bad = 0
    for start in range(0, grid.n_sc, int(window_sc)):
        stop = min(grid.n_sc, start + int(window_sc))
        mask = (pilot_local >= start) & (pilot_local < stop)
        pilot_idx = np.where(mask)[0]
        if len(pilot_idx) < 2:
            bad += 1
            continue
        A = C[pilot_local[pilot_idx], :]
        conds.append(float(np.linalg.cond(A)))
    if not conds:
        return float("inf"), float("inf"), 1.0, bad
    return float(np.mean(conds)), float(np.max(conds)), float(bad / max(math.ceil(grid.n_sc / int(window_sc)), 1)), bad


def optimized_qpsk_design(grid: ResourceGrid, window_sc: int, n_candidates: int, seed: int) -> PrecoderDesign:
    phases = np.asarray([1.0, 1.0j, -1.0, -1.0j], dtype=np.complex128)
    candidates: List[np.ndarray] = []
    idx = np.arange(grid.pilot_count)
    candidates.append(np.where(idx % 2 == 0, 1.0, -1.0).astype(np.complex128))
    candidates.append(phases[idx % 4])
    rng = np.random.default_rng(int(seed))
    for _ in range(int(n_candidates)):
        candidates.append(phases[rng.integers(0, len(phases), size=grid.pilot_count)])

    best = None
    best_score = (float("inf"), float("inf"), float("inf"))
    for cand in candidates:
        cond_mean, cond_max, bad_frac, _ = local_condition_for_pilot_phases(grid, cand, window_sc)
        score = (bad_frac, cond_max, cond_mean)
        if score < best_score:
            best_score = score
            best = cand.copy()
    assert best is not None
    z = phase_by_pilot_bins(grid, best)
    return PrecoderDesign(
        label=f"OPT-QPSK pilots W={int(window_sc)}",
        family="optimized_qpsk_pilot_phase",
        C=two_tx_precoder_from_relative_phase(z),
        metadata={
            "optimized_for_window_sc": int(window_sc),
            "n_candidates": int(n_candidates),
            "design_bad_fraction": float(best_score[0]),
            "design_cond_max": float(best_score[1]),
            "design_cond_mean": float(best_score[2]),
        },
    )


def build_designs(grid: ResourceGrid, cdd_delays: List[int], window_sc: int, opt_candidates: int, seed: int) -> List[PrecoderDesign]:
    out: List[PrecoderDesign] = [constant_design(grid)]
    out.extend(cdd_design(grid, d) for d in cdd_delays)
    out.append(pilot_cycle_design(grid, np.asarray([1.0, -1.0]), "Pilot ALT +/-"))
    out.append(pilot_cycle_design(grid, np.asarray([1.0, 1.0j, -1.0, -1.0j]), "Pilot QPSK cycle"))
    out.append(rb_qpsk_design(grid, rb_group=1, label="RB QPSK cycle"))
    out.append(rb_qpsk_design(grid, rb_group=4, label="PRG4RB QPSK cycle"))
    out.append(optimized_qpsk_design(grid, window_sc, opt_candidates, seed))
    return out


def make_local_estimator(
    grid: ResourceGrid,
    C: np.ndarray,
    noise_var_ls: float,
    loading: float,
    window_sc: int,
) -> LocalEstimator:
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    window = max(1, min(int(window_sc), int(grid.n_sc)))
    blocks: List[Dict[str, object]] = []
    conds: List[float] = []
    pilot_counts: List[int] = []
    for start in range(0, grid.n_sc, window):
        stop = min(grid.n_sc, start + window)
        target_local = np.arange(start, stop, dtype=np.int64)
        mask = (pilot_local >= start) & (pilot_local < stop)
        pilot_idx = np.where(mask)[0]
        if len(pilot_idx) == 0:
            center = 0.5 * (grid.subcarrier_indices[start] + grid.subcarrier_indices[stop - 1])
            nearest = int(np.argmin(np.abs(grid.pilot_subcarriers.astype(np.float64) - center)))
            pilot_idx = np.asarray([nearest], dtype=np.int64)
        A = np.asarray(C[pilot_local[pilot_idx], :], dtype=np.complex128)
        S = A @ A.conj().T
        S += (float(noise_var_ls) + float(loading)) * np.eye(S.shape[0], dtype=np.complex128)
        M = A.conj().T @ np.linalg.solve(S, np.eye(S.shape[0], dtype=np.complex128))
        cond = float(np.linalg.cond(A)) if len(pilot_idx) >= C.shape[1] else float("inf")
        conds.append(cond)
        pilot_counts.append(int(len(pilot_idx)))
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
    bad = [p for p in pilot_counts if p < C.shape[1]]
    return LocalEstimator(
        label=f"Alg2 local W={window}",
        blocks=blocks,
        cond_mean=float(np.mean(finite)) if finite else float("inf"),
        cond_max=float(max(finite)) if finite else float("inf"),
        metadata={
            "window_sc": int(window),
            "n_blocks": int(len(blocks)),
            "bad_window_fraction": float(len(bad) / max(len(pilot_counts), 1)),
            "min_pilots_per_window": int(min(pilot_counts)) if pilot_counts else 0,
            "mean_pilots_per_window": float(np.mean(pilot_counts)) if pilot_counts else 0.0,
        },
    )


def apply_local_estimator(est: LocalEstimator, ls_obs: np.ndarray, grid: ResourceGrid, C: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    n_rx = int(ls_obs.shape[0])
    n_tx = int(C.shape[1])
    h_hat = np.empty((n_rx, n_tx, grid.n_sc), dtype=np.complex128)
    g_hat = np.empty((n_rx, grid.n_sc), dtype=np.complex128)
    for block in est.blocks:
        pilot_idx = block["pilot_idx"]
        target_local = block["target_local"]
        h_block = ls_obs[:, pilot_idx] @ block["matrix"].T
        h_hat[:, :, target_local] = h_block[:, :, None]
        g_hat[:, target_local] = h_block @ C[target_local, :].T
    return g_hat, h_hat


def plot_best_gains(rows: List[Dict[str, object]], path: Path, metric: str = "gain_vs_best_cdd_db", title_suffix: str = "best CDD Alg2") -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    candidates = [r for r in rows if str(r["family"]) not in ("cdd_linear_phase", "constant_rank1")]
    candidates = [r for r in candidates if np.isfinite(float(r.get(metric, "nan")))]
    candidates.sort(key=lambda r: float(r[metric]), reverse=True)
    selected = candidates[:30]
    if not selected:
        return
    labels = [
        f"{r['precoder']}\nDS={r['delay_spread_ns']} sf={r['dmrs_spacing_sc']} W={r['window_sc']} SNR={r['snr_db']}"
        for r in selected
    ]
    gains = [float(r[metric]) for r in selected]
    colors = ["#2a9d8f" if g > 0 else "#b56576" for g in gains]
    fig, ax = plt.subplots(figsize=(12, 8.8))
    y = np.arange(len(selected))
    ax.barh(y, gains, color=colors)
    ax.axvline(0.0, color="black", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel(f"Alg2 NMSE gain vs {title_suffix} in same setting (dB)")
    ax.set_title(f"Best deterministic frequency-precoder designs for Alg2 vs {title_suffix}")
    ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_family_summary(rows: List[Dict[str, object]], path: Path, metric: str = "gain_vs_best_cdd_db", title_suffix: str = "best CDD Alg2") -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    groups: Dict[str, List[float]] = {}
    for row in rows:
        fam = str(row["precoder"])
        if str(row["family"]) == "constant_rank1":
            continue
        val = float(row.get(metric, "nan"))
        if np.isfinite(val):
            groups.setdefault(fam, []).append(val)
    labels = sorted(groups, key=lambda k: np.nanmedian(groups[k]), reverse=True)
    data = [groups[k] for k in labels]
    fig, ax = plt.subplots(figsize=(11, 5.8))
    ax.boxplot(data, tick_labels=labels, showfliers=False)
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_ylabel(f"Gain vs {title_suffix} in same setting (dB)")
    ax.set_title("Distribution over DS/DMRS/SNR/window settings")
    ax.tick_params(axis="x", rotation=25, labelsize=8)
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out) / f"precoder_design_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    delay_spreads = parse_float_list(args.delay_spreads_ns)
    dmrs_spacings = parse_int_list(args.dmrs_spacings)
    snrs = parse_float_list(args.snrs)
    windows = parse_int_list(args.windows_sc)
    cdd_delays = parse_int_list(args.cdd_delays)
    rows: List[Dict[str, object]] = []
    alg1_rows: List[Dict[str, object]] = []

    for delay_spread_ns in delay_spreads:
        for dmrs_spacing in dmrs_spacings:
            resource = ResourceConfig(
                carrier_bandwidth_mhz=100.0,
                scs_khz=30,
                n_fft=4096,
                n_prbs=48,
                pdsch_n_symbols=10,
                dmrs_symbol_indices=[2, 7],
                dmrs_spacing_sc=int(dmrs_spacing),
                prg_size_rb=4,
            )
            grid = build_resource_grid(resource)
            sample_period_ns = 1e9 / (float(resource.n_fft) * float(resource.scs_khz) * 1e3)
            pdp = make_exponential_pdp(delay_spread_ns, sample_period_ns, max_delay_factor=8.0)
            channel = ChannelConfig(delay_spread_ns=delay_spread_ns, max_delay_factor=8.0)
            rho_h = freq_corr_from_pdp(pdp, grid.n_fft, grid.n_sc - 1)
            bc_h_09 = first_below(rho_h, 0.9)
            bc_h_05 = first_below(rho_h, 0.5)
            pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)

            for window_sc in windows:
                design_seed = stable_seed(delay_spread_ns, dmrs_spacing, window_sc, "design", base=int(args.seed))
                designs = build_designs(grid, cdd_delays, int(window_sc), int(args.opt_candidates), design_seed)
                for snr_db in snrs:
                    noise_var = 10.0 ** (-float(snr_db) / 10.0)
                    noise_var_ls = noise_var / float(len(resource.dmrs_symbol_indices))
                    cdd_designs = [design for design in designs if design.family == "cdd_linear_phase"]
                    direct_estimators = {
                        design.label: make_direct_estimator(
                            grid=grid,
                            pdp=pdp,
                            delays=[0, int(design.metadata["cdd_delay_samples"])],
                            noise_var_ls=noise_var_ls,
                            loading=float(args.loading),
                        )
                        for design in cdd_designs
                    }
                    estimators = {
                        design.label: make_local_estimator(
                            grid, design.C, noise_var_ls, float(args.loading), int(window_sc)
                        )
                        for design in designs
                    }
                    accum: Dict[str, Dict[str, object]] = {
                        design.label: {
                            "eff": [],
                            "branch": [],
                            "design": design,
                            "est": estimators[design.label],
                        }
                        for design in designs
                    }
                    alg1_accum: Dict[str, List[float]] = {design.label: [] for design in cdd_designs}

                    rng = np.random.default_rng(stable_seed(delay_spread_ns, dmrs_spacing, window_sc, snr_db, base=int(args.seed)))
                    for _ in range(int(args.trials)):
                        realization = generate_tdl_channel(rng, grid, channel, n_tx=2, n_rx=4)
                        H = realization.H
                        shared_noise = (
                            rng.normal(size=(4, len(pilot_local)))
                            + 1j * rng.normal(size=(4, len(pilot_local)))
                        ) * math.sqrt(noise_var_ls / 2.0)
                        for design in designs:
                            g = equivalent_channel(H, design.C)
                            ls_obs = g[:, pilot_local] + shared_noise
                            est = estimators[design.label]
                            g_hat, h_hat = apply_local_estimator(est, ls_obs, grid, design.C)
                            accum[design.label]["eff"].append(nmse(g, g_hat))
                            accum[design.label]["branch"].append(nmse(H, h_hat))
                            if design.family == "cdd_linear_phase":
                                direct = direct_estimators[design.label]
                                g_hat_alg1 = ls_obs @ direct.matrix.T
                                alg1_accum[design.label].append(nmse(g, g_hat_alg1))

                    cdd_vals = {
                        label: mean_or_nan(values["eff"])
                        for label, values in accum.items()
                        if values["design"].family == "cdd_linear_phase"
                    }
                    best_cdd_label = min(cdd_vals, key=cdd_vals.get)
                    best_cdd_nmse = cdd_vals[best_cdd_label]
                    alg1_cdd_vals = {
                        label: mean_or_nan(vals)
                        for label, vals in alg1_accum.items()
                    }
                    best_alg1_cdd_label = min(alg1_cdd_vals, key=alg1_cdd_vals.get)
                    best_alg1_cdd_nmse = alg1_cdd_vals[best_alg1_cdd_label]
                    for label, vals in alg1_accum.items():
                        alg1_row = {
                            "estimator": "Alg1 Direct RMMSE CDD",
                            "precoder": label,
                            "family": "cdd_linear_phase",
                            "delay_spread_ns": float(delay_spread_ns),
                            "dmrs_spacing_sc": int(dmrs_spacing),
                            "window_sc": int(window_sc),
                            "snr_db": float(snr_db),
                            "trials": int(args.trials),
                            "n_sc": int(grid.n_sc),
                            "pilot_count_combined": int(grid.pilot_count),
                            "bc_h_0p9_sc": "" if bc_h_09 is None else int(bc_h_09),
                            "bc_h_0p5_sc": "" if bc_h_05 is None else int(bc_h_05),
                            "ce_nmse_eff_mean": mean_or_nan(vals),
                            "ce_nmse_eff_median": median_or_nan(vals),
                            "best_alg1_cdd_precoder": best_alg1_cdd_label,
                            "best_alg1_cdd_nmse": best_alg1_cdd_nmse,
                            "gain_vs_best_alg1_cdd_db": gain_db(best_alg1_cdd_nmse, mean_or_nan(vals)),
                        }
                        cdd_delay = label.split("d=")[-1] if "d=" in label else ""
                        alg1_row["cdd_delay_samples"] = cdd_delay
                        alg1_rows.append(alg1_row)
                    for label, values in accum.items():
                        design: PrecoderDesign = values["design"]
                        est: LocalEstimator = values["est"]
                        eff = mean_or_nan(values["eff"])
                        branch = mean_or_nan(values["branch"])
                        row = {
                            "estimator": "Alg2 LocalWindow",
                            "precoder": label,
                            "family": design.family,
                            "delay_spread_ns": float(delay_spread_ns),
                            "dmrs_spacing_sc": int(dmrs_spacing),
                            "window_sc": int(window_sc),
                            "snr_db": float(snr_db),
                            "trials": int(args.trials),
                            "n_sc": int(grid.n_sc),
                            "pilot_count_combined": int(grid.pilot_count),
                            "bc_h_0p9_sc": "" if bc_h_09 is None else int(bc_h_09),
                            "bc_h_0p5_sc": "" if bc_h_05 is None else int(bc_h_05),
                            "ce_nmse_eff_mean": eff,
                            "ce_nmse_eff_median": median_or_nan(values["eff"]),
                            "ce_nmse_branch_mean": branch,
                            "ce_nmse_branch_median": median_or_nan(values["branch"]),
                            "best_cdd_precoder": best_cdd_label,
                            "best_cdd_nmse": best_cdd_nmse,
                            "gain_vs_best_cdd_db": gain_db(best_cdd_nmse, eff),
                            "best_alg1_cdd_precoder": best_alg1_cdd_label,
                            "best_alg1_cdd_nmse": best_alg1_cdd_nmse,
                            "gain_vs_best_alg1_cdd_db": gain_db(best_alg1_cdd_nmse, eff),
                            "cond_mean": est.cond_mean,
                            "cond_max": est.cond_max,
                        }
                        row.update(est.metadata)
                        row.update({f"design_{k}": v for k, v in design.metadata.items() if not isinstance(v, list)})
                        rows.append(row)

    write_csv(rows, out_dir / "precoder_design_alg2_nmse.csv")
    write_csv(alg1_rows, out_dir / "precoder_design_alg1_cdd_baselines.csv")
    winners = [
        row for row in rows
        if row["family"] not in ("cdd_linear_phase", "constant_rank1")
        and np.isfinite(float(row["gain_vs_best_cdd_db"]))
        and float(row["gain_vs_best_cdd_db"]) > 0.0
    ]
    winners.sort(key=lambda r: float(r["gain_vs_best_cdd_db"]), reverse=True)
    write_csv(winners, out_dir / "precoder_design_alg2_wins.csv")
    alg1_winners = [
        row for row in rows
        if row["family"] not in ("cdd_linear_phase", "constant_rank1")
        and np.isfinite(float(row["gain_vs_best_alg1_cdd_db"]))
        and float(row["gain_vs_best_alg1_cdd_db"]) > 0.0
    ]
    alg1_winners.sort(key=lambda r: float(r["gain_vs_best_alg1_cdd_db"]), reverse=True)
    write_csv(alg1_winners, out_dir / "precoder_design_alg2_wins_vs_alg1_cdd.csv")
    with (out_dir / "run_config.json").open("w", encoding="utf-8") as f:
        json.dump({
            "delay_spreads_ns": delay_spreads,
            "dmrs_spacings": dmrs_spacings,
            "snrs": snrs,
            "windows_sc": windows,
            "cdd_delays": cdd_delays,
            "trials": int(args.trials),
            "opt_candidates": int(args.opt_candidates),
            "seed": int(args.seed),
            "note": "Algorithm-2-only NMSE search over deterministic frequency-domain precoder designs.",
        }, f, indent=2)
    plot_best_gains(rows, fig_dir / "experiment13_precoder_design_best_gains.png")
    plot_family_summary(rows, fig_dir / "experiment13_precoder_design_family_summary.png")
    plot_best_gains(
        rows,
        fig_dir / "experiment13_precoder_design_gain_vs_alg1_cdd.png",
        metric="gain_vs_best_alg1_cdd_db",
        title_suffix="best Alg1 CDD",
    )
    plot_family_summary(
        rows,
        fig_dir / "experiment13_precoder_design_family_vs_alg1_cdd.png",
        metric="gain_vs_best_alg1_cdd_db",
        title_suffix="best Alg1 CDD",
    )
    print(out_dir)
    return out_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/precoder_design_alg2")
    parser.add_argument("--fig-dir", default="docs/figures")
    parser.add_argument("--delay-spreads-ns", default="1,5,30,100")
    parser.add_argument("--dmrs-spacings", default="6,12,24")
    parser.add_argument("--snrs", default="-6,0,8")
    parser.add_argument("--windows-sc", default="48,72,96,144,192,288")
    parser.add_argument("--cdd-delays", default="64,128,256,512")
    parser.add_argument("--trials", type=int, default=80)
    parser.add_argument("--opt-candidates", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260614)
    parser.add_argument("--loading", type=float, default=1e-8)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
