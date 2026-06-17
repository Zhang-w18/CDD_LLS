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
from typing import Dict, Iterable, List

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ChannelConfig, ResourceConfig
from cdd_lls.core.mcs import get_mcs
from cdd_lls.phy.channel_tdl import generate_tdl_channel, make_exponential_pdp
from cdd_lls.phy.estimators import nmse
from cdd_lls.phy.precoding import equivalent_channel
from cdd_lls.phy.resource_grid import ResourceGrid, build_resource_grid, local_indices_for_subcarriers
from tools.search_precoder_design_alg1 import make_alg1_full_cov_estimator
from tools.search_precoder_design_alg2 import (
    PrecoderDesign,
    cdd_design,
    constant_design,
    two_tx_precoder_from_relative_phase,
)


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


def mean_or_nan(values: Iterable[float]) -> float:
    vals = [float(x) for x in values if np.isfinite(float(x))]
    return float(np.mean(vals)) if vals else float("nan")


def pct_or_nan(values: Iterable[float], pct: float) -> float:
    vals = [float(x) for x in values if np.isfinite(float(x))]
    return float(np.percentile(vals, float(pct))) if vals else float("nan")


def gain_db(reference: float, candidate: float) -> float:
    if not np.isfinite(reference) or not np.isfinite(candidate) or reference <= 0.0 or candidate <= 0.0:
        return float("nan")
    return float(10.0 * np.log10(reference / candidate))


def two_tx_precoder_from_power_phase(power0: np.ndarray, z: np.ndarray) -> np.ndarray:
    p0 = np.asarray(power0, dtype=np.float64).reshape(-1)
    zz = np.asarray(z, dtype=np.complex128).reshape(-1)
    if p0.shape != zz.shape:
        raise ValueError("power0 and z must have the same length.")
    p0 = np.clip(p0, 0.0, 1.0)
    C = np.empty((len(zz), 2), dtype=np.complex128)
    C[:, 0] = np.sqrt(p0)
    C[:, 1] = np.sqrt(1.0 - p0) * zz
    return C


def rb_mpsk_design(grid: ResourceGrid, m: int, rb_group: int, step: int = 1) -> PrecoderDesign:
    m = int(m)
    rb_group = max(int(rb_group), 1)
    step = int(step)
    phases = np.exp(1j * 2.0 * np.pi * step * np.arange(m, dtype=np.float64) / float(m))
    z = np.empty(grid.n_sc, dtype=np.complex128)
    for local in range(grid.n_sc):
        rb = int(local // 12)
        idx = int((rb // rb_group) % m)
        z[local] = phases[idx]
    return PrecoderDesign(
        label=f"RB{rb_group} {m}PSK step{step}",
        family="rb_mpsk_cycling",
        C=two_tx_precoder_from_relative_phase(z),
        metadata={"rb_group": rb_group, "mpsk_order": m, "phase_step": step},
    )


def sinusoid_phase_design(grid: ResourceGrid, cycles: float, beta_pi: float) -> PrecoderDesign:
    x = (np.arange(grid.n_sc, dtype=np.float64) + 0.5) / float(grid.n_sc)
    beta = float(beta_pi) * np.pi
    z = np.exp(1j * beta * np.sin(2.0 * np.pi * float(cycles) * x))
    return PrecoderDesign(
        label=f"sin phase cyc={float(cycles):g} beta={float(beta_pi):g}pi",
        family="smooth_phase",
        C=two_tx_precoder_from_relative_phase(z),
        metadata={"cycles": float(cycles), "beta_pi": float(beta_pi)},
    )


def chirp_phase_design(grid: ResourceGrid, chirp_rate: float) -> PrecoderDesign:
    x = (np.arange(grid.n_sc, dtype=np.float64) - 0.5 * (grid.n_sc - 1)) / float(grid.n_sc)
    z = np.exp(1j * 2.0 * np.pi * float(chirp_rate) * x * x)
    return PrecoderDesign(
        label=f"quadratic phase q={float(chirp_rate):g}",
        family="smooth_chirp",
        C=two_tx_precoder_from_relative_phase(z),
        metadata={"chirp_rate": float(chirp_rate)},
    )


def hybrid_cdd_rb_design(grid: ResourceGrid, delay: int, rb_group: int, m: int = 4) -> PrecoderDesign:
    cdd_z = np.exp(-1j * 2.0 * np.pi * grid.subcarrier_indices.astype(np.float64) * float(delay) / float(grid.n_fft))
    rb = rb_mpsk_design(grid, m=m, rb_group=rb_group, step=1)
    z = cdd_z * (rb.C[:, 1] / rb.C[:, 0])
    return PrecoderDesign(
        label=f"hybrid d={int(delay)} + RB{int(rb_group)} {int(m)}PSK",
        family="hybrid_cdd_mpsk",
        C=two_tx_precoder_from_relative_phase(z),
        metadata={"cdd_delay_samples": int(delay), "rb_group": int(rb_group), "mpsk_order": int(m)},
    )


def power_cycle_design(grid: ResourceGrid, rb_group: int, powers: List[float], phase_m: int = 4) -> PrecoderDesign:
    rb_group = max(int(rb_group), 1)
    powers_arr = np.asarray(powers, dtype=np.float64)
    phases = np.exp(1j * 2.0 * np.pi * np.arange(int(phase_m), dtype=np.float64) / float(phase_m))
    p0 = np.empty(grid.n_sc, dtype=np.float64)
    z = np.empty(grid.n_sc, dtype=np.complex128)
    for local in range(grid.n_sc):
        rb = int(local // 12)
        group = int(rb // rb_group)
        p0[local] = powers_arr[group % len(powers_arr)]
        z[local] = phases[group % len(phases)]
    label_p = "-".join(f"{p:.2g}" for p in powers)
    return PrecoderDesign(
        label=f"RB{rb_group} power[{label_p}] {phase_m}PSK",
        family="power_phase_cycling",
        C=two_tx_precoder_from_power_phase(p0, z),
        metadata={"rb_group": int(rb_group), "powers": label_p, "phase_m": int(phase_m)},
    )


def random_anchor_phase_design(
    grid: ResourceGrid,
    anchor_rb: int,
    phase_levels: int,
    seed: int,
    label: str,
) -> PrecoderDesign:
    anchor_rb = max(int(anchor_rb), 1)
    rng = np.random.default_rng(int(seed))
    n_rb = int(math.ceil(grid.n_sc / 12.0))
    n_anchor = int(math.ceil(n_rb / float(anchor_rb))) + 1
    phases = 2.0 * np.pi * rng.integers(0, int(phase_levels), size=n_anchor) / float(phase_levels)
    rb_points = np.arange(n_anchor, dtype=np.float64) * float(anchor_rb)
    rb_index = (np.arange(grid.n_sc, dtype=np.float64) + 0.5) / 12.0
    unwrapped = np.unwrap(phases)
    interp = np.interp(rb_index, rb_points, unwrapped)
    z = np.exp(1j * interp)
    return PrecoderDesign(
        label=label,
        family="smooth_random_phase",
        C=two_tx_precoder_from_relative_phase(z),
        metadata={"anchor_rb": int(anchor_rb), "phase_levels": int(phase_levels), "design_seed": int(seed)},
    )


def build_designs(grid: ResourceGrid, cdd_delays: List[int], rb_groups: List[int], seed: int) -> List[PrecoderDesign]:
    designs: List[PrecoderDesign] = [constant_design(grid)]
    designs.extend(cdd_design(grid, d) for d in cdd_delays)

    for rb_group in rb_groups:
        for m in (4, 8, 16):
            designs.append(rb_mpsk_design(grid, m=m, rb_group=rb_group, step=1))
        designs.append(power_cycle_design(grid, rb_group=rb_group, powers=[0.1, 0.9], phase_m=4))
        designs.append(power_cycle_design(grid, rb_group=rb_group, powers=[0.25, 0.75], phase_m=4))
        designs.append(power_cycle_design(grid, rb_group=rb_group, powers=[0.0, 1.0], phase_m=2))

    for cycles in (1, 2, 4, 8, 12):
        for beta_pi in (0.5, 1.0, 1.5, 2.0):
            designs.append(sinusoid_phase_design(grid, cycles=cycles, beta_pi=beta_pi))

    for q in (1, 2, 4, 8, 16):
        designs.append(chirp_phase_design(grid, chirp_rate=q))

    for delay in (16, 32, 48):
        for rb_group in (4, 8, 12):
            designs.append(hybrid_cdd_rb_design(grid, delay=delay, rb_group=rb_group, m=4))
            designs.append(hybrid_cdd_rb_design(grid, delay=delay, rb_group=rb_group, m=8))

    for anchor_rb in (2, 4, 8, 12):
        for i in range(3):
            designs.append(random_anchor_phase_design(
                grid,
                anchor_rb=anchor_rb,
                phase_levels=16,
                seed=stable_seed("rand", anchor_rb, i, base=seed),
                label=f"smooth rand anchorRB{anchor_rb} #{i}",
            ))

    unique: Dict[str, PrecoderDesign] = {}
    for design in designs:
        unique.setdefault(design.label, design)
    return list(unique.values())


def corr_effective_rank_and_bw(gamma_vectors: List[np.ndarray], max_delta: int = 144) -> tuple[float, float, float]:
    if len(gamma_vectors) < 3:
        return float("nan"), float("nan"), float("nan")
    X = np.vstack([np.asarray(x, dtype=np.float64).reshape(1, -1) for x in gamma_vectors])
    X = X - np.mean(X, axis=0, keepdims=True)
    std = np.std(X, axis=0, keepdims=True)
    X = X / np.maximum(std, 1e-12)
    C = (X.T @ X) / float(max(X.shape[0] - 1, 1))
    tr = float(np.trace(C).real)
    tr2 = float(np.sum(np.abs(C) ** 2).real)
    eff_rank = tr * tr / tr2 if tr2 > 0.0 else float("nan")
    corrs = []
    for delta in range(1, min(int(max_delta), X.shape[1] - 1) + 1):
        vals = X[:, :-delta] * X[:, delta:]
        corrs.append(float(np.mean(vals)))
    abs_corr = np.abs(np.asarray(corrs, dtype=np.float64))
    below_05 = np.where(abs_corr < 0.5)[0]
    below_02 = np.where(abs_corr < 0.2)[0]
    bw05 = float(below_05[0] + 1) if below_05.size else float("nan")
    bw02 = float(below_02[0] + 1) if below_02.size else float("nan")
    return float(eff_rank), bw05, bw02


def plot_tradeoff(rows: List[Dict[str, object]], path: Path) -> None:
    if not rows:
        return
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    candidates = [
        r for r in rows
        if str(r["family"]) not in ("cdd_linear_phase", "constant_rank1")
        and np.isfinite(float(r.get("ce_gain_vs_best_cdd_db", "nan")))
        and np.isfinite(float(r.get("mi_p05_delta_vs_best_cdd", "nan")))
    ]
    if not candidates:
        return
    families = sorted(set(str(r["family"]) for r in candidates))
    colors = {fam: plt.cm.tab20(i % 20) for i, fam in enumerate(families)}
    fig, ax = plt.subplots(figsize=(10, 6.2))
    for fam in families:
        pts = [r for r in candidates if str(r["family"]) == fam]
        x = [float(r["ce_gain_vs_best_cdd_db"]) for r in pts]
        y = [float(r["mi_p05_delta_vs_best_cdd"]) for r in pts]
        ax.scatter(x, y, s=22, alpha=0.62, label=fam, color=colors[fam])
    ax.axvline(0.0, color="black", linewidth=0.9)
    ax.axhline(0.0, color="black", linewidth=0.9)
    ax.set_xlabel("CE NMSE gain vs best CDD (dB)")
    ax.set_ylabel("5%-tile capped-MI margin vs best CDD (bit/RE)")
    ax.set_title("Precoder CE-vs-diversity tradeoff candidates")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_family_summary(rows: List[Dict[str, object]], path: Path) -> None:
    if not rows:
        return
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    groups: Dict[str, List[float]] = {}
    for row in rows:
        if str(row["family"]) in ("cdd_linear_phase", "constant_rank1"):
            continue
        if int(row["mcs_index"]) != 8 or float(row["snr_db"]) not in (4.0, 6.0):
            continue
        val = float(row.get("composite_score", "nan"))
        if np.isfinite(val):
            groups.setdefault(str(row["family"]), []).append(val)
    if not groups:
        return
    labels = sorted(groups, key=lambda k: np.median(groups[k]), reverse=True)
    data = [groups[k] for k in labels]
    fig, ax = plt.subplots(figsize=(11, 5.8))
    ax.boxplot(data, tick_labels=labels, showfliers=False)
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_ylabel("Composite score")
    ax.set_title("Experiment 16 family summary, MCS8 SNR 4/6 dB")
    ax.tick_params(axis="x", rotation=25, labelsize=8)
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out) / f"precoder_diversity_tradeoff_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    delay_spreads = parse_float_list(args.delay_spreads_ns)
    snrs = parse_float_list(args.snrs)
    mcs_indexes = parse_int_list(args.mcs_indexes)
    cdd_delays = parse_int_list(args.cdd_delays)
    baseline_cdd_delays = set(parse_int_list(args.baseline_cdd_delays))
    rb_groups = parse_int_list(args.rb_groups)
    rows: List[Dict[str, object]] = []

    resource = ResourceConfig(
        carrier_bandwidth_mhz=100.0,
        scs_khz=30,
        n_fft=4096,
        n_prbs=int(args.n_prbs),
        pdsch_n_symbols=10,
        dmrs_symbol_indices=[2, 7],
        dmrs_spacing_sc=int(args.dmrs_spacing_sc),
        prg_size_rb=4,
    )
    grid = build_resource_grid(resource)
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    data_unique_local = np.unique(local_indices_for_subcarriers(grid, grid.data_subcarrier_indices))
    designs = build_designs(grid, cdd_delays, rb_groups, int(args.seed))

    for delay_spread_ns in delay_spreads:
        sample_period_ns = 1e9 / (float(resource.n_fft) * float(resource.scs_khz) * 1e3)
        pdp = make_exponential_pdp(delay_spread_ns, sample_period_ns, max_delay_factor=8.0)
        channel = ChannelConfig(delay_spread_ns=float(delay_spread_ns), max_delay_factor=8.0)

        for snr_db in snrs:
            noise_var = 10.0 ** (-float(snr_db) / 10.0)
            noise_var_ls = noise_var / float(len(resource.dmrs_symbol_indices))
            estimators = {
                design.label: make_alg1_full_cov_estimator(
                    grid=grid,
                    pdp=pdp,
                    C=design.C,
                    noise_var_ls=noise_var_ls,
                    loading=float(args.loading),
                )
                for design in designs
            }
            accum = {
                design.label: {
                    "design": design,
                    "ce_nmse": [],
                    "slot_capacity": [],
                    "gamma_p05_db": [],
                    "gamma_p10_db": [],
                    "gamma_vectors": [],
                    "rpp_loaded_cond": estimators[design.label].metadata["rpp_loaded_cond"],
                }
                for design in designs
            }
            rng = np.random.default_rng(stable_seed(delay_spread_ns, snr_db, base=int(args.seed)))
            for trial in range(1, int(args.trials) + 1):
                realization = generate_tdl_channel(rng, grid, channel, n_tx=2, n_rx=4)
                H = realization.H
                pilot_noise = (
                    rng.normal(size=(4, len(pilot_local)))
                    + 1j * rng.normal(size=(4, len(pilot_local)))
                ) * math.sqrt(noise_var_ls / 2.0)
                for design in designs:
                    est = estimators[design.label]
                    g = equivalent_channel(H, design.C)
                    ls_obs = g[:, pilot_local] + pilot_noise
                    g_hat = ls_obs @ est.matrix.T
                    ce = nmse(g, g_hat)
                    gamma = np.sum(np.abs(g[:, data_unique_local]) ** 2, axis=0) / float(noise_var)
                    gamma_db = 10.0 * np.log10(np.maximum(gamma, 1e-12))
                    capacity = np.log2(1.0 + gamma)
                    acc = accum[design.label]
                    acc["ce_nmse"].append(float(ce))
                    acc["slot_capacity"].append(float(np.mean(capacity)))
                    acc["gamma_p05_db"].append(float(np.percentile(gamma_db, 5.0)))
                    acc["gamma_p10_db"].append(float(np.percentile(gamma_db, 10.0)))
                    if trial <= int(args.diversity_trials):
                        acc["gamma_vectors"].append(gamma.copy())
                if int(args.progress_every) > 0 and trial % int(args.progress_every) == 0:
                    print(
                        f"[progress] DS={float(delay_spread_ns):g} ns, "
                        f"SNR={float(snr_db):g} dB, trial={trial}/{int(args.trials)}",
                        flush=True,
                    )

            base_rows: List[Dict[str, object]] = []
            for design in designs:
                acc = accum[design.label]
                eff_rank, corr_bw_05, corr_bw_02 = corr_effective_rank_and_bw(
                    acc["gamma_vectors"],
                    max_delta=int(args.max_corr_delta),
                )
                cdd_delay = design.metadata.get("cdd_delay_samples", None)
                is_baseline_cdd = (
                    design.family == "cdd_linear_phase"
                    and cdd_delay is not None
                    and int(cdd_delay) in baseline_cdd_delays
                )
                row_base = {
                    "precoder": design.label,
                    "family": design.family,
                    "delay_spread_ns": float(delay_spread_ns),
                    "snr_db": float(snr_db),
                    "trials": int(args.trials),
                    "n_sc": int(grid.n_sc),
                    "n_prbs": int(args.n_prbs),
                    "dmrs_spacing_sc": int(args.dmrs_spacing_sc),
                    "dmrs_overhead": float(grid.dmrs_overhead),
                    "ce_nmse_mean": mean_or_nan(acc["ce_nmse"]),
                    "ce_nmse_median": pct_or_nan(acc["ce_nmse"], 50.0),
                    "slot_capacity_mean": mean_or_nan(acc["slot_capacity"]),
                    "slot_capacity_p05": pct_or_nan(acc["slot_capacity"], 5.0),
                    "gamma_p05_db_mean": mean_or_nan(acc["gamma_p05_db"]),
                    "gamma_p10_db_mean": mean_or_nan(acc["gamma_p10_db"]),
                    "gamma_corr_eff_rank": eff_rank,
                    "gamma_corr_bw_0p5_sc": corr_bw_05,
                    "gamma_corr_bw_0p2_sc": corr_bw_02,
                    "rpp_loaded_cond": float(acc["rpp_loaded_cond"]),
                    "is_baseline_cdd": int(is_baseline_cdd),
                }
                row_base.update({f"design_{k}": v for k, v in design.metadata.items() if not isinstance(v, list)})
                base_rows.append(row_base)

            for mcs_index in mcs_indexes:
                mcs = get_mcs("nr_256qam", int(mcs_index), None, None)
                spectral_eff = float(mcs.qm) * float(mcs.code_rate)
                tmp = []
                for row_base in base_rows:
                    design = accum[row_base["precoder"]]["design"]
                    acc = accum[row_base["precoder"]]
                    capped = [
                        min(float(x), float(mcs.qm))
                        for x in acc["slot_capacity"]
                    ]
                    row = dict(row_base)
                    row.update({
                        "mcs_table": "nr_256qam",
                        "mcs_index": int(mcs_index),
                        "qm": int(mcs.qm),
                        "code_rate": float(mcs.code_rate),
                        "required_bpre": spectral_eff,
                        "slot_mi_capped_mean": mean_or_nan(capped),
                        "slot_mi_capped_p05": pct_or_nan(capped, 5.0),
                        "slot_mi_margin_p05": pct_or_nan([x - spectral_eff for x in capped], 5.0),
                    })
                    tmp.append(row)

                baseline = [r for r in tmp if int(r["is_baseline_cdd"]) == 1]
                best_cdd_ce = min(baseline, key=lambda r: float(r["ce_nmse_mean"]))
                best_cdd_mi = max(baseline, key=lambda r: float(r["slot_mi_margin_p05"]))
                best_cdd_rank = max(baseline, key=lambda r: float(r["gamma_corr_eff_rank"]))
                for row in tmp:
                    row["best_cdd_ce_precoder"] = best_cdd_ce["precoder"]
                    row["best_cdd_ce_nmse"] = best_cdd_ce["ce_nmse_mean"]
                    row["ce_gain_vs_best_cdd_db"] = gain_db(float(best_cdd_ce["ce_nmse_mean"]), float(row["ce_nmse_mean"]))
                    row["best_cdd_mi_precoder"] = best_cdd_mi["precoder"]
                    row["best_cdd_mi_p05"] = best_cdd_mi["slot_mi_margin_p05"]
                    row["mi_p05_delta_vs_best_cdd"] = float(row["slot_mi_margin_p05"]) - float(best_cdd_mi["slot_mi_margin_p05"])
                    row["best_cdd_rank_precoder"] = best_cdd_rank["precoder"]
                    row["rank_ratio_vs_best_cdd"] = (
                        float(row["gamma_corr_eff_rank"]) / float(best_cdd_rank["gamma_corr_eff_rank"])
                        if float(best_cdd_rank["gamma_corr_eff_rank"]) > 0.0 else float("nan")
                    )
                    row["corr_bw_delta_vs_best_cdd"] = (
                        float(row["gamma_corr_bw_0p5_sc"]) - float(best_cdd_mi["gamma_corr_bw_0p5_sc"])
                        if np.isfinite(float(row["gamma_corr_bw_0p5_sc"])) and np.isfinite(float(best_cdd_mi["gamma_corr_bw_0p5_sc"]))
                        else float("nan")
                    )
                    row["qualified_ce_and_diversity"] = int(
                        str(row["family"]) not in ("cdd_linear_phase", "constant_rank1")
                        and float(row["ce_gain_vs_best_cdd_db"]) > float(args.min_ce_gain_db)
                        and float(row["mi_p05_delta_vs_best_cdd"]) >= -float(args.max_mi_loss)
                        and float(row["rank_ratio_vs_best_cdd"]) >= float(args.min_rank_ratio)
                    )
                    row["composite_score"] = (
                        float(row["ce_gain_vs_best_cdd_db"])
                        + 3.0 * float(row["mi_p05_delta_vs_best_cdd"])
                        + 0.5 * (float(row["rank_ratio_vs_best_cdd"]) - 1.0)
                    )
                    rows.append(row)

            write_csv(rows, out_dir / "precoder_diversity_tradeoff_partial.csv")
            print(
                f"[done] DS={float(delay_spread_ns):g} ns, SNR={float(snr_db):g} dB, "
                f"designs={len(designs)}",
                flush=True,
            )

    write_csv(rows, out_dir / "precoder_diversity_tradeoff.csv")
    qualified = [
        row for row in rows
        if int(row.get("qualified_ce_and_diversity", 0)) == 1
    ]
    qualified.sort(key=lambda r: float(r["composite_score"]), reverse=True)
    write_csv(qualified, out_dir / "qualified_precoders.csv")
    with (out_dir / "run_config.json").open("w", encoding="utf-8") as f:
        json.dump(vars(args), f, indent=2, ensure_ascii=False)
    plot_tradeoff(rows, fig_dir / "experiment16_precoder_ce_diversity_tradeoff.png")
    plot_family_summary(rows, fig_dir / "experiment16_precoder_family_summary.png")
    print(out_dir)
    return out_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/precoder_diversity_tradeoff")
    parser.add_argument("--fig-dir", default="docs/figures")
    parser.add_argument("--delay-spreads-ns", default="30,100")
    parser.add_argument("--snrs", default="4,6,12,18")
    parser.add_argument("--mcs-indexes", default="8,14,20")
    parser.add_argument("--trials", type=int, default=120)
    parser.add_argument("--diversity-trials", type=int, default=120)
    parser.add_argument("--n-prbs", type=int, default=48)
    parser.add_argument("--dmrs-spacing-sc", type=int, default=6)
    parser.add_argument("--cdd-delays", default="32,64,128,256")
    parser.add_argument("--baseline-cdd-delays", default="64,128")
    parser.add_argument("--rb-groups", default="1,2,4,6,8,12,16,24")
    parser.add_argument("--loading", type=float, default=1e-8)
    parser.add_argument("--seed", type=int, default=20260614)
    parser.add_argument("--progress-every", type=int, default=30)
    parser.add_argument("--max-corr-delta", type=int, default=144)
    parser.add_argument("--min-ce-gain-db", type=float, default=0.3)
    parser.add_argument("--max-mi-loss", type=float, default=0.02)
    parser.add_argument("--min-rank-ratio", type=float, default=0.9)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
