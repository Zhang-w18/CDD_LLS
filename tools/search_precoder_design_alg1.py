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
from typing import Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ChannelConfig, ResourceConfig
from cdd_lls.phy.channel_tdl import generate_tdl_channel, make_exponential_pdp
from cdd_lls.phy.estimators import covariance_matrix, nmse
from cdd_lls.phy.precoding import equivalent_channel
from cdd_lls.phy.resource_grid import ResourceGrid, build_resource_grid, local_indices_for_subcarriers
from tools.run_nmse_algorithm_comparison import first_below, freq_corr_from_pdp, mean_or_nan, median_or_nan
from tools.search_precoder_design_alg2 import (
    PrecoderDesign,
    cdd_design,
    constant_design,
    pilot_cycle_design,
    rb_qpsk_design,
)


@dataclass(frozen=True)
class Alg1Estimator:
    label: str
    matrix: np.ndarray
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


def full_precoded_covariance(
    grid: ResourceGrid,
    pdp: np.ndarray,
    C: np.ndarray,
    target_local: np.ndarray,
    pilot_local: np.ndarray,
) -> np.ndarray:
    target_local = np.asarray(target_local, dtype=np.int64)
    pilot_local = np.asarray(pilot_local, dtype=np.int64)
    target_k = grid.subcarrier_indices[target_local]
    pilot_k = grid.subcarrier_indices[pilot_local]
    base = covariance_matrix(target_k, pilot_k, pdp, grid.n_fft)
    inner = np.asarray(C[target_local, :], dtype=np.complex128) @ np.asarray(C[pilot_local, :], dtype=np.complex128).conj().T
    return base * inner


def make_alg1_full_cov_estimator(
    grid: ResourceGrid,
    pdp: np.ndarray,
    C: np.ndarray,
    noise_var_ls: float,
    loading: float,
) -> Alg1Estimator:
    target_local = np.arange(grid.n_sc, dtype=np.int64)
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    R_pp = full_precoded_covariance(grid, pdp, C, pilot_local, pilot_local)
    R_tp = full_precoded_covariance(grid, pdp, C, target_local, pilot_local)
    A = R_pp + (float(noise_var_ls) + float(loading)) * np.eye(len(pilot_local), dtype=np.complex128)
    W = R_tp @ np.linalg.solve(A, np.eye(len(pilot_local), dtype=np.complex128))
    cond = float(np.linalg.cond(A))
    return Alg1Estimator(
        label="Alg1 full-cov RMMSE",
        matrix=W,
        metadata={"rpp_loaded_cond": cond},
    )


def build_designs(grid: ResourceGrid, cdd_delays: List[int], rb_groups: List[int]) -> List[PrecoderDesign]:
    designs: List[PrecoderDesign] = [constant_design(grid)]
    designs.extend(cdd_design(grid, d) for d in cdd_delays)
    designs.append(pilot_cycle_design(grid, np.asarray([1.0, -1.0]), "Pilot ALT +/-"))
    designs.append(pilot_cycle_design(grid, np.asarray([1.0, 1.0j, -1.0, -1.0j]), "Pilot QPSK cycle"))
    for group in rb_groups:
        designs.append(rb_qpsk_design(grid, rb_group=int(group), label=f"RB{int(group)} QPSK cycle"))
    return designs


def plot_best_gains(rows: List[Dict[str, object]], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    candidates = [
        r for r in rows
        if str(r["family"]) not in ("cdd_linear_phase", "constant_rank1")
        and np.isfinite(float(r.get("gain_vs_best_cdd_db", "nan")))
    ]
    candidates.sort(key=lambda r: float(r["gain_vs_best_cdd_db"]), reverse=True)
    selected = candidates[:30]
    if not selected:
        return
    labels = [
        f"{r['precoder']}\nDS={r['delay_spread_ns']} sf={r['dmrs_spacing_sc']} SNR={r['snr_db']}"
        for r in selected
    ]
    gains = [float(r["gain_vs_best_cdd_db"]) for r in selected]
    colors = ["#2a9d8f" if g > 0 else "#b56576" for g in gains]
    fig, ax = plt.subplots(figsize=(12, 8.8))
    y = np.arange(len(selected))
    ax.barh(y, gains, color=colors)
    ax.axvline(0.0, color="black", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Alg1 full-cov RMMSE NMSE gain vs best CDD in same setting (dB)")
    ax.set_title("Best deterministic frequency precoders for Alg1")
    ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_family_summary(rows: List[Dict[str, object]], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    groups: Dict[str, List[float]] = {}
    for row in rows:
        if str(row["family"]) == "constant_rank1":
            continue
        val = float(row.get("gain_vs_best_cdd_db", "nan"))
        if np.isfinite(val):
            groups.setdefault(str(row["precoder"]), []).append(val)
    labels = sorted(groups, key=lambda k: np.nanmedian(groups[k]), reverse=True)
    data = [groups[k] for k in labels]
    fig, ax = plt.subplots(figsize=(11.5, 5.8))
    ax.boxplot(data, tick_labels=labels, showfliers=False)
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_ylabel("Gain vs best CDD in same setting (dB)")
    ax.set_title("Alg1 full-cov RMMSE: distribution over DS/DMRS/SNR")
    ax.tick_params(axis="x", rotation=25, labelsize=8)
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out) / f"alg1_precoder_design_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    delay_spreads = parse_float_list(args.delay_spreads_ns)
    dmrs_spacings = parse_int_list(args.dmrs_spacings)
    snrs = parse_float_list(args.snrs)
    cdd_delays = parse_int_list(args.cdd_delays)
    rb_groups = parse_int_list(args.rb_groups)
    rows: List[Dict[str, object]] = []

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
            designs = build_designs(grid, cdd_delays, rb_groups)

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
                accum: Dict[str, Dict[str, object]] = {
                    design.label: {"eff": [], "design": design, "est": estimators[design.label]}
                    for design in designs
                }
                rng = np.random.default_rng(stable_seed(delay_spread_ns, dmrs_spacing, snr_db, base=int(args.seed)))
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
                        g_hat = ls_obs @ estimators[design.label].matrix.T
                        accum[design.label]["eff"].append(nmse(g, g_hat))

                cdd_vals = {
                    label: mean_or_nan(values["eff"])
                    for label, values in accum.items()
                    if values["design"].family == "cdd_linear_phase"
                }
                best_cdd_label = min(cdd_vals, key=cdd_vals.get)
                best_cdd_nmse = cdd_vals[best_cdd_label]
                for label, values in accum.items():
                    design: PrecoderDesign = values["design"]
                    est: Alg1Estimator = values["est"]
                    eff = mean_or_nan(values["eff"])
                    row = {
                        "estimator": "Alg1 full-cov RMMSE",
                        "precoder": label,
                        "family": design.family,
                        "delay_spread_ns": float(delay_spread_ns),
                        "dmrs_spacing_sc": int(dmrs_spacing),
                        "snr_db": float(snr_db),
                        "trials": int(args.trials),
                        "n_sc": int(grid.n_sc),
                        "pilot_count_combined": int(grid.pilot_count),
                        "bc_h_0p9_sc": "" if bc_h_09 is None else int(bc_h_09),
                        "bc_h_0p5_sc": "" if bc_h_05 is None else int(bc_h_05),
                        "ce_nmse_eff_mean": eff,
                        "ce_nmse_eff_median": median_or_nan(values["eff"]),
                        "best_cdd_precoder": best_cdd_label,
                        "best_cdd_nmse": best_cdd_nmse,
                        "gain_vs_best_cdd_db": gain_db(best_cdd_nmse, eff),
                        "rpp_loaded_cond": est.metadata["rpp_loaded_cond"],
                    }
                    row.update({f"design_{k}": v for k, v in design.metadata.items() if not isinstance(v, list)})
                    rows.append(row)

    write_csv(rows, out_dir / "alg1_precoder_design_nmse.csv")
    winners = [
        row for row in rows
        if row["family"] not in ("cdd_linear_phase", "constant_rank1")
        and np.isfinite(float(row["gain_vs_best_cdd_db"]))
        and float(row["gain_vs_best_cdd_db"]) > 0.0
    ]
    winners.sort(key=lambda r: float(r["gain_vs_best_cdd_db"]), reverse=True)
    write_csv(winners, out_dir / "alg1_precoder_design_wins.csv")
    with (out_dir / "run_config.json").open("w", encoding="utf-8") as f:
        json.dump({
            "delay_spreads_ns": delay_spreads,
            "dmrs_spacings": dmrs_spacings,
            "snrs": snrs,
            "cdd_delays": cdd_delays,
            "rb_groups": rb_groups,
            "trials": int(args.trials),
            "seed": int(args.seed),
            "note": "Algorithm-1 full-covariance RMMSE comparison over deterministic frequency-domain precoder designs.",
        }, f, indent=2)
    plot_best_gains(rows, fig_dir / "experiment14_alg1_precoder_design_best_gains.png")
    plot_family_summary(rows, fig_dir / "experiment14_alg1_precoder_design_family_summary.png")
    print(out_dir)
    return out_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/precoder_design_alg1")
    parser.add_argument("--fig-dir", default="docs/figures")
    parser.add_argument("--delay-spreads-ns", default="1,5,30,100")
    parser.add_argument("--dmrs-spacings", default="6,12,24")
    parser.add_argument("--snrs", default="-6,0,8")
    parser.add_argument("--cdd-delays", default="64,128,256,512")
    parser.add_argument("--rb-groups", default="1,2,4,8,12")
    parser.add_argument("--trials", type=int, default=80)
    parser.add_argument("--seed", type=int, default=20260614)
    parser.add_argument("--loading", type=float, default=1e-8)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
