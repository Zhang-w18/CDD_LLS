from __future__ import annotations

import argparse
import csv
import json
import math
import os
from datetime import datetime
from pathlib import Path
import sys
from typing import Dict, Iterable, List

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ChannelConfig, ResourceConfig
from cdd_lls.phy.channel_tdl import generate_tdl_channel, make_exponential_pdp
from cdd_lls.phy.estimators import nmse, shifted_pdp
from cdd_lls.phy.precoding import cdd_equivalent_from_branches
from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
from tools.run_nmse_algorithm_comparison import (
    apply_basis,
    apply_block,
    apply_port_dmrs,
    first_below,
    freq_corr_from_pdp,
    make_basis_estimator,
    make_block_estimator,
    make_direct_estimator,
    make_port_estimators,
    median_or_nan,
    mean_or_nan,
    write_csv,
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


def add_result(accum: Dict[str, Dict[str, object]], label: str, family: str, eff: float, branch: float, meta: Dict[str, object]) -> None:
    if label not in accum:
        accum[label] = {"family": family, "eff": [], "branch": [], "meta": meta}
    accum[label]["eff"].append(float(eff))
    accum[label]["branch"].append(float(branch))


def gain_db(reference: float, candidate: float) -> float:
    if not np.isfinite(reference) or not np.isfinite(candidate) or reference <= 0.0 or candidate <= 0.0:
        return float("nan")
    return float(10.0 * np.log10(reference / candidate))


def plot_best_gains(rows: List[Dict[str, object]], path: Path, min_gain_db: float = -3.0) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    candidates = [r for r in rows if str(r["algorithm"]) != "Alg1 Direct RMMSE WB 576sc"]
    candidates = [r for r in candidates if np.isfinite(float(r.get("gain_vs_alg1_db", "nan")))]
    candidates.sort(key=lambda r: float(r["gain_vs_alg1_db"]), reverse=True)
    selected = candidates[:30]
    if not selected:
        return
    labels = [
        f"{r['algorithm']}\nDS={r['delay_spread_ns']}ns d={r['cdd_delay_samples']} sf={r['dmrs_spacing_sc']} SNR={r['snr_db']}"
        for r in selected
    ]
    gains = [float(r["gain_vs_alg1_db"]) for r in selected]
    colors = ["#2a9d8f" if g > 0 else "#b56576" for g in gains]
    fig, ax = plt.subplots(figsize=(12, 8.5))
    y = np.arange(len(selected))
    ax.barh(y, gains, color=colors)
    ax.axvline(0.0, color="black", linewidth=0.8)
    ax.axvline(min_gain_db, color="gray", linewidth=0.8, linestyle=":")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("NMSE gain vs Alg1 Direct RMMSE (dB)")
    ax.set_title("Best same-overhead NMSE gains over Alg1")
    ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_alg3_gain_maps(rows: List[Dict[str, object]], snr_db: float, path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    target = [
        r for r in rows
        if r["algorithm"] == "Alg3 Port-DMRS RMMSE equal-total"
        and abs(float(r["snr_db"]) - float(snr_db)) < 1e-9
    ]
    if not target:
        return
    spacings = sorted({int(r["dmrs_spacing_sc"]) for r in target})
    delay_spreads = sorted({float(r["delay_spread_ns"]) for r in target})
    cdd_delays = sorted({int(r["cdd_delay_samples"]) for r in target})
    fig, axes = plt.subplots(1, len(spacings), figsize=(4.8 * len(spacings), 4.4), sharey=True)
    if len(spacings) == 1:
        axes = [axes]
    vmax = max(abs(float(r["gain_vs_alg1_db"])) for r in target if np.isfinite(float(r["gain_vs_alg1_db"])))
    vmax = max(0.5, min(8.0, vmax))
    for ax, sf in zip(axes, spacings):
        mat = np.full((len(delay_spreads), len(cdd_delays)), np.nan)
        for r in target:
            if int(r["dmrs_spacing_sc"]) != sf:
                continue
            i = delay_spreads.index(float(r["delay_spread_ns"]))
            j = cdd_delays.index(int(r["cdd_delay_samples"]))
            mat[i, j] = float(r["gain_vs_alg1_db"])
        im = ax.imshow(mat, origin="lower", aspect="auto", cmap="RdBu", vmin=-vmax, vmax=vmax)
        ax.set_title(f"DMRS spacing={sf} sc")
        ax.set_xticks(np.arange(len(cdd_delays)))
        ax.set_xticklabels(cdd_delays)
        ax.set_yticks(np.arange(len(delay_spreads)))
        ax.set_yticklabels([f"{x:g}" for x in delay_spreads])
        ax.set_xlabel("CDD delay (samples)")
        ax.grid(False)
        for i in range(len(delay_spreads)):
            for j in range(len(cdd_delays)):
                val = mat[i, j]
                if np.isfinite(val):
                    ax.text(j, i, f"{val:.1f}", ha="center", va="center", fontsize=7)
    axes[0].set_ylabel("Delay spread (ns)")
    cbar = fig.colorbar(im, ax=axes, shrink=0.92)
    cbar.set_label("Alg3 equal-total NMSE gain vs Alg1 (dB)")
    fig.suptitle(f"Same-overhead Alg3 gain map @ SNR={snr_db:g} dB")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out) / f"search_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    delay_spreads = parse_float_list(args.delay_spreads_ns)
    cdd_delays = parse_int_list(args.cdd_delays)
    dmrs_spacings = parse_int_list(args.dmrs_spacings)
    snrs = parse_float_list(args.snrs)
    basis_thresholds = parse_float_list(args.basis_thresholds)

    all_rows: List[Dict[str, object]] = []
    win_rows: List[Dict[str, object]] = []

    for delay_spread_ns in delay_spreads:
        for cdd_delay in cdd_delays:
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
                delays = [0, int(cdd_delay)]
                channel = ChannelConfig(delay_spread_ns=delay_spread_ns, max_delay_factor=8.0)
                rho_h = freq_corr_from_pdp(pdp, grid.n_fft, grid.n_sc - 1)
                pdp_cdd = shifted_pdp(pdp, delays)
                rho_g = freq_corr_from_pdp(pdp_cdd, grid.n_fft, grid.n_sc - 1)
                bc_h_09 = first_below(rho_h, 0.9)
                bc_h_05 = first_below(rho_h, 0.5)
                bc_g_09 = first_below(rho_g, 0.9)
                bc_g_05 = first_below(rho_g, 0.5)
                pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)

                for snr_db in snrs:
                    noise_var = 10.0 ** (-float(snr_db) / 10.0)
                    noise_var_ls = noise_var / float(len(resource.dmrs_symbol_indices))
                    direct = make_direct_estimator(grid, pdp, delays, noise_var_ls, float(args.loading))
                    basis_ests = [
                        make_basis_estimator(grid, pdp, delays, noise_var_ls, float(args.loading), thr)
                        for thr in basis_thresholds
                    ]
                    block = make_block_estimator(grid, pdp, delays, noise_var_ls, float(args.loading), 0.9)
                    port_total = make_port_estimators(
                        grid,
                        pdp,
                        n_tx=2,
                        noise_var_ls=noise_var_ls,
                        loading=float(args.loading),
                        mode="equal_total_overhead",
                    )

                    accum: Dict[str, Dict[str, object]] = {}
                    rng = np.random.default_rng(stable_seed(delay_spread_ns, cdd_delay, dmrs_spacing, snr_db, base=int(args.seed)))
                    for _ in range(int(args.trials)):
                        realization = generate_tdl_channel(rng, grid, channel, n_tx=2, n_rx=4)
                        H = realization.H
                        g = cdd_equivalent_from_branches(H, grid, delays)
                        noise = (
                            rng.normal(size=(4, len(pilot_local)))
                            + 1j * rng.normal(size=(4, len(pilot_local)))
                        ) * math.sqrt(noise_var_ls / 2.0)
                        ls_cdd = g[:, pilot_local] + noise

                        g_hat = ls_cdd @ direct.matrix.T
                        add_result(accum, direct.label, direct.family, nmse(g, g_hat), float("nan"), direct.metadata)

                        for best in basis_ests:
                            g_hat, h_hat = apply_basis(best, ls_cdd, grid, delays)
                            meta = dict(best.metadata)
                            meta.update({"cond_number": best.cond_number, "effective_rank": best.effective_rank})
                            add_result(
                                accum,
                                best.label,
                                "alg2_basis_reduced_support",
                                nmse(g, g_hat),
                                nmse(H, h_hat),
                                meta,
                            )

                        g_hat, h_hat = apply_block(block, ls_cdd, grid, delays)
                        meta = dict(block.metadata)
                        meta.update({"cond_number": block.cond_max, "cond_mean": block.cond_mean})
                        add_result(
                            accum,
                            block.label,
                            "alg2_coherence_block",
                            nmse(g, g_hat),
                            nmse(H, h_hat),
                            meta,
                        )

                        g_hat, h_hat = apply_port_dmrs(port_total, H, grid, delays, noise_var_ls, rng, mode="equal_total_overhead")
                        add_result(
                            accum,
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

                    direct_nmse = mean_or_nan(accum[direct.label]["eff"])
                    scenario_rows: List[Dict[str, object]] = []
                    for label, values in accum.items():
                        meta = dict(values["meta"])
                        eff = mean_or_nan(values["eff"])
                        branch = mean_or_nan(values["branch"])
                        row = {
                            "algorithm": label,
                            "family": values["family"],
                            "delay_spread_ns": float(delay_spread_ns),
                            "cdd_delay_samples": int(cdd_delay),
                            "cdd_delay_ns": float(cdd_delay) * sample_period_ns,
                            "dmrs_spacing_sc": int(dmrs_spacing),
                            "snr_db": float(snr_db),
                            "trials": int(args.trials),
                            "n_sc": int(grid.n_sc),
                            "n_prbs": int(resource.n_prbs),
                            "pilot_count_combined": int(grid.pilot_count),
                            "dmrs_symbols": int(len(resource.dmrs_symbol_indices)),
                            "total_dmrs_re_alg1": int(grid.pilot_count * len(resource.dmrs_symbol_indices)),
                            "dmrs_overhead_alg1": float(grid.n_dmrs_re) / float(grid.n_sc * grid.n_symbols),
                            "bc_h_0p9_sc": "" if bc_h_09 is None else int(bc_h_09),
                            "bc_h_0p5_sc": "" if bc_h_05 is None else int(bc_h_05),
                            "bc_g_0p9_sc": "" if bc_g_09 is None else int(bc_g_09),
                            "bc_g_0p5_sc": "" if bc_g_05 is None else int(bc_g_05),
                            "ce_nmse_eff_mean": eff,
                            "ce_nmse_eff_median": median_or_nan(values["eff"]),
                            "ce_nmse_branch_mean": branch,
                            "ce_nmse_branch_median": median_or_nan(values["branch"]),
                            "gain_vs_alg1_db": gain_db(direct_nmse, eff),
                        }
                        row.update(meta)
                        scenario_rows.append(row)
                    all_rows.extend(scenario_rows)
                    win_rows.extend([
                        r for r in scenario_rows
                        if r["algorithm"] != direct.label and np.isfinite(float(r["gain_vs_alg1_db"])) and float(r["gain_vs_alg1_db"]) > 0.0
                    ])

    write_csv(all_rows, out_dir / "same_overhead_nmse_search.csv")
    win_rows.sort(key=lambda r: float(r["gain_vs_alg1_db"]), reverse=True)
    write_csv(win_rows, out_dir / "same_overhead_wins.csv")
    with (out_dir / "run_config.json").open("w", encoding="utf-8") as f:
        json.dump({
            "delay_spreads_ns": delay_spreads,
            "cdd_delays": cdd_delays,
            "dmrs_spacings": dmrs_spacings,
            "snrs": snrs,
            "basis_thresholds": basis_thresholds,
            "trials": int(args.trials),
            "seed": int(args.seed),
            "same_overhead_definition": "Alg3 uses equal-total FDM split, sum_m |P_m| = |P|; no equal-per-port extra-overhead curve is included.",
        }, f, indent=2)

    plot_best_gains(all_rows, fig_dir / "experiment10_best_same_overhead_gains.png")
    for snr in snrs:
        if abs(snr - round(snr)) < 1e-9:
            suffix = str(int(round(snr))).replace("-", "m")
        else:
            suffix = str(snr).replace("-", "m").replace(".", "p")
        plot_alg3_gain_maps(all_rows, snr, fig_dir / f"experiment10_alg3_gain_map_snr{suffix}.png")

    print(out_dir)
    return out_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/nmse_same_overhead_search")
    parser.add_argument("--fig-dir", default="docs/figures")
    parser.add_argument("--delay-spreads-ns", default="1,5,10,30,100")
    parser.add_argument("--cdd-delays", default="16,32,64,128,256")
    parser.add_argument("--dmrs-spacings", default="2,4,6,12")
    parser.add_argument("--snrs", default="-6,0,8")
    parser.add_argument("--basis-thresholds", default="0.90,0.95,0.99")
    parser.add_argument("--trials", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20260612)
    parser.add_argument("--loading", type=float, default=1e-8)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
