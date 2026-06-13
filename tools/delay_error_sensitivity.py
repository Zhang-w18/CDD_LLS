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
    apply_port_dmrs,
    first_below,
    freq_corr_from_pdp,
    make_basis_estimator,
    make_direct_estimator,
    make_port_estimators,
    mean_or_nan,
    median_or_nan,
    write_csv,
)


def parse_int_list(text: str) -> List[int]:
    return [int(float(x.strip())) for x in str(text).split(",") if x.strip()]


def stable_seed(*values: object, base: int) -> int:
    h = int(base) & 0x7FFFFFFF
    for value in values:
        for ch in str(value):
            h = (h * 131 + ord(ch)) & 0x7FFFFFFF
    return h


def add(accum: Dict[str, Dict[str, object]], algorithm: str, eff: float, branch: float, meta: Dict[str, object]) -> None:
    if algorithm not in accum:
        accum[algorithm] = {"eff": [], "branch": [], "meta": meta}
    accum[algorithm]["eff"].append(float(eff))
    accum[algorithm]["branch"].append(float(branch))


def gain_db(reference: float, candidate: float) -> float:
    if not np.isfinite(reference) or not np.isfinite(candidate) or reference <= 0.0 or candidate <= 0.0:
        return float("nan")
    return float(10.0 * np.log10(reference / candidate))


def plot_delay_error(rows: List[Dict[str, object]], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cases = list(dict.fromkeys(str(r["case_id"]) for r in rows))
    algs = list(dict.fromkeys(str(r["algorithm"]) for r in rows))
    fig, axes = plt.subplots(len(cases), 1, figsize=(10.5, max(3.4 * len(cases), 5.0)), sharex=True)
    if len(cases) == 1:
        axes = [axes]
    for ax, case_id in zip(axes, cases):
        for alg in algs:
            pts = [r for r in rows if r["case_id"] == case_id and r["algorithm"] == alg]
            pts = sorted(pts, key=lambda r: int(r["delay_error_samples"]))
            if not pts:
                continue
            ax.semilogy(
                [int(r["delay_error_samples"]) for r in pts],
                [float(r["ce_nmse_eff_mean"]) for r in pts],
                marker="o",
                linewidth=1.3,
                markersize=3.5,
                label=alg,
            )
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_title(case_id, loc="left", fontsize=10)
        ax.set_ylabel("Effective NMSE")
        ax.grid(True, which="both", alpha=0.28)
    axes[-1].set_xlabel("Assumed CDD delay error on port 1 (samples)")
    axes[0].legend(fontsize=8, ncol=3, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out) / f"delay_error_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        {
            "case_id": "smooth_ds1_d128_sf6_snr8",
            "delay_spread_ns": 1.0,
            "cdd_delay": 128,
            "dmrs_spacing_sc": 6,
            "snr_db": 8.0,
        },
        {
            "case_id": "smooth_ds5_d64_sf4_snrm6",
            "delay_spread_ns": 5.0,
            "cdd_delay": 64,
            "dmrs_spacing_sc": 4,
            "snr_db": -6.0,
        },
        {
            "case_id": "tdl30_d64_sf6_snrm6",
            "delay_spread_ns": 30.0,
            "cdd_delay": 64,
            "dmrs_spacing_sc": 6,
            "snr_db": -6.0,
        },
        {
            "case_id": "tdl100_d32_sf6_snr8",
            "delay_spread_ns": 100.0,
            "cdd_delay": 32,
            "dmrs_spacing_sc": 6,
            "snr_db": 8.0,
        },
    ]
    delay_errors = parse_int_list(args.delay_errors)
    rows: List[Dict[str, object]] = []

    for case in cases:
        resource = ResourceConfig(
            carrier_bandwidth_mhz=100.0,
            scs_khz=30,
            n_fft=4096,
            n_prbs=48,
            pdsch_n_symbols=10,
            dmrs_symbol_indices=[2, 7],
            dmrs_spacing_sc=int(case["dmrs_spacing_sc"]),
            prg_size_rb=4,
        )
        grid = build_resource_grid(resource)
        sample_period_ns = 1e9 / (float(resource.n_fft) * float(resource.scs_khz) * 1e3)
        pdp = make_exponential_pdp(float(case["delay_spread_ns"]), sample_period_ns, max_delay_factor=8.0)
        true_delays = [0, int(case["cdd_delay"])]
        rho_g = freq_corr_from_pdp(shifted_pdp(pdp, true_delays), grid.n_fft, grid.n_sc - 1)
        rho_h = freq_corr_from_pdp(pdp, grid.n_fft, grid.n_sc - 1)
        bc_h_05 = first_below(rho_h, 0.5)
        bc_g_05 = first_below(rho_g, 0.5)
        noise_var = 10.0 ** (-float(case["snr_db"]) / 10.0)
        noise_var_ls = noise_var / float(len(resource.dmrs_symbol_indices))
        channel = ChannelConfig(delay_spread_ns=float(case["delay_spread_ns"]), max_delay_factor=8.0)
        pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)

        for err in delay_errors:
            assumed_delays = [0, int(case["cdd_delay"]) + int(err)]
            accum: Dict[str, Dict[str, object]] = {}
            direct = make_direct_estimator(grid, pdp, assumed_delays, noise_var_ls, float(args.loading))
            basis = make_basis_estimator(grid, pdp, assumed_delays, noise_var_ls, float(args.loading), 0.99)
            port_total = make_port_estimators(
                grid,
                pdp,
                n_tx=2,
                noise_var_ls=noise_var_ls,
                loading=float(args.loading),
                mode="equal_total_overhead",
            )
            rng = np.random.default_rng(
                stable_seed(case["case_id"], err, base=int(args.seed))
            )
            for _ in range(int(args.trials)):
                realization = generate_tdl_channel(rng, grid, channel, n_tx=2, n_rx=4)
                H = realization.H
                true_g = cdd_equivalent_from_branches(H, grid, true_delays)
                noise = (
                    rng.normal(size=(4, len(pilot_local)))
                    + 1j * rng.normal(size=(4, len(pilot_local)))
                ) * math.sqrt(noise_var_ls / 2.0)
                ls_cdd = true_g[:, pilot_local] + noise

                g_hat = ls_cdd @ direct.matrix.T
                add(accum, "Alg1 Direct RMMSE", nmse(true_g, g_hat), float("nan"), direct.metadata)

                g_hat, h_hat = apply_basis(basis, ls_cdd, grid, assumed_delays)
                add(
                    accum,
                    "Alg2 Basis E99",
                    nmse(true_g, g_hat),
                    nmse(H, h_hat),
                    {
                        "support_len": int(len(basis.support)),
                        "cond_number": float(basis.cond_number),
                        "effective_rank": float(basis.effective_rank),
                    },
                )

                g_hat, h_hat = apply_port_dmrs(
                    port_total,
                    H,
                    grid,
                    assumed_delays,
                    noise_var_ls,
                    rng,
                    mode="equal_total_overhead",
                )
                add(
                    accum,
                    "Alg3 Port-DMRS equal-total",
                    nmse(true_g, g_hat),
                    nmse(H, h_hat),
                    {
                        "pilot_count_per_port": int(len(grid.pilot_subcarriers[0::2])),
                        "total_pilot_re_relative": 1.0,
                    },
                )

            ref = mean_or_nan(accum["Alg1 Direct RMMSE"]["eff"])
            for algorithm, values in accum.items():
                eff = mean_or_nan(values["eff"])
                branch = mean_or_nan(values["branch"])
                row = {
                    "case_id": case["case_id"],
                    "algorithm": algorithm,
                    "delay_spread_ns": float(case["delay_spread_ns"]),
                    "true_cdd_delay_samples": int(case["cdd_delay"]),
                    "assumed_cdd_delay_samples": int(case["cdd_delay"]) + int(err),
                    "delay_error_samples": int(err),
                    "delay_error_ns": float(err) * sample_period_ns,
                    "dmrs_spacing_sc": int(case["dmrs_spacing_sc"]),
                    "snr_db": float(case["snr_db"]),
                    "trials": int(args.trials),
                    "pilot_count_combined": int(grid.pilot_count),
                    "bc_h_0p5_sc": "" if bc_h_05 is None else int(bc_h_05),
                    "bc_g_0p5_sc": "" if bc_g_05 is None else int(bc_g_05),
                    "ce_nmse_eff_mean": eff,
                    "ce_nmse_eff_median": median_or_nan(values["eff"]),
                    "ce_nmse_branch_mean": branch,
                    "ce_nmse_branch_median": median_or_nan(values["branch"]),
                    "gain_vs_alg1_same_error_db": gain_db(ref, eff),
                }
                row.update(values["meta"])
                rows.append(row)

    write_csv(rows, out_dir / "delay_error_sensitivity.csv")
    with (out_dir / "run_config.json").open("w", encoding="utf-8") as f:
        json.dump({"cases": cases, "delay_errors": delay_errors, "trials": int(args.trials)}, f, indent=2)
    plot_delay_error(rows, fig_dir / "experiment11_delay_error_sensitivity.png")
    print(out_dir)
    return out_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/delay_error_sensitivity")
    parser.add_argument("--fig-dir", default="docs/figures")
    parser.add_argument("--delay-errors", default="-8,-4,-2,-1,0,1,2,4,8")
    parser.add_argument("--trials", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260612)
    parser.add_argument("--loading", type=float, default=1e-8)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
