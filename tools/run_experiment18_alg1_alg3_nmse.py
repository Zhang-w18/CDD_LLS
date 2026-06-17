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
from cdd_lls.phy.estimators import nmse
from cdd_lls.phy.precoding import cdd_equivalent_from_branches
from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
from tools.run_experiment17_cdd_diversity_bler import stable_seed, write_csv
from tools.run_nmse_algorithm_comparison import make_direct_estimator, make_port_estimators, mean_or_nan


def parse_float_list(text: str) -> List[float]:
    return [float(x.strip()) for x in str(text).split(",") if x.strip()]


def parse_int_list(text: str) -> List[int]:
    return [int(float(x.strip())) for x in str(text).split(",") if x.strip()]


def resource(args: argparse.Namespace) -> ResourceConfig:
    return ResourceConfig(
        carrier_bandwidth_mhz=100.0,
        scs_khz=30,
        n_fft=int(args.n_fft),
        n_prbs=int(args.n_prbs),
        pdsch_n_symbols=int(args.pdsch_symbols),
        dmrs_symbol_indices=parse_int_list(args.dmrs_symbols),
        dmrs_spacing_sc=int(args.dmrs_spacing_sc),
        prg_size_rb=4,
    )


def run_nmse(args: argparse.Namespace, out_dir: Path) -> List[Dict[str, object]]:
    cdd_delays = parse_int_list(args.cdd_delays)
    snrs = parse_float_list(args.snrs)
    res = resource(args)
    grid = build_resource_grid(res)
    sample_period_ns = 1e9 / (float(res.n_fft) * float(res.scs_khz) * 1e3)
    pdp = make_exponential_pdp(float(args.delay_spread_ns), sample_period_ns, float(args.max_delay_factor))
    channel_cfg = ChannelConfig(delay_spread_ns=float(args.delay_spread_ns), max_delay_factor=float(args.max_delay_factor))
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    rows: List[Dict[str, object]] = []

    for snr_db in snrs:
        noise_var = 10.0 ** (-float(snr_db) / 10.0)
        noise_var_ls = noise_var / float(len(res.dmrs_symbol_indices))
        direct_estimators = {
            int(d): make_direct_estimator(
                grid=grid,
                pdp=pdp,
                delays=[0, int(d)],
                noise_var_ls=noise_var_ls,
                loading=float(args.loading),
            )
            for d in cdd_delays
        }
        port_estimators = make_port_estimators(
            grid=grid,
            pdp=pdp,
            n_tx=2,
            noise_var_ls=noise_var_ls,
            loading=float(args.loading),
            mode="equal_total_overhead",
        )
        acc = {
            (int(d), alg): []
            for d in cdd_delays
            for alg in ("Alg1 Direct RMMSE WB", "Alg3 Port-DMRS RMMSE equal-total")
        }
        branch_nmse_values: List[float] = []
        print(f"[alg1-alg3-nmse] SNR={float(snr_db):g} dB, delays={cdd_delays}", flush=True)
        rng = np.random.default_rng(stable_seed("alg1-alg3-nmse", args.delay_spread_ns, snr_db, base=int(args.seed)))
        for trial in range(1, int(args.trials) + 1):
            realization = generate_tdl_channel(rng, grid, channel_cfg, n_tx=2, n_rx=4)

            # Alg1 observes the CDD-combined channel on the combined pilot comb.
            noise_alg1 = (
                rng.normal(size=(4, len(pilot_local)))
                + 1j * rng.normal(size=(4, len(pilot_local)))
            ) * math.sqrt(noise_var_ls / 2.0)

            # Alg3 estimates physical port channels once, then recombines for each CDD delay.
            h_hat = np.empty_like(realization.H)
            for m, est in enumerate(port_estimators):
                pilot_k = np.asarray(est.metadata["pilot_k"], dtype=np.int64)
                pl = local_indices_for_subcarriers(grid, pilot_k)
                noise = (
                    rng.normal(size=(4, len(pl)))
                    + 1j * rng.normal(size=(4, len(pl)))
                ) * math.sqrt(noise_var_ls / 2.0)
                obs = realization.H[:, m, pl] + noise
                h_hat[:, m, :] = obs @ est.matrix.T
            branch_nmse_values.append(nmse(realization.H, h_hat))

            for d in cdd_delays:
                delays = [0, int(d)]
                true_g = cdd_equivalent_from_branches(realization.H, grid, delays)
                ls_cdd = true_g[:, pilot_local] + noise_alg1
                g_hat_alg1 = ls_cdd @ direct_estimators[int(d)].matrix.T
                g_hat_alg3 = cdd_equivalent_from_branches(h_hat, grid, delays)
                acc[(int(d), "Alg1 Direct RMMSE WB")].append(nmse(true_g, g_hat_alg1))
                acc[(int(d), "Alg3 Port-DMRS RMMSE equal-total")].append(nmse(true_g, g_hat_alg3))

            if int(args.progress_every) > 0 and trial % int(args.progress_every) == 0:
                print(
                    f"[alg1-alg3-nmse-progress] SNR={float(snr_db):g}, "
                    f"trial={trial}/{int(args.trials)}",
                    flush=True,
                )

        for d in cdd_delays:
            for alg in ("Alg1 Direct RMMSE WB", "Alg3 Port-DMRS RMMSE equal-total"):
                row = {
                    "algorithm": alg,
                    "cdd_delay_samples": int(d),
                    "cdd_delay_ns": float(d) * sample_period_ns,
                    "delay_spread_ns": float(args.delay_spread_ns),
                    "dmrs_spacing_sc": int(args.dmrs_spacing_sc),
                    "snr_db": float(snr_db),
                    "trials": int(args.trials),
                    "ce_nmse_eff_mean": mean_or_nan(acc[(int(d), alg)]),
                    "processing_bandwidth_sc": int(grid.n_sc),
                    "pilot_count_combined": int(grid.pilot_count),
                    "alg3_pilot_count_per_port": int(len(grid.pilot_subcarriers[0::2])),
                }
                if alg.startswith("Alg3"):
                    row["ce_nmse_branch_mean"] = mean_or_nan(branch_nmse_values)
                rows.append(row)
        write_csv(rows, out_dir / "alg1_alg3_nmse_by_cdd_partial.csv")

    write_csv(rows, out_dir / "alg1_alg3_nmse_by_cdd.csv")
    return rows


def plot_nmse_by_snr(rows: List[Dict[str, object]], path: Path, selected_delays: List[int]) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.9), sharey=True)
    for ax, alg in zip(axes, ["Alg1 Direct RMMSE WB", "Alg3 Port-DMRS RMMSE equal-total"]):
        for d in selected_delays:
            pts = sorted(
                [r for r in rows if int(r["cdd_delay_samples"]) == int(d) and str(r["algorithm"]) == alg],
                key=lambda r: float(r["snr_db"]),
            )
            if not pts:
                continue
            ax.semilogy(
                [float(r["snr_db"]) for r in pts],
                [float(r["ce_nmse_eff_mean"]) for r in pts],
                marker="o",
                linewidth=1.25,
                label=f"d={d}",
            )
        ax.set_title(alg, fontsize=10)
        ax.set_xlabel("SNR (dB)")
        ax.grid(True, which="both", alpha=0.3)
    axes[0].set_ylabel("Equivalent-channel NMSE")
    axes[-1].legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_nmse_vs_cdd(rows: List[Dict[str, object]], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    available_snrs = sorted({float(r["snr_db"]) for r in rows})
    selected_snrs = []
    if available_snrs:
        selected_snrs.append(available_snrs[0])
    if 5.0 in available_snrs:
        selected_snrs.append(5.0)
    elif available_snrs:
        selected_snrs.append(available_snrs[-1])
    selected_snrs = selected_snrs[:2]

    fig, axes = plt.subplots(1, len(selected_snrs), figsize=(6.2 * len(selected_snrs), 4.9), sharey=True)
    if len(selected_snrs) == 1:
        axes = [axes]
    for ax, snr in zip(axes, selected_snrs):
        for alg in ["Alg1 Direct RMMSE WB", "Alg3 Port-DMRS RMMSE equal-total"]:
            pts = sorted(
                [r for r in rows if str(r["algorithm"]) == alg and abs(float(r["snr_db"]) - snr) < 1e-9],
                key=lambda r: int(r["cdd_delay_samples"]),
            )
            if not pts:
                continue
            ax.semilogy(
                [int(r["cdd_delay_samples"]) for r in pts],
                [float(r["ce_nmse_eff_mean"]) for r in pts],
                marker="o",
                linewidth=1.25,
                label=alg,
            )
        ax.set_title(f"SNR={snr:g} dB", fontsize=10)
        ax.set_xlabel("CDD delay (samples)")
        ax.grid(True, which="both", alpha=0.3)
    axes[0].set_ylabel("Equivalent-channel NMSE")
    axes[-1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out) / f"experiment18_alg1_alg3_nmse_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    rows = run_nmse(args, out_dir)
    plot_nmse_by_snr(
        rows,
        fig_dir / "experiment18_alg1_alg3_nmse_d128_d512.png",
        parse_int_list(args.selected_plot_delays),
    )
    plot_nmse_vs_cdd(rows, fig_dir / "experiment18_alg1_alg3_nmse_vs_cdd.png")

    with (out_dir / "run_config.json").open("w", encoding="utf-8") as f:
        json.dump({
            **vars(args),
            "notes": [
                "Alg1 uses LMMSE on the CDD-combined equivalent channel with CDD-shifted PDP.",
                "Alg3 uses per-port LMMSE on physical branch channels and recombines each CDD delay from the same branch estimate.",
                "Both target the full 576-sc allocation; equal-total Alg3 uses half the combined pilots per port.",
            ],
        }, f, indent=2, ensure_ascii=False)
    print(out_dir)
    return out_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/experiment18_alg1_alg3_nmse")
    parser.add_argument("--fig-dir", default="docs/figures")
    parser.add_argument("--delay-spread-ns", type=float, default=10.0)
    parser.add_argument("--max-delay-factor", type=float, default=8.0)
    parser.add_argument("--cdd-delays", default="0,64,128,256,512,768,1024,1536,2048")
    parser.add_argument("--selected-plot-delays", default="128,512")
    parser.add_argument("--dmrs-spacing-sc", type=int, default=24)
    parser.add_argument("--snrs", default="2,3,4,5,6")
    parser.add_argument("--n-prbs", type=int, default=48)
    parser.add_argument("--n-fft", type=int, default=4096)
    parser.add_argument("--pdsch-symbols", type=int, default=10)
    parser.add_argument("--dmrs-symbols", default="2,7")
    parser.add_argument("--trials", type=int, default=200)
    parser.add_argument("--loading", type=float, default=1e-8)
    parser.add_argument("--seed", type=int, default=20260616)
    parser.add_argument("--progress-every", type=int, default=50)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
