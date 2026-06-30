from __future__ import annotations

import argparse
import csv
import json
import math
import os
from datetime import datetime
from pathlib import Path
import sys
from typing import Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ChannelConfig, ResourceConfig
from cdd_lls.core.mcs import build_tb_layout, get_mcs
from cdd_lls.phy.channel_tdl import generate_tdl_channel, make_exponential_pdp
from cdd_lls.phy.estimators import nmse
from cdd_lls.phy.ldpc import SionnaLDPCAdapter
from cdd_lls.phy.precoding import cdd_equivalent_from_branches
from cdd_lls.phy.qam import qam_modulate
from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
from cdd_lls.sim.stats import interpolate_target_snr
from tools.run_experiment17_cdd_diversity_bler import (
    decode_same_tb_batch,
    llr_for_candidate,
    stable_seed,
    write_csv,
)
from tools.run_nmse_algorithm_comparison import (
    apply_port_dmrs,
    make_direct_estimator,
    make_port_estimators,
)


ALG1 = "Alg1 Direct RMMSE WB"
ALG3 = "Alg3 Port-DMRS RMMSE equal-total"


def parse_float_list(text: str) -> List[float]:
    return [float(x.strip()) for x in str(text).split(",") if x.strip()]


def parse_int_list(text: str) -> List[int]:
    return [int(float(x.strip())) for x in str(text).split(",") if x.strip()]


def read_csv(path: Path) -> List[Dict[str, object]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return [dict(row) for row in csv.DictReader(f)]


def result_key(row: Dict[str, object]) -> tuple[float, int, int, int, float, str]:
    return (
        float(row["delay_spread_ns"]),
        int(float(row["cdd_delay_samples"])),
        int(float(row["dmrs_spacing_sc"])),
        int(float(row["mcs_index"])),
        float(row["snr_db"]),
        str(row["algorithm"]),
    )


def make_resource(args: argparse.Namespace) -> ResourceConfig:
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


def target_summary(rows: List[Dict[str, object]], targets: List[float]) -> List[Dict[str, object]]:
    out: List[Dict[str, object]] = []
    delay_spreads = sorted({float(row["delay_spread_ns"]) for row in rows})
    for delay_spread in delay_spreads:
        spacings = sorted({
            int(float(row["dmrs_spacing_sc"]))
            for row in rows
            if float(row["delay_spread_ns"]) == delay_spread
        })
        for spacing in spacings:
            for algorithm in (ALG1, ALG3):
                points = [
                    row for row in rows
                    if float(row["delay_spread_ns"]) == delay_spread
                    and int(float(row["dmrs_spacing_sc"])) == spacing
                    and str(row["algorithm"]) == algorithm
                ]
                for target in targets:
                    out.append({
                        "delay_spread_ns": delay_spread,
                        "cdd_delay_samples": int(float(points[0]["cdd_delay_samples"])),
                        "dmrs_spacing_sc": spacing,
                        "algorithm": algorithm,
                        "target_bler": float(target),
                        "snr_at_target_bler_db": interpolate_target_snr(points, target_bler=float(target)),
                        "min_trials_per_snr": min(int(float(point["n_trials"])) for point in points),
                        "max_trials_per_snr": max(int(float(point["n_trials"])) for point in points),
                        "mcs_index": int(float(points[0]["mcs_index"])),
                        "qm": int(float(points[0]["qm"])),
                        "code_rate": float(points[0]["code_rate"]),
                    })
    return out


def plot_bler(rows: List[Dict[str, object]], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator

    colors = {
        6: "#0072B2",
        12: "#D55E00",
        18: "#C44E52",
        24: "#009E73",
    }
    styles = {
        ALG1: {"linestyle": (0, (5, 2))},
        ALG3: {"linestyle": "-"},
    }
    ds_markers = {
        1.0: "o",
        10.0: "s",
        100.0: "^",
    }
    max_trials = max(int(float(row["n_trials"])) for row in rows)
    zero_floor = 0.5 / float(max_trials)
    fig, ax = plt.subplots(figsize=(11.5, 6.8))
    for spacing in sorted({int(float(row["dmrs_spacing_sc"])) for row in rows}):
        for delay_spread in sorted({
            float(row["delay_spread_ns"])
            for row in rows
            if int(float(row["dmrs_spacing_sc"])) == spacing
        }):
            for algorithm in (ALG1, ALG3):
                points = sorted(
                    [
                        row for row in rows
                        if int(float(row["dmrs_spacing_sc"])) == spacing
                        and float(row["delay_spread_ns"]) == delay_spread
                        and str(row["algorithm"]) == algorithm
                    ],
                    key=lambda row: float(row["snr_db"]),
                )
                if not points:
                    continue
                ys = [
                    max(float(row["bler"]), 0.5 / float(row["n_trials"]))
                    for row in points
                ]
                label_alg = "Alg1" if algorithm == ALG1 else "Alg3"
                ax.semilogy(
                    [float(row["snr_db"]) for row in points],
                    ys,
                    color=colors.get(spacing, "#333333"),
                    linewidth=1.8,
                    markersize=5.6,
                    linestyle=styles[algorithm]["linestyle"],
                    marker=ds_markers.get(delay_spread, "o"),
                    markerfacecolor=colors.get(spacing, "#333333"),
                    markeredgewidth=1.2,
                    label=f"sp={spacing} | DS={delay_spread:g} | {label_alg}",
                )
    ax.axhline(0.01, color="#555555", linewidth=1.0, linestyle=":", label="BLER=0.01")
    ax.set_ylim(max(zero_floor / 1.6, 1e-4), 1.08)
    ax.set_xlabel("SNR (dB)")
    ax.set_ylabel("TB BLER")
    ax.set_title("Alg1 vs Alg3, MCS20 / 256QAM (CDD=512)")
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    spacing24_alg1 = [
        row for row in rows
        if int(float(row["dmrs_spacing_sc"])) == 24 and str(row["algorithm"]) == ALG1
    ]
    if spacing24_alg1 and all(abs(float(row["bler"]) - 1.0) < 1e-12 for row in spacing24_alg1):
        ax.text(
            0.985,
            0.955,
            "All spacing=24 Alg1 curves overlap at BLER=1",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8.3,
            bbox={"facecolor": "white", "edgecolor": "#999999", "alpha": 0.9, "pad": 3.0},
        )
    ax.grid(True, which="major", axis="x", color="#777777", linewidth=0.7, alpha=0.36)
    ax.grid(True, which="minor", axis="x", color="#999999", linewidth=0.55, linestyle=":", alpha=0.24)
    ax.grid(True, which="both", axis="y", alpha=0.24)
    ax.legend(ncol=2, fontsize=7.4, loc="best", handlelength=3.2, markerscale=1.05)
    fig.tight_layout()
    fig.savefig(path, dpi=190)
    plt.close(fig)


def run(args: argparse.Namespace) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out) / f"experiment11_alg3_gain_bler_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    delay_spreads = parse_float_list(args.delay_spreads_ns)
    snrs = parse_float_list(args.snrs)
    resource = make_resource(args)
    grid = build_resource_grid(resource)
    mcs = get_mcs("nr_256qam", int(args.mcs_index), None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(
        tb.cb_k_values,
        tb.cb_e_values,
        num_iter=int(args.ldpc_iterations),
        llr_clip=float(args.llr_clip),
    )
    sample_period_ns = 1e9 / (float(resource.n_fft) * float(resource.scs_khz) * 1e3)
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    delays = [0, int(args.cdd_delay)]
    rows: List[Dict[str, object]] = []
    for path_text in args.reuse_csv:
        path = Path(path_text)
        if not path.exists():
            raise FileNotFoundError(f"Reuse CSV does not exist: {path}")
        rows.extend(read_csv(path))
    rows_by_key = {result_key(row): row for row in rows}

    for delay_spread in delay_spreads:
        pdp = make_exponential_pdp(
            float(delay_spread),
            sample_period_ns,
            max_delay_factor=float(args.max_delay_factor),
        )
        channel_cfg = ChannelConfig(
            delay_spread_ns=float(delay_spread),
            max_delay_factor=float(args.max_delay_factor),
        )
        for snr_db in snrs:
            scenario_keys = {
                algorithm: (
                    float(delay_spread),
                    int(args.cdd_delay),
                    int(args.dmrs_spacing_sc),
                    int(args.mcs_index),
                    float(snr_db),
                    algorithm,
                )
                for algorithm in (ALG1, ALG3)
            }
            existing = {
                algorithm: rows_by_key.get(key)
                for algorithm, key in scenario_keys.items()
            }
            existing_trials = {
                algorithm: 0 if row is None else int(float(row["n_trials"]))
                for algorithm, row in existing.items()
            }
            if all(count >= int(args.trials) for count in existing_trials.values()):
                print(f"[exp11-reuse] DS={delay_spread:g} ns, SNR={snr_db:g} dB", flush=True)
                continue
            nonzero_counts = {count for count in existing_trials.values() if count > 0}
            if len(nonzero_counts) > 1 or (nonzero_counts and any(count == 0 for count in existing_trials.values())):
                raise ValueError(
                    f"Inconsistent reusable trial counts for DS={delay_spread}, SNR={snr_db}: {existing_trials}"
                )
            start_trial = next(iter(nonzero_counts), 0)
            noise_var = 10.0 ** (-float(snr_db) / 10.0)
            noise_var_ls = noise_var / float(len(resource.dmrs_symbol_indices))
            direct = make_direct_estimator(
                grid,
                pdp,
                delays,
                noise_var_ls,
                float(args.loading),
            )
            port_total = make_port_estimators(
                grid,
                pdp,
                n_tx=2,
                noise_var_ls=noise_var_ls,
                loading=float(args.loading),
                mode="equal_total_overhead",
            )
            accum = {}
            for algorithm in (ALG1, ALG3):
                prior = existing[algorithm]
                prior_trials = existing_trials[algorithm]
                accum[algorithm] = {
                    "tb_errors": 0 if prior is None else int(float(prior["tb_errors"])),
                    "cb_errors": 0 if prior is None else int(float(prior["cb_errors"])),
                    "ce_nmse_sum": 0.0 if prior is None else float(prior["ce_nmse_eff"]) * prior_trials,
                    "ce_nmse_count": prior_trials,
                }
            print(
                f"[exp11-bler] DS={delay_spread:g} ns, SNR={snr_db:g} dB, "
                f"MCS={int(args.mcs_index)}, trials={start_trial + 1}..{int(args.trials)}",
                flush=True,
            )
            for trial in range(start_trial + 1, int(args.trials) + 1):
                rng = np.random.default_rng(
                    stable_seed(
                        "exp11-alg3-gain-bler",
                        delay_spread,
                        args.cdd_delay,
                        args.dmrs_spacing_sc,
                        args.mcs_index,
                        snr_db,
                        trial,
                        base=int(args.seed),
                    )
                )
                realization = generate_tdl_channel(rng, grid, channel_cfg, n_tx=2, n_rx=4)
                true_g = cdd_equivalent_from_branches(realization.H, grid, delays)
                payload = [rng.integers(0, 2, size=int(k), dtype=np.int8) for k in tb.cb_k_values]
                symbols = qam_modulate(np.concatenate(adapter.encode(payload)), int(mcs.qm))
                data_noise = (
                    rng.normal(size=(4, len(data_local)))
                    + 1j * rng.normal(size=(4, len(data_local)))
                ) * math.sqrt(noise_var / 2.0)

                pilot_noise_alg1 = (
                    rng.normal(size=(4, len(pilot_local)))
                    + 1j * rng.normal(size=(4, len(pilot_local)))
                ) * math.sqrt(noise_var_ls / 2.0)
                ls_cdd = true_g[:, pilot_local] + pilot_noise_alg1
                g_hat_alg1 = ls_cdd @ direct.matrix.T
                g_hat_alg3, _ = apply_port_dmrs(
                    port_total,
                    realization.H,
                    grid,
                    delays,
                    noise_var_ls,
                    rng,
                    mode="equal_total_overhead",
                )

                g_hats = [g_hat_alg1, g_hat_alg3]
                llrs = [
                    llr_for_candidate(
                        true_g,
                        g_hat,
                        data_local,
                        symbols,
                        data_noise,
                        noise_var,
                        int(mcs.qm),
                    )
                    for g_hat in g_hats
                ]
                for algorithm, decoded, g_hat in zip(
                    (ALG1, ALG3),
                    decode_same_tb_batch(adapter, llrs, payload),
                    g_hats,
                ):
                    accum[algorithm]["tb_errors"] += int(not decoded.tb_success)
                    accum[algorithm]["cb_errors"] += int(sum(1 for ok in decoded.cb_success if not ok))
                    accum[algorithm]["ce_nmse_sum"] += nmse(true_g, g_hat)
                    accum[algorithm]["ce_nmse_count"] += 1

                if int(args.progress_every) > 0 and trial % int(args.progress_every) == 0:
                    status = ", ".join(
                        f"{algorithm.split()[0]}={accum[algorithm]['tb_errors'] / float(trial):.4f}"
                        for algorithm in (ALG1, ALG3)
                    )
                    print(
                        f"[exp11-progress] DS={delay_spread:g}, SNR={snr_db:g}, "
                        f"trial={trial}/{int(args.trials)} | {status}",
                        flush=True,
                    )

            scenario_key_set = set(scenario_keys.values())
            rows = [row for row in rows if result_key(row) not in scenario_key_set]
            for algorithm in (ALG1, ALG3):
                item = accum[algorithm]
                tb_errors = int(item["tb_errors"])
                cb_errors = int(item["cb_errors"])
                rows.append({
                    "algorithm": algorithm,
                    "delay_spread_ns": float(delay_spread),
                    "cdd_delay_samples": int(args.cdd_delay),
                    "dmrs_spacing_sc": int(args.dmrs_spacing_sc),
                    "snr_db": float(snr_db),
                    "n_trials": int(args.trials),
                    "tb_errors": tb_errors,
                    "cb_errors": cb_errors,
                    "bler": float(tb_errors) / float(args.trials),
                    "cb_bler": float(cb_errors) / float(args.trials * len(tb.cb_k_values)),
                    "bler_resolution": 1.0 / float(args.trials),
                    "ce_nmse_eff": float(item["ce_nmse_sum"]) / float(item["ce_nmse_count"]),
                    "n_tx": 2,
                    "n_rx": 4,
                    "n_prbs": int(args.n_prbs),
                    "n_sc": int(grid.n_sc),
                    "n_fft": int(args.n_fft),
                    "pdsch_symbols": int(args.pdsch_symbols),
                    "dmrs_symbols": int(len(resource.dmrs_symbol_indices)),
                    "pilot_count_combined": int(grid.pilot_count),
                    "alg3_pilot_count_per_port": int(len(grid.pilot_subcarriers[0::2])),
                    "same_total_dmrs_overhead": True,
                    "mcs_table": "nr_256qam",
                    "mcs_index": int(args.mcs_index),
                    "qm": int(mcs.qm),
                    "code_rate": float(mcs.code_rate),
                    "actual_tbs_bits": int(tb.tb_size),
                    "coded_bits": int(tb.coded_bits),
                    "n_cbs": int(len(tb.cb_k_values)),
                    "ldpc_iterations": int(args.ldpc_iterations),
                })
            rows_by_key = {result_key(row): row for row in rows}
            write_csv(rows, out_dir / "alg1_alg3_high_mcs_bler_partial.csv")

    targets = target_summary(rows, parse_float_list(args.target_blers))
    write_csv(rows, out_dir / "alg1_alg3_high_mcs_bler.csv")
    write_csv(targets, out_dir / "alg1_alg3_high_mcs_target_snr.csv")
    plot_bler(rows, fig_dir / "experiment11_alg1_alg3_high_mcs_bler.png")
    merged_delay_spreads = sorted({float(row["delay_spread_ns"]) for row in rows})
    merged_snrs_by_ds_spacing = {
        f"DS={delay_spread:g},spacing={spacing}": sorted({
            float(row["snr_db"])
            for row in rows
            if float(row["delay_spread_ns"]) == delay_spread
            and int(float(row["dmrs_spacing_sc"])) == spacing
        })
        for delay_spread in merged_delay_spreads
        for spacing in sorted({
            int(float(row["dmrs_spacing_sc"]))
            for row in rows
            if float(row["delay_spread_ns"]) == delay_spread
        })
    }
    with (out_dir / "run_config.json").open("w", encoding="utf-8") as f:
        json.dump(
            {
                **vars(args),
                "mcs_qm": int(mcs.qm),
                "mcs_code_rate": float(mcs.code_rate),
                "actual_tbs_bits": int(tb.tb_size),
                "coded_bits": int(tb.coded_bits),
                "n_cbs": int(len(tb.cb_k_values)),
                "merged_delay_spreads_ns": merged_delay_spreads,
                "merged_snrs_by_ds_spacing": merged_snrs_by_ds_spacing,
                "merged_scenario_points": int(len(rows) // 2),
                "merged_algorithm_rows": int(len(rows)),
                "reuse_csv": list(args.reuse_csv),
                "randomization": (
                    "Each trial regenerates TDL channel, payload, and noise. Alg1 and Alg3 share "
                    "the underlying channel, payload, and data noise within a trial; their DMRS "
                    "observation noises follow their respective pilot designs."
                ),
            },
            f,
            indent=2,
        )
    print(out_dir)
    return out_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/experiment11_alg3_gain_bler")
    parser.add_argument("--fig-dir", default="docs/figures")
    parser.add_argument("--delay-spreads-ns", default="1,10,100")
    parser.add_argument("--max-delay-factor", type=float, default=8.0)
    parser.add_argument("--cdd-delay", type=int, default=512)
    parser.add_argument("--dmrs-spacing-sc", type=int, default=24)
    parser.add_argument("--snrs", default="8,10,12,14,16,18,20")
    parser.add_argument("--target-blers", default="0.1,0.01")
    parser.add_argument("--n-prbs", type=int, default=48)
    parser.add_argument("--n-fft", type=int, default=4096)
    parser.add_argument("--pdsch-symbols", type=int, default=10)
    parser.add_argument("--dmrs-symbols", default="2,7")
    parser.add_argument("--mcs-index", type=int, default=20)
    parser.add_argument("--trials", type=int, default=200)
    parser.add_argument("--ldpc-iterations", type=int, default=8)
    parser.add_argument("--llr-clip", type=float, default=50.0)
    parser.add_argument("--loading", type=float, default=1e-8)
    parser.add_argument("--seed", type=int, default=20260619)
    parser.add_argument("--progress-every", type=int, default=50)
    parser.add_argument("--reuse-csv", action="append", default=[])
    run(parser.parse_args())


if __name__ == "__main__":
    main()
