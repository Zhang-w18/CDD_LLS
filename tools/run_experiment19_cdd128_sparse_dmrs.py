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
from cdd_lls.core.mcs import build_tb_layout, get_mcs
from cdd_lls.phy.channel_tdl import generate_tdl_channel, make_exponential_pdp
from cdd_lls.phy.estimators import nmse, shifted_pdp
from cdd_lls.phy.ldpc import SionnaLDPCAdapter
from cdd_lls.phy.precoding import cdd_equivalent_from_branches
from cdd_lls.phy.qam import qam_modulate
from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
from cdd_lls.sim.stats import interpolate_target_snr
from tools.run_experiment17_cdd_diversity_bler import decode_same_tb_batch, llr_for_candidate, stable_seed
from tools.run_nmse_algorithm_comparison import (
    first_below,
    freq_corr_from_pdp,
    make_direct_estimator,
    make_port_estimators,
    mean_or_nan,
)


ALG1 = "Alg1 Direct RMMSE WB"
ALG3 = "Alg3 Port-DMRS RMMSE equal-total"
IDEAL = "Ideal CSI"


def parse_float_list(text: str) -> List[float]:
    return [float(x.strip()) for x in str(text).split(",") if x.strip()]


def parse_int_list(text: str) -> List[int]:
    return [int(float(x.strip())) for x in str(text).split(",") if x.strip()]


def write_csv(rows: List[Dict[str, object]], path: Path) -> None:
    if not rows:
        return
    fields = sorted(set().union(*(row.keys() for row in rows)))
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def resource_for_spacing(args: argparse.Namespace, dmrs_spacing: int) -> ResourceConfig:
    return ResourceConfig(
        carrier_bandwidth_mhz=100.0,
        scs_khz=30,
        n_fft=int(args.n_fft),
        n_prbs=int(args.n_prbs),
        pdsch_n_symbols=int(args.pdsch_symbols),
        dmrs_symbol_indices=parse_int_list(args.dmrs_symbols),
        dmrs_spacing_sc=int(dmrs_spacing),
        prg_size_rb=4,
    )


def reference_tbs_bits(grid, mcs) -> int:
    total_re = int(grid.n_sc) * int(grid.n_symbols)
    return int(build_tb_layout(total_re, mcs).tb_size)


def common_metadata(args: argparse.Namespace, grid, resource, spacing: int, sample_period_ns: float) -> Dict[str, object]:
    total_re = int(grid.n_sc) * int(grid.n_symbols)
    return {
        "delay_spread_ns": float(args.delay_spread_ns),
        "cdd_delay_samples": int(args.cdd_delay),
        "cdd_delay_ns": float(args.cdd_delay) * float(sample_period_ns),
        "dmrs_spacing_sc": int(spacing),
        "dmrs_symbols": int(len(resource.dmrs_symbol_indices)),
        "n_prbs": int(args.n_prbs),
        "n_sc": int(grid.n_sc),
        "n_fft": int(args.n_fft),
        "pdsch_symbols": int(args.pdsch_symbols),
        "pilot_count_combined": int(grid.pilot_count),
        "alg3_pilot_count_per_port": int(len(grid.pilot_subcarriers[0::2])),
        "total_re": int(total_re),
        "data_re": int(grid.n_data_re),
        "total_dmrs_re": int(grid.n_dmrs_re),
        "dmrs_overhead": float(grid.dmrs_overhead),
        "throughput_factor": float(grid.n_data_re) / float(total_re),
        "same_total_dmrs_overhead_alg1_alg3": True,
        "processing_bandwidth_sc": int(grid.n_sc),
    }


def run_nmse_sweep(args: argparse.Namespace, out_dir: Path) -> tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    spacings = parse_int_list(args.dmrs_spacings)
    snrs = parse_float_list(args.snrs)
    rows: List[Dict[str, object]] = []
    gain_rows: List[Dict[str, object]] = []
    channel_cfg = ChannelConfig(
        delay_spread_ns=float(args.delay_spread_ns),
        max_delay_factor=float(args.max_delay_factor),
    )

    for spacing in spacings:
        resource = resource_for_spacing(args, int(spacing))
        grid = build_resource_grid(resource)
        sample_period_ns = 1e9 / (float(resource.n_fft) * float(resource.scs_khz) * 1e3)
        pdp = make_exponential_pdp(float(args.delay_spread_ns), sample_period_ns, float(args.max_delay_factor))
        pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
        base = common_metadata(args, grid, resource, int(spacing), sample_period_ns)
        delays = [0, int(args.cdd_delay)]

        for snr_db in snrs:
            noise_var = 10.0 ** (-float(snr_db) / 10.0)
            noise_var_ls = noise_var / float(len(resource.dmrs_symbol_indices))
            direct = make_direct_estimator(
                grid=grid,
                pdp=pdp,
                delays=delays,
                noise_var_ls=noise_var_ls,
                loading=float(args.loading),
            )
            port_total = make_port_estimators(
                grid=grid,
                pdp=pdp,
                n_tx=2,
                noise_var_ls=noise_var_ls,
                loading=float(args.loading),
                mode="equal_total_overhead",
            )
            acc = {ALG1: [], ALG3: []}
            branch_acc: List[float] = []
            rng = np.random.default_rng(
                stable_seed("exp19-nmse", args.delay_spread_ns, args.cdd_delay, spacing, snr_db, base=int(args.seed))
            )
            print(
                f"[exp19-nmse] spacing={int(spacing)}, SNR={float(snr_db):g} dB, "
                f"trials={int(args.trials_nmse)}",
                flush=True,
            )
            for trial in range(1, int(args.trials_nmse) + 1):
                realization = generate_tdl_channel(rng, grid, channel_cfg, n_tx=2, n_rx=4)
                true_g = cdd_equivalent_from_branches(realization.H, grid, delays)

                noise_alg1 = (
                    rng.normal(size=(4, len(pilot_local)))
                    + 1j * rng.normal(size=(4, len(pilot_local)))
                ) * math.sqrt(noise_var_ls / 2.0)
                ls_cdd = true_g[:, pilot_local] + noise_alg1
                g_hat_alg1 = ls_cdd @ direct.matrix.T

                h_hat = np.empty_like(realization.H)
                for m, est in enumerate(port_total):
                    pilot_k = np.asarray(est.metadata["pilot_k"], dtype=np.int64)
                    pl = local_indices_for_subcarriers(grid, pilot_k)
                    noise = (
                        rng.normal(size=(4, len(pl)))
                        + 1j * rng.normal(size=(4, len(pl)))
                    ) * math.sqrt(noise_var_ls / 2.0)
                    obs = realization.H[:, m, pl] + noise
                    h_hat[:, m, :] = obs @ est.matrix.T
                g_hat_alg3 = cdd_equivalent_from_branches(h_hat, grid, delays)

                acc[ALG1].append(nmse(true_g, g_hat_alg1))
                acc[ALG3].append(nmse(true_g, g_hat_alg3))
                branch_acc.append(nmse(realization.H, h_hat))

                if int(args.progress_every) > 0 and trial % int(args.progress_every) == 0:
                    print(
                        f"[exp19-nmse-progress] spacing={int(spacing)}, SNR={float(snr_db):g}, "
                        f"trial={trial}/{int(args.trials_nmse)}",
                        flush=True,
                    )

            alg1_nmse = mean_or_nan(acc[ALG1])
            alg3_nmse = mean_or_nan(acc[ALG3])
            for alg, value in ((ALG1, alg1_nmse), (ALG3, alg3_nmse)):
                row = dict(base)
                row.update({
                    "phase": "nmse_spacing_sweep",
                    "algorithm": alg,
                    "snr_db": float(snr_db),
                    "trials": int(args.trials_nmse),
                    "ce_nmse_eff_mean": float(value),
                })
                if alg == ALG3:
                    row["ce_nmse_branch_mean"] = mean_or_nan(branch_acc)
                rows.append(row)
            gain_rows.append({
                **base,
                "snr_db": float(snr_db),
                "trials": int(args.trials_nmse),
                "alg1_nmse": float(alg1_nmse),
                "alg3_nmse": float(alg3_nmse),
                "alg3_branch_nmse": mean_or_nan(branch_acc),
                "alg3_nmse_gain_db": float(10.0 * math.log10(float(alg1_nmse) / float(alg3_nmse))),
            })
            write_csv(rows, out_dir / "nmse_spacing_sweep_partial.csv")
            write_csv(gain_rows, out_dir / "nmse_gain_by_spacing_partial.csv")

    write_csv(rows, out_dir / "nmse_spacing_sweep.csv")
    write_csv(gain_rows, out_dir / "nmse_gain_by_spacing.csv")
    return rows, gain_rows


def run_bler_sweep(args: argparse.Namespace, out_dir: Path) -> tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    spacings = parse_int_list(args.dmrs_spacings)
    snrs = parse_float_list(args.snrs)
    mcs = get_mcs("nr_256qam", int(args.mcs_index), None, None)
    channel_cfg = ChannelConfig(
        delay_spread_ns=float(args.delay_spread_ns),
        max_delay_factor=float(args.max_delay_factor),
    )
    rows: List[Dict[str, object]] = []
    target_rows: List[Dict[str, object]] = []

    for spacing in spacings:
        resource = resource_for_spacing(args, int(spacing))
        grid = build_resource_grid(resource)
        tb = build_tb_layout(grid.n_data_re, mcs)
        adapter = SionnaLDPCAdapter(
            tb.cb_k_values,
            tb.cb_e_values,
            num_iter=int(args.ldpc_iterations),
            llr_clip=float(args.llr_clip),
        )
        sample_period_ns = 1e9 / (float(resource.n_fft) * float(resource.scs_khz) * 1e3)
        pdp = make_exponential_pdp(float(args.delay_spread_ns), sample_period_ns, float(args.max_delay_factor))
        data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
        pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
        base = common_metadata(args, grid, resource, int(spacing), sample_period_ns)
        ref_tbs = reference_tbs_bits(grid, mcs)
        base.update({
            "mcs_table": "nr_256qam",
            "mcs_index": int(args.mcs_index),
            "qm": int(mcs.qm),
            "code_rate": float(mcs.code_rate),
            "actual_tbs_bits": int(tb.tb_size),
            "reference_tbs_no_dmrs_bits": int(ref_tbs),
            "coded_bits": int(tb.coded_bits),
            "n_cbs": int(len(tb.cb_k_values)),
        })
        delays = [0, int(args.cdd_delay)]

        for snr_db in snrs:
            noise_var = 10.0 ** (-float(snr_db) / 10.0)
            noise_var_ls = noise_var / float(len(resource.dmrs_symbol_indices))
            direct = make_direct_estimator(
                grid=grid,
                pdp=pdp,
                delays=delays,
                noise_var_ls=noise_var_ls,
                loading=float(args.loading),
            )
            port_total = make_port_estimators(
                grid=grid,
                pdp=pdp,
                n_tx=2,
                noise_var_ls=noise_var_ls,
                loading=float(args.loading),
                mode="equal_total_overhead",
            )
            acc = {
                alg: {"tb_errors": 0, "cb_errors": 0, "ce_nmse": []}
                for alg in (IDEAL, ALG1, ALG3)
            }
            print(
                f"[exp19-bler] spacing={int(spacing)}, SNR={float(snr_db):g} dB, "
                f"trials={int(args.trials_bler)}",
                flush=True,
            )
            for trial in range(1, int(args.trials_bler) + 1):
                rng = np.random.default_rng(
                    stable_seed(
                        "exp19-bler",
                        args.delay_spread_ns,
                        args.cdd_delay,
                        spacing,
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

                h_hat = np.empty_like(realization.H)
                for m, est in enumerate(port_total):
                    pilot_k = np.asarray(est.metadata["pilot_k"], dtype=np.int64)
                    pl = local_indices_for_subcarriers(grid, pilot_k)
                    noise = (
                        rng.normal(size=(4, len(pl)))
                        + 1j * rng.normal(size=(4, len(pl)))
                    ) * math.sqrt(noise_var_ls / 2.0)
                    obs = realization.H[:, m, pl] + noise
                    h_hat[:, m, :] = obs @ est.matrix.T
                g_hat_alg3 = cdd_equivalent_from_branches(h_hat, grid, delays)

                g_hats = [true_g, g_hat_alg1, g_hat_alg3]
                llrs = [
                    llr_for_candidate(true_g, g_hat, data_local, symbols, data_noise, noise_var, int(mcs.qm))
                    for g_hat in g_hats
                ]
                for alg, dec, g_hat in zip((IDEAL, ALG1, ALG3), decode_same_tb_batch(adapter, llrs, payload), g_hats):
                    acc[alg]["tb_errors"] += int(not dec.tb_success)
                    acc[alg]["cb_errors"] += int(sum(1 for ok in dec.cb_success if not ok))
                    acc[alg]["ce_nmse"].append(nmse(true_g, g_hat))

                if int(args.progress_every) > 0 and trial % int(args.progress_every) == 0:
                    status = ", ".join(
                        f"{alg.split()[0]}:{acc[alg]['tb_errors'] / float(trial):.3f}"
                        for alg in (IDEAL, ALG1, ALG3)
                    )
                    print(
                        f"[exp19-bler-progress] spacing={int(spacing)}, SNR={float(snr_db):g}, "
                        f"trial={trial}/{int(args.trials_bler)} | {status}",
                        flush=True,
                    )

            for alg in (IDEAL, ALG1, ALG3):
                tb_errors = int(acc[alg]["tb_errors"])
                bler = float(tb_errors) / float(args.trials_bler)
                row = dict(base)
                row.update({
                    "phase": "bler_spacing_sweep",
                    "algorithm": alg,
                    "snr_db": float(snr_db),
                    "n_trials": int(args.trials_bler),
                    "tb_errors": tb_errors,
                    "cb_errors": int(acc[alg]["cb_errors"]),
                    "bler": float(bler),
                    "cb_bler": float(acc[alg]["cb_errors"]) / float(args.trials_bler * len(tb.cb_k_values)),
                    "bler_resolution": float(1.0 / float(args.trials_bler)),
                    "ce_nmse_eff": mean_or_nan(acc[alg]["ce_nmse"]),
                    "throughput_bits_per_slot": float(ref_tbs) * (1.0 - bler) * float(grid.n_data_re) / float(grid.n_sc * grid.n_symbols),
                    "raw_goodput_bits_per_slot": float(tb.tb_size) * (1.0 - bler),
                })
                rows.append(row)
            write_csv(rows, out_dir / "bler_spacing_sweep_partial.csv")

    for spacing in spacings:
        for alg in (IDEAL, ALG1, ALG3):
            pts = [r for r in rows if int(r["dmrs_spacing_sc"]) == int(spacing) and str(r["algorithm"]) == alg]
            first = pts[0]
            for target in parse_float_list(args.target_blers):
                target_rows.append({
                    "dmrs_spacing_sc": int(spacing),
                    "algorithm": alg,
                    "target_bler": float(target),
                    "snr_at_target_bler_db": interpolate_target_snr(pts, target_bler=float(target)),
                    "delay_spread_ns": first["delay_spread_ns"],
                    "cdd_delay_samples": first["cdd_delay_samples"],
                    "n_trials_per_snr": int(args.trials_bler),
                    "pilot_count_combined": first["pilot_count_combined"],
                    "alg3_pilot_count_per_port": first["alg3_pilot_count_per_port"],
                    "dmrs_overhead": first["dmrs_overhead"],
                })

    write_csv(rows, out_dir / "bler_spacing_sweep.csv")
    write_csv(target_rows, out_dir / "bler_target_snr_by_spacing.csv")
    return rows, target_rows


def make_bler_gain_rows(rows: List[Dict[str, object]], target_rows: List[Dict[str, object]], out_dir: Path) -> List[Dict[str, object]]:
    gain_rows: List[Dict[str, object]] = []
    spacings = sorted({int(r["dmrs_spacing_sc"]) for r in rows})
    snrs = sorted({float(r["snr_db"]) for r in rows})
    for spacing in spacings:
        for snr in snrs:
            alg1 = next(r for r in rows if int(r["dmrs_spacing_sc"]) == spacing and abs(float(r["snr_db"]) - snr) < 1e-9 and r["algorithm"] == ALG1)
            alg3 = next(r for r in rows if int(r["dmrs_spacing_sc"]) == spacing and abs(float(r["snr_db"]) - snr) < 1e-9 and r["algorithm"] == ALG3)
            gain_rows.append({
                "dmrs_spacing_sc": int(spacing),
                "snr_db": float(snr),
                "alg1_bler": float(alg1["bler"]),
                "alg3_bler": float(alg3["bler"]),
                "alg3_bler_reduction": float(alg1["bler"]) - float(alg3["bler"]),
                "alg1_ce_nmse_eff_in_bler_run": float(alg1["ce_nmse_eff"]),
                "alg3_ce_nmse_eff_in_bler_run": float(alg3["ce_nmse_eff"]),
                "alg3_nmse_gain_db_in_bler_run": float(10.0 * math.log10(float(alg1["ce_nmse_eff"]) / float(alg3["ce_nmse_eff"]))),
                "pilot_count_combined": int(alg1["pilot_count_combined"]),
                "alg3_pilot_count_per_port": int(alg1["alg3_pilot_count_per_port"]),
                "dmrs_overhead": float(alg1["dmrs_overhead"]),
            })

    for target in sorted({float(r["target_bler"]) for r in target_rows}):
        for spacing in spacings:
            alg1_t = next(r for r in target_rows if int(r["dmrs_spacing_sc"]) == spacing and r["algorithm"] == ALG1 and abs(float(r["target_bler"]) - target) < 1e-9)
            alg3_t = next(r for r in target_rows if int(r["dmrs_spacing_sc"]) == spacing and r["algorithm"] == ALG3 and abs(float(r["target_bler"]) - target) < 1e-9)
            gain_rows.append({
                "dmrs_spacing_sc": int(spacing),
                "target_bler": float(target),
                "alg1_target_snr_db": alg1_t["snr_at_target_bler_db"],
                "alg3_target_snr_db": alg3_t["snr_at_target_bler_db"],
                "alg3_target_snr_gain_db": (
                    float(alg1_t["snr_at_target_bler_db"]) - float(alg3_t["snr_at_target_bler_db"])
                    if np.isfinite(float(alg1_t["snr_at_target_bler_db"])) and np.isfinite(float(alg3_t["snr_at_target_bler_db"]))
                    else float("nan")
                ),
                "pilot_count_combined": int(alg1_t["pilot_count_combined"]),
                "alg3_pilot_count_per_port": int(alg1_t["alg3_pilot_count_per_port"]),
                "dmrs_overhead": float(alg1_t["dmrs_overhead"]),
            })
    write_csv(gain_rows, out_dir / "bler_gain_by_spacing.csv")
    return gain_rows


def coherence_rows(args: argparse.Namespace, out_dir: Path) -> List[Dict[str, object]]:
    resource = resource_for_spacing(args, parse_int_list(args.dmrs_spacings)[0])
    grid = build_resource_grid(resource)
    sample_period_ns = 1e9 / (float(resource.n_fft) * float(resource.scs_khz) * 1e3)
    pdp = make_exponential_pdp(float(args.delay_spread_ns), sample_period_ns, float(args.max_delay_factor))
    pdp_eff = shifted_pdp(pdp, [0, int(args.cdd_delay)])
    rho_h = freq_corr_from_pdp(pdp, grid.n_fft, grid.n_sc - 1)
    rho_g = freq_corr_from_pdp(pdp_eff, grid.n_fft, grid.n_sc - 1)
    rows = []
    for label, rho in (("physical_branch_H", rho_h), ("cdd_equivalent_g", rho_g)):
        row = {
            "channel": label,
            "delay_spread_ns": float(args.delay_spread_ns),
            "cdd_delay_samples": int(args.cdd_delay),
            "cdd_delay_ns": float(args.cdd_delay) * sample_period_ns,
            "bc_0p9_re": "" if first_below(rho, 0.9) is None else int(first_below(rho, 0.9)),
            "bc_0p7_re": "" if first_below(rho, 0.7) is None else int(first_below(rho, 0.7)),
            "bc_0p5_re": "" if first_below(rho, 0.5) is None else int(first_below(rho, 0.5)),
            "bc_0p2_re": "" if first_below(rho, 0.2) is None else int(first_below(rho, 0.2)),
            "bc_0p1_re": "" if first_below(rho, 0.1) is None else int(first_below(rho, 0.1)),
        }
        rows.append(row)
    write_csv(rows, out_dir / "coherence_bandwidth_cdd128.csv")
    return rows


def plot_nmse_vs_spacing(rows: List[Dict[str, object]], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    snrs = sorted({float(r["snr_db"]) for r in rows})
    selected = []
    for target in (2.0, 5.0, 8.0):
        if target in snrs:
            selected.append(target)
    if not selected:
        selected = snrs[:3]

    fig, axes = plt.subplots(1, len(selected), figsize=(5.2 * len(selected), 4.7), sharey=True)
    if len(selected) == 1:
        axes = [axes]
    for ax, snr in zip(axes, selected):
        for alg in (ALG1, ALG3):
            pts = sorted([r for r in rows if r["algorithm"] == alg and abs(float(r["snr_db"]) - snr) < 1e-9], key=lambda r: int(r["dmrs_spacing_sc"]))
            ax.semilogy([int(r["dmrs_spacing_sc"]) for r in pts], [float(r["ce_nmse_eff_mean"]) for r in pts], marker="o", linewidth=1.35, label=alg)
        ax.set_title(f"SNR={snr:g} dB", fontsize=10)
        ax.set_xlabel("DMRS spacing (sc)")
        ax.grid(True, which="both", alpha=0.3)
    axes[0].set_ylabel("Equivalent-channel NMSE")
    axes[-1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_nmse_gain(gain_rows: List[Dict[str, object]], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    snrs = sorted({float(r["snr_db"]) for r in gain_rows})
    selected = [x for x in (2.0, 5.0, 8.0) if x in snrs] or snrs[:3]
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    for snr in selected:
        pts = sorted([r for r in gain_rows if abs(float(r["snr_db"]) - snr) < 1e-9], key=lambda r: int(r["dmrs_spacing_sc"]))
        ax.plot([int(r["dmrs_spacing_sc"]) for r in pts], [float(r["alg3_nmse_gain_db"]) for r in pts], marker="o", linewidth=1.45, label=f"SNR={snr:g} dB")
    ax.axhline(0.0, color="black", linestyle=":", linewidth=0.9)
    ax.set_xlabel("DMRS spacing (sc)")
    ax.set_ylabel("Alg3 NMSE gain vs Alg1 (dB)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_bler(rows: List[Dict[str, object]], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    spacings = sorted({int(r["dmrs_spacing_sc"]) for r in rows})
    ncols = min(3, len(spacings))
    nrows = int(math.ceil(len(spacings) / float(ncols)))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.0 * ncols, 4.2 * nrows), sharey=True)
    axes = np.asarray(axes).reshape(-1)
    colors = {IDEAL: "#4c78a8", ALG1: "#f58518", ALG3: "#54a24b"}
    for ax, spacing in zip(axes, spacings):
        for alg in (IDEAL, ALG1, ALG3):
            pts = sorted([r for r in rows if int(r["dmrs_spacing_sc"]) == spacing and r["algorithm"] == alg], key=lambda r: float(r["snr_db"]))
            floor = 0.5 / max(float(pts[0]["n_trials"]), 1.0)
            ax.semilogy(
                [float(r["snr_db"]) for r in pts],
                [max(float(r["bler"]), floor) for r in pts],
                marker="o",
                linewidth=1.25,
                color=colors[alg],
                label=alg,
            )
        ax.axhline(0.1, color="black", linestyle=":", linewidth=0.8)
        ax.axhline(0.01, color="black", linestyle="--", linewidth=0.8)
        first = next(r for r in rows if int(r["dmrs_spacing_sc"]) == spacing)
        ax.set_title(
            f"spacing={spacing} sc, combined pilots={int(first['pilot_count_combined'])}",
            fontsize=10,
        )
        ax.set_xlabel("SNR (dB)")
        ax.grid(True, which="both", alpha=0.3)
    for ax in axes[len(spacings):]:
        ax.axis("off")
    axes[0].set_ylabel("TB BLER")
    axes[0].set_ylim(5e-3, 1.0)
    axes[min(len(spacings), len(axes)) - 1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_bler_gain(rows: List[Dict[str, object]], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    snrs = sorted({float(r["snr_db"]) for r in rows})
    selected = [x for x in (5.0, 6.0, 7.0) if x in snrs] or snrs[-3:]
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for snr in selected:
        pts = sorted([r for r in rows if abs(float(r["snr_db"]) - snr) < 1e-9], key=lambda r: int(r["dmrs_spacing_sc"]))
        ax.plot([int(r["dmrs_spacing_sc"]) for r in pts], [float(r["alg3_bler_reduction"]) for r in pts], marker="o", linewidth=1.35, label=f"SNR={snr:g} dB")
    ax.axhline(0.0, color="black", linestyle=":", linewidth=0.9)
    ax.set_xlabel("DMRS spacing (sc)")
    ax.set_ylabel("BLER reduction: Alg1 - Alg3")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out) / f"experiment19_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    coh = coherence_rows(args, out_dir)
    nmse_rows, nmse_gain_rows = run_nmse_sweep(args, out_dir)
    bler_rows, target_rows = run_bler_sweep(args, out_dir)
    bler_gain_rows = make_bler_gain_rows(bler_rows, target_rows, out_dir)

    plot_nmse_vs_spacing(nmse_rows, fig_dir / "experiment19_cdd128_sparse_dmrs_nmse.png")
    plot_nmse_gain(nmse_gain_rows, fig_dir / "experiment19_cdd128_sparse_dmrs_nmse_gain.png")
    plot_bler(bler_rows, fig_dir / "experiment19_cdd128_sparse_dmrs_bler.png")
    plot_bler_gain(
        [r for r in bler_gain_rows if "alg3_bler_reduction" in r],
        fig_dir / "experiment19_cdd128_sparse_dmrs_bler_gain.png",
    )

    with (out_dir / "run_config.json").open("w", encoding="utf-8") as f:
        json.dump({
            **vars(args),
            "coherence_summary": coh,
            "notes": [
                "CDD delay is fixed at the ideal-CSI-favored d=128 point from Experiment 17/18.",
                "DMRS spacing is swept above 24 to stress Alg1 equivalent-channel interpolation.",
                "Alg1 and Alg3 use the same total DMRS overhead at each spacing; Alg3 splits combined pilot positions across two physical ports.",
            ],
        }, f, indent=2, ensure_ascii=False)

    print(out_dir)
    return out_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/experiment19_cdd128_sparse_dmrs")
    parser.add_argument("--fig-dir", default="docs/figures")
    parser.add_argument("--delay-spread-ns", type=float, default=10.0)
    parser.add_argument("--max-delay-factor", type=float, default=8.0)
    parser.add_argument("--cdd-delay", type=int, default=128)
    parser.add_argument("--dmrs-spacings", default="24,36,48,72,96")
    parser.add_argument("--snrs", default="2,3,4,5,6,7,8")
    parser.add_argument("--target-blers", default="0.1,0.01")
    parser.add_argument("--n-prbs", type=int, default=48)
    parser.add_argument("--n-fft", type=int, default=4096)
    parser.add_argument("--pdsch-symbols", type=int, default=10)
    parser.add_argument("--dmrs-symbols", default="2,7")
    parser.add_argument("--mcs-index", type=int, default=8)
    parser.add_argument("--trials-bler", type=int, default=100)
    parser.add_argument("--trials-nmse", type=int, default=300)
    parser.add_argument("--ldpc-iterations", type=int, default=8)
    parser.add_argument("--llr-clip", type=float, default=50.0)
    parser.add_argument("--loading", type=float, default=1e-8)
    parser.add_argument("--seed", type=int, default=20260616)
    parser.add_argument("--progress-every", type=int, default=50)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
