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
from cdd_lls.phy.precoding import cdd_equivalent_from_branches
from cdd_lls.phy.qam import qam_modulate
from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
from cdd_lls.sim.stats import interpolate_target_snr
from tools.run_experiment17_cdd_diversity_bler import (
    decode_same_tb_batch,
    equalize_mrc,
    llr_for_candidate,
    stable_seed,
    write_csv,
)
from tools.run_nmse_algorithm_comparison import (
    apply_basis,
    first_below,
    freq_corr_from_pdp,
    make_basis_estimator,
    make_direct_estimator,
    make_port_estimators,
    mean_or_nan,
)


def parse_float_list(text: str) -> List[float]:
    return [float(x.strip()) for x in str(text).split(",") if x.strip()]


def parse_int_list(text: str) -> List[int]:
    return [int(float(x.strip())) for x in str(text).split(",") if x.strip()]


def resource(args: argparse.Namespace, dmrs_spacing: int) -> ResourceConfig:
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


def target_rows(rows: List[Dict[str, object]], keys: List[str], target: float = 0.1) -> List[Dict[str, object]]:
    groups: Dict[tuple, List[Dict[str, object]]] = {}
    for row in rows:
        groups.setdefault(tuple(row[k] for k in keys), []).append(row)
    out: List[Dict[str, object]] = []
    for key, items in groups.items():
        first = items[0]
        item = {k: v for k, v in zip(keys, key)}
        item.update({
            "target_bler": float(target),
            "snr_at_target_bler_db": interpolate_target_snr(items, target_bler=float(target)),
            "delay_spread_ns": first.get("delay_spread_ns"),
            "dmrs_spacing_sc": first.get("dmrs_spacing_sc"),
            "n_trials": first.get("n_trials"),
        })
        out.append(item)
    return out


def coherence_table(args: argparse.Namespace, cdd_delays: List[int], out_dir: Path) -> List[Dict[str, object]]:
    res = resource(args, int(args.dmrs_spacing_sc))
    grid = build_resource_grid(res)
    sample_period_ns = 1e9 / (float(res.n_fft) * float(res.scs_khz) * 1e3)
    pdp = make_exponential_pdp(float(args.delay_spread_ns), sample_period_ns, float(args.max_delay_factor))
    rows: List[Dict[str, object]] = []
    for d in cdd_delays:
        pdp_eff = shifted_pdp(pdp, [0, int(d)])
        rho = freq_corr_from_pdp(pdp_eff, grid.n_fft, grid.n_sc - 1)
        row = {
            "cdd_delay_samples": int(d),
            "cdd_delay_ns": float(d) * sample_period_ns,
            "delay_spread_ns": float(args.delay_spread_ns),
            "n_fft": int(grid.n_fft),
            "n_sc": int(grid.n_sc),
            "bc_0p9_re": "" if first_below(rho, 0.9) is None else int(first_below(rho, 0.9)),
            "bc_0p7_re": "" if first_below(rho, 0.7) is None else int(first_below(rho, 0.7)),
            "bc_0p5_re": "" if first_below(rho, 0.5) is None else int(first_below(rho, 0.5)),
            "bc_0p2_re": "" if first_below(rho, 0.2) is None else int(first_below(rho, 0.2)),
            "bc_0p1_re": "" if first_below(rho, 0.1) is None else int(first_below(rho, 0.1)),
            "min_abs_corr_within_allocation": float(np.min(rho)),
        }
        rows.append(row)
    write_csv(rows, out_dir / "cdd_coherence_bandwidth.csv")
    return rows


def run_extended_ideal_bler(args: argparse.Namespace, cdd_delays: List[int], out_dir: Path) -> List[Dict[str, object]]:
    snrs = parse_float_list(args.snrs)
    res = resource(args, int(args.dmrs_spacing_sc))
    grid = build_resource_grid(res)
    mcs = get_mcs("nr_256qam", int(args.mcs_index), None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    from cdd_lls.phy.ldpc import SionnaLDPCAdapter

    adapter = SionnaLDPCAdapter(
        tb.cb_k_values,
        tb.cb_e_values,
        num_iter=int(args.ldpc_iterations),
        llr_clip=float(args.llr_clip),
    )
    channel_cfg = ChannelConfig(
        delay_spread_ns=float(args.delay_spread_ns),
        max_delay_factor=float(args.max_delay_factor),
    )
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    rows: List[Dict[str, object]] = []
    for snr_db in snrs:
        noise_var = 10.0 ** (-float(snr_db) / 10.0)
        acc = {int(d): {"tb_errors": 0, "cb_errors": 0} for d in cdd_delays}
        print(f"[ideal-ext] SNR={float(snr_db):g} dB, delays={cdd_delays}", flush=True)
        for trial in range(1, int(args.trials_bler) + 1):
            rng = np.random.default_rng(stable_seed("ideal-ext", args.delay_spread_ns, snr_db, trial, base=int(args.seed)))
            realization = generate_tdl_channel(rng, grid, channel_cfg, n_tx=2, n_rx=4)
            payload = [rng.integers(0, 2, size=int(k), dtype=np.int8) for k in tb.cb_k_values]
            symbols = qam_modulate(np.concatenate(adapter.encode(payload)), int(mcs.qm))
            data_noise = (
                rng.normal(size=(4, len(data_local)))
                + 1j * rng.normal(size=(4, len(data_local)))
            ) * math.sqrt(noise_var / 2.0)
            llrs = []
            for d in cdd_delays:
                true_g = cdd_equivalent_from_branches(realization.H, grid, [0, int(d)])
                llrs.append(
                    llr_for_candidate(true_g, true_g, data_local, symbols, data_noise, noise_var, int(mcs.qm))
                )
            for d, dec in zip(cdd_delays, decode_same_tb_batch(adapter, llrs, payload)):
                acc[int(d)]["tb_errors"] += int(not dec.tb_success)
                acc[int(d)]["cb_errors"] += int(sum(1 for ok in dec.cb_success if not ok))
            if int(args.progress_every) > 0 and trial % int(args.progress_every) == 0:
                status = ", ".join(f"d{d}:{acc[int(d)]['tb_errors'] / float(trial):.3f}" for d in cdd_delays)
                print(f"[ideal-ext-progress] SNR={float(snr_db):g}, trial={trial}/{int(args.trials_bler)} | {status}", flush=True)
        for d in cdd_delays:
            tb_errors = int(acc[int(d)]["tb_errors"])
            cb_errors = int(acc[int(d)]["cb_errors"])
            rows.append({
                "phase": "extended_ideal_cdd_sweep",
                "algorithm": "Ideal CSI",
                "cdd_delay_samples": int(d),
                "delay_spread_ns": float(args.delay_spread_ns),
                "dmrs_spacing_sc": int(args.dmrs_spacing_sc),
                "snr_db": float(snr_db),
                "n_trials": int(args.trials_bler),
                "tb_errors": tb_errors,
                "cb_errors": cb_errors,
                "bler": float(tb_errors) / float(args.trials_bler),
                "cb_bler": float(cb_errors) / float(args.trials_bler * len(tb.cb_k_values)),
            })
        write_csv(rows, out_dir / "extended_ideal_cdd_bler_partial.csv")
    write_csv(rows, out_dir / "extended_ideal_cdd_bler.csv")
    write_csv(target_rows(rows, ["cdd_delay_samples"], 0.1), out_dir / "extended_ideal_cdd_10pct_snr.csv")
    return rows


def run_bler_128_512(args: argparse.Namespace, delays_for_bler: List[int], out_dir: Path) -> List[Dict[str, object]]:
    snrs = parse_float_list(args.snrs)
    res = resource(args, int(args.dmrs_spacing_sc))
    grid = build_resource_grid(res)
    mcs = get_mcs("nr_256qam", int(args.mcs_index), None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    from cdd_lls.phy.ldpc import SionnaLDPCAdapter

    adapter = SionnaLDPCAdapter(
        tb.cb_k_values,
        tb.cb_e_values,
        num_iter=int(args.ldpc_iterations),
        llr_clip=float(args.llr_clip),
    )
    channel_cfg = ChannelConfig(
        delay_spread_ns=float(args.delay_spread_ns),
        max_delay_factor=float(args.max_delay_factor),
    )
    sample_period_ns = 1e9 / (float(res.n_fft) * float(res.scs_khz) * 1e3)
    pdp = make_exponential_pdp(float(args.delay_spread_ns), sample_period_ns, float(args.max_delay_factor))
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    rows: List[Dict[str, object]] = []
    for snr_db in snrs:
        noise_var = 10.0 ** (-float(snr_db) / 10.0)
        noise_var_ls = noise_var / float(len(res.dmrs_symbol_indices))
        direct = {
            int(d): make_direct_estimator(grid, pdp, [0, int(d)], noise_var_ls, float(args.loading))
            for d in delays_for_bler
        }
        port_total = make_port_estimators(
            grid, pdp, n_tx=2, noise_var_ls=noise_var_ls, loading=float(args.loading), mode="equal_total_overhead"
        )
        algs = ["Ideal CSI", "Alg1 Direct RMMSE WB", "Alg3 Port-DMRS RMMSE equal-total"]
        acc = {
            (int(d), alg): {"tb_errors": 0, "cb_errors": 0, "ce_nmse": []}
            for d in delays_for_bler
            for alg in algs
        }
        print(f"[bler-compare] SNR={float(snr_db):g} dB, delays={delays_for_bler}", flush=True)
        for trial in range(1, int(args.trials_bler) + 1):
            rng = np.random.default_rng(stable_seed("bler-128-512", args.delay_spread_ns, snr_db, trial, base=int(args.seed)))
            realization = generate_tdl_channel(rng, grid, channel_cfg, n_tx=2, n_rx=4)
            payload = [rng.integers(0, 2, size=int(k), dtype=np.int8) for k in tb.cb_k_values]
            symbols = qam_modulate(np.concatenate(adapter.encode(payload)), int(mcs.qm))
            data_noise = (
                rng.normal(size=(4, len(data_local)))
                + 1j * rng.normal(size=(4, len(data_local)))
            ) * math.sqrt(noise_var / 2.0)
            pilot_noise = (
                rng.normal(size=(4, len(pilot_local)))
                + 1j * rng.normal(size=(4, len(pilot_local)))
            ) * math.sqrt(noise_var_ls / 2.0)
            llrs = []
            keys = []
            g_refs: Dict[int, np.ndarray] = {}
            for d in delays_for_bler:
                true_g = cdd_equivalent_from_branches(realization.H, grid, [0, int(d)])
                g_refs[int(d)] = true_g
                ls_cdd = true_g[:, pilot_local] + pilot_noise
                g_hat_alg1 = ls_cdd @ direct[int(d)].matrix.T
                llrs.append(llr_for_candidate(true_g, true_g, data_local, symbols, data_noise, noise_var, int(mcs.qm)))
                keys.append((int(d), "Ideal CSI", true_g))
                llrs.append(llr_for_candidate(true_g, g_hat_alg1, data_local, symbols, data_noise, noise_var, int(mcs.qm)))
                keys.append((int(d), "Alg1 Direct RMMSE WB", g_hat_alg1))
            # Estimate the physical branches once, then recombine with each CDD delay.
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
            for d in delays_for_bler:
                g_hat_alg3 = cdd_equivalent_from_branches(h_hat, grid, [0, int(d)])
                true_g = g_refs[int(d)]
                llrs.append(llr_for_candidate(true_g, g_hat_alg3, data_local, symbols, data_noise, noise_var, int(mcs.qm)))
                keys.append((int(d), "Alg3 Port-DMRS RMMSE equal-total", g_hat_alg3))
            for (d, alg, g_hat), dec in zip(keys, decode_same_tb_batch(adapter, llrs, payload)):
                acc[(int(d), alg)]["tb_errors"] += int(not dec.tb_success)
                acc[(int(d), alg)]["cb_errors"] += int(sum(1 for ok in dec.cb_success if not ok))
                acc[(int(d), alg)]["ce_nmse"].append(nmse(g_refs[int(d)], g_hat))
            if int(args.progress_every) > 0 and trial % int(args.progress_every) == 0:
                status = ", ".join(
                    f"d{d}-{alg.split()[0]}:{acc[(int(d), alg)]['tb_errors'] / float(trial):.3f}"
                    for d in delays_for_bler
                    for alg in algs
                )
                print(f"[bler-compare-progress] SNR={float(snr_db):g}, trial={trial}/{int(args.trials_bler)} | {status}", flush=True)
        for d in delays_for_bler:
            for alg in algs:
                item = acc[(int(d), alg)]
                rows.append({
                    "phase": "cdd128_vs_cdd512_bler",
                    "algorithm": alg,
                    "cdd_delay_samples": int(d),
                    "delay_spread_ns": float(args.delay_spread_ns),
                    "dmrs_spacing_sc": int(args.dmrs_spacing_sc),
                    "snr_db": float(snr_db),
                    "n_trials": int(args.trials_bler),
                    "tb_errors": int(item["tb_errors"]),
                    "cb_errors": int(item["cb_errors"]),
                    "bler": float(item["tb_errors"]) / float(args.trials_bler),
                    "cb_bler": float(item["cb_errors"]) / float(args.trials_bler * len(tb.cb_k_values)),
                    "ce_nmse_eff": mean_or_nan(item["ce_nmse"]),
                    "pilot_count_combined": int(grid.pilot_count),
                    "alg3_pilot_count_per_port": int(len(grid.pilot_subcarriers[0::2])),
                })
        write_csv(rows, out_dir / "cdd128_cdd512_bler_partial.csv")
    write_csv(rows, out_dir / "cdd128_cdd512_bler.csv")
    return rows


def run_nmse_alg1_alg2(args: argparse.Namespace, cdd_delays: List[int], out_dir: Path) -> List[Dict[str, object]]:
    snrs = parse_float_list(args.snrs)
    res = resource(args, int(args.dmrs_spacing_sc))
    grid = build_resource_grid(res)
    sample_period_ns = 1e9 / (float(res.n_fft) * float(res.scs_khz) * 1e3)
    pdp = make_exponential_pdp(float(args.delay_spread_ns), sample_period_ns, float(args.max_delay_factor))
    channel_cfg = ChannelConfig(delay_spread_ns=float(args.delay_spread_ns), max_delay_factor=float(args.max_delay_factor))
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    rows: List[Dict[str, object]] = []
    for snr_db in snrs:
        noise_var = 10.0 ** (-float(snr_db) / 10.0)
        noise_var_ls = noise_var / float(len(res.dmrs_symbol_indices))
        estimators = {}
        for d in cdd_delays:
            delays = [0, int(d)]
            estimators[(int(d), "Alg1 Direct RMMSE WB")] = make_direct_estimator(
                grid, pdp, delays, noise_var_ls, float(args.loading)
            )
            estimators[(int(d), "Alg2 Basis LMMSE E99")] = make_basis_estimator(
                grid, pdp, delays, noise_var_ls, float(args.loading), 0.99
            )
        acc = {(int(d), alg): [] for d in cdd_delays for alg in ("Alg1 Direct RMMSE WB", "Alg2 Basis LMMSE E99")}
        branch_acc = {(int(d), "Alg2 Basis LMMSE E99"): [] for d in cdd_delays}
        print(f"[nmse] SNR={float(snr_db):g} dB, delays={cdd_delays}", flush=True)
        rng = np.random.default_rng(stable_seed("nmse", args.delay_spread_ns, snr_db, base=int(args.seed)))
        for _ in range(int(args.trials_nmse)):
            realization = generate_tdl_channel(rng, grid, channel_cfg, n_tx=2, n_rx=4)
            noise = (
                rng.normal(size=(4, len(pilot_local)))
                + 1j * rng.normal(size=(4, len(pilot_local)))
            ) * math.sqrt(noise_var_ls / 2.0)
            for d in cdd_delays:
                delays = [0, int(d)]
                true_g = cdd_equivalent_from_branches(realization.H, grid, delays)
                ls_cdd = true_g[:, pilot_local] + noise
                direct = estimators[(int(d), "Alg1 Direct RMMSE WB")]
                g_hat = ls_cdd @ direct.matrix.T
                acc[(int(d), "Alg1 Direct RMMSE WB")].append(nmse(true_g, g_hat))
                basis = estimators[(int(d), "Alg2 Basis LMMSE E99")]
                g_hat_basis, h_hat = apply_basis(basis, ls_cdd, grid, delays)
                acc[(int(d), "Alg2 Basis LMMSE E99")].append(nmse(true_g, g_hat_basis))
                branch_acc[(int(d), "Alg2 Basis LMMSE E99")].append(nmse(realization.H, h_hat))
        for d in cdd_delays:
            for alg in ("Alg1 Direct RMMSE WB", "Alg2 Basis LMMSE E99"):
                row = {
                    "algorithm": alg,
                    "cdd_delay_samples": int(d),
                    "cdd_delay_ns": float(d) * sample_period_ns,
                    "delay_spread_ns": float(args.delay_spread_ns),
                    "dmrs_spacing_sc": int(args.dmrs_spacing_sc),
                    "snr_db": float(snr_db),
                    "trials": int(args.trials_nmse),
                    "ce_nmse_eff_mean": mean_or_nan(acc[(int(d), alg)]),
                    "pilot_count_combined": int(grid.pilot_count),
                    "processing_bandwidth_sc": int(grid.n_sc),
                }
                if alg == "Alg2 Basis LMMSE E99":
                    est = estimators[(int(d), alg)]
                    row.update({
                        "ce_nmse_branch_mean": mean_or_nan(branch_acc[(int(d), alg)]),
                        "support_len": int(est.metadata["support_len"]),
                        "n_unknown_taps": int(est.metadata["n_unknown_taps"]),
                        "cond_number": float(est.cond_number),
                        "effective_rank": float(est.effective_rank),
                    })
                rows.append(row)
        write_csv(rows, out_dir / "alg1_alg2_nmse_by_cdd_partial.csv")
    write_csv(rows, out_dir / "alg1_alg2_nmse_by_cdd.csv")
    return rows


def plot_extended_ideal(rows: List[Dict[str, object]], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    delays = sorted({int(r["cdd_delay_samples"]) for r in rows})
    fig, ax = plt.subplots(figsize=(9.5, 5.4))
    for d in delays:
        pts = sorted([r for r in rows if int(r["cdd_delay_samples"]) == d], key=lambda r: float(r["snr_db"]))
        if not pts:
            continue
        floor = 0.5 / max(float(pts[0]["n_trials"]), 1.0)
        ax.semilogy([float(r["snr_db"]) for r in pts], [max(float(r["bler"]), floor) for r in pts], marker="o", linewidth=1.2, label=f"d={d}")
    ax.axhline(0.1, color="black", linestyle=":", linewidth=0.9)
    ax.set_xlabel("SNR (dB)")
    ax.set_ylabel("TB BLER")
    ax.set_ylim(5e-3, 1.0)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7, ncol=3)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_bler_128_512(rows: List[Dict[str, object]], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    styles = {
        "Ideal CSI": ("#4c78a8", "-"),
        "Alg1 Direct RMMSE WB": ("#f58518", "-"),
        "Alg3 Port-DMRS RMMSE equal-total": ("#54a24b", "-"),
    }
    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    for d in sorted({int(r["cdd_delay_samples"]) for r in rows}):
        for alg, (color, line) in styles.items():
            pts = sorted([r for r in rows if int(r["cdd_delay_samples"]) == d and str(r["algorithm"]) == alg], key=lambda r: float(r["snr_db"]))
            if not pts:
                continue
            floor = 0.5 / max(float(pts[0]["n_trials"]), 1.0)
            marker = "o" if d == 128 else "s"
            linestyle = line if d == 128 else "--"
            ax.semilogy([float(r["snr_db"]) for r in pts], [max(float(r["bler"]), floor) for r in pts], marker=marker, linestyle=linestyle, color=color, linewidth=1.25, label=f"{alg}, d={d}")
    ax.axhline(0.1, color="black", linestyle=":", linewidth=0.9)
    ax.axhline(0.01, color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel("SNR (dB)")
    ax.set_ylabel("TB BLER")
    ax.set_ylim(5e-3, 1.0)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_nmse(rows: List[Dict[str, object]], path: Path, delays: List[int]) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.9), sharey=True)
    for ax, alg in zip(axes, ["Alg1 Direct RMMSE WB", "Alg2 Basis LMMSE E99"]):
        for d in delays:
            pts = sorted([r for r in rows if int(r["cdd_delay_samples"]) == int(d) and str(r["algorithm"]) == alg], key=lambda r: float(r["snr_db"]))
            if not pts:
                continue
            ax.semilogy([float(r["snr_db"]) for r in pts], [float(r["ce_nmse_eff_mean"]) for r in pts], marker="o", linewidth=1.25, label=f"d={d}")
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
    if 5.0 in available_snrs and 5.0 not in selected_snrs:
        selected_snrs.append(5.0)
    elif available_snrs and available_snrs[-1] not in selected_snrs:
        selected_snrs.append(available_snrs[-1])
    while len(selected_snrs) < 2 and available_snrs:
        selected_snrs.append(available_snrs[-1])

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.9), sharey=True)
    for ax, snr in zip(axes, selected_snrs[:2]):
        for alg in ["Alg1 Direct RMMSE WB", "Alg2 Basis LMMSE E99"]:
            pts = sorted([r for r in rows if str(r["algorithm"]) == alg and abs(float(r["snr_db"]) - snr) < 1e-9], key=lambda r: int(r["cdd_delay_samples"]))
            if not pts:
                continue
            ax.semilogy([int(r["cdd_delay_samples"]) for r in pts], [float(r["ce_nmse_eff_mean"]) for r in pts], marker="o", linewidth=1.25, label=alg)
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
    out_dir = Path(args.out) / f"experiment18_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    extended_delays = parse_int_list(args.extended_cdd_delays)
    bler_compare_delays = parse_int_list(args.bler_compare_cdd_delays)
    nmse_delays = parse_int_list(args.nmse_cdd_delays)

    coh_rows = coherence_table(args, extended_delays, out_dir)
    ideal_rows = run_extended_ideal_bler(args, extended_delays, out_dir)
    compare_rows = run_bler_128_512(args, bler_compare_delays, out_dir)
    nmse_rows = run_nmse_alg1_alg2(args, nmse_delays, out_dir)

    plot_extended_ideal(ideal_rows, fig_dir / "experiment18_extended_ideal_cdd_bler.png")
    plot_bler_128_512(compare_rows, fig_dir / "experiment18_cdd128_cdd512_bler.png")
    plot_nmse(nmse_rows, fig_dir / "experiment18_alg1_alg2_nmse_d128_d512.png", bler_compare_delays)
    plot_nmse_vs_cdd(nmse_rows, fig_dir / "experiment18_alg1_alg2_nmse_vs_cdd.png")

    with (out_dir / "run_config.json").open("w", encoding="utf-8") as f:
        json.dump({
            **vars(args),
            "notes": [
                "Alg1 Direct RMMSE uses shifted effective CDD PDP and full 576-sc target bandwidth.",
                "Alg2 is Basis LMMSE E99 in this follow-up.",
                "Alg3 Port-DMRS uses the same lmmse_matrix() per port and full 576-sc target bandwidth; equal-total mode splits the combined pilots across ports.",
            ],
        }, f, indent=2, ensure_ascii=False)

    print(out_dir)
    return out_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/experiment18_cdd_delay_followup")
    parser.add_argument("--fig-dir", default="docs/figures")
    parser.add_argument("--delay-spread-ns", type=float, default=10.0)
    parser.add_argument("--max-delay-factor", type=float, default=8.0)
    parser.add_argument("--extended-cdd-delays", default="0,64,128,256,512,768,1024,1536,2048")
    parser.add_argument("--bler-compare-cdd-delays", default="128,512")
    parser.add_argument("--nmse-cdd-delays", default="0,64,128,256,512,768,1024,1536,2048")
    parser.add_argument("--dmrs-spacing-sc", type=int, default=24)
    parser.add_argument("--snrs", default="2,3,4,5,6")
    parser.add_argument("--n-prbs", type=int, default=48)
    parser.add_argument("--n-fft", type=int, default=4096)
    parser.add_argument("--pdsch-symbols", type=int, default=10)
    parser.add_argument("--dmrs-symbols", default="2,7")
    parser.add_argument("--mcs-index", type=int, default=8)
    parser.add_argument("--trials-bler", type=int, default=100)
    parser.add_argument("--trials-nmse", type=int, default=200)
    parser.add_argument("--ldpc-iterations", type=int, default=8)
    parser.add_argument("--llr-clip", type=float, default=50.0)
    parser.add_argument("--loading", type=float, default=1e-8)
    parser.add_argument("--seed", type=int, default=20260616)
    parser.add_argument("--progress-every", type=int, default=25)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
