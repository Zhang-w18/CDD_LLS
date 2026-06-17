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
from cdd_lls.phy.estimators import nmse
from cdd_lls.phy.ldpc import SionnaLDPCAdapter, summarize_decode
from cdd_lls.phy.precoding import cdd_equivalent_from_branches
from cdd_lls.phy.qam import qam_demapper_maxlog, qam_modulate
from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
from cdd_lls.sim.stats import interpolate_target_snr
from tools.run_nmse_algorithm_comparison import (
    apply_port_dmrs,
    make_direct_estimator,
    make_port_estimators,
)


def parse_float_list(text: str) -> List[float]:
    return [float(x.strip()) for x in str(text).split(",") if x.strip()]


def parse_int_list(text: str) -> List[int]:
    return [int(float(x.strip())) for x in str(text).split(",") if x.strip()]


def stable_seed(*values: object, base: int = 0) -> int:
    h = int(base) & 0xFFFFFFFF
    for value in values:
        for ch in str(value):
            h ^= ord(ch)
            h = (h * 16777619) & 0xFFFFFFFF
    return int(h)


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


def split_llrs(llr: np.ndarray, cb_e_values: List[int]) -> List[np.ndarray]:
    arr = np.asarray(llr, dtype=np.float64).reshape(-1)
    out: List[np.ndarray] = []
    cursor = 0
    for e in cb_e_values:
        e = int(e)
        out.append(arr[cursor:cursor + e].copy())
        cursor += e
    if cursor != arr.size:
        raise ValueError(f"CB E values sum to {cursor}, LLR stream has {arr.size}.")
    return out


def decode_same_tb_batch(
    adapter: SionnaLDPCAdapter,
    llrs_by_candidate: List[np.ndarray],
    reference_payload_bits_by_cb: List[np.ndarray],
):
    if not llrs_by_candidate:
        return []
    split_by_candidate = [
        split_llrs(llr, adapter.cb_e_values)
        for llr in llrs_by_candidate
    ]
    decoded_by_candidate: List[List[np.ndarray]] = [
        [] for _ in llrs_by_candidate
    ]
    for cb_idx, _ in enumerate(adapter.cb_e_values):
        batch = np.stack(
            [parts[cb_idx] for parts in split_by_candidate],
            axis=0,
        ).astype(np.float32)
        if adapter.llr_clip > 0:
            batch = np.clip(batch, -adapter.llr_clip, adapter.llr_clip)
        b_hat = adapter.decoder_by_cb[int(cb_idx)](
            adapter.tf.constant(batch, dtype=adapter.tf.float32)
        )
        decoded_batch = np.rint(b_hat.numpy()).astype(np.int8)
        for cand_idx in range(len(llrs_by_candidate)):
            decoded_by_candidate[cand_idx].append(decoded_batch[cand_idx])
    return [
        summarize_decode(decoded, reference_payload_bits_by_cb)
        for decoded in decoded_by_candidate
    ]


def equalize_mrc(y: np.ndarray, g_hat: np.ndarray, noise_var: float) -> tuple[np.ndarray, np.ndarray]:
    gh = np.asarray(g_hat, dtype=np.complex128)
    yy = np.asarray(y, dtype=np.complex128)
    denom = np.sum(np.abs(gh) ** 2, axis=0)
    denom = np.maximum(denom, 1e-10)
    z = np.sum(np.conj(gh) * yy, axis=0) / denom
    no_eff = float(noise_var) / denom
    return z, no_eff


def resource_for_spacing(args: argparse.Namespace, dmrs_spacing_sc: int) -> ResourceConfig:
    return ResourceConfig(
        carrier_bandwidth_mhz=100.0,
        scs_khz=30,
        n_fft=int(args.n_fft),
        n_prbs=int(args.n_prbs),
        pdsch_n_symbols=int(args.pdsch_symbols),
        dmrs_symbol_indices=parse_int_list(args.dmrs_symbols),
        dmrs_spacing_sc=int(dmrs_spacing_sc),
        prg_size_rb=4,
    )


def throughput_factor(grid) -> float:
    total_re = int(grid.n_sc) * int(grid.n_symbols)
    return float(grid.n_data_re) / float(total_re)


def reference_tbs_bits(grid, mcs) -> int:
    total_re = int(grid.n_sc) * int(grid.n_symbols)
    return int(build_tb_layout(total_re, mcs).tb_size)


def make_base_metadata(args, resource, grid, mcs, tb, ref_tbs: int) -> Dict[str, object]:
    total_re = int(grid.n_sc) * int(grid.n_symbols)
    factor = throughput_factor(grid)
    return {
        "delay_spread_ns": float(args.delay_spread_ns),
        "n_tx": 2,
        "n_rx": 4,
        "n_prbs": int(args.n_prbs),
        "n_sc": int(grid.n_sc),
        "n_fft": int(args.n_fft),
        "pdsch_symbols": int(args.pdsch_symbols),
        "dmrs_symbols": int(len(resource.dmrs_symbol_indices)),
        "dmrs_spacing_sc": int(resource.dmrs_spacing_sc),
        "pilot_count_combined": int(grid.pilot_count),
        "total_dmrs_re": int(grid.n_dmrs_re),
        "total_re": int(total_re),
        "data_re": int(grid.n_data_re),
        "dmrs_overhead": float(grid.dmrs_overhead),
        "throughput_factor": float(factor),
        "reference_tbs_no_dmrs_bits": int(ref_tbs),
        "actual_tbs_bits": int(tb.tb_size),
        "coded_bits": int(tb.coded_bits),
        "n_cbs": int(len(tb.cb_k_values)),
        "mcs_table": "nr_256qam",
        "mcs_index": int(args.mcs_index),
        "qm": int(mcs.qm),
        "code_rate": float(mcs.code_rate),
        "ldpc_iterations": int(args.ldpc_iterations),
    }


def finalize_row(row: Dict[str, object], tb_errors: int, cb_errors: int, n_trials: int, n_cbs: int) -> Dict[str, object]:
    out = dict(row)
    bler = float(tb_errors) / float(n_trials)
    cb_bler = float(cb_errors) / float(max(n_trials * n_cbs, 1))
    ref_tbs = float(out["reference_tbs_no_dmrs_bits"])
    factor = float(out["throughput_factor"])
    actual_tbs = float(out["actual_tbs_bits"])
    out.update({
        "n_trials": int(n_trials),
        "tb_errors": int(tb_errors),
        "cb_errors": int(cb_errors),
        "bler": float(bler),
        "cb_bler": float(cb_bler),
        "bler_resolution": float(1.0 / float(n_trials)),
        "throughput_bits_per_slot": float(ref_tbs * (1.0 - bler) * factor),
        "raw_goodput_bits_per_slot": float(actual_tbs * (1.0 - bler)),
    })
    return out


def cdd_label(delay: int) -> str:
    return f"CDD d={int(delay)}"


def apply_data_channel(true_g: np.ndarray, data_local: np.ndarray, symbols: np.ndarray, noise: np.ndarray) -> np.ndarray:
    return true_g[:, data_local] * symbols[None, :] + noise


def llr_for_candidate(
    true_g: np.ndarray,
    g_hat: np.ndarray,
    data_local: np.ndarray,
    symbols: np.ndarray,
    data_noise: np.ndarray,
    noise_var: float,
    qm: int,
) -> np.ndarray:
    y = apply_data_channel(true_g, data_local, symbols, data_noise)
    z, no_eff = equalize_mrc(y, g_hat[:, data_local], noise_var)
    return qam_demapper_maxlog(z, no_eff, qm)


def run_ideal_cdd_scan(args: argparse.Namespace, out_dir: Path) -> tuple[List[Dict[str, object]], List[Dict[str, object]], int]:
    cdd_delays = parse_int_list(args.ideal_cdd_delays)
    snrs = parse_float_list(args.ideal_snrs)
    resource = resource_for_spacing(args, int(args.ideal_dmrs_spacing_sc))
    grid = build_resource_grid(resource)
    mcs = get_mcs("nr_256qam", int(args.mcs_index), None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    ref_tbs = reference_tbs_bits(grid, mcs)
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
    accum_debug: List[Dict[str, object]] = []

    for snr_db in snrs:
        noise_var = 10.0 ** (-float(snr_db) / 10.0)
        accum = {
            int(delay): {"tb_errors": 0, "cb_errors": 0}
            for delay in cdd_delays
        }
        print(
            f"[ideal-start] DS={float(args.delay_spread_ns):g} ns, "
            f"SNR={float(snr_db):g} dB, trials={int(args.trials)}, "
            f"CDD={cdd_delays}",
            flush=True,
        )
        for trial in range(1, int(args.trials) + 1):
            rng = np.random.default_rng(stable_seed("ideal", args.delay_spread_ns, snr_db, trial, base=int(args.seed)))
            realization = generate_tdl_channel(rng, grid, channel_cfg, n_tx=2, n_rx=4)
            payload = [
                rng.integers(0, 2, size=int(k), dtype=np.int8)
                for k in tb.cb_k_values
            ]
            coded = adapter.encode(payload)
            symbols = qam_modulate(np.concatenate(coded), int(mcs.qm))
            data_noise = (
                rng.normal(size=(4, len(data_local)))
                + 1j * rng.normal(size=(4, len(data_local)))
            ) * math.sqrt(noise_var / 2.0)
            llrs = []
            trial_delays = []
            for delay in cdd_delays:
                true_g = cdd_equivalent_from_branches(realization.H, grid, [0, int(delay)])
                llrs.append(
                    llr_for_candidate(
                        true_g=true_g,
                        g_hat=true_g,
                        data_local=data_local,
                        symbols=symbols,
                        data_noise=data_noise,
                        noise_var=noise_var,
                        qm=int(mcs.qm),
                    )
                )
                trial_delays.append(int(delay))
            dec_results = decode_same_tb_batch(adapter, llrs, payload)
            for delay, dec in zip(trial_delays, dec_results):
                accum[int(delay)]["tb_errors"] += int(not dec.tb_success)
                accum[int(delay)]["cb_errors"] += int(sum(1 for ok in dec.cb_success if not ok))
            if int(args.progress_every) > 0 and trial % int(args.progress_every) == 0:
                status = ", ".join(
                    f"d{delay}:BLER={accum[int(delay)]['tb_errors'] / float(trial):.3f}"
                    for delay in cdd_delays
                )
                print(
                    f"[ideal-progress] SNR={float(snr_db):g} dB, "
                    f"trial={trial}/{int(args.trials)} | {status}",
                    flush=True,
                )
        base = make_base_metadata(args, resource, grid, mcs, tb, ref_tbs)
        for delay in cdd_delays:
            row = dict(base)
            row.update({
                "phase": "ideal_csi_cdd_scan",
                "algorithm": "Ideal CSI",
                "cdd_delay_samples": int(delay),
                "cdd_delay_ns": float(delay) * float(realization.sample_period_ns),
                "precoder": cdd_label(delay),
                "snr_db": float(snr_db),
                "ce_nmse_eff": 0.0,
            })
            rows.append(
                finalize_row(
                    row,
                    tb_errors=int(accum[int(delay)]["tb_errors"]),
                    cb_errors=int(accum[int(delay)]["cb_errors"]),
                    n_trials=int(args.trials),
                    n_cbs=int(len(tb.cb_k_values)),
                )
            )
        write_csv(rows, out_dir / "ideal_csi_cdd_scan_partial.csv")
        accum_debug.append({"snr_db": float(snr_db), "accum": accum})

    target_rows = target_summary(
        rows,
        group_keys=["cdd_delay_samples", "precoder"],
        targets=parse_float_list(args.target_blers),
    )
    selected = select_cdd_delay(rows, target_rows, cdd_delays, parse_float_list(args.selection_targets))
    return rows, target_rows, selected


def select_cdd_delay(
    rows: List[Dict[str, object]],
    target_rows: List[Dict[str, object]],
    cdd_delays: List[int],
    selection_targets: List[float],
) -> int:
    for target in selection_targets:
        candidates = [
            r for r in target_rows
            if abs(float(r["target_bler"]) - float(target)) < 1e-12
            and np.isfinite(float(r["snr_at_target_bler_db"]))
        ]
        if candidates:
            best = min(candidates, key=lambda r: float(r["snr_at_target_bler_db"]))
            return int(best["cdd_delay_samples"])

    scores = []
    for delay in cdd_delays:
        pts = [r for r in rows if int(r["cdd_delay_samples"]) == int(delay)]
        if not pts:
            continue
        mean_bler = float(np.mean([float(r["bler"]) for r in pts]))
        min_bler = float(np.min([float(r["bler"]) for r in pts]))
        scores.append((mean_bler, min_bler, int(delay)))
    if not scores:
        return int(cdd_delays[0])
    scores.sort()
    return int(scores[0][2])


def run_algorithm_comparison(args: argparse.Namespace, out_dir: Path, selected_cdd_delay: int) -> tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    dmrs_spacings = parse_int_list(args.compare_dmrs_spacings)
    snrs = parse_float_list(args.compare_snrs)
    mcs = get_mcs("nr_256qam", int(args.mcs_index), None, None)
    rows: List[Dict[str, object]] = []

    for dmrs_spacing in dmrs_spacings:
        resource = resource_for_spacing(args, int(dmrs_spacing))
        grid = build_resource_grid(resource)
        tb = build_tb_layout(grid.n_data_re, mcs)
        ref_tbs = reference_tbs_bits(grid, mcs)
        adapter = SionnaLDPCAdapter(
            tb.cb_k_values,
            tb.cb_e_values,
            num_iter=int(args.ldpc_iterations),
            llr_clip=float(args.llr_clip),
        )
        sample_period_ns = 1e9 / (float(resource.n_fft) * float(resource.scs_khz) * 1e3)
        pdp = make_exponential_pdp(float(args.delay_spread_ns), sample_period_ns, float(args.max_delay_factor))
        channel_cfg = ChannelConfig(
            delay_spread_ns=float(args.delay_spread_ns),
            max_delay_factor=float(args.max_delay_factor),
        )
        data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
        pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
        delays = [0, int(selected_cdd_delay)]

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
            algorithms = [
                "Ideal CSI",
                "Alg1 Direct RMMSE WB",
                "Alg3 Port-DMRS RMMSE equal-total",
            ]
            accum: Dict[str, Dict[str, object]] = {
                alg: {"tb_errors": 0, "cb_errors": 0, "ce_nmse": []}
                for alg in algorithms
            }
            print(
                f"[compare-start] CDD={int(selected_cdd_delay)}, "
                f"DMRS spacing={int(dmrs_spacing)}, SNR={float(snr_db):g} dB, "
                f"trials={int(args.trials)}",
                flush=True,
            )
            for trial in range(1, int(args.trials) + 1):
                rng = np.random.default_rng(
                    stable_seed("compare", args.delay_spread_ns, selected_cdd_delay, dmrs_spacing, snr_db, trial, base=int(args.seed))
                )
                realization = generate_tdl_channel(rng, grid, channel_cfg, n_tx=2, n_rx=4)
                true_g = cdd_equivalent_from_branches(realization.H, grid, delays)
                payload = [
                    rng.integers(0, 2, size=int(k), dtype=np.int8)
                    for k in tb.cb_k_values
                ]
                coded = adapter.encode(payload)
                symbols = qam_modulate(np.concatenate(coded), int(mcs.qm))
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

                g_hats = [
                    true_g,
                    g_hat_alg1,
                    g_hat_alg3,
                ]
                llrs = [
                    llr_for_candidate(
                        true_g=true_g,
                        g_hat=g_hat,
                        data_local=data_local,
                        symbols=symbols,
                        data_noise=data_noise,
                        noise_var=noise_var,
                        qm=int(mcs.qm),
                    )
                    for g_hat in g_hats
                ]
                dec_results = decode_same_tb_batch(adapter, llrs, payload)
                for alg, dec, g_hat in zip(algorithms, dec_results, g_hats):
                    accum[alg]["tb_errors"] += int(not dec.tb_success)
                    accum[alg]["cb_errors"] += int(sum(1 for ok in dec.cb_success if not ok))
                    accum[alg]["ce_nmse"].append(nmse(true_g, g_hat))
                if int(args.progress_every) > 0 and trial % int(args.progress_every) == 0:
                    status = ", ".join(
                        f"{alg}:BLER={int(accum[alg]['tb_errors']) / float(trial):.3f}"
                        for alg in algorithms
                    )
                    print(
                        f"[compare-progress] sf={int(dmrs_spacing)}, "
                        f"SNR={float(snr_db):g} dB, trial={trial}/{int(args.trials)} | {status}",
                        flush=True,
                    )
            base = make_base_metadata(args, resource, grid, mcs, tb, ref_tbs)
            for alg in algorithms:
                row = dict(base)
                row.update({
                    "phase": "algorithm_comparison",
                    "algorithm": alg,
                    "cdd_delay_samples": int(selected_cdd_delay),
                    "cdd_delay_ns": float(selected_cdd_delay) * float(sample_period_ns),
                    "precoder": cdd_label(selected_cdd_delay),
                    "snr_db": float(snr_db),
                    "ce_nmse_eff": mean_or_nan(accum[alg]["ce_nmse"]),
                    "alg3_pilot_count_per_port": int(len(grid.pilot_subcarriers[0::2])),
                    "same_total_dmrs_overhead_alg1_alg3": True,
                })
                rows.append(
                    finalize_row(
                        row,
                        tb_errors=int(accum[alg]["tb_errors"]),
                        cb_errors=int(accum[alg]["cb_errors"]),
                        n_trials=int(args.trials),
                        n_cbs=int(len(tb.cb_k_values)),
                    )
                )
            write_csv(rows, out_dir / "algorithm_comparison_partial.csv")

    target_rows = target_summary(
        rows,
        group_keys=["dmrs_spacing_sc", "algorithm"],
        targets=parse_float_list(args.target_blers),
    )
    return rows, target_rows


def target_summary(
    rows: List[Dict[str, object]],
    group_keys: List[str],
    targets: List[float],
) -> List[Dict[str, object]]:
    groups: Dict[tuple, List[Dict[str, object]]] = {}
    for row in rows:
        key = tuple(row[k] for k in group_keys)
        groups.setdefault(key, []).append(row)

    out: List[Dict[str, object]] = []
    for key, items in groups.items():
        first = items[0]
        for target in targets:
            item = {k: v for k, v in zip(group_keys, key)}
            item.update({
                "target_bler": float(target),
                "snr_at_target_bler_db": interpolate_target_snr(items, target_bler=float(target)),
                "delay_spread_ns": first.get("delay_spread_ns"),
                "cdd_delay_samples": first.get("cdd_delay_samples"),
                "dmrs_spacing_sc": first.get("dmrs_spacing_sc"),
                "algorithm": first.get("algorithm"),
                "precoder": first.get("precoder"),
                "n_trials_per_snr": first.get("n_trials"),
            })
            out.append(item)
    return out


def choose_throughput_snr(rows: List[Dict[str, object]], args: argparse.Namespace) -> float:
    if str(args.throughput_snr).lower() != "auto":
        return float(args.throughput_snr)
    ideal_rows = [r for r in rows if str(r["algorithm"]) == "Ideal CSI"]
    if not ideal_rows:
        return float(parse_float_list(args.compare_snrs)[0])
    target = float(parse_float_list(args.target_blers)[0])
    best = min(ideal_rows, key=lambda r: abs(float(r["bler"]) - target))
    return float(best["snr_db"])


def plot_ideal_bler(rows: List[Dict[str, object]], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    delays = sorted({int(r["cdd_delay_samples"]) for r in rows})
    for delay in delays:
        pts = sorted([r for r in rows if int(r["cdd_delay_samples"]) == delay], key=lambda r: float(r["snr_db"]))
        if not pts:
            continue
        floor = 0.5 / max(float(pts[0]["n_trials"]), 1.0)
        ys = [max(float(r["bler"]), floor) for r in pts]
        ax.semilogy([float(r["snr_db"]) for r in pts], ys, marker="o", linewidth=1.35, label=cdd_label(delay))
    ax.axhline(0.1, color="black", linestyle=":", linewidth=0.9)
    ax.axhline(0.01, color="black", linestyle="--", linewidth=0.9)
    ax.set_xlabel("SNR (dB)")
    ax.set_ylabel("TB BLER")
    ax.set_ylim(5e-3, 1.0)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_compare_bler(rows: List[Dict[str, object]], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    spacings = sorted({int(r["dmrs_spacing_sc"]) for r in rows})
    fig, axes = plt.subplots(1, len(spacings), figsize=(5.0 * len(spacings), 4.8), sharey=True)
    if len(spacings) == 1:
        axes = [axes]
    for ax, spacing in zip(axes, spacings):
        sub = [r for r in rows if int(r["dmrs_spacing_sc"]) == spacing]
        algs = []
        for row in sub:
            alg = str(row["algorithm"])
            if alg not in algs:
                algs.append(alg)
        for alg in algs:
            pts = sorted([r for r in sub if str(r["algorithm"]) == alg], key=lambda r: float(r["snr_db"]))
            if not pts:
                continue
            floor = 0.5 / max(float(pts[0]["n_trials"]), 1.0)
            ys = [max(float(r["bler"]), floor) for r in pts]
            ax.semilogy([float(r["snr_db"]) for r in pts], ys, marker="o", linewidth=1.35, label=alg)
        ax.axhline(0.1, color="black", linestyle=":", linewidth=0.8)
        ax.axhline(0.01, color="black", linestyle="--", linewidth=0.8)
        ax.set_title(f"DMRS spacing={spacing} sc", fontsize=10)
        ax.set_xlabel("SNR (dB)")
        ax.grid(True, which="both", alpha=0.3)
    axes[0].set_ylabel("TB BLER")
    axes[0].set_ylim(5e-3, 1.0)
    axes[-1].legend(fontsize=7, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_throughput_bar(rows: List[Dict[str, object]], snr_db: float, path: Path) -> List[Dict[str, object]]:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    selected = [
        r for r in rows
        if abs(float(r["snr_db"]) - float(snr_db)) < 1e-9
    ]
    selected = sorted(selected, key=lambda r: (int(r["dmrs_spacing_sc"]), str(r["algorithm"])))
    if not selected:
        return []
    labels = [f"sf{int(r['dmrs_spacing_sc'])}\n{r['algorithm']}" for r in selected]
    vals = [float(r["throughput_bits_per_slot"]) for r in selected]
    colors = [
        "#4c78a8" if str(r["algorithm"]).startswith("Ideal") else
        "#f58518" if str(r["algorithm"]).startswith("Alg1") else
        "#54a24b"
        for r in selected
    ]
    fig, ax = plt.subplots(figsize=(max(9.0, 0.72 * len(selected)), 5.2))
    x = np.arange(len(selected))
    ax.bar(x, vals, color=colors)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Throughput bits/slot")
    ax.set_title(f"DMRS-overhead-adjusted throughput @ SNR={float(snr_db):g} dB", fontsize=11)
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return selected


def run(args: argparse.Namespace) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out) / f"experiment17_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    ideal_rows, ideal_target_rows, selected_cdd_delay = run_ideal_cdd_scan(args, out_dir)
    compare_rows, compare_target_rows = run_algorithm_comparison(args, out_dir, selected_cdd_delay)
    throughput_snr = choose_throughput_snr(compare_rows, args)
    throughput_rows = plot_throughput_bar(
        compare_rows,
        throughput_snr,
        fig_dir / "experiment17_throughput_bar.png",
    )

    write_csv(ideal_rows, out_dir / "ideal_csi_cdd_scan.csv")
    write_csv(ideal_target_rows, out_dir / "ideal_csi_cdd_target_snr.csv")
    write_csv(compare_rows, out_dir / "algorithm_comparison.csv")
    write_csv(compare_target_rows, out_dir / "algorithm_comparison_target_snr.csv")
    write_csv(throughput_rows, out_dir / "throughput_bar_rows.csv")

    plot_ideal_bler(ideal_rows, fig_dir / "experiment17_ideal_cdd_bler.png")
    plot_compare_bler(compare_rows, fig_dir / "experiment17_alg1_alg3_bler.png")

    with (out_dir / "run_config.json").open("w", encoding="utf-8") as f:
        json.dump({
            **vars(args),
            "selected_cdd_delay_samples": int(selected_cdd_delay),
            "throughput_bar_snr_db": float(throughput_snr),
            "throughput_definition": "reference_tbs_no_dmrs_bits * (1 - bler) * data_re / total_re",
            "notes": [
                "The factor data_re/total_re includes data REs on DMRS symbols that are not occupied by pilot REs.",
                "Alg3 uses equal-total-overhead port DMRS: two ports split the same combined pilot positions used by Alg1.",
            ],
        }, f, indent=2, ensure_ascii=False)

    print(out_dir)
    return out_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/experiment17_cdd_diversity_bler")
    parser.add_argument("--fig-dir", default="docs/figures")
    parser.add_argument("--delay-spread-ns", type=float, default=10.0)
    parser.add_argument("--max-delay-factor", type=float, default=8.0)
    parser.add_argument("--ideal-cdd-delays", default="0,64,128,256,512")
    parser.add_argument("--ideal-dmrs-spacing-sc", type=int, default=24)
    parser.add_argument("--ideal-snrs", default="2,3,4,5,6")
    parser.add_argument("--compare-dmrs-spacings", default="24,12,6")
    parser.add_argument("--compare-snrs", default="2,3,4,5,6")
    parser.add_argument("--target-blers", default="0.1,0.01")
    parser.add_argument("--selection-targets", default="0.01,0.1")
    parser.add_argument("--throughput-snr", default="auto")
    parser.add_argument("--n-prbs", type=int, default=48)
    parser.add_argument("--n-fft", type=int, default=4096)
    parser.add_argument("--pdsch-symbols", type=int, default=10)
    parser.add_argument("--dmrs-symbols", default="2,7")
    parser.add_argument("--mcs-index", type=int, default=8)
    parser.add_argument("--trials", type=int, default=100)
    parser.add_argument("--ldpc-iterations", type=int, default=8)
    parser.add_argument("--llr-clip", type=float, default=50.0)
    parser.add_argument("--loading", type=float, default=1e-8)
    parser.add_argument("--seed", type=int, default=20260615)
    parser.add_argument("--progress-every", type=int, default=25)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
