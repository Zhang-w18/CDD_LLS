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
from cdd_lls.phy.precoding import equivalent_channel
from cdd_lls.phy.qam import qam_demapper_maxlog, qam_modulate
from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
from cdd_lls.sim.stats import interpolate_target_snr, snr_values
from tools.search_precoder_design_alg1 import make_alg1_full_cov_estimator
from tools.search_precoder_design_alg2 import PrecoderDesign, cdd_design, rb_qpsk_design
from tools.search_precoder_diversity_tradeoff import (
    hybrid_cdd_rb_design,
    power_cycle_design,
    rb_mpsk_design,
    sinusoid_phase_design,
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


def gain_db(reference: float, candidate: float) -> float:
    if not np.isfinite(reference) or not np.isfinite(candidate) or reference <= 0.0 or candidate <= 0.0:
        return float("nan")
    return float(10.0 * np.log10(reference / candidate))


def split_llrs(llr: np.ndarray, cb_e_values: List[int]) -> List[np.ndarray]:
    arr = np.asarray(llr, dtype=np.float64).reshape(-1)
    out = []
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
    """Decode multiple candidate LLR streams sharing the same TB layout."""
    if not llrs_by_candidate:
        return []
    split_by_candidate = [
        split_llrs(llr, adapter.cb_e_values) for llr in llrs_by_candidate
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


def build_designs(
    grid,
    cdd_delays: List[int],
    rb_group: int,
    candidate_set: str = "rb12",
) -> List[PrecoderDesign]:
    designs = [cdd_design(grid, d) for d in cdd_delays]
    designs.append(rb_qpsk_design(grid, rb_group=int(rb_group), label=f"RB{int(rb_group)} QPSK cycle"))
    if str(candidate_set).lower() in ("diversity16", "experiment16", "exp16"):
        designs.extend([
            rb_mpsk_design(grid, m=16, rb_group=1, step=1),
            rb_mpsk_design(grid, m=8, rb_group=1, step=1),
            rb_mpsk_design(grid, m=8, rb_group=2, step=1),
            hybrid_cdd_rb_design(grid, delay=32, rb_group=12, m=8),
            hybrid_cdd_rb_design(grid, delay=32, rb_group=4, m=8),
            hybrid_cdd_rb_design(grid, delay=48, rb_group=8, m=4),
            power_cycle_design(grid, rb_group=2, powers=[0.25, 0.75], phase_m=4),
            sinusoid_phase_design(grid, cycles=2, beta_pi=2.0),
        ])
    unique: Dict[str, PrecoderDesign] = {}
    for design in designs:
        unique.setdefault(design.label, design)
    return list(unique.values())


def plot_metric(rows: List[Dict[str, object]], metric: str, path: Path, ylabel: str) -> None:
    if not rows:
        return
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    groups: Dict[str, List[Dict[str, object]]] = {}
    for row in rows:
        key = f"DS={row['delay_spread_ns']:g}ns | {row['precoder']}"
        groups.setdefault(key, []).append(row)
    fig, ax = plt.subplots(figsize=(10, 5.8))
    for label, items in sorted(groups.items()):
        pts = sorted(items, key=lambda r: float(r["snr_db"]))
        ys = [max(float(r[metric]), 1e-4 if metric == "bler" else 1e-9) for r in pts]
        if metric in ("bler", "ce_nmse_eff"):
            ax.semilogy([float(r["snr_db"]) for r in pts], ys, marker="o", linewidth=1.3, label=label)
        else:
            ax.plot([float(r["snr_db"]) for r in pts], ys, marker="o", linewidth=1.3, label=label)
    ax.grid(True, which="both", alpha=0.3)
    ax.set_xlabel("SNR (dB)")
    ax.set_ylabel(ylabel)
    if metric == "bler":
        ax.set_ylim(1e-4, 1.0)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out) / f"bler_precoder_alg1_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    delay_spreads = parse_float_list(args.delay_spreads_ns)
    cdd_delays = parse_int_list(args.cdd_delays)
    snrs = (
        parse_float_list(args.snrs)
        if str(args.snrs).strip()
        else snr_values(parse_float_list(args.snr_range_db))
    )
    rows: List[Dict[str, object]] = []
    trial_rows: List[Dict[str, object]] = []

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
    mcs = get_mcs("nr_256qam", int(args.mcs_index), None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(
        tb.cb_k_values,
        tb.cb_e_values,
        num_iter=int(args.ldpc_iterations),
        llr_clip=float(args.llr_clip),
    )
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)

    for delay_spread_ns in delay_spreads:
        sample_period_ns = 1e9 / (float(resource.n_fft) * float(resource.scs_khz) * 1e3)
        pdp = make_exponential_pdp(delay_spread_ns, sample_period_ns, max_delay_factor=8.0)
        channel_cfg = ChannelConfig(delay_spread_ns=float(delay_spread_ns), max_delay_factor=8.0)
        designs = build_designs(grid, cdd_delays, int(args.rb_group), str(args.candidate_set))

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
                    "tb_errors": 0,
                    "cb_errors": 0,
                    "goodput_bits": 0,
                    "ce_nmse": [],
                    "design": design,
                    "rpp_loaded_cond": estimators[design.label].metadata["rpp_loaded_cond"],
                }
                for design in designs
            }
            base_seed = stable_seed(delay_spread_ns, snr_db, base=int(args.seed))
            print(
                f"[start] DS={float(delay_spread_ns):g} ns, "
                f"SNR={float(snr_db):g} dB, trials={int(args.trials)}",
                flush=True,
            )
            for trial in range(1, int(args.trials) + 1):
                rng = np.random.default_rng(stable_seed(base_seed, trial))
                realization = generate_tdl_channel(rng, grid, channel_cfg, n_tx=2, n_rx=4)
                H = realization.H
                payload = [
                    rng.integers(0, 2, size=int(k), dtype=np.int8)
                    for k in tb.cb_k_values
                ]
                coded = adapter.encode(payload)
                coded_cw = np.concatenate(coded)
                symbols = qam_modulate(coded_cw, int(mcs.qm))
                pilot_noise = (
                    rng.normal(size=(4, len(pilot_local)))
                    + 1j * rng.normal(size=(4, len(pilot_local)))
                ) * math.sqrt(noise_var_ls / 2.0)
                data_noise = (
                    rng.normal(size=(4, len(data_local)))
                    + 1j * rng.normal(size=(4, len(data_local)))
                ) * math.sqrt(noise_var / 2.0)

                trial_designs = []
                trial_llrs = []
                trial_nmse = []
                for design in designs:
                    g = equivalent_channel(H, design.C)
                    ls_obs = g[:, pilot_local] + pilot_noise
                    g_hat = ls_obs @ estimators[design.label].matrix.T
                    y = g[:, data_local] * symbols[None, :] + data_noise
                    z, no_eff = equalize_mrc(y, g_hat[:, data_local], noise_var)
                    llr_cw = qam_demapper_maxlog(z, no_eff, int(mcs.qm))
                    trial_designs.append(design)
                    trial_llrs.append(llr_cw)
                    trial_nmse.append(nmse(g, g_hat))

                dec_results = decode_same_tb_batch(adapter, trial_llrs, payload)
                for design, dec, ce in zip(trial_designs, dec_results, trial_nmse):
                    acc = accum[design.label]
                    acc["tb_errors"] += int(not dec.tb_success)
                    acc["cb_errors"] += int(sum(1 for ok in dec.cb_success if not ok))
                    acc["goodput_bits"] += int(dec.goodput_bits)
                    acc["ce_nmse"].append(ce)
                    if bool(args.save_trial_metrics):
                        trial_rows.append({
                            "delay_spread_ns": float(delay_spread_ns),
                            "snr_db": float(snr_db),
                            "trial": int(trial),
                            "precoder": design.label,
                            "tb_error": int(not dec.tb_success),
                            "cb_errors": int(sum(1 for ok in dec.cb_success if not ok)),
                            "ce_nmse_eff": float(ce),
                        })
                if int(args.progress_every) > 0 and trial % int(args.progress_every) == 0:
                    short = ", ".join(
                        f"{label}:BLER={values['tb_errors'] / trial:.3f}"
                        for label, values in sorted(accum.items())
                    )
                    print(
                        f"[progress] DS={float(delay_spread_ns):g} ns, "
                        f"SNR={float(snr_db):g} dB, trial={trial}/{int(args.trials)} | {short}",
                        flush=True,
                    )

            cdd_blers = {
                label: float(values["tb_errors"]) / float(args.trials)
                for label, values in accum.items()
                if values["design"].family == "cdd_linear_phase"
            }
            cdd_nmse = {
                label: mean_or_nan(values["ce_nmse"])
                for label, values in accum.items()
                if values["design"].family == "cdd_linear_phase"
            }
            best_cdd_bler_label = min(cdd_blers, key=cdd_blers.get)
            best_cdd_nmse_label = min(cdd_nmse, key=cdd_nmse.get)
            for label, values in accum.items():
                design = values["design"]
                n_cb_trials = int(int(args.trials) * len(tb.cb_k_values))
                bler = float(values["tb_errors"]) / float(args.trials)
                cb_bler = float(values["cb_errors"]) / float(n_cb_trials)
                ce_nmse = mean_or_nan(values["ce_nmse"])
                rows.append({
                    "scenario_id": f"ds{delay_spread_ns:g}",
                    "precoder": label,
                    "family": design.family,
                    "delay_spread_ns": float(delay_spread_ns),
                    "snr_db": float(snr_db),
                    "n_trials": int(args.trials),
                    "tb_errors": int(values["tb_errors"]),
                    "cb_errors": int(values["cb_errors"]),
                    "bler": bler,
                    "cb_bler": cb_bler,
                    "ce_nmse_eff": ce_nmse,
                    "best_cdd_bler_precoder": best_cdd_bler_label,
                    "best_cdd_bler": cdd_blers[best_cdd_bler_label],
                    "best_cdd_nmse_precoder": best_cdd_nmse_label,
                    "best_cdd_nmse": cdd_nmse[best_cdd_nmse_label],
                    "bler_delta_vs_best_cdd": float(bler - cdd_blers[best_cdd_bler_label]),
                    "ce_nmse_gain_vs_best_cdd_db": gain_db(cdd_nmse[best_cdd_nmse_label], ce_nmse),
                    "rpp_loaded_cond": float(values["rpp_loaded_cond"]),
                    "goodput_bits_per_slot": float(values["goodput_bits"]) / float(args.trials),
                    "goodput_se_per_re": float(values["goodput_bits"]) / float(max(int(args.trials) * grid.n_data_re, 1)),
                    "n_tx": 2,
                    "n_rx": 4,
                    "n_prbs": int(args.n_prbs),
                    "n_sc": int(grid.n_sc),
                    "dmrs_spacing_sc": int(args.dmrs_spacing_sc),
                    "dmrs_symbols": int(len(resource.dmrs_symbol_indices)),
                    "dmrs_overhead": float(grid.dmrs_overhead),
                    "data_re": int(grid.n_data_re),
                    "mcs_table": "nr_256qam",
                    "mcs_index": int(args.mcs_index),
                    "qm": int(mcs.qm),
                    "code_rate": float(mcs.code_rate),
                    "tbs_bits": int(tb.tb_size),
                    "coded_bits": int(tb.coded_bits),
                    "n_cbs": int(len(tb.cb_k_values)),
                    "ldpc_iterations": int(args.ldpc_iterations),
                })
            write_csv(rows, out_dir / "bler_summary_partial.csv")
            write_csv(trial_rows, out_dir / "trial_metrics_partial.csv")
            short = ", ".join(
                f"{label}:BLER={values['tb_errors'] / float(args.trials):.3f},"
                f"NMSE={mean_or_nan(values['ce_nmse']):.3e}"
                for label, values in sorted(accum.items())
            )
            print(
                f"[done] DS={float(delay_spread_ns):g} ns, "
                f"SNR={float(snr_db):g} dB | {short}",
                flush=True,
            )

    write_csv(rows, out_dir / "bler_summary.csv")
    write_csv(trial_rows, out_dir / "trial_metrics.csv")
    target_rows = []
    for delay_spread_ns in delay_spreads:
        for precoder in sorted({str(r["precoder"]) for r in rows if float(r["delay_spread_ns"]) == float(delay_spread_ns)}):
            pts = [r for r in rows if float(r["delay_spread_ns"]) == float(delay_spread_ns) and str(r["precoder"]) == precoder]
            target_rows.append({
                "delay_spread_ns": float(delay_spread_ns),
                "precoder": precoder,
                "target_bler": float(args.target_bler),
                "snr_at_target_bler_db": interpolate_target_snr(pts, target_bler=float(args.target_bler)),
            })
    write_csv(target_rows, out_dir / "summary_10pct_bler_snr.csv")
    with (out_dir / "run_config.json").open("w", encoding="utf-8") as f:
        json.dump(vars(args), f, indent=2, ensure_ascii=False)
    plot_metric(rows, "bler", fig_dir / "experiment15_alg1_precoder_bler.png", "TB BLER")
    plot_metric(rows, "ce_nmse_eff", fig_dir / "experiment15_alg1_precoder_ce_nmse.png", "Equivalent-channel NMSE")
    print(out_dir)
    return out_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/bler_precoder_alg1")
    parser.add_argument("--fig-dir", default="docs/figures")
    parser.add_argument("--delay-spreads-ns", default="30,100")
    parser.add_argument("--snr-range-db", default="-2,8,2")
    parser.add_argument("--snrs", default="")
    parser.add_argument("--trials", type=int, default=20)
    parser.add_argument("--n-prbs", type=int, default=48)
    parser.add_argument("--dmrs-spacing-sc", type=int, default=6)
    parser.add_argument("--mcs-index", type=int, default=8)
    parser.add_argument("--cdd-delays", default="64,128")
    parser.add_argument("--rb-group", type=int, default=12)
    parser.add_argument("--candidate-set", default="rb12")
    parser.add_argument("--ldpc-iterations", type=int, default=8)
    parser.add_argument("--llr-clip", type=float, default=50.0)
    parser.add_argument("--loading", type=float, default=1e-8)
    parser.add_argument("--seed", type=int, default=20260614)
    parser.add_argument("--target-bler", type=float, default=0.10)
    parser.add_argument("--save-trial-metrics", action="store_true")
    parser.add_argument("--progress-every", type=int, default=25)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
