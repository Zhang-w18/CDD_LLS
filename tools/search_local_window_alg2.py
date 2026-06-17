from __future__ import annotations

import argparse
import csv
import json
import math
import os
from datetime import datetime
from pathlib import Path
import sys
from typing import Dict, Iterable, List, Optional

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ChannelConfig, ResourceConfig
from cdd_lls.phy.channel_tdl import generate_tdl_channel, make_exponential_pdp
from cdd_lls.phy.estimators import nmse, shifted_pdp
from cdd_lls.phy.precoding import cdd_equivalent_from_branches
from cdd_lls.phy.resource_grid import ResourceGrid, build_resource_grid, local_indices_for_subcarriers
from tools.run_nmse_algorithm_comparison import (
    BlockEstimator,
    apply_block,
    first_below,
    freq_corr_from_pdp,
    make_direct_estimator,
    mean_or_nan,
    median_or_nan,
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


def gain_db(reference: float, candidate: float) -> float:
    if not np.isfinite(reference) or not np.isfinite(candidate) or reference <= 0.0 or candidate <= 0.0:
        return float("nan")
    return float(10.0 * np.log10(reference / candidate))


def cdd_phase_span_rad(window_sc: int, cdd_delay: int, n_fft: int) -> float:
    return float(2.0 * np.pi * abs(int(cdd_delay)) * max(int(window_sc) - 1, 0) / float(n_fft))


def make_local_window_estimator(
    grid: ResourceGrid,
    delays: List[int],
    noise_var_ls: float,
    loading: float,
    window_sc: int,
) -> BlockEstimator:
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    delays_arr = np.asarray(delays, dtype=np.float64)
    n_tx = len(delays)
    window = max(1, min(int(window_sc), int(grid.n_sc)))
    blocks: List[Dict[str, object]] = []
    conds: List[float] = []
    pilot_counts: List[int] = []

    for start in range(0, grid.n_sc, window):
        stop = min(grid.n_sc, start + window)
        target_local = np.arange(start, stop, dtype=np.int64)
        mask = (pilot_local >= start) & (pilot_local < stop)
        pilot_idx = np.where(mask)[0]
        if len(pilot_idx) == 0:
            center = 0.5 * (grid.subcarrier_indices[start] + grid.subcarrier_indices[stop - 1])
            nearest = int(np.argmin(np.abs(grid.pilot_subcarriers.astype(np.float64) - center)))
            pilot_idx = np.asarray([nearest], dtype=np.int64)
        pilot_k = grid.pilot_subcarriers[pilot_idx].astype(np.float64)
        A = np.exp(
            -1j * 2.0 * np.pi * pilot_k[:, None] * delays_arr[None, :] / float(grid.n_fft)
        ) / math.sqrt(float(n_tx))
        S = A @ A.conj().T
        S += (float(noise_var_ls) + float(loading)) * np.eye(S.shape[0], dtype=np.complex128)
        M = A.conj().T @ np.linalg.solve(S, np.eye(S.shape[0], dtype=np.complex128))
        cond = float(np.linalg.cond(A)) if len(pilot_k) >= n_tx else float("inf")
        conds.append(cond)
        pilot_counts.append(int(len(pilot_idx)))
        blocks.append({
            "start": int(start),
            "stop": int(stop),
            "target_local": target_local,
            "pilot_idx": pilot_idx,
            "matrix": M,
            "cond": cond,
            "pilot_count": int(len(pilot_idx)),
        })

    finite = [x for x in conds if np.isfinite(x)]
    bad = [p for p in pilot_counts if p < n_tx]
    return BlockEstimator(
        label=f"Alg2 LocalWindow W{window}sc",
        blocks=blocks,
        block_width_sc=int(window),
        cond_max=float(max(finite)) if finite else float("inf"),
        cond_mean=float(np.mean(finite)) if finite else float("inf"),
        metadata={
            "window_sc": int(window),
            "n_blocks": int(len(blocks)),
            "cond_max": float(max(finite)) if finite else float("inf"),
            "cond_mean": float(np.mean(finite)) if finite else float("inf"),
            "min_pilots_per_window": int(min(pilot_counts)) if pilot_counts else 0,
            "mean_pilots_per_window": float(np.mean(pilot_counts)) if pilot_counts else 0.0,
            "bad_window_fraction": float(len(bad) / max(len(pilot_counts), 1)),
        },
    )


def plot_best(rows: List[Dict[str, object]], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    candidates = [
        r for r in rows
        if str(r["algorithm"]).startswith("Alg2 LocalWindow")
        and int(float(r.get("window_sc", 0))) < int(float(r.get("n_sc", 0)))
    ]
    candidates = [r for r in candidates if np.isfinite(float(r.get("gain_vs_alg1_db", "nan")))]
    candidates.sort(key=lambda r: float(r["gain_vs_alg1_db"]), reverse=True)
    selected = candidates[:30]
    if not selected:
        return
    labels = [
        f"W={r['window_sc']} DS={r['delay_spread_ns']} d={r['cdd_delay_samples']} sf={r['dmrs_spacing_sc']} SNR={r['snr_db']}"
        for r in selected
    ]
    gains = [float(r["gain_vs_alg1_db"]) for r in selected]
    colors = ["#2a9d8f" if g > 0 else "#b56576" for g in gains]
    fig, ax = plt.subplots(figsize=(12, 8.4))
    y = np.arange(len(selected))
    ax.barh(y, gains, color=colors)
    ax.axvline(0.0, color="black", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("NMSE gain vs matched Alg1 Direct RMMSE (dB)")
    ax.set_title("Best strict local-window Alg2 cases (W < allocation)")
    ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_window_tradeoff(rows: List[Dict[str, object]], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/cdd_lls_mplconfig")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    available = []
    for r in rows:
        key = (
            float(r["delay_spread_ns"]),
            int(float(r["cdd_delay_samples"])),
            int(float(r["dmrs_spacing_sc"])),
            float(r["snr_db"]),
        )
        if key not in available:
            available.append(key)
    preferred = [
        (1.0, 64, 24, 8.0),
        (1.0, 512, 24, 8.0),
        (5.0, 128, 12, 8.0),
        (10.0, 256, 12, 8.0),
        (30.0, 512, 24, 8.0),
        (30.0, 512, 12, 8.0),
        (100.0, 32, 24, -6.0),
        (100.0, 512, 24, 8.0),
    ]
    selected_keys = [key for key in preferred if key in available]
    for key in available:
        if len(selected_keys) >= 8:
            break
        if key not in selected_keys:
            selected_keys.append(key)
    if not selected_keys:
        return

    fig, axes = plt.subplots(len(selected_keys), 1, figsize=(10.5, max(2.5 * len(selected_keys), 5.0)), sharex=True)
    if len(selected_keys) == 1:
        axes = [axes]

    for ax, key in zip(axes, selected_keys):
        ds, cdd, sf, snr = key
        pts = [
            r for r in rows
            if float(r["delay_spread_ns"]) == ds
            and int(r["cdd_delay_samples"]) == cdd
            and int(r["dmrs_spacing_sc"]) == sf
            and float(r["snr_db"]) == snr
            and str(r["algorithm"]).startswith("Alg2 LocalWindow")
        ]
        pts.sort(key=lambda r: int(r["window_sc"]))
        if not pts:
            continue
        x = [int(r["window_sc"]) for r in pts]
        y = [float(r["gain_vs_alg1_db"]) for r in pts]
        ax.plot(x, y, marker="o", linewidth=1.2, markersize=3.2, label="gain")
        ax.axhline(0.0, color="black", linewidth=0.8)
        bc_h_09 = pts[0].get("bc_h_0p9_sc", "")
        bc_h_05 = pts[0].get("bc_h_0p5_sc", "")
        for bc, style, name in [(bc_h_09, ":", "BcH0.9"), (bc_h_05, "--", "BcH0.5")]:
            if str(bc) != "":
                ax.axvline(float(bc), color="gray", linestyle=style, linewidth=0.9)
                ax.text(float(bc), ax.get_ylim()[1], name, fontsize=7, va="top", ha="right", rotation=90)
        ax.set_title(f"DS={ds:g} ns, d={cdd}, sf={sf}, SNR={snr:g} dB", loc="left", fontsize=9)
        ax.set_ylabel("gain dB")
        ax.grid(True, alpha=0.25)
    axes[-1].set_xlabel("local window width (subcarriers)")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out) / f"local_window_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    delay_spreads = parse_float_list(args.delay_spreads_ns)
    cdd_delays = parse_int_list(args.cdd_delays)
    dmrs_spacings = parse_int_list(args.dmrs_spacings)
    snrs = parse_float_list(args.snrs)
    windows = parse_int_list(args.windows_sc)
    rows: List[Dict[str, object]] = []
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
                    block_ests = [
                        make_local_window_estimator(grid, delays, noise_var_ls, float(args.loading), w)
                        for w in windows
                    ]
                    accum: Dict[str, Dict[str, object]] = {
                        direct.label: {"eff": [], "branch": [], "meta": direct.metadata, "family": direct.family}
                    }
                    for est in block_ests:
                        accum[est.label] = {"eff": [], "branch": [], "meta": est.metadata, "family": "alg2_local_window"}

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
                        accum[direct.label]["eff"].append(nmse(g, g_hat))
                        accum[direct.label]["branch"].append(float("nan"))

                        for est in block_ests:
                            g_hat, h_hat = apply_block(est, ls_cdd, grid, delays)
                            accum[est.label]["eff"].append(nmse(g, g_hat))
                            accum[est.label]["branch"].append(nmse(H, h_hat))

                    direct_nmse = mean_or_nan(accum[direct.label]["eff"])
                    for label, values in accum.items():
                        meta = dict(values["meta"])
                        eff = mean_or_nan(values["eff"])
                        branch = mean_or_nan(values["branch"])
                        window_sc = meta.get("window_sc", "")
                        phase_span = "" if window_sc == "" else cdd_phase_span_rad(int(window_sc), int(cdd_delay), grid.n_fft)
                        flatness_rho = ""
                        if window_sc != "":
                            idx = min(max(int(window_sc) - 1, 0), len(rho_h) - 1)
                            flatness_rho = float(rho_h[idx])
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
                            "pilot_count_combined": int(grid.pilot_count),
                            "bc_h_0p9_sc": "" if bc_h_09 is None else int(bc_h_09),
                            "bc_h_0p5_sc": "" if bc_h_05 is None else int(bc_h_05),
                            "bc_g_0p9_sc": "" if bc_g_09 is None else int(bc_g_09),
                            "bc_g_0p5_sc": "" if bc_g_05 is None else int(bc_g_05),
                            "ce_nmse_eff_mean": eff,
                            "ce_nmse_eff_median": median_or_nan(values["eff"]),
                            "ce_nmse_branch_mean": branch,
                            "ce_nmse_branch_median": median_or_nan(values["branch"]),
                            "gain_vs_alg1_db": gain_db(direct_nmse, eff),
                            "window_phase_span_rad": phase_span,
                            "window_phase_span_pi": "" if phase_span == "" else float(phase_span) / math.pi,
                            "rho_h_at_window_edge": flatness_rho,
                        }
                        row.update(meta)
                        rows.append(row)
                        if label != direct.label and np.isfinite(float(row["gain_vs_alg1_db"])) and float(row["gain_vs_alg1_db"]) > 0.0:
                            win_rows.append(row)

    write_csv(rows, out_dir / "local_window_nmse_search.csv")
    win_rows.sort(key=lambda r: float(r["gain_vs_alg1_db"]), reverse=True)
    write_csv(win_rows, out_dir / "local_window_wins.csv")
    with (out_dir / "run_config.json").open("w", encoding="utf-8") as f:
        json.dump({
            "delay_spreads_ns": delay_spreads,
            "cdd_delays": cdd_delays,
            "dmrs_spacings": dmrs_spacings,
            "snrs": snrs,
            "windows_sc": windows,
            "trials": int(args.trials),
            "seed": int(args.seed),
            "note": "Targeted NMSE-only search for local-window deterministic CDD-aware estimator.",
        }, f, indent=2)

    plot_best(rows, fig_dir / "experiment12_local_window_best_gains.png")
    plot_window_tradeoff(rows, fig_dir / "experiment12_local_window_tradeoff.png")
    print(out_dir)
    return out_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/local_window_alg2_search")
    parser.add_argument("--fig-dir", default="docs/figures")
    parser.add_argument("--delay-spreads-ns", default="1,5,10,30,100")
    parser.add_argument("--cdd-delays", default="32,64,128,256,512")
    parser.add_argument("--dmrs-spacings", default="4,6,12,24")
    parser.add_argument("--snrs", default="-6,0,8")
    parser.add_argument("--windows-sc", default="12,18,24,36,48,72,96,144,192,288,576")
    parser.add_argument("--trials", type=int, default=80)
    parser.add_argument("--seed", type=int, default=20260614)
    parser.add_argument("--loading", type=float, default=1e-8)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
