from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ChannelConfig, ChannelEstimationConfig, ResourceConfig
from cdd_lls.phy.channel_tdl import make_exponential_pdp
from cdd_lls.phy.estimators import shifted_pdp
from cdd_lls.phy.resource_grid import build_resource_grid


def support_indices(pdp: np.ndarray, threshold: float = 0.99) -> np.ndarray:
    c = np.cumsum(np.asarray(pdp, dtype=np.float64))
    cutoff = int(np.searchsorted(c, float(threshold), side="left")) + 1
    return np.arange(max(1, min(cutoff, len(pdp))), dtype=np.int64)


def build_phi(grid, delays: List[int], support: np.ndarray) -> np.ndarray:
    n_tx = len(delays)
    pilot_k = grid.pilot_subcarriers.astype(np.float64)
    phi = np.empty((len(pilot_k), n_tx * len(support)), dtype=np.complex128)
    col = 0
    for d in delays:
        shifted = support.astype(np.float64) + float(d)
        phi[:, col:col + len(support)] = np.exp(
            -1j * 2.0 * np.pi * pilot_k[:, None] * shifted[None, :] / float(grid.n_fft)
        ) / math.sqrt(float(n_tx))
        col += len(support)
    return phi


def svd_metrics(a: np.ndarray, eps_scale: float = 1e-12) -> Dict[str, float]:
    s = np.linalg.svd(np.asarray(a, dtype=np.complex128), compute_uv=False)
    if len(s) == 0:
        return {
            "rank": 0,
            "nullity": 0,
            "cond_nonzero": float("nan"),
            "s_max": float("nan"),
            "s_min_nonzero": float("nan"),
            "effective_rank": float("nan"),
        }
    tol = float(s[0]) * max(a.shape) * np.finfo(float).eps
    rank = int(np.sum(s > tol))
    nullity = int(a.shape[1] - rank)
    s_min = float(s[rank - 1]) if rank > 0 else float("nan")
    cond = float(s[0] / s_min) if rank > 0 and s_min > 0 else float("inf")
    eff = float((np.sum(s) ** 2) / max(np.sum(s ** 2), eps_scale))
    return {
        "rank": rank,
        "nullity": nullity,
        "cond_nonzero": cond,
        "s_max": float(s[0]),
        "s_min_nonzero": s_min,
        "effective_rank": eff,
    }


def freq_corr_from_pdp(pdp: np.ndarray, n_fft: int, max_delta_sc: int) -> np.ndarray:
    p = np.asarray(pdp, dtype=np.float64)
    taps = np.arange(len(p), dtype=np.float64)
    deltas = np.arange(int(max_delta_sc) + 1, dtype=np.float64)
    r = np.exp(-1j * 2.0 * np.pi * deltas[:, None] * taps[None, :] / float(n_fft)) @ p
    return np.abs(r / max(float(np.sum(p)), 1e-30))


def first_below(rho: np.ndarray, threshold: float) -> int | None:
    idx = np.where(np.asarray(rho) <= float(threshold))[0]
    if len(idx) == 0:
        return None
    return int(idx[0])


def bc_label(rho: np.ndarray, threshold: float, scs_khz: int) -> str:
    sc = first_below(rho, threshold)
    if sc is None:
        return "> allocation"
    return f"{sc} sc / {sc * scs_khz / 1000.0:.2f} MHz"


def coherence_row(label: str, pdp: np.ndarray, grid, scs_khz: int, max_delta_sc: int) -> Dict[str, object]:
    rho = freq_corr_from_pdp(pdp, grid.n_fft, max_delta_sc)
    row: Dict[str, object] = {
        "label": label,
        "pdp_len": int(len(pdp)),
        "rho_edge": float(rho[-1]),
        "rho_min_in_allocation": float(np.min(rho[: max_delta_sc + 1])),
        "rho_mean_in_allocation": float(np.mean(rho[: max_delta_sc + 1])),
    }
    for thr in (0.9, 0.5):
        sc = first_below(rho, thr)
        row[f"bc_{thr}_sc"] = "" if sc is None else int(sc)
        row[f"bc_{thr}_khz"] = "" if sc is None else float(sc * scs_khz)
        row[f"bc_{thr}_rb"] = "" if sc is None else float(sc / 12.0)
    return row


def write_csv(rows: List[Dict[str, object]], path: Path) -> None:
    if not rows:
        return
    fields = sorted(set().union(*(row.keys() for row in rows)))
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot_corr(curves: Dict[str, np.ndarray], scs_khz: int, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 5.2))
    colors = {}
    for label, rho in curves.items():
        x_mhz = np.arange(len(rho), dtype=np.float64) * float(scs_khz) / 1000.0
        line, = ax.plot(x_mhz, rho, linewidth=1.5, label=label)
        colors[label] = line.get_color()
    for thr, style in ((0.9, "--"), (0.5, ":")):
        ax.axhline(thr, color="gray", linestyle=style, linewidth=0.9)
        for label, rho in curves.items():
            sc = first_below(rho, thr)
            if sc is not None:
                ax.axvline(
                    sc * scs_khz / 1000.0,
                    color=colors[label],
                    linestyle=style,
                    linewidth=0.9,
                    alpha=0.45,
                )
    lines = []
    for label, rho in curves.items():
        lines.append(
            f"{label}: Bc0.9={bc_label(rho, 0.9, scs_khz)}, "
            f"Bc0.5={bc_label(rho, 0.5, scs_khz)}"
        )
    ax.text(
        0.98,
        0.96,
        "\n".join(lines),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        bbox={"facecolor": "white", "alpha": 0.78, "edgecolor": "none"},
    )
    ax.set_xlabel("Subcarrier separation (MHz)")
    ax.set_ylabel("|normalized covariance|")
    ax.set_ylim(0.0, 1.03)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_corr_summary(curves_by_case: Dict[str, Dict[str, np.ndarray]], scs_khz: int, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(curves_by_case)
    fig, axes = plt.subplots(n, 1, figsize=(10.5, max(3.0 * n, 5.0)), sharex=True)
    if n == 1:
        axes = [axes]
    for ax, (case_id, curves) in zip(axes, curves_by_case.items()):
        colors = {}
        for label, rho in curves.items():
            x_mhz = np.arange(len(rho), dtype=np.float64) * float(scs_khz) / 1000.0
            line, = ax.plot(x_mhz, rho, linewidth=1.6, label=label)
            colors[label] = line.get_color()
        for thr, style in ((0.9, "--"), (0.5, ":")):
            ax.axhline(thr, color="gray", linestyle=style, linewidth=0.8)
            for label, rho in curves.items():
                sc = first_below(rho, thr)
                if sc is not None:
                    ax.axvline(
                        sc * scs_khz / 1000.0,
                        color=colors[label],
                        linestyle=style,
                        linewidth=0.9,
                        alpha=0.45,
                    )
        lines = []
        for label, rho in curves.items():
            lines.append(
                f"{label}: Bc0.9={bc_label(rho, 0.9, scs_khz)}, "
                f"Bc0.5={bc_label(rho, 0.5, scs_khz)}"
            )
        ax.text(
            0.985,
            0.94,
            "\n".join(lines),
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.80, "edgecolor": "none"},
        )
        ax.set_title(case_id, loc="left", fontsize=10)
        ax.set_ylabel(r"$|R[\Delta k]|/R[0]$")
        ax.set_ylim(0.0, 1.03)
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8, loc="lower left")
    axes[-1].set_xlabel("Subcarrier separation (MHz), SCS=30 kHz")
    fig.suptitle("Experiment 8: frequency covariance magnitude and coherence bandwidth", y=0.998)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_original_pdp(delay_spreads_ns: List[float], sample_period_ns: float, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9.5, 5.4))
    for ds in sorted(set(float(x) for x in delay_spreads_ns)):
        pdp = make_exponential_pdp(
            delay_spread_ns=ds,
            sample_period_ns=sample_period_ns,
            max_delay_factor=8.0,
        )
        delays_ns = np.arange(len(pdp), dtype=np.float64) * float(sample_period_ns)
        support = support_indices(pdp, threshold=0.99)
        support_edge_ns = delays_ns[int(support[-1])]
        ax.plot(
            delays_ns,
            pdp,
            marker="o",
            markersize=2.8,
            linewidth=1.4,
            label=f"Original PDP, DS={ds:g} ns, 99% support={len(support)} taps",
        )
        ax.axvline(support_edge_ns, linestyle="--", linewidth=0.9, alpha=0.45)
        ax.text(
            support_edge_ns,
            max(pdp) * 0.65,
            f"99% @ {support_edge_ns:.1f} ns",
            rotation=90,
            va="center",
            ha="right",
            fontsize=8,
        )
    ax.set_yscale("log")
    ax.set_xlabel(f"Tap delay (ns), sample period={sample_period_ns:.3f} ns")
    ax.set_ylabel("Normalized tap power")
    ax.set_title("Experiment 8: original exponential TDL PDP before CDD")
    ax.grid(True, which="both", alpha=0.28)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/diagnostics/reconstruction_covariance")
    parser.add_argument("--fig-dir", default="docs/figures")
    args = parser.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        {"case_id": "ds30_d16_sf6", "delay_spread_ns": 30.0, "cdd_delay": 16, "dmrs_spacing_sc": 6},
        {"case_id": "ds30_d32_sf6", "delay_spread_ns": 30.0, "cdd_delay": 32, "dmrs_spacing_sc": 6},
        {"case_id": "ds100_d16_sf6", "delay_spread_ns": 100.0, "cdd_delay": 16, "dmrs_spacing_sc": 6},
        {"case_id": "ds100_d32_sf6", "delay_spread_ns": 100.0, "cdd_delay": 32, "dmrs_spacing_sc": 6},
        {"case_id": "ds30_d64_sf2", "delay_spread_ns": 30.0, "cdd_delay": 64, "dmrs_spacing_sc": 2},
    ]

    cond_rows: List[Dict[str, object]] = []
    coh_rows: List[Dict[str, object]] = []
    curves_by_case: Dict[str, Dict[str, np.ndarray]] = {}

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
        pdp = make_exponential_pdp(
            delay_spread_ns=float(case["delay_spread_ns"]),
            sample_period_ns=sample_period_ns,
            max_delay_factor=8.0,
        )
        support = support_indices(pdp, threshold=0.99)
        delays = [0, int(case["cdd_delay"])]
        phi = build_phi(grid, delays, support)
        p_support = pdp[support]
        weighted_phi = phi * np.sqrt(np.tile(p_support, len(delays)))[None, :]
        obs_cov = weighted_phi @ weighted_phi.conj().T
        noise_var_ls_at_4db = (10.0 ** (-4.0 / 10.0)) / 2.0
        reg_obs_cov = obs_cov + noise_var_ls_at_4db * np.eye(obs_cov.shape[0], dtype=np.complex128)

        row = {
            "case_id": case["case_id"],
            "delay_spread_ns": float(case["delay_spread_ns"]),
            "cdd_delay_samples": int(case["cdd_delay"]),
            "cdd_delay_ns": float(case["cdd_delay"]) * sample_period_ns,
            "dmrs_spacing_sc": int(case["dmrs_spacing_sc"]),
            "n_pilots_freq": int(grid.pilot_count),
            "support_len_99pct": int(len(support)),
            "n_unknown_taps": int(phi.shape[1]),
            "sample_period_ns": float(sample_period_ns),
        }
        for prefix, mat in (
            ("phi", phi),
            ("weighted_phi", weighted_phi),
            ("obs_cov", obs_cov),
            ("reg_obs_cov_4db", reg_obs_cov),
        ):
            metrics = svd_metrics(mat)
            for k, v in metrics.items():
                row[f"{prefix}_{k}"] = v
        cond_rows.append(row)

        max_delta = grid.n_sc - 1
        pdp_cdd = shifted_pdp(pdp, delays)
        orig_label = f"{case['case_id']}_orig_branch"
        cdd_label = f"{case['case_id']}_cdd_equiv"
        orig_row = coherence_row(orig_label, pdp, grid, int(resource.scs_khz), max_delta)
        cdd_row = coherence_row(cdd_label, pdp_cdd, grid, int(resource.scs_khz), max_delta)
        for r in (orig_row, cdd_row):
            r.update({
                "case_id": case["case_id"],
                "delay_spread_ns": float(case["delay_spread_ns"]),
                "cdd_delay_samples": int(case["cdd_delay"]),
                "dmrs_spacing_sc": int(case["dmrs_spacing_sc"]),
                "allocation_sc": int(grid.n_sc),
                "allocation_mhz": float(grid.n_sc * resource.scs_khz / 1000.0),
            })
            coh_rows.append(r)
        curves_by_case[case["case_id"]] = {
            "orig branch": freq_corr_from_pdp(pdp, grid.n_fft, max_delta),
            "CDD equiv": freq_corr_from_pdp(pdp_cdd, grid.n_fft, max_delta),
        }

    write_csv(cond_rows, out_dir / "reconstruction_conditioning.csv")
    write_csv(coh_rows, out_dir / "coherence_bandwidth.csv")

    for case_id, curves in curves_by_case.items():
        plot_corr(curves, 30, out_dir / f"{case_id}_covariance_correlation.png")
    plot_corr_summary(
        curves_by_case,
        30,
        fig_dir / "experiment8_covariance_correlation_with_bc.png",
    )
    plot_original_pdp(
        [float(case["delay_spread_ns"]) for case in cases],
        sample_period_ns=1e9 / (4096.0 * 30.0 * 1e3),
        path=fig_dir / "experiment8_original_pdp.png",
    )

    print(f"Wrote diagnostics to {out_dir}")


if __name__ == "__main__":
    main()
