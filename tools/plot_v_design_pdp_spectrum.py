from __future__ import annotations

import argparse
import csv
import html
import json
import math
from datetime import datetime
from pathlib import Path
import sys
from typing import Dict, List, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ResourceConfig
from cdd_lls.phy.channel_tdl import make_exponential_pdp
from cdd_lls.phy.estimators import covariance_matrix
from cdd_lls.phy.resource_grid import build_resource_grid
from tools.run_v_design_piecewise_tradeoff import (
    VDesign,
    build_designs,
    cdd_design,
    effective_covariance_from_design,
    parse_int_list,
)


def write_csv(rows: List[Dict[str, object]], path: Path) -> None:
    if not rows:
        return
    fields = sorted(set().union(*(row.keys() for row in rows)))
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def find_design(designs: Sequence[VDesign], label: str) -> VDesign:
    for design in designs:
        if design.label == label:
            return design
    raise ValueError(f"Design not found: {label}")


def full_delay_power_from_frequency_cov(
    Rf: np.ndarray,
    active_subcarriers: np.ndarray,
    n_fft: int,
) -> Tuple[np.ndarray, float]:
    """Return diag(F^H Rf F) on the full nFFT delay grid.

    The full delay covariance is not explicitly formed. Its diagonal is the
    nFFT-point IDFT of frequency-covariance lag sums. The off-diagonal energy
    ratio is computed from Frobenius norms using the unitary full DFT.
    """
    Rf = np.asarray(Rf, dtype=np.complex128)
    active = np.asarray(active_subcarriers, dtype=np.int64)
    lag_sums = np.zeros(int(n_fft), dtype=np.complex128)
    for ia, ka in enumerate(active):
        diffs = (int(ka) - active) % int(n_fft)
        np.add.at(lag_sums, diffs, Rf[ia, :])
    diag_delay = np.fft.ifft(lag_sums)
    power = np.maximum(np.real(diag_delay), 0.0)
    total_frob = float(np.sum(np.abs(Rf) ** 2))
    diag_energy = float(np.sum(np.abs(diag_delay) ** 2))
    offdiag_ratio = float(max(total_frob - diag_energy, 0.0) / max(total_frob, 1e-30))
    return power, offdiag_ratio


def spectrum_rows_for_design(
    *,
    label: str,
    display_label: str,
    family: str,
    Rf: np.ndarray,
    active_subcarriers: np.ndarray,
    n_fft: int,
    sample_period_ns: float,
    min_signed_bin: int,
    max_signed_bin: int,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    power, offdiag_ratio = full_delay_power_from_frequency_cov(Rf, active_subcarriers, n_fft)
    total = float(np.sum(power))
    if total > 0.0:
        power_norm = power / total
    else:
        power_norm = power
    cdf = np.cumsum(power_norm)
    e90_idx = int(np.searchsorted(cdf, 0.90, side="left")) if power_norm.size else 0
    e95_idx = int(np.searchsorted(cdf, 0.95, side="left")) if power_norm.size else 0
    e99_idx = int(np.searchsorted(cdf, 0.99, side="left")) if power_norm.size else 0
    peak_idx = int(np.argmax(power_norm)) if power_norm.size else 0
    tail_64 = float(np.sum(power_norm[64:])) if power_norm.size > 64 else 0.0
    tail_128 = float(np.sum(power_norm[128:])) if power_norm.size > 128 else 0.0
    rows: List[Dict[str, object]] = []
    for tau in range(int(power_norm.size)):
        signed_tau = tau if tau <= int(n_fft) // 2 else tau - int(n_fft)
        if signed_tau < int(min_signed_bin) or signed_tau > int(max_signed_bin):
            continue
        rows.append({
            "label": display_label,
            "design_label": label,
            "family": family,
            "delay_bin": int(tau),
            "signed_delay_bin": int(signed_tau),
            "signed_delay_ns": float(signed_tau * sample_period_ns),
            "delay_ns": float(tau * sample_period_ns),
            "power_norm": float(power_norm[tau]),
            "power_db": float(10.0 * np.log10(max(power_norm[tau], 1e-16))),
            "offdiag_energy_ratio": offdiag_ratio,
        })
    summary = {
        "label": display_label,
        "design_label": label,
        "family": family,
        "peak_delay_bin": peak_idx,
        "peak_delay_ns": float(peak_idx * sample_period_ns),
        "pdp_e90_delay_bin": e90_idx,
        "pdp_e90_delay_ns": float(e90_idx * sample_period_ns),
        "pdp_e95_delay_bin": e95_idx,
        "pdp_e95_delay_ns": float(e95_idx * sample_period_ns),
        "pdp_e99_delay_bin": e99_idx,
        "pdp_e99_delay_ns": float(e99_idx * sample_period_ns),
        "tail_power_ge_64_bins": tail_64,
        "tail_power_ge_128_bins": tail_128,
        "offdiag_energy_ratio": offdiag_ratio,
        "total_power": total,
    }
    return rows, summary


def spectrum_rows_from_pdp(
    *,
    display_label: str,
    pdp: np.ndarray,
    sample_period_ns: float,
    min_signed_bin: int,
    max_signed_bin: int,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    power = np.asarray(pdp, dtype=np.float64)
    power_norm = power / max(float(np.sum(power)), 1e-30)
    rows: List[Dict[str, object]] = []
    for tau, value in enumerate(power_norm):
        if tau < int(min_signed_bin) or tau > int(max_signed_bin):
            continue
        rows.append({
            "label": display_label,
            "design_label": "input_exponential_pdp",
            "family": "physical_channel",
            "delay_bin": int(tau),
            "signed_delay_bin": int(tau),
            "signed_delay_ns": float(tau * sample_period_ns),
            "delay_ns": float(tau * sample_period_ns),
            "power_norm": float(value),
            "power_db": float(10.0 * np.log10(max(value, 1e-16))),
            "offdiag_energy_ratio": 0.0,
        })
    cdf = np.cumsum(power_norm)
    e90_idx = int(np.searchsorted(cdf, 0.90, side="left")) if power_norm.size else 0
    e95_idx = int(np.searchsorted(cdf, 0.95, side="left")) if power_norm.size else 0
    e99_idx = int(np.searchsorted(cdf, 0.99, side="left")) if power_norm.size else 0
    peak_idx = int(np.argmax(power_norm)) if power_norm.size else 0
    summary = {
        "label": display_label,
        "design_label": "input_exponential_pdp",
        "family": "physical_channel",
        "peak_delay_bin": peak_idx,
        "peak_delay_ns": float(peak_idx * sample_period_ns),
        "pdp_e90_delay_bin": e90_idx,
        "pdp_e90_delay_ns": float(e90_idx * sample_period_ns),
        "pdp_e95_delay_bin": e95_idx,
        "pdp_e95_delay_ns": float(e95_idx * sample_period_ns),
        "pdp_e99_delay_bin": e99_idx,
        "pdp_e99_delay_ns": float(e99_idx * sample_period_ns),
        "tail_power_ge_64_bins": float(np.sum(power_norm[64:])) if power_norm.size > 64 else 0.0,
        "tail_power_ge_128_bins": float(np.sum(power_norm[128:])) if power_norm.size > 128 else 0.0,
        "offdiag_energy_ratio": 0.0,
        "total_power": float(np.sum(power)),
    }
    return rows, summary


def svg_pdp_plot(rows: Sequence[Dict[str, object]], path: Path) -> None:
    groups: Dict[str, List[Dict[str, object]]] = {}
    for row in rows:
        groups.setdefault(str(row["label"]), []).append(row)
    width, height = 1120, 680
    left, right, top, bottom = 84, 190, 54, 78
    plot_w = width - left - right
    plot_h = height - top - bottom
    xs = [float(r["signed_delay_bin"]) for r in rows]
    ys = [float(r["power_db"]) for r in rows]
    x0, x1 = min(xs), max(xs)
    y0, y1 = max(min(ys), -90.0), min(max(ys), 0.0)
    y0 = min(y0, -60.0)
    y1 = 0.0

    def sx(x: float) -> float:
        return left + (x - x0) / max(x1 - x0, 1e-12) * plot_w

    def sy(y: float) -> float:
        yy = max(min(y, y1), y0)
        return top + plot_h - (yy - y0) / max(y1 - y0, 1e-12) * plot_h

    colors = {
        "Original input PDP": "#222222",
        "CDD C6 step=7": "#26547c",
        "PW P6 cyclic slope=6": "#06a77d",
        "PW P8 cyclic slope=16": "#7b2cbf",
        "PW reset R1 slope=2": "#d62828",
    }
    fallback = ["#ef476f", "#f77f00", "#118ab2", "#6a994e"]
    grid = []
    for i in range(7):
        x = left + plot_w * i / 6
        val = x0 + (x1 - x0) * i / 6
        grid.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top+plot_h}" stroke="#e6e8eb"/>')
        grid.append(f'<text x="{x:.2f}" y="{top+plot_h+24}" font-size="12" text-anchor="middle">{val:.0f}</text>')
    for i in range(7):
        y = top + plot_h - plot_h * i / 6
        val = y0 + (y1 - y0) * i / 6
        grid.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}" stroke="#e6e8eb"/>')
        grid.append(f'<text x="{left-8}" y="{y+4:.2f}" font-size="12" text-anchor="end">{val:.0f}</text>')

    elements = []
    legend = []
    for idx, (label, pts) in enumerate(groups.items()):
        pts = sorted(pts, key=lambda r: int(r["signed_delay_bin"]))
        color = colors.get(label, fallback[idx % len(fallback)])
        coords = [
            f'{sx(float(r["signed_delay_bin"])):.2f},{sy(float(r["power_db"])):.2f}'
            for r in pts
        ]
        elements.append(f'<polyline points="{" ".join(coords)}" fill="none" stroke="{color}" stroke-width="2"/>')
        lx = left + plot_w + 22
        ly = top + 20 + idx * 22
        legend.append(f'<line x1="{lx}" y1="{ly:.2f}" x2="{lx+18}" y2="{ly:.2f}" stroke="{color}" stroke-width="2.5"/>')
        legend.append(f'<text x="{lx+24}" y="{ly+4:.2f}" font-size="12">{html.escape(label)}</text>')

    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="white"/>
<text x="{width/2}" y="30" text-anchor="middle" font-size="18" font-family="Arial, sans-serif">Active-band effective PDP spectrum</text>
<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#333" stroke-width="1.2"/>
<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#333" stroke-width="1.2"/>
{''.join(grid)}
{''.join(elements)}
{''.join(legend)}
<text x="{left+plot_w/2}" y="{height-24}" text-anchor="middle" font-size="14">Signed delay bin on nFFT grid (8.138 ns/bin)</text>
<text x="22" y="{top+plot_h/2}" text-anchor="middle" font-size="14" transform="rotate(-90 22 {top+plot_h/2})">Normalized delay power (dB)</text>
</svg>
'''
    path.write_text(svg, encoding="utf-8")


def run(args: argparse.Namespace) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out) / f"v_design_pdp_spectrum_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    resource = ResourceConfig(
        carrier_bandwidth_mhz=100.0,
        scs_khz=30,
        n_fft=int(args.n_fft),
        n_prbs=int(args.n_prbs),
        pdsch_n_symbols=10,
        dmrs_symbol_indices=[2, 7],
        dmrs_spacing_sc=int(args.dmrs_spacing_sc),
        prg_size_rb=4,
    )
    grid = build_resource_grid(resource)
    sample_period_ns = 1e9 / (float(resource.n_fft) * float(resource.scs_khz) * 1e3)
    pdp = make_exponential_pdp(
        delay_spread_ns=float(args.delay_spread_ns),
        sample_period_ns=sample_period_ns,
        max_delay_factor=float(args.max_delay_factor),
    )
    base_cov = covariance_matrix(grid.subcarrier_indices, grid.subcarrier_indices, pdp, grid.n_fft)
    designs = build_designs(
        grid=grid,
        n_tx=int(args.n_tx),
        n_segments=int(args.n_segments),
        cdd_steps=parse_int_list(args.cdd_steps),
        piecewise_slopes=parse_int_list(args.piecewise_slopes),
        random_per_slope=int(args.random_per_slope),
        seed=int(args.seed),
    )
    selected: List[Tuple[str, VDesign]] = [
        ("CDD C6 step=7", cdd_design(grid, int(args.n_tx), 7)),
        ("PW P6 cyclic slope=6", find_design(designs, "PW cont cyclic slope=6")),
        ("PW P8 cyclic slope=16", find_design(designs, "PW cont cyclic slope=16")),
        ("PW reset R1 slope=2", find_design(designs, "PW reset DFTseg slope=2")),
    ]
    rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []
    pdp_rows, pdp_summary = spectrum_rows_from_pdp(
        display_label="Original input PDP",
        pdp=pdp,
        sample_period_ns=sample_period_ns,
        min_signed_bin=-int(args.max_negative_bins),
        max_signed_bin=int(args.max_positive_bins),
    )
    rows.extend(pdp_rows)
    summary_rows.append(pdp_summary)
    for display_label, design in selected:
        Rf = effective_covariance_from_design(base_cov, design)
        part_rows, summary = spectrum_rows_for_design(
            label=design.label,
            display_label=display_label,
            family=design.family,
            Rf=Rf,
            active_subcarriers=grid.subcarrier_indices,
            n_fft=grid.n_fft,
            sample_period_ns=sample_period_ns,
            min_signed_bin=-int(args.max_negative_bins),
            max_signed_bin=int(args.max_positive_bins),
        )
        rows.extend(part_rows)
        summary_rows.append(summary)
    write_csv(rows, out_dir / "effective_pdp_spectrum_selected.csv")
    write_csv(summary_rows, out_dir / "effective_pdp_spectrum_summary.csv")
    svg_pdp_plot(rows, fig_dir / "experiment18_v_design_effective_pdp_spectrum.svg")
    with (out_dir / "run_config.json").open("w", encoding="utf-8") as f:
        json.dump(vars(args), f, indent=2, ensure_ascii=False)
    print(out_dir)
    return out_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/v_design_pdp_spectrum")
    parser.add_argument("--fig-dir", default="docs/figures")
    parser.add_argument("--n-tx", type=int, default=8)
    parser.add_argument("--n-segments", type=int, default=8)
    parser.add_argument("--n-prbs", type=int, default=48)
    parser.add_argument("--n-fft", type=int, default=4096)
    parser.add_argument("--dmrs-spacing-sc", type=int, default=24)
    parser.add_argument("--delay-spread-ns", type=float, default=10.0)
    parser.add_argument("--max-delay-factor", type=float, default=8.0)
    parser.add_argument("--cdd-steps", default="1,2,3,4,5,6,7,8,10,12,16,24,32,48,64,96,128,192,256")
    parser.add_argument("--piecewise-slopes", default="0,1,2,4,6,8,12,16,24,32,48,64")
    parser.add_argument("--random-per-slope", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20260615)
    parser.add_argument("--max-positive-bins", type=int, default=220)
    parser.add_argument("--max-negative-bins", type=int, default=220)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
