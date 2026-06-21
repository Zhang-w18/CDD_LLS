from __future__ import annotations

import argparse
import csv
import json
import math
import os
from datetime import datetime
from pathlib import Path
import sys
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ChannelConfig, ResourceConfig
from cdd_lls.core.mcs import build_tb_layout, get_mcs
from cdd_lls.phy.channel_tdl import generate_tdl_channel, make_exponential_pdp
from cdd_lls.phy.estimators import covariance_matrix, nmse
from cdd_lls.phy.ldpc import SionnaLDPCAdapter
from cdd_lls.phy.precoding import equivalent_channel
from cdd_lls.phy.qam import qam_demapper_maxlog, qam_modulate
from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
from tools.run_v_design_piecewise_tradeoff import (
    VDesign,
    build_designs,
    decode_same_tb_batch,
    equalize_mrc,
    finite_float,
    make_alg1_estimator,
    parse_float_list,
    stable_seed,
    write_csv,
)


def read_pareto_ids(path: Path, include_families: Sequence[str]) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    allowed = set(include_families)
    out = [r for r in rows if r.get("family") in allowed and str(r.get("pareto_id", "")).strip()]
    order = {"cdd_linear_phase": 0, "piecewise_linear_continuous": 1}
    return sorted(out, key=lambda r: (order.get(r["family"], 99), str(r["pareto_id"])[0], int(str(r["pareto_id"])[1:])))


def select_designs(all_designs: Sequence[VDesign], pareto_rows: Sequence[Dict[str, str]]) -> List[Tuple[str, VDesign, Dict[str, str]]]:
    by_label = {d.label: d for d in all_designs}
    selected: List[Tuple[str, VDesign, Dict[str, str]]] = []
    missing = []
    for row in pareto_rows:
        label = str(row["precoder"])
        design = by_label.get(label)
        if design is None:
            missing.append(label)
            continue
        selected.append((str(row["pareto_id"]), design, row))
    if missing:
        raise ValueError(f"Pareto designs not found in generated design set: {missing}")
    return selected


def mean_or_nan(values: Iterable[float]) -> float:
    vals = [float(x) for x in values if np.isfinite(float(x))]
    return float(np.mean(vals)) if vals else float("nan")


def interp_target_snr(rows: Sequence[Dict[str, object]], target: float) -> float:
    pts = sorted(
        [(float(r["snr_db"]), float(r["bler"])) for r in rows if np.isfinite(float(r["bler"]))],
        key=lambda x: x[0],
    )
    if not pts:
        return float("nan")
    for (s0, b0), (s1, b1) in zip(pts[:-1], pts[1:]):
        if (b0 - target) == 0.0:
            return float(s0)
        if (b0 - target) * (b1 - target) <= 0.0 and b0 != b1:
            t = (target - b0) / (b1 - b0)
            return float(s0 + t * (s1 - s0))
    return float("nan")


def compute_target_rows(rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    target_rows = []
    for pid in sorted({str(r["pareto_id"]) for r in rows}, key=lambda x: (x[0], int(x[1:]))):
        pts = [r for r in rows if str(r["pareto_id"]) == pid]
        first = pts[0]
        target_rows.append({
            "pareto_id": pid,
            "precoder": first["precoder"],
            "family": first["family"],
            "snr_at_bler_0p1_db": interp_target_snr(pts, 0.10),
            "snr_at_bler_0p01_db": interp_target_snr(pts, 0.01),
            "min_observed_bler": min(float(r["bler"]) for r in pts),
            "coherence_bw_abs_0p5_sc": first["coherence_bw_abs_0p5_sc"],
            "diversity_log10_product_norm": first["diversity_log10_product_norm"],
        })
    return target_rows


def finite_or_nan(value: object) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return out if np.isfinite(out) else float("nan")


def best_target(target_rows: Sequence[Dict[str, object]], family: str, field: str) -> Tuple[float, str, str]:
    best = (float("nan"), "", "")
    for row in target_rows:
        if str(row.get("family")) != family:
            continue
        snr = finite_or_nan(row.get(field))
        if not np.isfinite(snr):
            continue
        if not np.isfinite(best[0]) or snr < best[0]:
            best = (snr, str(row["pareto_id"]), str(row["precoder"]))
    return best


def compute_piecewise_gain_rows(target_rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    cdd_10 = best_target(target_rows, "cdd_linear_phase", "snr_at_bler_0p1_db")
    cdd_1 = best_target(target_rows, "cdd_linear_phase", "snr_at_bler_0p01_db")
    gain_rows = []
    for row in target_rows:
        if str(row.get("family")) != "piecewise_linear_continuous":
            continue
        p10 = finite_or_nan(row.get("snr_at_bler_0p1_db"))
        p1 = finite_or_nan(row.get("snr_at_bler_0p01_db"))
        gain10 = cdd_10[0] - p10 if np.isfinite(cdd_10[0]) and np.isfinite(p10) else float("nan")
        gain1 = cdd_1[0] - p1 if np.isfinite(cdd_1[0]) and np.isfinite(p1) else float("nan")
        gain_rows.append({
            "pareto_id": row["pareto_id"],
            "precoder": row["precoder"],
            "coherence_bw_abs_0p5_sc": row["coherence_bw_abs_0p5_sc"],
            "diversity_log10_product_norm": row["diversity_log10_product_norm"],
            "piecewise_snr_at_bler_0p1_db": p10,
            "best_cdd_pareto_id_at_0p1": cdd_10[1],
            "best_cdd_precoder_at_0p1": cdd_10[2],
            "best_cdd_snr_at_bler_0p1_db": cdd_10[0],
            "snr_gain_vs_best_cdd_at_bler_0p1_db": gain10,
            "piecewise_snr_at_bler_0p01_db": p1,
            "best_cdd_pareto_id_at_0p01": cdd_1[1],
            "best_cdd_precoder_at_0p01": cdd_1[2],
            "best_cdd_snr_at_bler_0p01_db": cdd_1[0],
            "snr_gain_vs_best_cdd_at_bler_0p01_db": gain1,
            "min_observed_bler": row["min_observed_bler"],
        })
    return gain_rows


def write_analysis_outputs(rows: Sequence[Dict[str, object]], out_dir: Path, fig_dir: Path, fig_prefix: str = "experiment19") -> None:
    write_csv(list(rows), out_dir / "pareto_link_curves.csv")
    target_rows = compute_target_rows(rows)
    write_csv(target_rows, out_dir / "pareto_link_targets.csv")
    write_csv(compute_piecewise_gain_rows(target_rows), out_dir / "pareto_link_piecewise_gain_vs_best_cdd.csv")
    svg_line_plot(rows, metric="ce_nmse_mean", family="cdd_linear_phase", title="CDD Pareto Alg1 NMSE", ylabel="NMSE", path=fig_dir / f"{fig_prefix}_cdd_pareto_nmse.svg", log_y=True)
    svg_line_plot(rows, metric="ce_nmse_mean", family="piecewise_linear_continuous", title="Piecewise-linear Pareto Alg1 NMSE", ylabel="NMSE", path=fig_dir / f"{fig_prefix}_piecewise_pareto_nmse.svg", log_y=True)
    svg_line_plot(rows, metric="bler", family="cdd_linear_phase", title="CDD Pareto BLER", ylabel="BLER", path=fig_dir / f"{fig_prefix}_cdd_pareto_bler.svg", log_y=True, threshold_lines=(0.1, 0.01))
    svg_line_plot(rows, metric="bler", family="piecewise_linear_continuous", title="Piecewise-linear Pareto BLER", ylabel="BLER", path=fig_dir / f"{fig_prefix}_piecewise_pareto_bler.svg", log_y=True, threshold_lines=(0.1, 0.01))
    svg_all_pareto_plot(rows, metric="ce_nmse_mean", title="All Pareto Alg1 NMSE", ylabel="NMSE", path=fig_dir / f"{fig_prefix}_all_pareto_nmse.svg", log_y=True)
    svg_all_pareto_plot(rows, metric="bler", title="All Pareto BLER", ylabel="BLER", path=fig_dir / f"{fig_prefix}_all_pareto_bler.svg", log_y=True, threshold_lines=(0.1, 0.01))


def svg_line_plot(
    rows: Sequence[Dict[str, object]],
    *,
    metric: str,
    family: str,
    title: str,
    ylabel: str,
    path: Path,
    log_y: bool,
    threshold_lines: Sequence[float] = (),
) -> None:
    fam_rows = [r for r in rows if str(r["family"]) == family]
    if not fam_rows:
        return
    groups: Dict[str, List[Dict[str, object]]] = {}
    for row in fam_rows:
        groups.setdefault(str(row["pareto_id"]), []).append(row)
    width, height = 1060, 660
    left, right, top, bottom = 78, 126, 54, 78
    plot_w = width - left - right
    plot_h = height - top - bottom
    xs = [float(r["snr_db"]) for r in fam_rows]
    y_floor = 1e-4 if metric == "bler" else 1e-8
    raw_ys = [max(float(r[metric]), y_floor) for r in fam_rows]
    ys = [math.log10(y) for y in raw_ys] if log_y else raw_ys
    if threshold_lines:
        ys.extend(math.log10(x) if log_y else x for x in threshold_lines)
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys), max(ys)
    if abs(x1 - x0) < 1e-12:
        x0 -= 1.0
        x1 += 1.0
    if abs(y1 - y0) < 1e-12:
        y0 -= 0.5
        y1 += 0.5
    xpad = 0.04 * (x1 - x0)
    ypad = 0.08 * (y1 - y0)
    x0 -= xpad
    x1 += xpad
    y0 -= ypad
    y1 += ypad

    def sx(x: float) -> float:
        return left + (x - x0) / (x1 - x0) * plot_w

    def sy(y: float) -> float:
        yy = math.log10(max(y, y_floor)) if log_y else y
        return top + plot_h - (yy - y0) / (y1 - y0) * plot_h

    colors = [
        "#26547c", "#ef476f", "#06a77d", "#7b2cbf", "#f77f00", "#118ab2",
        "#9d4edd", "#2a9d8f", "#d62828", "#6a994e", "#003049", "#bc6c25",
    ]
    grid = []
    unique_xs = sorted(set(xs))
    x_ticks = unique_xs if len(unique_xs) <= 8 else [x0 + (x1 - x0) * i / 5 for i in range(6)]
    for val in x_ticks:
        x = sx(val)
        grid.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top+plot_h}" stroke="#e6e8eb"/>')
        grid.append(f'<text x="{x:.2f}" y="{top+plot_h+24}" font-size="12" text-anchor="middle">{val:.3g}</text>')
    if log_y:
        major_exponents = range(math.ceil(y0), math.floor(y1) + 1)
        for exponent in major_exponents:
            y = top + plot_h - (exponent - y0) / max(y1 - y0, 1e-12) * plot_h
            grid.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}" stroke="#d9dde2"/>')
            grid.append(f'<text x="{left-8}" y="{y+4:.2f}" font-size="12" text-anchor="end">1e{exponent}</text>')
        for exponent in range(math.floor(y0), math.ceil(y1) + 1):
            for multiplier in (2, 5):
                value = math.log10(multiplier) + exponent
                if y0 < value < y1:
                    y = top + plot_h - (value - y0) / max(y1 - y0, 1e-12) * plot_h
                    grid.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}" stroke="#eef0f2"/>')
    else:
        for i in range(6):
            y = top + plot_h - plot_h * i / 5
            val = y0 + (y1 - y0) * i / 5
            grid.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}" stroke="#e6e8eb"/>')
            grid.append(f'<text x="{left-8}" y="{y+4:.2f}" font-size="12" text-anchor="end">{val:.3g}</text>')
    elements = []
    for target in threshold_lines:
        y = sy(float(target))
        elements.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}" stroke="#111" stroke-width="1" stroke-dasharray="5 4"/>')
        elements.append(f'<text x="{left+plot_w+6}" y="{y+4:.2f}" font-size="11" fill="#111">{target:g}</text>')
    legend = []
    for idx, (pid, pts) in enumerate(sorted(groups.items(), key=lambda kv: int(kv[0][1:]))):
        pts = sorted(pts, key=lambda r: float(r["snr_db"]))
        color = colors[idx % len(colors)]
        coords = []
        for row in pts:
            x = sx(float(row["snr_db"]))
            y = sy(float(row[metric]))
            coords.append(f"{x:.2f},{y:.2f}")
            elements.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.6" fill="{color}"><title>{pid} {row["precoder"]}</title></circle>')
        points_attr = " ".join(coords)
        elements.append(f'<polyline points="{points_attr}" fill="none" stroke="{color}" stroke-width="1.7"/>')
        last = pts[-1]
        elements.append(f'<text x="{sx(float(last["snr_db"]))+5:.2f}" y="{sy(float(last[metric]))+4:.2f}" font-size="10" fill="{color}" font-weight="700">{pid}</text>')
        lx = left + plot_w + 18
        ly = top + 18 + idx * 18
        legend.append(f'<line x1="{lx}" y1="{ly:.2f}" x2="{lx+16}" y2="{ly:.2f}" stroke="{color}" stroke-width="2"/>')
        legend.append(f'<text x="{lx+20}" y="{ly+4:.2f}" font-size="11">{pid}</text>')
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="white"/>
<text x="{width/2}" y="28" text-anchor="middle" font-size="18" font-family="Arial, sans-serif">{title}</text>
<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#333" stroke-width="1.2"/>
<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#333" stroke-width="1.2"/>
{''.join(grid)}
{''.join(elements)}
{''.join(legend)}
<text x="{left+plot_w/2}" y="{height-24}" text-anchor="middle" font-size="14">SNR (dB)</text>
<text x="20" y="{top+plot_h/2}" text-anchor="middle" font-size="14" transform="rotate(-90 20 {top+plot_h/2})">{ylabel}</text>
</svg>
'''
    path.write_text(svg, encoding="utf-8")


def svg_all_pareto_plot(
    rows: Sequence[Dict[str, object]],
    *,
    metric: str,
    title: str,
    ylabel: str,
    path: Path,
    log_y: bool,
    threshold_lines: Sequence[float] = (),
    label_endpoints: bool = True,
    direct_labels: bool = False,
) -> None:
    families = {"cdd_linear_phase", "piecewise_linear_continuous"}
    plot_rows = [r for r in rows if str(r["family"]) in families]
    if not plot_rows:
        return
    groups: Dict[str, List[Dict[str, object]]] = {}
    for row in plot_rows:
        groups.setdefault(str(row["pareto_id"]), []).append(row)
    width, height = 1240, 720
    left, right, top, bottom = 84, 250, 56, 82
    plot_w = width - left - right
    plot_h = height - top - bottom
    xs = [float(r["snr_db"]) for r in plot_rows]
    y_floor = 1e-8

    def display_value(row: Dict[str, object]) -> float:
        value = float(row[metric])
        if metric.startswith("bler") and value <= 0.0:
            return 0.5 / max(float(row.get("trials", 1)), 1.0)
        return max(value, y_floor)

    raw_ys = [display_value(r) for r in plot_rows]
    ys = [math.log10(y) for y in raw_ys] if log_y else raw_ys
    if threshold_lines:
        ys.extend(math.log10(x) if log_y else x for x in threshold_lines)
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys), max(ys)
    xpad = 0.04 * (x1 - x0)
    ypad = 0.08 * (y1 - y0)
    x0 -= xpad
    x1 += xpad
    y0 -= ypad
    y1 += ypad

    def sx(x: float) -> float:
        return left + (x - x0) / max(x1 - x0, 1e-12) * plot_w

    def sy(y: float) -> float:
        yy = math.log10(max(y, y_floor)) if log_y else y
        return top + plot_h - (yy - y0) / max(y1 - y0, 1e-12) * plot_h

    cdd_colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#6F4E7C", "#222222"]
    pw_colors = ["#A0006D", "#7A3E00", "#117733", "#332288", "#AA4499", "#44AA99", "#882255", "#999933", "#661100", "#004488", "#EE7733", "#228833"]
    cdd_markers = ["circle", "square", "triangle", "diamond", "down_triangle", "cross", "x", "hexagon"]
    pw_markers = ["circle", "square", "triangle", "diamond"]
    family_order = {"cdd_linear_phase": 0, "piecewise_linear_continuous": 1}

    def marker_svg(marker: str, x: float, y: float, color: str, title_text: str, filled: bool) -> str:
        fill = color if filled else "white"
        common = f'stroke="{color}" stroke-width="1.6" fill="{fill}"'
        if marker == "circle":
            shape = f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4" {common}/>'
        elif marker == "square":
            shape = f'<rect x="{x-3.8:.2f}" y="{y-3.8:.2f}" width="7.6" height="7.6" {common}/>'
        elif marker == "triangle":
            shape = f'<polygon points="{x:.2f},{y-4.6:.2f} {x-4.3:.2f},{y+3.6:.2f} {x+4.3:.2f},{y+3.6:.2f}" {common}/>'
        elif marker == "down_triangle":
            shape = f'<polygon points="{x-4.3:.2f},{y-3.6:.2f} {x+4.3:.2f},{y-3.6:.2f} {x:.2f},{y+4.6:.2f}" {common}/>'
        elif marker == "diamond":
            shape = f'<polygon points="{x:.2f},{y-4.5:.2f} {x-4.5:.2f},{y:.2f} {x:.2f},{y+4.5:.2f} {x+4.5:.2f},{y:.2f}" {common}/>'
        elif marker == "hexagon":
            shape = f'<polygon points="{x-4.2:.2f},{y:.2f} {x-2.1:.2f},{y-3.7:.2f} {x+2.1:.2f},{y-3.7:.2f} {x+4.2:.2f},{y:.2f} {x+2.1:.2f},{y+3.7:.2f} {x-2.1:.2f},{y+3.7:.2f}" {common}/>'
        elif marker == "cross":
            shape = f'<path d="M {x-4:.2f} {y:.2f} H {x+4:.2f} M {x:.2f} {y-4:.2f} V {y+4:.2f}" stroke="{color}" stroke-width="2" fill="none"/>'
        else:
            shape = f'<path d="M {x-3.5:.2f} {y-3.5:.2f} L {x+3.5:.2f} {y+3.5:.2f} M {x+3.5:.2f} {y-3.5:.2f} L {x-3.5:.2f} {y+3.5:.2f}" stroke="{color}" stroke-width="2" fill="none"/>'
        return f'<g><title>{title_text}</title>{shape}</g>'

    def sort_key(pid: str) -> Tuple[int, int]:
        fam = str(groups[pid][0]["family"])
        return family_order.get(fam, 99), int(pid[1:])

    grid = []
    unique_xs = sorted(set(xs))
    x_ticks = unique_xs if len(unique_xs) <= 8 else [x0 + (x1 - x0) * i / 5 for i in range(6)]
    for val in x_ticks:
        x = sx(val)
        grid.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top+plot_h}" stroke="#e6e8eb"/>')
        grid.append(f'<text x="{x:.2f}" y="{top+plot_h+24}" font-size="12" text-anchor="middle">{val:.3g}</text>')
    if log_y:
        major_exponents = range(math.ceil(y0), math.floor(y1) + 1)
        for exponent in major_exponents:
            y = top + plot_h - (exponent - y0) / max(y1 - y0, 1e-12) * plot_h
            grid.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}" stroke="#d9dde2"/>')
            grid.append(f'<text x="{left-8}" y="{y+4:.2f}" font-size="12" text-anchor="end">1e{exponent}</text>')
        for exponent in range(math.floor(y0), math.ceil(y1) + 1):
            for multiplier in (2, 5):
                value = math.log10(multiplier) + exponent
                if y0 < value < y1:
                    y = top + plot_h - (value - y0) / max(y1 - y0, 1e-12) * plot_h
                    grid.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}" stroke="#eef0f2"/>')
    else:
        for i in range(6):
            y = top + plot_h - plot_h * i / 5
            val = y0 + (y1 - y0) * i / 5
            grid.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}" stroke="#e6e8eb"/>')
            grid.append(f'<text x="{left-8}" y="{y+4:.2f}" font-size="12" text-anchor="end">{val:.3g}</text>')
    elements = []
    for target in threshold_lines:
        y = sy(float(target))
        elements.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}" stroke="#111" stroke-width="1" stroke-dasharray="5 4"/>')
        elements.append(f'<text x="{left+plot_w+8}" y="{y+4:.2f}" font-size="11" fill="#111">{target:g}</text>')
    legend = []
    direct_label_specs = []
    cdd_idx = 0
    pw_idx = 0
    for idx, pid in enumerate(sorted(groups, key=sort_key)):
        pts = sorted(groups[pid], key=lambda r: float(r["snr_db"]))
        fam = str(pts[0]["family"])
        numeric_id = int(pid[1:])
        if fam == "cdd_linear_phase":
            color = cdd_colors[numeric_id % len(cdd_colors)]
            marker = cdd_markers[numeric_id % len(cdd_markers)]
            dash = "5 3"
            filled = False
            cdd_idx += 1
        else:
            color = pw_colors[pw_idx % len(pw_colors)]
            marker = pw_markers[pw_idx % len(pw_markers)]
            dash = ""
            filled = True
            pw_idx += 1
        coords = []
        for row in pts:
            x = sx(float(row["snr_db"]))
            y = sy(display_value(row))
            coords.append(f"{x:.2f},{y:.2f}")
            elements.append(marker_svg(marker, x, y, color, f'{pid} {row["precoder"]}', filled))
        points_attr = " ".join(coords)
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        elements.append(f'<polyline points="{points_attr}" fill="none" stroke="{color}" stroke-width="2.1"{dash_attr}/>')
        if label_endpoints:
            last = pts[-1]
            elements.append(f'<text x="{sx(float(last["snr_db"]))+5:.2f}" y="{sy(display_value(last))+4:.2f}" font-size="10" fill="{color}" font-weight="700">{pid}</text>')
        if direct_labels:
            target_log = math.log10(0.03) if metric.startswith("bler") else sum(ys) / len(ys)
            anchor = min(
                pts,
                key=lambda row: abs(
                    (math.log10(display_value(row)) if log_y else display_value(row)) - target_log
                ),
            )
            direct_label_specs.append({
                "pid": pid,
                "color": color,
                "x": sx(float(anchor["snr_db"])),
                "y": sy(display_value(anchor)),
            })
        lx = left + plot_w + 26
        ly = top + 18 + idx * 18
        legend.append(f'<line x1="{lx}" y1="{ly:.2f}" x2="{lx+18}" y2="{ly:.2f}" stroke="{color}" stroke-width="2"{dash_attr}/>')
        legend.append(marker_svg(marker, lx+9, ly, color, pid, filled))
        legend.append(f'<text x="{lx+24}" y="{ly+4:.2f}" font-size="11">{pid}</text>')

    if direct_label_specs:
        ordered = sorted(direct_label_specs, key=lambda item: float(item["y"]))
        gap = 15.0
        label_ys = []
        for item in ordered:
            desired = min(max(float(item["y"]), top + 10.0), top + plot_h - 8.0)
            label_ys.append(max(desired, label_ys[-1] + gap) if label_ys else desired)
        overflow = label_ys[-1] - (top + plot_h - 8.0)
        if overflow > 0:
            label_ys = [value - overflow for value in label_ys]
        for item, label_y in zip(ordered, label_ys):
            anchor_x = float(item["x"])
            anchor_y = float(item["y"])
            label_x = min(anchor_x + 58.0, left + plot_w - 30.0)
            end_x = label_x - 6.0
            end_y = label_y - 3.0
            dx = end_x - anchor_x
            dy = end_y - anchor_y
            norm = max(math.hypot(dx, dy), 1e-12)
            ux, uy = dx / norm, dy / norm
            base_x, base_y = end_x - 6.0 * ux, end_y - 6.0 * uy
            perp_x, perp_y = -uy * 3.0, ux * 3.0
            color = str(item["color"])
            elements.append(f'<line x1="{anchor_x:.2f}" y1="{anchor_y:.2f}" x2="{base_x:.2f}" y2="{base_y:.2f}" stroke="{color}" stroke-width="1.2"/>')
            elements.append(f'<polygon points="{end_x:.2f},{end_y:.2f} {base_x+perp_x:.2f},{base_y+perp_y:.2f} {base_x-perp_x:.2f},{base_y-perp_y:.2f}" fill="{color}"/>')
            elements.append(f'<text x="{label_x:.2f}" y="{label_y:.2f}" font-size="11" fill="{color}" font-weight="700" paint-order="stroke" stroke="white" stroke-width="3">{item["pid"]}</text>')
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="white"/>
<text x="{width/2}" y="30" text-anchor="middle" font-size="18" font-family="Arial, sans-serif">{title}</text>
<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#333" stroke-width="1.2"/>
<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#333" stroke-width="1.2"/>
{''.join(grid)}
{''.join(elements)}
{''.join(legend)}
<text x="{left+plot_w/2}" y="{height-26}" text-anchor="middle" font-size="14">SNR (dB)</text>
<text x="22" y="{top+plot_h/2}" text-anchor="middle" font-size="14" transform="rotate(-90 22 {top+plot_h/2})">{ylabel}</text>
<text x="{left+plot_w+26}" y="{height-28}" font-size="11" fill="#555">dashed: CDD, solid: piecewise</text>
</svg>
'''
    path.write_text(svg, encoding="utf-8")


def run(args: argparse.Namespace) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out) / f"v_design_pareto_link_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

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
    sample_period_ns = 1e9 / (float(resource.n_fft) * float(resource.scs_khz) * 1e3)
    pdp = make_exponential_pdp(
        delay_spread_ns=float(args.delay_spread_ns),
        sample_period_ns=sample_period_ns,
        max_delay_factor=float(args.max_delay_factor),
    )
    base_cov = covariance_matrix(grid.subcarrier_indices, grid.subcarrier_indices, pdp, grid.n_fft)
    all_designs = build_designs(
        grid=grid,
        n_tx=int(args.n_tx),
        n_segments=int(args.n_segments),
        cdd_steps=[1,2,3,4,5,6,7,8,10,12,16,24,32,48,64,96,128,192,256],
        piecewise_slopes=[0,1,2,4,6,8,12,16,24,32,48,64],
        random_per_slope=3,
        seed=int(args.design_seed),
    )
    pareto_rows = read_pareto_ids(
        Path(args.pareto_csv),
        include_families=("cdd_linear_phase", "piecewise_linear_continuous"),
    )
    selected = select_designs(all_designs, pareto_rows)
    snrs = parse_float_list(args.snrs)
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    mcs = get_mcs("nr_256qam", int(args.mcs_index), None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(
        tb.cb_k_values,
        tb.cb_e_values,
        num_iter=int(args.ldpc_iterations),
        llr_clip=float(args.llr_clip),
    )
    channel = ChannelConfig(delay_spread_ns=float(args.delay_spread_ns), max_delay_factor=float(args.max_delay_factor))
    rows: List[Dict[str, object]] = []

    for snr_db in snrs:
        noise_var = 10.0 ** (-float(snr_db) / 10.0)
        noise_var_ls = noise_var / float(len(resource.dmrs_symbol_indices))
        estimators = {
            pid: make_alg1_estimator(grid, pdp, design, noise_var_ls, float(args.loading))
            for pid, design, _ in selected
        }
        accum = {
            pid: {
                "tb_errors": 0,
                "cb_errors": 0,
                "goodput_bits": 0,
                "ce_nmse": [],
                "design": design,
                "pareto": row,
            }
            for pid, design, row in selected
        }
        base_seed = stable_seed("pareto-link", snr_db, base=int(args.seed))
        for trial in range(1, int(args.trials) + 1):
            rng = np.random.default_rng(stable_seed(base_seed, trial, base=int(args.seed)))
            realization = generate_tdl_channel(rng, grid, channel, n_tx=int(args.n_tx), n_rx=int(args.n_rx))
            H = realization.H
            payload = [rng.integers(0, 2, size=int(k), dtype=np.int8) for k in tb.cb_k_values]
            coded_cw = np.concatenate(adapter.encode(payload))
            symbols = qam_modulate(coded_cw, int(mcs.qm))
            pilot_noise = (
                rng.normal(size=(int(args.n_rx), len(pilot_local)))
                + 1j * rng.normal(size=(int(args.n_rx), len(pilot_local)))
            ) * math.sqrt(noise_var_ls / 2.0)
            data_noise = (
                rng.normal(size=(int(args.n_rx), len(data_local)))
                + 1j * rng.normal(size=(int(args.n_rx), len(data_local)))
            ) * math.sqrt(noise_var / 2.0)
            llrs = []
            trial_ids = []
            trial_nmse = []
            for pid, design, _ in selected:
                g = equivalent_channel(H, design.C)
                ls_obs = g[:, pilot_local] + pilot_noise
                g_hat = ls_obs @ estimators[pid].matrix.T
                y = g[:, data_local] * symbols[None, :] + data_noise
                z, no_eff = equalize_mrc(y, g_hat[:, data_local], noise_var)
                llrs.append(qam_demapper_maxlog(z, no_eff, int(mcs.qm)))
                trial_ids.append(pid)
                trial_nmse.append(nmse(g, g_hat))
            dec_results = decode_same_tb_batch(adapter, llrs, payload)
            for pid, dec, ce in zip(trial_ids, dec_results, trial_nmse):
                acc = accum[pid]
                acc["tb_errors"] += int(not dec.tb_success)
                acc["cb_errors"] += int(sum(1 for ok in dec.cb_success if not ok))
                acc["goodput_bits"] += int(dec.goodput_bits)
                acc["ce_nmse"].append(float(ce))
            if int(args.progress_every) > 0 and trial % int(args.progress_every) == 0:
                print(f"[progress] SNR={snr_db:g} dB trial={trial}/{int(args.trials)}", flush=True)

        for pid, acc in accum.items():
            design = acc["design"]
            prow = acc["pareto"]
            bler = float(acc["tb_errors"]) / float(args.trials)
            rows.append({
                "pareto_id": pid,
                "precoder": design.label,
                "family": design.family,
                "snr_db": float(snr_db),
                "trials": int(args.trials),
                "tb_errors": int(acc["tb_errors"]),
                "cb_errors": int(acc["cb_errors"]),
                "bler": bler,
                "ce_nmse_mean": mean_or_nan(acc["ce_nmse"]),
                "ce_nmse_median": float(np.median(acc["ce_nmse"])),
                "goodput_bits_per_slot": float(acc["goodput_bits"]) / float(args.trials),
                "goodput_se_per_re": float(acc["goodput_bits"]) / float(max(int(args.trials) * grid.n_data_re, 1)),
                "coherence_bw_abs_0p5_sc": prow.get("coherence_bw_abs_0p5_sc", ""),
                "diversity_log10_product_norm": prow.get("diversity_log10_product_norm", ""),
                "condition_number": prow.get("condition_number", ""),
                "rpp_loaded_cond": float(estimators[pid].cond),
                "n_tx": int(args.n_tx),
                "n_rx": int(args.n_rx),
                "delay_spread_ns": float(args.delay_spread_ns),
                "dmrs_spacing_sc": int(args.dmrs_spacing_sc),
                "mcs_index": int(args.mcs_index),
                "qm": int(mcs.qm),
                "code_rate": float(mcs.code_rate),
                "tbs_bits": int(tb.tb_size),
                "coded_bits": int(tb.coded_bits),
            })
        write_csv(rows, out_dir / "pareto_link_curves_partial.csv")
        print(f"[done] SNR={snr_db:g} dB", flush=True)

    write_analysis_outputs(rows, out_dir, fig_dir, str(args.fig_prefix))

    with (out_dir / "run_config.json").open("w", encoding="utf-8") as f:
        json.dump(vars(args), f, indent=2, ensure_ascii=False)
    print(out_dir)
    return out_dir


def postprocess_from_partial(args: argparse.Namespace) -> Path:
    partial_arg = Path(args.postprocess_from)
    if partial_arg.is_dir():
        out_dir = partial_arg
        partial_path = out_dir / "pareto_link_curves_partial.csv"
    else:
        partial_path = partial_arg
        out_dir = partial_path.parent
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    with partial_path.open("r", encoding="utf-8", newline="") as f:
        rows: List[Dict[str, object]] = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"No rows found in {partial_path}")
    write_analysis_outputs(rows, out_dir, fig_dir, str(args.fig_prefix))

    with (out_dir / "postprocess_config.json").open("w", encoding="utf-8") as f:
        json.dump(vars(args), f, indent=2, ensure_ascii=False)
    def unique_sorted(field: str) -> List[object]:
        vals = sorted({str(r[field]) for r in rows}, key=lambda x: float(x) if x.replace(".", "", 1).isdigit() else x)
        return [float(v) if v.replace(".", "", 1).isdigit() else v for v in vals]

    completed = {
        "source_csv": str(partial_path),
        "n_rows": len(rows),
        "pareto_ids": sorted({str(r["pareto_id"]) for r in rows}, key=lambda x: (x[0], int(x[1:]))),
        "snr_db": unique_sorted("snr_db"),
        "trials_per_snr": unique_sorted("trials"),
        "n_tx": unique_sorted("n_tx"),
        "n_rx": unique_sorted("n_rx"),
        "delay_spread_ns": unique_sorted("delay_spread_ns"),
        "dmrs_spacing_sc": unique_sorted("dmrs_spacing_sc"),
        "mcs_index": unique_sorted("mcs_index"),
        "note": "Recovered from completed SNR rows after the long run was interrupted before finishing later SNR points.",
    }
    with (out_dir / "completed_run_summary.json").open("w", encoding="utf-8") as f:
        json.dump(completed, f, indent=2, ensure_ascii=False)
    print(out_dir)
    return out_dir


def merge_csvs(args: argparse.Namespace) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.merge_output_dir) if args.merge_output_dir else Path(args.out) / f"v_design_pareto_link_merged_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    rows_by_key: Dict[Tuple[str, float], Dict[str, object]] = {}
    sources = [Path(x.strip()) for x in str(args.merge_csvs).split(",") if x.strip()]
    for source in sources:
        with source.open("r", encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f):
                rows_by_key[(str(row["pareto_id"]), float(row["snr_db"]))] = row
    family_order = {"cdd_linear_phase": 0, "piecewise_linear_continuous": 1}
    rows = sorted(
        rows_by_key.values(),
        key=lambda r: (family_order.get(str(r["family"]), 99), int(str(r["pareto_id"])[1:]), float(r["snr_db"])),
    )
    write_analysis_outputs(rows, out_dir, fig_dir, str(args.fig_prefix))
    completed = {
        "source_csvs": [str(s) for s in sources],
        "n_rows": len(rows),
        "pareto_ids": sorted({str(r["pareto_id"]) for r in rows}, key=lambda x: (x[0], int(x[1:]))),
        "snr_db": sorted({float(r["snr_db"]) for r in rows}),
        "trials_per_snr": sorted({float(r["trials"]) for r in rows}),
        "n_tx": sorted({float(r["n_tx"]) for r in rows}),
        "n_rx": sorted({float(r["n_rx"]) for r in rows}),
        "delay_spread_ns": sorted({float(r["delay_spread_ns"]) for r in rows}),
        "dmrs_spacing_sc": sorted({float(r["dmrs_spacing_sc"]) for r in rows}),
        "mcs_index": sorted({float(r["mcs_index"]) for r in rows}),
    }
    with (out_dir / "completed_run_summary.json").open("w", encoding="utf-8") as f:
        json.dump(completed, f, indent=2, ensure_ascii=False)
    print(out_dir)
    return out_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pareto-csv", default="outputs/v_design_piecewise_tradeoff/v_design_piecewise_tradeoff_20260615_231257/pareto_front_labeled.csv")
    parser.add_argument("--out", default="outputs/v_design_pareto_link_curves")
    parser.add_argument("--fig-dir", default="docs/figures")
    parser.add_argument("--fig-prefix", default="experiment19")
    parser.add_argument("--postprocess-from", default="")
    parser.add_argument("--merge-csvs", default="")
    parser.add_argument("--merge-output-dir", default="")
    parser.add_argument("--n-tx", type=int, default=8)
    parser.add_argument("--n-rx", type=int, default=1)
    parser.add_argument("--n-segments", type=int, default=8)
    parser.add_argument("--n-prbs", type=int, default=48)
    parser.add_argument("--dmrs-spacing-sc", type=int, default=6)
    parser.add_argument("--delay-spread-ns", type=float, default=10.0)
    parser.add_argument("--max-delay-factor", type=float, default=8.0)
    parser.add_argument("--snrs", default="8,10,12,14,16,18,20,22")
    parser.add_argument("--trials", type=int, default=300)
    parser.add_argument("--mcs-index", type=int, default=8)
    parser.add_argument("--ldpc-iterations", type=int, default=8)
    parser.add_argument("--llr-clip", type=float, default=50.0)
    parser.add_argument("--loading", type=float, default=1e-8)
    parser.add_argument("--seed", type=int, default=20260616)
    parser.add_argument("--design-seed", type=int, default=20260615)
    parser.add_argument("--progress-every", type=int, default=50)
    args = parser.parse_args()
    if args.merge_csvs:
        merge_csvs(args)
    elif args.postprocess_from:
        postprocess_from_partial(args)
    else:
        run(args)


if __name__ == "__main__":
    main()
