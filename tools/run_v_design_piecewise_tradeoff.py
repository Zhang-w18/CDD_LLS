from __future__ import annotations

import argparse
import csv
import json
import math
import os
import html
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import sys
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ChannelConfig, ResourceConfig
from cdd_lls.core.mcs import build_tb_layout, get_mcs
from cdd_lls.phy.channel_tdl import generate_tdl_channel, make_exponential_pdp
from cdd_lls.phy.estimators import covariance_matrix, nmse
from cdd_lls.phy.precoding import equivalent_channel
from cdd_lls.phy.qam import qam_demapper_maxlog, qam_modulate
from cdd_lls.phy.resource_grid import ResourceGrid, build_resource_grid, local_indices_for_subcarriers
from tools.search_precoder_design_alg1 import make_alg1_full_cov_estimator


@dataclass(frozen=True)
class VDesign:
    label: str
    family: str
    V: np.ndarray
    phase: np.ndarray
    metadata: Dict[str, object]

    @property
    def n_tx(self) -> int:
        return int(self.V.shape[1])

    @property
    def C(self) -> np.ndarray:
        return self.V / math.sqrt(float(self.n_tx))


@dataclass(frozen=True)
class LinearEstimator:
    matrix: np.ndarray
    cond: float


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


def mean_or_nan(values: Iterable[float]) -> float:
    vals = [float(x) for x in values if np.isfinite(float(x))]
    return float(np.mean(vals)) if vals else float("nan")


def median_or_nan(values: Iterable[float]) -> float:
    vals = [float(x) for x in values if np.isfinite(float(x))]
    return float(np.median(vals)) if vals else float("nan")


def gain_db(reference: float, candidate: float) -> float:
    if not np.isfinite(reference) or not np.isfinite(candidate) or reference <= 0.0 or candidate <= 0.0:
        return float("nan")
    return float(10.0 * np.log10(reference / candidate))


def write_csv(rows: List[Dict[str, object]], path: Path) -> None:
    if not rows:
        return
    fields = sorted(set().union(*(row.keys() for row in rows)))
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def first_below(values: np.ndarray, threshold: float, start: int = 1) -> Optional[int]:
    arr = np.asarray(values, dtype=np.float64)
    idx = np.where(arr[int(start):] <= float(threshold))[0]
    if idx.size == 0:
        return None
    return int(idx[0] + int(start))


def finite_float(value: object, default: float = float("nan")) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float(default)
    return out if np.isfinite(out) else float(default)


def segment_edges(n_sc: int, n_segments: int) -> np.ndarray:
    edges = np.linspace(0, int(n_sc), int(n_segments) + 1)
    return np.rint(edges).astype(np.int64)


def phase_from_piecewise_delays(
    grid: ResourceGrid,
    delay_samples_by_segment: np.ndarray,
    *,
    reset_offsets: Optional[np.ndarray] = None,
    continuous: bool,
) -> np.ndarray:
    delays = np.asarray(delay_samples_by_segment, dtype=np.float64)
    n_segments, n_tx = delays.shape
    edges = segment_edges(grid.n_sc, n_segments)
    phase = np.empty((grid.n_sc, n_tx), dtype=np.float64)
    if continuous:
        start_phase = np.zeros(n_tx, dtype=np.float64)
        for s in range(n_segments):
            start = int(edges[s])
            stop = int(edges[s + 1])
            local = np.arange(stop - start, dtype=np.float64)
            phase[start:stop, :] = (
                start_phase[None, :]
                - 2.0 * np.pi * local[:, None] * delays[s, :][None, :] / float(grid.n_fft)
            )
            start_phase = (
                start_phase
                - 2.0 * np.pi * float(stop - start) * delays[s, :] / float(grid.n_fft)
            )
        return phase

    offsets = np.zeros_like(delays) if reset_offsets is None else np.asarray(reset_offsets, dtype=np.float64)
    if offsets.shape != delays.shape:
        raise ValueError("reset_offsets must have the same shape as delay_samples_by_segment.")
    for s in range(n_segments):
        start = int(edges[s])
        stop = int(edges[s + 1])
        local = np.arange(stop - start, dtype=np.float64)
        phase[start:stop, :] = (
            offsets[s, :][None, :]
            - 2.0 * np.pi * local[:, None] * delays[s, :][None, :] / float(grid.n_fft)
        )
    return phase


def design_from_phase(label: str, family: str, phase: np.ndarray, metadata: Dict[str, object]) -> VDesign:
    return VDesign(
        label=label,
        family=family,
        V=np.exp(1j * np.asarray(phase, dtype=np.float64)).astype(np.complex128),
        phase=np.asarray(phase, dtype=np.float64),
        metadata=dict(metadata),
    )


def cdd_design(grid: ResourceGrid, n_tx: int, step: int) -> VDesign:
    delays = int(step) * np.arange(int(n_tx), dtype=np.float64)
    phase = -2.0 * np.pi * grid.subcarrier_indices.astype(np.float64)[:, None] * delays[None, :] / float(grid.n_fft)
    return design_from_phase(
        label=f"CDD step={int(step)}",
        family="cdd_linear_phase",
        phase=phase,
        metadata={
            "cdd_step_samples": int(step),
            "delay_min_samples": float(np.min(delays)),
            "delay_max_samples": float(np.max(delays)),
        },
    )


def no_cdd_design(grid: ResourceGrid, n_tx: int) -> VDesign:
    phase = np.zeros((grid.n_sc, int(n_tx)), dtype=np.float64)
    return design_from_phase(
        label="No CDD",
        family="constant_rank1",
        phase=phase,
        metadata={"cdd_step_samples": 0},
    )


def dft_segment_offsets(n_segments: int, n_tx: int) -> np.ndarray:
    seg = np.arange(int(n_segments), dtype=np.float64)[:, None]
    tx = np.arange(int(n_tx), dtype=np.float64)[None, :]
    return 2.0 * np.pi * seg * tx / float(n_tx)


def balanced_random_permutation_indices(
    n_segments: int,
    n_tx: int,
    rng: np.random.Generator,
) -> np.ndarray:
    n_segments = int(n_segments)
    n_tx = int(n_tx)
    if n_segments > n_tx:
        raise ValueError(
            "Cannot assign a distinct local delay to every segment when n_segments > n_tx."
        )
    symbol_order = rng.permutation(n_tx)
    tx_offsets = rng.permutation(n_tx)
    segment_offsets = rng.permutation(n_tx)[:n_segments]
    latin_indices = (segment_offsets[:, None] + tx_offsets[None, :]) % n_tx
    return np.asarray(symbol_order[latin_indices], dtype=np.float64)


def centered_delay_alphabet(n_tx: int, slope_step: int) -> np.ndarray:
    n_tx = int(n_tx)
    return float(slope_step) * (
        np.arange(n_tx, dtype=np.float64) - 0.5 * float(n_tx - 1)
    )


def piecewise_reset_design(grid: ResourceGrid, n_tx: int, n_segments: int, slope_step: int) -> VDesign:
    idx = np.arange(int(n_tx), dtype=np.float64)
    delays = float(slope_step) * np.tile(idx[None, :], (int(n_segments), 1))
    offsets = dft_segment_offsets(n_segments, n_tx)
    phase = phase_from_piecewise_delays(
        grid,
        delays,
        reset_offsets=offsets,
        continuous=False,
    )
    return design_from_phase(
        label=f"PW reset DFTseg slope={int(slope_step)}",
        family="piecewise_linear_reset",
        phase=phase,
        metadata={
            "n_segments": int(n_segments),
            "local_slope_step_samples": int(slope_step),
            "phase_continuous": 0,
            "segment_phase_code": "dft",
        },
    )


def piecewise_continuous_perm_design(
    grid: ResourceGrid,
    n_tx: int,
    n_segments: int,
    slope_step: int,
    pattern: str,
    seed: int,
) -> VDesign:
    n_tx = int(n_tx)
    n_segments = int(n_segments)
    base = np.arange(n_tx, dtype=np.float64)
    if pattern == "cyclic":
        idx = np.vstack([np.roll(base, s) for s in range(n_segments)])
    elif pattern == "zigzag":
        idx = np.vstack([base if s % 2 == 0 else base[::-1] for s in range(n_segments)])
    elif pattern == "random":
        rng = np.random.default_rng(int(seed))
        idx = balanced_random_permutation_indices(n_segments, n_tx, rng)
    else:
        raise ValueError(f"Unknown piecewise pattern={pattern}.")
    delay_alphabet = centered_delay_alphabet(n_tx, slope_step)
    delays = delay_alphabet[np.asarray(idx, dtype=np.int64)]
    phase = phase_from_piecewise_delays(
        grid,
        delays,
        continuous=True,
    )
    suffix = f" #{int(seed) % 1000}" if pattern == "random" else ""
    return design_from_phase(
        label=f"PW cont {pattern} slope={int(slope_step)}{suffix}",
        family="piecewise_linear_continuous",
        phase=phase,
        metadata={
            "n_segments": n_segments,
            "local_slope_step_samples": int(slope_step),
            "local_delay_min_samples": float(np.min(delay_alphabet)),
            "local_delay_max_samples": float(np.max(delay_alphabet)),
            "local_delay_centered": 1,
            "phase_continuous": 1,
            "segment_delay_pattern": pattern,
            "per_tx_delays_unique_across_segments": int(
                pattern == "random" or (pattern == "cyclic" and n_segments <= n_tx)
            ),
            "design_seed": int(seed),
        },
    )


def build_designs(
    grid: ResourceGrid,
    n_tx: int,
    n_segments: int,
    cdd_steps: Sequence[int],
    piecewise_slopes: Sequence[int],
    random_per_slope: int,
    seed: int,
) -> List[VDesign]:
    designs: List[VDesign] = [no_cdd_design(grid, n_tx)]
    designs.extend(cdd_design(grid, n_tx, int(step)) for step in cdd_steps)
    for slope in piecewise_slopes:
        designs.append(piecewise_reset_design(grid, n_tx, n_segments, int(slope)))
        if int(slope) > 0:
            designs.append(piecewise_continuous_perm_design(
                grid,
                n_tx,
                n_segments,
                int(slope),
                pattern="cyclic",
                seed=stable_seed("cyclic", slope, base=seed),
            ))
            designs.append(piecewise_continuous_perm_design(
                grid,
                n_tx,
                n_segments,
                int(slope),
                pattern="zigzag",
                seed=stable_seed("zigzag", slope, base=seed),
            ))
            for i in range(int(random_per_slope)):
                designs.append(piecewise_continuous_perm_design(
                    grid,
                    n_tx,
                    n_segments,
                    int(slope),
                    pattern="random",
                    seed=stable_seed("random", slope, i, base=seed),
                ))
    unique: Dict[str, VDesign] = {}
    for design in designs:
        unique.setdefault(design.label, design)
    return list(unique.values())


def matrix_metrics(design: VDesign, grid: ResourceGrid, n_segments: int) -> Dict[str, object]:
    V = np.asarray(design.V, dtype=np.complex128)
    K = int(V.shape[0])
    N = int(V.shape[1])
    gram = V.conj().T @ V
    eig = np.linalg.eigvalsh(gram)
    eig = np.maximum(np.real(eig), 0.0)
    eig_norm = eig / float(K)
    eps = 1e-15
    logdet_norm = float(np.sum(np.log(np.maximum(eig_norm, eps))))
    log10_prod_norm = float(logdet_norm / np.log(10.0))
    prod_norm = float(np.exp(logdet_norm)) if logdet_norm > -700 else 0.0
    min_eig_norm = float(np.min(eig_norm)) if eig_norm.size else float("nan")
    max_eig_norm = float(np.max(eig_norm)) if eig_norm.size else float("nan")
    cond = float(math.sqrt(max_eig_norm / max(min_eig_norm, eps))) if max_eig_norm > 0.0 else float("inf")
    col_norm = np.sqrt(np.maximum(np.real(np.diag(gram)), eps))
    coh = np.abs(gram / (col_norm[:, None] * col_norm[None, :]))
    if N > 1:
        mutual = float(np.max(coh[~np.eye(N, dtype=bool)]))
    else:
        mutual = 0.0

    dphi = np.angle(V[1:, :] * V[:-1, :].conj())
    group_delay_samples = -float(grid.n_fft) * dphi / (2.0 * np.pi)
    gd_spread = float(np.max(group_delay_samples) - np.min(group_delay_samples))
    gd_rms = float(np.std(group_delay_samples))
    rough1 = float(np.mean(np.abs(dphi) ** 2))
    d2 = np.angle(V[2:, :] * V[1:-1, :].conj()) - np.angle(V[1:-1, :] * V[:-2, :].conj())
    rough2 = float(np.mean(d2 ** 2)) if d2.size else 0.0
    edges = segment_edges(K, n_segments)
    boundary_steps = []
    for b in edges[1:-1]:
        if 0 < int(b) < K:
            boundary_steps.extend(np.abs(np.angle(V[int(b), :] * V[int(b) - 1, :].conj())).tolist())
    boundary_jump_max = float(np.max(boundary_steps)) if boundary_steps else 0.0
    boundary_jump_rms = float(np.sqrt(np.mean(np.asarray(boundary_steps) ** 2))) if boundary_steps else 0.0

    return {
        "diversity_logdet_norm": logdet_norm,
        "diversity_log10_product_norm": log10_prod_norm,
        "diversity_product_norm": prod_norm,
        "min_eig_norm": min_eig_norm,
        "max_eig_norm": max_eig_norm,
        "condition_number": cond,
        "mutual_coherence": mutual,
        "group_delay_spread_samples": gd_spread,
        "group_delay_rms_samples": gd_rms,
        "roughness_first": rough1,
        "roughness_second": rough2,
        "boundary_phase_step_max_rad": boundary_jump_max,
        "boundary_phase_step_rms_rad": boundary_jump_rms,
    }


def covariance_metrics(
    design: VDesign,
    base_cov: np.ndarray,
    max_delta: int,
) -> Tuple[Dict[str, object], np.ndarray, np.ndarray]:
    C = design.C
    inner = C @ C.conj().T
    R = np.asarray(base_cov, dtype=np.complex128) * inner
    diag = np.sqrt(np.maximum(np.real(np.diag(R)), 1e-15))
    Rn = R / (diag[:, None] * diag[None, :])
    max_delta = min(int(max_delta), Rn.shape[0] - 1)
    rho_abs = np.empty(max_delta + 1, dtype=np.float64)
    rho_complex_abs = np.empty(max_delta + 1, dtype=np.float64)
    rho_power = np.empty(max_delta + 1, dtype=np.float64)
    for delta in range(max_delta + 1):
        vals = np.diag(Rn, k=delta)
        rho_abs[delta] = float(np.mean(np.abs(vals)))
        rho_complex_abs[delta] = float(abs(np.mean(vals)))
        rho_power[delta] = float(np.mean(np.abs(vals) ** 2))
    bw_abs_09 = first_below(rho_abs, 0.9)
    bw_abs_05 = first_below(rho_abs, 0.5)
    bw_complex_05 = first_below(rho_complex_abs, 0.5)
    area_abs = float(np.sum(rho_abs))
    area_power = float(np.sum(rho_power))
    eff_rank = float((np.trace(R).real ** 2) / max(np.sum(np.abs(R) ** 2).real, 1e-15))
    return (
        {
            "coherence_bw_abs_0p9_sc": "" if bw_abs_09 is None else int(bw_abs_09),
            "coherence_bw_abs_0p5_sc": "" if bw_abs_05 is None else int(bw_abs_05),
            "coherence_bw_complex_0p5_sc": "" if bw_complex_05 is None else int(bw_complex_05),
            "coherence_area_abs_sc": area_abs,
            "coherence_area_power_sc": area_power,
            "covariance_effective_rank": eff_rank,
        },
        rho_abs,
        rho_complex_abs,
    )


def effective_covariance_from_design(base_cov: np.ndarray, design: VDesign) -> np.ndarray:
    return np.asarray(base_cov, dtype=np.complex128) * (design.C @ design.C.conj().T)


def make_alg1_estimator(grid: ResourceGrid, pdp: np.ndarray, design: VDesign, noise_var_ls: float, loading: float) -> LinearEstimator:
    est = make_alg1_full_cov_estimator(
        grid=grid,
        pdp=pdp,
        C=design.C,
        noise_var_ls=float(noise_var_ls),
        loading=float(loading),
    )
    return LinearEstimator(
        matrix=est.matrix,
        cond=float(est.metadata["rpp_loaded_cond"]),
    )


def pareto_front(rows: List[Dict[str, object]], x_key: str, y_key: str) -> List[Dict[str, object]]:
    pts = [
        row for row in rows
        if np.isfinite(finite_float(row.get(x_key))) and np.isfinite(finite_float(row.get(y_key)))
    ]
    front: List[Dict[str, object]] = []
    for row in pts:
        x = finite_float(row[x_key])
        y = finite_float(row[y_key])
        dominated = False
        for other in pts:
            if other is row:
                continue
            ox = finite_float(other[x_key])
            oy = finite_float(other[y_key])
            if ox >= x and oy >= y and (ox > x or oy > y):
                dominated = True
                break
        if not dominated:
            front.append(row)
    return sorted(front, key=lambda r: (finite_float(r[x_key]), finite_float(r[y_key])))


def pareto_id_prefix(family: str) -> str:
    if str(family) == "cdd_linear_phase":
        return "C"
    if str(family) == "piecewise_linear_continuous":
        return "P"
    if str(family) == "piecewise_linear_reset":
        return "R"
    if str(family) == "piecewise_balanced_signed":
        return "N"
    return "F"


def label_pareto_rows(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    labeled: List[Dict[str, object]] = []
    family_order = ["cdd_linear_phase", "piecewise_linear_continuous", "piecewise_linear_reset"]
    for family in family_order:
        pts = [r for r in rows if str(r.get("family")) == family]
        pts = sorted(
            pts,
            key=lambda r: (
                -finite_float(r.get("coherence_bw_abs_0p5_sc"), default=-1.0),
                -finite_float(r.get("diversity_log10_product_norm"), default=-1e9),
            ),
        )
        prefix = pareto_id_prefix(family)
        for idx, row in enumerate(pts, start=1):
            out = dict(row)
            out["pareto_id"] = f"{prefix}{idx}"
            labeled.append(out)
    for row in rows:
        if str(row.get("family")) not in family_order:
            labeled.append(dict(row))
    return labeled


def best_cdd_y_at_or_above_x(cdd_rows: List[Dict[str, object]], x: float, y_key: str) -> float:
    vals = [
        float(row[y_key]) for row in cdd_rows
        if np.isfinite(finite_float(row.get("coherence_bw_abs_0p5_sc")))
        and finite_float(row["coherence_bw_abs_0p5_sc"]) >= float(x)
    ]
    return max(vals) if vals else float("-inf")


def choose_piecewise_winner(rows: List[Dict[str, object]]) -> Tuple[Optional[str], Optional[str], float]:
    cdd_rows = [r for r in rows if str(r["family"]) == "cdd_linear_phase"]
    candidate_rows = [
        r for r in rows
        if str(r["family"]).startswith("piecewise")
        and np.isfinite(finite_float(r.get("coherence_bw_abs_0p5_sc")))
    ]
    best_label = None
    best_ref = None
    best_margin = float("-inf")
    for row in candidate_rows:
        x = finite_float(row["coherence_bw_abs_0p5_sc"])
        y = finite_float(row["diversity_log10_product_norm"])
        cdd_y = best_cdd_y_at_or_above_x(cdd_rows, x, "diversity_log10_product_norm")
        margin = y - cdd_y
        if margin > best_margin:
            best_margin = margin
            best_label = str(row["precoder"])
            eligible = [
                r for r in cdd_rows
                if np.isfinite(finite_float(r.get("coherence_bw_abs_0p5_sc")))
                and finite_float(r["coherence_bw_abs_0p5_sc"]) >= x
            ]
            if eligible:
                ref = max(eligible, key=lambda r: finite_float(r["diversity_log10_product_norm"]))
            else:
                ref = min(cdd_rows, key=lambda r: abs(finite_float(r["coherence_bw_abs_0p5_sc"], default=1e9) - x))
            best_ref = str(ref["precoder"])
    return best_label, best_ref, float(best_margin)


def _nice_bounds(values: Sequence[float], pad: float = 0.06, floor_zero: bool = False) -> Tuple[float, float]:
    vals = np.asarray([v for v in values if np.isfinite(float(v))], dtype=np.float64)
    if vals.size == 0:
        return 0.0, 1.0
    lo = float(np.min(vals))
    hi = float(np.max(vals))
    if abs(hi - lo) < 1e-12:
        delta = max(abs(hi), 1.0) * 0.1
        lo -= delta
        hi += delta
    span = hi - lo
    lo -= pad * span
    hi += pad * span
    if floor_zero:
        lo = min(0.0, lo)
    return lo, hi


def _svg_polyline(points: List[Tuple[float, float]], color: str, width: float = 2.0, dashed: bool = False) -> str:
    if not points:
        return ""
    data = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
    dash = ' stroke-dasharray="5 4"' if dashed else ""
    return f'<polyline points="{data}" fill="none" stroke="{color}" stroke-width="{width}"{dash}/>'


def _svg_marker(x: float, y: float, color: str, shape: str, label: str, size: float = 5.2) -> str:
    title = f"<title>{html.escape(label)}</title>"
    if shape == "square":
        return f'<rect x="{x-size:.2f}" y="{y-size:.2f}" width="{2*size:.2f}" height="{2*size:.2f}" fill="{color}" opacity="0.72">{title}</rect>'
    if shape == "triangle":
        pts = f"{x:.2f},{y-size:.2f} {x-size:.2f},{y+size:.2f} {x+size:.2f},{y+size:.2f}"
        return f'<polygon points="{pts}" fill="{color}" opacity="0.72">{title}</polygon>'
    if shape == "diamond":
        pts = f"{x:.2f},{y-size:.2f} {x-size:.2f},{y:.2f} {x:.2f},{y+size:.2f} {x+size:.2f},{y:.2f}"
        return f'<polygon points="{pts}" fill="{color}" opacity="0.72">{title}</polygon>'
    if shape == "cross":
        return (
            f'<g opacity="0.8">{title}'
            f'<line x1="{x-size:.2f}" y1="{y-size:.2f}" x2="{x+size:.2f}" y2="{y+size:.2f}" stroke="{color}" stroke-width="1.7"/>'
            f'<line x1="{x-size:.2f}" y1="{y+size:.2f}" x2="{x+size:.2f}" y2="{y-size:.2f}" stroke="{color}" stroke-width="1.7"/>'
            "</g>"
        )
    return f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{size:.2f}" fill="{color}" opacity="0.72">{title}</circle>'


def _write_svg_plot(
    path: Path,
    *,
    title: str,
    xlabel: str,
    ylabel: str,
    x_bounds: Tuple[float, float],
    y_bounds: Tuple[float, float],
    elements: List[str],
    legend: List[Tuple[str, str, str]],
    width: int = 1040,
    height: int = 640,
    top_margin: int = 58,
) -> None:
    left, right, top, bottom = 86, 28, int(top_margin), 82
    plot_w = width - left - right
    plot_h = height - top - bottom
    x0, x1 = x_bounds
    y0, y1 = y_bounds
    grid = []
    for i in range(6):
        tx = left + plot_w * i / 5.0
        val = x0 + (x1 - x0) * i / 5.0
        grid.append(f'<line x1="{tx:.2f}" y1="{top}" x2="{tx:.2f}" y2="{top+plot_h}" stroke="#e6e8eb" stroke-width="1"/>')
        grid.append(f'<text x="{tx:.2f}" y="{top+plot_h+24}" font-size="12" text-anchor="middle" fill="#333">{val:.3g}</text>')
    for i in range(6):
        ty = top + plot_h - plot_h * i / 5.0
        val = y0 + (y1 - y0) * i / 5.0
        grid.append(f'<line x1="{left}" y1="{ty:.2f}" x2="{left+plot_w}" y2="{ty:.2f}" stroke="#e6e8eb" stroke-width="1"/>')
        grid.append(f'<text x="{left-10}" y="{ty+4:.2f}" font-size="12" text-anchor="end" fill="#333">{val:.3g}</text>')
    legend_items = []
    lx = left + 10
    ly = 52 if top >= 100 else top + 12
    for idx, (name, color, shape) in enumerate(legend[:10]):
        x = lx + (idx % 2) * 310
        y = ly + (idx // 2) * 20
        legend_items.append(_svg_marker(x, y - 4, color, shape, name, size=4.2))
        legend_items.append(f'<text x="{x+12}" y="{y:.2f}" font-size="12" fill="#333">{html.escape(name)}</text>')
    body = "\n".join(grid + elements + legend_items)
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="white"/>
<text x="{width/2:.2f}" y="28" font-size="18" font-family="Arial, sans-serif" text-anchor="middle" fill="#111">{html.escape(title)}</text>
<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#333" stroke-width="1.2"/>
<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#333" stroke-width="1.2"/>
{body}
<text x="{left+plot_w/2:.2f}" y="{height-24}" font-size="14" font-family="Arial, sans-serif" text-anchor="middle" fill="#111">{html.escape(xlabel)}</text>
<text x="20" y="{top+plot_h/2:.2f}" font-size="14" font-family="Arial, sans-serif" text-anchor="middle" fill="#111" transform="rotate(-90 20 {top+plot_h/2:.2f})">{html.escape(ylabel)}</text>
</svg>
'''
    path.write_text(svg, encoding="utf-8")


def plot_tradeoff(
    rows: List[Dict[str, object]],
    pareto_rows: List[Dict[str, object]],
    path: Path,
    title: str = "8TX/1RX V-matrix diversity vs CE-coherence tradeoff",
) -> None:
    families = sorted(set(str(r["family"]) for r in rows))
    colors = {
        "constant_rank1": "#8d99ae",
        "cdd_linear_phase": "#26547c",
        "piecewise_linear_reset": "#ef476f",
        "piecewise_linear_continuous": "#06a77d",
        "piecewise_linear_continuous_legacy": "#8a9a9a",
        "piecewise_balanced_signed": "#d97706",
    }
    shapes = {
        "constant_rank1": "cross",
        "cdd_linear_phase": "circle",
        "piecewise_linear_reset": "square",
        "piecewise_linear_continuous": "triangle",
        "piecewise_linear_continuous_legacy": "cross",
        "piecewise_balanced_signed": "diamond",
    }
    pts_all = [
        (finite_float(r.get("coherence_bw_abs_0p5_sc")), finite_float(r.get("diversity_log10_product_norm")))
        for r in rows
    ]
    xs = [x for x, y in pts_all if np.isfinite(x) and np.isfinite(y)]
    ys = [y for x, y in pts_all if np.isfinite(x) and np.isfinite(y)]
    x_bounds = _nice_bounds(xs, floor_zero=True)
    y_bounds = _nice_bounds(ys)
    left, right, top, bottom = 86, 28, 120, 82
    width, height = 1040, 700
    plot_w = width - left - right
    plot_h = height - top - bottom

    def sx(x: float) -> float:
        return left + (x - x_bounds[0]) / max(x_bounds[1] - x_bounds[0], 1e-12) * plot_w

    def sy(y: float) -> float:
        return top + plot_h - (y - y_bounds[0]) / max(y_bounds[1] - y_bounds[0], 1e-12) * plot_h

    elements: List[str] = []
    for family in families:
        pts = [
            r for r in rows
            if str(r["family"]) == family
            and np.isfinite(finite_float(r.get("coherence_bw_abs_0p5_sc")))
        ]
        for row in pts:
            elements.append(_svg_marker(
                sx(finite_float(row["coherence_bw_abs_0p5_sc"])),
                sy(finite_float(row["diversity_log10_product_norm"])),
                colors.get(family, "#444"),
                shapes.get(family, "circle"),
                str(row["precoder"]),
            ))
    for family, color, dashed in (
        ("cdd_linear_phase", "#102a43", False),
        ("piecewise_linear_reset", "#c9184a", False),
        ("piecewise_linear_continuous", "#007f5f", False),
        ("piecewise_linear_continuous_legacy", "#697979", True),
        ("piecewise_balanced_signed", "#b45309", False),
    ):
        front = [
            r for r in pareto_rows
            if str(r["family"]) == family
            and np.isfinite(finite_float(r.get("coherence_bw_abs_0p5_sc")))
        ]
        front = sorted(front, key=lambda r: finite_float(r["coherence_bw_abs_0p5_sc"]))
        elements.append(_svg_polyline(
            [(sx(finite_float(r["coherence_bw_abs_0p5_sc"])), sy(finite_float(r["diversity_log10_product_norm"]))) for r in front],
            color,
            width=2.2,
            dashed=dashed,
        ))
        for row in front:
            label = str(row.get("pareto_id", ""))
            if not label:
                continue
            x = sx(finite_float(row["coherence_bw_abs_0p5_sc"]))
            y = sy(finite_float(row["diversity_log10_product_norm"]))
            elements.append(
                f'<text x="{x+6:.2f}" y="{y-6:.2f}" font-size="11" '
                f'font-family="Arial, sans-serif" fill="{color}" font-weight="700">'
                f'{html.escape(label)}</text>'
            )
    family_labels = {
        "constant_rank1": "No CDD",
        "cdd_linear_phase": "CDD",
        "piecewise_linear_reset": "PW reset diagnostic",
        "piecewise_linear_continuous": "PW continuous",
        "piecewise_linear_continuous_legacy": "Legacy PW random",
        "piecewise_balanced_signed": "New signed balanced PW",
    }
    legend = [
        (family_labels.get(fam, fam), colors.get(fam, "#444"), shapes.get(fam, "circle"))
        for fam in families
    ]
    _write_svg_plot(
        path,
        title=title,
        xlabel="CE gain proxy: coherence BW at mean |rho| <= 0.5 (subcarriers)",
        ylabel="Diversity gain: log10 det((V^H V)/K)",
        x_bounds=x_bounds,
        y_bounds=y_bounds,
        elements=elements,
        legend=legend,
        width=width,
        height=height,
        top_margin=top,
    )


def plot_covariance_curve(curve_rows: List[Dict[str, object]], path: Path) -> None:
    groups: Dict[str, List[Dict[str, object]]] = {}
    for row in curve_rows:
        groups.setdefault(str(row["precoder"]), []).append(row)
    colors = ["#ef476f", "#26547c", "#06a77d", "#7b2cbf"]
    xs = [int(r["delta_sc"]) for r in curve_rows]
    x_bounds = _nice_bounds(xs, floor_zero=True)
    y_bounds = (0.0, 1.05)
    left, right, top, bottom = 86, 28, 58, 82
    width, height = 1000, 600
    plot_w = width - left - right
    plot_h = height - top - bottom

    def sx(x: float) -> float:
        return left + (x - x_bounds[0]) / max(x_bounds[1] - x_bounds[0], 1e-12) * plot_w

    def sy(y: float) -> float:
        return top + plot_h - (y - y_bounds[0]) / max(y_bounds[1] - y_bounds[0], 1e-12) * plot_h

    elements = []
    elements.append(_svg_polyline([(sx(x_bounds[0]), sy(0.5)), (sx(x_bounds[1]), sy(0.5))], "#111", width=1.0, dashed=True))
    legend: List[Tuple[str, str, str]] = []
    for idx, (label, pts) in enumerate(sorted(groups.items())):
        pts = sorted(pts, key=lambda r: int(r["delta_sc"]))
        color = colors[idx % len(colors)]
        elements.append(_svg_polyline([(sx(int(r["delta_sc"])), sy(float(r["rho_abs_mean"]))) for r in pts], color, width=2.0))
        elements.append(_svg_polyline([(sx(int(r["delta_sc"])), sy(float(r["rho_complex_abs"]))) for r in pts], color, width=1.2, dashed=True))
        legend.append((label, color, "circle"))
    _write_svg_plot(
        path,
        title="Selected V-matrix covariance function",
        xlabel="Frequency separation Delta k (subcarriers)",
        ylabel="Normalized covariance",
        x_bounds=x_bounds,
        y_bounds=y_bounds,
        elements=elements,
        legend=legend,
        width=width,
        height=height,
    )


def plot_metric_curve(rows: List[Dict[str, object]], metric: str, path: Path, ylabel: str) -> None:
    if not rows:
        return
    groups: Dict[str, List[Dict[str, object]]] = {}
    for row in rows:
        groups.setdefault(str(row["precoder"]), []).append(row)
    log_y = metric in ("ce_nmse_mean", "bler")
    xs = [float(r["snr_db"]) for r in rows]
    raw_y = [max(float(r[metric]), 1e-5) for r in rows]
    ys = [math.log10(y) for y in raw_y] if log_y else raw_y
    x_bounds = _nice_bounds(xs)
    y_bounds = _nice_bounds(ys)
    left, right, top, bottom = 86, 28, 58, 82
    width, height = 900, 560
    plot_w = width - left - right
    plot_h = height - top - bottom

    def sx(x: float) -> float:
        return left + (x - x_bounds[0]) / max(x_bounds[1] - x_bounds[0], 1e-12) * plot_w

    def sy(y: float) -> float:
        return top + plot_h - (y - y_bounds[0]) / max(y_bounds[1] - y_bounds[0], 1e-12) * plot_h

    colors = ["#ef476f", "#26547c", "#06a77d", "#7b2cbf"]
    elements: List[str] = []
    legend: List[Tuple[str, str, str]] = []
    for idx, (label, pts) in enumerate(sorted(groups.items())):
        pts = sorted(pts, key=lambda r: float(r["snr_db"]))
        color = colors[idx % len(colors)]
        coords = []
        for row in pts:
            y = max(float(row[metric]), 1e-5)
            yy = math.log10(y) if log_y else y
            coords.append((sx(float(row["snr_db"])), sy(yy)))
            elements.append(_svg_marker(coords[-1][0], coords[-1][1], color, "circle", label, size=4.4))
        elements.append(_svg_polyline(coords, color, width=1.8))
        legend.append((label, color, "circle"))
    ylab = f"log10 {ylabel}" if log_y else ylabel
    _write_svg_plot(
        path,
        title=ylabel,
        xlabel="SNR (dB)",
        ylabel=ylab,
        x_bounds=x_bounds,
        y_bounds=y_bounds,
        elements=elements,
        legend=legend,
        width=width,
        height=height,
    )


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


def decode_same_tb_batch(adapter, llrs_by_candidate: List[np.ndarray], reference_payload_bits_by_cb: List[np.ndarray]):
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
    from cdd_lls.phy.ldpc import summarize_decode

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


def run_nmse_compare(
    args: argparse.Namespace,
    grid: ResourceGrid,
    pdp: np.ndarray,
    base_cov: np.ndarray,
    designs: List[VDesign],
    labels: Sequence[str],
    out_dir: Path,
    fig_dir: Path,
) -> List[Dict[str, object]]:
    selected = [d for d in designs if d.label in set(labels)]
    if len(selected) < 2:
        return []
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    channel = ChannelConfig(delay_spread_ns=float(args.delay_spread_ns), max_delay_factor=float(args.max_delay_factor))
    rows: List[Dict[str, object]] = []
    for snr_db in parse_float_list(args.nmse_snrs):
        noise_var = 10.0 ** (-float(snr_db) / 10.0)
        noise_var_ls = noise_var / float(len(args.dmrs_symbol_indices))
        est_by_label = {
            d.label: make_alg1_estimator(
                grid,
                pdp,
                d,
                noise_var_ls,
                float(args.loading),
            )
            for d in selected
        }
        accum = {d.label: [] for d in selected}
        rng = np.random.default_rng(stable_seed("nmse", snr_db, base=int(args.seed)))
        for trial in range(int(args.nmse_trials)):
            realization = generate_tdl_channel(rng, grid, channel, n_tx=int(args.n_tx), n_rx=int(args.n_rx))
            H = realization.H
            shared_noise = (
                rng.normal(size=(int(args.n_rx), len(pilot_local)))
                + 1j * rng.normal(size=(int(args.n_rx), len(pilot_local)))
            ) * math.sqrt(noise_var_ls / 2.0)
            for design in selected:
                g = equivalent_channel(H, design.C)
                ls_obs = g[:, pilot_local] + shared_noise
                g_hat = ls_obs @ est_by_label[design.label].matrix.T
                accum[design.label].append(nmse(g, g_hat))
        cdd_vals = {
            d.label: mean_or_nan(accum[d.label])
            for d in selected
            if d.family == "cdd_linear_phase"
        }
        best_cdd_nmse = min(cdd_vals.values()) if cdd_vals else float("nan")
        for design in selected:
            mean_nmse = mean_or_nan(accum[design.label])
            rows.append({
                "precoder": design.label,
                "family": design.family,
                "snr_db": float(snr_db),
                "delay_spread_ns": float(args.delay_spread_ns),
                "trials": int(args.nmse_trials),
                "n_tx": int(args.n_tx),
                "n_rx": int(args.n_rx),
                "ce_nmse_mean": mean_nmse,
                "ce_nmse_median": median_or_nan(accum[design.label]),
                "gain_vs_cdd_db": gain_db(best_cdd_nmse, mean_nmse),
                "rpp_loaded_cond": float(est_by_label[design.label].cond),
            })
    write_csv(rows, out_dir / "nmse_comparison.csv")
    plot_metric_curve(rows, "ce_nmse_mean", fig_dir / "experiment18_v_design_nmse.svg", "Alg1 RMMSE NMSE")
    return rows


def run_bler_compare(
    args: argparse.Namespace,
    grid: ResourceGrid,
    pdp: np.ndarray,
    base_cov: np.ndarray,
    designs: List[VDesign],
    labels: Sequence[str],
    out_dir: Path,
    fig_dir: Path,
) -> List[Dict[str, object]]:
    selected = [d for d in designs if d.label in set(labels)]
    if len(selected) < 2 or int(args.bler_trials) <= 0:
        return []
    from cdd_lls.phy.ldpc import SionnaLDPCAdapter

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
    for snr_db in parse_float_list(args.bler_snrs):
        noise_var = 10.0 ** (-float(snr_db) / 10.0)
        noise_var_ls = noise_var / float(len(args.dmrs_symbol_indices))
        est_by_label = {
            d.label: make_alg1_estimator(
                grid,
                pdp,
                d,
                noise_var_ls,
                float(args.loading),
            )
            for d in selected
        }
        accum = {
            d.label: {
                "tb_errors": 0,
                "cb_errors": 0,
                "goodput_bits": 0,
                "ce_nmse": [],
            }
            for d in selected
        }
        base_seed = stable_seed("bler", snr_db, base=int(args.seed))
        for trial in range(1, int(args.bler_trials) + 1):
            rng = np.random.default_rng(stable_seed(base_seed, trial, base=int(args.seed)))
            realization = generate_tdl_channel(rng, grid, channel, n_tx=int(args.n_tx), n_rx=int(args.n_rx))
            H = realization.H
            payload = [
                rng.integers(0, 2, size=int(k), dtype=np.int8)
                for k in tb.cb_k_values
            ]
            coded = adapter.encode(payload)
            coded_cw = np.concatenate(coded)
            symbols = qam_modulate(coded_cw, int(mcs.qm))
            pilot_noise = (
                rng.normal(size=(int(args.n_rx), len(pilot_local)))
                + 1j * rng.normal(size=(int(args.n_rx), len(pilot_local)))
            ) * math.sqrt(noise_var_ls / 2.0)
            data_noise = (
                rng.normal(size=(int(args.n_rx), len(data_local)))
                + 1j * rng.normal(size=(int(args.n_rx), len(data_local)))
            ) * math.sqrt(noise_var / 2.0)
            trial_llrs = []
            trial_nmse = []
            trial_designs = []
            for design in selected:
                g = equivalent_channel(H, design.C)
                ls_obs = g[:, pilot_local] + pilot_noise
                g_hat = ls_obs @ est_by_label[design.label].matrix.T
                y = g[:, data_local] * symbols[None, :] + data_noise
                z, no_eff = equalize_mrc(y, g_hat[:, data_local], noise_var)
                trial_llrs.append(qam_demapper_maxlog(z, no_eff, int(mcs.qm)))
                trial_nmse.append(nmse(g, g_hat))
                trial_designs.append(design)
            dec_results = decode_same_tb_batch(adapter, trial_llrs, payload)
            for design, dec, ce in zip(trial_designs, dec_results, trial_nmse):
                acc = accum[design.label]
                acc["tb_errors"] += int(not dec.tb_success)
                acc["cb_errors"] += int(sum(1 for ok in dec.cb_success if not ok))
                acc["goodput_bits"] += int(dec.goodput_bits)
                acc["ce_nmse"].append(float(ce))
        cdd_blers = {
            d.label: float(accum[d.label]["tb_errors"]) / float(args.bler_trials)
            for d in selected
            if d.family == "cdd_linear_phase"
        }
        best_cdd_bler = min(cdd_blers.values()) if cdd_blers else float("nan")
        for design in selected:
            acc = accum[design.label]
            bler = float(acc["tb_errors"]) / float(args.bler_trials)
            rows.append({
                "precoder": design.label,
                "family": design.family,
                "snr_db": float(snr_db),
                "delay_spread_ns": float(args.delay_spread_ns),
                "trials": int(args.bler_trials),
                "tb_errors": int(acc["tb_errors"]),
                "cb_errors": int(acc["cb_errors"]),
                "bler": bler,
                "bler_delta_vs_cdd": float(bler - best_cdd_bler),
                "ce_nmse_mean": mean_or_nan(acc["ce_nmse"]),
                "goodput_bits_per_slot": float(acc["goodput_bits"]) / float(args.bler_trials),
                "goodput_se_per_re": float(acc["goodput_bits"]) / float(max(int(args.bler_trials) * grid.n_data_re, 1)),
                "n_tx": int(args.n_tx),
                "n_rx": int(args.n_rx),
                "mcs_index": int(args.mcs_index),
                "qm": int(mcs.qm),
                "code_rate": float(mcs.code_rate),
                "tbs_bits": int(tb.tb_size),
                "coded_bits": int(tb.coded_bits),
                "n_cbs": int(len(tb.cb_k_values)),
            })
    write_csv(rows, out_dir / "bler_comparison.csv")
    plot_metric_curve(rows, "bler", fig_dir / "experiment18_v_design_bler.svg", "TB BLER")
    return rows


def run(args: argparse.Namespace) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out) / f"v_design_piecewise_tradeoff_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    dmrs_symbols = parse_int_list(args.dmrs_symbols)
    args.dmrs_symbol_indices = dmrs_symbols
    resource = ResourceConfig(
        carrier_bandwidth_mhz=float(args.carrier_bandwidth_mhz),
        scs_khz=int(args.scs_khz),
        n_fft=int(args.n_fft),
        n_prbs=int(args.n_prbs),
        pdsch_n_symbols=int(args.pdsch_n_symbols),
        dmrs_symbol_indices=dmrs_symbols,
        dmrs_spacing_sc=int(args.dmrs_spacing_sc),
        prg_size_rb=4,
    )
    grid = build_resource_grid(resource)
    if grid.n_sc % int(args.n_segments) != 0:
        raise ValueError("n_sc must be divisible by n_segments for the DFT segment code experiment.")

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

    metric_rows: List[Dict[str, object]] = []
    cov_curves: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for design in designs:
        mm = matrix_metrics(design, grid, int(args.n_segments))
        cm, rho_abs, rho_complex_abs = covariance_metrics(design, base_cov, int(args.max_corr_delta))
        cov_curves[design.label] = (rho_abs, rho_complex_abs)
        row = {
            "precoder": design.label,
            "family": design.family,
            "n_tx": int(args.n_tx),
            "n_rx": int(args.n_rx),
            "n_sc": int(grid.n_sc),
            "n_prbs": int(args.n_prbs),
            "n_segments": int(args.n_segments),
            "segment_len_sc": int(grid.n_sc // int(args.n_segments)),
            "delay_spread_ns": float(args.delay_spread_ns),
            "sample_period_ns": float(sample_period_ns),
            "pdp_taps": int(len(pdp)),
            "dmrs_spacing_sc": int(args.dmrs_spacing_sc),
            "dmrs_symbols": int(len(dmrs_symbols)),
            "pilot_count_combined": int(grid.pilot_count),
        }
        row.update(mm)
        row.update(cm)
        row.update({f"design_{k}": v for k, v in design.metadata.items()})
        metric_rows.append(row)

    cdd_front = pareto_front([r for r in metric_rows if r["family"] == "cdd_linear_phase"], "coherence_bw_abs_0p5_sc", "diversity_log10_product_norm")
    pw_reset_front = pareto_front([r for r in metric_rows if r["family"] == "piecewise_linear_reset"], "coherence_bw_abs_0p5_sc", "diversity_log10_product_norm")
    pw_cont_front = pareto_front([r for r in metric_rows if r["family"] == "piecewise_linear_continuous"], "coherence_bw_abs_0p5_sc", "diversity_log10_product_norm")
    front_rows = label_pareto_rows(cdd_front + pw_cont_front + pw_reset_front)

    winner_label, ref_cdd_label, margin = choose_piecewise_winner(metric_rows)
    if winner_label is None:
        piecewise_rows = [r for r in metric_rows if str(r["family"]).startswith("piecewise")]
        winner_label = str(piecewise_rows[0]["precoder"] if piecewise_rows else metric_rows[0]["precoder"])
    if ref_cdd_label is None:
        cdd_rows = [r for r in metric_rows if r["family"] == "cdd_linear_phase"]
        ref_cdd_label = str(max(cdd_rows, key=lambda r: finite_float(r["diversity_log10_product_norm"]))["precoder"]) if cdd_rows else winner_label

    curve_rows: List[Dict[str, object]] = []
    for label in [winner_label, ref_cdd_label]:
        rho_abs, rho_complex_abs = cov_curves[label]
        for delta in range(len(rho_abs)):
            curve_rows.append({
                "precoder": label,
                "delta_sc": int(delta),
                "rho_abs_mean": float(rho_abs[delta]),
                "rho_complex_abs": float(rho_complex_abs[delta]),
            })

    write_csv(metric_rows, out_dir / "v_tradeoff_metrics.csv")
    write_csv(front_rows, out_dir / "pareto_front.csv")
    write_csv(front_rows, out_dir / "pareto_front_labeled.csv")
    write_csv(curve_rows, out_dir / "covariance_function_selected.csv")
    plot_tradeoff(metric_rows, front_rows, fig_dir / "experiment18_v_design_tradeoff_scatter.svg")
    plot_covariance_curve(curve_rows, fig_dir / "experiment18_v_design_selected_covariance.svg")

    labels_for_link = [winner_label, ref_cdd_label]
    nmse_rows = run_nmse_compare(args, grid, pdp, base_cov, designs, labels_for_link, out_dir, fig_dir)
    bler_rows: List[Dict[str, object]] = []
    try:
        bler_rows = run_bler_compare(args, grid, pdp, base_cov, designs, labels_for_link, out_dir, fig_dir)
    except Exception as exc:
        with (out_dir / "bler_error.txt").open("w", encoding="utf-8") as f:
            f.write(repr(exc) + "\n")
        print(f"[warn] BLER comparison skipped/failed: {exc!r}", flush=True)

    selected_rows = {
        str(row["precoder"]): row
        for row in metric_rows
        if str(row["precoder"]) in set(labels_for_link)
    }
    summary = {
        "out_dir": str(out_dir),
        "winner_piecewise": winner_label,
        "reference_cdd": ref_cdd_label,
        "piecewise_margin_log10_product_vs_cdd_at_same_or_better_bw": margin,
        "selected_metric_rows": selected_rows,
        "nmse_rows": nmse_rows,
        "bler_rows": bler_rows,
        "figures": {
            "tradeoff": str(fig_dir / "experiment18_v_design_tradeoff_scatter.svg"),
            "covariance": str(fig_dir / "experiment18_v_design_selected_covariance.svg"),
            "nmse": str(fig_dir / "experiment18_v_design_nmse.svg"),
            "bler": str(fig_dir / "experiment18_v_design_bler.svg"),
        },
    }
    with (out_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    with (out_dir / "run_config.json").open("w", encoding="utf-8") as f:
        cfg = vars(args).copy()
        cfg["dmrs_symbol_indices"] = dmrs_symbols
        json.dump(cfg, f, indent=2, ensure_ascii=False)
    print(out_dir)
    return out_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/v_design_piecewise_tradeoff")
    parser.add_argument("--fig-dir", default="docs/figures")
    parser.add_argument("--n-tx", type=int, default=8)
    parser.add_argument("--n-rx", type=int, default=1)
    parser.add_argument("--n-segments", type=int, default=8)
    parser.add_argument("--carrier-bandwidth-mhz", type=float, default=100.0)
    parser.add_argument("--scs-khz", type=int, default=30)
    parser.add_argument("--n-fft", type=int, default=4096)
    parser.add_argument("--n-prbs", type=int, default=48)
    parser.add_argument("--pdsch-n-symbols", type=int, default=10)
    parser.add_argument("--dmrs-symbols", default="2,7")
    parser.add_argument("--dmrs-spacing-sc", type=int, default=6)
    parser.add_argument("--delay-spread-ns", type=float, default=5.0)
    parser.add_argument("--max-delay-factor", type=float, default=8.0)
    parser.add_argument("--cdd-steps", default="1,2,3,4,5,6,7,8,10,12,16,24,32,48,64,96,128,192,256")
    parser.add_argument("--piecewise-slopes", default="0,1,2,4,6,8,12,16,24,32,48,64")
    parser.add_argument("--random-per-slope", type=int, default=3)
    parser.add_argument("--max-corr-delta", type=int, default=192)
    parser.add_argument("--nmse-snrs", default="-4,0,4,8")
    parser.add_argument("--nmse-trials", type=int, default=80)
    parser.add_argument("--bler-snrs", default="0,4,8")
    parser.add_argument("--bler-trials", type=int, default=20)
    parser.add_argument("--mcs-index", type=int, default=8)
    parser.add_argument("--ldpc-iterations", type=int, default=8)
    parser.add_argument("--llr-clip", type=float, default=50.0)
    parser.add_argument("--loading", type=float, default=1e-8)
    parser.add_argument("--seed", type=int, default=20260615)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
