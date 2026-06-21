from __future__ import annotations

import argparse
import csv
import itertools
import json
from datetime import datetime
from pathlib import Path
import sys
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ResourceConfig
from cdd_lls.phy.channel_tdl import make_exponential_pdp
from cdd_lls.phy.estimators import covariance_matrix
from cdd_lls.phy.resource_grid import ResourceGrid, build_resource_grid
from tools.run_v_design_piecewise_tradeoff import (
    VDesign,
    balanced_random_permutation_indices,
    cdd_design,
    covariance_metrics,
    design_from_phase,
    finite_float,
    matrix_metrics,
    pareto_front,
    phase_from_piecewise_delays,
    plot_tradeoff,
    stable_seed,
    write_csv,
)


LEGACY_DIR = Path(
    "outputs/v_design_piecewise_tradeoff/"
    "v_design_piecewise_tradeoff_20260615_231257"
)


def parse_int_list(text: str) -> List[int]:
    return [int(item.strip()) for item in str(text).split(",") if item.strip()]


def format_number(value: float) -> str:
    value = float(value)
    if value.is_integer():
        return str(int(value))
    return f"{value:g}".replace(".", "p")


def alphabet_key(values: Sequence[float]) -> Tuple[float, ...]:
    return tuple(round(float(value), 9) for value in values)


def build_alphabet_catalog(
    uniform_steps: Sequence[int],
    magnitude_grid: Sequence[int],
) -> List[Dict[str, object]]:
    catalog_by_key: Dict[Tuple[float, ...], Dict[str, object]] = {}
    for step in uniform_steps:
        values = float(step) * (np.arange(8, dtype=np.float64) - 3.5)
        key = alphabet_key(values)
        catalog_by_key.setdefault(key, {
            "source": "uniform_centered",
            "source_parameter": f"step={int(step)}",
            "values": list(key),
        })
    for magnitudes in itertools.combinations(sorted(set(int(x) for x in magnitude_grid)), 4):
        values = tuple(float(x) for x in (-magnitudes[3], -magnitudes[2], -magnitudes[1], -magnitudes[0], *magnitudes))
        key = alphabet_key(values)
        catalog_by_key.setdefault(key, {
            "source": "irregular_symmetric_integer",
            "source_parameter": "positive=" + ",".join(str(x) for x in magnitudes),
            "values": list(key),
        })
    catalog = sorted(
        catalog_by_key.values(),
        key=lambda row: (
            max(abs(float(x)) for x in row["values"]),
            tuple(float(x) for x in row["values"]),
        ),
    )
    for index, row in enumerate(catalog, start=1):
        row["alphabet_id"] = f"A{index:04d}"
    return catalog


def cyclic_indices(n_segments: int, n_tx: int) -> np.ndarray:
    if int(n_segments) > int(n_tx):
        raise ValueError("n_segments must not exceed n_tx for unique per-Tx delays.")
    base = np.arange(int(n_tx), dtype=np.int64)
    return np.vstack([np.roll(base, segment) for segment in range(int(n_segments))])


def build_named_design(
    grid: ResourceGrid,
    alphabet_row: Dict[str, object],
    n_segments: int,
    layout_kind: str,
    layout_index: int,
    base_seed: int,
) -> Tuple[VDesign, Dict[str, object]]:
    alphabet = np.asarray(alphabet_row["values"], dtype=np.float64)
    n_tx = int(len(alphabet))
    alphabet_id = str(alphabet_row["alphabet_id"])
    magnitudes = [format_number(value) for value in alphabet if value > 0]
    magnitude_code = "_".join(magnitudes)
    if layout_kind == "cyclic":
        layout_id = "LC"
        layout_seed = ""
        layout_code = layout_id
        indices = cyclic_indices(n_segments, n_tx)
    elif layout_kind == "randomized_cyclic":
        layout_id = f"LR{int(layout_index):02d}"
        seed = stable_seed("balanced-signed", alphabet_id, layout_index, base=int(base_seed))
        layout_seed = int(seed)
        layout_code = f"{layout_id}-S{int(seed)}"
        indices = balanced_random_permutation_indices(
            n_segments,
            n_tx,
            np.random.default_rng(seed),
        ).astype(np.int64)
    else:
        raise ValueError(f"Unsupported layout_kind={layout_kind}")
    delays = alphabet[indices]
    name = f"NS-{alphabet_id}-M{magnitude_code}-{layout_code}"
    phase = phase_from_piecewise_delays(grid, delays, continuous=True)
    design = design_from_phase(
        label=name,
        family="piecewise_balanced_signed",
        phase=phase,
        metadata={
            "alphabet_id": alphabet_id,
            "alphabet_source": str(alphabet_row["source"]),
            "alphabet_source_parameter": str(alphabet_row["source_parameter"]),
            "layout_id": layout_id,
            "layout_kind": layout_kind,
            "layout_seed": layout_seed,
            "n_segments": int(n_segments),
            "phase_continuous": 1,
            "per_tx_delays_unique_across_segments": 1,
            "signed_centered_delays": 1,
        },
    )
    manifest = {
        "design_name": name,
        "alphabet_id": alphabet_id,
        "alphabet_source": alphabet_row["source"],
        "alphabet_source_parameter": alphabet_row["source_parameter"],
        "delay_alphabet_samples": json.dumps(alphabet.tolist(), separators=(",", ":")),
        "layout_id": layout_id,
        "layout_kind": layout_kind,
        "layout_seed": layout_seed,
        "segment_tx_delay_samples": json.dumps(delays.tolist(), separators=(",", ":")),
        "matrix_orientation": "rows=segments, columns=Tx0..Tx7",
        "phase_recurrence": "phi[k+1]=phi[k]-2*pi*d[s(k),n]/Nfft",
    }
    return design, manifest


def read_csv(path: Path) -> List[Dict[str, object]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def label_front(
    rows: Sequence[Dict[str, object]],
    prefix: str,
) -> List[Dict[str, object]]:
    ordered = sorted(
        rows,
        key=lambda row: (
            -finite_float(row.get("coherence_bw_abs_0p5_sc"), default=-1.0),
            -finite_float(row.get("diversity_log10_product_norm"), default=-1e9),
        ),
    )
    labeled = []
    for index, row in enumerate(ordered, start=1):
        item = dict(row)
        item["pareto_id"] = f"{prefix}{index}"
        labeled.append(item)
    return labeled


def label_cdd_front(rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    ids_by_step = {1: "C0", 2: "C1", 3: "C2", 4: "C3", 5: "C4", 6: "C5", 7: "C6", 64: "C7"}
    labeled = []
    for row in rows:
        item = dict(row)
        step = int(str(item["precoder"]).split("=")[-1])
        item["pareto_id"] = ids_by_step.get(step, f"Cstep{step}")
        labeled.append(item)
    return labeled


def add_cdd_comparison(
    front_rows: Sequence[Dict[str, object]],
    cdd_rows: Sequence[Dict[str, object]],
) -> List[Dict[str, object]]:
    output = []
    for row in front_rows:
        bw = finite_float(row["coherence_bw_abs_0p5_sc"])
        eligible = [
            ref for ref in cdd_rows
            if finite_float(ref.get("coherence_bw_abs_0p5_sc"), default=-1.0) >= bw
        ]
        item = dict(row)
        if eligible:
            ref = max(eligible, key=lambda value: finite_float(value["diversity_log10_product_norm"]))
            item["reference_cdd"] = ref["precoder"]
            item["reference_cdd_log10_product"] = finite_float(ref["diversity_log10_product_norm"])
            item["log10_product_gain_vs_cdd"] = (
                finite_float(row["diversity_log10_product_norm"])
                - finite_float(ref["diversity_log10_product_norm"])
            )
        else:
            item["reference_cdd"] = ""
            item["reference_cdd_log10_product"] = ""
            item["log10_product_gain_vs_cdd"] = ""
        output.append(item)
    return output


def apply_band_limited_coherence_width(
    row: Dict[str, object],
    n_sc: int,
) -> None:
    raw = row.get("coherence_bw_abs_0p5_sc", "")
    if raw == "" or not np.isfinite(finite_float(raw)):
        row["coherence_bw_abs_0p5_sc"] = int(n_sc)
        row["coherence_bw_abs_0p5_right_censored"] = 1
    else:
        row["coherence_bw_abs_0p5_right_censored"] = 0


def select_representative_front(
    rows: Sequence[Dict[str, object]],
    target_bandwidths: Sequence[int],
) -> List[Dict[str, object]]:
    selected: List[Dict[str, object]] = []
    seen = set()
    for target in target_bandwidths:
        row = min(
            rows,
            key=lambda item: abs(
                finite_float(item["coherence_bw_abs_0p5_sc"]) - float(target)
            ),
        )
        pareto_id = str(row["pareto_id"])
        if pareto_id not in seen:
            selected.append(dict(row))
            seen.add(pareto_id)
    return selected


def make_front_for_plot(
    rows: Sequence[Dict[str, object]],
    labeled_new_ids: Iterable[str],
) -> List[Dict[str, object]]:
    allowed = set(str(value) for value in labeled_new_ids)
    output = []
    for row in rows:
        item = dict(row)
        family = str(item.get("family"))
        if family in {"piecewise_linear_continuous_legacy", "piecewise_linear_reset"}:
            item["pareto_id"] = ""
        elif family == "piecewise_balanced_signed" and str(item.get("pareto_id")) not in allowed:
            item["pareto_id"] = ""
        output.append(item)
    return output


def run(args: argparse.Namespace) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path(args.out) / f"v_design_balanced_slope_scan_{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_path = Path(args.figure)
    figure_path.parent.mkdir(parents=True, exist_ok=True)

    resource = ResourceConfig(
        carrier_bandwidth_mhz=100.0,
        scs_khz=30,
        n_fft=4096,
        n_prbs=48,
        pdsch_n_symbols=10,
        dmrs_symbol_indices=[2, 7],
        dmrs_spacing_sc=6,
        prg_size_rb=4,
    )
    grid = build_resource_grid(resource)
    sample_period_ns = 1e9 / (resource.n_fft * resource.scs_khz * 1e3)
    pdp = make_exponential_pdp(
        delay_spread_ns=5.0,
        sample_period_ns=sample_period_ns,
        max_delay_factor=8.0,
    )
    base_covariance = covariance_matrix(
        grid.subcarrier_indices,
        grid.subcarrier_indices,
        pdp,
        grid.n_fft,
    )

    catalog = build_alphabet_catalog(
        uniform_steps=parse_int_list(args.uniform_steps),
        magnitude_grid=parse_int_list(args.magnitude_grid),
    )
    if int(args.max_alphabets) > 0:
        catalog = catalog[: int(args.max_alphabets)]
    layouts = [("cyclic", 0)] + [
        ("randomized_cyclic", index)
        for index in range(int(args.random_layouts_per_alphabet))
    ]

    metric_rows: List[Dict[str, object]] = []
    manifests: List[Dict[str, object]] = []
    total = len(catalog) * len(layouts)
    completed = 0
    for alphabet_row in catalog:
        for layout_kind, layout_index in layouts:
            design, manifest = build_named_design(
                grid,
                alphabet_row,
                int(args.n_segments),
                layout_kind,
                layout_index,
                int(args.seed),
            )
            matrix_row = matrix_metrics(design, grid, int(args.n_segments))
            covariance_row, _, _ = covariance_metrics(
                design,
                base_covariance,
                int(args.max_corr_delta),
            )
            row: Dict[str, object] = {
                "precoder": design.label,
                "family": design.family,
                "n_tx": design.n_tx,
                "n_rx": 1,
                "n_sc": grid.n_sc,
                "n_segments": int(args.n_segments),
                "segment_len_sc": grid.n_sc // int(args.n_segments),
                "n_fft": grid.n_fft,
                "scs_khz": grid.scs_khz,
                "delay_spread_ns": 5.0,
                "sample_period_ns": sample_period_ns,
                "dmrs_spacing_sc": 6,
            }
            row.update(matrix_row)
            row.update(covariance_row)
            apply_band_limited_coherence_width(row, grid.n_sc)
            row.update({f"design_{key}": value for key, value in design.metadata.items()})
            metric_rows.append(row)
            manifests.append(manifest)
            completed += 1
            if int(args.progress_every) > 0 and completed % int(args.progress_every) == 0:
                print(f"[progress] {completed}/{total}", flush=True)

    front = pareto_front(
        metric_rows,
        "coherence_bw_abs_0p5_sc",
        "diversity_log10_product_norm",
    )
    labeled_front = label_front(front, "N")

    legacy_dir = Path(args.legacy_dir)
    legacy_metrics = read_csv(legacy_dir / "v_tradeoff_metrics.csv")
    legacy_pareto = read_csv(legacy_dir / "pareto_front_labeled.csv")
    cdd_rows: List[Dict[str, object]] = []
    for step in parse_int_list(args.cdd_steps):
        design = cdd_design(grid, 8, step)
        matrix_row = matrix_metrics(design, grid, int(args.n_segments))
        covariance_row, _, _ = covariance_metrics(
            design,
            base_covariance,
            int(args.max_corr_delta),
        )
        row: Dict[str, object] = {
            "precoder": design.label,
            "family": design.family,
            "n_tx": 8,
            "n_rx": 1,
            "n_sc": grid.n_sc,
            "n_segments": int(args.n_segments),
            "segment_len_sc": grid.n_sc // int(args.n_segments),
            "n_fft": grid.n_fft,
            "scs_khz": grid.scs_khz,
            "delay_spread_ns": 5.0,
            "sample_period_ns": sample_period_ns,
            "dmrs_spacing_sc": 6,
        }
        row.update(matrix_row)
        row.update(covariance_row)
        apply_band_limited_coherence_width(row, grid.n_sc)
        cdd_rows.append(row)
    cdd_front = label_cdd_front(
        pareto_front(
            cdd_rows,
            "coherence_bw_abs_0p5_sc",
            "diversity_log10_product_norm",
        )
    )
    labeled_front = add_cdd_comparison(labeled_front, cdd_rows)

    combined_rows: List[Dict[str, object]] = []
    for row in legacy_metrics:
        if row.get("family") == "cdd_linear_phase":
            continue
        item = dict(row)
        if item.get("family") == "piecewise_linear_continuous":
            item["family"] = "piecewise_linear_continuous_legacy"
        combined_rows.append(item)
    combined_rows.extend(cdd_rows)
    combined_rows.extend(metric_rows)

    combined_front: List[Dict[str, object]] = []
    legacy_counter = 0
    for row in legacy_pareto:
        if row.get("family") == "cdd_linear_phase":
            continue
        item = dict(row)
        if item.get("family") == "piecewise_linear_continuous":
            item["family"] = "piecewise_linear_continuous_legacy"
            legacy_counter += 1
            item["pareto_id"] = f"L{legacy_counter}"
        combined_front.append(item)
    combined_front.extend(cdd_front)
    combined_front.extend(labeled_front)
    representative_front = select_representative_front(
        labeled_front,
        [576, 315, 230, 155, 104, 78, 63, 45, 5],
    )
    plot_front = make_front_for_plot(
        combined_front,
        [row["pareto_id"] for row in representative_front],
    )

    write_csv(metric_rows, output_dir / "new_signed_balanced_metrics.csv")
    write_csv(labeled_front, output_dir / "new_signed_balanced_pareto.csv")
    write_csv(representative_front, output_dir / "new_signed_balanced_pareto_representatives.csv")
    write_csv(cdd_rows, output_dir / "recomputed_cdd_metrics.csv")
    write_csv(cdd_front, output_dir / "recomputed_cdd_pareto.csv")
    write_csv(manifests, output_dir / "new_signed_balanced_implementation_manifest.csv")
    with (output_dir / "new_signed_balanced_implementation_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifests, handle, indent=2, ensure_ascii=False)
    write_csv(combined_rows, output_dir / "combined_legacy_and_new_metrics.csv")
    write_csv(combined_front, output_dir / "combined_legacy_and_new_pareto.csv")
    with (output_dir / "run_config.json").open("w", encoding="utf-8") as handle:
        json.dump(vars(args), handle, indent=2, ensure_ascii=False)

    plot_tradeoff(combined_rows, plot_front, figure_path)
    zoom_rows = [
        row for row in combined_rows
        if finite_float(row.get("coherence_bw_abs_0p5_sc"), default=1e9) <= 180.0
        and finite_float(row.get("diversity_log10_product_norm"), default=-1e9) >= -6.0
    ]
    zoom_front = [
        row for row in plot_front
        if finite_float(row.get("coherence_bw_abs_0p5_sc"), default=1e9) <= 180.0
        and finite_float(row.get("diversity_log10_product_norm"), default=-1e9) >= -6.0
    ]
    zoom_path = figure_path.with_name(f"{figure_path.stem}_zoom{figure_path.suffix}")
    plot_tradeoff(
        zoom_rows,
        zoom_front,
        zoom_path,
        title="8TX/1RX V-matrix tradeoff: decision-region zoom",
    )
    print(output_dir)
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/v_design_balanced_slope_scan")
    parser.add_argument("--figure", default="docs/figures/experiment18_v_design_tradeoff_scatter.svg")
    parser.add_argument("--legacy-dir", default=str(LEGACY_DIR))
    parser.add_argument("--n-segments", type=int, default=8)
    parser.add_argument("--uniform-steps", default="1,2,3,4,6,8,12,16,24,32,48,64")
    parser.add_argument("--magnitude-grid", default="1,2,3,4,6,8,12,16,24,32,48,64")
    parser.add_argument("--random-layouts-per-alphabet", type=int, default=3)
    parser.add_argument("--cdd-steps", default="1,2,3,4,5,6,7,8,10,12,16,24,32,48,64,96,128,192,256")
    parser.add_argument("--max-alphabets", type=int, default=0)
    parser.add_argument("--max-corr-delta", type=int, default=575)
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260619)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
