"""Plot the result-028 A100 subset with transparent PRG and CDD baselines.

This is a read-only analysis of the merged result-028 point CSV files.  It does
not run new link trials.  The three figures deliberately omit titles for direct
placement in presentation slides.
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys
from pathlib import Path
from typing import Sequence

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".runtime/mpl"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
FINAL = ROOT / "outputs/experiment028_csi_curves/20260803_main/curve_augmentation_v2/a100/final"
TRANSPARENT_FINAL = ROOT / "outputs/experiment028_csi_curves/20260803_main/transparent_prg_baselines/a100/final"
TRANSPARENT_PRG6_EXTENSION_FINAL = ROOT / "outputs/experiment028_csi_curves/20260803_main/transparent_prg_6rb_1pct_extension/a100/final"
TRANSPARENT_CDD_FINAL = ROOT / "outputs/experiment028_csi_curves/20260803_main/transparent_cdd_physical_covariance/a100/final"
SMALL_DELAY_FINAL = ROOT / "outputs/experiment028_csi_curves/20260803_main/transparent_cdd_small_delay/a100/final"
SMALL_DELAY_MATCHED_FINAL = ROOT / "outputs/experiment028_csi_curves/20260803_main/small_delay_cdd_matched/a100/final"
SMALL_DELAY_IDEAL_FINAL = ROOT / "outputs/experiment028_csi_curves/20260803_main/small_delay_cdd_ideal/a100/final"
FIGURES = ROOT / "docs/figures/result-028"
SECTION37_FINAL = ROOT / "outputs/experiment028_csi_curves/20260803_main/section37_receiver_comparison/a100/final"
MANIFEST = ROOT / "outputs/experiment028_csi_curves/20260803_main/a100_comb6/manifest/source_manifest.json"

DELAY_SELECTED = (
    ("delay set 1", "A100_B0_QC"),
    ("delay set 2", "A100_AP_RMS_T1"),
    ("delay set 3", "A100_AP_TU_NT"),
    ("delay set 4", "A100_S0_SIDON"),
    ("delay set 5", "A100_MEFF_T2_06"),
)
TRANSPARENT_SELECTED = (
    ("transparent PRG 4 RB", "A100_PRG_DFT8_4RB"),
    ("transparent PRG 6 RB", "A100_PRG_DFT8_6RB"),
)
TRANSPARENT_CDD_SELECTED = (
    ("delay set 1, transparent CDD", "A100_B0_QC_TRANSPARENT_CDD"),
    ("delay set 2, transparent CDD", "A100_AP_RMS_T1_TRANSPARENT_CDD"),
    ("delay set 3, transparent CDD", "A100_AP_TU_NT_TRANSPARENT_CDD"),
    ("delay set 4, transparent CDD", "A100_S0_SIDON_TRANSPARENT_CDD"),
    ("delay set 5, transparent CDD", "A100_MEFF_T2_06_TRANSPARENT_CDD"),
)
SMALL_DELAY_SELECTED = (
    ("small-delay transparent CDD", "A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD"),
)
SMALL_DELAY_MATCHED_SELECTED = (
    ("small-delay CDD, matched CE", "A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD"),
)
SMALL_DELAY_IDEAL_SELECTED = (
    ("small-delay CDD", "A100_SMALL_CDD_QSTEP0P25_IDEAL_CDD"),
)
SELECTED = (*DELAY_SELECTED, *TRANSPARENT_SELECTED)
ESTIMATED_SELECTED = (
    *SELECTED,
    *TRANSPARENT_CDD_SELECTED,
    *SMALL_DELAY_SELECTED,
    *SMALL_DELAY_MATCHED_SELECTED,
)
IDEAL_SELECTED = (*SELECTED, *SMALL_DELAY_IDEAL_SELECTED)
SECTION37_SELECTED = (
    ("large-delay CDD, transparent", "A100_AP_RMS_T1_TRANSPARENT_CDD"),
    ("large-delay CDD, non-transparent", "A100_AP_RMS_T1"),
    ("small-delay CDD, transparent", "A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD"),
    ("small-delay CDD, non-transparent", "A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD"),
    ("precoder cycling 6 RB, transparent", "A100_PRG_DFT8_6RB"),
)
SNR_MAX_DB = 16.0
X_MIN_DB = 13.75
X_MAX_DB = 16.25
ESTIMATED_X_MAX_DB = 21.25
SECTION37_X_MIN_DB = 13.75
SECTION37_X_MAX_DB = 20.25
SECTION37_SNR_MAX_BY_CANDIDATE = {
    "A100_PRG_DFT8_6RB": 17.0,
    "A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD": 18.5,
}
PLOT_FONT = 18
TICK_FONT = 16
LINE_WIDTH = 3.0
MARKER_SIZE = 9.0


def fft_sample_delays(
    delay_grid_coordinates: Sequence[float],
    n_fft: int = 4096,
    phase_denominator: int = 576,
) -> list[float]:
    """Express /K grid delays as equivalent physical delays in FFT samples."""
    scale = float(n_fft) / float(phase_denominator)
    return [float(value) * scale for value in delay_grid_coordinates]


def _read_rows(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def subset_rows(
    rows: Sequence[dict],
    selected_curves: Sequence[tuple[str, str]] = SELECTED,
    require_plot_included: bool = False,
    snr_max_db: float = SNR_MAX_DB,
) -> list[dict]:
    wanted = {candidate_id for _, candidate_id in selected_curves}
    selected = [
        dict(row)
        for row in rows
        if str(row["candidate_id"]) in wanted
        and float(row["snr_db"]) <= snr_max_db
        and (
            not require_plot_included
            or str(row.get("plot_included", "")).strip().lower() == "true"
        )
    ]
    selected.sort(key=lambda row: (str(row["candidate_id"]), float(row["snr_db"])))
    return selected


def _matplotlib_linestyle(value: object) -> object:
    if (
        isinstance(value, list)
        and len(value) == 2
        and isinstance(value[1], list)
    ):
        return (value[0], tuple(value[1]))
    return value


def _plot_metric(
    rows: Sequence[dict],
    styles: dict,
    selected_curves: Sequence[tuple[str, str]],
    metric: str,
    ylabel: str,
    output: Path,
    log_y: bool,
    x_min_db: float = X_MIN_DB,
    x_max_db: float = X_MAX_DB,
    legend_loc: str | None = None,
    legend_bbox_to_anchor: tuple[float, float] | None = None,
) -> None:
    by_id: dict[str, list[dict]] = {}
    for row in rows:
        by_id.setdefault(str(row["candidate_id"]), []).append(dict(row))
    figure, axis = plt.subplots(figsize=(12.8, 7.2))
    for label, candidate_id in selected_curves:
        values = sorted(by_id[candidate_id], key=lambda row: float(row["snr_db"]))
        x = [float(row["snr_db"]) for row in values]
        if metric == "bler":
            y = [
                float(row[metric])
                if float(row[metric]) > 0.0
                else 0.5 / float(row["trials"])
                for row in values
            ]
        else:
            y = [float(row[metric]) for row in values]
        style = styles[candidate_id]
        axis.plot(
            x,
            y,
            label=label,
            color=style["color"],
            linestyle=_matplotlib_linestyle(style["linestyle"]),
            marker=style["marker"],
            linewidth=LINE_WIDTH,
            markersize=MARKER_SIZE,
        )
    if log_y:
        axis.set_yscale("log")
    axis.set_xlim(x_min_db, x_max_db)
    axis.xaxis.set_major_locator(MultipleLocator(0.5))
    axis.xaxis.set_minor_locator(MultipleLocator(0.1))
    axis.set_xlabel("SNR (dB)", fontsize=PLOT_FONT)
    axis.set_ylabel(ylabel, fontsize=PLOT_FONT)
    axis.tick_params(axis="both", which="both", labelsize=TICK_FONT)
    axis.grid(True, axis="y", which="both", alpha=0.28)
    axis.grid(True, axis="x", which="major", alpha=0.40, linewidth=1.0)
    axis.grid(True, axis="x", which="minor", alpha=0.20, linewidth=0.7)
    legend_options = {
        "loc": "best",
        "ncol": 2 if len(selected_curves) <= 8 else 3,
        "fontsize": TICK_FONT if len(selected_curves) <= 8 else 12,
        "frameon": True,
    }
    if metric == "ce_nmse_mean_db":
        legend_options.update({"loc": "center right", "bbox_to_anchor": (0.99, 0.63)})
    if legend_loc is not None:
        legend_options.update({"loc": legend_loc})
        legend_options.pop("bbox_to_anchor", None)
    if legend_bbox_to_anchor is not None:
        legend_options["bbox_to_anchor"] = legend_bbox_to_anchor
    axis.legend(**legend_options)
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(figure)


def _numerical_rank(singular_values: np.ndarray, relative_tolerance: float = 1e-10) -> int:
    return int(np.sum(singular_values > singular_values[0] * relative_tolerance))


def _write_rows(rows: Sequence[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted(set().union(*(row.keys() for row in rows))))
        writer.writeheader()
        writer.writerows(rows)


def _plot_data_rows(
    rows: Sequence[dict], metric: str, selected_curves: Sequence[tuple[str, str]]
) -> list[dict]:
    labels = {candidate_id: label for label, candidate_id in selected_curves}
    output = []
    for row in rows:
        item = {
            "label": labels[str(row["candidate_id"])],
            "candidate_id": row["candidate_id"],
            "snr_db": row["snr_db"],
            "trials": row["trials"],
            metric: row[metric],
        }
        if metric == "bler":
            item.update(
                {
                    "tb_errors": row["tb_errors"],
                    "bler_wilson95_lo": row.get("bler_wilson95_lo", ""),
                    "bler_wilson95_hi": row.get("bler_wilson95_hi", ""),
                }
            )
        else:
            item.update(
                {
                    "ce_nmse_mean": row.get("ce_nmse_mean", ""),
                    "ce_nmse_standard_error": row.get("ce_nmse_standard_error", ""),
                    "ce_nmse_ci95_lo": row.get("ce_nmse_ci95_lo", ""),
                    "ce_nmse_ci95_hi": row.get("ce_nmse_ci95_hi", ""),
                }
            )
        output.append(item)
    return sorted(output, key=lambda row: (str(row["candidate_id"]), float(row["snr_db"])))


def export_diagnostics(estimated_rows: Sequence[dict]) -> None:
    from cdd_lls.phy.resource_grid import local_indices_for_subcarriers
    from tools import run_plan027_bler as bler027

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    all_candidates = bler027.scenario_candidates(manifest, "A100")
    wanted = {candidate_id for _, candidate_id in DELAY_SELECTED}
    candidates = [row for row in all_candidates if str(row["candidate_id"]) in wanted]
    _, grid, precoders, covariances = bler027.build_scene(manifest, candidates, "A100")
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    ls_noise_variance = (8.0 / (10.0 ** (SNR_MAX_DB / 10.0))) / 2.0
    measured = {
        str(row["candidate_id"]): float(row["ce_nmse_mean_db"])
        for row in estimated_rows
        if math.isclose(float(row["snr_db"]), SNR_MAX_DB)
    }
    by_id = {str(row["candidate_id"]): row for row in candidates}
    diagnostic_rows = []
    delay_rows = []
    for label, candidate_id in DELAY_SELECTED:
        candidate = by_id[candidate_id]
        covariance = np.asarray(covariances[candidate_id], dtype=np.complex128)
        r_pp = covariance[np.ix_(pilot_local, pilot_local)]
        r_dd = covariance[np.ix_(data_local, data_local)]
        r_dp = covariance[np.ix_(data_local, pilot_local)]
        full_singular = np.linalg.svd(covariance, compute_uv=False)
        pilot_singular = np.linalg.svd(r_pp, compute_uv=False)
        full_rank = _numerical_rank(full_singular)
        pilot_covariance_rank = _numerical_rank(pilot_singular)
        pilot_matrix = np.asarray(precoders[candidate_id].C)[pilot_local]
        pilot_matrix_singular = np.linalg.svd(pilot_matrix, compute_uv=False)
        pilot_matrix_rank = _numerical_rank(pilot_matrix_singular)
        pilot_matrix_condition = float(
            candidate.get(
                "pilot_condition_number",
                pilot_matrix_singular[0] / max(pilot_matrix_singular[-1], 1e-300),
            )
        )
        noise_free_error = r_dd - r_dp @ np.linalg.pinv(r_pp, rcond=1e-10) @ r_dp.conj().T
        finite_error = r_dd - r_dp @ np.linalg.solve(
            r_pp + ls_noise_variance * np.eye(len(pilot_local)), r_dp.conj().T
        )
        denominator = float(np.real(np.trace(r_dd)))
        noise_free_nmse = max(float(np.real(np.trace(noise_free_error))) / denominator, 1e-30)
        finite_nmse = float(np.real(np.trace(finite_error))) / denominator
        diagnostic_rows.append(
            {
                "candidate_label": label,
                "candidate_id": candidate_id,
                "pilot_matrix_rank": pilot_matrix_rank,
                "pilot_matrix_condition": pilot_matrix_condition,
                "pilot_residues": json.dumps(
                    candidate.get("residues", []), separators=(",", ":")
                ),
                "full_covariance_rank": full_rank,
                "pilot_covariance_rank": pilot_covariance_rank,
                "rank_gap": full_rank - pilot_covariance_rank,
                "noise_free_nmse_floor_db": 10.0 * math.log10(noise_free_nmse),
                "theoretical_nmse_16db": 10.0 * math.log10(finite_nmse),
                "monte_carlo_nmse_16db": measured[candidate_id],
                "rank_relative_tolerance": 1e-10,
            }
        )
        q = [float(value) for value in candidate["delay_grid_coordinates"]]
        ns = [value / (576.0 * 30e3) * 1e9 for value in q]
        delay_rows.append(
            {
                "candidate_label": label,
                "candidate_id": candidate_id,
                "delay_grid_coordinates": json.dumps(q, separators=(",", ":")),
                "delay_ns": json.dumps(ns, separators=(",", ":")),
                "delay_fft_samples": json.dumps(
                    fft_sample_delays(q), separators=(",", ":")
                ),
            }
        )
    _write_rows(diagnostic_rows, FINAL / "a100_subset_matrix_diagnostics.csv")
    _write_rows(delay_rows, FINAL / "a100_subset_delay_table.csv")


def main() -> None:
    estimated = subset_rows(_read_rows(FINAL / "estimated_csi_bler_points.csv"))
    estimated.extend(subset_rows(_read_rows(TRANSPARENT_FINAL / "estimated_csi_bler_points.csv")))
    estimated.extend(
        subset_rows(
            _read_rows(TRANSPARENT_PRG6_EXTENSION_FINAL / "estimated_csi_bler_points.csv"),
            (("transparent PRG 6 RB", "A100_PRG_DFT8_6RB"),),
            snr_max_db=18.0,
        )
    )
    estimated.extend(
        subset_rows(
            _read_rows(TRANSPARENT_CDD_FINAL / "estimated_csi_bler_points.csv"),
            TRANSPARENT_CDD_SELECTED,
        )
    )
    estimated.extend(
        subset_rows(
            _read_rows(SMALL_DELAY_FINAL / "estimated_csi_bler_points.csv"),
            SMALL_DELAY_SELECTED,
            snr_max_db=21.0,
        )
    )
    estimated.extend(
        subset_rows(
            _read_rows(SMALL_DELAY_MATCHED_FINAL / "estimated_csi_bler_points.csv"),
            SMALL_DELAY_MATCHED_SELECTED,
            snr_max_db=20.0,
        )
    )
    ideal = subset_rows(_read_rows(FINAL / "ideal_csi_bler_points.csv"), require_plot_included=True)
    ideal.extend(subset_rows(_read_rows(TRANSPARENT_FINAL / "ideal_csi_bler_points.csv")))
    ideal.extend(
        subset_rows(
            _read_rows(SMALL_DELAY_IDEAL_FINAL / "ideal_csi_bler_points.csv"),
            SMALL_DELAY_IDEAL_SELECTED,
            snr_max_db=19.0,
        )
    )
    styles = json.loads((FINAL / "curve_styles.json").read_text(encoding="utf-8"))
    styles.update(
        json.loads((TRANSPARENT_FINAL / "curve_styles.json").read_text(encoding="utf-8"))
    )
    styles.update(
        json.loads(
            (TRANSPARENT_PRG6_EXTENSION_FINAL / "curve_styles.json").read_text(encoding="utf-8")
        )
    )
    styles.update(
        json.loads((TRANSPARENT_CDD_FINAL / "curve_styles.json").read_text(encoding="utf-8"))
    )
    styles.update(
        json.loads((SMALL_DELAY_FINAL / "curve_styles.json").read_text(encoding="utf-8"))
    )
    styles.update(
        json.loads((SMALL_DELAY_MATCHED_FINAL / "curve_styles.json").read_text(encoding="utf-8"))
    )
    styles.update(
        json.loads((SMALL_DELAY_IDEAL_FINAL / "curve_styles.json").read_text(encoding="utf-8"))
    )
    _write_rows(
        _plot_data_rows(estimated, "bler", ESTIMATED_SELECTED),
        SMALL_DELAY_FINAL / "a100_subset_with_small_delay_transparent_cdd_estimated_points.csv",
    )
    _write_rows(
        _plot_data_rows(ideal, "bler", IDEAL_SELECTED),
        TRANSPARENT_FINAL / "a100_subset_with_transparent_ideal_points.csv",
    )
    _write_rows(
        _plot_data_rows(estimated, "ce_nmse_mean_db", ESTIMATED_SELECTED),
        SMALL_DELAY_FINAL / "a100_subset_with_small_delay_transparent_cdd_ce_points.csv",
    )
    _plot_metric(
        estimated,
        styles,
        ESTIMATED_SELECTED,
        "bler",
        "Estimated-CSI BLER",
        FIGURES / "a100_comb6_subset_estimated_csi_bler_snr_le_21.png",
        True,
        x_max_db=ESTIMATED_X_MAX_DB,
    )
    _plot_metric(
        ideal,
        styles,
        IDEAL_SELECTED,
        "bler",
        "Ideal-CSI BLER",
        FIGURES / "a100_comb6_subset_ideal_csi_bler_snr_le_19.png",
        True,
        x_max_db=19.25,
    )
    _plot_metric(
        estimated,
        styles,
        ESTIMATED_SELECTED,
        "ce_nmse_mean_db",
        "CE NMSE (dB)",
        FIGURES / "a100_comb6_subset_ce_nmse_snr_le_21.png",
        False,
        x_max_db=ESTIMATED_X_MAX_DB,
    )
    section37 = subset_rows(
        estimated,
        SECTION37_SELECTED,
        snr_max_db=20.0,
    )
    section37 = [
        row
        for row in section37
        if float(row["snr_db"])
        <= SECTION37_SNR_MAX_BY_CANDIDATE.get(str(row["candidate_id"]), 20.0)
    ]
    _write_rows(
        _plot_data_rows(section37, "bler", SECTION37_SELECTED),
        SECTION37_FINAL / "a100_section37_estimated_csi_bler_points.csv",
    )
    _write_rows(
        _plot_data_rows(section37, "ce_nmse_mean_db", SECTION37_SELECTED),
        SECTION37_FINAL / "a100_section37_ce_nmse_points.csv",
    )
    _plot_metric(
        section37,
        styles,
        SECTION37_SELECTED,
        "bler",
        "Estimated-CSI BLER",
        FIGURES / "a100_comb6_section37_estimated_csi_bler_snr_14_20.png",
        True,
        x_min_db=SECTION37_X_MIN_DB,
        x_max_db=SECTION37_X_MAX_DB,
        legend_loc="upper right",
        legend_bbox_to_anchor=(0.99, 0.88),
    )
    _plot_metric(
        section37,
        styles,
        SECTION37_SELECTED,
        "ce_nmse_mean_db",
        "CE NMSE (dB)",
        FIGURES / "a100_comb6_section37_ce_nmse_snr_14_20.png",
        False,
        x_min_db=SECTION37_X_MIN_DB,
        x_max_db=SECTION37_X_MAX_DB,
    )
    export_diagnostics(estimated)


if __name__ == "__main__":
    main()
