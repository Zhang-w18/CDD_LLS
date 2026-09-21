"""Merge, audit, and plot the formal plan-032 mobility experiment."""

from __future__ import annotations

import json
import math
import shutil
import sys
import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import j0

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.run_plan025_delay_matched_tdl import read_csv_rows, write_csv_rows


RUN_ROOT = ROOT / "outputs/experiment032_tdl_mobility/20260906_main"
FORMAL_DIRS = [
    RUN_ROOT / "formal",
    RUN_ROOT / "formal_shard1",
    RUN_ROOT / "formal_shard2",
    RUN_ROOT / "formal_shard3",
]
FINAL = RUN_ROOT / "final"
FIGURES = ROOT / "docs/figures/result-032"
EXPECTED_SNR = [14.0 + 0.25 * index for index in range(25)] + [22.0, 24.0]
PRG = "A100_PRG_DFT8_6RB"
MRT = "A100_AGED_CSI_MRT_PRG6_SLOTS10"
TARGETS = (0.1, 0.01)


LABELS = {
    "A100_B0_QC": "B0_QC",
    "A100_AP_RMS_T1": "AP_RMS_T1",
    "A100_AP_TEPS_T1": "AP_TEPS_T1",
    "A100_AP_TU_NT": "AP_TU_NT",
    "A100_AP_TU_NTM1": "AP_TU_NTM1",
    "A100_AP_TALIAS_NT": "AP_TALIAS_NT",
    "A100_S0_SIDON": "S0_SIDON",
    "A100_AP_T2_01": "AP_T2_01",
    "A100_MEFF_T2_04": "MEFF_T2_04",
    "A100_MEFF_T2_06": "MEFF_T2_06",
    PRG: "transparent PRG6",
    MRT: "aged-CSI MRT PRG6",
}


def _write_rows(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_csv_rows(rows, path)


def _read_formal() -> tuple[list[dict], list[dict]]:
    points = []
    intervals = []
    snr_sets = []
    for directory in FORMAL_DIRS:
        current = [
            dict(row)
            for row in read_csv_rows(directory / "final/estimated_csi_bler_points.csv")
        ]
        current_intervals = [dict(row) for row in read_csv_rows(directory / "intervals.csv")]
        if not current or not current_intervals:
            raise RuntimeError(f"Incomplete formal directory: {directory}")
        points.extend(current)
        intervals.extend(current_intervals)
        snr_sets.append({float(row["snr_db"]) for row in current})
    for left in range(len(snr_sets)):
        for right in range(left + 1, len(snr_sets)):
            overlap = snr_sets[left] & snr_sets[right]
            if overlap:
                raise RuntimeError(f"Formal shard SNR overlap: {sorted(overlap)}")
    candidate_ids = set(LABELS)
    actual_snr = sorted(set().union(*snr_sets))
    if actual_snr != EXPECTED_SNR:
        raise RuntimeError(f"Formal SNR union changed: {actual_snr}")
    counts = defaultdict(int)
    for row in points:
        counts[(str(row["candidate_id"]), float(row["snr_db"]))] += 1
        if int(row["trials"]) < 1000:
            raise RuntimeError("Formal public-grid budget is below 1000 trials")
    expected_keys = {(candidate_id, snr) for candidate_id in candidate_ids for snr in EXPECTED_SNR}
    if set(counts) != expected_keys or set(counts.values()) != {1}:
        raise RuntimeError("Formal candidate/SNR grid is incomplete or duplicated")
    points.sort(key=lambda row: (str(row["candidate_id"]), float(row["snr_db"])))
    intervals.sort(
        key=lambda row: (
            str(row["candidate_id"]), float(row["snr_db"]), int(row["trial_start"])
        )
    )
    return points, intervals


def _target_crossing(rows: list[dict], target: float) -> dict | None:
    ordered = sorted(rows, key=lambda row: float(row["snr_db"]))
    for left, right in zip(ordered, ordered[1:]):
        p0 = float(left["bler"])
        p1 = float(right["bler"])
        if p0 >= target and p1 <= target and p0 > p1:
            x0 = float(left["snr_db"])
            x1 = float(right["snr_db"])
            y0 = math.log10(max(p0, 0.5 / int(left["trials"])))
            y1 = math.log10(max(p1, 0.5 / int(right["trials"])))
            value = x0 + (math.log10(target) - y0) * (x1 - x0) / (y1 - y0)
            return {
                "target_snr_db": value,
                "snr_lo_db": x0,
                "snr_hi_db": x1,
                "bler_lo_snr": p0,
                "bler_hi_snr": p1,
            }
    return None


def _static_targets(candidate_ids: set[str]) -> dict[tuple[str, float], float]:
    paths = [
        ROOT
        / "outputs/experiment028_csi_curves/20260803_main/curve_augmentation_v2/a100/final/estimated_csi_bler_points.csv",
        ROOT
        / "outputs/experiment028_csi_curves/20260803_main/transparent_prg_baselines/a100/final/estimated_csi_bler_points.csv",
        ROOT
        / "outputs/experiment028_csi_curves/20260803_main/transparent_prg_6rb_1pct_extension/a100/final/estimated_csi_bler_points.csv",
    ]
    rows = []
    for path in paths:
        rows.extend(dict(row) for row in read_csv_rows(path))
    grouped = defaultdict(list)
    for row in rows:
        candidate_id = str(row["candidate_id"])
        if candidate_id in candidate_ids:
            grouped[candidate_id].append(row)
    output = {}
    for candidate_id, values in grouped.items():
        unique = {}
        for row in values:
            key = float(row["snr_db"])
            if key not in unique or int(row["trials"]) > int(unique[key]["trials"]):
                unique[key] = row
        for target in TARGETS:
            crossing = _target_crossing(list(unique.values()), target)
            if crossing:
                output[(candidate_id, target)] = float(crossing["target_snr_db"])
    return output


def _target_tables(points: list[dict]) -> tuple[list[dict], list[dict]]:
    grouped = defaultdict(list)
    for row in points:
        grouped[str(row["candidate_id"])].append(row)
    current = {}
    summary = []
    static = _static_targets(set(LABELS) - {MRT})
    for candidate_id in LABELS:
        for target in TARGETS:
            crossing = _target_crossing(grouped[candidate_id], target)
            key = (candidate_id, target)
            current[key] = None if crossing is None else float(crossing["target_snr_db"])
            summary.append(
                {
                    "candidate_id": candidate_id,
                    "label": LABELS[candidate_id],
                    "target_bler": target,
                    "bracketed": crossing is not None,
                    "snr_lo_db": "" if crossing is None else crossing["snr_lo_db"],
                    "snr_hi_db": "" if crossing is None else crossing["snr_hi_db"],
                    "bler_lo_snr": "" if crossing is None else crossing["bler_lo_snr"],
                    "bler_hi_snr": "" if crossing is None else crossing["bler_hi_snr"],
                    "target_snr_db": "" if crossing is None else crossing["target_snr_db"],
                    "static_target_snr_db": static.get(key, ""),
                    "mobility_shift_db": (
                        ""
                        if crossing is None or key not in static
                        else float(crossing["target_snr_db"]) - static[key]
                    ),
                }
            )
    comparisons = []
    for candidate_id in LABELS:
        if candidate_id in {PRG, MRT}:
            continue
        for target in TARGETS:
            value = current[(candidate_id, target)]
            for baseline in (PRG, MRT):
                reference = current[(baseline, target)]
                comparisons.append(
                    {
                        "candidate_id": candidate_id,
                        "label": LABELS[candidate_id],
                        "baseline_id": baseline,
                        "baseline_label": LABELS[baseline],
                        "target_bler": target,
                        "both_bracketed": value is not None and reference is not None,
                        "candidate_target_snr_db": "" if value is None else value,
                        "baseline_target_snr_db": "" if reference is None else reference,
                        "candidate_minus_baseline_db": (
                            "" if value is None or reference is None else value - reference
                        ),
                    }
                )
    return summary, comparisons


def _load_flag_map(intervals: list[dict]) -> dict[tuple[str, float], np.ndarray]:
    grouped = defaultdict(list)
    for row in intervals:
        grouped[(str(row["candidate_id"]), float(row["snr_db"]))].append(row)
    output = {}
    for key, values in grouped.items():
        values.sort(key=lambda row: int(row["trial_start"]))
        arrays = []
        for row in values:
            path = Path(str(row["error_flags"]))
            path = path if path.is_absolute() else ROOT / path
            array = np.load(path).astype(bool)
            if len(array) != int(row["trials"]):
                raise RuntimeError(f"Flag length mismatch: {path}")
            arrays.append(array)
        output[key] = np.concatenate(arrays)
    return output


def _interpolate_pair(x0: float, x1: float, p0: float, p1: float,
                      target: float, trials: int) -> float | None:
    floor = 0.5 / float(trials)
    p0 = max(float(p0), floor)
    p1 = max(float(p1), floor)
    # The measured point estimate is already guarded by _target_crossing().
    # Bootstrap draws use that fixed measured bracket and may fluctuate just
    # outside the target; retaining monotone draws avoids conditioning the CI
    # on the random event that every resample still straddles the target.
    if p0 <= p1:
        return None
    y0 = math.log10(p0)
    y1 = math.log10(p1)
    return x0 + (math.log10(target) - y0) * (x1 - x0) / (y1 - y0)


def _add_bootstrap_intervals(
    comparisons: list[dict], summary: list[dict], intervals: list[dict], repeats: int = 1000
) -> None:
    flags = _load_flag_map(intervals)
    brackets = {
        (str(row["candidate_id"]), float(row["target_bler"])): row
        for row in summary
        if bool(row["bracketed"])
    }
    rng = np.random.default_rng(20260906)
    for comparison in comparisons:
        comparison["bootstrap_repeats"] = repeats
        comparison["bootstrap_valid_repeats"] = 0
        comparison["difference_ci95_lo_db"] = ""
        comparison["difference_ci95_hi_db"] = ""
        if not bool(comparison["both_bracketed"]):
            continue
        target = float(comparison["target_bler"])
        candidate_id = str(comparison["candidate_id"])
        baseline_id = str(comparison["baseline_id"])
        definitions = {}
        for identifier in (candidate_id, baseline_id):
            bracket = brackets[(identifier, target)]
            x0 = float(bracket["snr_lo_db"])
            x1 = float(bracket["snr_hi_db"])
            definitions[identifier] = (x0, x1)
        samples = []
        for _ in range(repeats):
            resample_by_snr = {}
            estimates = {}
            for identifier in (candidate_id, baseline_id):
                x0, x1 = definitions[identifier]
                probabilities = []
                for snr in (x0, x1):
                    array = flags[(identifier, snr)]
                    if snr not in resample_by_snr:
                        resample_by_snr[snr] = rng.integers(0, len(array), size=len(array))
                    probabilities.append(float(np.mean(array[resample_by_snr[snr]])))
                estimates[identifier] = _interpolate_pair(
                    x0, x1, probabilities[0], probabilities[1], target, len(array)
                )
            if estimates[candidate_id] is not None and estimates[baseline_id] is not None:
                samples.append(estimates[candidate_id] - estimates[baseline_id])
        comparison["bootstrap_valid_repeats"] = len(samples)
        if len(samples) >= int(0.8 * repeats):
            comparison["difference_ci95_lo_db"] = float(np.percentile(samples, 2.5))
            comparison["difference_ci95_hi_db"] = float(np.percentile(samples, 97.5))


def _channel_diagnostics() -> dict:
    rows = []
    for directory in FORMAL_DIRS:
        path = directory / "channel_correlation_intervals.jsonl"
        for line in path.read_text(encoding="utf-8").splitlines():
            if line.strip():
                rows.append(json.loads(line))
    old = np.asarray([row["old_to_current_first_complex_correlation"] for row in rows])
    edge = np.asarray([row["current_symbol0_to9_complex_correlation"] for row in rows])
    weights = np.asarray(
        [int(row["trial_end"]) - int(row["trial_start"]) + 1 for row in rows],
        dtype=np.float64,
    )
    old_mean = float(np.average(old[:, 0], weights=weights))
    edge_mean = float(np.average(edge[:, 0], weights=weights))
    symbol_duration = (1.0 + 288.0 / 4096.0) / 30000.0
    age = 140.0 * symbol_duration
    maximum_doppler = (60.0 / 3.6) * 3.5e9 / 299792458.0
    return {
        "formal_interval_records": len(rows),
        "formal_total_trials_across_snr_intervals": int(np.sum(weights)),
        "maximum_doppler_hz": maximum_doppler,
        "actual_csi_age_ms": age * 1e3,
        "theoretical_j0_old_to_current": float(j0(2.0 * math.pi * maximum_doppler * age)),
        "empirical_old_to_current_weighted_mean_real": old_mean,
        "empirical_old_to_current_weighted_std_real_across_intervals": float(
            np.sqrt(np.average((old[:, 0] - old_mean) ** 2, weights=weights))
        ),
        "theoretical_j0_current_symbol0_to9": float(
            j0(2.0 * math.pi * maximum_doppler * 9.0 * symbol_duration)
        ),
        "empirical_current_symbol0_to9_weighted_mean_real": edge_mean,
        "empirical_current_symbol0_to9_weighted_std_real_across_intervals": float(
            np.sqrt(np.average((edge[:, 0] - edge_mean) ** 2, weights=weights))
        ),
    }


def _styles(candidate_ids: list[str]) -> dict[str, dict]:
    colors = list(plt.get_cmap("tab20").colors)
    markers = ["o", "s", "^", "v", "D", "P", "X", "<", ">", "h", "*", "d"]
    styles = {}
    for index, candidate_id in enumerate(candidate_ids):
        styles[candidate_id] = {
            "color": colors[index % len(colors)],
            "linestyle": "--" if candidate_id == MRT else "-." if candidate_id == PRG else "-",
            "marker": markers[index % len(markers)],
        }
    return styles


def _plot(points: list[dict], field: str, ylabel: str, filename: str, log_y: bool) -> None:
    grouped = defaultdict(list)
    for row in points:
        if float(row["snr_db"]) <= 20.0:
            grouped[str(row["candidate_id"])].append(row)
    order = list(LABELS)
    styles = _styles(order)
    figure, axis = plt.subplots(figsize=(15.0, 9.0))
    for candidate_id in order:
        rows = sorted(grouped[candidate_id], key=lambda row: float(row["snr_db"]))
        x_values = np.asarray([float(row["snr_db"]) for row in rows])
        y_values = np.asarray([float(row[field]) for row in rows])
        style = styles[candidate_id]
        if log_y:
            trials = np.asarray([float(row["trials"]) for row in rows])
            plotting_values = np.maximum(y_values, 0.5 / trials)
            axis.plot(
                x_values,
                plotting_values,
                label=LABELS[candidate_id],
                color=style["color"],
                linestyle=style["linestyle"],
                marker=style["marker"],
                linewidth=2.5,
                markersize=8,
                markevery=1,
                zorder=3,
            )
        else:
            axis.plot(
                x_values,
                y_values,
                label=LABELS[candidate_id],
                linewidth=2.5,
                markersize=8,
                markevery=2,
                **style,
            )
    if log_y:
        axis.set_yscale("log")
        # The experiment targets 1% BLER. Keep the lower bracket visible while
        # excluding the few-error 0.1% region whose Monte Carlo granularity
        # would dominate the visual impression on a logarithmic axis.
        axis.set_ylim(5e-3, 1.1)
        axis.axhline(1e-2, color="0.25", linewidth=1.8, linestyle=":", zorder=0)
    axis.set_xlabel("SNR (dB)", fontsize=17)
    axis.set_ylabel(ylabel, fontsize=17)
    axis.tick_params(labelsize=14)
    axis.grid(True, which="both", alpha=0.28)
    axis.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.14),
        ncol=4,
        fontsize=13,
        frameon=False,
    )
    figure.subplots_adjust(
        bottom=0.29,
        left=0.10,
        right=0.98,
        top=0.98,
    )
    output = FINAL / filename
    figure.savefig(output, dpi=180)
    FIGURES.mkdir(parents=True, exist_ok=True)
    shutil.copy2(output, FIGURES / filename)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--skip-ce-plot",
        action="store_true",
        help="Update statistics and the BLER figure without rewriting the CE-NMSE figure.",
    )
    args = parser.parse_args()
    FINAL.mkdir(parents=True, exist_ok=True)
    points, intervals = _read_formal()
    _write_rows(points, FINAL / "estimated_csi_bler_points.csv")
    _write_rows(intervals, FINAL / "formal_intervals.csv")
    summary, comparisons = _target_tables(points)
    _add_bootstrap_intervals(comparisons, summary, intervals)
    _write_rows(summary, FINAL / "target_summary.csv")
    _write_rows(comparisons, FINAL / "baseline_comparisons.csv")
    diagnostics = _channel_diagnostics()
    (FINAL / "channel_correlation_summary.json").write_text(
        json.dumps(diagnostics, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    styles = _styles(list(LABELS))
    (FINAL / "curve_styles.json").write_text(
        json.dumps(styles, indent=2, ensure_ascii=False, default=list) + "\n",
        encoding="utf-8",
    )
    _plot(points, "bler", "Estimated-CSI BLER", "a100_v60_estimated_csi_bler.png", True)
    if not args.skip_ce_plot:
        _plot(points, "ce_nmse_mean_db", "CE NMSE (dB)", "a100_v60_ce_nmse.png", False)
    analysis = {
        "schema": "result032-analysis-v1",
        "formal_points": len(points),
        "candidate_count": len(LABELS),
        "snr_count": len(EXPECTED_SNR),
        "minimum_trials_per_point": min(int(row["trials"]) for row in points),
        "candidate_trials": sum(int(row["trials"]) for row in points),
        "channel_diagnostics": diagnostics,
        "target_summary": summary,
        "baseline_comparisons": comparisons,
    }
    (FINAL / "analysis_summary.json").write_text(
        json.dumps(analysis, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(
        f"[analyze] points={len(points)} "
        f"candidate_trials={sum(int(row['trials']) for row in points)} "
        f"output={FINAL}",
        flush=True,
    )


if __name__ == "__main__":
    main()
