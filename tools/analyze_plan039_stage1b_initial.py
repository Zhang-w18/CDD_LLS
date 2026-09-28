"""Merge and plot the currently completed plan-039 stage-1B batches.

This script is intentionally limited to batch directories that contain a
``batch_receipt.json``.  Incomplete/stale directories are therefore excluded.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


VARIANTS = {
    "b0__plan039_common_reference_pdp": ("B0_QC", "#1f77b4", "o", "-"),
    "sidon_selected__plan039_common_reference_pdp": ("SIDON_SELECTED", "#d62728", "s", "--"),
    "cycling__plan039_prg_common_reference_pdp": ("PRG6_CYCLING", "#2ca02c", "^", "-."),
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def wilson(errors: int, trials: int) -> tuple[float, float]:
    if trials == 0:
        return math.nan, math.nan
    z = 1.959963984540054
    p = errors / trials
    den = 1.0 + z * z / trials
    center = (p + z * z / (2.0 * trials)) / den
    half = z * math.sqrt(p * (1.0 - p) / trials + z * z / (4.0 * trials * trials)) / den
    return center - half, center + half


def completed_batch_rows(base: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for receipt in sorted((base / "batches").glob("*/batch_receipt.json")):
        summary = receipt.parent / "summary.csv"
        if not summary.is_file():
            raise FileNotFoundError(f"Receipt has no summary: {receipt}")
        rows.extend(read_csv(summary))
    return [row for row in rows if row["variant_id"] in VARIANTS]


def aggregate(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    groups: dict[tuple[str, float], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(row["variant_id"], float(row["snr_db"]))].append(row)
    output: list[dict[str, object]] = []
    for (variant, snr), items in sorted(groups.items(), key=lambda item: (item[0][0], item[0][1])):
        intervals = sorted((int(x["absolute_trial_start"]), int(x["absolute_trial_stop"])) for x in items)
        if any(stop <= start for start, stop in intervals):
            raise ValueError(f"Invalid interval for {variant} at {snr}: {intervals}")
        if any(intervals[i][0] < intervals[i - 1][1] for i in range(1, len(intervals))):
            raise ValueError(f"Overlapping intervals for {variant} at {snr}: {intervals}")
        trials = sum(int(x["n_trials"]) for x in items)
        errors = sum(int(x["tb_errors"]) for x in items)
        nmse_sum = sum(float(x["ce_nmse_sum"]) for x in items)
        low, high = wilson(errors, trials)
        output.append({
            "variant_id": variant,
            "label": VARIANTS[variant][0],
            "snr_db": snr,
            "trials": trials,
            "tb_errors": errors,
            "bler": errors / trials,
            "wilson95_low": low,
            "wilson95_high": high,
            "ce_nmse_mean_linear": nmse_sum / trials,
            "ce_nmse_mean_db": 10.0 * math.log10(nmse_sum / trials),
            "absolute_trial_intervals": ";".join(f"[{a},{b})" for a, b in intervals),
        })
    return output


def coverage(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output = []
    for variant in VARIANTS:
        selected = [row for row in rows if row["variant_id"] == variant]
        trials = [int(row["trials"]) for row in selected]
        output.append({
            "variant_id": variant,
            "label": VARIANTS[variant][0],
            "snr_point_count": len(selected),
            "snr_min_db": min(float(row["snr_db"]) for row in selected),
            "snr_max_db": max(float(row["snr_db"]) for row in selected),
            "trials_per_point_min": min(trials),
            "trials_per_point_max": max(trials),
            "total_trials": sum(trials),
            "points_with_1000_trials": sum(value == 1000 for value in trials),
            "points_with_2000_trials": sum(value == 2000 for value in trials),
        })
    return output


def crossings(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for variant in VARIANTS:
        selected = sorted((row for row in rows if row["variant_id"] == variant),
                          key=lambda row: float(row["snr_db"]))
        for target in (0.1, 0.01):
            brackets = []
            for left, right in zip(selected, selected[1:]):
                y0, y1 = float(left["bler"]), float(right["bler"])
                if y0 >= target >= y1 and y0 > 0.0 and y1 > 0.0:
                    x0, x1 = float(left["snr_db"]), float(right["snr_db"])
                    crossing = x0 + (math.log(target) - math.log(y0)) / (math.log(y1) - math.log(y0)) * (x1 - x0)
                    brackets.append((x0, x1, crossing))
            output.append({
                "variant_id": variant,
                "label": VARIANTS[variant][0],
                "target_bler": target,
                "status": "current_raw_bracket" if brackets else "no_current_bracket",
                "crossing_snr_db": brackets[0][2] if brackets else "",
                "bracket_low_snr_db": brackets[0][0] if brackets else "",
                "bracket_high_snr_db": brackets[0][1] if brackets else "",
                "bracket_count": len(brackets),
            })
    return output


def plot(rows: list[dict[str, object]], output: Path, preview: bool = False) -> None:
    size = (5.12, 6.4) if preview else (10.0, 8.0)
    fig, axes = plt.subplots(2, 1, figsize=size, sharex=True, constrained_layout=True)
    for variant, (label, color, marker, line) in VARIANTS.items():
        selected = sorted((r for r in rows if r["variant_id"] == variant), key=lambda r: float(r["snr_db"]))
        x = [float(r["snr_db"]) for r in selected]
        y = [max(float(r["bler"]), 0.5 / int(r["trials"])) for r in selected]
        axes[0].plot(x, y, label=label, color=color, marker=marker, linestyle=line,
                     linewidth=2.5, markersize=8, markevery=max(1, len(x) // 12))
        axes[1].plot(x, [float(r["ce_nmse_mean_db"]) for r in selected], label=label,
                     color=color, marker=marker, linestyle=line, linewidth=2.5,
                     markersize=8, markevery=max(1, len(x) // 12))
    axes[0].set_yscale("log")
    axes[0].set_ylabel("TB BLER", fontsize=16)
    axes[0].set_title("Plan-039 stage-1B: currently completed trials", fontsize=16)
    axes[0].grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=14)
    axes[1].set_xlabel("SNR (dB)", fontsize=16)
    axes[1].set_ylabel("CE NMSE (dB)", fontsize=16)
    axes[1].grid(True, alpha=0.3)
    for axis in axes:
        axis.tick_params(labelsize=14)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    args = parser.parse_args()
    merged = aggregate(completed_batch_rows(args.input))
    counts = coverage(merged)
    args.output.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output / "stage1b_current_points.csv", merged)
    write_csv(args.output / "stage1b_trial_coverage.csv", counts)
    write_csv(args.output / "stage1b_current_crossings.csv", crossings(merged))
    style = {key: {"label": value[0], "color": value[1], "marker": value[2], "linestyle": value[3]}
             for key, value in VARIANTS.items()}
    (args.output / "plot_style.json").write_text(json.dumps(style, indent=2) + "\n", encoding="utf-8")
    plot(merged, args.figure_dir / "stage1b_current_bler_nmse.png")
    plot(merged, args.figure_dir / "stage1b_current_bler_nmse_preview_13cm.png", preview=True)
    print(json.dumps(counts, indent=2))


if __name__ == "__main__":
    main()
