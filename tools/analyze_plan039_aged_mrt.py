"""Validate, summarize, and plot the plan-039 aged-MRT formal curve."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt


VARIANT_ID = "aged_mrt_prg6__60kmh__csi_age40ms"
STYLE = {
    "label": "Aged MRT, 60 km/h, CSI age 40 ms",
    "color": "#9467bd",
    "marker": "D",
    "linestyle": "-",
}
STAGE1B_STYLES = {
    "b0__plan039_common_reference_pdp": ("B0_QC (3 km/h)", "#1f77b4", "o", "-"),
    "sidon_selected__plan039_common_reference_pdp": ("SIDON_SELECTED (3 km/h)", "#d62728", "s", "--"),
    "cycling__plan039_prg_common_reference_pdp": ("PRG6_CYCLING (3 km/h)", "#2ca02c", "^", "-."),
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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def crossing(rows: list[dict[str, object]], target: float) -> dict[str, object]:
    brackets: list[tuple[float, float, float]] = []
    for left, right in zip(rows, rows[1:]):
        x0, x1 = float(left["snr_db"]), float(right["snr_db"])
        y0, y1 = float(left["bler"]), float(right["bler"])
        if y0 >= target >= y1 and y0 > 0.0 and y1 > 0.0:
            value = x0 + (math.log(target) - math.log(y0)) / (math.log(y1) - math.log(y0)) * (x1 - x0)
            brackets.append((x0, x1, value))
    return {
        "target_bler": target,
        "status": "descriptive_raw_bracket" if brackets else "no_raw_bracket",
        "crossing_snr_db": brackets[0][2] if brackets else "",
        "bracket_low_snr_db": brackets[0][0] if brackets else "",
        "bracket_high_snr_db": brackets[0][1] if brackets else "",
        "bracket_count": len(brackets),
    }


def validate_and_summarize(input_dir: Path) -> tuple[list[dict[str, object]], dict[str, object]]:
    summary_path = input_dir / "aged_mrt_formal_summary.csv"
    trials_path = input_dir / "aged_mrt_formal_trial_metrics.csv"
    freeze_path = input_dir / "aged_mrt_formal_freeze.json"
    report_path = input_dir / "aged_mrt_formal_report.json"
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    report = json.loads(report_path.read_text(encoding="utf-8"))
    rows = read_csv(summary_path)
    trial_rows = read_csv(trials_path)

    expected_grid = [7.5 + 0.5 * index for index in range(26)]
    actual_grid = [float(row["snr_db"]) for row in rows]
    if report.get("status") != "COMPLETE" or actual_grid != expected_grid:
        raise ValueError("Formal report is incomplete or the frozen 26-point SNR grid does not match.")
    if any(row["variant_id"] != VARIANT_ID for row in rows):
        raise ValueError("Unexpected variant in aged-MRT summary.")
    if any(int(row["n_trials"]) != 1000 for row in rows) or len(trial_rows) != 26000:
        raise ValueError("Expected exactly 1,000 trials per SNR and 26,000 trial rows.")
    if sha256(summary_path) != report["summary_sha256"] or sha256(trials_path) != report["trial_metrics_sha256"]:
        raise ValueError("Formal summary/trial-metrics digest does not match the completion report.")

    output: list[dict[str, object]] = []
    for row in rows:
        trials = int(row["n_trials"])
        errors = int(row["tb_errors"])
        nmse = float(row["ce_nmse_sum"]) / trials
        output.append({
            "variant_id": row["variant_id"],
            "snr_db": float(row["snr_db"]),
            "trials": trials,
            "tb_errors": errors,
            "bler": errors / trials,
            "wilson95_low": float(row["wilson95_low"]),
            "wilson95_high": float(row["wilson95_high"]),
            "ce_nmse_mean_linear": nmse,
            "ce_nmse_mean_db": 10.0 * math.log10(nmse),
            "csi_age_ms": float(row["precoder_csi_age_ms"]),
            "replay_error_max": float(row["aged_csi_replay_error_max"]),
            "precoder_raw_power_min": float(row["precoder_raw_power_min"]),
            "precoder_raw_power_max": float(row["precoder_raw_power_max"]),
        })

    audit = {
        "schema": "plan039-aged-mrt-analysis-audit-v1",
        "status": "PASS",
        "formal_report_status": report["status"],
        "snr_point_count": len(output),
        "snr_min_db": output[0]["snr_db"],
        "snr_max_db": output[-1]["snr_db"],
        "total_trials": len(trial_rows),
        "total_tb_errors": sum(int(row["tb_errors"]) for row in output),
        "zero_error_point_count": sum(int(row["tb_errors"]) == 0 for row in output),
        "csi_age_ms_values": sorted({float(row["csi_age_ms"]) for row in output}),
        "replay_error_max": max(float(row["replay_error_max"]) for row in output),
        "precoder_raw_power_min": min(float(row["precoder_raw_power_min"]) for row in output),
        "precoder_raw_power_max": max(float(row["precoder_raw_power_max"]) for row in output),
        "prescan_run": report["prescan_run"],
        "adaptive_additions": report["adaptive_additions"],
        "summary_sha256": sha256(summary_path),
        "trial_metrics_sha256": sha256(trials_path),
        "freeze_sha256": sha256(freeze_path),
        "frozen_grid_sha256": freeze["snr_grid_sha256"],
    }
    return output, audit


def plot(rows: list[dict[str, object]], output: Path, preview: bool = False) -> None:
    size = (5.12, 6.4) if preview else (10.0, 8.0)
    fig, axes = plt.subplots(2, 1, figsize=size, sharex=True, constrained_layout=True)
    x = [float(row["snr_db"]) for row in rows]
    trials = [int(row["trials"]) for row in rows]
    raw_bler = [float(row["bler"]) for row in rows]
    plotted_bler = [max(value, 0.5 / count) for value, count in zip(raw_bler, trials)]
    low = [max(float(row["wilson95_low"]), 0.5 / count) for row, count in zip(rows, trials)]
    high = [float(row["wilson95_high"]) for row in rows]

    axes[0].fill_between(x, low, high, color=STYLE["color"], alpha=0.18, label="Wilson 95% interval")
    axes[0].plot(x, plotted_bler, color=STYLE["color"], marker=STYLE["marker"],
                 linestyle=STYLE["linestyle"], linewidth=2.5, markersize=8, label=STYLE["label"])
    zero_x = [snr for snr, value in zip(x, raw_bler) if value == 0.0]
    zero_y = [0.5 / trials[x.index(snr)] for snr in zero_x]
    if zero_x:
        axes[0].scatter(zero_x, zero_y, facecolors="none", edgecolors=STYLE["color"], marker="o",
                        s=80, linewidths=2.0, label="0 errors; plotted at 0.5/trials")
    axes[0].axhline(0.1, color="#666666", linewidth=1.5, linestyle=":")
    axes[0].axhline(0.01, color="#666666", linewidth=1.5, linestyle=":")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("TB BLER", fontsize=16)
    axes[0].set_title("Plan-039 aged-MRT closed-loop supplement", fontsize=16)
    axes[0].grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=12, loc="upper right")

    axes[1].plot(x, [float(row["ce_nmse_mean_db"]) for row in rows], color=STYLE["color"],
                 marker=STYLE["marker"], linestyle=STYLE["linestyle"], linewidth=2.5, markersize=8)
    axes[1].set_xlabel("SNR (dB)", fontsize=16)
    axes[1].set_ylabel("CE NMSE (dB)", fontsize=16)
    axes[1].grid(True, alpha=0.3)
    for axis in axes:
        axis.tick_params(labelsize=14)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_combined(
    stage1b_rows: list[dict[str, str]],
    aged_rows: list[dict[str, object]],
    output: Path,
    preview: bool = False,
) -> None:
    size = (5.12, 6.4) if preview else (10.0, 8.0)
    fig, axes = plt.subplots(2, 1, figsize=size, sharex=True, constrained_layout=True)
    for variant, (label, color, marker, linestyle) in STAGE1B_STYLES.items():
        selected = sorted(
            (row for row in stage1b_rows if row["variant_id"] == variant),
            key=lambda row: float(row["snr_db"]),
        )
        x = [float(row["snr_db"]) for row in selected]
        y = [max(float(row["bler"]), 0.5 / int(row["trials"])) for row in selected]
        axes[0].plot(x, y, label=label, color=color, marker=marker, linestyle=linestyle,
                     linewidth=2.5, markersize=8, markevery=max(1, len(x) // 12))
        axes[1].plot(x, [float(row["ce_nmse_mean_db"]) for row in selected], label=label,
                     color=color, marker=marker, linestyle=linestyle, linewidth=2.5,
                     markersize=8, markevery=max(1, len(x) // 12))

    aged_x = [float(row["snr_db"]) for row in aged_rows]
    aged_y = [max(float(row["bler"]), 0.5 / int(row["trials"])) for row in aged_rows]
    axes[0].plot(aged_x, aged_y, label=STYLE["label"], color=STYLE["color"],
                 marker=STYLE["marker"], linestyle=STYLE["linestyle"], linewidth=2.5,
                 markersize=8, markevery=2)
    axes[1].plot(aged_x, [float(row["ce_nmse_mean_db"]) for row in aged_rows],
                 label=STYLE["label"], color=STYLE["color"], marker=STYLE["marker"],
                 linestyle=STYLE["linestyle"], linewidth=2.5, markersize=8, markevery=2)

    axes[0].set_yscale("log")
    axes[0].set_ylabel("TB BLER", fontsize=16)
    axes[0].set_title("Plan-039: four measured curves (different scenario groups)", fontsize=16)
    axes[0].grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=11, loc="lower left")
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
    parser.add_argument("--stage1b-points", type=Path)
    args = parser.parse_args()

    rows, audit = validate_and_summarize(args.input)
    crossing_rows = [crossing(rows, target) for target in (0.1, 0.01)]
    args.output.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output / "aged_mrt_points.csv", rows)
    write_csv(args.output / "aged_mrt_crossings.csv", crossing_rows)
    (args.output / "aged_mrt_analysis_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )
    (args.output / "aged_mrt_plot_style.json").write_text(
        json.dumps({VARIANT_ID: STYLE}, indent=2) + "\n", encoding="utf-8"
    )
    plot(rows, args.figure_dir / "aged_mrt_bler_nmse.png")
    plot(rows, args.figure_dir / "aged_mrt_bler_nmse_preview_13cm.png", preview=True)
    if args.stage1b_points is not None:
        stage1b_rows = read_csv(args.stage1b_points)
        missing = set(STAGE1B_STYLES) - {row["variant_id"] for row in stage1b_rows}
        if missing:
            raise ValueError(f"Stage-1B points are missing variants: {sorted(missing)}")
        plot_combined(stage1b_rows, rows, args.figure_dir / "four_curves_bler_nmse.png")
        plot_combined(stage1b_rows, rows, args.figure_dir / "four_curves_bler_nmse_preview_13cm.png",
                      preview=True)
    print(json.dumps({"audit": audit, "crossings": crossing_rows}, indent=2))


if __name__ == "__main__":
    main()
