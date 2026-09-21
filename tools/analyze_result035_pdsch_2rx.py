"""Audit and plot plan-035 formal estimated/ideal BLER and CE NMSE curves."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.run_plan025_delay_matched_tdl import read_csv_rows, wilson, write_csv_rows


STYLES = {
    "B0_QC": {"color": "#1f77b4", "marker": "o"},
    "S0_SIDON": {"color": "#ff7f0e", "marker": "s"},
    "TRANSPARENT_PRG6": {"color": "#2ca02c", "marker": "^"},
    "AGED_MRT_PRG6": {"color": "#d62728", "marker": "D"},
    "SMALL_CDD_QSTEP0P25_TRANSPARENT": {"color": "#9467bd", "marker": "v"},
    "SMALL_CDD_QSTEP0P25_MATCHED": {"color": "#8c564b", "marker": "P"},
}


def short_id(candidate_id: str) -> str:
    for key in sorted(STYLES, key=len, reverse=True):
        if candidate_id.endswith(key):
            return key
    raise ValueError(f"Unknown candidate id: {candidate_id}")


def _read(path: Path) -> list[dict]:
    if not path.exists():
        raise FileNotFoundError(path)
    return [dict(row) for row in read_csv_rows(path)]


def _write(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_csv_rows(rows, path)


def audit_scene(scene: Path) -> dict:
    adaptive = json.loads((scene / "adaptive_status.json").read_text(encoding="utf-8"))
    if adaptive["status"] not in {"complete", "complete_with_caps"}:
        raise RuntimeError(f"Adaptive run is not complete: {adaptive['status']}")
    intervals = _read(scene / "intervals.csv")
    groups: dict[tuple[str, str, float], list[dict]] = {}
    checked_arrays = 0
    for row in intervals:
        key = (row["candidate_id"], row["receiver"], float(row["snr_db"]))
        groups.setdefault(key, []).append(row)
        flags = np.load(ROOT / row["error_flags"])
        if flags.shape != (int(row["trials"]),) or int(np.sum(flags)) != int(row["tb_errors"]):
            raise RuntimeError(f"Error-flag mismatch: {row['error_flags']}")
        checked_arrays += 1
        if row["receiver"] == "estimated":
            nmse = np.load(ROOT / row["ce_nmse_trial_path"])
            if nmse.shape != flags.shape or not np.all(np.isfinite(nmse)):
                raise RuntimeError(f"NMSE array mismatch: {row['ce_nmse_trial_path']}")
            if not math.isclose(float(np.sum(nmse)), float(row["ce_nmse_sum"]), rel_tol=1e-10, abs_tol=1e-10):
                raise RuntimeError(f"NMSE sum mismatch: {row['ce_nmse_trial_path']}")
            checked_arrays += 1
    for key, values in groups.items():
        values.sort(key=lambda row: int(row["trial_start"]))
        expected = 1
        for row in values:
            if int(row["trial_start"]) != expected:
                raise RuntimeError(f"Trial gap/overlap: {key}")
            expected = int(row["trial_end"]) + 1
    estimated = _read(scene / "final" / "estimated_csi_bler_points.csv")
    ideal = _read(scene / "final" / "ideal_csi_bler_points.csv")
    if len(estimated) != 150 or len(ideal) != 150 or len(groups) != 300:
        raise RuntimeError("Expected 6 candidates x 25 SNR x 2 receivers.")
    runtime = json.loads((scene / "runtime_diagnostics.json").read_text(encoding="utf-8"))
    if runtime["execution_device"] != "CPU" or runtime["logical_gpus"]:
        raise RuntimeError("Formal run was not CPU-only.")
    payload = {
        "status": "passed" if adaptive["status"] == "complete" else "passed_with_caps",
        "scene": scene.name,
        "curve_points": {"estimated": len(estimated), "ideal": len(ideal)},
        "interval_rows": len(intervals),
        "checked_arrays": checked_arrays,
        "adaptive_status": adaptive["status"],
        "capped_endpoint_count": len(adaptive.get("capped_endpoints", [])),
        "execution_device": runtime["execution_device"],
        "peak_rss_bytes": runtime["process_peak_rss_bytes_sampled"],
    }
    (scene / "analysis").mkdir(parents=True, exist_ok=True)
    (scene / "analysis" / "data_audit.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload


def _plot_bler(rows: list[dict], scene_name: str, receiver: str, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
    plot_rows = []
    for candidate_id in sorted({row["candidate_id"] for row in rows}):
        points = sorted((row for row in rows if row["candidate_id"] == candidate_id), key=lambda row: float(row["snr_db"]))
        key = short_id(candidate_id)
        x = [float(row["snr_db"]) for row in points]
        y = [max(float(row["bler"]), 0.5 / int(row["trials"])) for row in points]
        ax.semilogy(x, y, color=STYLES[key]["color"], marker=STYLES[key]["marker"], linewidth=2.5, markersize=8, markevery=2, label=key)
        for row, plotted in zip(points, y):
            plot_rows.append({**row, "plot_bler": plotted})
    ax.axhline(0.1, color="black", linestyle="--", linewidth=1.5)
    ax.axhline(0.01, color="black", linestyle=":", linewidth=1.5)
    ax.set_ylim(0.008, 1.1)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.25))
    ax.set(xlabel="SNR per Rx branch (dB)", ylabel="BLER", title=f"{scene_name}: {receiver} CSI")
    ax.grid(True, which="major", axis="both", alpha=0.35, linewidth=0.9)
    ax.grid(True, which="minor", axis="x", alpha=0.22, linewidth=0.6)
    ax.tick_params(labelsize=14)
    ax.xaxis.label.set_size(16); ax.yaxis.label.set_size(16); ax.title.set_size(16)
    ax.legend(fontsize=10, ncol=2)
    fig.savefig(output, dpi=300)
    plt.close(fig)
    _write(output.with_suffix(".csv"), plot_rows)


def _plot_nmse(rows: list[dict], scene_name: str, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
    plot_rows = []
    for candidate_id in sorted({row["candidate_id"] for row in rows}):
        points = sorted((row for row in rows if row["candidate_id"] == candidate_id), key=lambda row: float(row["snr_db"]))
        key = short_id(candidate_id)
        x = [float(row["snr_db"]) for row in points]
        y = [float(row["ce_nmse_mean_db"]) for row in points]
        ax.plot(x, y, color=STYLES[key]["color"], marker=STYLES[key]["marker"], linewidth=2.5, markersize=8, markevery=2, label=key)
        plot_rows.extend(points)
    ax.set(xlabel="SNR per Rx branch (dB)", ylabel="Data-RE CE NMSE (dB)", title=f"{scene_name}: estimated CSI CE NMSE")
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=14)
    ax.xaxis.label.set_size(16); ax.yaxis.label.set_size(16); ax.title.set_size(16)
    ax.legend(fontsize=10, ncol=2)
    fig.savefig(output, dpi=300)
    plt.close(fig)
    _write(output.with_suffix(".csv"), plot_rows)


def plot_scene(scene: Path) -> list[Path]:
    analysis = scene / "analysis"
    analysis.mkdir(parents=True, exist_ok=True)
    estimated = _read(scene / "final" / "estimated_csi_bler_points.csv")
    ideal = _read(scene / "final" / "ideal_csi_bler_points.csv")
    outputs = [analysis / "estimated_csi_bler.png", analysis / "ideal_csi_bler.png", analysis / "estimated_csi_ce_nmse.png"]
    _plot_bler(estimated, scene.name, "estimated", outputs[0])
    _plot_bler(ideal, scene.name, "ideal", outputs[1])
    _plot_nmse(estimated, scene.name, outputs[2])
    (analysis / "style_map.json").write_text(json.dumps(STYLES, indent=2) + "\n", encoding="utf-8")
    dimensions = {}
    from PIL import Image
    for path in outputs:
        with Image.open(path) as image:
            dimensions[path.name] = {"width": image.width, "height": image.height, "format": image.format, "bytes": path.stat().st_size}
    (analysis / "figure_audit.json").write_text(json.dumps(dimensions, indent=2) + "\n", encoding="utf-8")
    return outputs


def _crossing(points: list[dict], target: float) -> dict:
    ordered = sorted(points, key=lambda row: float(row["snr_db"]))
    for left, right in zip(ordered, ordered[1:]):
        p0, p1 = float(left["bler"]), float(right["bler"])
        if p0 >= target and p1 <= target and (p0 > target or p1 < target):
            if p0 <= 0.0 or p1 <= 0.0 or math.isclose(p0, p1):
                value = float("nan")
            else:
                value = float(left["snr_db"]) + (math.log(target) - math.log(p0)) * (float(right["snr_db"]) - float(left["snr_db"])) / (math.log(p1) - math.log(p0))
            return {"target_bler": target, "crossing_snr_db": value, "bracket_low_snr_db": float(left["snr_db"]), "bracket_high_snr_db": float(right["snr_db"]), "bracket_low_bler": p0, "bracket_high_bler": p1, "bracket_low_errors": int(left["tb_errors"]), "bracket_high_errors": int(right["tb_errors"]), "bracket_low_trials": int(left["trials"]), "bracket_high_trials": int(right["trials"])}
    raise RuntimeError(f"No bracket for target {target}")


def summarize_crossings(scene: Path) -> list[dict]:
    rows = _read(scene / "final" / "estimated_csi_bler_points.csv") + _read(scene / "final" / "ideal_csi_bler_points.csv")
    output = []
    for (candidate_id, receiver) in sorted({(row["candidate_id"], row["receiver"]) for row in rows}):
        points = [row for row in rows if row["candidate_id"] == candidate_id and row["receiver"] == receiver]
        for target in (0.10, 0.01):
            output.append({"scenario_id": scene.name, "candidate_id": candidate_id, "candidate_short": short_id(candidate_id), "receiver": receiver, **_crossing(points, target)})
    _write(scene / "analysis" / "target_crossings.csv", output)
    return output


def summarize_differences(scene: Path, crossings: list[dict]) -> list[dict]:
    lookup = {(row["candidate_short"], row["receiver"], float(row["target_bler"])): float(row["crossing_snr_db"]) for row in crossings}
    comparisons = []
    pairs = [("S0_SIDON", "B0_QC", "S0_SIDON minus B0_QC"), ("SMALL_CDD_QSTEP0P25_MATCHED", "SMALL_CDD_QSTEP0P25_TRANSPARENT", "small matched minus transparent"), ("AGED_MRT_PRG6", "TRANSPARENT_PRG6", "aged MRT minus transparent PRG6")]
    for receiver in ("estimated", "ideal"):
        for target in (0.10, 0.01):
            for candidate, reference, label in pairs:
                comparisons.append({"scenario_id": scene.name, "comparison": label, "receiver": receiver, "target_bler": target, "candidate_snr_db": lookup[(candidate, receiver, target)], "reference_snr_db": lookup[(reference, receiver, target)], "delta_snr_db": lookup[(candidate, receiver, target)] - lookup[(reference, receiver, target)]})
    for candidate in STYLES:
        for target in (0.10, 0.01):
            comparisons.append({"scenario_id": scene.name, "comparison": f"{candidate} estimated minus ideal", "receiver": "estimated-ideal", "target_bler": target, "candidate_snr_db": lookup[(candidate, "estimated", target)], "reference_snr_db": lookup[(candidate, "ideal", target)], "delta_snr_db": lookup[(candidate, "estimated", target)] - lookup[(candidate, "ideal", target)]})
    _write(scene / "analysis" / "target_snr_differences.csv", comparisons)
    return comparisons


def paired_bootstrap(scene: Path, crossings: list[dict], comparisons: list[dict], replicates: int = 500) -> None:
    intervals = _read(scene / "intervals.csv")
    curve_keys = sorted({(row["candidate_id"], row["receiver"]) for row in intervals})
    snrs = sorted({float(row["snr_db"]) for row in intervals})
    bootstrap_bler: dict[tuple[str, str, float], np.ndarray] = {}
    rng = np.random.default_rng(20260918)
    for snr_db in snrs:
        columns = []
        for candidate_id, receiver in curve_keys:
            values = [row for row in intervals if row["candidate_id"] == candidate_id and row["receiver"] == receiver and math.isclose(float(row["snr_db"]), snr_db, abs_tol=1e-12)]
            values.sort(key=lambda row: int(row["trial_start"]))
            columns.append(np.concatenate([np.load(ROOT / row["error_flags"]).astype(np.float64) for row in values]))
        matrix = np.stack(columns, axis=1)
        estimates = np.empty((replicates, len(curve_keys)), dtype=np.float64)
        for start in range(0, replicates, 25):
            used = min(25, replicates - start)
            indices = rng.integers(0, len(matrix), size=(used, len(matrix)))
            estimates[start:start + used] = np.mean(matrix[indices], axis=1)
        for index, (candidate_id, receiver) in enumerate(curve_keys):
            bootstrap_bler[(candidate_id, receiver, snr_db)] = estimates[:, index]
    crossing_samples: dict[tuple[str, str, float], np.ndarray] = {}
    for candidate_id, receiver in curve_keys:
        key = short_id(candidate_id)
        for target in (0.10, 0.01):
            samples = np.full(replicates, np.nan)
            for rep in range(replicates):
                points = [{"snr_db": snr, "bler": bootstrap_bler[(candidate_id, receiver, snr)][rep], "tb_errors": 0, "trials": 1} for snr in snrs]
                try:
                    samples[rep] = _crossing(points, target)["crossing_snr_db"]
                except RuntimeError:
                    pass
            crossing_samples[(key, receiver, target)] = samples
    crossing_rows = []
    for row in crossings:
        samples = crossing_samples[(row["candidate_short"], row["receiver"], float(row["target_bler"]))]
        finite = samples[np.isfinite(samples)]
        crossing_rows.append({**row, "bootstrap_replicates": replicates, "bootstrap_finite": len(finite), "crossing_ci95_lo_db": float(np.quantile(finite, 0.025)), "crossing_ci95_hi_db": float(np.quantile(finite, 0.975))})
    _write(scene / "analysis" / "target_crossings_bootstrap.csv", crossing_rows)
    difference_rows = []
    for row in comparisons:
        target = float(row["target_bler"])
        if row["receiver"] == "estimated-ideal":
            candidate = row["comparison"].removesuffix(" estimated minus ideal")
            samples = crossing_samples[(candidate, "estimated", target)] - crossing_samples[(candidate, "ideal", target)]
        else:
            labels = {"S0_SIDON minus B0_QC": ("S0_SIDON", "B0_QC"), "small matched minus transparent": ("SMALL_CDD_QSTEP0P25_MATCHED", "SMALL_CDD_QSTEP0P25_TRANSPARENT"), "aged MRT minus transparent PRG6": ("AGED_MRT_PRG6", "TRANSPARENT_PRG6")}
            candidate, reference = labels[row["comparison"]]
            samples = crossing_samples[(candidate, row["receiver"], target)] - crossing_samples[(reference, row["receiver"], target)]
        finite = samples[np.isfinite(samples)]
        difference_rows.append({**row, "bootstrap_replicates": replicates, "bootstrap_finite": len(finite), "delta_ci95_lo_db": float(np.quantile(finite, 0.025)), "delta_ci95_hi_db": float(np.quantile(finite, 0.975))})
    _write(scene / "analysis" / "target_snr_differences_bootstrap.csv", difference_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scene", required=True, type=Path)
    args = parser.parse_args()
    scene = args.scene.resolve()
    audit = audit_scene(scene)
    outputs = plot_scene(scene)
    crossings = summarize_crossings(scene)
    differences = summarize_differences(scene, crossings)
    paired_bootstrap(scene, crossings, differences)
    print(json.dumps({"audit": audit, "figures": [str(path) for path in outputs], "crossings": len(crossings), "differences": len(differences)}, indent=2))


if __name__ == "__main__":
    main()
