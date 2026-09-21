"""Audit and analyze plan-037 formal PDSCH Tx-scaling results."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.run_plan025_delay_matched_tdl import read_csv_rows, write_csv_rows


STYLES = {
    "B0_QC": {"color": "#1f77b4", "marker": "o", "label": "B0 QC"},
    "S0_SIDON": {"color": "#ff7f0e", "marker": "s", "label": "Strict Sidon"},
    "TRANSPARENT_PRG6": {"color": "#2ca02c", "marker": "^", "label": "PRG6 cycling"},
}


def short_id(candidate_id: str) -> str:
    for key in sorted(STYLES, key=len, reverse=True):
        if candidate_id.endswith(key):
            return key
    raise ValueError(f"Unknown plan-037 candidate: {candidate_id}")


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
        key = (str(row["candidate_id"]), str(row["receiver"]), float(row["snr_db"]))
        groups.setdefault(key, []).append(row)
        flags = np.load(ROOT / row["error_flags"])
        if flags.shape != (int(row["trials"]),) or int(np.sum(flags)) != int(row["tb_errors"]):
            raise RuntimeError(f"Error-flag mismatch: {row['error_flags']}")
        checked_arrays += 1
        if row["receiver"] == "estimated":
            nmse = np.load(ROOT / row["ce_nmse_trial_path"])
            if nmse.shape != flags.shape or not np.all(np.isfinite(nmse)):
                raise RuntimeError(f"NMSE mismatch: {row['ce_nmse_trial_path']}")
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
    candidate_ids = sorted({str(row["candidate_id"]) for row in estimated})
    snrs = sorted({float(row["snr_db"]) for row in estimated})
    expected_points = len(candidate_ids) * len(snrs)
    if len(candidate_ids) != 3 or len(estimated) != expected_points or len(ideal) != expected_points:
        raise RuntimeError("Formal point table dimensions do not match 3 candidates x common SNR grid.")
    runtime = json.loads((scene / "runtime_diagnostics.json").read_text(encoding="utf-8"))
    payload = {
        "status": adaptive["status"], "scene": scene.name,
        "candidate_count": len(candidate_ids), "snr_point_count": len(snrs),
        "curve_points": {"estimated": len(estimated), "ideal": len(ideal)},
        "interval_rows": len(intervals), "checked_arrays": checked_arrays,
        "capped_endpoint_count": len(adaptive.get("capped_endpoints", [])),
        "execution_device": runtime["execution_device"],
        "logical_gpus": runtime["logical_gpus"],
        "peak_rss_bytes": runtime["process_peak_rss_bytes_sampled"],
    }
    analysis = scene / "analysis"
    analysis.mkdir(parents=True, exist_ok=True)
    (analysis / "data_audit.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload


def _plot_bler(rows: list[dict], scene_name: str, receiver: str, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 5.3), constrained_layout=True)
    exported = []
    for candidate_id in sorted({str(row["candidate_id"]) for row in rows}):
        points = sorted((row for row in rows if row["candidate_id"] == candidate_id), key=lambda row: float(row["snr_db"]))
        key = short_id(candidate_id); style = STYLES[key]
        x = np.asarray([float(row["snr_db"]) for row in points])
        y = np.asarray([max(float(row["bler"]), 0.5 / int(row["trials"])) for row in points])
        ax.semilogy(x, y, color=style["color"], marker=style["marker"], linewidth=2.5, markersize=8, markevery=2, label=style["label"])
        exported.extend({**row, "plot_bler": float(value)} for row, value in zip(points, y))
    ax.axhline(0.1, color="black", linestyle="--", linewidth=1.5, label="10% BLER")
    ax.axhline(0.01, color="black", linestyle=":", linewidth=1.5, label="1% BLER")
    ax.set_ylim(8e-5, 1.1)
    ax.set_xlabel("SNR per Rx branch (dB)", fontsize=16)
    ax.set_ylabel("BLER", fontsize=16)
    ax.set_title(f"4T1R, 60 km/h: {receiver} CSI", fontsize=16)
    ax.tick_params(labelsize=14); ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=11, ncol=2)
    fig.savefig(output, dpi=300); plt.close(fig)
    _write(output.with_suffix(".csv"), exported)


def _plot_nmse(rows: list[dict], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 5.3), constrained_layout=True)
    exported = []
    for candidate_id in sorted({str(row["candidate_id"]) for row in rows}):
        points = sorted((row for row in rows if row["candidate_id"] == candidate_id), key=lambda row: float(row["snr_db"]))
        key = short_id(candidate_id); style = STYLES[key]
        x = np.asarray([float(row["snr_db"]) for row in points])
        y = np.asarray([float(row["ce_nmse_mean_db"]) for row in points])
        lo = np.asarray([10.0 * math.log10(max(float(row["ce_nmse_ci95_lo"]), np.finfo(float).tiny)) for row in points])
        hi = np.asarray([10.0 * math.log10(float(row["ce_nmse_ci95_hi"])) for row in points])
        ax.plot(x, y, color=style["color"], marker=style["marker"], linewidth=2.5, markersize=8, markevery=2, label=style["label"])
        ax.fill_between(x, lo, hi, color=style["color"], alpha=0.15)
        exported.extend({**row, "ce_nmse_ci95_lo_db": float(a), "ce_nmse_ci95_hi_db": float(b)} for row, a, b in zip(points, lo, hi))
    ax.set_xlabel("SNR per Rx branch (dB)", fontsize=16)
    ax.set_ylabel("Data-RE CE NMSE (dB)", fontsize=16)
    ax.set_title("4T1R, 60 km/h: estimated-CSI CE NMSE", fontsize=16)
    ax.tick_params(labelsize=14); ax.grid(True, alpha=0.3); ax.legend(fontsize=11)
    fig.savefig(output, dpi=300); plt.close(fig)
    _write(output.with_suffix(".csv"), exported)


def plot_scene(scene: Path) -> list[Path]:
    analysis = scene / "analysis"
    estimated = _read(scene / "final" / "estimated_csi_bler_points.csv")
    ideal = _read(scene / "final" / "ideal_csi_bler_points.csv")
    outputs = [analysis / "estimated_csi_bler.png", analysis / "ideal_csi_bler.png", analysis / "estimated_csi_ce_nmse.png"]
    _plot_bler(estimated, scene.name, "estimated", outputs[0])
    _plot_bler(ideal, scene.name, "ideal", outputs[1])
    _plot_nmse(estimated, outputs[2])
    (analysis / "style_map.json").write_text(json.dumps(STYLES, indent=2) + "\n", encoding="utf-8")
    from PIL import Image
    dimensions = {}
    for path in outputs:
        with Image.open(path) as image:
            dimensions[path.name] = {"width": image.width, "height": image.height, "format": image.format, "bytes": path.stat().st_size}
    (analysis / "figure_audit.json").write_text(json.dumps(dimensions, indent=2) + "\n", encoding="utf-8")
    return outputs


def crossing(points: list[dict], target: float) -> dict:
    ordered = sorted(points, key=lambda row: float(row["snr_db"]))
    for left, right in zip(ordered, ordered[1:]):
        p0, p1 = float(left["bler"]), float(right["bler"])
        if p0 >= target and p1 <= target and (p0 > target or p1 < target):
            if p0 <= 0.0 or p1 <= 0.0 or math.isclose(p0, p1):
                raise RuntimeError(f"Invalid log-BLER bracket for target {target}: {p0}, {p1}")
            value = float(left["snr_db"]) + (math.log(target) - math.log(p0)) * (float(right["snr_db"]) - float(left["snr_db"])) / (math.log(p1) - math.log(p0))
            return {
                "target_bler": target, "crossing_snr_db": value,
                "bracket_low_snr_db": float(left["snr_db"]), "bracket_high_snr_db": float(right["snr_db"]),
                "bracket_low_bler": p0, "bracket_high_bler": p1,
                "bracket_low_errors": int(left["tb_errors"]), "bracket_high_errors": int(right["tb_errors"]),
                "bracket_low_trials": int(left["trials"]), "bracket_high_trials": int(right["trials"]),
            }
    raise RuntimeError(f"No bracket for target {target}")


def _load_trial_matrix(scene: Path, snr_db: float, curve_keys: list[tuple[str, str]]) -> np.ndarray:
    intervals = _read(scene / "intervals.csv")
    columns = []
    for candidate_id, receiver in curve_keys:
        rows = [row for row in intervals if row["candidate_id"] == candidate_id and row["receiver"] == receiver and math.isclose(float(row["snr_db"]), snr_db, abs_tol=1e-12)]
        rows.sort(key=lambda row: int(row["trial_start"]))
        columns.append(np.concatenate([np.load(ROOT / row["error_flags"]).astype(np.float64) for row in rows]))
    lengths = {len(column) for column in columns}
    if len(lengths) != 1:
        raise RuntimeError(f"Paired trial lengths differ at {snr_db:g} dB.")
    return np.stack(columns, axis=1)


def bootstrap(scene: Path, replicates: int) -> tuple[list[dict], list[dict], list[dict]]:
    rows = _read(scene / "final" / "estimated_csi_bler_points.csv") + _read(scene / "final" / "ideal_csi_bler_points.csv")
    curve_keys = sorted({(str(row["candidate_id"]), str(row["receiver"])) for row in rows})
    snrs = sorted({float(row["snr_db"]) for row in rows})
    point_lookup = {(str(row["candidate_id"]), str(row["receiver"]), float(row["snr_db"])): row for row in rows}
    rng = np.random.default_rng(20260921)
    boot_bler: dict[tuple[str, str, float], np.ndarray] = {}
    for snr_db in snrs:
        matrix = _load_trial_matrix(scene, snr_db, curve_keys)
        estimates = np.empty((replicates, len(curve_keys)))
        for start in range(0, replicates, 25):
            used = min(25, replicates - start)
            indices = rng.integers(0, len(matrix), size=(used, len(matrix)))
            estimates[start:start + used] = np.mean(matrix[indices], axis=1)
        for index, key in enumerate(curve_keys):
            boot_bler[(*key, snr_db)] = estimates[:, index]
    samples: dict[tuple[str, str, float], np.ndarray] = {}
    crossing_rows = []
    for candidate_id, receiver in curve_keys:
        for target in (0.10, 0.01):
            point_rows = [point_lookup[(candidate_id, receiver, snr)] for snr in snrs]
            estimate = crossing(point_rows, target)
            values = np.full(replicates, np.nan)
            for rep in range(replicates):
                boot_points = [{"snr_db": snr, "bler": boot_bler[(candidate_id, receiver, snr)][rep], "tb_errors": 0, "trials": 1} for snr in snrs]
                try:
                    values[rep] = crossing(boot_points, target)["crossing_snr_db"]
                except RuntimeError:
                    pass
            samples[(short_id(candidate_id), receiver, target)] = values
            finite = values[np.isfinite(values)]; valid = len(finite) >= 0.8 * replicates
            crossing_rows.append({
                "scenario_id": scene.name, "candidate_id": candidate_id, "candidate_short": short_id(candidate_id), "receiver": receiver,
                **estimate, "bootstrap_replicates": replicates, "bootstrap_finite": len(finite),
                "crossing_ci95_lo_db": float(np.quantile(finite, 0.025)) if valid else "",
                "crossing_ci95_hi_db": float(np.quantile(finite, 0.975)) if valid else "",
                "bootstrap_status": "valid" if valid else "insufficient_finite_replicates",
            })
    gains = []
    penalties = []
    point_cross = {(row["candidate_short"], row["receiver"], float(row["target_bler"])): float(row["crossing_snr_db"]) for row in crossing_rows}
    for receiver in ("estimated", "ideal"):
        for target in (0.10, 0.01):
            for candidate in ("B0_QC", "S0_SIDON"):
                values = samples[("TRANSPARENT_PRG6", receiver, target)] - samples[(candidate, receiver, target)]
                finite = values[np.isfinite(values)]; valid = len(finite) >= 0.8 * replicates
                gains.append({
                    "scenario_id": scene.name, "candidate_short": candidate, "reference_short": "TRANSPARENT_PRG6", "receiver": receiver, "target_bler": target,
                    "gain_definition": "SNR_PRG6_minus_SNR_candidate", "gain_db": point_cross[("TRANSPARENT_PRG6", receiver, target)] - point_cross[(candidate, receiver, target)],
                    "bootstrap_replicates": replicates, "bootstrap_finite": len(finite),
                    "gain_ci95_lo_db": float(np.quantile(finite, 0.025)) if valid else "", "gain_ci95_hi_db": float(np.quantile(finite, 0.975)) if valid else "",
                    "bootstrap_status": "valid" if valid else "insufficient_finite_replicates",
                })
    for candidate in STYLES:
        for target in (0.10, 0.01):
            values = samples[(candidate, "estimated", target)] - samples[(candidate, "ideal", target)]
            finite = values[np.isfinite(values)]; valid = len(finite) >= 0.8 * replicates
            penalties.append({
                "scenario_id": scene.name, "candidate_short": candidate, "target_bler": target,
                "penalty_definition": "SNR_estimated_minus_SNR_ideal", "penalty_db": point_cross[(candidate, "estimated", target)] - point_cross[(candidate, "ideal", target)],
                "bootstrap_replicates": replicates, "bootstrap_finite": len(finite),
                "penalty_ci95_lo_db": float(np.quantile(finite, 0.025)) if valid else "", "penalty_ci95_hi_db": float(np.quantile(finite, 0.975)) if valid else "",
                "bootstrap_status": "valid" if valid else "insufficient_finite_replicates",
            })
    analysis = scene / "analysis"
    _write(analysis / "target_crossings_bootstrap.csv", crossing_rows)
    _write(analysis / "target_gains_vs_prg6_bootstrap.csv", gains)
    _write(analysis / "estimated_minus_ideal_penalties_bootstrap.csv", penalties)
    return crossing_rows, gains, penalties


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scene", required=True, type=Path)
    parser.add_argument("--bootstrap-replicates", type=int, default=1000)
    args = parser.parse_args()
    if args.bootstrap_replicates < 1000:
        raise ValueError("Plan-037 requires at least 1000 bootstrap replicates.")
    scene = args.scene.resolve()
    audit = audit_scene(scene)
    figures = plot_scene(scene)
    crossings, gains, penalties = bootstrap(scene, args.bootstrap_replicates)
    print(json.dumps({"audit": audit, "figures": [str(path) for path in figures], "crossings": len(crossings), "gains": len(gains), "penalties": len(penalties)}, indent=2))


if __name__ == "__main__":
    main()
