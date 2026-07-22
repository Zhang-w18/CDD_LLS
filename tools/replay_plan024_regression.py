"""Capture and compare the plan-024 compatibility regression for plan-025."""

from __future__ import annotations

import argparse
import csv
import importlib.metadata
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parents[1]
SNRS = "14.25,14.5,14.75,16.25,17.25"


def run_logged(command: List[str], log_path: Path) -> int:
    environment = os.environ.copy()
    environment.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
    environment.setdefault("MPLCONFIGDIR", str(ROOT / "outputs" / "experiment025_sionna_tdl_rmmse" / "20260722_main" / "mplcache"))
    completed = subprocess.run(
        command,
        cwd=ROOT,
        env=environment,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(
        completed.stdout + ("\n[stderr]\n" + completed.stderr if completed.stderr else ""),
        encoding="utf-8",
    )
    return int(completed.returncode)


def package_receipt() -> Dict[str, object]:
    packages = {}
    for name in ("numpy", "scipy", "matplotlib", "PyYAML", "tensorflow", "sionna", "pytest"):
        packages[name] = importlib.metadata.version(name)
    try:
        import tensorflow as tf
        devices = [str(device) for device in tf.config.list_physical_devices()]
    except Exception as exc:  # pragma: no cover - diagnostic path
        devices = [f"ERROR: {exc!r}"]
    return {
        "python": sys.version,
        "python_executable": sys.executable,
        "packages": packages,
        "tensorflow_devices": devices,
    }


def read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def compare_link(before: Path, after: Path) -> Dict[str, object]:
    left = {(row["id"], row["snr_db"]): row for row in read_rows(before)}
    right = {(row["id"], row["snr_db"]): row for row in read_rows(after)}
    result: Dict[str, object] = {"keys_equal": left.keys() == right.keys(), "rows": []}
    passed = bool(result["keys_equal"])
    for key in sorted(left.keys() & right.keys()):
        a, b = left[key], right[key]
        exact_fields = ("trials", "tb_errors")
        float_fields = ("bler", "ce_nmse_mean", "ce_nmse_mean_dB", "estimator_cond")
        exact = all(a[field] == b[field] for field in exact_fields)
        differences = {field: abs(float(a[field]) - float(b[field])) for field in float_fields}
        row_passed = exact and max(differences.values(), default=0.0) <= 1e-12
        passed = passed and row_passed
        result["rows"].append({
            "id": key[0],
            "snr_db": key[1],
            "exact_counts": exact,
            "absolute_differences": differences,
            "passed": row_passed,
        })
    result["passed"] = passed
    return result


def compare_segment(before: Path, after: Path) -> Dict[str, object]:
    if not before.exists() or not after.exists():
        return {"passed": None, "status": "blocked_missing_plan023_source_data"}
    left = {(row["bandwidth_prb"], row["id"]): row for row in read_rows(before)}
    right = {(row["bandwidth_prb"], row["id"]): row for row in read_rows(after)}
    passed = left.keys() == right.keys()
    max_difference = 0.0
    exact_fields = ("pilots_per_segment_min", "pilots_per_segment_max", "diagnostic_h1_pass")
    numeric_fields = ("r1_nmse_16dB", "r4_nmse_16dB", "gain_r4_vs_r1_dB")
    for key in left.keys() & right.keys():
        passed = passed and all(left[key][field] == right[key][field] for field in exact_fields)
        for field in numeric_fields:
            max_difference = max(max_difference, abs(float(left[key][field]) - float(right[key][field])))
    passed = passed and max_difference <= 1e-12
    return {"passed": passed, "max_absolute_difference": max_difference}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=("before", "after"), required=True)
    parser.add_argument(
        "--root", type=Path,
        default=ROOT / "outputs" / "experiment025_sionna_tdl_rmmse" / "20260722_main",
    )
    parser.add_argument("--skip-run", action="store_true")
    args = parser.parse_args()
    phase_dir = args.root / f"regression_{args.phase}"
    phase_dir.mkdir(parents=True, exist_ok=True)
    commands = {
        "tests": [sys.executable, "-m", "pytest", "-q"],
        "segment": [
            sys.executable, "tools/run_experiment024_segment_sidon_qc.py",
            "--stage", "segment", "--out", str(phase_dir),
        ],
        "link": [
            sys.executable, "tools/run_experiment024_segment_sidon_qc.py",
            "--stage", "link", "--snrs", SNRS, "--trials", "100",
            "--batch-size", "1", "--out", str(phase_dir),
        ],
    }
    statuses: Dict[str, object] = {}
    if not args.skip_run:
        statuses["tests_exit_code"] = run_logged(commands["tests"], phase_dir / "tests.log")
        source = ROOT / "outputs" / "track_b_pilot_scan" / "20260709_main" / "candidates_48prb.csv"
        if source.exists():
            statuses["segment_exit_code"] = run_logged(commands["segment"], phase_dir / "segment.log")
        else:
            statuses["segment_exit_code"] = None
            statuses["segment_status"] = f"blocked: missing {source.relative_to(ROOT)}"
            (phase_dir / "segment.log").write_text(str(statuses["segment_status"]), encoding="utf-8")
        statuses["link_exit_code"] = run_logged(commands["link"], phase_dir / "link.log")
    (phase_dir / "environment.json").write_text(
        json.dumps(package_receipt(), indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (phase_dir / "commands.json").write_text(
        json.dumps(commands, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    statuses["phase"] = args.phase
    if args.phase == "after":
        before = args.root / "regression_before"
        link_before = before / "sidon_qc_bler.csv"
        link_after = phase_dir / "sidon_qc_bler.csv"
        statuses["link_comparison"] = (
            compare_link(link_before, link_after)
            if link_before.exists() and link_after.exists()
            else {"passed": None, "status": "missing regression CSV"}
        )
        statuses["segment_comparison"] = compare_segment(
            before / "segment_nmse.csv", phase_dir / "segment_nmse.csv"
        )
    (phase_dir / "regression_summary.json").write_text(
        json.dumps(statuses, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(statuses, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
