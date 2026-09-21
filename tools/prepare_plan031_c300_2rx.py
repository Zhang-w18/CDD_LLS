"""Prepare plan-031 C300 4Tx/2Rx smoke, prescan, and formal configs."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import importlib.metadata
import json
import math
import platform
from pathlib import Path
import subprocess
import sys
from typing import Any

import numpy as np
import yaml

import analyze_plan031 as common


ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = (
    ROOT
    / "outputs"
    / "experiment031_pdcch_cdd"
    / "20260911_c300_4tx_2rx_2sym"
)
ALS = (1, 2, 4)
MODES = ("estimated", "ideal")
TARGETS = (0.10, 0.01)


def _source_config(al: int, mode: str) -> Path:
    infix = "_ideal" if mode == "ideal" else ""
    return ROOT / "configs" / f"pdcch_result031_c300{infix}_al{al}_formal.yaml"


def _destination_config(al: int, mode: str, stage: str) -> Path:
    return (
        ROOT
        / "configs"
        / f"pdcch_result031_c300_2rx_{mode}_al{al}_{stage}.yaml"
    )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_yaml(path: Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle) or {}
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a YAML mapping: {path}")
    return payload


def _write_yaml(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(
        yaml.safe_dump(payload, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )


def _base_config(al: int, mode: str) -> dict[str, Any]:
    config = copy.deepcopy(_read_yaml(_source_config(al, mode)))
    config["seed"] = 20260911
    config["batch_size"] = 50
    config["antenna"]["n_rx"] = 2
    for candidate in config["candidates"]:
        candidate.pop("snr_points_db", None)
    expected_ce = "ideal" if mode == "ideal" else "frequency_lmmse"
    actual_ce = str(config["receiver"]["channel_estimation"]).lower()
    if actual_ce != expected_ce:
        raise ValueError(
            f"AL{al} {mode}: expected channel_estimation={expected_ce}, got {actual_ce}"
        )
    return config


def _normalized_for_diff(payload: dict[str, Any]) -> dict[str, Any]:
    normalized = copy.deepcopy(payload)
    for key in (
        "output_dir",
        "seed",
        "batch_size",
        "random_stream_namespace",
        "simulation",
    ):
        normalized.pop(key, None)
    normalized["antenna"].pop("n_rx", None)
    for candidate in normalized["candidates"]:
        candidate.pop("snr_points_db", None)
    return normalized


def _assert_only_allowed_differences(
    source: dict[str, Any], destination: dict[str, Any], label: str
) -> None:
    if _normalized_for_diff(source) != _normalized_for_diff(destination):
        raise ValueError(f"Unexpected non-control config difference: {label}")
    if int(source["antenna"]["n_rx"]) != 1:
        raise ValueError(f"Source n_rx is not 1: {label}")
    if int(destination["antenna"]["n_rx"]) != 2:
        raise ValueError(f"Destination n_rx is not 2: {label}")


def prepare_initial() -> None:
    receipt: dict[str, Any] = {
        "schema": "plan031-c300-2rx-config-diff-v1",
        "seed": 20260911,
        "unique_physical_difference": "antenna.n_rx: 1 -> 2",
        "allowed_control_differences": [
            "output_dir",
            "seed",
            "batch_size",
            "random_stream_namespace",
            "simulation",
            "candidates[*].snr_points_db",
        ],
        "comparisons": [],
    }
    stage_settings = {
        "smoke": {
            "snr_points_db": [-12.0, 18.0],
            "min_trials_per_snr": 20,
            "target_errors": 1_000_000,
            "max_trials_per_snr": 20,
            "target_bler": 0.01,
            "save_trial_error_flags": True,
            "resume": True,
        },
        "prescan": {
            "snr_points_db": [float(value) for value in range(-12, 20, 2)],
            "min_trials_per_snr": 300,
            "target_errors": 1_000_000,
            "max_trials_per_snr": 300,
            "target_bler": 0.01,
            "save_trial_error_flags": True,
            "resume": True,
        },
    }
    for al in ALS:
        for mode in MODES:
            source_path = _source_config(al, mode)
            source = _read_yaml(source_path)
            for stage, simulation in stage_settings.items():
                config = _base_config(al, mode)
                config["output_dir"] = str(RUN_ROOT / stage / mode / f"al{al}")
                config["random_stream_namespace"] = f"c300-2rx-{mode}-{stage}-v1"
                config["simulation"] = copy.deepcopy(simulation)
                if stage == "smoke":
                    config["batch_size"] = 20
                destination_path = _destination_config(al, mode, stage)
                _assert_only_allowed_differences(
                    source, config, f"AL{al} {mode} {stage}"
                )
                _write_yaml(destination_path, config)
                receipt["comparisons"].append(
                    {
                        "aggregation_level": al,
                        "csi_mode": mode,
                        "stage": stage,
                        "source_config": source_path.relative_to(ROOT).as_posix(),
                        "source_sha256": _sha256(source_path),
                        "new_config": destination_path.relative_to(ROOT).as_posix(),
                        "new_sha256": _sha256(destination_path),
                        "recursive_comparison_after_allowed_fields_removed": "equal",
                    }
                )
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    (RUN_ROOT / "config_diff_receipt.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    packages = {}
    for name in ("numpy", "PyYAML", "tensorflow", "sionna", "matplotlib"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    git_head = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    git_dirty = bool(
        subprocess.run(
            ["git", "status", "--short"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
        ).stdout.strip()
    )
    environment = {
        "schema": "plan031-c300-2rx-environment-v1",
        "python": sys.version,
        "platform": platform.platform(),
        "packages": packages,
        "git_head": git_head,
        "working_tree_dirty": git_dirty,
    }
    (RUN_ROOT / "environment_receipt.json").write_text(
        json.dumps(environment, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _read_candidate_points(al: int, mode: str, candidate_id: str) -> list[dict[str, str]]:
    path = RUN_ROOT / "prescan" / mode / f"al{al}" / candidate_id / "bler_points.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return sorted(csv.DictReader(handle), key=lambda row: float(row["snr_db"]))


def _raw_bracket(rows: list[dict[str, str]], target: float) -> tuple[float, float]:
    for left, right in zip(rows[:-1], rows[1:]):
        if float(left["bler"]) >= target >= float(right["bler"]):
            return float(left["snr_db"]), float(right["snr_db"])
    raise ValueError("target is not double-sided bracketed")


def _target_estimate(rows: list[dict[str, str]], target: float) -> float:
    snr = np.asarray([float(row["snr_db"]) for row in rows], dtype=np.float64)
    probabilities = np.asarray(
        [(int(row["errors"]) + 0.5) / (int(row["trials"]) + 1.0) for row in rows],
        dtype=np.float64,
    )
    weights = np.asarray([int(row["trials"]) for row in rows], dtype=np.float64)
    fitted = common.isotonic_decreasing(probabilities, weights)
    estimate, _ = common.interpolate_target(snr, fitted, target)
    return float(estimate)


def _quarter_db_grid(estimate: float) -> list[float]:
    low = math.floor((estimate - 0.75) * 4.0) / 4.0
    high = math.ceil((estimate + 0.75) * 4.0) / 4.0
    count = int(round((high - low) * 4.0)) + 1
    return [round(low + 0.25 * index, 2) for index in range(count)]


def prepare_formal() -> None:
    receipt: dict[str, Any] = {
        "schema": "plan031-c300-2rx-prescan-freeze-v1",
        "seed": 20260911,
        "targets": list(TARGETS),
        "batch_size": 50,
        "formal_grid_rule": "union of 0.25 dB grids spanning estimate +/- 0.75 dB",
        "config_sha256": {},
        "config_diff_comparisons": [],
        "candidates": [],
    }
    failures: list[str] = []
    pending: list[tuple[int, str, Path, dict[str, Any], Path]] = []
    for al in ALS:
        for mode in MODES:
            source_path = _source_config(al, mode)
            source = _read_yaml(source_path)
            config = _base_config(al, mode)
            for candidate in config["candidates"]:
                candidate_id = str(candidate["candidate_id"])
                rows = _read_candidate_points(al, mode, candidate_id)
                grids: set[float] = set()
                target_receipts = []
                for target in TARGETS:
                    try:
                        bracket = _raw_bracket(rows, target)
                        estimate = _target_estimate(rows, target)
                    except ValueError:
                        failures.append(f"AL{al}:{mode}:{candidate_id}:{target:g}")
                        continue
                    grids.update(_quarter_db_grid(estimate))
                    target_receipts.append(
                        {
                            "target_bler": target,
                            "prescan_bracket_db": list(bracket),
                            "prescan_estimate_db": estimate,
                        }
                    )
                candidate["snr_points_db"] = sorted(grids)
                receipt["candidates"].append(
                    {
                        "aggregation_level": al,
                        "csi_mode": mode,
                        "candidate_id": candidate_id,
                        "targets": target_receipts,
                        "formal_snr_points_db": candidate["snr_points_db"],
                    }
                )
            config["output_dir"] = str(RUN_ROOT / "formal" / mode / f"al{al}")
            config["random_stream_namespace"] = f"c300-2rx-{mode}-formal-v1"
            config["simulation"] = {
                "snr_points_db": [0.0],
                "min_trials_per_snr": 10_000,
                "target_errors": 200,
                "max_trials_per_snr": 50_000,
                "target_bler": 0.01,
                "save_trial_error_flags": True,
                "resume": True,
                "progress_every_batches": 50,
            }
            _assert_only_allowed_differences(source, config, f"AL{al} {mode} formal")
            path = _destination_config(al, mode, "formal")
            pending.append((al, mode, path, config, source_path))
    if failures:
        raise RuntimeError(
            "Prescan extension is required at 2 dB spacing for: " + ", ".join(failures)
        )
    for al, mode, path, config, source_path in pending:
        _write_yaml(path, config)
        receipt["config_sha256"][path.relative_to(ROOT).as_posix()] = _sha256(path)
        receipt["config_diff_comparisons"].append(
            {
                "aggregation_level": al,
                "csi_mode": mode,
                "source_config": source_path.relative_to(ROOT).as_posix(),
                "source_sha256": _sha256(source_path),
                "new_config": path.relative_to(ROOT).as_posix(),
                "new_sha256": _sha256(path),
                "recursive_comparison_after_allowed_fields_removed": "equal",
            }
        )
    (RUN_ROOT / "prescan_freeze_receipt.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("initial", "formal"), required=True)
    args = parser.parse_args()
    if args.stage == "initial":
        prepare_initial()
    else:
        prepare_formal()


if __name__ == "__main__":
    main()
