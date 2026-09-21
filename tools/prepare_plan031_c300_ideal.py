"""Create the plan-031 C300 two-symbol ideal-CSI prescan and formal configs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import yaml

import analyze_plan031 as common


ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "outputs" / "experiment031_pdcch_cdd" / "20260910_c300_4tx_2sym_ideal_csi"
TARGETS = (0.10, 0.01)
ALS = (1, 2, 4)


def _source_config(al: int) -> Path:
    return ROOT / "configs" / f"pdcch_result031_c300_al{al}_prescan.yaml"


def _destination_config(al: int, stage: str) -> Path:
    return ROOT / "configs" / f"pdcch_result031_c300_ideal_al{al}_{stage}.yaml"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _ideal_candidates(source: list[dict]) -> list[dict]:
    candidates: list[dict] = []
    small_template: dict | None = None
    for candidate in source:
        item = dict(candidate)
        candidate_id = str(item["candidate_id"])
        if candidate_id.startswith("C300_SMALL_CDD_QSTEP0P25_"):
            if small_template is None:
                small_template = item
            continue
        candidates.append(item)
    if small_template is None:
        raise ValueError("The estimated-CSI source config has no small-CDD waveform.")
    small_template["candidate_id"] = "C300_SMALL_CDD_QSTEP0P25_IDEAL_CSI"
    small_template["receiver_covariance_mode"] = "matched_effective"
    small_template["label"] = "C300 small CDD ideal CSI"
    small_template["style"] = {
        **dict(small_template.get("style", {})),
        "color": "#17becf",
        "linestyle": "-.",
        "marker": "h",
    }
    candidates.append(small_template)
    return candidates


def _base_config(al: int) -> dict:
    with open(_source_config(al), "r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    config["seed"] = 20260910
    config["receiver"]["channel_estimation"] = "ideal"
    config["candidates"] = _ideal_candidates(config["candidates"])
    return config


def _write_yaml(path: Path, payload: dict) -> None:
    path.write_text(
        yaml.safe_dump(payload, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )


def prepare_prescan() -> None:
    for al in ALS:
        config = _base_config(al)
        config["output_dir"] = str(RUN_ROOT / "prescan" / f"al{al}")
        config["random_stream_namespace"] = "c300-ideal-prescan-v1"
        config["simulation"] = {
            "snr_points_db": [float(value) for value in range(-12, 20, 2)],
            "min_trials_per_snr": 300,
            "target_errors": 1_000_000,
            "max_trials_per_snr": 300,
            "target_bler": 0.01,
            "save_trial_error_flags": True,
            "resume": True,
        }
        _write_yaml(_destination_config(al, "prescan"), config)


def _read_candidate_points(al: int, candidate_id: str) -> list[dict[str, object]]:
    path = RUN_ROOT / "prescan" / f"al{al}" / candidate_id / "bler_points.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    with open(path, "r", encoding="utf-8", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle)]
    return sorted(rows, key=lambda row: float(row["snr_db"]))


def _raw_bracket(rows: list[dict[str, object]], target: float) -> tuple[float, float]:
    for left, right in zip(rows[:-1], rows[1:]):
        if float(left["bler"]) >= target >= float(right["bler"]):
            return float(left["snr_db"]), float(right["snr_db"])
    raise ValueError("target is not double-sided bracketed")


def _target_estimate(rows: list[dict[str, object]], target: float) -> float:
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
    receipt: dict[str, object] = {
        "schema": "plan031-c300-ideal-prescan-freeze-v1",
        "seed": 20260910,
        "targets": list(TARGETS),
        "formal_grid_rule": "union of 0.25 dB grids spanning estimate +/- 0.75 dB",
        "config_sha256": {},
        "candidates": [],
    }
    failures: list[str] = []
    for al in ALS:
        config = _base_config(al)
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            rows = _read_candidate_points(al, candidate_id)
            grids: set[float] = set()
            target_receipts = []
            for target in TARGETS:
                try:
                    bracket = _raw_bracket(rows, target)
                    estimate = _target_estimate(rows, target)
                except ValueError:
                    failures.append(f"AL{al}:{candidate_id}:{target:g}")
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
                    "candidate_id": candidate_id,
                    "targets": target_receipts,
                    "formal_snr_points_db": candidate["snr_points_db"],
                }
            )
        if failures:
            continue
        config["output_dir"] = str(RUN_ROOT / "formal" / f"al{al}")
        config["random_stream_namespace"] = "c300-ideal-formal-v1"
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
        destination = _destination_config(al, "formal")
        _write_yaml(destination, config)
        receipt["config_sha256"][destination.relative_to(ROOT).as_posix()] = _sha256(destination)
    if failures:
        raise RuntimeError(
            "Prescan extension is required at 2 dB spacing for: " + ", ".join(failures)
        )
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    (RUN_ROOT / "prescan_freeze_receipt.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("prescan", "formal"), required=True)
    args = parser.parse_args()
    if args.stage == "prescan":
        prepare_prescan()
    else:
        prepare_formal()


if __name__ == "__main__":
    main()
