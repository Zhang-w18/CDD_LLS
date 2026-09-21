"""Freeze 0.5 dB/3000-trial coarse-confirmation configs from the completed prescan."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search"
PRESCAN_ROOT = RUN_ROOT / "prescan"
SCREENING_CSV = RUN_ROOT / "screening/candidate_screening.csv"
RECEIPT_PATH = RUN_ROOT / "coarse_confirmation_config_receipt.json"
TARGETS = (0.10, 0.01)
TRIALS_PER_POINT = 3000
SOURCE_CONFIGS = {
    "a100_al2": ROOT / "configs/pdcch_result031_strict_sidon_a100_al2_prescan.yaml",
    "a100_al4": ROOT / "configs/pdcch_result031_strict_sidon_a100_al4_prescan.yaml",
    "a100_al8": ROOT / "configs/pdcch_result031_strict_sidon_a100_al8_prescan.yaml",
    "c300_al1": ROOT / "configs/pdcch_result031_strict_sidon_c300_al1_prescan.yaml",
    "c300_al2": ROOT / "configs/pdcch_result031_strict_sidon_c300_al2_prescan.yaml",
    "c300_al4": ROOT / "configs/pdcch_result031_strict_sidon_c300_al4_prescan.yaml",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _first_raw_bracket(rows: list[dict[str, str]], target: float) -> tuple[float, float]:
    ordered = sorted(rows, key=lambda row: float(row["snr_db"]))
    for left, right in zip(ordered[:-1], ordered[1:]):
        if float(left["bler"]) >= target >= float(right["bler"]):
            lo = float(left["snr_db"])
            hi = float(right["snr_db"])
            if abs((hi - lo) - 1.0) > 1e-12:
                raise RuntimeError(f"Expected a 1 dB prescan bracket, got {lo:g}..{hi:g}")
            return lo, hi
    raise RuntimeError(f"No raw bracket for target {target:g}")


def _grid(rows: list[dict[str, str]]) -> tuple[list[float], dict[str, list[float]]]:
    brackets = {}
    points = set()
    for target in TARGETS:
        lo, hi = _first_raw_bracket(rows, target)
        brackets[f"{target:g}"] = [lo, hi]
        points.update((lo, lo + 0.5, hi))
    return sorted(points), brackets


def main() -> None:
    screening = _read_csv(SCREENING_CSV)
    retained = {
        str(row["candidate_id"])
        for row in screening
        if str(row["screened_before_formal"]).lower() != "true"
    }
    if len(retained) != 43:
        raise RuntimeError(f"Expected 43 retained candidates, found {len(retained)}")

    receipt_rows = []
    config_paths = []
    total_points = 0
    for scene, source_path in SOURCE_CONFIGS.items():
        config = yaml.safe_load(source_path.read_text(encoding="utf-8")) or {}
        candidates = []
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            if candidate_id not in retained:
                continue
            csv_path = PRESCAN_ROOT / scene / candidate_id / "bler_points.csv"
            rows = _read_csv(csv_path)
            grid, brackets = _grid(rows)
            frozen = dict(candidate)
            frozen["receiver_covariance_mode"] = "matched_effective"
            frozen["snr_points_db"] = grid
            candidates.append(frozen)
            total_points += len(grid)
            receipt_rows.append(
                {
                    "scene": scene,
                    "candidate_id": candidate_id,
                    "prescan_brackets_db": brackets,
                    "coarse_snr_points_db": grid,
                    "prescan_csv": csv_path.relative_to(ROOT).as_posix(),
                    "prescan_csv_sha256": _sha256(csv_path),
                }
            )
        config["output_dir"] = (
            f"outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search/"
            f"coarse_confirmation/{scene}"
        )
        config["simulation"] = {
            "snr_points_db": [0.0],
            "min_trials_per_snr": TRIALS_PER_POINT,
            "target_errors": 1_000_000,
            "max_trials_per_snr": TRIALS_PER_POINT,
            "target_bler": 0.01,
            "save_trial_error_flags": True,
            "resume": True,
            "progress_every_batches": 1,
        }
        config["candidates"] = candidates
        destination = ROOT / "configs" / f"pdcch_result031_strict_sidon_{scene}_coarse.yaml"
        destination.write_text(
            yaml.safe_dump(config, allow_unicode=True, sort_keys=False), encoding="utf-8"
        )
        config_paths.append(destination)

    receipt = {
        "schema": "plan031-strict-sidon-coarse-config-v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "selection_input": SCREENING_CSV.relative_to(ROOT).as_posix(),
        "selection_input_sha256": _sha256(SCREENING_CSV),
        "candidate_count": len(receipt_rows),
        "point_count": total_points,
        "candidate_trials": total_points * TRIALS_PER_POINT,
        "trials_per_point": TRIALS_PER_POINT,
        "grid_rule": "for each retained candidate and each 10%/1% target, use the first descending 1 dB raw prescan bracket endpoints and its 0.5 dB midpoint; take the union",
        "performance_selection_from_prescan": False,
        "configs": [
            {
                "path": path.relative_to(ROOT).as_posix(),
                "sha256": _sha256(path),
            }
            for path in config_paths
        ],
        "candidates": receipt_rows,
        "reproduce_command": "python tools/prepare_plan031_strict_sidon_coarse_configs.py",
    }
    RECEIPT_PATH.write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
