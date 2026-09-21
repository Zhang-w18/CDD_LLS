"""Freeze final 0.25 dB strict-Sidon configs from the coarse analysis."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search"
ANALYSIS_ROOT = RUN_ROOT / "coarse_analysis"
ANALYSIS_RECEIPT = ANALYSIS_ROOT / "analysis_receipt.json"
TARGET_CSV = ANALYSIS_ROOT / "target_snr_coarse.csv"
RECEIPT_PATH = RUN_ROOT / "fine_config_receipt.json"
MIN_TRIALS = 10_000
TARGET_ERRORS = 200
MAX_TRIALS = 50_000
SOURCE_CONFIGS = {
    scene: ROOT / "configs" / f"pdcch_result031_strict_sidon_{scene}_coarse.yaml"
    for scene in (
        "a100_al2",
        "a100_al4",
        "a100_al8",
        "c300_al1",
        "c300_al2",
        "c300_al4",
    )
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _quarter_grid(row: dict[str, str]) -> tuple[list[float], str]:
    status = row["status"]
    if status == "estimated":
        low = float(row["raw_bracket_lo_db"])
        high = float(row["raw_bracket_hi_db"])
        steps = round((high - low) / 0.25)
        if steps < 1 or abs(low + 0.25 * steps - high) > 1e-12:
            raise RuntimeError(f"Coarse bracket is not on the 0.25 dB grid: {low:g}..{high:g}")
        return (
            [low + 0.25 * index for index in range(steps + 1)],
            "fill_raw_bracket_at_0p25dB",
        )
    bound = row["bound"].split()
    if len(bound) != 3 or bound[2] != "dB":
        raise RuntimeError(f"Invalid censored target bound: {row['bound']!r}")
    edge = float(bound[1])
    if status == "censored_below_target_at_min_snr" and bound[0] == "<":
        return [edge - 0.5, edge - 0.25, edge], "extend_0p5dB_below_grid"
    if status == "censored_above_target_at_max_snr" and bound[0] == ">":
        return [edge + 0.25 * index for index in range(5)], "extend_1dB_above_grid"
    raise RuntimeError(f"Unsupported target status: {status}")


def main() -> None:
    analysis = json.loads(ANALYSIS_RECEIPT.read_text(encoding="utf-8"))
    selected_by_scene = {
        str(scene): [str(value) for value in values]
        for scene, values in analysis["preliminary_candidates_by_scene"].items()
    }
    if sum(len(values) for values in selected_by_scene.values()) != 18:
        raise RuntimeError("Expected exactly 18 coarse-selected candidates")
    targets = _read_csv(TARGET_CSV)
    target_lookup = {
        (row["scene"], row["candidate_id"], float(row["target_bler"])): row
        for row in targets
    }

    receipt_rows = []
    config_paths = []
    total_points = 0
    for scene, source_path in SOURCE_CONFIGS.items():
        config = yaml.safe_load(source_path.read_text(encoding="utf-8")) or {}
        source_candidates = {
            str(candidate["candidate_id"]): candidate for candidate in config["candidates"]
        }
        candidates = []
        for candidate_id in selected_by_scene[scene]:
            if candidate_id not in source_candidates:
                raise RuntimeError(f"Missing {candidate_id} from {source_path.name}")
            points: set[float] = set()
            target_rules = {}
            for target in (0.10, 0.01):
                row = target_lookup[(scene, candidate_id, target)]
                target_points, rule = _quarter_grid(row)
                points.update(round(value, 8) for value in target_points)
                target_rules[f"{target:g}"] = {
                    "coarse_status": row["status"],
                    "coarse_bound": row["bound"] or None,
                    "rule": rule,
                    "points_db": target_points,
                }
            frozen = dict(source_candidates[candidate_id])
            frozen["snr_points_db"] = sorted(points)
            candidates.append(frozen)
            total_points += len(points)
            receipt_rows.append(
                {
                    "scene": scene,
                    "candidate_id": candidate_id,
                    "fine_snr_points_db": sorted(points),
                    "target_rules": target_rules,
                }
            )

        config["random_stream_namespace"] = "plan031_strict_sidon_final_fine_v1"
        config["output_dir"] = (
            "outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search/"
            f"fine/{scene}"
        )
        config["simulation"] = {
            "snr_points_db": [0.0],
            "min_trials_per_snr": MIN_TRIALS,
            "target_errors": TARGET_ERRORS,
            "max_trials_per_snr": MAX_TRIALS,
            "target_bler": 0.01,
            "save_trial_error_flags": True,
            "resume": True,
            "progress_every_batches": 10,
        }
        config["candidates"] = candidates
        destination = ROOT / "configs" / f"pdcch_result031_strict_sidon_{scene}_fine.yaml"
        destination.write_text(
            yaml.safe_dump(config, allow_unicode=True, sort_keys=False), encoding="utf-8"
        )
        config_paths.append(destination)

    receipt = {
        "schema": "plan031-strict-sidon-fine-config-v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "selection_input": ANALYSIS_RECEIPT.relative_to(ROOT).as_posix(),
        "selection_input_sha256": _sha256(ANALYSIS_RECEIPT),
        "target_input": TARGET_CSV.relative_to(ROOT).as_posix(),
        "target_input_sha256": _sha256(TARGET_CSV),
        "candidate_count": len(receipt_rows),
        "point_count": total_points,
        "minimum_candidate_trials": total_points * MIN_TRIALS,
        "maximum_candidate_trials": total_points * MAX_TRIALS,
        "min_trials_per_snr": MIN_TRIALS,
        "target_errors": TARGET_ERRORS,
        "max_trials_per_snr": MAX_TRIALS,
        "grid_rule": (
            "fill each observed raw coarse bracket at 0.25 dB spacing; "
            "for censored targets extend 0.5 dB below or 1.0 dB above the coarse edge"
        ),
        "random_stream_namespace": "plan031_strict_sidon_final_fine_v1",
        "configs": [
            {"path": path.relative_to(ROOT).as_posix(), "sha256": _sha256(path)}
            for path in config_paths
        ],
        "candidates": receipt_rows,
        "reproduce_command": "python tools/prepare_plan031_strict_sidon_fine_configs.py",
    }
    RECEIPT_PATH.write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
