"""Freeze the minimal boundary/bracket supplement for plan-031 strict Sidon."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search"
TARGET_INPUT = RUN_ROOT / "fine_analysis/target_snr.csv"
RECEIPT = RUN_ROOT / "fine_supplement_config_receipt.json"
POINTS = {
    "a100_al4": {"A100_SIDON8_AL4_02": [1.25, 1.5]},
    "a100_al8": {
        "A100_SIDON8_AL8_01": [-4.25],
        "A100_SIDON8_AL8_02": [-2.25, -2.0],
        "A100_SIDON8_AL8_05": [-2.25, -2.0],
        "A100_SIDON8_AL8_06": [-2.25, -2.0],
    },
    "c300_al2": {
        "C300_SIDON4_AL2_03": [2.5, 2.75, 6.25, 6.5],
        "C300_SIDON4_AL2_06": [2.5, 2.75, 6.25, 6.5],
    },
    "c300_al4": {
        "C300_SIDON4_AL4_01": [1.25, 1.5],
        "C300_SIDON4_AL4_02": [-1.5, -1.25],
        "C300_SIDON4_AL4_03": [1.25, 1.5],
    },
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if not TARGET_INPUT.exists():
        raise FileNotFoundError("Run the fine analysis audit before freezing the supplement")
    configs = []
    rows = []
    for scene, candidate_points in POINTS.items():
        source = ROOT / "configs" / f"pdcch_result031_strict_sidon_{scene}_fine.yaml"
        config = yaml.safe_load(source.read_text(encoding="utf-8")) or {}
        source_candidates = {
            str(candidate["candidate_id"]): candidate for candidate in config["candidates"]
        }
        candidates = []
        for candidate_id, snr_points in candidate_points.items():
            candidate = dict(source_candidates[candidate_id])
            overlap = set(float(value) for value in candidate["snr_points_db"]) & set(snr_points)
            if overlap:
                raise RuntimeError(f"Supplement repeats frozen fine points: {candidate_id} {overlap}")
            candidate["snr_points_db"] = snr_points
            candidates.append(candidate)
            rows.append(
                {"scene": scene, "candidate_id": candidate_id, "snr_points_db": snr_points}
            )
        config["output_dir"] = (
            "outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search/"
            f"fine_supplement/{scene}"
        )
        config["candidates"] = candidates
        destination = (
            ROOT / "configs" / f"pdcch_result031_strict_sidon_{scene}_fine_supplement.yaml"
        )
        destination.write_text(
            yaml.safe_dump(config, allow_unicode=True, sort_keys=False), encoding="utf-8"
        )
        configs.append(destination)

    receipt = {
        "schema": "plan031-strict-sidon-fine-supplement-v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "trigger": (
            "8 boundary-limited bootstrap targets, 3 unbracketed targets, and one "
            "1.5 dB raw bracket in the completed 149-point fine scan"
        ),
        "target_input": TARGET_INPUT.relative_to(ROOT).as_posix(),
        "target_input_sha256": _sha256(TARGET_INPUT),
        "candidate_count": len(rows),
        "point_count": sum(len(row["snr_points_db"]) for row in rows),
        "min_trials_per_snr": 10_000,
        "target_errors": 200,
        "max_trials_per_snr": 50_000,
        "minimum_candidate_trials": 230_000,
        "maximum_candidate_trials": 1_150_000,
        "configs": [
            {"path": path.relative_to(ROOT).as_posix(), "sha256": _sha256(path)}
            for path in configs
        ],
        "candidates": rows,
        "reproduce_command": "python tools/prepare_plan031_strict_sidon_fine_supplement.py",
    }
    RECEIPT.write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
