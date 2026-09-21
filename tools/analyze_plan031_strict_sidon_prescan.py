"""Record the frozen pre-formal screening decision for plan-031 strict Sidon sets."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search"
PRESCAN_ROOT = RUN_ROOT / "prescan"
OUTPUT_ROOT = RUN_ROOT / "screening"
TARGET_1PCT = 0.01
TARGET_10PCT = 0.10
CE_FLOOR_GAP_LIMIT_DB = 3.0
HIGH_SNR_TAIL_POINTS = 5


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return sorted(csv.DictReader(handle), key=lambda row: float(row["snr_db"]))


def _raw_bracket(rows: list[dict[str, str]], target: float) -> bool:
    return any(
        float(left["bler"]) >= target >= float(right["bler"])
        for left, right in zip(rows[:-1], rows[1:])
    )


def _summarize(path: Path) -> dict[str, object]:
    rows = _read_rows(path)
    if not rows:
        raise ValueError(f"Empty prescan CSV: {path}")
    candidate_id = str(rows[0]["candidate_id"])
    if any(str(row["candidate_id"]) != candidate_id for row in rows):
        raise ValueError(f"Mixed candidates in {path}")
    minimum = min(rows, key=lambda row: float(row["bler"]))
    last = rows[-1]
    ce_floor_db = float(last["ce_floor_nmse_db"])
    max_snr_ce_db = float(last["ce_nmse_db"])
    tail_median = float(
        np.median([float(row["bler"]) for row in rows[-HIGH_SNR_TAIL_POINTS:]])
    )
    bracket_1pct = _raw_bracket(rows, TARGET_1PCT)
    bracket_10pct = _raw_bracket(rows, TARGET_10PCT)
    minimum_wilson_low = float(minimum["bler_wilson95_lo"])
    ce_gap_db = max_snr_ce_db - ce_floor_db
    tail_above_minimum = tail_median > float(minimum["bler"])
    screened = bool(
        not bracket_1pct
        and minimum_wilson_low > TARGET_1PCT
        and ce_gap_db <= CE_FLOOR_GAP_LIMIT_DB
        and tail_above_minimum
    )
    return {
        "scene": path.parents[1].name,
        "candidate_id": candidate_id,
        "prescan_points": len(rows),
        "prescan_trials": sum(int(row["trials"]) for row in rows),
        "snr_min_db": float(rows[0]["snr_db"]),
        "snr_max_db": float(last["snr_db"]),
        "raw_10pct_bracket": bracket_10pct,
        "raw_1pct_bracket": bracket_1pct,
        "minimum_bler": float(minimum["bler"]),
        "minimum_bler_snr_db": float(minimum["snr_db"]),
        "minimum_point_wilson95_lo": minimum_wilson_low,
        "minimum_point_wilson95_hi": float(minimum["bler_wilson95_hi"]),
        "high_snr_tail_points": HIGH_SNR_TAIL_POINTS,
        "high_snr_tail_median_bler": tail_median,
        "tail_median_above_observed_minimum": tail_above_minimum,
        "zero_noise_ce_floor_nmse_db": ce_floor_db,
        "max_snr_ce_nmse_db": max_snr_ce_db,
        "max_snr_ce_gap_to_floor_db": ce_gap_db,
        "screened_before_formal": screened,
        "decision": (
            "exclude_from_fine-SNR formal run: 1% is unbracketed, the minimum-point "
            "Wilson 95% lower bound exceeds 1%, high-SNR CE is within 3 dB of the "
            "analytic zero-noise floor, and tail BLER has flattened or increased"
            if screened
            else "retain for fine-SNR formal run"
        ),
        "source_csv": path.relative_to(ROOT).as_posix(),
        "source_csv_sha256": _sha256(path),
    }


def main() -> None:
    paths = sorted(PRESCAN_ROOT.glob("*/*/bler_points.csv"))
    rows = [_summarize(path) for path in paths]
    if len(rows) != 48:
        raise RuntimeError(f"Expected 48 candidate CSVs, found {len(rows)}")
    screened = [row for row in rows if bool(row["screened_before_formal"])]
    expected_screened = {
        "C300_SIDON4_AL1_05",
        "C300_SIDON4_AL1_06",
        "C300_SIDON4_AL1_08",
        "C300_SIDON4_AL2_05",
        "C300_SIDON4_AL4_05",
    }
    actual_screened = {str(row["candidate_id"]) for row in screened}
    if actual_screened != expected_screened:
        raise RuntimeError(
            f"Screening result changed: expected {sorted(expected_screened)}, "
            f"got {sorted(actual_screened)}"
        )

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT_ROOT / "candidate_screening.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    receipt = {
        "schema": "plan031-strict-sidon-prescan-screening-v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "researcher_decision_date": "2026-09-10",
        "purpose": "avoid fine-SNR formal runs for candidates with a demonstrated CE/BLER floor; retain the search for the best strict Sidon candidate",
        "criteria": {
            "all_required": True,
            "raw_1pct_bracket": False,
            "minimum_point_wilson95_lower_bound_strictly_above": TARGET_1PCT,
            "max_snr_ce_gap_to_analytic_zero_noise_floor_db_at_most": CE_FLOOR_GAP_LIMIT_DB,
            "high_snr_tail_median_bler_above_observed_minimum": True,
            "high_snr_tail_points": HIGH_SNR_TAIL_POINTS,
        },
        "candidate_count": len(rows),
        "retained_count": len(rows) - len(screened),
        "screened_count": len(screened),
        "screened_candidate_ids": sorted(actual_screened),
        "candidate_csv": csv_path.relative_to(ROOT).as_posix(),
        "candidate_csv_sha256": _sha256(csv_path),
        "reproduce_command": "python tools/analyze_plan031_strict_sidon_prescan.py",
        "interpretation_limit": "CE floor is supporting evidence, not by itself a proof of a BLER floor; the exclusion requires all recorded BLER and CE conditions.",
    }
    (OUTPUT_ROOT / "screening_receipt.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
