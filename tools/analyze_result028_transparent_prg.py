"""Validate and summarize the formal plan-028 transparent PRG baseline run."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools import run_bler_curves as runner


RUN = ROOT / "outputs/experiment028_csi_curves/20260803_main/transparent_prg_baselines/a100"
FINAL = RUN / "final"
EXPECTED_SNRS = (14.0, 14.25, 14.5, 14.75, 15.0, 15.5, 15.75, 16.0)
CANDIDATES = ("A100_PRG_DFT8_4RB", "A100_PRG_DFT8_6RB")


def read_rows(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def write_rows(rows: Sequence[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted(set().union(*(row.keys() for row in rows)))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resolve_repo(path: str) -> Path:
    value = Path(path)
    return value if value.is_absolute() else ROOT / value


def validate_intervals(rows: Sequence[dict]) -> dict[tuple[str, float], list[tuple[int, int]]]:
    grouped: dict[tuple[str, float], list[tuple[int, int]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["candidate_id"]), float(row["snr_db"]))].append(
            (int(row["trial_start"]), int(row["trial_end"]))
        )
    for key, intervals in grouped.items():
        expected = 1
        for start, end in sorted(intervals):
            if start != expected or end < start:
                raise RuntimeError(f"Non-contiguous intervals for {key}: expected {expected}, got {start}-{end}")
            expected = end + 1
    return grouped


def validate_mapping(rows: Sequence[dict]) -> dict:
    candidate = {
        "candidate_id": "A100_PRG_DFT8_4RB",
        "prg_size_rb": 4,
        "mapping": "cycle8_random_tail4",
    }
    unique_files = {}
    unique_tuples = set()
    checked_trials = 0
    for row in rows:
        if str(row["candidate_id"]) != candidate["candidate_id"]:
            continue
        key = (float(row["snr_db"]), int(row["trial_start"]), int(row["trial_end"]))
        if key in unique_files:
            continue
        path = resolve_repo(str(row["prg_tail_order_path"]))
        digest = sha256(path)
        if digest != str(row["prg_tail_order_sha256"]):
            raise RuntimeError(f"PRG tail mapping hash mismatch: {path}")
        actual = np.load(path).astype(np.int64)
        expected = np.asarray(
            [
                runner.transparent_prg_order(candidate, 20260727, key[0], trial)[8:12]
                for trial in range(key[1], key[2] + 1)
            ],
            dtype=np.int64,
        )
        np.testing.assert_array_equal(actual, expected)
        if any(len(set(row_values.tolist())) != 4 for row_values in actual):
            raise RuntimeError("A 4-RB tail mapping contains a repeated DFT vector.")
        unique_tuples.update(tuple(values.tolist()) for values in actual)
        checked_trials += len(actual)
        unique_files[key] = {"path": str(path.relative_to(ROOT)), "sha256": digest}
    return {
        "checked_intervals": len(unique_files),
        "checked_trials": checked_trials,
        "unique_ordered_tail_tuples": len(unique_tuples),
        "all_without_replacement": True,
        "all_replay_exact": True,
    }


def validate_ce_arrays(rows: Sequence[dict]) -> dict:
    checked = 0
    for row in rows:
        nmse = np.load(resolve_repo(str(row["ce_nmse_trial_path"]))).astype(np.float64)
        error = np.load(resolve_repo(str(row["ce_error_energy_path"]))).astype(np.float64)
        signal = np.load(resolve_repo(str(row["ce_true_energy_path"]))).astype(np.float64)
        if len(nmse) != int(row["trials"]) or np.any(~np.isfinite(nmse)):
            raise RuntimeError("CE NMSE trial array is missing or non-finite.")
        np.testing.assert_allclose(nmse, error / np.maximum(signal, 1e-30), rtol=1e-13, atol=1e-15)
        np.testing.assert_allclose(float(row["ce_nmse_sum"]), np.sum(nmse), rtol=1e-12)
        np.testing.assert_allclose(float(row["ce_nmse_sumsq"]), np.sum(nmse**2), rtol=1e-12)
        checked += len(nmse)
    return {"checked_candidate_trials": checked, "trial_ratio_and_sums_exact": True}


def main() -> None:
    resolved = json.loads((RUN / "resolved_run.json").read_text(encoding="utf-8"))
    transparent_manifest_path = resolve_repo(str(resolved["transparent_manifest_path"]))
    transparent_manifest_sha256 = sha256(transparent_manifest_path)
    if transparent_manifest_sha256 != str(resolved["transparent_manifest_sha256"]):
        raise RuntimeError("Transparent PRG manifest hash differs from resolved_run.json.")
    final_rows = {receiver: read_rows(FINAL / f"{receiver}_csi_bler_points.csv") for receiver in runner.RECEIVERS}
    supplemental = {
        receiver: read_rows(RUN / receiver / "supplemental_points.csv") for receiver in runner.RECEIVERS
    }
    intervals = {receiver: validate_intervals(rows) for receiver, rows in supplemental.items()}
    if intervals["estimated"] != intervals["ideal"]:
        raise RuntimeError("Estimated and ideal absolute trial intervals differ.")
    expected_keys = {(candidate, snr) for candidate in CANDIDATES for snr in EXPECTED_SNRS}
    for receiver, rows in final_rows.items():
        keys = {(str(row["candidate_id"]), float(row["snr_db"])) for row in rows}
        if keys != expected_keys:
            raise RuntimeError(f"{receiver} final grid differs from the frozen grid.")
    by_receiver_key = {
        receiver: {(str(row["candidate_id"]), float(row["snr_db"])): row for row in rows}
        for receiver, rows in final_rows.items()
    }
    trials_by_snr = {}
    for snr in EXPECTED_SNRS:
        trials = {
            int(by_receiver_key[receiver][(candidate, snr)]["trials"])
            for receiver in runner.RECEIVERS
            for candidate in CANDIDATES
        }
        if len(trials) != 1:
            raise RuntimeError(f"Paired trial counts differ at {snr:g} dB.")
        count = trials.pop()
        trials_by_snr[f"{snr:g}"] = count
        if count < 10000:
            raise RuntimeError(f"Minimum formal budget not reached at {snr:g} dB.")
        errors = [
            int(by_receiver_key[receiver][(candidate, snr)]["tb_errors"])
            for receiver in runner.RECEIVERS
            for candidate in CANDIDATES
        ]
        if any(value < 200 for value in errors) and count != 50000:
            raise RuntimeError(f"Stopping rule not satisfied at {snr:g} dB: trials={count}, errors={errors}")

    mapping_audit = validate_mapping(supplemental["estimated"])
    ce_audit = validate_ce_arrays(supplemental["estimated"])
    stability_rows = []
    reversals = []
    for receiver in runner.RECEIVERS:
        for candidate in CANDIDATES:
            values = [by_receiver_key[receiver][(candidate, snr)] for snr in EXPECTED_SNRS]
            previous = None
            for row in values:
                bler = float(row["bler"])
                reversal = previous is not None and bler > float(previous["bler"])
                overlap = ""
                if previous is not None:
                    overlap = not (
                        float(row["bler_wilson95_lo"]) > float(previous["bler_wilson95_hi"])
                        or float(previous["bler_wilson95_lo"]) > float(row["bler_wilson95_hi"])
                    )
                item = {
                    "receiver": receiver,
                    "candidate_id": candidate,
                    "snr_db": float(row["snr_db"]),
                    "trials": int(row["trials"]),
                    "tb_errors": int(row["tb_errors"]),
                    "bler": bler,
                    "bler_wilson95_lo": float(row["bler_wilson95_lo"]),
                    "bler_wilson95_hi": float(row["bler_wilson95_hi"]),
                    "wilson95_width": float(row["bler_wilson95_hi"]) - float(row["bler_wilson95_lo"]),
                    "errors_below_200": int(row["tb_errors"]) < 200,
                    "errors_below_30": int(row["tb_errors"]) < 30,
                    "reversal_from_previous_snr": reversal,
                    "wilson_interval_overlaps_previous": overlap,
                }
                stability_rows.append(item)
                if reversal:
                    reversals.append(item)
                previous = row
    ce_rows = []
    for candidate in CANDIDATES:
        for snr in EXPECTED_SNRS:
            row = by_receiver_key["estimated"][(candidate, snr)]
            mean = float(row["ce_nmse_mean"])
            lo = float(row["ce_nmse_ci95_lo"])
            hi = float(row["ce_nmse_ci95_hi"])
            ce_rows.append(
                {
                    "candidate_id": candidate,
                    "snr_db": snr,
                    "trials": int(row["trials"]),
                    "ce_nmse_mean": mean,
                    "ce_nmse_mean_db": float(row["ce_nmse_mean_db"]),
                    "ce_nmse_standard_error": float(row["ce_nmse_standard_error"]),
                    "ce_nmse_ci95_lo": lo,
                    "ce_nmse_ci95_hi": hi,
                    "ce_nmse_ci95_lo_db": 10.0 * math.log10(max(lo, 1e-300)),
                    "ce_nmse_ci95_hi_db": 10.0 * math.log10(max(hi, 1e-300)),
                }
            )
    write_rows(stability_rows, FINAL / "stability_audit.csv")
    write_rows(ce_rows, FINAL / "ce_nmse_uncertainty.csv")
    (FINAL / "prg_mapping_audit.json").write_text(
        json.dumps(mapping_audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    summary = {
        "schema": "result028-transparent-prg-analysis-v1",
        "source_manifest_sha256": str(resolved["manifest_sha256"]),
        "transparent_manifest_path": str(Path(resolved["transparent_manifest_path"])),
        "transparent_manifest_sha256": transparent_manifest_sha256,
        "resolved_config_sha256": str(resolved["config_sha256"]),
        "snr_grid_db": list(EXPECTED_SNRS),
        "candidates": list(CANDIDATES),
        "trials_by_snr": trials_by_snr,
        "common_paired_trials": sum(trials_by_snr.values()),
        "candidate_trials_per_receiver": len(CANDIDATES) * sum(trials_by_snr.values()),
        "paired_intervals_exact": True,
        "mapping_audit": mapping_audit,
        "ce_audit": ce_audit,
        "reversal_count": len(reversals),
        "reversals": reversals,
        "points_below_200_errors": sum(bool(row["errors_below_200"]) for row in stability_rows),
        "points_below_30_errors": sum(bool(row["errors_below_30"]) for row in stability_rows),
    }
    (FINAL / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
