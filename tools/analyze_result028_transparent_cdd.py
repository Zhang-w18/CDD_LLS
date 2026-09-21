"""Validate and summarize the formal plan-028 transparent CDD receiver run."""

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


RUN = ROOT / "outputs/experiment028_csi_curves/20260803_main/transparent_cdd_physical_covariance/a100"
FINAL = RUN / "final"
EXPECTED_SNRS = (14.0, 14.25, 14.5, 14.75, 15.0, 15.5, 15.75, 16.0)
CANDIDATES = (
    "A100_B0_QC_TRANSPARENT_CDD",
    "A100_AP_RMS_T1_TRANSPARENT_CDD",
    "A100_AP_TU_NT_TRANSPARENT_CDD",
    "A100_S0_SIDON_TRANSPARENT_CDD",
    "A100_MEFF_T2_06_TRANSPARENT_CDD",
)


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
                raise RuntimeError(
                    f"Non-contiguous intervals for {key}: expected {expected}, got {start}-{end}"
                )
            if end - start + 1 > 1000:
                raise RuntimeError(f"Resumable interval exceeds 1000 trials for {key}.")
            expected = end + 1
    return grouped


def validate_ce_arrays(rows: Sequence[dict]) -> dict:
    checked = 0
    for row in rows:
        nmse = np.load(resolve_repo(str(row["ce_nmse_trial_path"]))).astype(np.float64)
        error = np.load(resolve_repo(str(row["ce_error_energy_path"]))).astype(np.float64)
        signal = np.load(resolve_repo(str(row["ce_true_energy_path"]))).astype(np.float64)
        if len(nmse) != int(row["trials"]) or np.any(~np.isfinite(nmse)):
            raise RuntimeError("CE NMSE trial array is missing or non-finite.")
        np.testing.assert_allclose(
            nmse, error / np.maximum(signal, 1e-30), rtol=1e-13, atol=1e-15
        )
        np.testing.assert_allclose(float(row["ce_nmse_sum"]), np.sum(nmse), rtol=1e-12)
        np.testing.assert_allclose(float(row["ce_nmse_sumsq"]), np.sum(nmse**2), rtol=1e-12)
        checked += len(nmse)
    return {"checked_candidate_trials": checked, "trial_ratio_and_sums_exact": True}


def main() -> None:
    resolved = json.loads((RUN / "resolved_run.json").read_text(encoding="utf-8"))
    manifest_path = resolve_repo(str(resolved["transparent_manifest_path"]))
    manifest_digest = sha256(manifest_path)
    if manifest_digest != str(resolved["transparent_manifest_sha256"]):
        raise RuntimeError("Transparent CDD manifest hash differs from resolved_run.json.")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    receiver = manifest["receiver"]
    if receiver != {
        "knows_cdd_delays": False,
        "assumed_covariance": "8 * R_phy(static TDL-A 100 ns)",
        "assumed_covariance_power_factor": 8.0,
        "uses_candidate_specific_covariance": False,
        "averaged_dmrs_count": 2,
        "averaged_ls_noise_variance": "4/SNR_linear",
        "receiver_modes": ["estimated"],
    }:
        raise RuntimeError("Transparent CDD receiver definition changed.")
    if (RUN / "ideal").exists() or (FINAL / "ideal_csi_bler_points.csv").exists():
        raise RuntimeError("Transparent CDD run unexpectedly contains ideal-CSI data.")

    final_rows = read_rows(FINAL / "estimated_csi_bler_points.csv")
    supplemental = read_rows(RUN / "estimated" / "supplemental_points.csv")
    intervals = validate_intervals(supplemental)
    expected_keys = {(candidate, snr) for candidate in CANDIDATES for snr in EXPECTED_SNRS}
    by_key = {
        (str(row["candidate_id"]), float(row["snr_db"])): row for row in final_rows
    }
    if set(by_key) != expected_keys:
        raise RuntimeError("Transparent CDD final grid differs from the frozen grid.")
    trials_by_snr = {}
    for snr in EXPECTED_SNRS:
        trial_counts = {int(by_key[(candidate, snr)]["trials"]) for candidate in CANDIDATES}
        if len(trial_counts) != 1:
            raise RuntimeError(f"Paired trial counts differ at {snr:g} dB.")
        interval_lists = {tuple(intervals[(candidate, snr)]) for candidate in CANDIDATES}
        if len(interval_lists) != 1:
            raise RuntimeError(f"Paired absolute-trial intervals differ at {snr:g} dB.")
        count = trial_counts.pop()
        trials_by_snr[f"{snr:g}"] = count
        if count < 10000:
            raise RuntimeError(f"Minimum formal budget not reached at {snr:g} dB.")
        errors = [int(by_key[(candidate, snr)]["tb_errors"]) for candidate in CANDIDATES]
        if any(value < 200 for value in errors) and count != 50000:
            raise RuntimeError(
                f"Stopping rule not satisfied at {snr:g} dB: trials={count}, errors={errors}"
            )

    ce_audit = validate_ce_arrays(supplemental)
    stability_rows = []
    reversals = []
    ce_rows = []
    analytic_comparison = []
    for candidate in CANDIDATES:
        previous = None
        for snr in EXPECTED_SNRS:
            row = by_key[(candidate, snr)]
            bler = float(row["bler"])
            reversal = previous is not None and bler > float(previous["bler"])
            item = {
                "receiver": "estimated",
                "candidate_id": candidate,
                "snr_db": snr,
                "trials": int(row["trials"]),
                "tb_errors": int(row["tb_errors"]),
                "bler": bler,
                "bler_wilson95_lo": float(row["bler_wilson95_lo"]),
                "bler_wilson95_hi": float(row["bler_wilson95_hi"]),
                "wilson95_width": float(row["bler_wilson95_hi"])
                - float(row["bler_wilson95_lo"]),
                "errors_below_200": int(row["tb_errors"]) < 200,
                "reversal_from_previous_snr": reversal,
            }
            stability_rows.append(item)
            if reversal:
                reversals.append(item)
            previous = row
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
            diagnostic_name = (
                f"snr_{snr:+.2f}".replace("+", "p").replace("-", "m").replace(".", "p")
                + ".json"
            )
            diagnostic = json.loads(
                (RUN / "filter_diagnostics" / diagnostic_name).read_text(encoding="utf-8")
            )
            analytic = float(diagnostic["candidate_analytic_trace_nmse"][candidate])
            analytic_comparison.append(
                {
                    "candidate_id": candidate,
                    "snr_db": snr,
                    "monte_carlo_trial_ratio_mean": mean,
                    "analytic_trace_ratio": analytic,
                    "monte_carlo_db": 10.0 * math.log10(max(mean, 1e-300)),
                    "analytic_db": 10.0 * math.log10(max(analytic, 1e-300)),
                    "difference_db": 10.0 * math.log10(max(mean, 1e-300))
                    - 10.0 * math.log10(max(analytic, 1e-300)),
                }
            )

    write_rows(stability_rows, FINAL / "stability_audit.csv")
    write_rows(ce_rows, FINAL / "ce_nmse_uncertainty.csv")
    write_rows(analytic_comparison, FINAL / "analytic_vs_monte_carlo_ce.csv")
    mapping = manifest["source_candidate_mapping"]
    (FINAL / "source_mapping_audit.json").write_text(
        json.dumps(
            {
                "mapping": mapping,
                "unique_transparent_ids": len({row["candidate_id"] for row in mapping}) == 5,
                "unique_source_ids": len({row["source_candidate_id"] for row in mapping}) == 5,
                "receiver_knows_cdd_delays": False,
                "ideal_data_absent": True,
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    summary = {
        "schema": "result028-transparent-cdd-analysis-v1",
        "source_manifest_sha256": str(resolved["manifest_sha256"]),
        "transparent_manifest_path": str(Path(resolved["transparent_manifest_path"])),
        "transparent_manifest_sha256": manifest_digest,
        "resolved_config_sha256": str(resolved["config_sha256"]),
        "snr_grid_db": list(EXPECTED_SNRS),
        "candidates": list(CANDIDATES),
        "trials_by_snr": trials_by_snr,
        "common_paired_trials": sum(trials_by_snr.values()),
        "estimated_candidate_trials": len(CANDIDATES) * sum(trials_by_snr.values()),
        "paired_intervals_exact": True,
        "ideal_data_absent": True,
        "ce_audit": ce_audit,
        "reversal_count": len(reversals),
        "reversals": reversals,
        "points_below_200_errors": sum(bool(row["errors_below_200"]) for row in stability_rows),
        "maximum_abs_analytic_mc_difference_db": max(
            abs(float(row["difference_db"])) for row in analytic_comparison
        ),
    }
    (FINAL / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
