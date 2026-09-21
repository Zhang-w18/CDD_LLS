"""Validate and summarize the formal plan-028 small-delay transparent CDD run."""

from __future__ import annotations

import json
import math
from pathlib import Path

from analyze_result028_transparent_cdd import (
    ROOT,
    read_rows,
    resolve_repo,
    sha256,
    validate_ce_arrays,
    validate_intervals,
    write_rows,
)

RUN = ROOT / "outputs/experiment028_csi_curves/20260803_main/transparent_cdd_small_delay/a100"
FINAL = RUN / "final"
CANDIDATE = "A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD"
EXPECTED_DELAYS = (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75)
EXPECTED_SNRS = (
    14.0, 14.5, 15.0, 15.5, 15.75, 16.0, 16.25, 16.5, 17.0, 17.5, 18.0, 18.5,
    19.0, 19.25, 19.5, 19.75, 20.0, 20.25, 20.5, 21.0,
)


def _bracket(rows: list[dict], target: float) -> dict:
    ordered = sorted(rows, key=lambda row: float(row["snr_db"]))
    pair = next(
        (left, right)
        for left, right in zip(ordered, ordered[1:])
        if float(left["bler"]) >= target and float(right["bler"]) <= target
    ) if any(
        float(left["bler"]) >= target and float(right["bler"]) <= target
        for left, right in zip(ordered, ordered[1:])
    ) else None
    if pair is None:
        raise RuntimeError(f"BLER target {target:g} is not bracketed by adjacent sampled points.")
    lower, upper = pair
    return {
        "target_bler": target,
        "lower_snr_db": float(lower["snr_db"]),
        "lower_bler": float(lower["bler"]),
        "upper_snr_db": float(upper["snr_db"]),
        "upper_bler": float(upper["bler"]),
    }


def main() -> None:
    resolved = json.loads((RUN / "resolved_run.json").read_text(encoding="utf-8"))
    manifest_path = resolve_repo(str(resolved["transparent_manifest_path"]))
    manifest_digest = sha256(manifest_path)
    if manifest_digest != str(resolved["transparent_manifest_sha256"]):
        raise RuntimeError("Transparent CDD manifest hash differs from resolved_run.json.")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest["schema"] != "plan028-transparent-cdd-explicit-delay-manifest-v2":
        raise RuntimeError("Unexpected explicit-delay manifest schema.")
    mapping = manifest["source_candidate_mapping"]
    if len(mapping) != 1 or mapping[0]["candidate_id"] != CANDIDATE:
        raise RuntimeError("Small-delay candidate mapping changed.")
    if tuple(float(value) for value in mapping[0]["delay_grid_coordinates"]) != EXPECTED_DELAYS:
        raise RuntimeError("Frozen small-delay coordinates changed.")
    if mapping[0]["source_candidate_id"] or mapping[0]["transmitter_definition"] != "explicit_delay_grid_coordinates":
        raise RuntimeError("Small-delay transmitter is not the frozen explicit-delay definition.")
    receiver = manifest["receiver"]
    if receiver["knows_cdd_delays"] or receiver["uses_candidate_specific_covariance"]:
        raise RuntimeError("Receiver unexpectedly knows the CDD delays.")
    if receiver["assumed_covariance"] != "8 * R_phy(static TDL-A 100 ns)":
        raise RuntimeError("Receiver covariance assumption changed.")
    if (RUN / "ideal").exists() or (FINAL / "ideal_csi_bler_points.csv").exists():
        raise RuntimeError("Small-delay run unexpectedly contains ideal-CSI data.")

    rows = read_rows(FINAL / "estimated_csi_bler_points.csv")
    supplemental = read_rows(RUN / "estimated" / "supplemental_points.csv")
    by_snr = {float(row["snr_db"]): row for row in rows}
    if set(by_snr) != set(EXPECTED_SNRS) or any(row["candidate_id"] != CANDIDATE for row in rows):
        raise RuntimeError("Small-delay final grid differs from the frozen grid.")
    intervals = validate_intervals(supplemental)
    trials_by_snr = {}
    for snr in EXPECTED_SNRS:
        row = by_snr[snr]
        trials = int(row["trials"])
        errors = int(row["tb_errors"])
        if trials < 10000 or (errors < 200 and trials != 50000):
            raise RuntimeError(f"Stopping rule not satisfied at {snr:g} dB.")
        covered = sum(end - start + 1 for start, end in intervals[(CANDIDATE, snr)])
        if covered != trials:
            raise RuntimeError(f"Interval coverage differs from final trial count at {snr:g} dB.")
        trials_by_snr[f"{snr:g}"] = trials

    ce_audit = validate_ce_arrays(supplemental)
    stability_rows = []
    ce_rows = []
    analytic_rows = []
    reversals = []
    previous = None
    for snr in EXPECTED_SNRS:
        row = by_snr[snr]
        bler = float(row["bler"])
        reversal = previous is not None and bler > float(previous["bler"])
        stability = {
            "receiver": "estimated",
            "candidate_id": CANDIDATE,
            "snr_db": snr,
            "trials": int(row["trials"]),
            "tb_errors": int(row["tb_errors"]),
            "bler": bler,
            "bler_wilson95_lo": float(row["bler_wilson95_lo"]),
            "bler_wilson95_hi": float(row["bler_wilson95_hi"]),
            "reversal_from_previous_snr": reversal,
        }
        stability_rows.append(stability)
        if reversal:
            reversals.append(stability)
        previous = row
        mean = float(row["ce_nmse_mean"])
        mean_db = float(row["ce_nmse_mean_db"])
        if mean_db >= -15.0:
            raise RuntimeError(f"CE NMSE fails the -15 dB requirement at {snr:g} dB.")
        lo = float(row["ce_nmse_ci95_lo"])
        hi = float(row["ce_nmse_ci95_hi"])
        ce_rows.append({
            "candidate_id": CANDIDATE,
            "snr_db": snr,
            "trials": int(row["trials"]),
            "ce_nmse_mean": mean,
            "ce_nmse_mean_db": mean_db,
            "ce_nmse_standard_error": float(row["ce_nmse_standard_error"]),
            "ce_nmse_ci95_lo": lo,
            "ce_nmse_ci95_hi": hi,
            "ce_nmse_ci95_lo_db": 10.0 * math.log10(max(lo, 1e-300)),
            "ce_nmse_ci95_hi_db": 10.0 * math.log10(max(hi, 1e-300)),
        })
        diagnostic_name = f"snr_{snr:+.2f}".replace("+", "p").replace("-", "m").replace(".", "p") + ".json"
        diagnostic = json.loads((RUN / "filter_diagnostics" / diagnostic_name).read_text(encoding="utf-8"))
        analytic = float(diagnostic["candidate_analytic_trace_nmse"][CANDIDATE])
        analytic_db = 10.0 * math.log10(max(analytic, 1e-300))
        analytic_rows.append({
            "candidate_id": CANDIDATE,
            "snr_db": snr,
            "monte_carlo_trial_ratio_mean": mean,
            "analytic_trace_ratio": analytic,
            "monte_carlo_db": mean_db,
            "analytic_db": analytic_db,
            "difference_db": mean_db - analytic_db,
        })

    bracket_rows = [_bracket(rows, 0.1), _bracket(rows, 0.01)]
    write_rows(stability_rows, FINAL / "stability_audit.csv")
    write_rows(ce_rows, FINAL / "ce_nmse_uncertainty.csv")
    write_rows(analytic_rows, FINAL / "analytic_vs_monte_carlo_ce.csv")
    write_rows(bracket_rows, FINAL / "bler_bracket_audit.csv")
    summary = {
        "schema": "result028-transparent-cdd-small-delay-analysis-v1",
        "source_manifest_sha256": str(resolved["manifest_sha256"]),
        "transparent_manifest_sha256": manifest_digest,
        "resolved_config_sha256": str(resolved["config_sha256"]),
        "candidate": CANDIDATE,
        "delay_grid_coordinates": list(EXPECTED_DELAYS),
        "snr_grid_db": list(EXPECTED_SNRS),
        "trials_by_snr": trials_by_snr,
        "total_trials": sum(trials_by_snr.values()),
        "ideal_data_absent": True,
        "ce_all_below_minus_15_db": True,
        "ce_nmse_db_range": [min(row["ce_nmse_mean_db"] for row in ce_rows), max(row["ce_nmse_mean_db"] for row in ce_rows)],
        "bler_brackets": bracket_rows,
        "ce_audit": ce_audit,
        "reversal_count": len(reversals),
        "reversals": reversals,
        "maximum_abs_analytic_mc_difference_db": max(abs(row["difference_db"]) for row in analytic_rows),
    }
    (FINAL / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
