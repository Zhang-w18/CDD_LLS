"""Audit plan-028 section 16 matched-estimated and ideal small-delay CDD runs."""

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

MATCHED = ROOT / "outputs/experiment028_csi_curves/20260803_main/small_delay_cdd_matched/a100"
IDEAL = ROOT / "outputs/experiment028_csi_curves/20260803_main/small_delay_cdd_ideal/a100"
MATCHED_ID = "A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD"
IDEAL_ID = "A100_SMALL_CDD_QSTEP0P25_IDEAL_CDD"
DELAYS = (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75)
MATCHED_SNRS = tuple(14.0 + 0.5 * index for index in range(13))
IDEAL_SNRS = tuple(14.0 + 0.5 * index for index in range(11))
MATCHED_SOURCES = (
    MATCHED,
    *(MATCHED.parent.parent / f"small_delay_cdd_matched_shard{i}" / "a100" for i in range(1, 4)),
)
IDEAL_SOURCES = (
    IDEAL,
    *(IDEAL.parent.parent / f"small_delay_cdd_ideal_shard{i}" / "a100" for i in range(1, 4)),
)


def supplements(sources: tuple[Path, ...], receiver: str) -> list[dict]:
    rows: list[dict] = []
    for source in sources:
        path = source / receiver / "supplemental_points.csv"
        if path.exists():
            rows.extend(read_rows(path))
    return rows


def bracket(rows: list[dict], target: float) -> dict:
    ordered = sorted(rows, key=lambda row: float(row["snr_db"]))
    left, right = next(
        (left, right)
        for left, right in zip(ordered, ordered[1:])
        if float(left["bler"]) >= target >= float(right["bler"])
    )
    return {
        "target_bler": target,
        "lower_snr_db": float(left["snr_db"]),
        "lower_bler": float(left["bler"]),
        "upper_snr_db": float(right["snr_db"]),
        "upper_bler": float(right["bler"]),
    }


def audit_manifest(run: Path, receiver: str) -> tuple[dict, dict, str]:
    resolved = json.loads((run / "resolved_run.json").read_text(encoding="utf-8"))
    manifest_path = resolve_repo(resolved["transparent_manifest_path"])
    digest = sha256(manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if digest != resolved["transparent_manifest_sha256"]:
        raise RuntimeError("Section-16 manifest hash mismatch.")
    if manifest["schema"] != "plan028-small-delay-cdd-matched-csi-manifest-v1":
        raise RuntimeError("Unexpected section-16 manifest schema.")
    mapping = manifest["source_candidate_mapping"]
    if len(mapping) != 1 or tuple(mapping[0]["delay_grid_coordinates"]) != DELAYS:
        raise RuntimeError("Section-16 transmitter changed.")
    if receiver == "estimated":
        policy = manifest["receiver"]
        if not policy["knows_cdd_delays"] or policy["covariance_mode"] != "matched_effective":
            raise RuntimeError("Estimated receiver is not matched-effective.")
    elif manifest["receiver"] != {
        "receiver_modes": ["ideal"],
        "ideal_uses_true_effective_data_csi": True,
        "estimated_receiver_not_run": True,
    }:
        raise RuntimeError("Ideal receiver manifest changed.")
    return resolved, manifest, digest


def audit_run(
    run: Path,
    sources: tuple[Path, ...],
    receiver: str,
    candidate: str,
    snrs: tuple[float, ...],
) -> dict:
    resolved, _, manifest_digest = audit_manifest(run, receiver)
    final = run / "final"
    rows = read_rows(final / f"{receiver}_csi_bler_points.csv")
    if tuple(float(row["snr_db"]) for row in rows) != snrs:
        raise RuntimeError(f"{receiver} grid changed.")
    raw = supplements(sources, receiver)
    intervals = validate_intervals(raw)
    interval_count = len(raw)
    for row in rows:
        snr = float(row["snr_db"])
        trials = int(row["trials"])
        errors = int(row["tb_errors"])
        if trials < 10000 or (errors < 200 and trials != 50000):
            raise RuntimeError(f"Stopping rule failed for {receiver} at {snr:g} dB.")
        covered = sum(end - start + 1 for start, end in intervals[(candidate, snr)])
        if covered != trials:
            raise RuntimeError(f"Interval coverage failed for {receiver} at {snr:g} dB.")
    stability = []
    previous = None
    for row in rows:
        reversal = previous is not None and float(row["bler"]) > float(previous["bler"])
        stability.append({
            "receiver": receiver,
            "candidate_id": candidate,
            "snr_db": float(row["snr_db"]),
            "trials": int(row["trials"]),
            "tb_errors": int(row["tb_errors"]),
            "bler": float(row["bler"]),
            "bler_wilson95_lo": float(row["bler_wilson95_lo"]),
            "bler_wilson95_hi": float(row["bler_wilson95_hi"]),
            "reversal_from_previous_snr": reversal,
        })
        previous = row
    write_rows(stability, final / "stability_audit.csv")
    brackets = [bracket(rows, 0.1), bracket(rows, 0.01)]
    write_rows(brackets, final / "bler_bracket_audit.csv")
    summary = {
        "source_manifest_sha256": resolved["manifest_sha256"],
        "receiver_manifest_sha256": manifest_digest,
        "resolved_config_sha256": resolved["config_sha256"],
        "receiver": receiver,
        "candidate": candidate,
        "snr_grid_db": list(snrs),
        "total_trials": sum(int(row["trials"]) for row in rows),
        "interval_count": interval_count,
        "reversal_count": sum(bool(row["reversal_from_previous_snr"]) for row in stability),
        "bler_brackets": brackets,
    }
    if receiver == "estimated":
        ce_audit = validate_ce_arrays(raw)
        ce_rows = []
        analytic_rows = []
        for row in rows:
            snr = float(row["snr_db"])
            lo = float(row["ce_nmse_ci95_lo"])
            hi = float(row["ce_nmse_ci95_hi"])
            ce_rows.append({
                "candidate_id": candidate,
                "snr_db": snr,
                "trials": int(row["trials"]),
                "ce_nmse_mean": float(row["ce_nmse_mean"]),
                "ce_nmse_mean_db": float(row["ce_nmse_mean_db"]),
                "ce_nmse_standard_error": float(row["ce_nmse_standard_error"]),
                "ce_nmse_ci95_lo": lo,
                "ce_nmse_ci95_hi": hi,
                "ce_nmse_ci95_lo_db": 10.0 * math.log10(lo),
                "ce_nmse_ci95_hi_db": 10.0 * math.log10(hi),
            })
            name = f"snr_{snr:+.2f}".replace("+", "p").replace("-", "m").replace(".", "p") + ".json"
            diagnostic = json.loads((run / "filter_diagnostics" / name).read_text(encoding="utf-8"))
            analytic = float(diagnostic["candidate_analytic_trace_nmse"][candidate])
            analytic_rows.append({
                "candidate_id": candidate,
                "snr_db": snr,
                "monte_carlo_db": float(row["ce_nmse_mean_db"]),
                "analytic_db": 10.0 * math.log10(analytic),
                "difference_db": float(row["ce_nmse_mean_db"]) - 10.0 * math.log10(analytic),
            })
        write_rows(ce_rows, final / "ce_nmse_uncertainty.csv")
        write_rows(analytic_rows, final / "analytic_vs_monte_carlo_ce.csv")
        summary.update({
            "ce_audit": ce_audit,
            "ce_nmse_db_range": [
                min(row["ce_nmse_mean_db"] for row in ce_rows),
                max(row["ce_nmse_mean_db"] for row in ce_rows),
            ],
            "maximum_abs_analytic_mc_difference_db": max(
                abs(row["difference_db"]) for row in analytic_rows
            ),
        })
    (final / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    output = {
        "schema": "result028-small-delay-cdd-matched-ideal-analysis-v1",
        "matched": audit_run(MATCHED, MATCHED_SOURCES, "estimated", MATCHED_ID, MATCHED_SNRS),
        "ideal": audit_run(IDEAL, IDEAL_SOURCES, "ideal", IDEAL_ID, IDEAL_SNRS),
    }
    print(json.dumps(output, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
