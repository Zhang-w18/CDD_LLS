"""Validate and summarize the plan-028 section 15 PRG-6 estimated-CSI extension."""

from __future__ import annotations

import json
import math
from pathlib import Path

from analyze_result028_transparent_prg import (
    ROOT,
    read_rows,
    resolve_repo,
    sha256,
    validate_ce_arrays,
    validate_intervals,
    write_rows,
)


RUN = ROOT / "outputs/experiment028_csi_curves/20260803_main/transparent_prg_6rb_1pct_extension/a100"
FINAL = RUN / "final"
ORIGINAL_FINAL = ROOT / "outputs/experiment028_csi_curves/20260803_main/transparent_prg_baselines/a100/final"
CANDIDATE = "A100_PRG_DFT8_6RB"
EXPECTED_SNRS = (16.25, 16.5, 16.75, 17.0, 17.25, 17.5, 17.75)
SUPPLEMENT_RUNS = (
    "transparent_prg_6rb_1pct_extension",
    "transparent_prg_6rb_1pct_extension_low",
    "transparent_prg_6rb_1pct_extension_high",
    "transparent_prg_6rb_1pct_extension_shard1",
    "transparent_prg_6rb_1pct_extension_shard2",
    "transparent_prg_6rb_1pct_extension_shard3",
    "transparent_prg_6rb_1pct_extension_shard4",
    "transparent_prg_6rb_1pct_extension_final1",
    "transparent_prg_6rb_1pct_extension_final2",
    "transparent_prg_6rb_1pct_extension_final3",
)


def _one_percent_bracket(rows: list[dict]) -> dict:
    ordered = sorted(rows, key=lambda row: float(row["snr_db"]))
    for lower, upper in zip(ordered, ordered[1:]):
        if float(lower["bler"]) >= 0.01 and float(upper["bler"]) <= 0.01:
            return {
                "target_bler": 0.01,
                "lower_snr_db": float(lower["snr_db"]),
                "lower_bler": float(lower["bler"]),
                "upper_snr_db": float(upper["snr_db"]),
                "upper_bler": float(upper["bler"]),
            }
    raise RuntimeError("PRG-6 estimated-CSI curve does not bracket 1% BLER.")


def main() -> None:
    resolved = json.loads((RUN / "resolved_run.json").read_text(encoding="utf-8"))
    manifest_path = resolve_repo(str(resolved["transparent_manifest_path"]))
    manifest_digest = sha256(manifest_path)
    if manifest_digest != str(resolved["transparent_manifest_sha256"]):
        raise RuntimeError("Transparent PRG manifest hash differs from resolved_run.json.")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    candidates = manifest.get("candidates", [])
    if len(candidates) != 1 or candidates[0].get("candidate_id") != CANDIDATE:
        raise RuntimeError("Section 15 manifest must contain only transparent PRG 6 RB.")
    if int(candidates[0].get("prg_size_rb", 0)) != 6 or candidates[0].get("mapping") != "cycle_all":
        raise RuntimeError("Section 15 PRG-6 transmitter definition changed.")
    if manifest.get("plan") != "research/plan-028-comb6三类CSI-TDL-A300ns.md#15-a100-8tx1rx-第-36-节两条曲线的-snr-延伸":
        raise RuntimeError("Section 15 plan provenance is missing.")
    if (RUN / "ideal").exists() or (FINAL / "ideal_csi_bler_points.csv").exists():
        raise RuntimeError("Section 15 extension unexpectedly contains ideal-CSI data.")

    rows = read_rows(FINAL / "estimated_csi_bler_points.csv")
    supplemental = []
    for run_name in SUPPLEMENT_RUNS:
        path = (
            ROOT
            / "outputs/experiment028_csi_curves/20260803_main"
            / run_name
            / "a100/estimated/supplemental_points.csv"
        )
        if path.exists():
            supplemental.extend(read_rows(path))
    by_snr = {float(row["snr_db"]): row for row in rows}
    if set(by_snr) != set(EXPECTED_SNRS) or any(row["candidate_id"] != CANDIDATE for row in rows):
        raise RuntimeError("Section 15 PRG-6 final grid differs from the frozen grid.")
    intervals = validate_intervals(supplemental)
    stability_rows = []
    ce_rows = []
    previous = None
    reversals = []
    for snr in EXPECTED_SNRS:
        row = by_snr[snr]
        trials = int(row["trials"])
        errors = int(row["tb_errors"])
        if trials < 10000 or (errors < 200 and trials != 50000):
            raise RuntimeError(f"Stopping rule not satisfied at {snr:g} dB.")
        covered = sum(end - start + 1 for start, end in intervals[(CANDIDATE, snr)])
        if covered != trials:
            raise RuntimeError(f"Interval coverage differs at {snr:g} dB.")
        bler = float(row["bler"])
        reversal = previous is not None and bler > float(previous["bler"])
        audit = {
            "candidate_id": CANDIDATE,
            "receiver": "estimated",
            "snr_db": snr,
            "trials": trials,
            "tb_errors": errors,
            "bler": bler,
            "bler_wilson95_lo": float(row["bler_wilson95_lo"]),
            "bler_wilson95_hi": float(row["bler_wilson95_hi"]),
            "reversal_from_previous_snr": reversal,
        }
        stability_rows.append(audit)
        if reversal:
            reversals.append(audit)
        previous = row
        lo = float(row["ce_nmse_ci95_lo"])
        hi = float(row["ce_nmse_ci95_hi"])
        ce_rows.append(
            {
                "candidate_id": CANDIDATE,
                "snr_db": snr,
                "trials": trials,
                "ce_nmse_mean": float(row["ce_nmse_mean"]),
                "ce_nmse_mean_db": float(row["ce_nmse_mean_db"]),
                "ce_nmse_standard_error": float(row["ce_nmse_standard_error"]),
                "ce_nmse_ci95_lo": lo,
                "ce_nmse_ci95_hi": hi,
                "ce_nmse_ci95_lo_db": 10.0 * math.log10(max(lo, 1e-300)),
                "ce_nmse_ci95_hi_db": 10.0 * math.log10(max(hi, 1e-300)),
            }
        )

    original = [
        row
        for row in read_rows(ORIGINAL_FINAL / "estimated_csi_bler_points.csv")
        if row["candidate_id"] == CANDIDATE and math.isclose(float(row["snr_db"]), 16.0)
    ]
    if len(original) != 1:
        raise RuntimeError("Original PRG-6 16 dB anchor is missing.")
    bracket = _one_percent_bracket([original[0], *rows])
    ce_audit = validate_ce_arrays(supplemental)
    write_rows(stability_rows, FINAL / "stability_audit.csv")
    write_rows(ce_rows, FINAL / "ce_nmse_uncertainty.csv")
    write_rows([bracket], FINAL / "bler_bracket_audit.csv")
    summary = {
        "schema": "result028-transparent-prg6-1pct-extension-analysis-v1",
        "candidate": CANDIDATE,
        "snr_grid_db": list(EXPECTED_SNRS),
        "source_manifest_sha256": str(resolved["manifest_sha256"]),
        "transparent_manifest_sha256": manifest_digest,
        "resolved_config_sha256": str(resolved["config_sha256"]),
        "total_trials": sum(int(row["trials"]) for row in rows),
        "ideal_data_absent": True,
        "one_percent_bracket": bracket,
        "ce_nmse_db_range": [
            min(row["ce_nmse_mean_db"] for row in ce_rows),
            max(row["ce_nmse_mean_db"] for row in ce_rows),
        ],
        "ce_audit": ce_audit,
        "reversal_count": len(reversals),
        "reversals": reversals,
    }
    (FINAL / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
