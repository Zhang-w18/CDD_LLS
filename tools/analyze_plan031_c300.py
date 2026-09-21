"""Validate and analyze the plan-031 TDL-C/4Tx/2-symbol formal runs."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import platform
import subprocess

import numpy as np

import analyze_plan031 as common


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs" / "experiment031_pdcch_cdd" / "20260907_c300_4tx_2sym"
ANALYSIS = OUTPUT / "analysis"
FIGURES = ROOT / "docs" / "figures" / "result-031" / "c300_4tx_2sym"
CONFIGS = {
    al: ROOT / "configs" / f"pdcch_result031_c300_al{al}_bracket_completion.yaml"
    for al in (1, 2, 4)
}
TARGETS = (0.10, 0.01)
BASELINES = ("C300_B0_QC", "C300_PRG_DFT4_6REG")
MAX_BRACKET_WIDTH_DB = 0.2500001


def qualified_raw_bracket(rows: list[dict[str, object]], target: float) -> tuple[float, float]:
    brackets = []
    for left, right in zip(rows[:-1], rows[1:]):
        lo = float(left["snr_db"])
        hi = float(right["snr_db"])
        if float(left["bler"]) >= target >= float(right["bler"]) and hi - lo <= MAX_BRACKET_WIDTH_DB:
            brackets.append((lo, hi))
    if not brackets:
        raise ValueError("No adjacent raw bracket with width <= 0.25 dB.")
    probabilities = np.asarray(
        [(int(row["errors"]) + 0.5) / (int(row["trials"]) + 1.0) for row in rows]
    )
    trials = np.asarray([int(row["trials"]) for row in rows], dtype=np.float64)
    fitted = common.isotonic_decreasing(probabilities, trials)
    estimate, _ = common.interpolate_target(
        np.asarray([float(row["snr_db"]) for row in rows]), fitted, target
    )
    return min(brackets, key=lambda pair: abs(0.5 * (pair[0] + pair[1]) - estimate))


def analyze(configs: dict[int, dict], points: list[dict[str, object]], repeats: int):
    targets: list[dict[str, object]] = []
    replicate_map: dict[tuple[int, str, float], np.ndarray] = {}
    for al, config in configs.items():
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            rows = sorted(
                [
                    row for row in points
                    if int(row["aggregation_level"]) == al
                    and str(row["candidate_id"]) == candidate_id
                ],
                key=lambda row: float(row["snr_db"]),
            )
            for target in TARGETS:
                key = (al, candidate_id, target)
                try:
                    bracket_lo, bracket_hi = qualified_raw_bracket(rows, target)
                    estimate, ci_lo, ci_hi, replicates = common.bootstrap_targets(
                        rows,
                        target,
                        repeats,
                        common.stable_seed(20260907, al, candidate_id, target, "c300-bootstrap"),
                    )
                except ValueError as error:
                    replicate_map[key] = np.asarray([], dtype=np.float64)
                    targets.append(
                        {
                            "aggregation_level": al,
                            "candidate_id": candidate_id,
                            "target_bler": target,
                            "target_snr_db": None,
                            "target_snr_ci95_lo_db": None,
                            "target_snr_ci95_hi_db": None,
                            "raw_bracket_lo_db": None,
                            "raw_bracket_hi_db": None,
                            "raw_bracket_width_db": None,
                            "status": "unqualified",
                            "method": f"not estimated: {error}",
                            "ci_method": "not applicable",
                            "bootstrap_repeats_requested": repeats,
                            "bootstrap_repeats_valid": 0,
                        }
                    )
                    continue
                replicate_map[key] = replicates
                targets.append(
                    {
                        "aggregation_level": al,
                        "candidate_id": candidate_id,
                        "target_bler": target,
                        "target_snr_db": estimate,
                        "target_snr_ci95_lo_db": ci_lo,
                        "target_snr_ci95_hi_db": ci_hi,
                        "raw_bracket_lo_db": bracket_lo,
                        "raw_bracket_hi_db": bracket_hi,
                        "raw_bracket_width_db": bracket_hi - bracket_lo,
                        "status": "estimated",
                        "method": "Jeffreys-smoothed weighted decreasing isotonic curve; local log10-BLER interpolation; raw adjacent bracket <=0.25 dB",
                        "ci_method": "independent per-point Bernoulli bootstrap",
                        "bootstrap_repeats_requested": repeats,
                        "bootstrap_repeats_valid": len(replicates),
                    }
                )
    lookup = {
        (int(row["aggregation_level"]), str(row["candidate_id"]), float(row["target_bler"])): row
        for row in targets
    }
    gains: list[dict[str, object]] = []
    for al, config in configs.items():
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            for baseline in BASELINES:
                for target in TARGETS:
                    candidate_row = lookup[(al, candidate_id, target)]
                    baseline_row = lookup[(al, baseline, target)]
                    candidate_rep = replicate_map[(al, candidate_id, target)]
                    baseline_rep = replicate_map[(al, baseline, target)]
                    count = min(len(candidate_rep), len(baseline_rep))
                    if count == 0:
                        gain, low, high, status = None, None, None, "unavailable_unqualified_target"
                    else:
                        samples = baseline_rep[:count] - candidate_rep[:count]
                        low, high = (float(value) for value in np.quantile(samples, [0.025, 0.975]))
                        gain = float(baseline_row["target_snr_db"]) - float(candidate_row["target_snr_db"])
                        status = "estimated"
                    gains.append(
                        {
                            "aggregation_level": al,
                            "candidate_id": candidate_id,
                            "baseline_id": baseline,
                            "target_bler": target,
                            "gain_db": gain,
                            "gain_ci95_lo_db": low,
                            "gain_ci95_hi_db": high,
                            "status": status,
                            "positive_means_candidate_requires_less_snr": True,
                            "bootstrap_repeats_valid": count,
                        }
                    )
    return targets, gains


def environment_receipt() -> dict[str, object]:
    import matplotlib
    import numpy
    import scipy
    import sionna
    import tensorflow
    import yaml

    return {
        "platform": platform.platform(),
        "python": platform.python_version(),
        "packages": {
            "numpy": numpy.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
            "pyyaml": yaml.__version__,
            "tensorflow": tensorflow.__version__,
            "sionna": sionna.__version__,
        },
    }


def ce_diagnostic_rows(
    points: list[dict[str, object]],
    diagnostics: list[dict[str, object]],
    targets: list[dict[str, object]],
) -> list[dict[str, object]]:
    target_lookup = {
        (int(row["aggregation_level"]), str(row["candidate_id"])): row
        for row in targets
        if math.isclose(float(row["target_bler"]), 0.01)
    }
    rows = []
    for diagnostic in diagnostics:
        al = int(diagnostic["aggregation_level"])
        candidate_id = str(diagnostic["candidate_id"])
        candidate_points = [
            row for row in points
            if int(row["aggregation_level"]) == al and str(row["candidate_id"]) == candidate_id
        ]
        target = target_lookup[(al, candidate_id)]
        if target["status"] == "estimated":
            reference = min(
                candidate_points,
                key=lambda row: abs(float(row["snr_db"]) - float(target["target_snr_db"])),
            )
            selection = "raw SNR point nearest the qualified SNR@1% estimate"
        else:
            reference = min(candidate_points, key=lambda row: abs(float(row["bler"]) - 0.01))
            selection = "raw point nearest 1% BLER; target SNR unqualified"
        rows.append(
            {
                "aggregation_level": al,
                "candidate_id": candidate_id,
                "receiver_covariance_mode": diagnostic["receiver_covariance_mode"],
                "pilot_rank": diagnostic["pilot_rank"],
                "pilot_condition_number": diagnostic["pilot_condition_number"],
                "ce_floor_nmse_db": diagnostic["ce_floor_nmse_db"],
                "reference_snr_db": reference["snr_db"],
                "reference_bler": reference["bler"],
                "reference_ce_nmse_db": reference["ce_nmse_db"],
                "selection": selection,
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bootstrap-repeats", type=int, default=2000)
    args = parser.parse_args()
    common.CONFIGS = CONFIGS
    common.OUTPUT = OUTPUT
    common.ANALYSIS = ANALYSIS
    common.FIGURES = FIGURES
    configs, points, diagnostics = common.collect()
    targets, gains = analyze(configs, points, int(args.bootstrap_repeats))
    ce_rows = ce_diagnostic_rows(points, diagnostics, targets)
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    common.write_csv(ANALYSIS / "formal_points.csv", points)
    common.write_csv(ANALYSIS / "diagnostics.csv", diagnostics)
    common.write_csv(ANALYSIS / "target_snr.csv", targets)
    common.write_csv(ANALYSIS / "target_gains.csv", gains)
    common.write_csv(ANALYSIS / "ce_near_1pct.csv", ce_rows)
    common.plot(configs, points)
    git_head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True
    ).stdout.strip()
    metadata = {
        "schema": "plan031-c300-analysis-v1",
        "plan": "research/plan-031-PDCCH-CDD时延-BLER.md",
        "bootstrap_repeats": int(args.bootstrap_repeats),
        "git_head": git_head,
        "working_tree_note": "plan-031 C300 implementation and result files are uncommitted",
        "config_sha256": {
            path.relative_to(ROOT).as_posix(): common.sha256(path) for path in CONFIGS.values()
        },
        "formal_point_count": len(points),
        "total_trials": sum(int(row["trials"]) for row in points),
        "total_errors": sum(int(row["errors"]) for row in points),
        "target_rows": len(targets),
        "qualified_target_rows": sum(row["status"] == "estimated" for row in targets),
        "gain_rows": len(gains),
    }
    (ANALYSIS / "analysis_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (OUTPUT / "environment_receipt.json").write_text(
        json.dumps(environment_receipt(), indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
