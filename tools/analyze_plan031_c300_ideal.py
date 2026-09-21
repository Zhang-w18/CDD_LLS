"""Analyze the plan-031 C300 two-symbol ideal-CSI formal runs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import platform
import subprocess

import numpy as np

import analyze_plan031 as common
import analyze_plan031_c300 as c300


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs" / "experiment031_pdcch_cdd" / "20260910_c300_4tx_2sym_ideal_csi"
ANALYSIS = OUTPUT / "analysis"
FIGURES = ROOT / "docs" / "figures" / "result-031" / "c300_4tx_2sym_ideal_csi"
CONFIGS = {
    al: ROOT / "configs" / f"pdcch_result031_c300_ideal_al{al}_formal.yaml"
    for al in (1, 2, 4)
}
ESTIMATED_ANALYSIS = (
    ROOT
    / "outputs"
    / "experiment031_pdcch_cdd"
    / "20260907_c300_4tx_2sym"
    / "analysis"
)
TARGETS = (0.10, 0.01)
BASELINES = ("C300_B0_QC", "C300_PRG_DFT4_6REG")


def _read_csv(path: Path) -> list[dict[str, str]]:
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _analyze_targets(
    configs: dict[int, dict], points: list[dict[str, object]], repeats: int
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[tuple[int, str, float], np.ndarray]]:
    rows_out: list[dict[str, object]] = []
    replicate_map: dict[tuple[int, str, float], np.ndarray] = {}
    for al, config in configs.items():
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            rows = sorted(
                [
                    row
                    for row in points
                    if int(row["aggregation_level"]) == al
                    and str(row["candidate_id"]) == candidate_id
                ],
                key=lambda row: float(row["snr_db"]),
            )
            for target in TARGETS:
                key = (al, candidate_id, target)
                try:
                    bracket_lo, bracket_hi = c300.qualified_raw_bracket(rows, target)
                    estimate, ci_lo, ci_hi, replicates = common.bootstrap_targets(
                        rows,
                        target,
                        repeats,
                        common.stable_seed(
                            20260910, al, candidate_id, target, "c300-ideal-bootstrap"
                        ),
                    )
                    status = "estimated"
                    method = "Jeffreys-smoothed weighted decreasing isotonic curve; local log10-BLER interpolation; raw adjacent bracket <=0.25 dB"
                except ValueError as error:
                    bracket_lo = bracket_hi = estimate = ci_lo = ci_hi = None
                    replicates = np.asarray([], dtype=np.float64)
                    status = "unqualified"
                    method = f"not estimated: {error}"
                replicate_map[key] = replicates
                rows_out.append(
                    {
                        "aggregation_level": al,
                        "candidate_id": candidate_id,
                        "target_bler": target,
                        "target_snr_db": estimate,
                        "target_snr_ci95_lo_db": ci_lo,
                        "target_snr_ci95_hi_db": ci_hi,
                        "raw_bracket_lo_db": bracket_lo,
                        "raw_bracket_hi_db": bracket_hi,
                        "raw_bracket_width_db": (
                            None if bracket_lo is None else float(bracket_hi) - float(bracket_lo)
                        ),
                        "status": status,
                        "method": method,
                        "ci_method": (
                            "independent per-point Bernoulli bootstrap"
                            if status == "estimated"
                            else "not applicable"
                        ),
                        "bootstrap_repeats_requested": repeats,
                        "bootstrap_repeats_valid": len(replicates),
                    }
                )
    lookup = {
        (int(row["aggregation_level"]), str(row["candidate_id"]), float(row["target_bler"])): row
        for row in rows_out
    }
    gains: list[dict[str, object]] = []
    for al, config in configs.items():
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            for baseline in BASELINES:
                for target in TARGETS:
                    key = (al, candidate_id, target)
                    base_key = (al, baseline, target)
                    candidate_rep = replicate_map[key]
                    baseline_rep = replicate_map[base_key]
                    count = min(len(candidate_rep), len(baseline_rep))
                    if count:
                        samples = baseline_rep[:count] - candidate_rep[:count]
                        lo, hi = (float(value) for value in np.quantile(samples, [0.025, 0.975]))
                        gain = float(lookup[base_key]["target_snr_db"]) - float(
                            lookup[key]["target_snr_db"]
                        )
                        status = "estimated"
                    else:
                        gain = lo = hi = None
                        status = "unavailable_unqualified_target"
                    gains.append(
                        {
                            "aggregation_level": al,
                            "candidate_id": candidate_id,
                            "baseline_id": baseline,
                            "target_bler": target,
                            "gain_db": gain,
                            "gain_ci95_lo_db": lo,
                            "gain_ci95_hi_db": hi,
                            "status": status,
                            "positive_means_candidate_requires_less_snr": True,
                            "bootstrap_repeats_valid": count,
                        }
                    )
    return rows_out, gains, replicate_map


def _estimated_rows_by_candidate(points: list[dict[str, str]], al: int, candidate_id: str):
    return sorted(
        [
            row
            for row in points
            if int(row["aggregation_level"]) == al
            and str(row["candidate_id"]) == candidate_id
        ],
        key=lambda row: float(row["snr_db"]),
    )


def _ce_penalties(
    ideal_targets: list[dict[str, object]],
    ideal_replicates: dict[tuple[int, str, float], np.ndarray],
    repeats: int,
) -> list[dict[str, object]]:
    estimated_targets = _read_csv(ESTIMATED_ANALYSIS / "target_snr.csv")
    estimated_points = _read_csv(ESTIMATED_ANALYSIS / "formal_points.csv")
    estimated_lookup = {
        (int(row["aggregation_level"]), str(row["candidate_id"]), float(row["target_bler"])): row
        for row in estimated_targets
    }
    ideal_lookup = {
        (int(row["aggregation_level"]), str(row["candidate_id"]), float(row["target_bler"])): row
        for row in ideal_targets
    }
    rows_out: list[dict[str, object]] = []
    for ideal_key, ideal_row in ideal_lookup.items():
        al, ideal_id, target = ideal_key
        estimated_ids = (
            (
                "C300_SMALL_CDD_QSTEP0P25_MATCHED_CDD",
                "C300_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD",
            )
            if ideal_id == "C300_SMALL_CDD_QSTEP0P25_IDEAL_CSI"
            else (ideal_id,)
        )
        for estimated_id in estimated_ids:
            estimated_key = (al, estimated_id, target)
            estimated_row = estimated_lookup[estimated_key]
            if ideal_row["status"] != "estimated" or estimated_row["status"] != "estimated":
                rows_out.append(
                    {
                        "aggregation_level": al,
                        "ideal_candidate_id": ideal_id,
                        "estimated_candidate_id": estimated_id,
                        "target_bler": target,
                        "delta_ce_db": None,
                        "delta_ce_ci95_lo_db": None,
                        "delta_ce_ci95_hi_db": None,
                        "status": "unavailable_unqualified_target",
                        "positive_means_estimated_csi_requires_more_snr": True,
                        "bootstrap_repeats_valid": 0,
                    }
                )
                continue
            estimated_point_rows = _estimated_rows_by_candidate(
                estimated_points, al, estimated_id
            )
            _, _, _, estimated_rep = common.bootstrap_targets(
                estimated_point_rows,
                target,
                repeats,
                common.stable_seed(20260907, al, estimated_id, target, "c300-bootstrap"),
            )
            ideal_rep = ideal_replicates[ideal_key]
            count = min(len(estimated_rep), len(ideal_rep))
            samples = estimated_rep[:count] - ideal_rep[:count]
            lo, hi = (float(value) for value in np.quantile(samples, [0.025, 0.975]))
            rows_out.append(
                {
                    "aggregation_level": al,
                    "ideal_candidate_id": ideal_id,
                    "estimated_candidate_id": estimated_id,
                    "target_bler": target,
                    "delta_ce_db": float(estimated_row["target_snr_db"])
                    - float(ideal_row["target_snr_db"]),
                    "delta_ce_ci95_lo_db": lo,
                    "delta_ce_ci95_hi_db": hi,
                    "status": "estimated",
                    "positive_means_estimated_csi_requires_more_snr": True,
                    "bootstrap_repeats_valid": count,
                }
            )
    return rows_out


def _plot(configs: dict[int, dict], points: list[dict[str, object]]) -> None:
    import matplotlib.pyplot as plt

    FIGURES.mkdir(parents=True, exist_ok=True)
    for al, config in configs.items():
        fig, axis = plt.subplots(figsize=(15.0, 8.5))
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            rows = sorted(
                [
                    row
                    for row in points
                    if int(row["aggregation_level"]) == al
                    and str(row["candidate_id"]) == candidate_id
                ],
                key=lambda row: float(row["snr_db"]),
            )
            style = candidate.get("style", {})
            axis.semilogy(
                [float(row["snr_db"]) for row in rows],
                [
                    float(row["bler"])
                    if float(row["bler"]) > 0.0
                    else 0.5 / int(row["trials"])
                    for row in rows
                ],
                label=str(candidate.get("label", candidate_id)),
                color=style.get("color"),
                linestyle=style.get("linestyle", "-"),
                marker=style.get("marker", "o"),
                linewidth=2.5,
                markersize=8,
            )
        axis.axhline(0.10, color="black", linestyle=":", linewidth=1.8)
        axis.axhline(0.01, color="black", linestyle="--", linewidth=1.8)
        axis.set(xlabel="SNR (dB)", ylabel="Ideal-CSI DCI BLER")
        axis.set_ylim(bottom=0.005)
        axis.set_title(f"C300 two-symbol ideal-CSI PDCCH AL{al}", fontsize=17)
        axis.xaxis.label.set_size(16)
        axis.yaxis.label.set_size(16)
        axis.tick_params(labelsize=14)
        axis.grid(True, which="both", alpha=0.35)
        axis.legend(fontsize=16, loc="center left", bbox_to_anchor=(1.01, 0.5))
        fig.tight_layout()
        fig.savefig(FIGURES / f"al{al}_ideal_csi_bler.png", dpi=180)
        plt.close(fig)


def _environment_receipt() -> dict[str, object]:
    import matplotlib
    import numpy
    import scipy
    import sionna
    import tensorflow

    return {
        "platform": platform.platform(),
        "python": platform.python_version(),
        "packages": {
            "numpy": numpy.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
            "tensorflow": tensorflow.__version__,
            "sionna": sionna.__version__,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bootstrap-repeats", type=int, default=2000)
    args = parser.parse_args()
    common.CONFIGS = CONFIGS
    common.OUTPUT = OUTPUT
    common.ANALYSIS = ANALYSIS
    common.FIGURES = FIGURES
    configs, points, diagnostics = common.collect()
    targets, gains, replicate_map = _analyze_targets(
        configs, points, int(args.bootstrap_repeats)
    )
    penalties = _ce_penalties(targets, replicate_map, int(args.bootstrap_repeats))
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    common.write_csv(ANALYSIS / "formal_points.csv", points)
    common.write_csv(ANALYSIS / "diagnostics.csv", diagnostics)
    common.write_csv(ANALYSIS / "target_snr.csv", targets)
    common.write_csv(ANALYSIS / "target_gains.csv", gains)
    common.write_csv(ANALYSIS / "delta_ce.csv", penalties)
    _plot(configs, points)
    git_head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True
    ).stdout.strip()
    metadata = {
        "schema": "plan031-c300-ideal-analysis-v1",
        "plan": "research/plan-031-PDCCH-CDD时延-BLER.md#2026-09-10-补充规划c300-2-symbol-场景-ideal-csi-bler待执行",
        "bootstrap_repeats": int(args.bootstrap_repeats),
        "git_head": git_head,
        "working_tree_note": "plan-031 ideal-CSI implementation and results are uncommitted",
        "config_sha256": {
            path.relative_to(ROOT).as_posix(): common.sha256(path) for path in CONFIGS.values()
        },
        "formal_point_count": len(points),
        "total_trials": sum(int(row["trials"]) for row in points),
        "total_errors": sum(int(row["errors"]) for row in points),
        "max_abs_ce_nmse": max(abs(float(row["ce_nmse_db"]) + 300.0) for row in points),
        "target_rows": len(targets),
        "qualified_target_rows": sum(row["status"] == "estimated" for row in targets),
        "delta_ce_rows": len(penalties),
        "qualified_delta_ce_rows": sum(row["status"] == "estimated" for row in penalties),
    }
    (ANALYSIS / "analysis_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (OUTPUT / "environment_receipt.json").write_text(
        json.dumps(_environment_receipt(), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(metadata, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
