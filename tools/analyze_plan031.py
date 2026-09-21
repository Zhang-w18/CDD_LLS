"""Validate and analyze the formal plan-031 PDCCH CDD runs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
from typing import Iterable

import numpy as np
import yaml


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs" / "experiment031_pdcch_cdd" / "20260903_main"
ANALYSIS = OUTPUT / "analysis"
FIGURES = ROOT / "docs" / "figures" / "result-031"
CONFIGS = {
    2: ROOT / "configs" / "pdcch_result031_al2_formal.yaml",
    4: ROOT / "configs" / "pdcch_result031_al4_formal.yaml",
    8: ROOT / "configs" / "pdcch_result031_al8_formal.yaml",
}
TARGETS = (0.10, 0.01)
BASELINES = ("B0_QC", "A100_PRG_DFT8_6RB")


def stable_seed(*items: object) -> int:
    text = "|".join(str(item) for item in items)
    acc = 2166136261
    for char in text:
        acc ^= ord(char)
        acc = (acc * 16777619) % (2**32)
    return int(acc)


def read_csv(path: Path) -> list[dict[str, str]]:
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"Refusing to write empty CSV {path}.")
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def isotonic_decreasing(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    y = np.asarray(values, dtype=np.float64)
    w = np.asarray(weights, dtype=np.float64)
    blocks: list[list[float | int]] = []
    for index, (value, weight) in enumerate(zip(y, w)):
        blocks.append([index, index + 1, float(value), float(weight)])
        while len(blocks) >= 2 and float(blocks[-2][2]) < float(blocks[-1][2]):
            right = blocks.pop()
            left = blocks.pop()
            total_weight = float(left[3]) + float(right[3])
            mean = (
                float(left[2]) * float(left[3]) + float(right[2]) * float(right[3])
            ) / total_weight
            blocks.append([int(left[0]), int(right[1]), mean, total_weight])
    out = np.empty_like(y)
    for start, stop, value, _ in blocks:
        out[int(start) : int(stop)] = float(value)
    return out


def interpolate_target(snr: np.ndarray, probabilities: np.ndarray, target: float) -> tuple[float, int]:
    x = np.asarray(snr, dtype=np.float64)
    p = np.asarray(probabilities, dtype=np.float64)
    if p[0] < target or p[-1] > target:
        raise ValueError("Target is not double-sided bracketed.")
    for index in range(len(x) - 1):
        if p[index] >= target >= p[index + 1]:
            if abs(p[index] - p[index + 1]) <= 1e-15:
                return float(0.5 * (x[index] + x[index + 1])), index
            fraction = (
                math.log10(target) - math.log10(p[index])
            ) / (math.log10(p[index + 1]) - math.log10(p[index]))
            return float(x[index] + fraction * (x[index + 1] - x[index])), index
    raise ValueError("No local target bracket was found.")


def raw_bracket(rows: list[dict[str, object]], target: float) -> tuple[float, float]:
    for left, right in zip(rows[:-1], rows[1:]):
        p_left = float(left["bler"])
        p_right = float(right["bler"])
        if p_left >= target >= p_right:
            return float(left["snr_db"]), float(right["snr_db"])
    raise ValueError(f"Raw BLER does not provide an adjacent bracket for target {target}.")


def bootstrap_targets(
    rows: list[dict[str, object]], target: float, repeats: int, seed: int
) -> tuple[float, float, float, np.ndarray]:
    snr = np.asarray([float(row["snr_db"]) for row in rows], dtype=np.float64)
    trials = np.asarray([int(row["trials"]) for row in rows], dtype=np.int64)
    errors = np.asarray([int(row["errors"]) for row in rows], dtype=np.int64)
    probabilities = (errors + 0.5) / (trials + 1.0)
    fitted = isotonic_decreasing(probabilities, trials)
    estimate, _ = interpolate_target(snr, fitted, target)
    rng = np.random.default_rng(int(seed))
    replicates: list[float] = []
    raw_probability = errors / trials
    for _ in range(int(repeats)):
        sampled_errors = rng.binomial(trials, raw_probability)
        sampled = (sampled_errors + 0.5) / (trials + 1.0)
        sampled_fit = isotonic_decreasing(sampled, trials)
        try:
            value, _ = interpolate_target(snr, sampled_fit, target)
        except ValueError:
            continue
        replicates.append(value)
    if len(replicates) < int(0.95 * repeats):
        raise ValueError(
            f"Only {len(replicates)}/{repeats} bootstrap target replicates were bracketed."
        )
    array = np.asarray(replicates, dtype=np.float64)
    low, high = np.quantile(array, [0.025, 0.975])
    return estimate, float(low), float(high), array


def _flag_stem(candidate_id: str, snr_db: float) -> str:
    safe = "".join(char if char.isalnum() or char in "-_" else "_" for char in candidate_id)
    tag = str(float(snr_db)).replace("-", "m").replace(".", "p")
    return f"{safe}_snr_{tag}"


def collect() -> tuple[dict[int, dict], list[dict[str, object]], list[dict[str, object]]]:
    configs: dict[int, dict] = {}
    points: list[dict[str, object]] = []
    diagnostics: list[dict[str, object]] = []
    for aggregation_level, config_path in CONFIGS.items():
        with open(config_path, "r", encoding="utf-8") as handle:
            config = yaml.safe_load(handle) or {}
        configs[aggregation_level] = config
        default_snr = [float(value) for value in config["simulation"]["snr_points_db"]]
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            expected_snr = [
                float(value) for value in candidate.get("snr_points_db", default_snr)
            ]
            candidate_dir = ROOT / config["output_dir"] / candidate_id
            csv_path = candidate_dir / "bler_points.csv"
            metadata_path = candidate_dir / "run_metadata.json"
            if not csv_path.exists() or not metadata_path.exists():
                raise FileNotFoundError(f"Incomplete formal shard: {candidate_dir}")
            rows = read_csv(csv_path)
            actual_snr = [float(row["snr_db"]) for row in rows]
            if sorted(actual_snr) != sorted(expected_snr) or len(actual_snr) != len(set(actual_snr)):
                raise ValueError(f"{candidate_id} AL{aggregation_level}: incomplete or duplicate SNR grid.")
            rows = sorted(rows, key=lambda row: float(row["snr_db"]))
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            meta_candidate = metadata["candidates"][0]
            receiver_mode = str(
                meta_candidate.get(
                    "receiver_covariance_mode",
                    candidate.get(
                        "receiver_covariance_mode",
                        (
                            "physical_prg"
                            if str(candidate.get("scheme", "")).lower()
                            == "reg_bundle_dft_cycling"
                            else "matched_effective"
                        ),
                    ),
                )
            )
            violations = 0
            previous = None
            total_trials = 0
            total_errors = 0
            max_flag_error = 0
            max_nmse_error_db = 0.0
            for row in rows:
                normalized: dict[str, object] = dict(row)
                for key in (
                    "snr_db",
                    "bler",
                    "bler_wilson95_lo",
                    "bler_wilson95_hi",
                    "ce_nmse_db",
                    "pilot_condition_number",
                    "ce_floor_nmse_db",
                ):
                    normalized[key] = float(row[key])
                for key in (
                    "trials",
                    "errors",
                    "aggregation_level",
                    "occupied_rb",
                    "k_active",
                    "data_re",
                    "dmrs_re",
                    "coded_bits",
                    "pilot_rank",
                ):
                    normalized[key] = int(float(row[key]))
                normalized["source_csv"] = csv_path.relative_to(ROOT).as_posix()
                normalized["receiver_covariance_mode"] = receiver_mode
                points.append(normalized)
                if previous is not None and float(row["bler"]) > previous + 1e-15:
                    violations += 1
                previous = float(row["bler"])
                total_trials += int(row["trials"])
                total_errors += int(row["errors"])
                stem = _flag_stem(candidate_id, float(row["snr_db"]))
                flags = np.load(candidate_dir / "trial_error_flags" / f"{stem}_error_flags.npy")
                ce = np.load(candidate_dir / "trial_error_flags" / f"{stem}_ce_nmse.npy")
                max_flag_error = max(
                    max_flag_error,
                    abs(len(flags) - int(row["trials"])),
                    abs(int(np.sum(flags)) - int(row["errors"])),
                )
                ce_db = 10.0 * math.log10(max(float(np.mean(ce)), 1e-30))
                max_nmse_error_db = max(max_nmse_error_db, abs(ce_db - float(row["ce_nmse_db"])))
            diagnostics.append(
                {
                    "aggregation_level": aggregation_level,
                    "candidate_id": candidate_id,
                    "receiver_covariance_mode": receiver_mode,
                    "pilot_rank": int(meta_candidate["pilot_rank"]),
                    "pilot_condition_number": float(meta_candidate["pilot_condition_number"]),
                    "ce_floor_nmse_db": 10.0
                    * math.log10(max(float(meta_candidate["ce_floor_nmse"]), 1e-30)),
                    "formal_points": len(rows),
                    "total_trials": total_trials,
                    "total_errors": total_errors,
                    "adjacent_bler_increases": violations,
                    "max_flag_count_or_sum_error": max_flag_error,
                    "max_ce_recompute_error_db": max_nmse_error_db,
                    "reached_max_without_200_errors": sum(
                        int(row["trials"]) == int(config["simulation"]["max_trials_per_snr"])
                        and int(row["errors"]) < int(config["simulation"]["target_errors"])
                        for row in rows
                    ),
                }
            )
    return configs, points, diagnostics


def analyze(
    configs: dict[int, dict],
    points: list[dict[str, object]],
    bootstrap_repeats: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    target_rows: list[dict[str, object]] = []
    replicate_map: dict[tuple[int, str, float], np.ndarray] = {}
    for aggregation_level, config in configs.items():
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            rows = sorted(
                [
                    row
                    for row in points
                    if int(row["aggregation_level"]) == aggregation_level
                    and str(row["candidate_id"]) == candidate_id
                ],
                key=lambda row: float(row["snr_db"]),
            )
            for target in TARGETS:
                try:
                    bracket_low, bracket_high = raw_bracket(rows, target)
                except ValueError:
                    replicate_map[(aggregation_level, candidate_id, target)] = np.asarray(
                        [], dtype=np.float64
                    )
                    target_rows.append(
                        {
                            "aggregation_level": aggregation_level,
                            "candidate_id": candidate_id,
                            "target_bler": target,
                            "target_snr_db": None,
                            "target_snr_ci95_lo_db": None,
                            "target_snr_ci95_hi_db": None,
                            "raw_bracket_lo_db": None,
                            "raw_bracket_hi_db": None,
                            "raw_bracket_width_db": None,
                            "status": "unbracketed",
                            "method": "not estimated; no adjacent raw double-sided bracket",
                            "ci_method": "not applicable",
                            "bootstrap_repeats_requested": bootstrap_repeats,
                            "bootstrap_repeats_valid": 0,
                        }
                    )
                    continue
                estimate, low, high, replicates = bootstrap_targets(
                    rows,
                    target,
                    bootstrap_repeats,
                    stable_seed(20260904, aggregation_level, candidate_id, target, "bootstrap"),
                )
                replicate_map[(aggregation_level, candidate_id, target)] = replicates
                target_rows.append(
                    {
                        "aggregation_level": aggregation_level,
                        "candidate_id": candidate_id,
                        "target_bler": target,
                        "target_snr_db": estimate,
                        "target_snr_ci95_lo_db": low,
                        "target_snr_ci95_hi_db": high,
                        "raw_bracket_lo_db": bracket_low,
                        "raw_bracket_hi_db": bracket_high,
                        "raw_bracket_width_db": bracket_high - bracket_low,
                        "status": "estimated",
                        "method": "Jeffreys-smoothed weighted decreasing isotonic curve; local log10-BLER interpolation",
                        "ci_method": "independent per-point Bernoulli bootstrap",
                        "bootstrap_repeats_requested": bootstrap_repeats,
                        "bootstrap_repeats_valid": len(replicates),
                    }
                )
    target_lookup = {
        (int(row["aggregation_level"]), str(row["candidate_id"]), float(row["target_bler"])): row
        for row in target_rows
    }
    gain_rows: list[dict[str, object]] = []
    for aggregation_level, config in configs.items():
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            for baseline in BASELINES:
                for target in TARGETS:
                    base_row = target_lookup[(aggregation_level, baseline, target)]
                    candidate_row = target_lookup[(aggregation_level, candidate_id, target)]
                    base_rep = replicate_map[(aggregation_level, baseline, target)]
                    candidate_rep = replicate_map[(aggregation_level, candidate_id, target)]
                    count = min(len(base_rep), len(candidate_rep))
                    if count == 0:
                        gain_rows.append(
                            {
                                "aggregation_level": aggregation_level,
                                "candidate_id": candidate_id,
                                "baseline_id": baseline,
                                "target_bler": target,
                                "gain_db": None,
                                "gain_ci95_lo_db": None,
                                "gain_ci95_hi_db": None,
                                "status": "unavailable_unbracketed_target",
                                "positive_means_candidate_requires_less_snr": True,
                                "bootstrap_repeats_valid": 0,
                            }
                        )
                        continue
                    gains = base_rep[:count] - candidate_rep[:count]
                    low, high = np.quantile(gains, [0.025, 0.975])
                    gain_rows.append(
                        {
                            "aggregation_level": aggregation_level,
                            "candidate_id": candidate_id,
                            "baseline_id": baseline,
                            "target_bler": target,
                            "gain_db": float(base_row["target_snr_db"])
                            - float(candidate_row["target_snr_db"]),
                            "gain_ci95_lo_db": float(low),
                            "gain_ci95_hi_db": float(high),
                            "status": "estimated",
                            "positive_means_candidate_requires_less_snr": True,
                            "bootstrap_repeats_valid": count,
                        }
                    )
    return target_rows, gain_rows


def plot(configs: dict[int, dict], points: list[dict[str, object]]) -> None:
    import matplotlib.pyplot as plt

    FIGURES.mkdir(parents=True, exist_ok=True)
    for aggregation_level, config in configs.items():
        for metric, ylabel, stem, log_scale in (
            ("bler", "Estimated-CSI DCI BLER", f"al{aggregation_level}_bler", True),
            ("ce_nmse_db", "CE NMSE (dB)", f"al{aggregation_level}_ce_nmse", False),
        ):
            fig, axis = plt.subplots(figsize=(15.0, 8.5))
            for candidate in config["candidates"]:
                candidate_id = str(candidate["candidate_id"])
                rows = sorted(
                    [
                        row
                        for row in points
                        if int(row["aggregation_level"]) == aggregation_level
                        and str(row["candidate_id"]) == candidate_id
                    ],
                    key=lambda row: float(row["snr_db"]),
                )
                style = candidate.get("style", {})
                x = [float(row["snr_db"]) for row in rows]
                if metric == "bler":
                    y = [
                        float(row[metric])
                        if float(row[metric]) > 0.0
                        else 0.5 / int(row["trials"])
                        for row in rows
                    ]
                else:
                    y = [float(row[metric]) for row in rows]
                axis.plot(
                    x,
                    y,
                    label=str(candidate.get("label", candidate_id)),
                    color=style.get("color"),
                    linestyle=style.get("linestyle", "-"),
                    marker=style.get("marker", "o"),
                    linewidth=2.5,
                    markersize=8,
                )
            if log_scale:
                axis.set_yscale("log")
                axis.axhline(0.10, color="black", linestyle=":", linewidth=1.8)
                axis.axhline(0.01, color="black", linestyle="--", linewidth=1.8)
                axis.set_ylim(bottom=0.01)
            axis.set_xlabel("SNR (dB)", fontsize=16)
            axis.set_ylabel(ylabel, fontsize=16)
            axis.set_title(f"PDCCH AL{aggregation_level}", fontsize=17)
            axis.grid(True, which="both", alpha=0.35)
            axis.tick_params(labelsize=14)
            axis.legend(
                fontsize=16,
                ncol=1,
                loc="center left",
                bbox_to_anchor=(1.01, 0.5),
            )
            fig.tight_layout()
            fig.savefig(FIGURES / f"{stem}.png", dpi=180)
            plt.close(fig)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bootstrap-repeats", type=int, default=2000)
    args = parser.parse_args()
    configs, points, diagnostics = collect()
    targets, gains = analyze(configs, points, int(args.bootstrap_repeats))
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    write_csv(ANALYSIS / "formal_points.csv", points)
    write_csv(ANALYSIS / "diagnostics.csv", diagnostics)
    write_csv(ANALYSIS / "target_snr.csv", targets)
    write_csv(ANALYSIS / "target_gains.csv", gains)
    plot(configs, points)
    git_head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True
    ).stdout.strip()
    summary = {
        "schema": "plan031-analysis-v1",
        "plan": "research/plan-031-PDCCH-CDD时延-BLER.md",
        "bootstrap_repeats": int(args.bootstrap_repeats),
        "git_head": git_head,
        "working_tree_note": "plan-031 implementation and result files are uncommitted",
        "config_sha256": {
            path.relative_to(ROOT).as_posix(): sha256(path) for path in CONFIGS.values()
        },
        "formal_point_count": len(points),
        "total_trials": sum(int(row["trials"]) for row in points),
        "total_errors": sum(int(row["errors"]) for row in points),
        "target_rows": len(targets),
        "gain_rows": len(gains),
    }
    (ANALYSIS / "analysis_metadata.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
