"""Audit and analyze the final plan-031 strict-Sidon fine scan."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml

from analyze_plan031 import interpolate_target, isotonic_decreasing, raw_bracket, stable_seed


ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search"
DATA_ROOT = RUN_ROOT / "fine"
ANALYSIS_ROOT = RUN_ROOT / "fine_analysis"
FIGURE_ROOT = ROOT / "docs/figures/result-031"
CONFIGS = {
    scene: ROOT / "configs" / f"pdcch_result031_strict_sidon_{scene}_fine.yaml"
    for scene in (
        "a100_al2",
        "a100_al4",
        "a100_al8",
        "c300_al1",
        "c300_al2",
        "c300_al4",
    )
}
SUPPLEMENT_CONFIGS = {
    scene: ROOT
    / "configs"
    / f"pdcch_result031_strict_sidon_{scene}_fine_supplement.yaml"
    for scene in ("a100_al4", "a100_al8", "c300_al2", "c300_al4")
}
LEGACY_S0 = {
    "a100_al2": (
        ROOT / "outputs/experiment031_pdcch_cdd/20260903_main/formal/al2/S0_SIDON/bler_points.csv",
        "S0_SIDON",
    ),
    "c300_al1": (
        ROOT / "outputs/experiment031_pdcch_cdd/20260907_c300_4tx_2sym/formal/al1/C300_S0_SIDON/bler_points.csv",
        "C300_S0_SIDON",
    ),
}
TARGETS = (0.10, 0.01)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _flag_stem(candidate_id: str, snr_db: float) -> str:
    tag = str(float(snr_db)).replace("-", "m").replace(".", "p")
    return f"{candidate_id}_snr_{tag}"


def _normalize_rows(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    normalized_rows = []
    for source in rows:
        row: dict[str, object] = dict(source)
        for key in (
            "snr_db",
            "bler",
            "bler_wilson95_lo",
            "bler_wilson95_hi",
            "ce_nmse_db",
            "pilot_condition_number",
            "ce_floor_nmse_db",
        ):
            row[key] = float(source[key])
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
            row[key] = int(float(source[key]))
        normalized_rows.append(row)
    return sorted(normalized_rows, key=lambda row: float(row["snr_db"]))


def _bootstrap_target(
    rows: list[dict], target: float, repeats: int, seed: int
) -> tuple[float, float, float, np.ndarray]:
    snr = np.asarray([float(row["snr_db"]) for row in rows])
    trials = np.asarray([int(row["trials"]) for row in rows], dtype=np.int64)
    errors = np.asarray([int(row["errors"]) for row in rows], dtype=np.int64)
    fitted = isotonic_decreasing((errors + 0.5) / (trials + 1.0), trials)
    estimate, _ = interpolate_target(snr, fitted, target)
    rng = np.random.default_rng(seed)
    replicates = []
    raw = errors / trials
    for _ in range(repeats):
        sampled = rng.binomial(trials, raw)
        sampled_fit = isotonic_decreasing((sampled + 0.5) / (trials + 1.0), trials)
        try:
            value, _ = interpolate_target(snr, sampled_fit, target)
        except ValueError:
            continue
        replicates.append(value)
    array = np.asarray(replicates, dtype=np.float64)
    if not len(array):
        return estimate, math.nan, math.nan, array
    low, high = np.quantile(array, [0.025, 0.975])
    return estimate, float(low), float(high), array


def collect() -> tuple[dict[str, dict], dict[str, list[dict]], list[dict], list[dict]]:
    configs: dict[str, dict] = {}
    scene_candidates: dict[str, list[dict]] = {}
    all_points: list[dict] = []
    diagnostics: list[dict] = []
    supplement_points: dict[tuple[str, str], list[float]] = {}
    for scene, config_path in SUPPLEMENT_CONFIGS.items():
        config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
        for candidate in config["candidates"]:
            supplement_points[(scene, str(candidate["candidate_id"]))] = sorted(
                float(value) for value in candidate["snr_points_db"]
            )
    for scene, config_path in CONFIGS.items():
        config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
        configs[scene] = config
        default_snr = [float(value) for value in config["simulation"]["snr_points_db"]]
        scene_candidates[scene] = []
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            expected_snr = sorted(
                float(value) for value in candidate.get("snr_points_db", default_snr)
            )
            candidate_dir = DATA_ROOT / scene / candidate_id
            csv_path = candidate_dir / "bler_points.csv"
            metadata_path = candidate_dir / "run_metadata.json"
            if not csv_path.exists() or not metadata_path.exists():
                raise FileNotFoundError(f"Incomplete final candidate output: {candidate_dir}")
            rows = _normalize_rows(_read_csv(csv_path))
            actual_snr = [float(row["snr_db"]) for row in rows]
            if actual_snr != expected_snr or len(actual_snr) != len(set(actual_snr)):
                raise RuntimeError(f"Fine grid mismatch for {candidate_id}")
            source_by_snr = {
                float(row["snr_db"]): (candidate_dir, csv_path, "fine") for row in rows
            }
            extra_expected = supplement_points.get((scene, candidate_id), [])
            if extra_expected:
                supplement_dir = RUN_ROOT / "fine_supplement" / scene / candidate_id
                supplement_csv = supplement_dir / "bler_points.csv"
                supplement_metadata = supplement_dir / "run_metadata.json"
                if not supplement_csv.exists() or not supplement_metadata.exists():
                    raise FileNotFoundError(
                        f"Incomplete fine supplement output: {supplement_dir}"
                    )
                extra_rows = _normalize_rows(_read_csv(supplement_csv))
                extra_actual = [float(row["snr_db"]) for row in extra_rows]
                if extra_actual != extra_expected:
                    raise RuntimeError(f"Supplement grid mismatch for {candidate_id}")
                if set(actual_snr) & set(extra_actual):
                    raise RuntimeError(f"Fine supplement duplicates an existing SNR for {candidate_id}")
                rows = sorted(rows + extra_rows, key=lambda row: float(row["snr_db"]))
                source_by_snr.update(
                    {
                        float(row["snr_db"]): (
                            supplement_dir,
                            supplement_csv,
                            "fine_supplement",
                        )
                        for row in extra_rows
                    }
                )
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            meta_candidate = metadata["candidates"][0]
            max_evidence_delta = 0
            max_ce_delta_db = 0.0
            nmse_rows = []
            for row in rows:
                row_dir, row_csv, source_stage = source_by_snr[float(row["snr_db"])]
                stem = _flag_stem(candidate_id, float(row["snr_db"]))
                trial_dir = row_dir / "trial_error_flags"
                flags = np.load(trial_dir / f"{stem}_error_flags.npy")
                ce = np.load(trial_dir / f"{stem}_ce_nmse.npy")
                trials = int(row["trials"])
                max_evidence_delta = max(
                    max_evidence_delta,
                    abs(len(flags) - trials),
                    abs(len(ce) - trials),
                    abs(int(np.count_nonzero(flags)) - int(row["errors"])),
                )
                ce_mean_db = 10.0 * math.log10(max(float(np.mean(ce)), 1e-30))
                max_ce_delta_db = max(
                    max_ce_delta_db, abs(ce_mean_db - float(row["ce_nmse_db"]))
                )
                quantiles = np.quantile(ce, [0.025, 0.10, 0.50, 0.90, 0.975])
                nmse_rows.append(
                    {
                        "scene": scene,
                        "candidate_id": candidate_id,
                        "snr_db": row["snr_db"],
                        "trials": trials,
                        "source_stage": source_stage,
                        "mean_linear": float(np.mean(ce)),
                        "mean_db": ce_mean_db,
                        "q02p5_db": 10.0 * math.log10(max(float(quantiles[0]), 1e-30)),
                        "q10_db": 10.0 * math.log10(max(float(quantiles[1]), 1e-30)),
                        "median_db": 10.0 * math.log10(max(float(quantiles[2]), 1e-30)),
                        "q90_db": 10.0 * math.log10(max(float(quantiles[3]), 1e-30)),
                        "q97p5_db": 10.0 * math.log10(max(float(quantiles[4]), 1e-30)),
                    }
                )
                all_points.append(
                    {
                        "scene": scene,
                        "source_stage": source_stage,
                        "source_csv": row_csv.relative_to(ROOT).as_posix(),
                        **row,
                    }
                )
            if max_evidence_delta != 0 or max_ce_delta_db > 1e-10:
                raise RuntimeError(f"Per-trial evidence mismatch for {candidate_id}")
            bler = [float(row["bler"]) for row in rows]
            floor_db = 10.0 * math.log10(
                max(float(meta_candidate["ce_floor_nmse"]), 1e-30)
            )
            diagnostics.append(
                {
                    "scene": scene,
                    "candidate_id": candidate_id,
                    "pilot_rank": int(meta_candidate["pilot_rank"]),
                    "pilot_condition_number": float(meta_candidate["pilot_condition_number"]),
                    "ce_floor_nmse_db": floor_db,
                    "highest_snr_ce_nmse_db": float(rows[-1]["ce_nmse_db"]),
                    "highest_snr_ce_minus_floor_db": float(rows[-1]["ce_nmse_db"]) - floor_db,
                    "fine_points": len(rows),
                    "total_trials": sum(int(row["trials"]) for row in rows),
                    "total_errors": sum(int(row["errors"]) for row in rows),
                    "adjacent_raw_bler_increases": sum(
                        right > left + 1e-15 for left, right in zip(bler[:-1], bler[1:])
                    ),
                    "max_trials_points_below_200_errors": sum(
                        int(row["trials"]) == 50_000 and int(row["errors"]) < 200
                        for row in rows
                    ),
                    "max_flag_length_or_sum_delta": max_evidence_delta,
                    "max_ce_recompute_delta_db": max_ce_delta_db,
                }
            )
            scene_candidates[scene].append(
                {"candidate_id": candidate_id, "config": candidate, "rows": rows, "nmse": nmse_rows}
            )
    if sum(len(values) for values in scene_candidates.values()) != 18:
        raise RuntimeError("Expected 18 final candidates")
    if len(all_points) != 172:
        raise RuntimeError("Expected 149 fine plus 23 supplement candidate/SNR points")
    return configs, scene_candidates, all_points, diagnostics


def analyze_targets(
    scene_candidates: dict[str, list[dict]], repeats: int
) -> tuple[list[dict], dict[tuple[str, str, float], np.ndarray], list[dict]]:
    target_rows = []
    samples: dict[tuple[str, str, float], np.ndarray] = {}
    bracket_rows = []
    for scene, candidates in sorted(scene_candidates.items()):
        for candidate in candidates:
            candidate_id = candidate["candidate_id"]
            rows = candidate["rows"]
            for target in TARGETS:
                try:
                    low, high = raw_bracket(rows, target)
                except ValueError:
                    samples[(scene, candidate_id, target)] = np.asarray([], dtype=np.float64)
                    target_rows.append(
                        {
                            "scene": scene,
                            "candidate_id": candidate_id,
                            "target_bler": target,
                            "status": "unbracketed",
                            "target_snr_db": None,
                            "ci95_lo_db": None,
                            "ci95_hi_db": None,
                            "raw_bracket_lo_db": None,
                            "raw_bracket_hi_db": None,
                            "raw_bracket_width_db": None,
                            "bootstrap_valid": 0,
                            "bootstrap_valid_fraction": 0.0,
                        }
                    )
                    continue
                estimate, ci_low, ci_high, replicate = _bootstrap_target(
                    rows,
                    target,
                    repeats,
                    stable_seed(20260910, scene, candidate_id, target, "strict-sidon-fine"),
                )
                if high - low > 0.2500000001:
                    status = "wide_raw_bracket"
                elif len(replicate) < int(0.95 * repeats):
                    status = "boundary_limited_bootstrap"
                else:
                    status = "estimated"
                samples[(scene, candidate_id, target)] = (
                    replicate if status == "estimated" else np.asarray([], dtype=np.float64)
                )
                target_rows.append(
                    {
                        "scene": scene,
                        "candidate_id": candidate_id,
                        "target_bler": target,
                        "status": status,
                        "target_snr_db": estimate,
                        "ci95_lo_db": ci_low,
                        "ci95_hi_db": ci_high,
                        "raw_bracket_lo_db": low,
                        "raw_bracket_hi_db": high,
                        "raw_bracket_width_db": high - low,
                        "bootstrap_valid": len(replicate),
                        "bootstrap_valid_fraction": len(replicate) / repeats,
                    }
                )
                by_snr = {float(row["snr_db"]): row for row in rows}
                for side, snr in (("low", low), ("high", high)):
                    row = by_snr[snr]
                    bracket_rows.append(
                        {
                            "scene": scene,
                            "candidate_id": candidate_id,
                            "target_bler": target,
                            "side": side,
                            "snr_db": snr,
                            "errors": row["errors"],
                            "trials": row["trials"],
                            "bler": row["bler"],
                            "wilson95_lo": row["bler_wilson95_lo"],
                            "wilson95_hi": row["bler_wilson95_hi"],
                        }
                    )
    return target_rows, samples, bracket_rows


def _legacy_s0_samples(scene: str, target: float, repeats: int) -> tuple[float, np.ndarray]:
    path, candidate_id = LEGACY_S0[scene]
    rows = _normalize_rows(_read_csv(path))
    estimate, _, _, replicate = _bootstrap_target(
        rows,
        target,
        repeats,
        stable_seed(20260910, scene, candidate_id, target, "legacy-s0-comparison"),
    )
    if len(replicate) < int(0.95 * repeats):
        raise RuntimeError(f"Legacy S0 target is boundary-limited for {scene}/{target:g}")
    return estimate, replicate


def rank_candidates(
    scene_candidates: dict[str, list[dict]],
    targets: list[dict],
    samples: dict[tuple[str, str, float], np.ndarray],
    repeats: int,
) -> tuple[list[dict], list[dict], list[dict]]:
    lookup = {
        (row["scene"], row["candidate_id"], row["target_bler"]): row for row in targets
    }
    winner_rows = []
    pairwise_rows = []
    s0_rows = []
    for scene, candidates in sorted(scene_candidates.items()):
        ids = [candidate["candidate_id"] for candidate in candidates]
        estimated = [candidate_id for candidate_id in ids if lookup[(scene, candidate_id, 0.01)]["status"] == "estimated"]
        if not estimated:
            raise RuntimeError(f"No bracketed 1% target in {scene}")
        ordered = sorted(estimated, key=lambda value: lookup[(scene, value, 0.01)]["target_snr_db"])
        leader = ordered[0]
        common = min(len(samples[(scene, candidate_id, 0.01)]) for candidate_id in ordered)
        matrix = np.vstack([samples[(scene, candidate_id, 0.01)][:common] for candidate_id in ordered])
        winners = np.argmin(matrix, axis=0)
        co_best = []
        for rank, candidate_id in enumerate(ordered, start=1):
            difference = samples[(scene, candidate_id, 0.01)][:common] - samples[(scene, leader, 0.01)][:common]
            low, high = np.quantile(difference, [0.025, 0.975])
            indistinguishable = bool(low <= 0.0)
            if indistinguishable:
                co_best.append(candidate_id)
            pairwise_rows.append(
                {
                    "scene": scene,
                    "candidate_id": candidate_id,
                    "leader_id": leader,
                    "rank_by_1pct_point_estimate": rank,
                    "candidate_minus_leader_snr_db": float(lookup[(scene, candidate_id, 0.01)]["target_snr_db"]) - float(lookup[(scene, leader, 0.01)]["target_snr_db"]),
                    "difference_ci95_lo_db": float(low),
                    "difference_ci95_hi_db": float(high),
                    "bootstrap_probability_best": float(np.mean(winners == rank - 1)),
                    "statistically_indistinguishable_from_leader_95pct": indistinguishable,
                }
            )
        winner_rows.append(
            {
                "scene": scene,
                "point_estimate_leader": leader,
                "leader_snr_at_1pct_db": lookup[(scene, leader, 0.01)]["target_snr_db"],
                "leader_ci95_lo_db": lookup[(scene, leader, 0.01)]["ci95_lo_db"],
                "leader_ci95_hi_db": lookup[(scene, leader, 0.01)]["ci95_hi_db"],
                "co_best_95pct": ";".join(co_best),
                "co_best_count": len(co_best),
                "candidate_count": len(ids),
                "unbracketed_1pct_count": len(ids) - len(estimated),
            }
        )

        for target in TARGETS:
            if scene in LEGACY_S0:
                s0_estimate, s0_rep = _legacy_s0_samples(scene, target, repeats)
                s0_source = LEGACY_S0[scene][0].relative_to(ROOT).as_posix()
                s0_id = LEGACY_S0[scene][1]
            else:
                s0_id = next(value for value in ids if value.endswith("_01"))
                s0_estimate = float(lookup[(scene, s0_id, target)]["target_snr_db"])
                s0_rep = samples[(scene, s0_id, target)]
                s0_source = "final fine scan candidate 01"
            for candidate_id in ids:
                candidate_row = lookup[(scene, candidate_id, target)]
                candidate_rep = samples[(scene, candidate_id, target)]
                count = min(len(s0_rep), len(candidate_rep))
                if candidate_row["status"] != "estimated" or count == 0:
                    s0_rows.append(
                        {
                            "scene": scene,
                            "candidate_id": candidate_id,
                            "s0_id": s0_id,
                            "target_bler": target,
                            "gain_vs_s0_db": None,
                            "ci95_lo_db": None,
                            "ci95_hi_db": None,
                            "status": "unavailable_unbracketed",
                            "s0_source": s0_source,
                        }
                    )
                    continue
                gain = s0_rep[:count] - candidate_rep[:count]
                low, high = np.quantile(gain, [0.025, 0.975])
                s0_rows.append(
                    {
                        "scene": scene,
                        "candidate_id": candidate_id,
                        "s0_id": s0_id,
                        "target_bler": target,
                        "gain_vs_s0_db": s0_estimate - float(candidate_row["target_snr_db"]),
                        "ci95_lo_db": float(low),
                        "ci95_hi_db": float(high),
                        "status": "estimated",
                        "s0_source": s0_source,
                    }
                )
    return winner_rows, pairwise_rows, s0_rows


def _segments(rows: list[dict]) -> list[list[dict]]:
    segments: list[list[dict]] = []
    for row in rows:
        if not segments or float(row["snr_db"]) - float(segments[-1][-1]["snr_db"]) > 0.26:
            segments.append([])
        segments[-1].append(row)
    return segments


def plot(scene_candidates: dict[str, list[dict]]) -> dict[str, dict[str, object]]:
    markers = ("o", "s", "^", "v", "D", "P", "X", "h")
    colors = plt.get_cmap("tab10")(np.linspace(0.0, 0.9, 8))
    styles = {
        f"{index + 1:02d}": {
            "color": colors[index].tolist(),
            "marker": markers[index],
            "linestyle": "-",
            "linewidth": 2.5,
            "markersize": 8,
        }
        for index in range(8)
    }
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    for scene, candidates in sorted(scene_candidates.items()):
        fig, (bler_axis, nmse_axis) = plt.subplots(2, 1, figsize=(14.5, 10.5), sharex=True)
        for candidate in candidates:
            candidate_id = candidate["candidate_id"]
            suffix = candidate_id.rsplit("_", 1)[1]
            style = styles[suffix]
            label = f"candidate {suffix}" + (" (S0)" if suffix == "01" else "")
            first = True
            for segment in _segments(candidate["rows"]):
                x = np.asarray([float(row["snr_db"]) for row in segment])
                trials = np.asarray([int(row["trials"]) for row in segment])
                bler = np.asarray([float(row["bler"]) for row in segment])
                display = np.where(bler > 0.0, bler, 0.5 / trials)
                bler_axis.plot(x, display, label=label if first else None, **style)
                nmse_axis.plot(
                    x,
                    [float(row["ce_nmse_db"]) for row in segment],
                    label=label if first else None,
                    **style,
                )
                first = False
        bler_axis.set_yscale("log")
        bler_axis.axhline(0.10, color="black", linestyle=":", linewidth=1.8)
        bler_axis.axhline(0.01, color="black", linestyle="--", linewidth=1.8)
        bler_axis.set_ylabel("Estimated-CSI DCI BLER", fontsize=16)
        bler_axis.set_title(f"Plan-031 strict Sidon final fine scan: {scene}", fontsize=17)
        nmse_axis.set_xlabel("SNR (dB)", fontsize=16)
        nmse_axis.set_ylabel("Data-RE CE NMSE (dB)", fontsize=16)
        for axis in (bler_axis, nmse_axis):
            axis.grid(True, which="both", alpha=0.3)
            axis.tick_params(labelsize=14)
        handles, labels = bler_axis.get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncol=2, fontsize=14)
        fig.tight_layout(rect=(0.0, 0.08, 1.0, 1.0))
        fig.savefig(FIGURE_ROOT / f"strict_sidon_{scene}_bler_nmse.png", dpi=180)
        plt.close(fig)
    return styles


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bootstrap-repeats", type=int, default=10_000)
    args = parser.parse_args()
    repeats = int(args.bootstrap_repeats)
    configs, candidates, points, diagnostics = collect()
    targets, samples, brackets = analyze_targets(candidates, repeats)
    winners, pairwise, gains = rank_candidates(candidates, targets, samples, repeats)
    styles = plot(candidates)
    nmse_rows = [row for scene in candidates.values() for candidate in scene for row in candidate["nmse"]]

    ANALYSIS_ROOT.mkdir(parents=True, exist_ok=True)
    _write_csv(ANALYSIS_ROOT / "formal_points.csv", points)
    _write_csv(ANALYSIS_ROOT / "nmse_trial_summary.csv", nmse_rows)
    _write_csv(ANALYSIS_ROOT / "diagnostics.csv", diagnostics)
    _write_csv(ANALYSIS_ROOT / "target_snr.csv", targets)
    _write_csv(ANALYSIS_ROOT / "target_bracket_points.csv", brackets)
    _write_csv(ANALYSIS_ROOT / "winner_summary.csv", winners)
    _write_csv(ANALYSIS_ROOT / "pairwise_1pct.csv", pairwise)
    _write_csv(ANALYSIS_ROOT / "gain_vs_s0.csv", gains)
    (ANALYSIS_ROOT / "curve_styles.json").write_text(
        json.dumps(styles, indent=2) + "\n", encoding="utf-8"
    )
    git_head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True
    ).stdout.strip()
    metadata = {
        "schema": "plan031-strict-sidon-fine-analysis-v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "plan": "research/plan-031-PDCCH-CDD时延-BLER.md",
        "git_head": git_head,
        "working_tree_note": "plan-031 strict-Sidon implementation and result update are uncommitted",
        "bootstrap_repeats": repeats,
        "candidate_count": 18,
        "point_count": len(points),
        "total_trials": sum(int(row["trials"]) for row in points),
        "total_errors": sum(int(row["errors"]) for row in points),
        "estimated_target_count": sum(row["status"] == "estimated" for row in targets),
        "unqualified_target_count": sum(row["status"] != "estimated" for row in targets),
        "maximum_raw_bracket_width_db": max(
            float(row["raw_bracket_width_db"])
            for row in targets
            if row["status"] == "estimated"
        ),
        "config_sha256": {
            path.relative_to(ROOT).as_posix(): _sha256(path) for path in CONFIGS.values()
        },
        "supplement_config_sha256": {
            path.relative_to(ROOT).as_posix(): _sha256(path)
            for path in SUPPLEMENT_CONFIGS.values()
        },
        "fine_config_receipt_sha256": _sha256(RUN_ROOT / "fine_config_receipt.json"),
        "fine_supplement_config_receipt_sha256": _sha256(
            RUN_ROOT / "fine_supplement_config_receipt.json"
        ),
        "source_sha256": {
            "tools/analyze_plan031_strict_sidon_fine.py": _sha256(Path(__file__)),
            "tools/run_plan031_strict_sidon_prescan.py": _sha256(
                ROOT / "tools/run_plan031_strict_sidon_prescan.py"
            ),
            "cdd_lls/sim/pdcch_cdd.py": _sha256(ROOT / "cdd_lls/sim/pdcch_cdd.py"),
        },
        "winner_summary": winners,
        "reproduce_command": "python tools/analyze_plan031_strict_sidon_fine.py --bootstrap-repeats 10000",
    }
    (ANALYSIS_ROOT / "analysis_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
