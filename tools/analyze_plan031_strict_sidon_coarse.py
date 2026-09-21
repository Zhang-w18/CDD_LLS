"""Validate, plot, and rank the plan-031 3000-trial coarse-confirmation run."""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml

from analyze_plan031 import interpolate_target, isotonic_decreasing, stable_seed


ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search"
DATA_ROOT = RUN_ROOT / "coarse_confirmation"
OUTPUT_ROOT = RUN_ROOT / "coarse_analysis"
CONFIG_PATHS = sorted(ROOT.glob("configs/pdcch_result031_strict_sidon_*_coarse.yaml"))
TARGETS = (0.10, 0.01)
BOOTSTRAP_REPEATS = 4000


def _snr_token(value: float) -> str:
    return f"{value:.1f}".replace("-", "m").replace(".", "p")


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _target(
    rows: list[dict[str, str]], target: float, seed: int
) -> tuple[dict[str, object], np.ndarray]:
    ordered = sorted(rows, key=lambda row: float(row["snr_db"]))
    snr = np.asarray([float(row["snr_db"]) for row in ordered])
    trials = np.asarray([int(row["trials"]) for row in ordered], dtype=np.int64)
    errors = np.asarray([int(row["errors"]) for row in ordered], dtype=np.int64)
    raw = errors / trials
    brackets = [
        (float(left["snr_db"]), float(right["snr_db"]))
        for left, right in zip(ordered[:-1], ordered[1:])
        if float(left["bler"]) >= target >= float(right["bler"])
    ]
    if not brackets:
        if np.all(raw > target):
            status = "censored_above_target_at_max_snr"
            bound = f"> {snr[-1]:g} dB"
        elif np.all(raw < target):
            status = "censored_below_target_at_min_snr"
            bound = f"< {snr[0]:g} dB"
        else:
            status = "unbracketed_nonmonotonic_or_disjoint_grid"
            bound = None
        return {
            "target_bler": target,
            "status": status,
            "target_snr_db": None,
            "ci95_lo_db": None,
            "ci95_hi_db": None,
            "bound": bound,
            "raw_bracket_lo_db": None,
            "raw_bracket_hi_db": None,
            "bootstrap_valid": 0,
        }, np.asarray([], dtype=np.float64)

    probabilities = (errors + 0.5) / (trials + 1.0)
    fitted = isotonic_decreasing(probabilities, trials.astype(np.float64))
    estimate, index = interpolate_target(snr, fitted, target)
    rng = np.random.default_rng(seed)
    samples = []
    for _ in range(BOOTSTRAP_REPEATS):
        sampled_errors = rng.binomial(trials, raw)
        sampled = (sampled_errors + 0.5) / (trials + 1.0)
        sampled_fit = isotonic_decreasing(sampled, trials.astype(np.float64))
        try:
            value, _ = interpolate_target(snr, sampled_fit, target)
        except ValueError:
            continue
        samples.append(value)
    array = np.asarray(samples, dtype=np.float64)
    if len(array) < BOOTSTRAP_REPEATS // 2:
        raise RuntimeError(
            f"Only {len(array)}/{BOOTSTRAP_REPEATS} target bootstrap samples were bracketed"
        )
    low, high = np.quantile(array, [0.025, 0.975])
    return {
        "target_bler": target,
        "status": "estimated",
        "target_snr_db": estimate,
        "ci95_lo_db": float(low),
        "ci95_hi_db": float(high),
        "bound": None,
        "raw_bracket_lo_db": float(snr[index]),
        "raw_bracket_hi_db": float(snr[index + 1]),
        "bootstrap_valid": len(array),
    }, array


def _load() -> tuple[dict[str, list[dict]], list[dict[str, object]]]:
    scenes: dict[str, list[dict]] = {}
    diagnostics = []
    for config_path in CONFIG_PATHS:
        config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
        scene = Path(config["output_dir"]).name
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            candidate_dir = DATA_ROOT / scene / candidate_id
            csv_path = candidate_dir / "bler_points.csv"
            metadata_path = candidate_dir / "run_metadata.json"
            if not csv_path.exists() or not metadata_path.exists():
                raise FileNotFoundError(f"Incomplete candidate output: {candidate_dir}")
            rows = sorted(_read_csv(csv_path), key=lambda row: float(row["snr_db"]))
            expected = sorted(float(value) for value in candidate["snr_points_db"])
            actual = [float(row["snr_db"]) for row in rows]
            if actual != expected or any(int(row["trials"]) != 3000 for row in rows):
                raise RuntimeError(f"Grid/trial mismatch for {candidate_id}")
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            meta_candidate = metadata["candidates"][0]
            max_error_count_delta = 0
            max_ce_mean_delta_db = 0.0
            for row in rows:
                snr = float(row["snr_db"])
                stem = f"{candidate_id}_snr_{_snr_token(snr)}"
                trial_dir = candidate_dir / "trial_error_flags"
                error_flags = np.load(trial_dir / f"{stem}_error_flags.npy")
                ce_nmse = np.load(trial_dir / f"{stem}_ce_nmse.npy")
                trials = int(row["trials"])
                if error_flags.shape != (trials,) or ce_nmse.shape != (trials,):
                    raise RuntimeError(f"Per-trial array length mismatch for {stem}")
                error_delta = abs(int(np.count_nonzero(error_flags)) - int(row["errors"]))
                ce_mean_db = 10.0 * math.log10(max(float(np.mean(ce_nmse)), 1e-30))
                ce_delta_db = abs(ce_mean_db - float(row["ce_nmse_db"]))
                max_error_count_delta = max(max_error_count_delta, error_delta)
                max_ce_mean_delta_db = max(max_ce_mean_delta_db, ce_delta_db)
            if max_error_count_delta != 0 or max_ce_mean_delta_db > 1e-10:
                raise RuntimeError(f"Per-trial evidence mismatch for {candidate_id}")
            scenes.setdefault(scene, []).append(
                {"candidate_id": candidate_id, "rows": rows}
            )
            diagnostics.append(
                {
                    "scene": scene,
                    "candidate_id": candidate_id,
                    "pilot_rank": int(meta_candidate["pilot_rank"]),
                    "pilot_condition_number": float(meta_candidate["pilot_condition_number"]),
                    "ce_floor_nmse_db": 10.0
                    * math.log10(max(float(meta_candidate["ce_floor_nmse"]), 1e-30)),
                    "points": len(rows),
                    "trials": sum(int(row["trials"]) for row in rows),
                    "max_error_count_delta": max_error_count_delta,
                    "max_ce_mean_delta_db": max_ce_mean_delta_db,
                }
            )
    if sum(len(values) for values in scenes.values()) != 43:
        raise RuntimeError("Expected 43 coarse-confirmation candidates")
    return scenes, diagnostics


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _plot(scene: str, candidates: list[dict], styles: dict[str, dict[str, object]]) -> None:
    fig, (bler_axis, nmse_axis) = plt.subplots(2, 1, figsize=(14.5, 11.0), sharex=True)
    for candidate in candidates:
        candidate_id = candidate["candidate_id"]
        rows = candidate["rows"]
        snr = np.asarray([float(row["snr_db"]) for row in rows])
        trials = np.asarray([int(row["trials"]) for row in rows])
        bler = np.asarray([float(row["bler"]) for row in rows])
        display = np.where(bler > 0.0, bler, 0.5 / trials)
        style = styles[candidate_id]
        bler_axis.semilogy(snr, display, label=candidate_id, **style)
        nmse_axis.plot(
            snr,
            [float(row["ce_nmse_db"]) for row in rows],
            label=candidate_id,
            **style,
        )
    bler_axis.axhline(0.10, color="black", linestyle=":", linewidth=1.8)
    bler_axis.axhline(0.01, color="black", linestyle="--", linewidth=1.8)
    bler_axis.set_ylabel("Estimated-CSI DCI BLER", fontsize=16)
    bler_axis.set_title(f"Plan-031 strict Sidon coarse confirmation: {scene}", fontsize=17)
    nmse_axis.set_xlabel("SNR (dB)", fontsize=16)
    nmse_axis.set_ylabel("Data-RE CE NMSE (dB)", fontsize=16)
    for axis in (bler_axis, nmse_axis):
        axis.grid(True, which="both", alpha=0.3)
        axis.tick_params(labelsize=14)
    handles, labels = bler_axis.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, fontsize=11)
    fig.tight_layout(rect=(0.0, 0.11, 1.0, 1.0))
    fig.savefig(OUTPUT_ROOT / f"{scene}_coarse_bler_nmse.png", dpi=180)
    plt.close(fig)


def main() -> None:
    scenes, diagnostics = _load()
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    target_rows = []
    sample_map: dict[tuple[str, str, float], np.ndarray] = {}
    for scene, candidates in sorted(scenes.items()):
        for candidate in candidates:
            candidate_id = candidate["candidate_id"]
            for target in TARGETS:
                result, samples = _target(
                    candidate["rows"],
                    target,
                    stable_seed(20260908, scene, candidate_id, target, "coarse-confirmation"),
                )
                target_rows.append(
                    {"scene": scene, "candidate_id": candidate_id, **result}
                )
                sample_map[(scene, candidate_id, target)] = samples

    preliminary: dict[str, list[str]] = {}
    rank_rows = []
    target_lookup = {
        (str(row["scene"]), str(row["candidate_id"]), float(row["target_bler"])): row
        for row in target_rows
    }
    for scene, candidates in sorted(scenes.items()):
        ids = [str(candidate["candidate_id"]) for candidate in candidates]
        estimated = [
            candidate_id
            for candidate_id in ids
            if target_lookup[(scene, candidate_id, 0.01)]["status"] == "estimated"
        ]
        ordered = sorted(
            estimated,
            key=lambda candidate_id: float(
                target_lookup[(scene, candidate_id, 0.01)]["target_snr_db"]
            ),
        )
        leader = ordered[0]
        leader_high = float(target_lookup[(scene, leader, 0.01)]["ci95_hi_db"])
        common = min(len(sample_map[(scene, candidate_id, 0.01)]) for candidate_id in ordered)
        matrix = np.vstack(
            [sample_map[(scene, candidate_id, 0.01)][:common] for candidate_id in ordered]
        )
        winners = np.argmin(matrix, axis=0)
        probabilities = {
            candidate_id: float(np.mean(winners == index))
            for index, candidate_id in enumerate(ordered)
        }
        plausible = [
            candidate_id
            for candidate_id in ordered
            if probabilities[candidate_id] >= 0.05
            or float(target_lookup[(scene, candidate_id, 0.01)]["ci95_lo_db"])
            <= leader_high
        ]
        for candidate_id in ids:
            row = target_lookup[(scene, candidate_id, 0.01)]
            if row["status"] == "censored_below_target_at_min_snr":
                plausible.append(candidate_id)
            elif (
                row["status"] == "censored_above_target_at_max_snr"
                and float(str(row["bound"]).split()[1]) <= leader_high
            ):
                plausible.append(candidate_id)
        plausible = sorted(set(plausible))
        preliminary[scene] = plausible
        for rank, candidate_id in enumerate(ordered, start=1):
            row = target_lookup[(scene, candidate_id, 0.01)]
            rank_rows.append(
                {
                    "scene": scene,
                    "candidate_id": candidate_id,
                    "rank_by_1pct_point_estimate": rank,
                    "snr_at_1pct_db": row["target_snr_db"],
                    "ci95_lo_db": row["ci95_lo_db"],
                    "ci95_hi_db": row["ci95_hi_db"],
                    "bootstrap_probability_best": probabilities[candidate_id],
                    "preliminary_plausible_best": candidate_id in plausible,
                }
            )

    markers = ("o", "s", "^", "v", "D", "P", "X", "h")
    colors = plt.get_cmap("tab10")(np.linspace(0.0, 0.9, 8))
    styles = {}
    for candidates in scenes.values():
        for candidate in candidates:
            candidate_id = candidate["candidate_id"]
            index = int(candidate_id.rsplit("_", 1)[1]) - 1
            styles[candidate_id] = {
                "color": colors[index].tolist(),
                "linestyle": "-",
                "marker": markers[index],
                "linewidth": 2.5,
                "markersize": 8,
            }
    for scene, candidates in sorted(scenes.items()):
        _plot(scene, candidates, styles)

    _write_csv(OUTPUT_ROOT / "target_snr_coarse.csv", target_rows)
    _write_csv(OUTPUT_ROOT / "candidate_ranking_coarse.csv", rank_rows)
    _write_csv(OUTPUT_ROOT / "diagnostics.csv", diagnostics)
    (OUTPUT_ROOT / "curve_styles.json").write_text(
        json.dumps(styles, indent=2) + "\n", encoding="utf-8"
    )
    receipt = {
        "schema": "plan031-strict-sidon-coarse-analysis-v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "candidate_count": 43,
        "point_count": 254,
        "candidate_trials": 762000,
        "bootstrap_repeats": BOOTSTRAP_REPEATS,
        "target_rows": len(target_rows),
        "unbracketed_target_rows": sum(row["status"] != "estimated" for row in target_rows),
        "preliminary_rule": "1% bootstrap P(best)>=5% or 1% CI overlaps the point-estimate leader CI; censored potentially-better candidates are retained",
        "preliminary_candidates_by_scene": preliminary,
        "status": "exploratory coarse ranking; not a final winner declaration",
        "reproduce_command": "python tools/analyze_plan031_strict_sidon_coarse.py",
    }
    (OUTPUT_ROOT / "analysis_receipt.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
