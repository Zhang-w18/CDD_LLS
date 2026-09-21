"""Analyze completed AL1/AL2 shards of the plan-031 C300 4Tx/2Rx run."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
from typing import Any

import numpy as np
import yaml

import analyze_plan031 as common
import analyze_plan031_c300 as c300


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "outputs"
    / "experiment031_pdcch_cdd"
    / "20260911_c300_4tx_2rx_2sym"
)
ANALYSIS = OUTPUT / "analysis" / "partial_al1_al2"
FIGURES = ROOT / "docs" / "figures" / "result-031" / "c300_4tx_2rx_2sym"
ALS = (1, 2)
MODES = ("estimated", "ideal")
TARGETS = (0.10, 0.01)
BASELINES = ("C300_B0_QC", "C300_PRG_DFT4_6REG")
HISTORICAL = {
    "estimated": (
        ROOT
        / "outputs"
        / "experiment031_pdcch_cdd"
        / "20260907_c300_4tx_2sym"
        / "analysis"
    ),
    "ideal": (
        ROOT
        / "outputs"
        / "experiment031_pdcch_cdd"
        / "20260910_c300_4tx_2sym_ideal_csi"
        / "analysis"
    ),
}


def _config_path(al: int, mode: str) -> Path:
    return (
        ROOT
        / "configs"
        / f"pdcch_result031_c300_2rx_{mode}_al{al}_formal.yaml"
    )


def _read_csv(path: Path) -> list[dict[str, str]]:
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _flag_stem(candidate_id: str, snr_db: float) -> str:
    safe = "".join(char if char.isalnum() or char in "-_" else "_" for char in candidate_id)
    tag = str(float(snr_db)).replace("-", "m").replace(".", "p")
    return f"{safe}_snr_{tag}"


def collect() -> tuple[dict[tuple[int, str], dict], list[dict[str, Any]], list[dict[str, Any]]]:
    configs: dict[tuple[int, str], dict] = {}
    points: list[dict[str, Any]] = []
    diagnostics: list[dict[str, Any]] = []
    for al in ALS:
        for mode in MODES:
            config_path = _config_path(al, mode)
            with open(config_path, "r", encoding="utf-8") as handle:
                config = yaml.safe_load(handle) or {}
            configs[(al, mode)] = config
            if int(config["antenna"]["n_tx"]) != 4 or int(config["antenna"]["n_rx"]) != 2:
                raise ValueError(f"AL{al} {mode}: expected 4Tx/2Rx")
            expected_ce = "ideal" if mode == "ideal" else "frequency_lmmse"
            if str(config["receiver"]["channel_estimation"]).lower() != expected_ce:
                raise ValueError(f"AL{al} {mode}: unexpected channel estimation mode")
            default_snr = [float(value) for value in config["simulation"]["snr_points_db"]]
            for candidate in config["candidates"]:
                candidate_id = str(candidate["candidate_id"])
                expected_snr = sorted(
                    float(value) for value in candidate.get("snr_points_db", default_snr)
                )
                candidate_dir = Path(config["output_dir"]) / candidate_id
                csv_path = candidate_dir / "bler_points.csv"
                metadata_path = candidate_dir / "run_metadata.json"
                if not csv_path.exists() or not metadata_path.exists():
                    raise FileNotFoundError(f"Incomplete formal candidate: {candidate_dir}")
                rows = sorted(_read_csv(csv_path), key=lambda row: float(row["snr_db"]))
                actual_snr = [float(row["snr_db"]) for row in rows]
                if actual_snr != expected_snr or len(actual_snr) != len(set(actual_snr)):
                    raise ValueError(f"AL{al} {mode} {candidate_id}: incomplete SNR grid")
                metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
                meta_candidate = metadata["candidates"][0]
                previous_bler = None
                increases = 0
                max_flag_error = 0
                max_nmse_error_db = 0.0
                max_ideal_nmse = 0.0
                for row in rows:
                    if int(row["n_tx"]) != 4 or int(row["n_rx"]) != 2:
                        raise ValueError(f"AL{al} {mode} {candidate_id}: CSV antenna mismatch")
                    if str(row["csi_mode"]).lower() != expected_ce:
                        raise ValueError(f"AL{al} {mode} {candidate_id}: CSV CSI mismatch")
                    normalized: dict[str, Any] = dict(row)
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
                        "n_tx",
                        "n_rx",
                    ):
                        normalized[key] = int(float(row[key]))
                    normalized["csi_mode_label"] = mode
                    normalized["source_csv"] = csv_path.relative_to(ROOT).as_posix()
                    points.append(normalized)
                    bler = float(row["bler"])
                    if previous_bler is not None and bler > previous_bler + 1e-15:
                        increases += 1
                    previous_bler = bler
                    stem = _flag_stem(candidate_id, float(row["snr_db"]))
                    flags = np.load(candidate_dir / "trial_error_flags" / f"{stem}_error_flags.npy")
                    nmse = np.load(candidate_dir / "trial_error_flags" / f"{stem}_ce_nmse.npy")
                    max_flag_error = max(
                        max_flag_error,
                        abs(len(flags) - int(row["trials"])),
                        abs(int(np.sum(flags)) - int(row["errors"])),
                    )
                    recomputed = 10.0 * math.log10(max(float(np.mean(nmse)), 1e-30))
                    max_nmse_error_db = max(
                        max_nmse_error_db, abs(recomputed - float(row["ce_nmse_db"]))
                    )
                    if mode == "ideal":
                        max_ideal_nmse = max(max_ideal_nmse, float(np.max(np.abs(nmse))))
                diagnostics.append(
                    {
                        "aggregation_level": al,
                        "csi_mode": mode,
                        "candidate_id": candidate_id,
                        "receiver_covariance_mode": meta_candidate["receiver_covariance_mode"],
                        "pilot_rank": int(meta_candidate["pilot_rank"]),
                        "pilot_condition_number": float(meta_candidate["pilot_condition_number"]),
                        "ce_floor_nmse_db": 10.0
                        * math.log10(max(float(meta_candidate["ce_floor_nmse"]), 1e-30)),
                        "formal_points": len(rows),
                        "total_trials": sum(int(row["trials"]) for row in rows),
                        "total_errors": sum(int(row["errors"]) for row in rows),
                        "adjacent_bler_increases": increases,
                        "max_flag_count_or_sum_error": max_flag_error,
                        "max_ce_recompute_error_db": max_nmse_error_db,
                        "max_ideal_ce_nmse_linear": max_ideal_nmse,
                        "reached_max_without_200_errors": sum(
                            int(row["trials"])
                            == int(config["simulation"]["max_trials_per_snr"])
                            and int(row["errors"])
                            < int(config["simulation"]["target_errors"])
                            for row in rows
                        ),
                    }
                )
    return configs, points, diagnostics


def analyze_targets(
    configs: dict[tuple[int, str], dict],
    points: list[dict[str, Any]],
    repeats: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[tuple[str, int, str, float], np.ndarray]]:
    targets: list[dict[str, Any]] = []
    replicates: dict[tuple[str, int, str, float], np.ndarray] = {}
    for (al, mode), config in configs.items():
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            rows = sorted(
                [
                    row
                    for row in points
                    if int(row["aggregation_level"]) == al
                    and str(row["csi_mode_label"]) == mode
                    and str(row["candidate_id"]) == candidate_id
                ],
                key=lambda row: float(row["snr_db"]),
            )
            for target in TARGETS:
                key = (mode, al, candidate_id, target)
                try:
                    bracket_lo, bracket_hi = c300.qualified_raw_bracket(rows, target)
                    estimate, ci_lo, ci_hi, samples = common.bootstrap_targets(
                        rows,
                        target,
                        repeats,
                        common.stable_seed(
                            20260911, mode, al, candidate_id, target, "2rx-bootstrap"
                        ),
                    )
                    status = "estimated"
                except ValueError as error:
                    bracket_lo = bracket_hi = estimate = ci_lo = ci_hi = None
                    samples = np.asarray([], dtype=np.float64)
                    status = f"unqualified: {error}"
                replicates[key] = samples
                targets.append(
                    {
                        "aggregation_level": al,
                        "csi_mode": mode,
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
                        "method": "Jeffreys-smoothed trial-weighted decreasing isotonic fit and local log10-BLER interpolation",
                        "ci_method": "independent per-point Bernoulli bootstrap",
                        "bootstrap_repeats_requested": repeats,
                        "bootstrap_repeats_valid": len(samples),
                    }
                )
    lookup = {
        (str(row["csi_mode"]), int(row["aggregation_level"]), str(row["candidate_id"]), float(row["target_bler"])): row
        for row in targets
    }
    gains: list[dict[str, Any]] = []
    for (al, mode), config in configs.items():
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            for baseline in BASELINES:
                for target in TARGETS:
                    key = (mode, al, candidate_id, target)
                    base_key = (mode, al, baseline, target)
                    count = min(len(replicates[key]), len(replicates[base_key]))
                    if count:
                        samples = replicates[base_key][:count] - replicates[key][:count]
                        low, high = (float(value) for value in np.quantile(samples, [0.025, 0.975]))
                        gain = float(lookup[base_key]["target_snr_db"]) - float(
                            lookup[key]["target_snr_db"]
                        )
                        status = "estimated"
                    else:
                        gain = low = high = None
                        status = "unavailable_unqualified_target"
                    gains.append(
                        {
                            "aggregation_level": al,
                            "csi_mode": mode,
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
    return targets, gains, replicates


def delta_ce(
    configs: dict[tuple[int, str], dict],
    targets: list[dict[str, Any]],
    replicates: dict[tuple[str, int, str, float], np.ndarray],
) -> list[dict[str, Any]]:
    lookup = {
        (str(row["csi_mode"]), int(row["aggregation_level"]), str(row["candidate_id"]), float(row["target_bler"])): row
        for row in targets
    }
    rows_out: list[dict[str, Any]] = []
    for al in ALS:
        for candidate in configs[(al, "estimated")]["candidates"]:
            estimated_id = str(candidate["candidate_id"])
            ideal_id = (
                "C300_SMALL_CDD_QSTEP0P25_IDEAL_CSI"
                if estimated_id.startswith("C300_SMALL_CDD_QSTEP0P25_")
                else estimated_id
            )
            for target in TARGETS:
                estimated_key = ("estimated", al, estimated_id, target)
                ideal_key = ("ideal", al, ideal_id, target)
                count = min(len(replicates[estimated_key]), len(replicates[ideal_key]))
                if count:
                    samples = replicates[estimated_key][:count] - replicates[ideal_key][:count]
                    low, high = (float(value) for value in np.quantile(samples, [0.025, 0.975]))
                    delta = float(lookup[estimated_key]["target_snr_db"]) - float(
                        lookup[ideal_key]["target_snr_db"]
                    )
                    status = "estimated"
                else:
                    delta = low = high = None
                    status = "unavailable_unqualified_target"
                rows_out.append(
                    {
                        "aggregation_level": al,
                        "estimated_candidate_id": estimated_id,
                        "ideal_candidate_id": ideal_id,
                        "target_bler": target,
                        "delta_ce_db": delta,
                        "delta_ce_ci95_lo_db": low,
                        "delta_ce_ci95_hi_db": high,
                        "status": status,
                        "positive_means_estimated_requires_more_snr": True,
                        "bootstrap_repeats_valid": count,
                    }
                )
    return rows_out


def delta_rx(
    configs: dict[tuple[int, str], dict],
    targets: list[dict[str, Any]],
    replicates: dict[tuple[str, int, str, float], np.ndarray],
    repeats: int,
) -> list[dict[str, Any]]:
    two_lookup = {
        (str(row["csi_mode"]), int(row["aggregation_level"]), str(row["candidate_id"]), float(row["target_bler"])): row
        for row in targets
    }
    historical_targets = {
        mode: _read_csv(path / "target_snr.csv") for mode, path in HISTORICAL.items()
    }
    historical_points = {
        mode: _read_csv(path / "formal_points.csv") for mode, path in HISTORICAL.items()
    }
    historical_lookup = {
        mode: {
            (int(row["aggregation_level"]), str(row["candidate_id"]), float(row["target_bler"])): row
            for row in rows
        }
        for mode, rows in historical_targets.items()
    }
    rows_out: list[dict[str, Any]] = []
    for (al, mode), config in configs.items():
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            for target in TARGETS:
                two_key = (mode, al, candidate_id, target)
                one_row = historical_lookup[mode][(al, candidate_id, target)]
                two_samples = replicates[two_key]
                if str(one_row["status"]) != "estimated" or len(two_samples) == 0:
                    delta = low = high = None
                    count = 0
                    status = "unavailable_unqualified_target"
                else:
                    one_points = sorted(
                        [
                            row
                            for row in historical_points[mode]
                            if int(row["aggregation_level"]) == al
                            and str(row["candidate_id"]) == candidate_id
                        ],
                        key=lambda row: float(row["snr_db"]),
                    )
                    _, _, _, one_samples = common.bootstrap_targets(
                        one_points,
                        target,
                        repeats,
                        common.stable_seed(
                            20260913, mode, al, candidate_id, target, "1rx-for-delta-rx"
                        ),
                    )
                    count = min(len(one_samples), len(two_samples))
                    samples = one_samples[:count] - two_samples[:count]
                    low, high = (float(value) for value in np.quantile(samples, [0.025, 0.975]))
                    delta = float(one_row["target_snr_db"]) - float(
                        two_lookup[two_key]["target_snr_db"]
                    )
                    status = "estimated"
                rows_out.append(
                    {
                        "aggregation_level": al,
                        "csi_mode": mode,
                        "candidate_id": candidate_id,
                        "target_bler": target,
                        "delta_rx_db": delta,
                        "delta_rx_ci95_lo_db": low,
                        "delta_rx_ci95_hi_db": high,
                        "status": status,
                        "positive_means_2rx_requires_less_snr": True,
                        "bootstrap_repeats_valid": count,
                    }
                )
    return rows_out


def plot(configs: dict[tuple[int, str], dict], points: list[dict[str, Any]]) -> None:
    import matplotlib.pyplot as plt

    FIGURES.mkdir(parents=True, exist_ok=True)
    style_map: dict[str, dict[str, Any]] = {}
    for config in configs.values():
        for candidate in config["candidates"]:
            style_map[str(candidate["candidate_id"])] = dict(candidate.get("style", {}))
    (ANALYSIS / "style_map.json").write_text(
        json.dumps(style_map, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    for al in ALS:
        for mode in MODES:
            config = configs[(al, mode)]
            fig, axis = plt.subplots(figsize=(15.0, 8.5))
            for candidate in config["candidates"]:
                candidate_id = str(candidate["candidate_id"])
                rows = sorted(
                    [
                        row
                        for row in points
                        if int(row["aggregation_level"]) == al
                        and str(row["csi_mode_label"]) == mode
                        and str(row["candidate_id"]) == candidate_id
                    ],
                    key=lambda row: float(row["snr_db"]),
                )
                style = style_map[candidate_id]
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
            axis.set_xlabel("SNR (dB)", fontsize=16)
            axis.set_ylabel(f"{mode.capitalize()}-CSI DCI BLER", fontsize=16)
            axis.set_title(f"C300 4Tx/2Rx two-symbol PDCCH AL{al}: {mode} CSI", fontsize=17)
            axis.tick_params(labelsize=14)
            axis.grid(True, which="both", alpha=0.35)
            axis.legend(fontsize=14, loc="center left", bbox_to_anchor=(1.01, 0.5))
            fig.tight_layout()
            fig.savefig(FIGURES / f"al{al}_{mode}_csi_bler.png", dpi=180)
            plt.close(fig)
        config = configs[(al, "estimated")]
        fig, axis = plt.subplots(figsize=(15.0, 8.5))
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            rows = sorted(
                [
                    row
                    for row in points
                    if int(row["aggregation_level"]) == al
                    and str(row["csi_mode_label"]) == "estimated"
                    and str(row["candidate_id"]) == candidate_id
                ],
                key=lambda row: float(row["snr_db"]),
            )
            style = style_map[candidate_id]
            axis.plot(
                [float(row["snr_db"]) for row in rows],
                [float(row["ce_nmse_db"]) for row in rows],
                label=str(candidate.get("label", candidate_id)),
                color=style.get("color"),
                linestyle=style.get("linestyle", "-"),
                marker=style.get("marker", "o"),
                linewidth=2.5,
                markersize=8,
            )
        axis.set_xlabel("SNR (dB)", fontsize=16)
        axis.set_ylabel("Estimated-CSI CE NMSE (dB)", fontsize=16)
        axis.set_title(f"C300 4Tx/2Rx two-symbol PDCCH AL{al}: CE NMSE", fontsize=17)
        axis.tick_params(labelsize=14)
        axis.grid(True, alpha=0.35)
        axis.legend(fontsize=14, loc="center left", bbox_to_anchor=(1.01, 0.5))
        fig.tight_layout()
        fig.savefig(FIGURES / f"al{al}_estimated_csi_nmse.png", dpi=180)
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bootstrap-repeats", type=int, default=2000)
    args = parser.parse_args()
    repeats = int(args.bootstrap_repeats)
    configs, points, diagnostics = collect()
    targets, gains, replicates = analyze_targets(configs, points, repeats)
    ce_rows = delta_ce(configs, targets, replicates)
    rx_rows = delta_rx(configs, targets, replicates, repeats)
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    common.write_csv(ANALYSIS / "formal_points.csv", points)
    common.write_csv(ANALYSIS / "diagnostics.csv", diagnostics)
    common.write_csv(ANALYSIS / "target_snr.csv", targets)
    common.write_csv(ANALYSIS / "target_gains.csv", gains)
    common.write_csv(ANALYSIS / "delta_ce.csv", ce_rows)
    common.write_csv(ANALYSIS / "delta_rx.csv", rx_rows)
    plot(configs, points)
    git_head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True
    ).stdout.strip()
    metadata = {
        "schema": "plan031-c300-2rx-partial-al1-al2-analysis-v1",
        "scope": "AL1 and AL2 only; AL4 formal run still in progress",
        "bootstrap_repeats": repeats,
        "git_head": git_head,
        "config_sha256": {
            path.relative_to(ROOT).as_posix(): common.sha256(path)
            for path in (_config_path(al, mode) for al in ALS for mode in MODES)
        },
        "formal_points": len(points),
        "total_trials": sum(int(row["trials"]) for row in points),
        "total_errors": sum(int(row["errors"]) for row in points),
        "target_rows": len(targets),
        "qualified_targets": sum(row["status"] == "estimated" for row in targets),
        "delta_ce_rows": len(ce_rows),
        "delta_rx_rows": len(rx_rows),
        "qualified_delta_rx": sum(row["status"] == "estimated" for row in rx_rows),
        "max_flag_count_or_sum_error": max(
            int(row["max_flag_count_or_sum_error"]) for row in diagnostics
        ),
        "max_ce_recompute_error_db": max(
            float(row["max_ce_recompute_error_db"]) for row in diagnostics
        ),
        "max_ideal_ce_nmse_linear": max(
            float(row["max_ideal_ce_nmse_linear"]) for row in diagnostics
        ),
    }
    (ANALYSIS / "analysis_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
