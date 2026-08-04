"""Plan-027 E3 matched-CE diagnostics with conditional DMRS density."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.design import FrequencyCEMetric
from cdd_lls.phy.resource_grid import local_indices_for_subcarriers
from tools.run_plan025_delay_matched_tdl import write_csv_rows
from tools.run_plan027_meff_design import (
    DEFAULT_OUTPUT_ROOT,
    DEFAULT_RUN_ID,
    N_TX,
    SCENARIOS,
    decode,
    physical_covariance,
    read_csv,
    resource_and_grid,
    save_json,
)


TRIGGER_DB = 0.5


def common_data_subcarriers() -> np.ndarray:
    common = None
    for spacing in (24, 12, 6):
        _, grid = resource_and_grid(spacing)
        coordinates = {tuple(value) for value in grid.data_coordinates.tolist()}
        common = coordinates if common is None else common & coordinates
    return np.asarray(
        [coordinate[1] for coordinate in sorted(common or set())],
        dtype=np.int64,
    )


def build_metrics(spread_ns: float) -> dict[int, FrequencyCEMetric]:
    _, covariance = physical_covariance(spread_ns)
    common = common_data_subcarriers()
    metrics = {}
    for spacing in (24, 12, 6):
        _, grid = resource_and_grid(spacing)
        pilots = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
        data = local_indices_for_subcarriers(grid, common)
        metrics[spacing] = FrequencyCEMetric(covariance, pilots, data, N_TX)
    return metrics


def ce_row(
    metric: FrequencyCEMetric,
    source: dict[str, str],
    coordinates: Sequence[float],
    spacing: int,
    target_probability: float,
    working_point: str,
    snr_db: float,
    reference_candidate_id: str,
) -> dict[str, object]:
    result = metric.evaluate_grid_coordinates(coordinates, snr_db)
    return {
        "scenario_id": source["scenario_id"],
        "candidate_id": source["candidate_id"],
        "family": source["family"],
        "target_outage_probability": target_probability,
        "working_point": working_point,
        "reference_candidate_id": reference_candidate_id,
        "dmrs_spacing_subcarriers": spacing,
        "pilot_count_per_symbol": 576 // spacing,
        "pilot_re_count": 2 * (576 // spacing),
        "common_data_re_count": 5568,
        "snr_db": snr_db,
        "ce_nmse": result["nmse"],
        "ce_nmse_db": result["nmse_db"],
        "pilot_rank": result["pilot_rank"],
        "pilot_condition_number": result["pilot_condition_number"],
        "pilot_minimum_singular_value": result["pilot_minimum_singular_value"],
        "lmmse_system_condition_number": result["condition_number"],
        "lmmse_system_minimum_singular_value": result["minimum_singular_value"],
        "numerical_jitter": result["numerical_jitter"],
        "r_pp_trace": result["r_pp_trace"],
        "r_pp_frobenius_norm": result["r_pp_frobenius_norm"],
        "r_dp_frobenius_norm": result["r_fp_frobenius_norm"],
    }


def target_snr(row: dict[str, str], probability: float) -> float:
    return float(row["outage10_snr_db" if probability == 0.10 else "outage1_snr_db"])


def baseline_id(row: dict[str, str], probability: float) -> str:
    return row[
        "ap_t1_base_10_candidate_id"
        if probability == 0.10
        else "ap_t1_base_1_candidate_id"
    ]


def is_positive(row: dict[str, str]) -> bool:
    return (
        float(row["gain_vs_ap_t1_base_10_db"]) > 0.0
        or float(row["gain_vs_ap_t1_base_1_db"]) > 0.0
    )


def run_ce(output: Path) -> None:
    source_rows = read_csv(output / "e1_outage" / "e1_outage_targets.csv")
    stage = output / "e3_ce_density"
    stage.mkdir(parents=True, exist_ok=True)
    ce_rows: list[dict[str, object]] = []
    trigger_rows: list[dict[str, object]] = []
    for scenario_id, spread_ns in SCENARIOS:
        scenario = [row for row in source_rows if row["scenario_id"] == scenario_id]
        by_id = {row["candidate_id"]: row for row in scenario}
        positives = [
            row
            for row in scenario
            if row["family"] in ("AP_T2_CTRL", "GEO_T1_CTRL", "MEFF_T2_CAND")
            and is_positive(row)
        ]
        diagnostics = {
            row["candidate_id"]: row
            for row in scenario
            if row["family"]
            in (
                "B0_QC",
                "S0_SIDON",
                "AP_RMS_T1",
                "AP_TEPS_T1",
                "AP_TU_NT",
                "AP_TALIAS_NT",
            )
        }
        diagnostics.update({row["candidate_id"]: row for row in positives})
        metrics = build_metrics(spread_ns)
        comb24_deltas: dict[str, list[float]] = {row["candidate_id"]: [] for row in positives}
        comb12_deltas: dict[str, list[float]] = {row["candidate_id"]: [] for row in positives}
        for probability in (0.10, 0.01):
            for row in diagnostics.values():
                coordinates = decode(row["delay_grid_coordinates"])
                own_snr = target_snr(row, probability)
                base = by_id[baseline_id(row, probability)]
                base_coordinates = decode(base["delay_grid_coordinates"])
                base_snr = target_snr(base, probability)
                ce_rows.append(
                    ce_row(
                        metrics[24],
                        row,
                        coordinates,
                        24,
                        probability,
                        "candidate_own_target",
                        own_snr,
                        base["candidate_id"],
                    )
                )
                candidate_at_base = ce_row(
                    metrics[24],
                    row,
                    coordinates,
                    24,
                    probability,
                    "common_baseline_target",
                    base_snr,
                    base["candidate_id"],
                )
                baseline_at_base = ce_row(
                    metrics[24],
                    base,
                    base_coordinates,
                    24,
                    probability,
                    "common_baseline_target",
                    base_snr,
                    base["candidate_id"],
                )
                ce_rows.extend([candidate_at_base, baseline_at_base])
                if row["candidate_id"] in comb24_deltas:
                    comb24_deltas[row["candidate_id"]].append(
                        float(candidate_at_base["ce_nmse_db"])
                        - float(baseline_at_base["ce_nmse_db"])
                    )
        for row in positives:
            candidate_id = row["candidate_id"]
            trigger12 = max(comb24_deltas[candidate_id]) > TRIGGER_DB
            trigger6 = False
            if trigger12:
                coordinates = decode(row["delay_grid_coordinates"])
                for probability in (0.10, 0.01):
                    base = by_id[baseline_id(row, probability)]
                    base_coordinates = decode(base["delay_grid_coordinates"])
                    base_snr = target_snr(base, probability)
                    candidate_at_base = ce_row(
                        metrics[12],
                        row,
                        coordinates,
                        12,
                        probability,
                        "common_baseline_target",
                        base_snr,
                        base["candidate_id"],
                    )
                    baseline_at_base = ce_row(
                        metrics[12],
                        base,
                        base_coordinates,
                        12,
                        probability,
                        "common_baseline_target",
                        base_snr,
                        base["candidate_id"],
                    )
                    ce_rows.extend([candidate_at_base, baseline_at_base])
                    comb12_deltas[candidate_id].append(
                        float(candidate_at_base["ce_nmse_db"])
                        - float(baseline_at_base["ce_nmse_db"])
                    )
                trigger6 = max(comb12_deltas[candidate_id]) > TRIGGER_DB
            if trigger6:
                coordinates = decode(row["delay_grid_coordinates"])
                for probability in (0.10, 0.01):
                    base = by_id[baseline_id(row, probability)]
                    base_coordinates = decode(base["delay_grid_coordinates"])
                    base_snr = target_snr(base, probability)
                    ce_rows.extend(
                        [
                            ce_row(
                                metrics[6],
                                row,
                                coordinates,
                                6,
                                probability,
                                "common_baseline_target",
                                base_snr,
                                base["candidate_id"],
                            ),
                            ce_row(
                                metrics[6],
                                base,
                                base_coordinates,
                                6,
                                probability,
                                "common_baseline_target",
                                base_snr,
                                base["candidate_id"],
                            ),
                        ]
                    )
            trigger_rows.append(
                {
                    "scenario_id": scenario_id,
                    "candidate_id": candidate_id,
                    "family": row["family"],
                    "comb24_max_ce_degradation_db": max(comb24_deltas[candidate_id]),
                    "trigger_threshold_db": TRIGGER_DB,
                    "comb12_triggered": trigger12,
                    "comb12_max_ce_degradation_db": (
                        max(comb12_deltas[candidate_id]) if comb12_deltas[candidate_id] else ""
                    ),
                    "comb6_triggered": trigger6,
                    "pilot_re_power_policy": "fixed_per_pilot_re",
                    "common_data_re_count": 5568,
                }
            )
        print(
            f"[CE] {scenario_id}: positives={len(positives)}, "
            f"comb12={sum(row['comb12_triggered'] for row in trigger_rows if row['scenario_id'] == scenario_id)}, "
            f"comb6={sum(row['comb6_triggered'] for row in trigger_rows if row['scenario_id'] == scenario_id)}",
            flush=True,
        )
    write_csv_rows(ce_rows, stage / "e3_ce_at_targets.csv")
    write_csv_rows(trigger_rows, stage / "e3_density_trigger.csv")
    save_json(
        {
            "trigger_threshold_db": TRIGGER_DB,
            "pilot_re_power_policy": "fixed_per_pilot_re",
            "pilot_spacings_subcarriers": [24, 12, 6],
            "pilot_counts_per_symbol": [24, 48, 96],
            "common_data_re_count": 5568,
            "interpretation": "CE recoverability with increased pilot energy and overhead; not net link gain",
        },
        stage / "resolved_experiment.json",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("ce",), required=True)
    parser.add_argument("--run-id", default=DEFAULT_RUN_ID)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    args = parser.parse_args()
    output = args.output_root / str(args.run_id)
    save_json(
        {
            "command": [sys.executable, *sys.argv],
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        output / "commands_ce.json",
    )
    run_ce(output)


if __name__ == "__main__":
    main()
