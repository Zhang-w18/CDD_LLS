"""Plan-027 dense-DMRS candidate redesign, outage, gate, and plots."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from pathlib import Path
from typing import Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.design import (
    FrequencyCEMetric,
    arithmetic_delay_sets,
    canonical_delay_set,
    continuous_arithmetic_delays_ns,
    delay_indices_to_ns,
    delay_ns_to_grid_coordinates,
    effective_moments_lag,
    thick_sidon_requirements,
    unique_random_delay_sets,
)
from cdd_lls.phy.resource_grid import local_indices_for_subcarriers
from tools.run_plan025_delay_matched_tdl import stable_seed, write_csv_rows
from tools.run_plan027_bler import curve_label, curve_style
from tools.run_plan027_meff_design import (
    DEFAULT_OUTPUT_ROOT,
    DEFAULT_RUN_ID,
    GEOMETRY_STATES,
    GRID_NS,
    K,
    MEFF_BANDS_DB,
    N_TX,
    QC,
    SEARCH_STATES,
    SCS_HZ,
    SEED,
    SIDON,
    add_equivalent_delay_metrics,
    batch_moments,
    continuous_candidate,
    decode,
    encode,
    geometry_selections,
    geometry_states,
    grid_candidate_row,
    iqam_table,
    outage_target,
    paired_bootstrap_gain,
    physical_covariance,
    read_csv,
    resource_and_grid,
    save_json,
    scenario_outage,
    select_meff_indices,
    support_for_scenario,
    write_large_candidate_csv,
)


SCENARIO_SPECS = {
    "A30": {
        "spread_ns": 30.0,
        "stage_root": "e5_dense_dmrs",
        "stage_label": "e5",
        "allowed_spacings": (12, 6),
    },
    "A100": {
        "spread_ns": 100.0,
        "stage_root": "e6_a100_dense_dmrs",
        "stage_label": "e6",
        "allowed_spacings": (6,),
    },
    "A300": {
        "spread_ns": 300.0,
        "stage_root": "a300_comb6",
        "stage_label": "p28",
        "allowed_spacings": (6,),
    },
}
SCENARIO_ID = "A30"
SPREAD_NS = 30.0
DENSE_STAGE_ROOT = "e5_dense_dmrs"
STAGE_LABEL = "e5"
DENSE_SPACINGS = (12, 6)
BASELINE_FAMILIES = (
    "B0_QC",
    "AP_RMS_T1",
    "AP_TEPS_T1",
    "AP_TU_NT",
    "AP_TALIAS_NT",
)
A100_ADDITIONAL_BASELINE_FAMILY = "AP_TU_NTM1"
SEARCH_FAMILIES = ("AP_T2_CTRL", "GEO_T1_CTRL", "MEFF_T2_CAND")
AP_TU_NT = (0, 72, 144, 216, 288, 360, 432, 504)
USEFUL_SYMBOL_NS = 1e9 / SCS_HZ
TARGETS = (("10pct", 0.10, "outage10_snr_db"), ("1pct", 0.01, "outage1_snr_db"))


def configure_scenario(scenario_id: str, spacing: int) -> None:
    global SCENARIO_ID, SPREAD_NS, DENSE_STAGE_ROOT, STAGE_LABEL
    spec = SCENARIO_SPECS[scenario_id]
    if int(spacing) not in spec["allowed_spacings"]:
        raise ValueError(f"{scenario_id} does not authorize comb-{int(spacing)}.")
    SCENARIO_ID = scenario_id
    SPREAD_NS = float(spec["spread_ns"])
    DENSE_STAGE_ROOT = str(spec["stage_root"])
    STAGE_LABEL = str(spec["stage_label"])


def dense_root(output: Path, spacing: int) -> Path:
    return output / DENSE_STAGE_ROOT / f"comb{int(spacing)}"


def candidate_id(suffix: str) -> str:
    return f"{SCENARIO_ID}_{suffix}"


def active_baseline_families() -> tuple[str, ...]:
    if SCENARIO_ID in ("A100", "A300"):
        return (
            *BASELINE_FAMILIES[:-1],
            A100_ADDITIONAL_BASELINE_FAMILY,
            BASELINE_FAMILIES[-1],
        )
    return BASELINE_FAMILIES


def manifest_filename(suffix: str) -> str:
    return f"{STAGE_LABEL}_link_{suffix}"


def dynamic_alias_delays(spacing: int) -> tuple[int, ...]:
    pilot_period = K // int(spacing)
    if pilot_period % N_TX != 0:
        raise ValueError("Pilot alias period must be divisible by N_tx.")
    step = pilot_period // N_TX
    return tuple(step * branch for branch in range(N_TX))


def useful_symbol_nt_minus_one_delays_ns() -> tuple[float, ...]:
    return continuous_arithmetic_delays_ns(
        USEFUL_SYMBOL_NS / (N_TX - 1),
        N_TX,
    )


def _grid_row(
    candidate_id: str,
    family: str,
    rule: str,
    state: Sequence[int],
    covariance: np.ndarray,
    required_pair: int,
    required_fold: int,
    spacing: int,
) -> dict[str, object]:
    moments = effective_moments_lag(state, covariance)
    return grid_candidate_row(
        SCENARIO_ID,
        candidate_id,
        family,
        rule,
        state,
        moments["m2_eff"],
        moments["m4_eff"],
        required_pair,
        required_fold,
        pilot_spacing=int(spacing),
    )


def run_search(
    output: Path,
    spacing: int,
    search_states_count: int,
    geometry_states_count: int,
) -> None:
    stage = dense_root(output, spacing) / "search"
    stage.mkdir(parents=True, exist_ok=True)
    frozen_path = stage / "frozen_candidates.csv"
    if frozen_path.exists():
        raise RuntimeError("Dense-DMRS search output already exists.")

    pilot_period = K // int(spacing)
    _, covariance = physical_covariance(SPREAD_NS)
    _, grid = resource_and_grid(spacing)
    support = support_for_scenario(SPREAD_NS)
    support_ns = float(support["support_width_s"]) * 1e9
    requirements = thick_sidon_requirements(
        support_ns,
        GRID_NS,
        K,
        N_TX,
        pilot_period,
    )
    required_pair = int(requirements["required_pair_gap_indices"])
    required_fold = int(requirements["required_fold_gap_indices"])

    fixed_rows = [
        _grid_row(
            candidate_id("B0_QC"),
            "B0_QC",
            "fixed_reference",
            QC,
            covariance,
            required_pair,
            required_fold,
            spacing,
        ),
        _grid_row(
            candidate_id("S0_SIDON"),
            "S0_SIDON",
            "fixed_reference",
            SIDON,
            covariance,
            required_pair,
            required_fold,
            spacing,
        ),
        _grid_row(
            candidate_id("AP_TU_NT"),
            "AP_TU_NT",
            "fixed_system_useful_symbol_nt",
            AP_TU_NT,
            covariance,
            required_pair,
            required_fold,
            spacing,
        ),
        _grid_row(
            candidate_id("AP_TALIAS_NT"),
            "AP_TALIAS_NT",
            f"current_comb{spacing}_pilot_alias_period_nt",
            dynamic_alias_delays(spacing),
            covariance,
            required_pair,
            required_fold,
            spacing,
        ),
    ]
    fixed_rows.extend(
        [
            continuous_candidate(
                SCENARIO_ID,
                "AP_RMS_T1",
                continuous_arithmetic_delays_ns(SPREAD_NS, N_TX),
                covariance,
                grid,
            ),
            continuous_candidate(
                SCENARIO_ID,
                "AP_TEPS_T1",
                continuous_arithmetic_delays_ns(support_ns, N_TX),
                covariance,
                grid,
            ),
        ]
    )
    if SCENARIO_ID in ("A100", "A300"):
        tu_ntm1 = continuous_candidate(
            SCENARIO_ID,
            A100_ADDITIONAL_BASELINE_FAMILY,
            useful_symbol_nt_minus_one_delays_ns(),
            covariance,
            grid,
        )
        tu_ntm1["selection_rule"] = "fixed_system_useful_symbol_nt_minus_1"
        fixed_rows.append(tu_ntm1)

    arithmetic = []
    for state in arithmetic_delay_sets(N_TX, K):
        residues = tuple(int(value) % pilot_period for value in state)
        if len(set(residues)) == N_TX:
            arithmetic.append(state)
    ap_m2, ap_m4 = batch_moments(arithmetic, covariance)
    ap_selected = select_meff_indices(ap_m2, ap_m4, limit=2)
    ap_rows = [
        grid_candidate_row(
            SCENARIO_ID,
            candidate_id(f"AP_T2_{position:02d}"),
            "AP_T2_CTRL",
            ";".join(sorted(rules)),
            arithmetic[index],
            ap_m2[index],
            ap_m4[index],
            required_pair,
            required_fold,
            pilot_spacing=int(spacing),
        )
        for position, (index, rules) in enumerate(ap_selected.items())
    ]

    geo_states = geometry_states(
        stable_seed(SEED, SCENARIO_ID, f"comb{spacing}", "geometry"),
        geometry_states_count,
        required_fold,
        pilot_period=pilot_period,
    )
    geo_pair, geo_fold, geo_hard, geo_selected = geometry_selections(
        geo_states,
        required_pair,
        required_fold,
        limit=3,
        pilot_period=pilot_period,
    )
    geo_m2, geo_m4 = batch_moments(geo_states, covariance)
    write_large_candidate_csv(
        stage / "geometry_all_candidates.csv",
        SCENARIO_ID,
        "GEO_T1_CTRL",
        geo_states,
        geo_m2,
        geo_m4,
        geo_pair,
        geo_fold,
        geo_hard,
        pilot_period=pilot_period,
    )
    geo_rows = [
        grid_candidate_row(
            SCENARIO_ID,
            candidate_id(f"GEO_T1_{position:02d}"),
            "GEO_T1_CTRL",
            ";".join(sorted(rules)),
            geo_states[index],
            geo_m2[index],
            geo_m4[index],
            required_pair,
            required_fold,
            pilot_spacing=int(spacing),
        )
        for position, (index, rules) in enumerate(geo_selected.items())
    ]

    initial: list[Sequence[int]] = [QC, SIDON, AP_TU_NT, dynamic_alias_delays(spacing)]
    initial.extend(arithmetic)
    for physical_spacing_ns in (SPREAD_NS, support_ns):
        projection = np.rint(
            delay_ns_to_grid_coordinates(
                continuous_arithmetic_delays_ns(physical_spacing_ns, N_TX),
                K,
                SCS_HZ,
            )
        ).astype(int) % K
        initial.append(tuple(int(value) for value in projection))
    initial.extend(geo_states[index] for index in geo_selected)
    meff_states = unique_random_delay_sets(
        stable_seed(SEED, SCENARIO_ID, f"comb{spacing}", "t2"),
        search_states_count,
        initial=initial,
        n_tx=N_TX,
        period=K,
        pilot_period=pilot_period,
    )
    meff_m2, meff_m4 = batch_moments(meff_states, covariance)
    write_large_candidate_csv(
        stage / "meff_all_candidates.csv",
        SCENARIO_ID,
        "MEFF_T2_CAND",
        meff_states,
        meff_m2,
        meff_m4,
        pilot_period=pilot_period,
    )
    meff_selected = select_meff_indices(meff_m2, meff_m4, limit=7)
    meff_rows = [
        grid_candidate_row(
            SCENARIO_ID,
            candidate_id(f"MEFF_T2_{position:02d}"),
            "MEFF_T2_CAND",
            ";".join(sorted(rules)),
            meff_states[index],
            meff_m2[index],
            meff_m4[index],
            required_pair,
            required_fold,
            pilot_spacing=int(spacing),
        )
        for position, (index, rules) in enumerate(meff_selected.items())
    ]

    rows = [*fixed_rows, *ap_rows, *geo_rows, *meff_rows]
    add_equivalent_delay_metrics(rows, SPREAD_NS)
    unique_rows: dict[tuple[str, str], dict[str, object]] = {}
    for row in rows:
        identity = (str(row["family"]), str(row["delay_grid_coordinates"]))
        unique_rows[identity] = row
    rows = list(unique_rows.values())
    if len(rows) > 18:
        raise RuntimeError("Dense-DMRS search froze more than 18 candidates.")
    write_csv_rows(rows, frozen_path)
    save_json(
        {
            "scenario_id": SCENARIO_ID,
            "dmrs_spacing_subcarriers": int(spacing),
            "pilot_period_indices": pilot_period,
            "support_width_ns": support_ns,
            "required_pair_gap_indices": required_pair,
            "required_fold_gap_indices": required_fold,
            "search_states": int(search_states_count),
            "geometry_states": int(geometry_states_count),
            "seed": SEED,
            "seed_derivation_includes_spacing": True,
            "meff_bands_db": list(MEFF_BANDS_DB),
            "candidate_count": len(rows),
            "candidate_freeze_precedes_outage": True,
        },
        stage / "resolved_experiment.json",
    )


def run_outage(output: Path, spacing: int, n_mc: int, repeats: int) -> None:
    root = dense_root(output, spacing)
    rows = read_csv(root / "search" / "frozen_candidates.csv")
    stage = root / "outage"
    stage.mkdir(parents=True, exist_ok=True)
    table = iqam_table(output)
    snr_grid = np.arange(0.0, 28.0 + 1e-9, 0.5)
    started = time.time()
    curves, thresholds = scenario_outage(
        rows,
        SPREAD_NS,
        table,
        n_mc,
        stable_seed(SEED, SCENARIO_ID, f"{STAGE_LABEL}_dense_outage"),
        snr_grid,
    )
    np.savez(
        stage / "paired_threshold_snr_db.npz",
        candidate_ids=np.asarray([row["candidate_id"] for row in rows]),
        threshold_snr_db=thresholds,
        snr_grid_db=snr_grid,
    )
    targets = {
        (index, target): outage_target(curves[index], snr_grid, target, n_mc)
        for index in range(len(rows))
        for target in (0.10, 0.01)
    }
    family_index = {
        row["family"]: index
        for index, row in enumerate(rows)
        if row["family"] in ("AP_RMS_T1", "AP_TEPS_T1")
    }
    bootstrap_rows = []
    target_rows = []
    for target in (0.10, 0.01):
        rms = family_index["AP_RMS_T1"]
        teps = family_index["AP_TEPS_T1"]
        baseline = rms if targets[(rms, target)] <= targets[(teps, target)] else teps
        for index, row in enumerate(rows):
            gains = paired_bootstrap_gain(
                thresholds[index],
                thresholds[baseline],
                curves[index],
                curves[baseline],
                snr_grid,
                target,
                repeats,
                stable_seed(
                    SEED,
                    SCENARIO_ID,
                    f"comb{spacing}",
                    row["candidate_id"],
                    target,
                    "bootstrap",
                ),
            )
            finite = gains[np.isfinite(gains)]
            if finite.size != int(repeats):
                raise RuntimeError("Dense-DMRS bootstrap produced non-finite values.")
            low, high = np.percentile(finite, [2.5, 97.5])
            bootstrap_rows.append(
                {
                    "scenario_id": SCENARIO_ID,
                    "dmrs_spacing_subcarriers": int(spacing),
                    "candidate_id": row["candidate_id"],
                    "target_outage_probability": target,
                    "baseline_candidate_id": rows[baseline]["candidate_id"],
                    "gain_snr_db": targets[(baseline, target)]
                    - targets[(index, target)],
                    "gain_ci95_low_db": float(low),
                    "gain_ci95_high_db": float(high),
                    "bootstrap_repeats": int(repeats),
                }
            )
    bootstrap_lookup = {
        (row["candidate_id"], float(row["target_outage_probability"])): row
        for row in bootstrap_rows
    }
    for index, row in enumerate(rows):
        ten = bootstrap_lookup[(row["candidate_id"], 0.10)]
        one = bootstrap_lookup[(row["candidate_id"], 0.01)]
        target_rows.append(
            {
                **row,
                "dmrs_spacing_subcarriers": int(spacing),
                "outage10_snr_db": targets[(index, 0.10)],
                "outage1_snr_db": targets[(index, 0.01)],
                "ap_t1_base_10_candidate_id": ten["baseline_candidate_id"],
                "ap_t1_base_1_candidate_id": one["baseline_candidate_id"],
                "gain_vs_ap_t1_base_10_db": ten["gain_snr_db"],
                "gain_vs_ap_t1_base_10_ci95_low_db": ten["gain_ci95_low_db"],
                "gain_vs_ap_t1_base_10_ci95_high_db": ten["gain_ci95_high_db"],
                "gain_vs_ap_t1_base_1_db": one["gain_snr_db"],
                "gain_vs_ap_t1_base_1_ci95_low_db": one["gain_ci95_low_db"],
                "gain_vs_ap_t1_base_1_ci95_high_db": one["gain_ci95_high_db"],
            }
        )
    curve_rows = [
        {
            "scenario_id": SCENARIO_ID,
            "dmrs_spacing_subcarriers": int(spacing),
            "candidate_id": row["candidate_id"],
            "family": row["family"],
            "snr_db": float(snr_db),
            "outage_probability": float(curves[index, snr_index]),
            "n_mc": int(n_mc),
        }
        for index, row in enumerate(rows)
        for snr_index, snr_db in enumerate(snr_grid)
    ]
    write_csv_rows(curve_rows, stage / "outage_curves.csv")
    write_csv_rows(target_rows, stage / "outage_targets.csv")
    write_csv_rows(bootstrap_rows, stage / "paired_bootstrap.csv")
    save_json(
        {
            "scenario_id": SCENARIO_ID,
            "dmrs_spacing_subcarriers": int(spacing),
            "n_mc": int(n_mc),
            "bootstrap_repeats": int(repeats),
            "snr_grid_db": [float(value) for value in snr_grid],
            "common_random_numbers_across_candidates_and_dense_spacings": True,
            "elapsed_seconds": time.time() - started,
        },
        stage / "resolved_experiment.json",
    )


def _decode_array(row: dict[str, str], field: str) -> list[float] | None:
    text = str(row.get(field, "")).strip()
    if not text:
        return None
    values = json.loads(text)
    return [float(value) for value in values]


def _manifest_candidates(rows: list[dict[str, str]], spacing: int) -> list[dict]:
    selected: dict[str, dict[str, object]] = {}
    for family in active_baseline_families():
        matches = [row for row in rows if row["family"] == family]
        if len(matches) != 1:
            raise RuntimeError(f"Expected one dense-DMRS {family} baseline.")
        selected[matches[0]["candidate_id"]] = {
            "row": matches[0],
            "roles": ["baseline"],
            "outage_winner_for": [],
        }
    s0 = [row for row in rows if row["family"] == "S0_SIDON"]
    if len(s0) != 1:
        raise RuntimeError("Expected one dense-DMRS S0 reference.")
    selected[s0[0]["candidate_id"]] = {
        "row": s0[0],
        "roles": ["historical_reference"],
        "outage_winner_for": ["10pct", "1pct"],
    }
    for family in SEARCH_FAMILIES:
        family_rows = [row for row in rows if row["family"] == family]
        for label, _, field in TARGETS:
            winner = min(family_rows, key=lambda row: float(row[field]))
            item = selected.setdefault(
                winner["candidate_id"],
                {
                    "row": winner,
                    "roles": ["outage_family_winner"],
                    "outage_winner_for": [],
                },
            )
            item["outage_winner_for"].append(label)

    _, covariance = physical_covariance(SPREAD_NS)
    _, grid = resource_and_grid(spacing)
    pilots = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    data = local_indices_for_subcarriers(grid, grid.data_coordinates[:, 1])
    metric = FrequencyCEMetric(covariance, pilots, data, N_TX)
    output = []
    for candidate_id, item in selected.items():
        row = item["row"]
        coordinates = _decode_array(row, "delay_grid_coordinates")
        delay_ns = _decode_array(row, "delay_ns")
        indices = _decode_array(row, "delay_indices")
        if coordinates is None or len(coordinates) != N_TX:
            raise RuntimeError(f"{candidate_id} has invalid delay coordinates.")
        ce10 = metric.evaluate_grid_coordinates(
            coordinates,
            float(row["outage10_snr_db"]),
        )
        ce1 = metric.evaluate_grid_coordinates(
            coordinates,
            float(row["outage1_snr_db"]),
        )
        output.append(
            {
                "curve_id": f"{SCENARIO_ID}:comb{spacing}:{candidate_id}",
                "scenario_id": SCENARIO_ID,
                "profile": "A",
                "delay_spread_ns": SPREAD_NS,
                "candidate_id": candidate_id,
                "family": row["family"],
                "roles": item["roles"],
                "outage_winner_for": sorted(set(item["outage_winner_for"])),
                "represented_families": [row["family"]],
                "represented_candidate_ids": [candidate_id],
                "outage_winner_roles": [
                    f"{row['family']}:{label}"
                    for label in sorted(set(item["outage_winner_for"]))
                ],
                "selection_rule": row["selection_rule"],
                "delay_grid_coordinates": coordinates,
                "delay_indices": (
                    [int(round(value)) for value in indices]
                    if indices is not None
                    else None
                ),
                "delay_ns": delay_ns or [],
                "grid_aligned": row["grid_aligned"] == "True",
                "residues": _decode_array(row, "residues"),
                "lifts": _decode_array(row, "lifts"),
                "pilot_rank": int(row["pilot_rank"]),
                "pilot_condition_number": float(row["pilot_condition_number"]),
                "outage10_snr_db": float(row["outage10_snr_db"]),
                "outage1_snr_db": float(row["outage1_snr_db"]),
                "ideal_outage10_ce_nmse_db": float(ce10["nmse_db"]),
                "ideal_outage1_ce_nmse_db": float(ce1["nmse_db"]),
                "phase_denominator_active_subcarriers": K,
                "phase_reference": "first_active_subcarrier",
                "precoder_normalized": False,
                "dmrs_spacing_subcarriers": int(spacing),
            }
        )
    family_order = (
        *active_baseline_families(),
        "S0_SIDON",
        *SEARCH_FAMILIES,
    )
    output.sort(
        key=lambda row: (
            family_order.index(str(row["family"])),
            str(row["candidate_id"]),
        )
    )
    deduplicated: dict[tuple[float, ...], dict] = {}
    for row in output:
        key = tuple(
            round(float(value), 12)
            for value in row["delay_grid_coordinates"]
        )
        if key not in deduplicated:
            deduplicated[key] = row
            continue
        kept = deduplicated[key]
        kept["roles"] = sorted(set([*kept["roles"], *row["roles"]]))
        kept["outage_winner_for"] = sorted(
            set([*kept["outage_winner_for"], *row["outage_winner_for"]])
        )
        kept["represented_families"] = [
            *kept["represented_families"],
            *[
                value
                for value in row["represented_families"]
                if value not in kept["represented_families"]
            ],
        ]
        kept["represented_candidate_ids"] = [
            *kept["represented_candidate_ids"],
            *[
                value
                for value in row["represented_candidate_ids"]
                if value not in kept["represented_candidate_ids"]
            ],
        ]
        kept["outage_winner_roles"] = sorted(
            set([*kept["outage_winner_roles"], *row["outage_winner_roles"]])
        )
    output = list(deduplicated.values())
    if not 8 <= len(output) <= 13:
        raise RuntimeError(f"Dense-DMRS manifest has {len(output)} curves.")
    return output


def run_manifest(output: Path, spacing: int) -> None:
    root = dense_root(output, spacing)
    rows = read_csv(root / "outage" / "outage_targets.csv")
    candidates = _manifest_candidates(rows, spacing)
    stage = root / "manifest"
    stage.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema": (
            "plan028-a300-comb6-manifest-v1"
            if SCENARIO_ID == "A300"
            else f"plan027-{STAGE_LABEL}-{SCENARIO_ID.lower()}-dense-dmrs-manifest-v1"
        ),
        "status": (
            "researcher_confirmed_in_plan028"
            if SCENARIO_ID == "A300"
            else "researcher_confirmed_in_plan027_addendum"
        ),
        "run_id": output.name,
        "selection_scope": (
            f"{SCENARIO_ID} {len(active_baseline_families())} baselines, "
            "S0, and deduplicated 10%/1% "
            "ideal-outage "
            "family winners after current-density redesign"
        ),
        "curve_count": len(candidates),
        "scene_counts": {SCENARIO_ID: len(candidates)},
        "scene_order": [SCENARIO_ID],
        "execution_order": [SCENARIO_ID],
        "dmrs_spacing_subcarriers": int(spacing),
        "output_subdir": f"{DENSE_STAGE_ROOT}/comb{spacing}/link",
        "physical_definition": {
            "K_active_subcarriers": K,
            "subcarrier_spacing_hz": SCS_HZ,
            "phase_reference": "first_active_subcarrier",
            "phase_denominator": K,
            "delay_input_field": "delay_grid_coordinates",
            "continuous_coordinates_must_not_be_rounded": True,
            "precoder_normalized": False,
        },
        "trial_budget": {
            "smoke_first_scene_only": True,
            "smoke_trials_per_snr": 20,
            "prescan_trials_per_snr": 400,
            "prescan_step_db": 0.5,
            "initial_bler_gap_from_outage_db": 5.0,
            "maximum_extension_each_target_db": 2.0,
            "refine_trials_per_snr": 3000,
            "refine_step_db": 0.25,
            "minimum_errors_in_half_to_twice_target_band": 30,
        },
        "researcher_authorization": (
            "2026-07-27 request to append A30 comb-12/6 estimated-CSI BLER "
            "to plan-027/result-027"
            if SCENARIO_ID == "A30"
            else "2026-07-30 request to append A100 comb-6 estimated-CSI BLER "
            "and 10%/1% scan plus scatter plots to plan-027/result-027"
            if SCENARIO_ID == "A100"
            else
            "2026-08-03 request to run A300 comb-6 comparison with A100-comb-6 "
            "system, search, outage, and link budgets under plan-028"
        ),
        "candidates": candidates,
    }
    manifest_path = stage / manifest_filename("manifest.json")
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    digest = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    (stage / manifest_filename("manifest.sha256")).write_text(
        f"{digest}  {manifest_path.name}\n",
        encoding="utf-8",
    )
    approval = {
        "approved": True,
        "manifest_sha256": digest,
        "approved_curve_count": len(candidates),
        "approved_scenario": SCENARIO_ID,
        "approved_dmrs_spacing_subcarriers": int(spacing),
        "approved_trial_budget": manifest["trial_budget"],
        "researcher_confirmation": manifest["researcher_authorization"],
    }
    (stage / manifest_filename("approval.json")).write_text(
        json.dumps(approval, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    csv_rows = []
    for row in candidates:
        encoded = dict(row)
        for field in (
            "roles",
            "outage_winner_for",
            "represented_families",
            "represented_candidate_ids",
            "outage_winner_roles",
            "delay_grid_coordinates",
            "delay_indices",
            "delay_ns",
            "residues",
            "lifts",
        ):
            encoded[field] = json.dumps(row[field], separators=(",", ":"))
        csv_rows.append(encoded)
    write_csv_rows(csv_rows, stage / manifest_filename("manifest.csv"))


def _read_target_summary(root: Path) -> list[dict[str, str]]:
    path = root / "link" / SCENARIO_ID / "final" / "target_summary.csv"
    return read_csv(path) if path.exists() else []


def _scatter_points(output: Path, spacing: int) -> list[dict[str, object]]:
    root = dense_root(output, spacing)
    manifest = json.loads(
        (
            root
            / "manifest"
            / manifest_filename("manifest.json")
        ).read_text(encoding="utf-8")
    )
    summaries = {
        (row["candidate_id"], row["target"]): row
        for row in _read_target_summary(root)
    }
    points = []
    for row in manifest["candidates"]:
        for label, target, outage_field in TARGETS:
            points.append(
                {
                    "scenario_id": SCENARIO_ID,
                    "dmrs_spacing_subcarriers": int(spacing),
                    "candidate_id": row["candidate_id"],
                    "family": row["family"],
                    "target": label,
                    "target_probability": target,
                    "point_kind": "ideal_outage",
                    "snr_db": float(row[outage_field]),
                    "ce_nmse_db": float(
                        row[
                            "ideal_outage10_ce_nmse_db"
                            if label == "10pct"
                            else "ideal_outage1_ce_nmse_db"
                        ]
                    ),
                    "ci95_low_db": "",
                    "ci95_high_db": "",
                    "closed": True,
                }
            )
            summary = summaries.get((row["candidate_id"], label))
            if summary is not None:
                points.append(
                    {
                        "scenario_id": SCENARIO_ID,
                        "dmrs_spacing_subcarriers": int(spacing),
                        "candidate_id": row["candidate_id"],
                        "family": row["family"],
                        "target": label,
                        "target_probability": target,
                        "point_kind": "estimated_csi_bler",
                        "snr_db": float(summary["target_snr_db"]),
                        "ce_nmse_db": float(summary["target_ce_nmse_db"]),
                        "ci95_low_db": float(summary["ci95_low_db"]),
                        "ci95_high_db": float(summary["ci95_high_db"]),
                        "closed": True,
                    }
                )
    return points


def _plot_scatter(points: list[dict[str, object]], path: Path, target: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator

    selected = [row for row in points if row["target"] == target]
    figure, axis = plt.subplots(figsize=(11.5, 7.0))
    identifiers = sorted({str(row["candidate_id"]) for row in selected})
    for candidate_id in identifiers:
        rows = [row for row in selected if row["candidate_id"] == candidate_id]
        ideal = next(row for row in rows if row["point_kind"] == "ideal_outage")
        solid = next(
            (row for row in rows if row["point_kind"] == "estimated_csi_bler"),
            None,
        )
        style = curve_style(candidate_id, str(ideal["family"]))
        if solid is not None:
            axis.plot(
                [float(ideal["snr_db"]), float(solid["snr_db"])],
                [float(ideal["ce_nmse_db"]), float(solid["ce_nmse_db"])],
                color=style["color"],
                linestyle=style["linestyle"],
                linewidth=2.5,
                alpha=0.9,
            )
        axis.scatter(
            [float(ideal["snr_db"])],
            [float(ideal["ce_nmse_db"])],
            facecolors="none",
            edgecolors=style["color"],
            marker=style["marker"],
            s=96,
            linewidths=2.5,
            label=curve_label(candidate_id),
        )
        if solid is not None:
            x = float(solid["snr_db"])
            axis.errorbar(
                [x],
                [float(solid["ce_nmse_db"])],
                xerr=[
                    [x - float(solid["ci95_low_db"])],
                    [float(solid["ci95_high_db"]) - x],
                ],
                fmt=style["marker"],
                color=style["color"],
                markerfacecolor=style["color"],
                markeredgecolor=style["color"],
                markersize=8,
                capsize=2.5,
                linewidth=2.5,
            )
    spacing = int(selected[0]["dmrs_spacing_subcarriers"])
    axis.set_xlabel("Target SNR (dB)", fontsize=16)
    axis.set_ylabel("matched CE NMSE at target SNR (dB)", fontsize=16)
    axis.set_title(
        f"{SCENARIO_ID} comb-{spacing} {target} outage–CE–BLER link"
        , fontsize=18
    )
    axis.tick_params(axis="both", which="both", labelsize=14)
    axis.xaxis.set_major_locator(MultipleLocator(1.0))
    axis.xaxis.set_minor_locator(MultipleLocator(0.25))
    axis.yaxis.set_major_locator(MultipleLocator(2.0))
    axis.yaxis.set_minor_locator(MultipleLocator(0.5))
    axis.grid(True, which="major", alpha=0.35)
    axis.grid(True, which="minor", alpha=0.15)
    axis.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=16)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def _write_delay_density_table(output: Path) -> None:
    """Write the selected physical-delay sets across available densities."""
    manifests = [(24, output / "e4_link_gate" / "e4_link_manifest.json")]
    for density in SCENARIO_SPECS[SCENARIO_ID]["allowed_spacings"]:
        manifests.append(
            (
                int(density),
                dense_root(output, int(density))
                / "manifest"
                / manifest_filename("manifest.json"),
            )
        )
    rows: list[dict[str, object]] = []
    for spacing, path in manifests:
        if not path.exists():
            continue
        manifest = json.loads(path.read_text(encoding="utf-8"))
        for candidate in manifest["candidates"]:
            if not str(candidate["candidate_id"]).startswith(f"{SCENARIO_ID}_"):
                continue
            represented = candidate.get("represented_families") or [
                candidate["family"]
            ]
            rows.append(
                {
                    "scenario_id": SCENARIO_ID,
                    "dmrs_spacing_subcarriers": spacing,
                    "candidate_id": candidate["candidate_id"],
                    "family": candidate["family"],
                    "represented_families": json.dumps(
                        represented, separators=(",", ":")
                    ),
                    "delay_grid_coordinates": json.dumps(
                        candidate.get("delay_grid_coordinates", []),
                        separators=(",", ":"),
                    ),
                    "delay_indices": json.dumps(
                        candidate.get("delay_indices", []),
                        separators=(",", ":"),
                    ),
                    "delay_ns": json.dumps(
                        candidate["delay_ns"], separators=(",", ":")
                    ),
                    "source_manifest": str(path.relative_to(output)),
                }
            )
    destination = output / DENSE_STAGE_ROOT / "final"
    destination.mkdir(parents=True, exist_ok=True)
    write_csv_rows(rows, destination / "selected_delay_sets_by_density.csv")


def run_analyze(output: Path, spacing: int) -> None:
    root = dense_root(output, spacing)
    final = root / "final"
    final.mkdir(parents=True, exist_ok=True)
    points = _scatter_points(output, spacing)
    write_csv_rows(points, final / "outage_ce_link_scatter_points.csv")
    for target in ("10pct", "1pct"):
        _plot_scatter(
            points,
            final
            / (
                f"{SCENARIO_ID.lower()}_comb{spacing}_{target}"
                "_outage_ce_link.png"
            ),
            target,
        )
    _write_delay_density_table(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage",
        choices=("search", "outage", "manifest", "analyze"),
        required=True,
    )
    parser.add_argument(
        "--scenario",
        choices=tuple(SCENARIO_SPECS),
        default="A30",
    )
    parser.add_argument("--spacing", type=int, choices=DENSE_SPACINGS, required=True)
    parser.add_argument("--run-id", default=DEFAULT_RUN_ID)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--search-states", type=int, default=SEARCH_STATES)
    parser.add_argument("--geometry-states", type=int, default=GEOMETRY_STATES)
    parser.add_argument("--outage-samples", type=int, default=200_000)
    parser.add_argument("--bootstrap-repeats", type=int, default=1_000)
    args = parser.parse_args()
    configure_scenario(str(args.scenario), int(args.spacing))
    output = args.output_root / str(args.run_id)
    if args.stage == "search":
        run_search(
            output,
            int(args.spacing),
            int(args.search_states),
            int(args.geometry_states),
        )
    elif args.stage == "outage":
        run_outage(
            output,
            int(args.spacing),
            int(args.outage_samples),
            int(args.bootstrap_repeats),
        )
    elif args.stage == "manifest":
        run_manifest(output, int(args.spacing))
    else:
        run_analyze(output, int(args.spacing))


if __name__ == "__main__":
    main()
