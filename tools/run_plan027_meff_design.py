"""Plan-027 Phase 0, E1/E2 search, freeze, and ideal-CSI outage."""

from __future__ import annotations

import argparse
import csv
import hashlib
import inspect
import json
import math
import platform
import sys
import time
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ChannelConfig, ResourceConfig
from cdd_lls.design import (
    FrequencyCEMetric,
    arithmetic_delay_sets,
    canonical_delay_set,
    continuous_arithmetic_delays_ns,
    delay_indices_to_ns,
    delay_ns_to_grid_coordinates,
    effective_moments_direct,
    effective_moments_lag,
    equivalent_cdd_delay_spread,
    fold_min_gap_indices,
    pair_sum_min_gap_indices,
    residues_and_lifts,
    tdl_effective_support,
    thick_sidon_requirements,
    unique_random_delay_sets,
)
from cdd_lls.phy.estimators import tdl_active_frequency_covariance
from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
from tools.run_plan025_delay_matched_tdl import environment_receipt, stable_seed, write_csv_rows
from tools.run_track_b_pilot_scan import (
    HIST_DB_MAX,
    HIST_DB_MIN,
    IQAM_DB_MAX,
    IQAM_DB_MIN,
    IQAM_DB_STEP,
    R_SE,
    build_iqam_table,
    build_shift_table,
)


K = 576
N_TX = 8
SCS_HZ = 30e3
PILOT_SPACING = 24
PILOT_PERIOD = K // PILOT_SPACING
GRID_NS = 1e9 / (K * SCS_HZ)
SEED = 20260727
SEARCH_STATES = 200_000
GEOMETRY_STATES = 200_000
N_MC = 200_000
BOOTSTRAP_REPEATS = 1_000
MEFF_BANDS_DB = (0.10, 0.25, 0.50)
DEFAULT_OUTPUT_ROOT = ROOT / "outputs" / "experiment027_meff_sidon"
DEFAULT_RUN_ID = "20260726_main"
QC = (0, 9, 18, 27, 36, 45, 54, 63)
SIDON = (0, 1, 3, 7, 12, 20, 30, 65)
AP_TU_NT = tuple(index * (K // N_TX) for index in range(N_TX))
AP_TALIAS_NT = tuple(
    index * (PILOT_PERIOD // N_TX) for index in range(N_TX)
)
SYSTEM_BASELINES = (
    ("AP_TU_NT", AP_TU_NT, "useful_symbol_period_over_n_tx"),
    (
        "AP_TALIAS_NT",
        AP_TALIAS_NT,
        "pilot_alias_free_period_over_n_tx",
    ),
)
SCENARIOS = (
    ("A1", 1.0),
    ("A5", 5.0),
    ("A10", 10.0),
    ("A30", 30.0),
    ("A100", 100.0),
)
EXPECTED_SUPPORT_NS = {
    "A1": 4.7966,
    "A5": 23.983,
    "A10": 47.966,
    "A30": 143.898,
    "A100": 479.660,
}


def encode(values: Sequence[object]) -> str:
    return json.dumps(list(values), separators=(",", ":"))


def decode(text: str) -> tuple[float, ...]:
    return tuple(float(value) for value in json.loads(text))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def save_json(value: object, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding="utf-8")


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resource_and_grid(dmrs_spacing: int = PILOT_SPACING):
    resource = ResourceConfig(
        carrier_bandwidth_mhz=100.0,
        scs_khz=30,
        n_fft=4096,
        n_prbs=48,
        pdsch_n_symbols=10,
        dmrs_symbol_indices=[2, 7],
        dmrs_spacing_sc=int(dmrs_spacing),
        dmrs_offset_sc=0,
        cyclic_prefix_length=288,
    )
    return resource, build_resource_grid(resource)


def channel_config(delay_spread_ns: float) -> ChannelConfig:
    return ChannelConfig(
        backend="sionna_tdl",
        model="3gpp_tr38901_tdl",
        tdl_profile="A",
        delay_spread_ns=float(delay_spread_ns),
        carrier_frequency_hz=3.5e9,
        ue_speed_kmh=0.0,
        num_sinusoids=20,
    )


def physical_covariance(delay_spread_ns: float) -> tuple[object, np.ndarray]:
    _, grid = resource_and_grid()
    covariance = tdl_active_frequency_covariance(grid, channel_config(delay_spread_ns))
    return grid, np.asarray(covariance, dtype=np.complex128)


def tdl_pdp(delay_spread_ns: float) -> tuple[np.ndarray, np.ndarray]:
    from sionna.phy.channel.tr38901 import TDL

    tdl = TDL(
        model="A",
        delay_spread=float(delay_spread_ns) * 1e-9,
        carrier_frequency=3.5e9,
        num_sinusoids=20,
        min_speed=0.0,
        max_speed=0.0,
        num_rx_ant=1,
        num_tx_ant=1,
        precision="double",
    )
    return (
        np.asarray(tdl.delays.numpy(), dtype=np.float64),
        np.asarray(tdl.mean_powers.numpy(), dtype=np.float64),
    )


def support_for_scenario(delay_spread_ns: float) -> dict[str, float | int]:
    delays_s, powers = tdl_pdp(delay_spread_ns)
    return tdl_effective_support(delays_s, powers, epsilon=0.01)


def add_equivalent_delay_metrics(
    rows: Sequence[dict[str, object]],
    delay_spread_ns: float,
) -> None:
    physical_delays_s, physical_powers = tdl_pdp(delay_spread_ns)
    for row in rows:
        metrics = equivalent_cdd_delay_spread(
            physical_delays_s,
            physical_powers,
            decode(str(row["delay_grid_coordinates"])),
            period=K,
            subcarrier_spacing_hz=SCS_HZ,
            contained_probability=0.99,
        )
        row.update(
            {
                "equivalent_circular_rms_delay_spread_ns": 1e9
                * float(metrics["equivalent_circular_rms_delay_spread_s"]),
                "equivalent_99pct_circular_support_width_ns": 1e9
                * float(metrics["equivalent_circular_support_width_s"]),
                "equivalent_delay_period_ns": 1e9
                * float(metrics["equivalent_delay_period_s"]),
            }
        )


def candidate_geometry(
    indices: Sequence[int],
    pilot_spacing: int = PILOT_SPACING,
) -> dict[str, object]:
    pilot_period = K // int(pilot_spacing)
    key = canonical_delay_set(indices, K)
    residues, lifts = residues_and_lifts(key, pilot_period)
    pilots = np.arange(0, K, int(pilot_spacing), dtype=np.float64)
    phase = np.exp(
        -2j
        * np.pi
        * pilots[:, None]
        * np.asarray(key, dtype=np.float64)[None, :]
        / float(K)
    )
    singular = np.linalg.svd(phase, compute_uv=False)
    return {
        "delay_indices": encode(key),
        "delay_grid_coordinates": encode(key),
        "delay_ns": encode(delay_indices_to_ns(key, K, SCS_HZ)),
        "residues": encode(residues),
        "lifts": encode(lifts),
        "canonical_key": encode(key),
        "pair_sum_gap_indices": pair_sum_min_gap_indices(key, K),
        "fold_gap_indices": fold_min_gap_indices(key, pilot_period),
        "grid_aligned": True,
        "pilot_rank": int(np.linalg.matrix_rank(phase, tol=1e-10)),
        "pilot_condition_number": float(singular[0] / singular[-1]),
    }


def continuous_candidate(
    scenario_id: str,
    family: str,
    delay_ns: Sequence[float],
    covariance: np.ndarray,
    grid,
) -> dict[str, object]:
    coordinates = delay_ns_to_grid_coordinates(delay_ns, K, SCS_HZ)
    moments = effective_moments_lag(coordinates, covariance)
    pilots = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    phase = np.exp(
        -2j * np.pi * pilots[:, None] * coordinates[None, :] / float(K)
    )
    singular = np.linalg.svd(phase, compute_uv=False)
    return {
        "scenario_id": scenario_id,
        "candidate_id": f"{scenario_id}_{family}",
        "family": family,
        "selection_rule": "direct_physical_spacing",
        "delay_indices": "",
        "delay_grid_coordinates": encode(coordinates),
        "delay_ns": encode(delay_ns),
        "residues": "",
        "lifts": "",
        "canonical_key": "",
        "pair_sum_gap_indices": "",
        "fold_gap_indices": "",
        "grid_aligned": False,
        "pilot_rank": int(np.linalg.matrix_rank(phase, tol=1e-10)),
        "pilot_condition_number": float(singular[0] / singular[-1]),
        **moments,
    }


def lag_weights(covariance: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    diagonal = float(np.real(np.mean(np.diag(covariance))))
    physical = np.abs(covariance[0, 1:] / diagonal)
    lags = np.arange(1, K, dtype=np.float64)
    multiplicity = 2.0 * (K - lags)
    return multiplicity * physical**2, multiplicity * physical**4


def batch_moments(
    states: Sequence[Sequence[int]],
    covariance: np.ndarray,
    batch_size: int = 4_000,
) -> tuple[np.ndarray, np.ndarray]:
    weights2, weights4 = lag_weights(covariance)
    m2 = np.empty(len(states), dtype=np.float64)
    m4 = np.empty(len(states), dtype=np.float64)
    for start in range(0, len(states), int(batch_size)):
        stop = min(start + int(batch_size), len(states))
        counts = np.zeros((stop - start, K), dtype=np.float64)
        rows = np.arange(stop - start)[:, None]
        indices = np.asarray(states[start:stop], dtype=np.int64)
        counts[rows, indices] = 1.0
        factor = np.fft.fft(counts, axis=1)[:, 1:] / float(N_TX)
        magnitude = np.abs(factor)
        m2[start:stop] = magnitude**2 @ weights2
        m4[start:stop] = magnitude**4 @ weights4
    return m2, m4


def pareto_front(m2: np.ndarray, m4: np.ndarray) -> list[int]:
    order = np.lexsort((m4, m2))
    best_m4 = math.inf
    front: list[int] = []
    for index in order:
        value = float(m4[index])
        if value < best_m4 - 1e-12:
            front.append(int(index))
            best_m4 = value
    return front


def select_meff_indices(m2: np.ndarray, m4: np.ndarray, limit: int = 7) -> dict[int, set[str]]:
    selected: dict[int, set[str]] = {}

    def add(index: int, rule: str) -> None:
        selected.setdefault(int(index), set()).add(rule)

    minimum = int(np.argmin(m2))
    add(minimum, "minimum_m2")
    for band_db in MEFF_BANDS_DB:
        allowed = 10.0 * np.log10(np.maximum(m2, 1e-300) / float(m2[minimum])) <= band_db
        candidates = np.flatnonzero(allowed)
        add(int(candidates[np.argmin(m4[candidates])]), f"minimum_m4_within_{band_db:.2f}db_m2")
    add(int(np.argmin(m4)), "minimum_m4_ablation")
    front = pareto_front(m2, m4)
    if front:
        positions = np.unique(
            np.round(np.linspace(0, len(front) - 1, min(8, len(front)))).astype(int)
        )
        for position in positions:
            add(front[int(position)], "pareto_representative")
    ordered = sorted(
        selected,
        key=lambda index: (
            0 if "minimum_m2" in selected[index] else 1,
            float(m2[index]),
            float(m4[index]),
            index,
        ),
    )
    return {index: selected[index] for index in ordered[: int(limit)]}


def random_residue_set(
    rng: np.random.Generator,
    minimum_gap: int,
    pilot_period: int = PILOT_PERIOD,
) -> np.ndarray:
    if int(minimum_gap) >= 3:
        offset = int(rng.integers(0, 3))
        return (offset + 3 * np.arange(N_TX, dtype=np.int64)) % int(pilot_period)
    for _ in range(2_000):
        residues = np.sort(rng.choice(int(pilot_period), size=N_TX, replace=False))
        gaps = np.diff(np.r_[residues, residues[0] + int(pilot_period)])
        if int(np.min(gaps)) >= int(minimum_gap):
            return residues
    raise RuntimeError(f"Could not generate residues with gap {minimum_gap}.")


def geometry_states(
    seed: int,
    count: int,
    required_fold: int,
    pilot_period: int = PILOT_PERIOD,
) -> list[tuple[int, ...]]:
    rng = np.random.default_rng(int(seed))
    states: dict[tuple[int, ...], None] = {}
    fixed: Iterable[Sequence[int]] = (QC, SIDON)
    for value in fixed:
        key = canonical_delay_set(value, K)
        if fold_min_gap_indices(key, int(pilot_period)) >= min(required_fold, 3):
            states[key] = None
    while len(states) < int(count):
        residues = random_residue_set(
            rng,
            min(int(required_fold), 3),
            int(pilot_period),
        )
        lifts = rng.integers(0, K // int(pilot_period), size=N_TX)
        states[
            canonical_delay_set(residues + int(pilot_period) * lifts, K)
        ] = None
    return list(states)[: int(count)]


def geometry_selections(
    states: Sequence[Sequence[int]],
    required_pair: int,
    required_fold: int,
    limit: int = 3,
    pilot_period: int = PILOT_PERIOD,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[int, set[str]]]:
    pair = np.asarray([pair_sum_min_gap_indices(value, K) for value in states], dtype=np.int16)
    fold = np.asarray(
        [fold_min_gap_indices(value, int(pilot_period)) for value in states],
        dtype=np.int16,
    )
    hard = (pair >= int(required_pair)) & (fold >= int(required_fold))
    selected: dict[int, set[str]] = {}

    def add(index: int, rule: str) -> None:
        selected.setdefault(int(index), set()).add(rule)

    score = np.minimum(pair / max(1, required_pair), fold / max(1, required_fold))
    add(int(np.argmax(score)), "maximum_balanced_geometry")
    add(int(np.lexsort((-fold, -pair))[0]), "maximum_pair_gap")
    add(int(np.lexsort((-pair, -fold))[0]), "maximum_fold_gap")
    hard_indices = np.flatnonzero(hard)
    if hard_indices.size:
        add(int(hard_indices[0]), "hard_thick_sidon_representative")
    ordered = sorted(
        selected,
        key=lambda index: (-float(score[index]), -int(pair[index]), -int(fold[index]), index),
    )
    return pair, fold, hard, {index: selected[index] for index in ordered[: int(limit)]}


def write_large_candidate_csv(
    path: Path,
    scenario_id: str,
    family: str,
    states: Sequence[Sequence[int]],
    m2: np.ndarray,
    m4: np.ndarray,
    pair: np.ndarray | None = None,
    fold: np.ndarray | None = None,
    hard: np.ndarray | None = None,
    pilot_period: int = PILOT_PERIOD,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists()
    fields = [
        "scenario_id",
        "candidate_id",
        "family",
        "delay_indices",
        "delay_grid_coordinates",
        "delay_ns",
        "residues",
        "lifts",
        "canonical_key",
        "grid_aligned",
        "pilot_rank",
        "pilot_condition_number",
        "m2_eff",
        "m4_eff",
        "pair_sum_gap_indices",
        "fold_gap_indices",
        "hard_thick_sidon",
    ]
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        if not exists:
            writer.writeheader()
        for index, state in enumerate(states):
            key = canonical_delay_set(state, K)
            residues, lifts = residues_and_lifts(key, int(pilot_period))
            writer.writerow(
                {
                    "scenario_id": scenario_id,
                    "candidate_id": f"{scenario_id}_{family}_{index:06d}",
                    "family": family,
                    "delay_indices": encode(key),
                    "delay_grid_coordinates": encode(key),
                    "delay_ns": encode(delay_indices_to_ns(key, K, SCS_HZ)),
                    "residues": encode(residues),
                    "lifts": encode(lifts),
                    "canonical_key": encode(key),
                    "grid_aligned": True,
                    "pilot_rank": N_TX,
                    "pilot_condition_number": 1.0,
                    "m2_eff": float(m2[index]),
                    "m4_eff": float(m4[index]),
                    "pair_sum_gap_indices": "" if pair is None else int(pair[index]),
                    "fold_gap_indices": "" if fold is None else int(fold[index]),
                    "hard_thick_sidon": "" if hard is None else bool(hard[index]),
                }
            )


def grid_candidate_row(
    scenario_id: str,
    candidate_id: str,
    family: str,
    selection_rule: str,
    state: Sequence[int],
    m2: float,
    m4: float,
    required_pair: int,
    required_fold: int,
    pilot_spacing: int = PILOT_SPACING,
) -> dict[str, object]:
    geometry = candidate_geometry(state, int(pilot_spacing))
    return {
        "scenario_id": scenario_id,
        "candidate_id": candidate_id,
        "family": family,
        "selection_rule": selection_rule,
        **geometry,
        "m2_eff": float(m2),
        "m4_eff": float(m4),
        "hard_thick_sidon": bool(
            int(geometry["pair_sum_gap_indices"]) >= int(required_pair)
            and int(geometry["fold_gap_indices"]) >= int(required_fold)
        ),
    }


def run_search(output: Path, search_states: int, geometry_count: int) -> None:
    e1 = output / "e1_search"
    e2 = output / "e2_sidon_map"
    e1.mkdir(parents=True, exist_ok=True)
    e2.mkdir(parents=True, exist_ok=True)
    all_path = e1 / "e1_all_candidates.csv"
    sidon_path = e2 / "e2_sidon_candidates.csv"
    if all_path.exists() or sidon_path.exists():
        raise RuntimeError("Search output already exists; use a new run-id.")
    frozen: list[dict[str, object]] = []
    pareto_rows: list[dict[str, object]] = []
    support_rows: list[dict[str, object]] = []
    applicability_rows: list[dict[str, object]] = []
    for scenario_id, spread_ns in SCENARIOS:
        print(f"[search] {scenario_id}: covariance/support", flush=True)
        grid, covariance = physical_covariance(spread_ns)
        support = support_for_scenario(spread_ns)
        support_ns = float(support["support_width_s"]) * 1e9
        requirements = thick_sidon_requirements(
            support_ns,
            GRID_NS,
            K,
            N_TX,
            PILOT_PERIOD,
        )
        expected = EXPECTED_SUPPORT_NS[scenario_id]
        if not math.isclose(support_ns, expected, rel_tol=0.0, abs_tol=0.002):
            raise RuntimeError(
                f"{scenario_id} T_epsilon={support_ns:.6f} ns differs from "
                f"the frozen expectation {expected:.6f} ns."
            )
        required_pair = int(requirements["required_pair_gap_indices"])
        required_fold = int(requirements["required_fold_gap_indices"])
        support_rows.append(
            {
                "scenario_id": scenario_id,
                "tdl_profile": "A",
                "rms_delay_spread_ns": spread_ns,
                "epsilon": 0.01,
                "support_start_ns": float(support["support_start_s"]) * 1e9,
                "support_stop_ns": float(support["support_stop_s"]) * 1e9,
                "support_width_ns": support_ns,
                "contained_power_probability": float(support["contained_power"]),
                "grid_resolution_ns": GRID_NS,
                "support_over_grid_resolution": support_ns / GRID_NS,
                "pair_packing_ratio": N_TX * (N_TX + 1) * required_pair / (2.0 * K),
                "fold_packing_ratio": N_TX * required_fold / float(PILOT_PERIOD),
                **requirements,
            }
        )

        fixed_rows: list[dict[str, object]] = []
        for family, state in (("B0_QC", QC), ("S0_SIDON", SIDON)):
            moments = effective_moments_lag(state, covariance)
            fixed_rows.append(
                grid_candidate_row(
                    scenario_id,
                    f"{scenario_id}_{family}",
                    family,
                    "fixed_reference",
                    state,
                    moments["m2_eff"],
                    moments["m4_eff"],
                    required_pair,
                    required_fold,
                )
            )
        for family, state, selection_rule in SYSTEM_BASELINES:
            moments = effective_moments_lag(state, covariance)
            fixed_rows.append(
                grid_candidate_row(
                    scenario_id,
                    f"{scenario_id}_{family}",
                    family,
                    selection_rule,
                    state,
                    moments["m2_eff"],
                    moments["m4_eff"],
                    required_pair,
                    required_fold,
                )
            )
        rms_ns = continuous_arithmetic_delays_ns(spread_ns, N_TX)
        teps_ns = continuous_arithmetic_delays_ns(support_ns, N_TX)
        fixed_rows.extend(
            [
                continuous_candidate(scenario_id, "AP_RMS_T1", rms_ns, covariance, grid),
                continuous_candidate(scenario_id, "AP_TEPS_T1", teps_ns, covariance, grid),
            ]
        )

        arithmetic = []
        for state in arithmetic_delay_sets(N_TX, K):
            residues, _ = residues_and_lifts(state, PILOT_PERIOD)
            if len(set(residues)) == N_TX:
                arithmetic.append(state)
        ap_m2, ap_m4 = batch_moments(arithmetic, covariance)
        ap_selected = select_meff_indices(ap_m2, ap_m4, limit=2)
        ap_rows = [
            grid_candidate_row(
                scenario_id,
                f"{scenario_id}_AP_T2_{position:02d}",
                "AP_T2_CTRL",
                ";".join(sorted(rules)),
                arithmetic[index],
                ap_m2[index],
                ap_m4[index],
                required_pair,
                required_fold,
            )
            for position, (index, rules) in enumerate(ap_selected.items())
        ]

        print(f"[search] {scenario_id}: {geometry_count} geometry states", flush=True)
        geo_states = geometry_states(
            stable_seed(SEED, scenario_id, "geometry"),
            geometry_count,
            required_fold,
        )
        geo_pair, geo_fold, geo_hard, geo_selected = geometry_selections(
            geo_states,
            required_pair,
            required_fold,
            limit=3,
        )
        geo_m2, geo_m4 = batch_moments(geo_states, covariance)
        write_large_candidate_csv(
            sidon_path,
            scenario_id,
            "GEO_T1_CTRL",
            geo_states,
            geo_m2,
            geo_m4,
            geo_pair,
            geo_fold,
            geo_hard,
        )
        geo_rows = [
            grid_candidate_row(
                scenario_id,
                f"{scenario_id}_GEO_T1_{position:02d}",
                "GEO_T1_CTRL",
                ";".join(sorted(rules)),
                geo_states[index],
                geo_m2[index],
                geo_m4[index],
                required_pair,
                required_fold,
            )
            for position, (index, rules) in enumerate(geo_selected.items())
        ]
        applicability_rows.append(
            {
                "scenario_id": scenario_id,
                "hard_feasible_by_packing_bound": requirements[
                    "hard_feasible_by_packing_bound"
                ],
                "geometry_states": geometry_count,
                "hard_candidates_found": int(np.sum(geo_hard)),
                "maximum_pair_gap_indices": int(np.max(geo_pair)),
                "maximum_fold_gap_indices": int(np.max(geo_fold)),
                "s0_strict_sidon": pair_sum_min_gap_indices(SIDON, K) >= 1,
                "s0_hard_thick_sidon": bool(
                    pair_sum_min_gap_indices(SIDON, K) >= required_pair
                    and fold_min_gap_indices(SIDON, PILOT_PERIOD) >= required_fold
                ),
            }
        )

        initial = [QC, SIDON, *arithmetic]
        for spacing in (spread_ns, support_ns):
            projection = np.rint(
                delay_ns_to_grid_coordinates(
                    continuous_arithmetic_delays_ns(spacing, N_TX), K, SCS_HZ
                )
            ).astype(int) % K
            initial.append(tuple(int(value) for value in projection))
        initial.extend(row["delay_indices"] for row in [])
        initial.extend(geo_states[index] for index in geo_selected)
        print(f"[search] {scenario_id}: {search_states} T2 covariance states", flush=True)
        t2_states = unique_random_delay_sets(
            stable_seed(SEED, scenario_id, "t2"),
            search_states,
            initial=initial,
            n_tx=N_TX,
            period=K,
            pilot_period=PILOT_PERIOD,
        )
        t2_m2, t2_m4 = batch_moments(t2_states, covariance)
        write_large_candidate_csv(
            all_path,
            scenario_id,
            "MEFF_T2_CAND",
            t2_states,
            t2_m2,
            t2_m4,
        )
        t2_selected = select_meff_indices(t2_m2, t2_m4, limit=7)
        t2_rows = [
            grid_candidate_row(
                scenario_id,
                f"{scenario_id}_MEFF_T2_{position:02d}",
                "MEFF_T2_CAND",
                ";".join(sorted(rules)),
                t2_states[index],
                t2_m2[index],
                t2_m4[index],
                required_pair,
                required_fold,
            )
            for position, (index, rules) in enumerate(t2_selected.items())
        ]
        front = pareto_front(t2_m2, t2_m4)
        for index in front:
            pareto_rows.append(
                {
                    "scenario_id": scenario_id,
                    "delay_indices": encode(t2_states[index]),
                    "m2_eff": float(t2_m2[index]),
                    "m4_eff": float(t2_m4[index]),
                }
            )

        scenario_frozen = fixed_rows + ap_rows + geo_rows + t2_rows
        add_equivalent_delay_metrics(scenario_frozen, spread_ns)
        unique_rows: dict[tuple[str, str], dict[str, object]] = {}
        for row in scenario_frozen:
            identity = (
                str(row["family"]),
                str(row["delay_grid_coordinates"]),
            )
            unique_rows[identity] = row
        if len(unique_rows) > 18:
            raise RuntimeError(f"{scenario_id} froze more than 18 candidates.")
        frozen.extend(unique_rows.values())
        print(f"[search] {scenario_id}: froze {len(unique_rows)} candidates", flush=True)

    write_csv_rows(frozen, e1 / "e1_frozen_candidates.csv")
    write_csv_rows(pareto_rows, e1 / "e1_meff_pareto.csv")
    write_csv_rows(support_rows, e2 / "e2_support_geometry.csv")
    write_csv_rows(applicability_rows, e2 / "e2_applicability.csv")
    save_json(
        {
            "search_states_per_scenario": int(search_states),
            "geometry_states_per_scenario": int(geometry_count),
            "seed": SEED,
            "meff_constraint_bands_db": list(MEFF_BANDS_DB),
            "candidate_freeze_precedes_outage": True,
        },
        e1 / "search_summary.json",
    )


def iqam_table(output: Path) -> np.ndarray:
    cache = output / "validation" / "iqam_16qam_table.npz"
    return build_iqam_table(cache)


def run_add_system_baselines(output: Path) -> None:
    source = output / "e1_search" / "e1_frozen_candidates.csv"
    existing = read_csv(source)
    merged: list[dict[str, object]] = []
    for scenario_id, spread_ns in SCENARIOS:
        grid, covariance = physical_covariance(spread_ns)
        support_ns = float(
            support_for_scenario(spread_ns)["support_width_s"]
        ) * 1e9
        requirements = thick_sidon_requirements(
            support_ns,
            GRID_NS,
            K,
            N_TX,
            PILOT_PERIOD,
        )
        required_pair = int(requirements["required_pair_gap_indices"])
        required_fold = int(requirements["required_fold_gap_indices"])
        system_rows = []
        for family, state, selection_rule in SYSTEM_BASELINES:
            moments = effective_moments_lag(state, covariance)
            system_rows.append(
                grid_candidate_row(
                    scenario_id,
                    f"{scenario_id}_{family}",
                    family,
                    selection_rule,
                    state,
                    moments["m2_eff"],
                    moments["m4_eff"],
                    required_pair,
                    required_fold,
                )
            )
        system_ids = {str(row["candidate_id"]) for row in system_rows}
        scenario_rows: list[dict[str, object]] = [
            dict(row)
            for row in existing
            if row["scenario_id"] == scenario_id
            and row["candidate_id"] not in system_ids
        ]
        scenario_rows.extend(system_rows)
        add_equivalent_delay_metrics(scenario_rows, spread_ns)
        merged.extend(scenario_rows)

    temporary = source.with_suffix(source.suffix + ".tmp")
    write_csv_rows(merged, temporary)
    temporary.replace(source)
    save_json(
        {
            "status": "prepared",
            "baseline_families": [family for family, _, _ in SYSTEM_BASELINES],
            "AP_TU_NT_delay_indices": list(AP_TU_NT),
            "AP_TALIAS_NT_delay_indices": list(AP_TALIAS_NT),
            "equivalent_delay_metrics": {
                "rms": "minimum weighted squared circular distance over Tu",
                "support": "shortest circular interval containing at least 99% power",
                "period_s": 1.0 / SCS_HZ,
            },
            "candidate_count": len(merged),
            "scenario_counts": {
                scenario_id: sum(
                    row["scenario_id"] == scenario_id for row in merged
                )
                for scenario_id, _ in SCENARIOS
            },
        },
        output / "e1_search" / "system_baseline_addendum.json",
    )


def outage_target(curve: np.ndarray, snr_grid: np.ndarray, target: float, n_mc: int) -> float:
    probability = np.maximum(np.asarray(curve, dtype=np.float64), 0.5 / float(n_mc))
    for index in range(len(probability) - 1):
        if probability[index] >= float(target) > probability[index + 1]:
            first = math.log10(probability[index])
            second = math.log10(probability[index + 1])
            ratio = (first - math.log10(float(target))) / (first - second)
            return float(snr_grid[index] + ratio * (snr_grid[index + 1] - snr_grid[index]))
    return float("nan")


def threshold_from_mi(mi: np.ndarray, snr_grid: np.ndarray) -> np.ndarray:
    crossed = mi >= float(R_SE)
    first = np.argmax(crossed, axis=1)
    missing = ~np.any(crossed, axis=1)
    first = np.maximum(first, 1)
    rows = np.arange(len(mi))
    lo = first - 1
    hi = first
    denominator = mi[rows, hi] - mi[rows, lo]
    ratio = np.divide(
        float(R_SE) - mi[rows, lo],
        denominator,
        out=np.zeros_like(denominator),
        where=np.abs(denominator) > 1e-14,
    )
    threshold = snr_grid[lo] + ratio * (snr_grid[hi] - snr_grid[lo])
    threshold[missing] = snr_grid[-1] + 0.5
    return threshold.astype(np.float32)


def scenario_outage(
    candidate_rows: Sequence[dict[str, str]],
    spread_ns: float,
    table: np.ndarray,
    n_mc: int,
    seed: int,
    snr_grid: np.ndarray,
    chunk: int = 400,
) -> tuple[np.ndarray, np.ndarray]:
    delays_s, powers = tdl_pdp(spread_ns)
    physical_phase = np.exp(
        -2j
        * np.pi
        * np.arange(K, dtype=np.float64)[None, :]
        * SCS_HZ
        * delays_s[:, None]
    )
    coordinates = [decode(row["delay_grid_coordinates"]) for row in candidate_rows]
    artificial = np.asarray(
        [
            np.exp(
                -2j
                * np.pi
                * np.asarray(value, dtype=np.float64)[:, None]
                * np.arange(K, dtype=np.float64)[None, :]
                / float(K)
            )
            for value in coordinates
        ],
        dtype=np.complex128,
    )
    sqrt_powers = np.sqrt(np.asarray(powers, dtype=np.float64))
    tmat = build_shift_table(table, snr_grid)
    n_bins = int(round((HIST_DB_MAX - HIST_DB_MIN) / IQAM_DB_STEP)) + 1
    counts = np.zeros((len(candidate_rows), len(snr_grid)), dtype=np.int64)
    thresholds = np.empty((len(candidate_rows), int(n_mc)), dtype=np.float32)
    rng = np.random.default_rng(int(seed))
    for start in range(0, int(n_mc), int(chunk)):
        current = min(int(chunk), int(n_mc) - start)
        taps = (
            rng.standard_normal((current, N_TX, len(powers)))
            + 1j * rng.standard_normal((current, N_TX, len(powers)))
        ) / math.sqrt(2.0)
        taps *= sqrt_powers[None, None, :]
        branch = np.einsum("bnl,lk->bnk", taps, physical_phase, optimize=True)
        effective = np.einsum("bnk,cnk->bck", branch, artificial, optimize=True)
        for candidate in range(len(candidate_rows)):
            x_db = 10.0 * np.log10(
                np.abs(effective[:, candidate, :]) ** 2 / float(N_TX) + 1e-300
            )
            np.clip(
                x_db,
                HIST_DB_MIN,
                HIST_DB_MAX - IQAM_DB_STEP - 1e-6,
                out=x_db,
            )
            position = (x_db - HIST_DB_MIN) / IQAM_DB_STEP
            lower = np.floor(position).astype(np.int64)
            fraction = position - lower
            offsets = np.arange(current, dtype=np.int64)[:, None] * n_bins
            histogram = np.bincount(
                (offsets + lower).ravel(),
                weights=(1.0 - fraction).ravel(),
                minlength=current * n_bins,
            )
            histogram += np.bincount(
                (offsets + lower + 1).ravel(),
                weights=fraction.ravel(),
                minlength=current * n_bins,
            )
            mi = histogram.reshape(current, n_bins) @ tmat / float(K)
            counts[candidate] += np.sum(mi < float(R_SE), axis=0)
            thresholds[candidate, start : start + current] = threshold_from_mi(mi, snr_grid)
        if start == 0 or (start + current) % 20_000 == 0:
            print(f"  outage samples {start + current}/{n_mc}", flush=True)
    return counts / float(n_mc), thresholds


def bracket_indices(curve: np.ndarray, target: float) -> tuple[int, int]:
    for index in range(len(curve) - 1):
        if curve[index] >= target > curve[index + 1]:
            return index, index + 1
    raise RuntimeError(f"Outage target {target} is not bracketed.")


def interpolate_probabilities(
    first: np.ndarray,
    second: np.ndarray,
    first_snr: float,
    second_snr: float,
    target: float,
    n_mc: int,
) -> np.ndarray:
    denominator_count = float(n_mc) + 1.0
    p1 = (np.asarray(first, dtype=np.float64) * float(n_mc) + 0.5) / denominator_count
    p2 = (np.asarray(second, dtype=np.float64) * float(n_mc) + 0.5) / denominator_count
    p1 = np.maximum(p1, 0.5 / denominator_count)
    p2 = np.maximum(p2, 0.5 / denominator_count)
    p2 = np.minimum(p2, p1 - 0.5 / denominator_count)
    p2 = np.maximum(p2, 0.5 / denominator_count)
    bracket_epsilon = 0.25 / denominator_count
    p1 = np.maximum(p1, float(target) + bracket_epsilon)
    p2 = np.minimum(p2, max(0.5 / denominator_count, float(target) - bracket_epsilon))
    denominator = np.log10(p1) - np.log10(p2)
    ratio = np.divide(
        np.log10(p1) - math.log10(target),
        denominator,
        out=np.full_like(p1, np.nan, dtype=np.float64),
        where=np.abs(denominator) > 1e-15,
    )
    return first_snr + ratio * (second_snr - first_snr)


def paired_bootstrap_gain(
    candidate_threshold: np.ndarray,
    baseline_threshold: np.ndarray,
    candidate_curve: np.ndarray,
    baseline_curve: np.ndarray,
    snr_grid: np.ndarray,
    target: float,
    repeats: int,
    seed: int,
) -> np.ndarray:
    c0, c1 = bracket_indices(candidate_curve, target)
    b0, b1 = bracket_indices(baseline_curve, target)
    bits = (
        (candidate_threshold > snr_grid[c0]).astype(np.int8)
        + 2 * (candidate_threshold > snr_grid[c1]).astype(np.int8)
        + 4 * (baseline_threshold > snr_grid[b0]).astype(np.int8)
        + 8 * (baseline_threshold > snr_grid[b1]).astype(np.int8)
    )
    categories = np.bincount(bits, minlength=16)
    rng = np.random.default_rng(int(seed))
    sampled = rng.multinomial(
        len(bits),
        categories / float(len(bits)),
        size=int(repeats),
    )
    codes = np.arange(16)
    cp0 = sampled[:, (codes & 1) != 0].sum(axis=1) / float(len(bits))
    cp1 = sampled[:, (codes & 2) != 0].sum(axis=1) / float(len(bits))
    bp0 = sampled[:, (codes & 4) != 0].sum(axis=1) / float(len(bits))
    bp1 = sampled[:, (codes & 8) != 0].sum(axis=1) / float(len(bits))
    candidate = interpolate_probabilities(
        cp0, cp1, snr_grid[c0], snr_grid[c1], target, len(bits)
    )
    baseline = interpolate_probabilities(
        bp0, bp1, snr_grid[b0], snr_grid[b1], target, len(bits)
    )
    return baseline - candidate


def run_outage(output: Path, n_mc: int, repeats: int) -> None:
    source = output / "e1_search" / "e1_frozen_candidates.csv"
    rows = read_csv(source)
    stage = output / "e1_outage"
    stage.mkdir(parents=True, exist_ok=True)
    table = iqam_table(output)
    snr_grid = np.arange(0.0, 28.0 + 1e-9, 0.5)
    curve_rows: list[dict[str, object]] = []
    target_rows: list[dict[str, object]] = []
    bootstrap_rows: list[dict[str, object]] = []
    for scenario_id, spread_ns in SCENARIOS:
        selected = [row for row in rows if row["scenario_id"] == scenario_id]
        print(f"[outage] {scenario_id}: {len(selected)} candidates, n={n_mc}", flush=True)
        started = time.time()
        curves, thresholds = scenario_outage(
            selected,
            spread_ns,
            table,
            n_mc,
            stable_seed(SEED, scenario_id, "outage"),
            snr_grid,
        )
        np.savez(
            stage / f"{scenario_id}_paired_threshold_snr_db.npz",
            candidate_ids=np.asarray([row["candidate_id"] for row in selected]),
            threshold_snr_db=thresholds,
            snr_grid_db=snr_grid,
        )
        targets = {
            (index, target): outage_target(curves[index], snr_grid, target, n_mc)
            for index in range(len(selected))
            for target in (0.10, 0.01)
        }
        family_index = {
            row["family"]: index
            for index, row in enumerate(selected)
            if row["family"] in ("AP_RMS_T1", "AP_TEPS_T1")
        }
        for target in (0.10, 0.01):
            rms = family_index["AP_RMS_T1"]
            teps = family_index["AP_TEPS_T1"]
            baseline = rms if targets[(rms, target)] <= targets[(teps, target)] else teps
            baseline_row = selected[baseline]
            for index, row in enumerate(selected):
                point_gain = targets[(baseline, target)] - targets[(index, target)]
                gains = paired_bootstrap_gain(
                    thresholds[index],
                    thresholds[baseline],
                    curves[index],
                    curves[baseline],
                    snr_grid,
                    target,
                    repeats,
                    stable_seed(SEED, scenario_id, row["candidate_id"], target, "bootstrap"),
                )
                finite = gains[np.isfinite(gains)]
                if finite.size != int(repeats):
                    raise RuntimeError("Paired bootstrap produced non-finite replicates.")
                low, high = np.percentile(finite, [2.5, 97.5])
                bootstrap_rows.append(
                    {
                        "scenario_id": scenario_id,
                        "candidate_id": row["candidate_id"],
                        "target_outage_probability": target,
                        "baseline_candidate_id": baseline_row["candidate_id"],
                        "gain_snr_db": point_gain,
                        "gain_ci95_low_db": float(low),
                        "gain_ci95_high_db": float(high),
                        "bootstrap_repeats": int(repeats),
                        "paired": True,
                    }
                )
        lookup = {
            (row["candidate_id"], float(row["target_outage_probability"])): row
            for row in bootstrap_rows
            if row["scenario_id"] == scenario_id
        }
        for index, row in enumerate(selected):
            for snr_db, probability in zip(snr_grid, curves[index]):
                curve_rows.append(
                    {
                        "scenario_id": scenario_id,
                        "candidate_id": row["candidate_id"],
                        "family": row["family"],
                        "snr_db": float(snr_db),
                        "outage_probability": float(probability),
                        "n_mc": int(n_mc),
                    }
                )
            ten = lookup[(row["candidate_id"], 0.10)]
            one = lookup[(row["candidate_id"], 0.01)]
            target_rows.append(
                {
                    **row,
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
                    "elapsed_scenario_seconds": time.time() - started,
                }
            )
        write_csv_rows(curve_rows, stage / "e1_outage_curves.csv")
        write_csv_rows(target_rows, stage / "e1_outage_targets.csv")
        write_csv_rows(bootstrap_rows, stage / "e1_paired_bootstrap.csv")
    save_json(
        {
            "n_mc": int(n_mc),
            "bootstrap_repeats": int(repeats),
            "snr_grid_db": [float(value) for value in snr_grid],
            "initial_plan_grid_db": [0.0, 24.0, 0.5],
            "preauthorized_extension_db": 4.0,
            "spectral_efficiency_bit_per_re": float(R_SE),
            "common_random_numbers": "same TDL taps for all candidates within a scenario",
            "target_method": "log10 outage interpolation",
            "bootstrap_probability_smoothing": "Jeffreys 0.5-count with minimum adjacent 0.5-count decrease",
        },
        stage / "resolved_experiment.json",
    )
    delay_rows = [
        {
            "scenario_id": row["scenario_id"],
            "candidate_id": row["candidate_id"],
            "family": row["family"],
            "delay_grid_coordinates": row["delay_grid_coordinates"],
            "equivalent_circular_rms_delay_spread_ns": row[
                "equivalent_circular_rms_delay_spread_ns"
            ],
            "equivalent_99pct_circular_support_width_ns": row[
                "equivalent_99pct_circular_support_width_ns"
            ],
            "equivalent_delay_period_ns": row["equivalent_delay_period_ns"],
        }
        for row in target_rows
    ]
    write_csv_rows(
        delay_rows,
        output / "final" / "all_candidate_equivalent_delay_spread.csv",
    )
    delay_lines = [
        "# Plan-027 全部冻结基线/候选的等效时延扩展",
        "",
        "圆周 RMS 按有用符号周期 $T_u$ 上的最小加权平方圆周距离定义；"
        "99%支撑为包含至少99%复合PDP功率的最短圆周区间。",
        "",
        "| 场景 | candidate | family | delay-grid coordinates | 圆周RMS(ns) | 99%圆周支撑(ns) |",
        "|---|---|---|---|---:|---:|",
    ]
    for row in delay_rows:
        delay_lines.append(
            f"| {row['scenario_id']} | {row['candidate_id']} | "
            f"{row['family']} | `{row['delay_grid_coordinates']}` | "
            f"{float(row['equivalent_circular_rms_delay_spread_ns']):.3f} | "
            f"{float(row['equivalent_99pct_circular_support_width_ns']):.3f} |"
        )
    (output / "final" / "all_candidate_equivalent_delay_spread.md").write_text(
        "\n".join(delay_lines) + "\n",
        encoding="utf-8",
    )


def source_hash(function) -> str:
    return hashlib.sha256(inspect.getsource(function).encode("utf-8")).hexdigest()


def run_validate(output: Path) -> None:
    validation = output / "validation"
    validation.mkdir(parents=True, exist_ok=True)
    grid, covariance = physical_covariance(100.0)
    support = support_for_scenario(100.0)
    support_ns = float(support["support_width_s"]) * 1e9
    requirements = thick_sidon_requirements(support_ns, GRID_NS)
    direct_lag_rows = []
    for coordinates in (QC, SIDON, (0.0, 0.1, 0.7, 3.2, 8.1, 13.0, 21.3, 40.4)):
        direct = effective_moments_direct(coordinates, covariance)
        lag = effective_moments_lag(coordinates, covariance)
        direct_lag_rows.append(
            {
                "delay_grid_coordinates": encode(coordinates),
                "direct_m2_eff": direct["m2_eff"],
                "lag_m2_eff": lag["m2_eff"],
                "direct_m4_eff": direct["m4_eff"],
                "lag_m4_eff": lag["m4_eff"],
                "m2_abs_error": abs(direct["m2_eff"] - lag["m2_eff"]),
                "m4_abs_error": abs(direct["m4_eff"] - lag["m4_eff"]),
            }
        )
    replay_first = unique_random_delay_sets(SEED, 200)
    replay_second = unique_random_delay_sets(SEED, 200)
    density_grids = []
    common_data = None
    for spacing in (24, 12, 6):
        _, density_grid = resource_and_grid(spacing)
        density_grids.append((spacing, density_grid))
        coordinates = {
            tuple(value)
            for value in density_grid.data_coordinates.tolist()
        }
        common_data = coordinates if common_data is None else common_data & coordinates
    density_rows = []
    common_subcarriers = np.asarray(
        [coordinate[1] for coordinate in sorted(common_data or set())],
        dtype=np.int64,
    )
    for spacing, density_grid in density_grids:
        pilots = local_indices_for_subcarriers(
            density_grid, density_grid.pilot_subcarriers
        )
        phase = np.exp(-2j * np.pi * pilots[:, None] * np.asarray(SIDON)[None, :] / K)
        ce = FrequencyCEMetric(
            covariance,
            pilots,
            local_indices_for_subcarriers(density_grid, common_subcarriers),
            N_TX,
        ).evaluate(SIDON, 16.0)
        density_rows.append(
            {
                "dmrs_spacing_subcarriers": spacing,
                "pilot_count_per_symbol": density_grid.pilot_count,
                "pilot_re_count": density_grid.pilot_re_count,
                "fold_period_indices": K // spacing,
                "pilot_rank": int(np.linalg.matrix_rank(phase, tol=1e-10)),
                "data_re_count": density_grid.n_data_re,
                "common_data_re_count": len(common_subcarriers),
                "common_data_ce_nmse_db": float(ce["nmse_db"]),
            }
        )
    covariance_eigenvalues = np.linalg.eigvalsh(covariance)
    checks = {
        "K": grid.n_sc,
        "pilot_count_per_symbol": grid.pilot_count,
        "pilot_re_count": grid.pilot_re_count,
        "data_re_count": grid.n_data_re,
        "grid_resolution_ns": GRID_NS,
        "a100_support_ns": support_ns,
        "a100_requirements": requirements,
        "covariance_hermitian": bool(np.allclose(covariance, covariance.conj().T, atol=1e-11)),
        "covariance_minimum_eigenvalue": float(np.min(covariance_eigenvalues)),
        "covariance_diagonal_min": float(np.min(np.real(np.diag(covariance)))),
        "covariance_diagonal_max": float(np.max(np.real(np.diag(covariance)))),
        "direct_lag_rows": direct_lag_rows,
        "continuous_ap_rms_grid_coordinates": delay_ns_to_grid_coordinates(
            continuous_arithmetic_delays_ns(5.0)
        ).tolist(),
        "continuous_ap_teps_grid_coordinates": delay_ns_to_grid_coordinates(
            continuous_arithmetic_delays_ns(23.983)
        ).tolist(),
        "search_replay_equal": replay_first == replay_second,
        "search_unique_count": len(set(replay_first)),
        "density_rows": density_rows,
        "common_data_re_count": len(common_data or set()),
    }
    checks["passed"] = bool(
        checks["K"] == K
        and checks["pilot_count_per_symbol"] == 24
        and checks["pilot_re_count"] == 48
        and checks["data_re_count"] == 5712
        and math.isclose(support_ns, EXPECTED_SUPPORT_NS["A100"], abs_tol=0.002)
        and requirements["required_pair_gap_indices"] == 17
        and requirements["required_fold_gap_indices"] == 9
        and requirements["hard_feasible_by_packing_bound"] is False
        and checks["covariance_hermitian"]
        and checks["covariance_minimum_eigenvalue"] >= -1e-9
        and max(row["m2_abs_error"] for row in direct_lag_rows) < 1e-7
        and max(row["m4_abs_error"] for row in direct_lag_rows) < 1e-7
        and checks["search_replay_equal"]
        and checks["search_unique_count"] == 200
        and [row["pilot_count_per_symbol"] for row in density_rows] == [24, 48, 96]
        and [row["fold_period_indices"] for row in density_rows] == [24, 48, 96]
    )
    save_json(checks, validation / "preflight.json")
    write_csv_rows(direct_lag_rows, validation / "meff_direct_lag.csv")
    write_csv_rows(density_rows, validation / "dmrs_density.csv")
    related = {
        "cdd_metrics.py": file_sha256(ROOT / "cdd_lls" / "design" / "cdd_metrics.py"),
        "cdd_search.py": file_sha256(ROOT / "cdd_lls" / "design" / "cdd_search.py"),
        "estimators.py": file_sha256(ROOT / "cdd_lls" / "phy" / "estimators.py"),
        "run_track_b_pilot_scan.py": file_sha256(
            ROOT / "tools" / "run_track_b_pilot_scan.py"
        ),
    }
    save_json(
        {
            "source_experiment": "026",
            "source_evidence": "outputs/experiment026_cdd_design/20260724_main/validation/preflight.json",
            "source_hashes": related,
            "items": {
                "grid_phase_and_delay_conversion": "revalidated",
                "B0_S0_geometry": "revalidated",
                "TDL_A_covariance": "revalidated_with_A100",
                "matched_CE_formula": "revalidated_by_tests",
                "16QAM_MI_and_outage_normalization": "reused_unchanged_function_hashes",
            },
        },
        validation / "reuse_receipt.json",
    )
    receipt = environment_receipt()
    receipt.update(
        {
            "command": [sys.executable, *sys.argv],
            "platform": platform.platform(),
            "python": sys.version,
        }
    )
    save_json(receipt, validation / "environment.json")
    if not checks["passed"]:
        raise RuntimeError("Plan-027 Phase 0 validation failed.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage",
        choices=("validate", "search", "add-system-baselines", "outage"),
        required=True,
    )
    parser.add_argument("--run-id", default=DEFAULT_RUN_ID)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--search-states", type=int, default=SEARCH_STATES)
    parser.add_argument("--geometry-states", type=int, default=GEOMETRY_STATES)
    parser.add_argument("--n-mc", type=int, default=N_MC)
    parser.add_argument("--bootstrap-repeats", type=int, default=BOOTSTRAP_REPEATS)
    args = parser.parse_args()
    output = args.output_root / str(args.run_id)
    output.mkdir(parents=True, exist_ok=True)
    save_json(
        {
            "command": [sys.executable, *sys.argv],
            "plan_sha256": file_sha256(ROOT / "research" / "plan-027-有效矩-Sidon-DMRS.md"),
            "script_sha256": file_sha256(Path(__file__)),
        },
        output / f"commands_{args.stage}.json",
    )
    if args.stage == "validate":
        run_validate(output)
    elif args.stage == "search":
        run_search(output, int(args.search_states), int(args.geometry_states))
    elif args.stage == "add-system-baselines":
        run_add_system_baselines(output)
    else:
        run_outage(output, int(args.n_mc), int(args.bootstrap_repeats))


if __name__ == "__main__":
    main()
