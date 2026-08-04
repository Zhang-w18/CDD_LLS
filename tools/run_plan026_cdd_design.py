"""Plan-026 E1/E2 CDD search, deterministic metrics, and QAM/BICM outage."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ChannelConfig, ResourceConfig
from cdd_lls.design import (
    CDDKernel,
    FrequencyCEMetric,
    arithmetic_delay_sets,
    build_cdd_kernel,
    canonical_delay_set,
    delay_indices_to_ns,
    fold_min_gap_indices,
    jcdd_metrics,
    pair_sum_min_gap_indices,
    residues_and_lifts,
    tdl_effective_support,
)
from cdd_lls.design.cdd_search import random_unique_residue_delays
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
PILOT_PERIOD = 24
SCS_HZ = 30e3
GRID_NS = 1e9 / (K * SCS_HZ)
SEARCH_SNRS_DB = (14.0, 16.0, 18.0)
SEARCH_WEIGHTS = (0.0, 0.25, 0.5, 0.75, 1.0)
SEED = 20260726
DEFAULT_OUTPUT_ROOT = ROOT / "outputs" / "experiment026_cdd_design"
DEFAULT_RUN_ID = "20260724_main"
QC = (0, 9, 18, 27, 36, 45, 54, 63)
SIDON = (0, 1, 3, 7, 12, 20, 30, 65)
UNIFORM_RESIDUES = (0, 3, 6, 9, 12, 15, 18, 21)
E2_SCENARIOS = (
    ("E2-A10", "A", 10.0),
    ("E2-A20", "A", 20.0),
    ("E2-A30", "A", 30.0),
    ("E2-C5", "C", 5.0),
    ("E2-C10", "C", 10.0),
)


def resource_and_grid():
    resource = ResourceConfig(
        carrier_bandwidth_mhz=100.0,
        scs_khz=30,
        n_fft=4096,
        n_prbs=48,
        pdsch_n_symbols=10,
        dmrs_symbol_indices=[2, 7],
        dmrs_spacing_sc=24,
        dmrs_offset_sc=0,
        cyclic_prefix_length=288,
    )
    return resource, build_resource_grid(resource)


def channel_config(profile: str, delay_spread_ns: float) -> ChannelConfig:
    return ChannelConfig(
        backend="sionna_tdl",
        model="3gpp_tr38901_tdl",
        tdl_profile=str(profile).upper(),
        delay_spread_ns=float(delay_spread_ns),
        carrier_frequency_hz=3.5e9,
        ue_speed_kmh=0.0,
        num_sinusoids=20,
    )


def physical_covariance(profile: str, delay_spread_ns: float) -> tuple[object, np.ndarray]:
    _, grid = resource_and_grid()
    channel = channel_config(profile, delay_spread_ns)
    covariance = tdl_active_frequency_covariance(grid, channel)
    return grid, np.asarray(covariance, dtype=np.complex128)


def tdl_pdp(profile: str, delay_spread_ns: float) -> tuple[np.ndarray, np.ndarray]:
    from sionna.phy.channel.tr38901 import TDL

    tdl = TDL(
        model=str(profile).upper(),
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


def save_json(value: object, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding="utf-8")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def candidate_fields(delays: Sequence[int]) -> dict[str, object]:
    canonical = canonical_delay_set(delays, K)
    residues, lifts = residues_and_lifts(canonical, PILOT_PERIOD)
    return {
        "delay_indices": list(canonical),
        "delay_ns": [float(value) for value in delay_indices_to_ns(canonical, K, SCS_HZ)],
        "residues": list(residues),
        "lifts": list(lifts),
        "pair_sum_min_gap_indices": pair_sum_min_gap_indices(canonical, K),
        "pair_sum_min_gap_ns": pair_sum_min_gap_indices(canonical, K) * GRID_NS,
        "fold_min_gap_indices": fold_min_gap_indices(canonical, PILOT_PERIOD),
        "fold_min_gap_ns": fold_min_gap_indices(canonical, PILOT_PERIOD) * GRID_NS,
    }


def encode_list(values: Sequence[object]) -> str:
    return json.dumps(list(values), separators=(",", ":"))


def decode_indices(text: str) -> tuple[int, ...]:
    return tuple(int(value) for value in json.loads(text))


def iqam_and_kernels(output: Path) -> tuple[np.ndarray, np.ndarray, dict[float, CDDKernel]]:
    cache = output / "validation" / "iqam_16qam_table.npz"
    table = build_iqam_table(cache)
    grid_db = np.arange(IQAM_DB_MIN, IQAM_DB_MAX + 1e-9, IQAM_DB_STEP)
    kernels: dict[float, CDDKernel] = {}
    kernel_dir = output / "validation"
    for snr_db in SEARCH_SNRS_DB:
        path = kernel_dir / f"jcdd_kernel_{snr_db:g}db.npz"
        if path.exists():
            data = np.load(path)
            kernel = CDDKernel(
                float(snr_db),
                np.asarray(data["rho_grid"], dtype=np.float64),
                np.asarray(data["covariance_values"], dtype=np.float64),
            )
        else:
            kernel = build_cdd_kernel(snr_db, grid_db, table)
            np.savez(
                path,
                snr_db=float(snr_db),
                rho_grid=kernel.rho_grid,
                covariance_values=kernel.covariance_values,
            )
        kernels[float(snr_db)] = kernel
    return grid_db, table, kernels


class MetricEvaluator:
    def __init__(
        self,
        frequency_covariance: np.ndarray,
        kernels: dict[float, CDDKernel],
        grid,
    ) -> None:
        pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
        data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
        self.frequency_covariance = frequency_covariance
        self.kernels = kernels
        self.ce = FrequencyCEMetric(frequency_covariance, pilot_local, data_local, N_TX)
        self.cache16: dict[tuple[int, ...], dict[str, float | int]] = {}

    def at16(self, delays: Sequence[int]) -> dict[str, float | int]:
        key = canonical_delay_set(delays, K)
        if key not in self.cache16:
            row: dict[str, float | int] = {}
            row.update(jcdd_metrics(key, self.frequency_covariance, self.kernels[16.0]))
            ce = self.ce.evaluate(key, 16.0)
            row.update({f"ce_{name}": value for name, value in ce.items()})
            self.cache16[key] = row
        return self.cache16[key]

    def all_metrics(self, delays: Sequence[int]) -> dict[str, object]:
        key = canonical_delay_set(delays, K)
        row: dict[str, object] = {}
        for snr_db in SEARCH_SNRS_DB:
            jcdd = jcdd_metrics(key, self.frequency_covariance, self.kernels[snr_db])
            ce = self.ce.evaluate(key, snr_db)
            suffix = f"{int(snr_db)}db"
            row.update({f"{name}_{suffix}": value for name, value in jcdd.items()})
            row.update({f"ce_{name}_{suffix}": value for name, value in ce.items()})
        floor = self.ce.evaluate(key, 80.0)
        row.update({f"ce_{name}_80db": value for name, value in floor.items()})
        row.update(candidate_fields(key))
        return row


def fixed_residue_delays(lifts: Sequence[int]) -> tuple[int, ...]:
    return tuple(
        int(residue + PILOT_PERIOD * (int(lift) % (K // PILOT_PERIOD)))
        for residue, lift in zip(UNIFORM_RESIDUES, lifts)
    )


def qc_in_uniform_residue_order() -> tuple[int, ...]:
    return tuple(next(value for value in QC if value % PILOT_PERIOD == residue) for residue in UNIFORM_RESIDUES)


def pareto_indices(rows: Sequence[dict[str, object]], x_key: str, y_key: str) -> list[int]:
    out = []
    for index, row in enumerate(rows):
        x = float(row[x_key])
        y = float(row[y_key])
        dominated = any(
            other_index != index
            and float(other[x_key]) <= x
            and float(other[y_key]) <= y
            and (
                float(other[x_key]) < x
                or float(other[y_key]) < y
            )
            for other_index, other in enumerate(rows)
        )
        if not dominated:
            out.append(index)
    return out


def trace_row(
    family: str,
    weight: float,
    start_index: int,
    trace: dict[str, object],
) -> dict[str, object]:
    return {
        "family": family,
        "weight_jcdd": float(weight),
        "start_index": int(start_index),
        "sweep": int(trace["sweep"]),
        "coordinate": int(trace["coordinate"]),
        "objective": float(trace["score"]),
        "accepted": bool(trace["accepted"]),
        "delay_indices": encode_list(trace["delay_indices"]),
    }


def coordinate_run(
    initial: Sequence[int],
    neighbor_function,
    objective,
    max_sweeps: int = 20,
) -> tuple[tuple[int, ...], list[dict[str, object]]]:
    current = tuple(int(value) for value in initial)
    score = float(objective(current))
    trace = [{
        "sweep": 0,
        "coordinate": -1,
        "score": score,
        "accepted": True,
        "delay_indices": list(current),
    }]
    stagnant = 0
    for sweep in range(1, int(max_sweeps) + 1):
        improved = False
        for coordinate in range(N_TX):
            best = current
            best_score = score
            for proposal in neighbor_function(current, coordinate):
                candidate_score = float(objective(proposal))
                if candidate_score < best_score - 1e-12:
                    best = tuple(proposal)
                    best_score = candidate_score
            accepted = best != current
            if accepted:
                current = best
                score = best_score
                improved = True
            trace.append({
                "sweep": sweep,
                "coordinate": coordinate,
                "score": score,
                "accepted": accepted,
                "delay_indices": list(current),
            })
        stagnant = 0 if improved else stagnant + 1
        if stagnant >= 2:
            break
    return canonical_delay_set(current, K), trace


def deterministic_starts() -> list[tuple[int, ...]]:
    starts = [qc_in_uniform_residue_order()]
    bit_reverse = np.asarray([0, 4, 2, 6, 1, 5, 3, 7], dtype=np.int64)
    for scale in (0, 1, 2, 3):
        starts.append(fixed_residue_delays(scale * bit_reverse))
    rng = np.random.default_rng(SEED)
    while len(starts) < 32:
        starts.append(fixed_residue_delays(rng.integers(0, 24, size=N_TX)))
    return starts[:32]


def run_e1_search(output: Path) -> None:
    stage = output / "e1_search"
    stage.mkdir(parents=True, exist_ok=True)
    grid, rf = physical_covariance("A", 5.0)
    _, _, kernels = iqam_and_kernels(output)
    evaluator = MetricEvaluator(rf, kernels, grid)
    qc_metric = evaluator.at16(QC)
    qc_j = float(qc_metric["jcdd"])
    qc_ce = float(qc_metric["ce_nmse_db"])

    catalog: dict[tuple[int, ...], dict[str, object]] = {}
    traces: list[dict[str, object]] = []

    def add(delays: Sequence[int], family: str, source: str) -> None:
        key = canonical_delay_set(delays, K)
        row = catalog.setdefault(key, {
            "candidate_id": "",
            "family": family,
            "sources": set(),
            "delay_indices": key,
        })
        row["sources"].add(source)
        if family not in str(row["family"]):
            row["family"] = f"{row['family']}+{family}"

    add(QC, "B0", "fixed")
    add(SIDON, "S0", "fixed")

    print("[E1] evaluating arithmetic baseline", flush=True)
    arithmetic_rows = []
    for delays in arithmetic_delay_sets(N_TX, K):
        metric = evaluator.at16(delays)
        if int(metric["ce_pilot_rank"]) < N_TX:
            continue
        arithmetic_rows.append({
            "delays": delays,
            "jcdd": float(metric["jcdd"]),
            "ce_nmse_db": float(metric["ce_nmse_db"]),
        })
    arithmetic_front = pareto_indices(arithmetic_rows, "ce_nmse_db", "jcdd")
    for index in arithmetic_front:
        add(arithmetic_rows[index]["delays"], "B1", "arithmetic_pareto")

    starts = deterministic_starts()
    for weight in SEARCH_WEIGHTS:
        def objective(delays):
            metric = evaluator.at16(delays)
            if int(metric["ce_pilot_rank"]) < N_TX:
                return 1e6
            jcdd_delta = 10.0 * math.log10(max(float(metric["jcdd"]), 1e-300) / qc_j)
            ce_delta = float(metric["ce_nmse_db"]) - qc_ce
            return float(weight) * jcdd_delta + (1.0 - float(weight)) * ce_delta

        def fixed_neighbors(current, coordinate):
            residue = UNIFORM_RESIDUES[coordinate]
            for lift in range(24):
                proposal = list(current)
                proposal[coordinate] = residue + PILOT_PERIOD * lift
                yield tuple(proposal)

        for start_index, initial in enumerate(starts):
            print(f"[E1-U] weight={weight:g} start={start_index+1}/32", flush=True)
            final, trace = coordinate_run(initial, fixed_neighbors, objective)
            add(final, "E1_U", f"weight={weight:g};start={start_index}")
            traces.extend(trace_row("E1_U", weight, start_index, item) for item in trace)

    fixed_rows = [
        {"delays": key, **evaluator.at16(key)}
        for key, row in catalog.items()
        if "E1_U" in str(row["family"])
    ]
    fixed_front = [fixed_rows[index] for index in pareto_indices(fixed_rows, "ce_nmse_db", "jcdd")]
    fixed_front.sort(key=lambda row: (float(row["jcdd"]), float(row["ce_nmse_db"])))
    relaxed_starts = [tuple(row["delays"]) for row in fixed_front[:32]]

    for weight in SEARCH_WEIGHTS:
        def objective_relaxed(delays):
            key = canonical_delay_set(delays, K)
            residues, _ = residues_and_lifts(key, PILOT_PERIOD)
            if len(set(residues)) < N_TX:
                return 1e6
            metric = evaluator.at16(key)
            jcdd_delta = 10.0 * math.log10(max(float(metric["jcdd"]), 1e-300) / qc_j)
            ce_delta = float(metric["ce_nmse_db"]) - qc_ce
            return float(weight) * jcdd_delta + (1.0 - float(weight)) * ce_delta

        def relaxed_neighbors(current, coordinate):
            for shift in (-2, -1, 0, 1, 2):
                proposal = list(current)
                proposal[coordinate] = (proposal[coordinate] + shift) % K
                residues = [value % PILOT_PERIOD for value in proposal]
                if len(set(residues)) == N_TX:
                    yield tuple(proposal)

        for start_index, initial in enumerate(relaxed_starts):
            print(
                f"[E1-R] weight={weight:g} start={start_index+1}/{len(relaxed_starts)}",
                flush=True,
            )
            final, trace = coordinate_run(initial, relaxed_neighbors, objective_relaxed)
            add(final, "E1_R", f"weight={weight:g};start={start_index}")
            traces.extend(trace_row("E1_R", weight, start_index, item) for item in trace)

    rng = np.random.default_rng(SEED + 10)
    unrestricted_starts = [
        random_unique_residue_delays(rng, N_TX, K, PILOT_PERIOD)
        for _ in range(32)
    ]
    for weight in SEARCH_WEIGHTS:
        def objective_unrestricted(delays):
            key = canonical_delay_set(delays, K)
            residues, _ = residues_and_lifts(key, PILOT_PERIOD)
            if len(set(residues)) < N_TX:
                return 1e6
            metric = evaluator.at16(key)
            jcdd_delta = 10.0 * math.log10(max(float(metric["jcdd"]), 1e-300) / qc_j)
            ce_delta = float(metric["ce_nmse_db"]) - qc_ce
            return float(weight) * jcdd_delta + (1.0 - float(weight)) * ce_delta

        neighbor_seed = stable_seed(SEED, "E1_O", weight)
        neighbor_rng = np.random.default_rng(neighbor_seed)

        def unrestricted_neighbors(current, coordinate):
            residues_other = {
                int(value) % PILOT_PERIOD
                for index, value in enumerate(current)
                if index != coordinate
            }
            yield tuple(current)
            accepted = 0
            for value in neighbor_rng.permutation(K):
                if int(value) % PILOT_PERIOD in residues_other:
                    continue
                proposal = list(current)
                proposal[coordinate] = int(value)
                yield tuple(proposal)
                accepted += 1
                if accepted >= 23:
                    break

        for start_index, initial in enumerate(unrestricted_starts):
            print(f"[E1-O] weight={weight:g} start={start_index+1}/32", flush=True)
            final, trace = coordinate_run(initial, unrestricted_neighbors, objective_unrestricted)
            add(final, "E1_O", f"weight={weight:g};start={start_index}")
            traces.extend(trace_row("E1_O", weight, start_index, item) for item in trace)

    rows: list[dict[str, object]] = []
    for index, (key, entry) in enumerate(sorted(catalog.items())):
        metrics = evaluator.all_metrics(key)
        rows.append({
            "candidate_id": f"E1_{index:04d}",
            "family": entry["family"],
            "sources": ";".join(sorted(entry["sources"])),
            "delay_indices": encode_list(key),
            "delay_ns": encode_list(metrics.pop("delay_ns")),
            "residues": encode_list(metrics.pop("residues")),
            "lifts": encode_list(metrics.pop("lifts")),
            **metrics,
        })
    front = pareto_indices(rows, "ce_nmse_db_16db", "jcdd_16db")
    for index, row in enumerate(rows):
        row["pareto_16db"] = index in front
        row["jcdd_ratio_vs_qc_db_16db"] = 10.0 * math.log10(
            max(float(row["jcdd_16db"]), 1e-300) / qc_j
        )
        row["ce_nmse_delta_vs_qc_db_16db"] = float(row["ce_nmse_db_16db"]) - qc_ce
    write_csv_rows(rows, stage / "e1_metrics.csv")
    write_csv_rows(traces, stage / "search_trace.csv")
    write_csv_rows([rows[index] for index in front], stage / "e1_pareto.csv")
    freeze_e1_outage(rows, stage / "frozen_outage_candidates.csv")
    plot_e1_scatter(rows, stage / "e1_jcdd_ce_scatter_16db.png")
    save_json({
        "candidate_count": len(rows),
        "pareto_count": len(front),
        "metric_cache_count": len(evaluator.cache16),
        "search_weights": list(SEARCH_WEIGHTS),
        "starts_per_weight": 32,
        "max_sweeps": 20,
    }, stage / "search_summary.json")


def representative_rows(rows: Sequence[dict[str, object]], count: int) -> list[dict[str, object]]:
    if len(rows) <= int(count):
        return list(rows)
    ordered = sorted(rows, key=lambda row: float(row["ce_nmse_db_16db"]))
    indices = np.unique(np.round(np.linspace(0, len(ordered) - 1, int(count))).astype(int))
    return [ordered[int(index)] for index in indices]


def freeze_e1_outage(rows: Sequence[dict[str, object]], path: Path) -> None:
    selected: dict[str, dict[str, object]] = {}
    for family in ("B1", "E1_U", "E1_R", "E1_O"):
        family_rows = [
            row for row in rows
            if family in str(row["family"]) and bool(row["pareto_16db"])
        ]
        for row in representative_rows(family_rows, 7):
            selected[str(row["candidate_id"])] = row
    for row in rows:
        if "B0" in str(row["family"]) or "S0" in str(row["family"]):
            selected[str(row["candidate_id"])] = row
    frozen = list(selected.values())
    if len(frozen) > 32:
        frozen = frozen[:32]
    write_csv_rows(frozen, path)


def plot_e1_scatter(rows: Sequence[dict[str, object]], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(path.parent / "mplcache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    styles = {
        "B0": ("#f77f00", "o"),
        "B1": ("#5f5f5f", "^"),
        "S0": ("#111111", "s"),
        "E1_U": ("#277da1", "o"),
        "E1_R": ("#43aa8b", "D"),
        "E1_O": ("#9b5de5", "x"),
    }
    display_floor_db = -40.0

    def display_y(row: dict[str, object]) -> float:
        return max(display_floor_db, float(row["jcdd_ratio_vs_qc_db_16db"]))

    def on_pareto(row: dict[str, object]) -> bool:
        value = row["pareto_16db"]
        return value if isinstance(value, bool) else str(value).strip().lower() == "true"

    figure, axis = plt.subplots(figsize=(8.6, 6.2))
    for family, (color, marker) in styles.items():
        selected = [row for row in rows if family in str(row["family"])]
        if not selected:
            continue
        axis.scatter(
            [float(row["ce_nmse_db_16db"]) for row in selected],
            [display_y(row) for row in selected],
            s=28 if family.startswith("E1") else 70,
            alpha=0.72,
            color=color,
            marker=marker,
            label=family,
        )
    front = sorted(
        (row for row in rows if on_pareto(row)),
        key=lambda row: float(row["ce_nmse_db_16db"]),
    )
    if front:
        axis.plot(
            [float(row["ce_nmse_db_16db"]) for row in front],
            [display_y(row) for row in front],
            color="#d62828",
            linewidth=1.4,
            label="all-candidate Pareto",
        )
    annotated_ids = {"E1_0000", "E1_0028", "E1_0121", "E1_0161"}
    for row in rows:
        family = str(row["family"])
        if row["candidate_id"] not in annotated_ids:
            continue
        label = {
            "E1_0000": "B1",
            "E1_0028": "best joint",
            "E1_0121": "S0",
            "E1_0161": "B0",
        }[str(row["candidate_id"])]
        axis.annotate(
            label,
            (
                float(row["ce_nmse_db_16db"]),
                display_y(row),
            ),
            xytext=(3, 4),
            textcoords="offset points",
            fontsize=6.5,
            color="#333333",
        )
    axis.axhline(0.0, color="#888888", linestyle="--", linewidth=0.8)
    axis.axhline(
        display_floor_db,
        color="#aaaaaa",
        linestyle=":",
        linewidth=0.8,
    )
    axis.text(
        0.01,
        0.015,
        r"$J_{\rm CDD}=0$ numerical values displayed at -40 dB",
        transform=axis.transAxes,
        fontsize=8,
        color="#555555",
    )
    axis.set_ylim(display_floor_db - 2.0, 5.0)
    axis.set_xlabel("Matched data-RE CE NMSE at 16 dB (dB)")
    axis.set_ylabel(r"$10\log_{10}(J_{\rm CDD}/J_{\rm CDD,QC})$ at 16 dB (dB)")
    axis.set_title("Plan-026 E1: covariance-adapted CDD search")
    axis.grid(True, alpha=0.25)
    axis.legend(fontsize=8, ncol=2)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def random_residue_set(rng: np.random.Generator, minimum_gap: int) -> np.ndarray:
    if int(minimum_gap) >= 3:
        offset = int(rng.integers(0, 3))
        return (offset + 3 * np.arange(N_TX, dtype=np.int64)) % PILOT_PERIOD
    for _ in range(1000):
        residues = np.sort(rng.choice(PILOT_PERIOD, size=N_TX, replace=False))
        gaps = np.diff(np.r_[residues, residues[0] + PILOT_PERIOD])
        if int(np.min(gaps)) >= int(minimum_gap):
            return residues
    raise RuntimeError(f"Could not draw residues with minimum gap {minimum_gap}.")


def generate_thick_pool(
    support_ns: float,
    seed: int,
    samples: int = 200_000,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    required_pair = int(math.floor(2.0 * float(support_ns) / GRID_NS) + 1)
    required_fold = int(math.floor(float(support_ns) / GRID_NS) + 1)
    upper_pair = K // 36
    upper_fold = PILOT_PERIOD // N_TX
    if required_pair > upper_pair or required_fold > upper_fold:
        return [], {
            "required_pair_gap_indices": required_pair,
            "required_fold_gap_indices": required_fold,
            "pair_gap_packing_upper_bound": upper_pair,
            "fold_gap_packing_upper_bound": upper_fold,
            "hard_feasible_by_packing_bound": False,
        }
    rng = np.random.default_rng(int(seed))
    best: dict[tuple[int, ...], dict[str, object]] = {}
    seeds = [QC, SIDON]
    bit_reverse = np.asarray([0, 4, 2, 6, 1, 5, 3, 7], dtype=np.int64)
    for scale in (0, 1, 2, 3):
        seeds.append(fixed_residue_delays(scale * bit_reverse))
    for delays in seeds:
        key = canonical_delay_set(delays, K)
        best[key] = candidate_fields(key)
    for _ in range(int(samples)):
        residues = random_residue_set(rng, required_fold)
        lifts = rng.integers(0, K // PILOT_PERIOD, size=N_TX)
        delays = canonical_delay_set(residues + PILOT_PERIOD * lifts, K)
        pair_gap = pair_sum_min_gap_indices(delays, K)
        fold_gap = fold_min_gap_indices(delays, PILOT_PERIOD)
        if (
            pair_gap >= max(1, required_pair - 2)
            or (pair_gap >= required_pair and fold_gap >= required_fold)
        ):
            best[delays] = candidate_fields(delays)
    rows = list(best.values())
    for row in rows:
        row["hard_thick_sidon"] = bool(
            int(row["pair_sum_min_gap_indices"]) >= required_pair
            and int(row["fold_min_gap_indices"]) >= required_fold
        )
        row["required_pair_gap_indices"] = required_pair
        row["required_fold_gap_indices"] = required_fold
    rows.sort(
        key=lambda row: (
            not bool(row["hard_thick_sidon"]),
            -int(row["pair_sum_min_gap_indices"]),
            -int(row["fold_min_gap_indices"]),
            tuple(row["delay_indices"]),
        )
    )
    return rows[:500], {
        "required_pair_gap_indices": required_pair,
        "required_fold_gap_indices": required_fold,
        "pair_gap_packing_upper_bound": upper_pair,
        "fold_gap_packing_upper_bound": upper_fold,
        "hard_feasible_by_packing_bound": True,
        "random_samples": int(samples),
        "retained_candidates": min(len(rows), 500),
        "hard_candidate_count": sum(bool(row["hard_thick_sidon"]) for row in rows),
    }


def choose_t1_arithmetic(support_ns: float) -> tuple[int, ...]:
    required_fold = int(math.floor(float(support_ns) / GRID_NS) + 1)
    candidates = []
    for delays in arithmetic_delay_sets(N_TX, K):
        if fold_min_gap_indices(delays, PILOT_PERIOD) < required_fold:
            continue
        factor = np.abs(np.asarray([
            np.mean(np.exp(-2j * np.pi * lag * np.asarray(delays) / K))
            for lag in range(1, K)
        ]))
        candidates.append((
            float(np.sum(factor**2)),
            float(np.sum(factor**4)),
            delays,
        ))
    if not candidates:
        raise RuntimeError("No T1 arithmetic baseline satisfies the fold constraint.")
    return min(candidates)[2]


def run_e2(output: Path) -> None:
    stage = output / "e2_thick_sidon"
    stage.mkdir(parents=True, exist_ok=True)
    _, _, kernels = iqam_and_kernels(output)
    scenario_rows: list[dict[str, object]] = []
    frozen_rows: list[dict[str, object]] = []
    all_metric_rows: list[dict[str, object]] = []
    for scenario_id, profile, spread_ns in E2_SCENARIOS:
        print(f"[E2] {scenario_id}", flush=True)
        delays_s, powers = tdl_pdp(profile, spread_ns)
        support = tdl_effective_support(delays_s, powers, epsilon=0.01)
        support_ns = float(support["support_width_s"]) * 1e9
        pool, search_summary = generate_thick_pool(
            support_ns,
            stable_seed(SEED, scenario_id, "geometry"),
        )
        baseline = choose_t1_arithmetic(support_ns)
        pool.extend([
            {**candidate_fields(QC), "hard_thick_sidon": False},
            {**candidate_fields(SIDON), "hard_thick_sidon": False},
            {**candidate_fields(baseline), "hard_thick_sidon": False},
        ])
        unique = {canonical_delay_set(row["delay_indices"], K): row for row in pool}
        grid, rf = physical_covariance(profile, spread_ns)
        evaluator = MetricEvaluator(rf, kernels, grid)
        rows: list[dict[str, object]] = []
        baseline_key = canonical_delay_set(baseline, K)
        for index, (key, geometry) in enumerate(sorted(unique.items())):
            metrics = evaluator.all_metrics(key)
            family = (
                "B0" if key == canonical_delay_set(QC, K)
                else "S0" if key == canonical_delay_set(SIDON, K)
                else "B1_T1" if key == baseline_key
                else "THICK" if bool(geometry.get("hard_thick_sidon", False))
                else "GEOMETRY_NEAR"
            )
            rows.append({
                "scenario_id": scenario_id,
                "profile": profile,
                "delay_spread_ns": spread_ns,
                "candidate_id": f"{scenario_id}_{index:04d}",
                "family": family,
                "delay_indices": encode_list(key),
                "delay_ns": encode_list(metrics.pop("delay_ns")),
                "residues": encode_list(metrics.pop("residues")),
                "lifts": encode_list(metrics.pop("lifts")),
                "hard_thick_sidon": bool(geometry.get("hard_thick_sidon", False)),
                **metrics,
            })
        front = pareto_indices(rows, "ce_nmse_db_16db", "jcdd_16db")
        for index, row in enumerate(rows):
            row["pareto_16db"] = index in front
        all_metric_rows.extend(rows)
        selected: dict[str, dict[str, object]] = {}
        for row in rows:
            if row["family"] in ("B0", "S0", "B1_T1"):
                selected[str(row["candidate_id"])] = row
        thick_front = [
            rows[index] for index in front
            if rows[index]["family"] in ("THICK", "GEOMETRY_NEAR")
        ]
        for row in representative_rows(thick_front, 5):
            selected[str(row["candidate_id"])] = row
        frozen_rows.extend(selected.values())
        write_csv_rows(rows, stage / f"{scenario_id}_metrics.csv")
        scenario_rows.append({
            "scenario_id": scenario_id,
            "profile": profile,
            "delay_spread_ns": spread_ns,
            "epsilon": 0.01,
            "support_start_ns": float(support["support_start_s"]) * 1e9,
            "support_stop_ns": float(support["support_stop_s"]) * 1e9,
            "support_width_ns": support_ns,
            "contained_power": float(support["contained_power"]),
            **search_summary,
        })
    write_csv_rows(scenario_rows, stage / "e2_scenarios.csv")
    write_csv_rows(all_metric_rows, stage / "e2_metrics.csv")
    write_csv_rows(frozen_rows, stage / "frozen_outage_candidates.csv")
    plot_e2_map(scenario_rows, stage / "e2_potential_map.png")


def plot_e2_map(rows: Sequence[dict[str, object]], path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(path.parent / "mplcache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = [str(row["scenario_id"]) for row in rows]
    supports = [float(row["support_width_ns"]) for row in rows]
    required_pair = [int(row["required_pair_gap_indices"]) for row in rows]
    hard_count = [int(row.get("hard_candidate_count", 0)) for row in rows]
    figure, axes = plt.subplots(1, 2, figsize=(10.2, 4.4))
    axes[0].bar(labels, supports, color="#277da1")
    axes[0].set_ylabel(r"$T_{\epsilon}$ (ns), $\epsilon=1\%$")
    axes[0].tick_params(axis="x", rotation=35)
    axes[0].grid(True, axis="y", alpha=0.25)
    axes[1].bar(labels, required_pair, color="#43aa8b", label="required pair gap")
    axes[1].scatter(labels, hard_count, color="#d62828", label="hard candidates found")
    axes[1].set_ylabel("Grid steps / candidate count")
    axes[1].tick_params(axis="x", rotation=35)
    axes[1].legend(fontsize=8)
    axes[1].grid(True, axis="y", alpha=0.25)
    figure.suptitle("Plan-026 E2 thick-Sidon geometry map")
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def candidate_outage(
    delays: Sequence[int],
    profile: str,
    delay_spread_ns: float,
    iqam_table: np.ndarray,
    n_mc: int,
    seed: int,
    snr_grid_db: np.ndarray,
    chunk: int = 800,
) -> np.ndarray:
    delays_s, powers = tdl_pdp(profile, delay_spread_ns)
    physical_phase = np.exp(
        -2j
        * np.pi
        * np.arange(K, dtype=np.float64)[None, :]
        * SCS_HZ
        * delays_s[:, None]
    )
    artificial = np.exp(
        -2j
        * np.pi
        * np.arange(K, dtype=np.float64)[None, :]
        * np.asarray(delays, dtype=np.float64)[:, None]
        / float(K)
    )
    phase = artificial[:, None, :] * physical_phase[None, :, :]
    sqrt_powers = np.sqrt(np.asarray(powers, dtype=np.float64))
    n_bins = int(round((HIST_DB_MAX - HIST_DB_MIN) / IQAM_DB_STEP)) + 1
    tmat = build_shift_table(iqam_table, snr_grid_db)
    counts = np.zeros(len(snr_grid_db), dtype=np.int64)
    rng = np.random.default_rng(int(seed))
    for start in range(0, int(n_mc), int(chunk)):
        current = min(int(chunk), int(n_mc) - start)
        taps = (
            rng.standard_normal((current, N_TX, len(powers)))
            + 1j * rng.standard_normal((current, N_TX, len(powers)))
        ) / math.sqrt(2.0)
        taps *= sqrt_powers[None, None, :]
        effective = np.einsum("bnl,nlk->bk", taps, phase, optimize=True)
        x_db = 10.0 * np.log10(np.abs(effective) ** 2 / N_TX + 1e-300)
        np.clip(x_db, HIST_DB_MIN, HIST_DB_MAX - IQAM_DB_STEP - 1e-6, out=x_db)
        position = (x_db - HIST_DB_MIN) / IQAM_DB_STEP
        lower = np.floor(position).astype(np.int64)
        fraction = position - lower
        row_offsets = np.arange(current, dtype=np.int64)[:, None] * n_bins
        histogram = np.bincount(
            (row_offsets + lower).ravel(),
            weights=(1.0 - fraction).ravel(),
            minlength=current * n_bins,
        )
        histogram += np.bincount(
            (row_offsets + lower + 1).ravel(),
            weights=fraction.ravel(),
            minlength=current * n_bins,
        )
        mi = histogram.reshape(current, n_bins) @ tmat / float(K)
        counts += np.sum(mi < R_SE, axis=0)
    return counts / float(n_mc)


def outage_snr(curve: np.ndarray, snr_grid_db: np.ndarray, target: float, n_mc: int) -> float:
    probability = np.maximum(np.asarray(curve, dtype=np.float64), 0.5 / float(n_mc))
    for index in range(len(probability) - 1):
        if probability[index] >= float(target) > probability[index + 1]:
            first = math.log10(probability[index])
            second = math.log10(probability[index + 1])
            ratio = (first - math.log10(float(target))) / (first - second)
            return float(
                snr_grid_db[index] + ratio * (snr_grid_db[index + 1] - snr_grid_db[index])
            )
    return float("nan")


def run_outage_stage(
    output: Path,
    source_csv: Path,
    stage_name: str,
    n_mc: int,
) -> None:
    rows = read_csv(source_csv)
    stage = output / stage_name
    stage.mkdir(parents=True, exist_ok=True)
    _, table, _ = iqam_and_kernels(output)
    snr_grid = np.arange(0.0, 24.0 + 1e-9, 0.5)
    curve_rows: list[dict[str, object]] = []
    target_rows: list[dict[str, object]] = []
    for index, row in enumerate(rows):
        scenario_id = str(row.get("scenario_id", "E1-A5"))
        profile = str(row.get("profile", "A"))
        spread = float(row.get("delay_spread_ns", 5.0))
        delays = decode_indices(str(row["delay_indices"]))
        print(
            f"[outage] {index+1}/{len(rows)} {scenario_id} {row['candidate_id']} n={n_mc}",
            flush=True,
        )
        started = time.time()
        curve = candidate_outage(
            delays,
            profile,
            spread,
            table,
            int(n_mc),
            stable_seed(SEED, scenario_id, "outage"),
            snr_grid,
        )
        for snr_db, probability in zip(snr_grid, curve):
            curve_rows.append({
                "scenario_id": scenario_id,
                "profile": profile,
                "delay_spread_ns": spread,
                "candidate_id": row["candidate_id"],
                "family": row["family"],
                "snr_db": float(snr_db),
                "outage": float(probability),
                "n_mc": int(n_mc),
            })
        target_rows.append({
            "scenario_id": scenario_id,
            "profile": profile,
            "delay_spread_ns": spread,
            "candidate_id": row["candidate_id"],
            "family": row["family"],
            "delay_indices": row["delay_indices"],
            "delay_ns": row["delay_ns"],
            "outage10_snr_db": outage_snr(curve, snr_grid, 0.10, int(n_mc)),
            "outage1_snr_db": outage_snr(curve, snr_grid, 0.01, int(n_mc)),
            "elapsed_seconds": time.time() - started,
        })
        write_csv_rows(curve_rows, stage / "outage_curves.csv")
        write_csv_rows(target_rows, stage / "outage_targets.csv")
    save_json({
        "n_mc": int(n_mc),
        "snr_grid_db": [float(value) for value in snr_grid],
        "spectral_efficiency_bit_per_re": R_SE,
        "common_random_numbers": "same analytic TDL tap seed within each scenario",
        "source_candidates": str(source_csv),
    }, stage / "resolved_experiment.json")


def run_validate(output: Path) -> None:
    validation = output / "validation"
    validation.mkdir(parents=True, exist_ok=True)
    grid, rf = physical_covariance("A", 5.0)
    _, _, kernels = iqam_and_kernels(output)
    evaluator = MetricEvaluator(rf, kernels, grid)
    qc = evaluator.all_metrics(QC)
    sidon = evaluator.all_metrics(SIDON)
    checks = {
        "K": int(grid.n_sc),
        "pilot_count": int(grid.pilot_count),
        "data_re": int(grid.n_data_re),
        "grid_ns": GRID_NS,
        "qc": qc,
        "sidon": sidon,
        "sidon_pair_gap_expected_1": pair_sum_min_gap_indices(SIDON, K),
        "uniform_fold_gap_expected_3": fold_min_gap_indices(UNIFORM_RESIDUES, PILOT_PERIOD),
        "passed": bool(
            grid.n_sc == K
            and grid.pilot_count == PILOT_PERIOD
            and grid.n_data_re == 5712
            and pair_sum_min_gap_indices(SIDON, K) == 1
            and fold_min_gap_indices(UNIFORM_RESIDUES, PILOT_PERIOD) == 3
            and np.isfinite(float(qc["jcdd_16db"]))
            and np.isfinite(float(qc["ce_nmse_db_16db"]))
        ),
    }
    save_json(checks, validation / "preflight.json")
    receipt = environment_receipt()
    receipt["command"] = [sys.executable, *sys.argv]
    receipt["platform"] = platform.platform()
    save_json(receipt, validation / "environment.json")
    if not checks["passed"]:
        raise RuntimeError("Plan-026 design preflight failed.")


def analyze(output: Path) -> None:
    final = output / "final"
    final.mkdir(parents=True, exist_ok=True)
    summary: dict[str, object] = {}
    e1_targets = output / "e1_outage" / "outage_targets.csv"
    if e1_targets.exists():
        rows = read_csv(e1_targets)
        best = min(rows, key=lambda row: float(row["outage10_snr_db"]))
        summary["e1"] = {
            "candidate_count": len(rows),
            "best_10pct": best,
        }
    e2_targets = output / "e2_outage" / "outage_targets.csv"
    if e2_targets.exists():
        rows = read_csv(e2_targets)
        by_scenario = {}
        for scenario in sorted({row["scenario_id"] for row in rows}):
            selected = [row for row in rows if row["scenario_id"] == scenario]
            by_scenario[scenario] = min(
                selected,
                key=lambda row: float(row["outage10_snr_db"]),
            )
        summary["e2"] = {
            "candidate_count": len(rows),
            "best_by_scenario": by_scenario,
        }
    save_json(summary, final / "design_summary.json")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage",
        choices=("validate", "e1-search", "e1-outage", "e2", "e2-outage", "analyze"),
        required=True,
    )
    parser.add_argument("--run-id", default=DEFAULT_RUN_ID)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--n-mc", type=int, default=200_000)
    args = parser.parse_args()
    output = args.output_root / str(args.run_id)
    output.mkdir(parents=True, exist_ok=True)
    save_json(
        {
            "command": [sys.executable, *sys.argv],
            "sha256_plan": hashlib.sha256(
                (ROOT / "research" / "plan-026.md").read_bytes()
            ).hexdigest(),
        },
        output / f"commands_{args.stage}.json",
    )
    if args.stage == "validate":
        run_validate(output)
    elif args.stage == "e1-search":
        run_e1_search(output)
    elif args.stage == "e1-outage":
        run_outage_stage(
            output,
            output / "e1_search" / "frozen_outage_candidates.csv",
            "e1_outage",
            int(args.n_mc),
        )
    elif args.stage == "e2":
        run_e2(output)
    elif args.stage == "e2-outage":
        run_outage_stage(
            output,
            output / "e2_thick_sidon" / "frozen_outage_candidates.csv",
            "e2_outage",
            int(args.n_mc),
        )
    else:
        analyze(output)


if __name__ == "__main__":
    main()
