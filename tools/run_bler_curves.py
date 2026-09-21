"""Configuration-driven, resumable BLER-curve runner for static TDL-A studies.

The runner can start a new curve set or append non-overlapping trial intervals to
existing point-level CSV data.  Estimated-CSI and ideal-CSI receivers use the
same absolute-trial seed derivation, so samples remain paired wherever their
trial intervals overlap.
"""

from __future__ import annotations

import argparse
import faulthandler
import hashlib
import json
import math
import re
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

try:
    import yaml
except ModuleNotFoundError as exc:  # pragma: no cover
    raise RuntimeError("PyYAML is required by tools/run_bler_curves.py") from exc

from cdd_lls.core.config import ChannelConfig
from cdd_lls.core.mcs import build_tb_layout, get_mcs
from cdd_lls.phy.channel_tdl import generate_sionna_tdl_channel
from cdd_lls.phy.estimators import (
    build_frequency_rmmse_filter,
    build_prg_frequency_rmmse_filter,
    tdl_active_frequency_covariance,
)
from cdd_lls.phy.ldpc import SionnaLDPCAdapter
from cdd_lls.phy.precoding import build_prg_dft_precoder_batch, equivalent_channel
from cdd_lls.phy.resource_grid import build_resource_grid
from cdd_lls.phy.qam import qam_demapper_maxlog, qam_modulate
from cdd_lls.phy.resource_grid import local_indices_for_subcarriers
from tools import run_plan027_bler as bler027
from tools import run_plan028_csi_curves as curves028
from tools.run_plan025_delay_matched_tdl import (
    read_csv_rows,
    stable_seed,
    wilson,
    write_csv_rows,
)
from tools.run_v_design_piecewise_tradeoff import decode_same_tb_batch


SCHEMA = "bler-curve-runner-v1"
RECEIVERS = ("estimated", "ideal")
TRANSPARENT_MODE = "transparent_prg_dft"
TRANSPARENT_CDD_MODE = "transparent_cdd_physical_covariance"


def _resolve(path: str | Path, base: Path = ROOT) -> Path:
    value = Path(path)
    return value if value.is_absolute() else base / value


def _repo_relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return str(path.resolve())


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _is_transparent_scene(scene: dict) -> bool:
    return str(scene.get("mode", "")).lower() == TRANSPARENT_MODE


def _is_transparent_cdd_scene(scene: dict) -> bool:
    return str(scene.get("mode", "")).lower() == TRANSPARENT_CDD_MODE


def _scene_n_rx(scene: dict) -> int:
    """Return the configured receive-antenna count with strict validation."""
    value = scene.get("n_rx", 1)
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError("scene.n_rx must be a positive integer")
    return int(value)


def _transparent_candidates(scene: dict) -> list[dict]:
    rows = []
    for source in scene.get("transparent_prg_baselines", []):
        row = dict(source)
        row.setdefault("family", "TRANSPARENT_PRG_DFT")
        row.setdefault("scenario_id", str(scene["scenario_id"]))
        rows.append(row)
    identifiers = [str(row.get("candidate_id", "")) for row in rows]
    if not rows or any(not value for value in identifiers) or len(identifiers) != len(set(identifiers)):
        raise ValueError("transparent_prg_baselines require unique non-empty candidate_id values")
    return rows


def _transparent_cdd_candidates(scene: dict, manifest: dict) -> list[dict]:
    """Materialize receiver-variant IDs while preserving the source CDD transmitter."""
    source_rows = {
        str(row["candidate_id"]): row
        for row in bler027.scenario_candidates(manifest, str(scene["scenario_id"]))
    }
    rows = []
    for definition in scene.get("transparent_cdd_baselines", []):
        source_id = str(definition.get("source_candidate_id", ""))
        if source_id:
            if source_id not in source_rows:
                raise ValueError(f"Unknown transparent CDD source_candidate_id: {source_id!r}")
            row = dict(source_rows[source_id])
            row.update(dict(definition))
            row["source_family"] = str(source_rows[source_id]["family"])
            row["transmitter_definition"] = "source_manifest_candidate"
        else:
            coordinates = definition.get("delay_grid_coordinates")
            if not isinstance(coordinates, list) or len(coordinates) != 8:
                raise ValueError(
                    "Explicit transparent CDD candidates require eight delay_grid_coordinates"
                )
            values = [float(value) for value in coordinates]
            if not all(math.isfinite(value) for value in values):
                raise ValueError("Explicit transparent CDD delays must be finite")
            if any(right < left for left, right in zip(values, values[1:])):
                raise ValueError("Explicit transparent CDD delays must be nondecreasing")
            row = dict(definition)
            row["delay_grid_coordinates"] = values
            row["scenario_id"] = str(scene["scenario_id"])
            row["dmrs_spacing_subcarriers"] = int(manifest["dmrs_spacing_subcarriers"])
            row["source_candidate_id"] = ""
            row["source_family"] = "EXPLICIT_DELAY_GRID"
            row["transmitter_definition"] = "explicit_delay_grid_coordinates"
        row.setdefault("family", "TRANSPARENT_CDD_PHYSCOV")
        rows.append(row)
    identifiers = [str(row.get("candidate_id", "")) for row in rows]
    sources = [str(row.get("source_candidate_id", "")) for row in rows]
    if not rows or any(not value for value in identifiers) or len(identifiers) != len(set(identifiers)):
        raise ValueError("transparent_cdd_baselines require unique non-empty candidate_id values")
    nonempty_sources = [value for value in sources if value]
    if len(nonempty_sources) != len(set(nonempty_sources)):
        raise ValueError("transparent_cdd_baselines require unique source_candidate_id values")
    return rows


def _write_rows_atomic(rows: Sequence[dict], path: Path) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    write_csv_rows(rows, temporary)
    temporary.replace(path)


def load_runner_config(path: Path) -> dict:
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if payload.get("schema") != SCHEMA:
        raise ValueError(f"config.schema must be {SCHEMA!r}")
    if not isinstance(payload.get("scenes"), list) or not payload["scenes"]:
        raise ValueError("config.scenes must be a non-empty list")
    batch_size = int(payload.get("batch_size", 100))
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    payload["batch_size"] = batch_size
    payload["seed"] = int(payload.get("seed", bler027.SEED))
    payload["max_passes"] = int(payload.get("max_passes", 3))
    payload["output_dir"] = str(_resolve(payload["output_dir"]))
    payload["config_path"] = str(path.resolve())
    for scene in payload["scenes"]:
        if not isinstance(scene, dict):
            raise ValueError("config.scenes entries must be mappings")
        scene["n_rx"] = _scene_n_rx(scene)
    return payload


def _load_manifest(scene: dict) -> tuple[dict, str]:
    bundle = scene.get("manifest") or {}
    path = _resolve(bundle["path"])
    digest = _sha256(path)
    if bundle.get("sha256_path"):
        recorded = _resolve(bundle["sha256_path"]).read_text(encoding="utf-8").split()[0]
        if recorded.lower() != digest.lower():
            raise RuntimeError(f"Manifest SHA-256 mismatch: {path}")
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if bundle.get("approval_path"):
        approval_path = _resolve(bundle["approval_path"])
        approval = json.loads(approval_path.read_text(encoding="utf-8"))
        if approval.get("approved") is not True:
            raise RuntimeError(f"Manifest approval is not affirmative: {approval_path}")
        if str(approval.get("manifest_sha256", "")).lower() != digest.lower():
            raise RuntimeError(f"Manifest approval hash mismatch: {approval_path}")
    physical = manifest.get("physical_definition", {})
    if int(physical.get("phase_denominator", -1)) != 576:
        raise RuntimeError("BLER runner requires phase_denominator=576")
    if physical.get("delay_input_field") != "delay_grid_coordinates":
        raise RuntimeError("BLER runner requires delay_grid_coordinates")
    if int(manifest.get("dmrs_spacing_subcarriers", 0)) <= 0:
        raise RuntimeError("Manifest DMRS spacing must be positive")
    scenario_id = str(scene["scenario_id"])
    try:
        float(scenario_id.removeprefix("A"))
    except ValueError as exc:
        raise ValueError("scenario_id must encode TDL-A delay spread as A<number>") from exc
    if not _is_transparent_scene(scene) and not _is_transparent_cdd_scene(scene):
        candidates = bler027.scenario_candidates(manifest, scenario_id)
        requested = scene.get("candidate_ids", "all")
        if requested != "all":
            wanted = {str(value) for value in requested}
            actual = {str(row["candidate_id"]) for row in candidates}
            if not wanted <= actual:
                raise ValueError(f"Unknown candidate_ids: {sorted(wanted-actual)}")
    return manifest, digest


def _selected_candidates(scene: dict, manifest: dict) -> list[dict]:
    if _is_transparent_scene(scene):
        return _transparent_candidates(scene)
    if _is_transparent_cdd_scene(scene):
        return _transparent_cdd_candidates(scene, manifest)
    candidates = bler027.scenario_candidates(manifest, str(scene["scenario_id"]))
    requested = scene.get("candidate_ids", "all")
    if requested == "all":
        return candidates
    wanted = {str(value) for value in requested}
    return [row for row in candidates if str(row["candidate_id"]) in wanted]


def _base_rows(scene: dict, receiver: str, candidates: Sequence[dict]) -> list[dict]:
    receiver_cfg = scene["receivers"][receiver]
    base_csv = receiver_cfg.get("base_csv")
    if base_csv:
        rows = [dict(row) for row in read_csv_rows(_resolve(base_csv))]
        wanted = {str(row["candidate_id"]) for row in candidates}
        rows = [row for row in rows if str(row["candidate_id"]) in wanted]
        if not rows:
            raise RuntimeError(f"No selected rows in base_csv={base_csv}")
        expected_n_rx = _scene_n_rx(scene)
        for row in rows:
            row_n_rx = int(row.get("n_rx", 1))
            if row_n_rx != expected_n_rx:
                raise RuntimeError(
                    f"base_csv n_rx={row_n_rx} does not match scene.n_rx={expected_n_rx}"
                )
            row["n_rx"] = row_n_rx
        return rows
    snrs = [float(value) for value in receiver_cfg.get("snr_db", [])]
    if not snrs:
        raise ValueError(f"{receiver} requires base_csv or snr_db")
    rows = []
    for candidate in candidates:
        for snr_db in snrs:
            rows.append(
                {
                    "scenario_id": scene["scenario_id"],
                    "candidate_id": candidate["candidate_id"],
                    "family": candidate["family"],
                    "n_rx": _scene_n_rx(scene),
                    "snr_db": snr_db,
                    "trials": 0,
                    "tb_errors": 0,
                    "bler": 0.0,
                }
            )
    return rows


def _scene_output(config: dict, scene: dict) -> Path:
    return Path(config["output_dir"]) / str(scene["scenario_id"]).lower()


def _supplement_path(config: dict, scene: dict, receiver: str) -> Path:
    return _scene_output(config, scene) / receiver / "supplemental_points.csv"


def _flag_path(
    config: dict,
    scene: dict,
    receiver: str,
    snr_db: float,
    candidate_id: str,
    trial_start: int,
    trial_end: int,
) -> Path:
    snr_key = f"{float(snr_db):+.2f}".replace("+", "p").replace("-", "m").replace(".", "p")
    safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", candidate_id)
    name = f"{safe_id}_t{trial_start:06d}_{trial_end:06d}.npy"
    return _scene_output(config, scene) / receiver / "error_flags" / snr_key / name


def _trial_array_path(
    config: dict,
    scene: dict,
    category: str,
    snr_db: float,
    name: str,
    trial_start: int,
    trial_end: int,
) -> Path:
    snr_key = f"{float(snr_db):+.2f}".replace("+", "p").replace("-", "m").replace(".", "p")
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
    filename = f"{safe}_t{trial_start:06d}_{trial_end:06d}.npy"
    return _scene_output(config, scene) / category / snr_key / filename


def _read_supplements(config: dict, scene: dict, receiver: str) -> list[dict]:
    path = _supplement_path(config, scene, receiver)
    rows = [dict(row) for row in read_csv_rows(path)] if path.exists() else []
    for source in scene.get("additional_supplement_csvs", {}).get(receiver, []):
        source_path = _resolve(source)
        if source_path.exists():
            rows.extend(dict(row) for row in read_csv_rows(source_path))
    return rows


def _filter_run_snr_tasks(scene: dict, tasks: Sequence[dict]) -> list[dict]:
    configured = scene.get("run_snr_db")
    if configured is None:
        return list(tasks)
    allowed = {float(value) for value in configured}
    return [task for task in tasks if float(task["snr_db"]) in allowed]


def _validate_intervals(base_rows: Sequence[dict], supplement_rows: Sequence[dict]) -> None:
    base = {
        (str(row["candidate_id"]), float(row["snr_db"])): int(row["trials"])
        for row in base_rows
    }
    grouped: dict[tuple[str, float], list[tuple[int, int]]] = defaultdict(list)
    for row in supplement_rows:
        key = (str(row["candidate_id"]), float(row["snr_db"]))
        grouped[key].append((int(row["trial_start"]), int(row["trial_end"])))
    for key, intervals in grouped.items():
        expected = base.get(key, 0) + 1
        for start, end in sorted(intervals):
            if start != expected or end < start:
                raise RuntimeError(
                    f"Non-contiguous or overlapping supplemental interval for {key}: "
                    f"expected start {expected}, got [{start}, {end}]"
                )
            expected = end + 1


def merge_rows(base_rows: Sequence[dict], supplement_rows: Sequence[dict]) -> list[dict]:
    _validate_intervals(base_rows, supplement_rows)
    additions: dict[tuple[str, float], list[dict]] = defaultdict(list)
    for row in supplement_rows:
        additions[(str(row["candidate_id"]), float(row["snr_db"]))].append(row)
    merged = []
    for source in base_rows:
        row = dict(source)
        key = (str(row["candidate_id"]), float(row["snr_db"]))
        extra = additions.get(key, [])
        base_trials = int(row["trials"])
        base_errors = int(row["tb_errors"])
        trials = base_trials + sum(int(item["trials"]) for item in extra)
        errors = base_errors + sum(int(item["tb_errors"]) for item in extra)
        low, high = wilson(errors, trials)
        row.update(
            {
                "base_trials": base_trials,
                "supplemental_trials": trials - base_trials,
                "trials": trials,
                "tb_errors": errors,
                "bler": errors / trials if trials > 0 else 0.0,
                "bler_wilson95_lo": low,
                "bler_wilson95_hi": high,
            }
        )
        has_nmse = "ce_nmse_mean" in row or any(
            str(item.get("ce_nmse_sum", "")).strip() for item in extra
        )
        if has_nmse:
            total_nmse = float(row.get("ce_nmse_mean", 0.0)) * base_trials
            total_nmse += sum(float(item.get("ce_nmse_sum", 0.0)) for item in extra)
            total_nmse_sq = float(row.get("ce_nmse_sumsq", 0.0))
            if not total_nmse_sq and base_trials > 0 and "ce_nmse_std" in row:
                base_mean = float(row["ce_nmse_mean"])
                base_std = float(row["ce_nmse_std"])
                total_nmse_sq = (base_std**2 * max(base_trials - 1, 0)) + base_trials * base_mean**2
            total_nmse_sq += sum(float(item.get("ce_nmse_sumsq", 0.0)) for item in extra)
            mean_nmse = total_nmse / trials
            variance = max((total_nmse_sq - trials * mean_nmse**2) / max(trials - 1, 1), 0.0)
            standard_error = math.sqrt(variance / trials)
            row["ce_nmse_sum"] = total_nmse
            row["ce_nmse_sumsq"] = total_nmse_sq
            row["ce_nmse_mean"] = mean_nmse
            row["ce_nmse_std"] = math.sqrt(variance)
            row["ce_nmse_standard_error"] = standard_error
            row["ce_nmse_ci95_lo"] = max(mean_nmse - 1.96 * standard_error, 0.0)
            row["ce_nmse_ci95_hi"] = mean_nmse + 1.96 * standard_error
            row["ce_nmse_mean_db"] = 10.0 * math.log10(max(mean_nmse, 1e-300))
        merged.append(row)
    merged.sort(key=lambda row: (str(row["candidate_id"]), float(row["snr_db"])))
    return merged


def _plot_keys(rows: Sequence[dict], stop_below_bler: float | None) -> set[tuple[str, float]]:
    if stop_below_bler is None:
        return {(str(row["candidate_id"]), float(row["snr_db"])) for row in rows}
    output = set()
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[str(row["candidate_id"])].append(row)
    for candidate_id, values in grouped.items():
        for row in sorted(values, key=lambda item: float(item["snr_db"])):
            output.add((candidate_id, float(row["snr_db"])))
            if float(row["bler"]) < float(stop_below_bler):
                break
    return output


def _candidate_policy(policy: dict, candidate_id: str) -> dict:
    """Return the receiver policy after applying an optional candidate override."""
    effective = {
        key: value for key, value in policy.items() if key != "candidate_overrides"
    }
    overrides = policy.get("candidate_overrides", {})
    if overrides:
        override = overrides.get(candidate_id, {})
        if not isinstance(override, dict):
            raise ValueError(f"candidate_overrides[{candidate_id!r}] must be a mapping")
        effective.update(override)
    return effective


def _policy_plot_keys(rows: Sequence[dict], policy: dict) -> set[tuple[str, float]]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[str(row["candidate_id"])].append(row)
    output = set()
    for candidate_id, values in grouped.items():
        effective = _candidate_policy(policy, candidate_id)
        stop = effective.get("stop_below_bler")
        output.update(_plot_keys(values, None if stop is None else float(stop)))
    return output


def build_tasks(
    rows: Sequence[dict], policy: dict, batch_size: int
) -> list[dict]:
    selected = _policy_plot_keys(rows, policy)
    tasks = []
    for row in rows:
        key = (str(row["candidate_id"]), float(row["snr_db"]))
        if key not in selected:
            continue
        effective = _candidate_policy(policy, key[0])
        target_errors = int(effective.get("target_errors", 0))
        maximum = int(effective["max_total_trials"])
        minimum = int(effective.get("min_total_trials", 0))
        fixed = effective.get("fixed_total_trials")
        current = int(row["trials"])
        errors = int(row["tb_errors"])
        if current % batch_size:
            raise ValueError(f"Existing trials must align to batch_size for {key}")
        if fixed is not None:
            target = int(fixed)
        elif target_errors > 0 and errors >= target_errors:
            target = current
        elif errors > 0:
            target = math.ceil((target_errors * current / errors) / batch_size) * batch_size
        else:
            target = maximum
        target = max(current, minimum, target)
        target = min(maximum, math.ceil(target / batch_size) * batch_size)
        max_interval = int(effective.get("max_interval_trials", 0))
        if max_interval > 0:
            target = min(target, current + max_interval)
        if target > current:
            tasks.append(
                {
                    "candidate_id": key[0],
                    "snr_db": key[1],
                    "trial_start": current + 1,
                    "trial_end": target,
                    "trials": target - current,
                }
            )
    return tasks


def transparent_prg_order(
    candidate: dict,
    base_seed: int,
    snr_db: float,
    absolute_trial_index: int,
) -> list[int]:
    """Return the reproducible PRG-to-DFT-vector order for one absolute trial."""
    prg_size_rb = int(candidate["prg_size_rb"])
    mapping = str(candidate["mapping"])
    if prg_size_rb == 6 and mapping == "cycle_all":
        return list(range(8))
    if prg_size_rb == 4 and mapping == "cycle8_random_tail4":
        derived = stable_seed(
            int(base_seed),
            str(candidate["candidate_id"]),
            float(snr_db),
            int(absolute_trial_index),
            "tail_prg_order",
        )
        rng = np.random.Generator(np.random.PCG64(derived))
        tail = rng.choice(8, size=4, replace=False).astype(int).tolist()
        return [*range(8), *tail]
    raise ValueError(
        f"Unsupported transparent PRG mapping {mapping!r} for prg_size_rb={prg_size_rb}"
    )


def build_paired_tasks(
    rows_by_receiver: dict[str, Sequence[dict]],
    policy: dict,
    batch_size: int,
) -> list[dict]:
    """Build common SNR trial intervals shared by all candidates and receivers."""
    by_key = {
        receiver: {
            (str(row["candidate_id"]), float(row["snr_db"])): row for row in rows
        }
        for receiver, rows in rows_by_receiver.items()
    }
    keys = set(next(iter(by_key.values())))
    if any(set(values) != keys for values in by_key.values()):
        raise RuntimeError("Paired receivers do not have identical candidate/SNR grids.")
    target_errors = int(policy["target_errors"])
    minimum = int(policy["min_total_trials"])
    maximum = int(policy["max_total_trials"])
    tasks = []
    for snr_db in sorted({key[1] for key in keys}):
        point_rows = [
            values[key]
            for values in by_key.values()
            for key in sorted(values)
            if math.isclose(key[1], snr_db, abs_tol=1e-12)
        ]
        current_values = {int(row["trials"]) for row in point_rows}
        if len(current_values) != 1:
            raise RuntimeError(f"Paired trial counts differ at SNR={snr_db:g} dB.")
        current = current_values.pop()
        if current % int(batch_size):
            raise ValueError("Existing paired trials must align to batch_size.")
        errors = [int(row["tb_errors"]) for row in point_rows]
        if current < minimum:
            target = minimum
        elif all(value >= target_errors for value in errors):
            target = current
        else:
            estimates = [
                maximum
                if value == 0
                else math.ceil(target_errors * current / value / batch_size) * batch_size
                for value in errors
            ]
            target = max([minimum, current, *estimates])
        target = min(maximum, math.ceil(target / batch_size) * batch_size)
        max_interval = int(policy.get("max_interval_trials", 0))
        if max_interval > 0:
            target = min(target, current + max_interval)
        if target > current:
            tasks.append(
                {
                    "snr_db": snr_db,
                    "trial_start": current + 1,
                    "trial_end": target,
                    "trials": target - current,
                }
            )
    return tasks


def _decode_flags(
    adapter: SionnaLDPCAdapter,
    equalized: np.ndarray,
    effective_noise: np.ndarray,
    qm: int,
    payload: Sequence[np.ndarray],
) -> np.ndarray:
    llrs = [
        qam_demapper_maxlog(equalized[index], effective_noise[index], qm)
        for index in range(len(equalized))
    ]
    decoded = decode_same_tb_batch(adapter, llrs, payload)
    return np.asarray([not item.tb_success for item in decoded], dtype=bool)


def _ce_trial_statistics(
    estimate_data: np.ndarray, true_data: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return trial-level CE energies and ratios summed over Rx and data REs."""
    estimate = np.asarray(estimate_data, dtype=np.complex128)
    truth = np.asarray(true_data, dtype=np.complex128)
    if estimate.shape != truth.shape or estimate.ndim != 3:
        raise ValueError("CE arrays must share shape [batch,n_rx,n_data_re]")
    error_energy = np.sum(np.abs(estimate - truth) ** 2, axis=(1, 2))
    true_energy = np.sum(np.abs(truth) ** 2, axis=(1, 2))
    ratios = error_energy / np.maximum(true_energy, 1e-30)
    return error_energy, true_energy, ratios


def _estimated_csi_mrc_equalize(
    estimate_data: np.ndarray, received: np.ndarray, noise_variance: float
) -> tuple[np.ndarray, np.ndarray]:
    """Apply scalar single-layer MRC across the configured Rx dimension."""
    estimate = np.asarray(estimate_data, dtype=np.complex128)
    samples = np.asarray(received, dtype=np.complex128)
    if estimate.shape != samples.shape or estimate.ndim != 3:
        raise ValueError("MRC arrays must share shape [batch,n_rx,n_data_re]")
    denominator = np.maximum(np.sum(np.abs(estimate) ** 2, axis=1), 1e-10)
    equalized = np.sum(np.conj(estimate) * samples, axis=1) / denominator
    return equalized, float(noise_variance) / denominator


def _build_transparent_scene(manifest: dict, scenario_id: str):
    spacing = int(manifest["dmrs_spacing_subcarriers"])
    if spacing != 6 or scenario_id != "A100":
        raise ValueError("Transparent PRG baselines are fixed to A100 comb-6.")
    resource = bler027.resource_config(spacing)
    channel = ChannelConfig(
        backend="sionna_tdl",
        model="3gpp_tr38901_tdl",
        tdl_profile="A",
        delay_spread_ns=100.0,
        carrier_frequency_hz=3.5e9,
        ue_speed_kmh=0.0,
        num_sinusoids=20,
    )
    grid = build_resource_grid(resource)
    covariance = tdl_active_frequency_covariance(grid, channel)
    return channel, grid, covariance


def _write_transparent_manifest(
    config: dict,
    scene: dict,
    source_manifest_sha256: str,
    candidates: Sequence[dict],
) -> tuple[Path, str]:
    output = _scene_output(config, scene)
    plan_section = int(scene.get("plan_section", 10))
    if plan_section == 14:
        status = "researcher_confirmed_plan028_section14"
        plan = "research/plan-028-comb6三类CSI-TDL-A300ns.md#14-a100-8tx4rx-首批-idealestimated-csi-bler-增补"
    elif plan_section == 15:
        status = "researcher_confirmed_plan028_section15"
        plan = "research/plan-028-comb6三类CSI-TDL-A300ns.md#15-a100-8tx1rx-第-36-节两条曲线的-snr-延伸"
    else:
        status = "researcher_confirmed_in_plan028_section10"
        plan = "research/plan-028-comb6三类CSI-TDL-A300ns.md#10-a100-透明-prg-precoder-cycling-基线"
    payload = {
        "schema": "plan028-transparent-prg-manifest-v1",
        "status": status,
        "plan": plan,
        "scenario_id": "A100",
        "source_manifest_sha256": source_manifest_sha256,
        "system": {
            "n_prbs": 48,
            "active_subcarriers": 576,
            "subcarrier_spacing_hz": 30000.0,
            "n_tx": 8,
            "n_rx": _scene_n_rx(scene),
            "layers": 1,
            "tdl_profile": "A",
            "delay_spread_ns": 100.0,
            "ue_speed_kmh": 0.0,
            "dmrs_spacing_subcarriers": 6,
            "dmrs_symbol_indices": [2, 7],
            "data_re": 5568,
            "mcs_table": "nr_256qam",
            "mcs_index": 8,
        },
        "precoder": {
            "codebook": "8x8 spatial DFT",
            "coefficient_magnitude": 1.0,
            "vector_power": 8.0,
            "phase_sign": "negative",
            "prg_reference": "first_active_rb",
        },
        "randomization": {
            "generator": "numpy.random.Generator(PCG64)",
            "base_seed": int(config["seed"]),
            "seed_fields": [
                "base_seed",
                "candidate_id",
                "snr_db",
                "absolute_trial_index",
                "tail_prg_order",
            ],
        },
        "candidates": [dict(row) for row in candidates],
        "paired_policy": dict(scene["paired_policy"]),
    }
    path = output / "manifest" / "transparent_prg_manifest.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    digest = _sha256(path)
    (path.with_suffix(".sha256")).write_text(f"{digest}  {path.name}\n", encoding="utf-8")
    return path, digest


def _run_transparent_paired_tasks(
    config: dict,
    scene: dict,
    manifest: dict,
    candidates: Sequence[dict],
    tasks: Sequence[dict],
) -> None:
    if not tasks:
        return
    scenario_id = str(scene["scenario_id"])
    n_rx = _scene_n_rx(scene)
    receiver_ids = tuple(
        receiver for receiver in RECEIVERS if receiver in scene.get("receivers", {})
    )
    channel, grid, base_covariance = _build_transparent_scene(manifest, scenario_id)
    by_id = {str(row["candidate_id"]): row for row in candidates}
    candidate_ids = sorted(by_id)
    mcs = get_mcs("nr_256qam", 8, None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(tb.cb_k_values, tb.cb_e_values, num_iter=8, llr_clip=50.0)
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    stored = {receiver: _read_supplements(config, scene, receiver) for receiver in receiver_ids}
    for task in sorted(tasks, key=lambda item: float(item["snr_db"])):
        snr_db = float(task["snr_db"])
        trial_start = int(task["trial_start"])
        trial_end = int(task["trial_end"])
        trial_count = int(task["trials"])
        flags = {
            receiver: {candidate_id: np.zeros(trial_count, dtype=bool) for candidate_id in candidate_ids}
            for receiver in receiver_ids
        }
        ce_error = {candidate_id: np.zeros(trial_count, dtype=np.float64) for candidate_id in candidate_ids}
        ce_signal = {candidate_id: np.zeros(trial_count, dtype=np.float64) for candidate_id in candidate_ids}
        ce_nmse = {candidate_id: np.zeros(trial_count, dtype=np.float64) for candidate_id in candidate_ids}
        tail_orders = np.full((trial_count, 4), -1, dtype=np.int8)
        snr_linear = 10.0 ** (snr_db / 10.0)
        noise_variance = 8.0 / snr_linear
        ls_noise_variance = noise_variance / 2.0
        filters = {
            candidate_id: build_prg_frequency_rmmse_filter(
                8.0 * base_covariance,
                pilot_local,
                12 * int(by_id[candidate_id]["prg_size_rb"]),
                ls_noise_variance,
                diagonal_loading=1e-10,
            )
            for candidate_id in candidate_ids
        }
        diagnostic_path = _scene_output(config, scene) / "filter_diagnostics" / (
            f"snr_{snr_db:+.2f}".replace("+", "p").replace("-", "m").replace(".", "p")
            + ".json"
        )
        diagnostic_path.parent.mkdir(parents=True, exist_ok=True)
        diagnostic_path.write_text(
            json.dumps(
                {
                    "snr_db": snr_db,
                    "n_rx": n_rx,
                    "noise_variance": noise_variance,
                    "averaged_ls_noise_variance": ls_noise_variance,
                    "effective_covariance_power_factor": 8.0,
                    "candidates": {
                        candidate_id: list(filters[candidate_id].prg_diagnostics)
                        for candidate_id in candidate_ids
                    },
                },
                indent=2,
                ensure_ascii=False,
            )
            + "\n",
            encoding="utf-8",
        )
        started = time.time()
        for absolute_start in range(trial_start, trial_end + 1, int(config["batch_size"])):
            current = min(int(config["batch_size"]), trial_end - absolute_start + 1)
            offset = absolute_start - trial_start
            absolute_indices = range(absolute_start, absolute_start + current)
            payload_rng = np.random.default_rng(
                stable_seed(config["seed"], scenario_id, snr_db, absolute_start, "payload")
            )
            payload = [payload_rng.integers(0, 2, size=int(k), dtype=np.int8) for k in tb.cb_k_values]
            symbols = qam_modulate(np.concatenate(adapter.encode(payload)), int(mcs.qm))
            realization = generate_sionna_tdl_channel(
                grid,
                channel,
                n_tx=8,
                n_rx=n_rx,
                batch_size=current,
                seed=stable_seed(config["seed"], scenario_id, snr_db, absolute_start, "channel"),
            )
            if absolute_start == trial_start:
                print(
                    f"[TRACE {scenario_id} transparent paired] snr={snr_db:g} channel_ready "
                    f"elapsed={time.time()-started:.1f}s",
                    flush=True,
                )
            noise_rng = np.random.default_rng(
                stable_seed(config["seed"], scenario_id, snr_db, absolute_start, "noise")
            )
            averaged_ls_noise = math.sqrt(ls_noise_variance / 2.0) * (
                noise_rng.normal(size=(current, n_rx, grid.pilot_count))
                + 1j * noise_rng.normal(size=(current, n_rx, grid.pilot_count))
            )
            data_noise = math.sqrt(noise_variance / 2.0) * (
                noise_rng.normal(size=(current, n_rx, grid.n_data_re))
                + 1j * noise_rng.normal(size=(current, n_rx, grid.n_data_re))
            )
            decode_batches = []
            for candidate_id in candidate_ids:
                candidate = by_id[candidate_id]
                orders = np.asarray(
                    [
                        transparent_prg_order(candidate, config["seed"], snr_db, absolute_index)
                        for absolute_index in absolute_indices
                    ],
                    dtype=np.int64,
                )
                if int(candidate["prg_size_rb"]) == 4:
                    tail_orders[offset : offset + current] = orders[:, 8:12]
                precoder = build_prg_dft_precoder_batch(
                    grid,
                    n_tx=8,
                    prg_size_rb=int(candidate["prg_size_rb"]),
                    prg_vector_indices=orders,
                    normalize=False,
                )
                effective = equivalent_channel(realization.H, precoder)
                true_data = effective[:, :, grid.data_symbol_indices, data_local]
                received = true_data * symbols[None, None, :] + data_noise
                pilot0 = effective[:, :, int(grid.pilot_symbol_indices[0]), pilot_local]
                pilot1 = effective[:, :, int(grid.pilot_symbol_indices[-1]), pilot_local]
                estimate_full = filters[candidate_id].estimate_full_band(
                    0.5 * (pilot0 + pilot1) + averaged_ls_noise
                )
                estimate_data = estimate_full[:, :, data_local]
                error_energy, signal_energy, ratios = _ce_trial_statistics(
                    estimate_data, true_data
                )
                ce_error[candidate_id][offset : offset + current] = error_energy
                ce_signal[candidate_id][offset : offset + current] = signal_energy
                ce_nmse[candidate_id][offset : offset + current] = ratios
                estimated_equalized, estimated_noise = _estimated_csi_mrc_equalize(
                    estimate_data, received, noise_variance
                )
                if "estimated" in receiver_ids:
                    decode_batches.append(
                        ("estimated", candidate_id, estimated_equalized, estimated_noise)
                    )
                if "ideal" in receiver_ids:
                    ideal_equalized, ideal_noise = curves028.ideal_csi_equalize(
                        true_data, received, noise_variance
                    )
                    decode_batches.append(("ideal", candidate_id, ideal_equalized, ideal_noise))
                if absolute_start == trial_start:
                    print(
                        f"[TRACE {scenario_id} transparent paired] snr={snr_db:g} "
                        f"candidate={candidate_id} front_end_ready elapsed={time.time()-started:.1f}s",
                        flush=True,
                    )
            if absolute_start == trial_start:
                print(
                    f"[TRACE {scenario_id} transparent paired] snr={snr_db:g} decode_start "
                    f"elapsed={time.time()-started:.1f}s",
                    flush=True,
                )
            combined_flags = _decode_flags(
                adapter,
                np.concatenate([item[2] for item in decode_batches], axis=0),
                np.concatenate([item[3] for item in decode_batches], axis=0),
                int(mcs.qm),
                payload,
            )
            if absolute_start == trial_start:
                print(
                    f"[TRACE {scenario_id} transparent paired] snr={snr_db:g} decode_ready "
                    f"elapsed={time.time()-started:.1f}s",
                    flush=True,
                )
            cursor = 0
            for receiver, candidate_id, _, _ in decode_batches:
                flags[receiver][candidate_id][offset : offset + current] = combined_flags[
                    cursor : cursor + current
                ]
                cursor += current
            completed = absolute_start + current - 1
            if absolute_start == trial_start or completed % 500 < current or completed == trial_end:
                print(
                    f"[BLER {scenario_id} transparent paired] snr={snr_db:g} "
                    f"absolute_trials={completed}/{trial_end} elapsed={time.time()-started:.1f}s",
                    flush=True,
                )
        has_random_tail = any(
            int(by_id[candidate_id]["prg_size_rb"]) == 4 for candidate_id in candidate_ids
        )
        mapping_path = None
        mapping_sha256 = ""
        if has_random_tail:
            mapping_path = _trial_array_path(
                config,
                scene,
                "prg_tail_orders",
                snr_db,
                "A100_PRG_DFT8_4RB",
                trial_start,
                trial_end,
            )
            mapping_path.parent.mkdir(parents=True, exist_ok=True)
            np.save(mapping_path, tail_orders)
            mapping_sha256 = _sha256(mapping_path)
        for candidate_id in candidate_ids:
            ce_paths = {}
            for label, values in (
                ("error_energy", ce_error[candidate_id]),
                ("true_energy", ce_signal[candidate_id]),
                ("nmse", ce_nmse[candidate_id]),
            ):
                path = _trial_array_path(
                    config,
                    scene,
                    "ce_trial_stats",
                    snr_db,
                    f"{candidate_id}_{label}",
                    trial_start,
                    trial_end,
                )
                path.parent.mkdir(parents=True, exist_ok=True)
                np.save(path, values)
                ce_paths[label] = _repo_relative(path)
            for receiver in receiver_ids:
                flag_path = _flag_path(
                    config, scene, receiver, snr_db, candidate_id, trial_start, trial_end
                )
                flag_path.parent.mkdir(parents=True, exist_ok=True)
                np.save(flag_path, flags[receiver][candidate_id])
                errors = int(np.sum(flags[receiver][candidate_id]))
                low, high = wilson(errors, trial_count)
                row = {
                    "scenario_id": scenario_id,
                    "candidate_id": candidate_id,
                    "family": by_id[candidate_id]["family"],
                    "receiver": receiver,
                    "n_rx": n_rx,
                    "snr_db": snr_db,
                    "trial_start": trial_start,
                    "trial_end": trial_end,
                    "trials": trial_count,
                    "tb_errors": errors,
                    "bler": errors / trial_count,
                    "bler_wilson95_lo": low,
                    "bler_wilson95_hi": high,
                    "ce_nmse_sum": float(np.sum(ce_nmse[candidate_id])) if receiver == "estimated" else "",
                    "ce_nmse_sumsq": float(np.sum(ce_nmse[candidate_id] ** 2)) if receiver == "estimated" else "",
                    "ce_error_energy_path": ce_paths["error_energy"] if receiver == "estimated" else "",
                    "ce_true_energy_path": ce_paths["true_energy"] if receiver == "estimated" else "",
                    "ce_nmse_trial_path": ce_paths["nmse"] if receiver == "estimated" else "",
                    "prg_tail_order_path": _repo_relative(mapping_path) if mapping_path else "",
                    "prg_tail_order_sha256": mapping_sha256,
                    "seed": int(config["seed"]),
                    "error_flags": _repo_relative(flag_path),
                }
                stored[receiver].append(row)
        for receiver in receiver_ids:
            stored[receiver].sort(
                key=lambda row: (
                    str(row["candidate_id"]),
                    float(row["snr_db"]),
                    int(row["trial_start"]),
                )
            )
            _write_rows_atomic(stored[receiver], _supplement_path(config, scene, receiver))
        progress_path = _scene_output(config, scene) / "run_progress.jsonl"
        with progress_path.open("a", encoding="utf-8") as handle:
            handle.write(
                json.dumps(
                    {
                        "snr_db": snr_db,
                        "trial_start": trial_start,
                        "trial_end": trial_end,
                        "trials": trial_count,
                        "elapsed_seconds": time.time() - started,
                        "n_rx": n_rx,
                        "errors": {
                            receiver: {
                                candidate_id: int(np.sum(flags[receiver][candidate_id]))
                                for candidate_id in candidate_ids
                            }
                            for receiver in receiver_ids
                        },
                        "prg_tail_order_sha256": mapping_sha256,
                    },
                    ensure_ascii=False,
                )
                + "\n"
            )


def mismatched_frequency_lmmse_trace_nmse(
    true_covariance: np.ndarray,
    filter_weights: np.ndarray,
    pilot_local_indices: np.ndarray,
    data_local_indices: np.ndarray,
    observation_noise_variance: float,
) -> float:
    """Return trace NMSE for a fixed (possibly mismatched) frequency-LMMSE filter."""
    covariance = np.asarray(true_covariance, dtype=np.complex128)
    pilots = np.asarray(pilot_local_indices, dtype=np.int64)
    data = np.asarray(data_local_indices, dtype=np.int64)
    weights = np.asarray(filter_weights, dtype=np.complex128)[data]
    r_dp = covariance[np.ix_(data, pilots)]
    r_pp = covariance[np.ix_(pilots, pilots)]
    system = r_pp + float(observation_noise_variance) * np.eye(len(pilots))
    signal_trace = float(np.real(np.sum(covariance[data, data])))
    cross_trace = np.sum(weights * r_dp.conj())
    filtered_trace = np.einsum(
        "ip,pq,iq->", weights, system, weights.conj(), optimize=True
    )
    numerator = max(
        float(np.real(signal_trace - cross_trace - cross_trace.conjugate() + filtered_trace)),
        0.0,
    )
    if signal_trace <= 0.0:
        raise ValueError("True data covariance must have positive trace.")
    return numerator / signal_trace


def _transparent_cdd_receiver_filters(
    covariance_mode: str,
    candidate_ids: Sequence[str],
    true_covariances: dict[str, np.ndarray],
    physical_covariance: np.ndarray,
    pilot_local: np.ndarray,
    ls_noise_variance: float,
) -> tuple[dict[str, object], dict[str, str]]:
    """Build either the transparent shared filter or candidate-matched filters."""
    mode = str(covariance_mode).strip().lower()
    if mode == "physical":
        assumed = 8.0 * physical_covariance
        shared = build_frequency_rmmse_filter(
            assumed,
            pilot_local,
            ls_noise_variance,
            diagonal_loading=1e-10,
        )
        return (
            {candidate_id: shared for candidate_id in candidate_ids},
            {candidate_id: "8 * physical TDL-A covariance" for candidate_id in candidate_ids},
        )
    if mode == "matched_effective":
        filters = {
            candidate_id: build_frequency_rmmse_filter(
                true_covariances[candidate_id],
                pilot_local,
                ls_noise_variance,
                diagonal_loading=1e-10,
            )
            for candidate_id in candidate_ids
        }
        return filters, {
            candidate_id: "candidate-specific true effective CDD covariance"
            for candidate_id in candidate_ids
        }
    raise ValueError(f"Unsupported transparent CDD covariance_mode: {covariance_mode!r}")


def _write_transparent_cdd_manifest(
    config: dict,
    scene: dict,
    source_manifest_sha256: str,
    candidates: Sequence[dict],
) -> tuple[Path, str]:
    output = _scene_output(config, scene)
    has_explicit = any(
        str(row["transmitter_definition"]) == "explicit_delay_grid_coordinates"
        for row in candidates
    )
    receiver_modes = sorted(str(value) for value in scene.get("receivers", {}))
    plan_section = int(scene.get("plan_section", 12 if has_explicit else 11))
    mapping = []
    for row in candidates:
        item = {
            "candidate_id": str(row["candidate_id"]),
            "source_candidate_id": str(row.get("source_candidate_id", "")),
            "label": str(row.get("label", row["candidate_id"])),
        }
        if str(row["transmitter_definition"]) == "explicit_delay_grid_coordinates":
            item.update(
                {
                    "transmitter_definition": "explicit_delay_grid_coordinates",
                    "delay_grid_coordinates": [
                        float(value) for value in row["delay_grid_coordinates"]
                    ],
                }
            )
        mapping.append(item)
    if plan_section == 16:
        schema = "plan028-small-delay-cdd-matched-csi-manifest-v1"
        status = "researcher_requested_plan028_section16"
        plan = "research/plan-028-comb6三类CSI-TDL-A300ns.md#16-a100-小时延-cdd-的-ideal-csi-与-matched-covariance-estimated-csi"
    elif plan_section == 14:
        schema = "plan028-transparent-cdd-4rx-subset-manifest-v1"
        status = "researcher_confirmed_plan028_section14"
        plan = "research/plan-028-comb6三类CSI-TDL-A300ns.md#14-a100-8tx4rx-首批-idealestimated-csi-bler-增补"
    elif plan_section == 13:
        schema = "plan028-transparent-cdd-small-delay-diversity-manifest-v3"
        status = "researcher_requested_plan028_section13"
        plan = (
            "research/plan-028-comb6三类CSI-TDL-A300ns.md#13-a100-小时延透明-cdd分集与信道估计联合筛选"
        )
    else:
        schema = (
            "plan028-transparent-cdd-explicit-delay-manifest-v2"
            if has_explicit
            else "plan028-transparent-cdd-physical-covariance-manifest-v1"
        )
        status = (
            "researcher_requested_plan028_section12"
            if has_explicit
            else "researcher_requested_plan028_section11"
        )
        plan = (
            "research/plan-028-comb6三类CSI-TDL-A300ns.md#12-a100-小时延透明-cddce-nmse-低于--15-db-且-estimated-bler-闭合"
            if has_explicit
            else "research/plan-028-comb6三类CSI-TDL-A300ns.md#11-a100-透明-cddue-使用底层物理信道协方差"
        )
    covariance_mode = str(scene.get("covariance_mode", "physical")).strip().lower()
    matched_effective = covariance_mode == "matched_effective"
    receiver_payload = {
        "knows_cdd_delays": matched_effective,
        "covariance_mode": covariance_mode,
        "assumed_covariance": (
            "candidate-specific true effective CDD covariance"
            if matched_effective
            else "8 * R_phy(static TDL-A 100 ns)"
        ),
        "assumed_covariance_power_factor": None if matched_effective else 8.0,
        "uses_candidate_specific_covariance": matched_effective,
        "averaged_dmrs_count": 2,
        "averaged_ls_noise_variance": "4/SNR_linear",
        "receiver_modes": receiver_modes,
    }
    if receiver_modes == ["ideal"]:
        receiver_payload = {
            "receiver_modes": receiver_modes,
            "ideal_uses_true_effective_data_csi": True,
            "estimated_receiver_not_run": True,
        }
    payload = {
        "schema": schema,
        "status": status,
        "plan": plan,
        "scenario_id": "A100",
        "system": {
            "n_tx": 8,
            "n_rx": _scene_n_rx(scene),
            "layers": 1,
            "snr_definition": "per-Rx-branch average receive SNR",
        },
        "source_manifest_sha256": source_manifest_sha256,
        "transmitter": {
            "definition": (
                "explicit or source CDD candidate"
                if has_explicit
                else "unchanged source CDD candidate"
            ),
            "phase_denominator": 576,
            "coefficient_magnitude": 1.0,
            "total_power": 8.0,
            "noise_variance": "8/SNR_linear",
        },
        "receiver": receiver_payload,
        "source_candidate_mapping": mapping,
        "paired_policy": dict(scene["paired_policy"]),
        "seed": int(config["seed"]),
        "absolute_trial_seeded": True,
    }
    path = output / "manifest" / "transparent_cdd_manifest.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    digest = _sha256(path)
    (path.with_suffix(".sha256")).write_text(f"{digest}  {path.name}\n", encoding="utf-8")
    return path, digest


def _run_transparent_cdd_paired_tasks(
    config: dict,
    scene: dict,
    manifest: dict,
    candidates: Sequence[dict],
    tasks: Sequence[dict],
) -> None:
    """Run transparent CDD transmitters for one configured receiver mode."""
    if not tasks:
        return
    scenario_id = str(scene["scenario_id"])
    receiver_modes = sorted(str(value) for value in scene.get("receivers", {}))
    if len(receiver_modes) != 1 or receiver_modes[0] not in RECEIVERS:
        raise ValueError("Transparent CDD execution requires exactly one receiver mode.")
    receiver = receiver_modes[0]
    n_rx = _scene_n_rx(scene)
    by_id = {str(row["candidate_id"]): row for row in candidates}
    candidate_ids = [str(row["candidate_id"]) for row in candidates]
    channel, grid, precoders, true_covariances = bler027.build_scene(
        manifest, candidates, scenario_id
    )
    base_covariance = tdl_active_frequency_covariance(grid, channel)
    covariance_mode = str(scene.get("covariance_mode", "physical")).strip().lower()
    mcs = get_mcs("nr_256qam", 8, None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(tb.cb_k_values, tb.cb_e_values, num_iter=8, llr_clip=50.0)
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    stored = _read_supplements(config, scene, receiver)
    for task in sorted(tasks, key=lambda item: float(item["snr_db"])):
        snr_db = float(task["snr_db"])
        trial_start = int(task["trial_start"])
        trial_end = int(task["trial_end"])
        trial_count = int(task["trials"])
        flags = {
            candidate_id: np.zeros(trial_count, dtype=bool) for candidate_id in candidate_ids
        }
        ce_error = {
            candidate_id: np.zeros(trial_count, dtype=np.float64) for candidate_id in candidate_ids
        }
        ce_signal = {
            candidate_id: np.zeros(trial_count, dtype=np.float64) for candidate_id in candidate_ids
        }
        ce_nmse = {
            candidate_id: np.zeros(trial_count, dtype=np.float64) for candidate_id in candidate_ids
        }
        snr_linear = 10.0 ** (snr_db / 10.0)
        noise_variance = 8.0 / snr_linear
        ls_noise_variance = noise_variance / 2.0
        receiver_filters, covariance_descriptions = _transparent_cdd_receiver_filters(
            covariance_mode,
            candidate_ids,
            true_covariances,
            base_covariance,
            pilot_local,
            ls_noise_variance,
        )
        diagnostic = {
            "snr_db": snr_db,
            "n_rx": n_rx,
            "noise_variance": noise_variance,
            "averaged_ls_noise_variance": ls_noise_variance,
            "covariance_mode": covariance_mode,
            "knows_cdd_delays": covariance_mode == "matched_effective",
            "candidate_assumed_covariance": covariance_descriptions,
            "candidate_filter_condition_number": {
                candidate_id: receiver_filters[candidate_id].condition_number
                for candidate_id in candidate_ids
            },
            "candidate_filter_minimum_singular_value": {
                candidate_id: receiver_filters[candidate_id].minimum_singular_value
                for candidate_id in candidate_ids
            },
            "candidate_filter_numerical_jitter": {
                candidate_id: receiver_filters[candidate_id].numerical_jitter
                for candidate_id in candidate_ids
            },
            "candidate_analytic_trace_nmse": {
                candidate_id: mismatched_frequency_lmmse_trace_nmse(
                    true_covariances[candidate_id],
                    receiver_filters[candidate_id].weights,
                    pilot_local,
                    data_local,
                    ls_noise_variance,
                )
                for candidate_id in candidate_ids
            },
        }
        diagnostic_path = _scene_output(config, scene) / "filter_diagnostics" / (
            f"snr_{snr_db:+.2f}".replace("+", "p").replace("-", "m").replace(".", "p")
            + ".json"
        )
        diagnostic_path.parent.mkdir(parents=True, exist_ok=True)
        diagnostic_path.write_text(
            json.dumps(diagnostic, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        started = time.time()
        for absolute_start in range(trial_start, trial_end + 1, int(config["batch_size"])):
            current = min(int(config["batch_size"]), trial_end - absolute_start + 1)
            offset = absolute_start - trial_start
            payload_rng = np.random.default_rng(
                stable_seed(config["seed"], scenario_id, snr_db, absolute_start, "payload")
            )
            payload = [
                payload_rng.integers(0, 2, size=int(k), dtype=np.int8) for k in tb.cb_k_values
            ]
            symbols = qam_modulate(np.concatenate(adapter.encode(payload)), int(mcs.qm))
            realization = generate_sionna_tdl_channel(
                grid,
                channel,
                n_tx=8,
                n_rx=n_rx,
                batch_size=current,
                seed=stable_seed(config["seed"], scenario_id, snr_db, absolute_start, "channel"),
            )
            noise_rng = np.random.default_rng(
                stable_seed(config["seed"], scenario_id, snr_db, absolute_start, "noise")
            )
            averaged_ls_noise = math.sqrt(ls_noise_variance / 2.0) * (
                noise_rng.normal(size=(current, n_rx, grid.pilot_count))
                + 1j * noise_rng.normal(size=(current, n_rx, grid.pilot_count))
            )
            data_noise = math.sqrt(noise_variance / 2.0) * (
                noise_rng.normal(size=(current, n_rx, grid.n_data_re))
                + 1j * noise_rng.normal(size=(current, n_rx, grid.n_data_re))
            )
            decode_batches = []
            for candidate_id in candidate_ids:
                effective = equivalent_channel(
                    realization.H, precoders[candidate_id].C
                )
                true_data = effective[:, :, grid.data_symbol_indices, data_local]
                received = true_data * symbols[None, None, :] + data_noise
                if receiver == "estimated":
                    pilot0 = effective[:, :, int(grid.pilot_symbol_indices[0]), pilot_local]
                    pilot1 = effective[:, :, int(grid.pilot_symbol_indices[-1]), pilot_local]
                    estimate_full = receiver_filters[candidate_id].estimate_full_band(
                        0.5 * (pilot0 + pilot1) + averaged_ls_noise
                    )
                    estimate_data = estimate_full[:, :, data_local]
                    error_energy, signal_energy, ratios = _ce_trial_statistics(
                        estimate_data, true_data
                    )
                    ce_error[candidate_id][offset : offset + current] = error_energy
                    ce_signal[candidate_id][offset : offset + current] = signal_energy
                    ce_nmse[candidate_id][offset : offset + current] = ratios
                    equalized, effective_noise = _estimated_csi_mrc_equalize(
                        estimate_data, received, noise_variance
                    )
                else:
                    equalized, effective_noise = curves028.ideal_csi_equalize(
                        true_data, received, noise_variance
                    )
                decode_batches.append((candidate_id, equalized, effective_noise))
            for candidate_id, equalized, effective_noise in decode_batches:
                flags[candidate_id][offset : offset + current] = _decode_flags(
                    adapter,
                    equalized,
                    effective_noise,
                    int(mcs.qm),
                    payload,
                )
            completed = absolute_start + current - 1
            if absolute_start == trial_start or completed % 500 < current or completed == trial_end:
                print(
                    f"[BLER {scenario_id} transparent CDD {receiver}] snr={snr_db:g} "
                    f"absolute_trials={completed}/{trial_end} elapsed={time.time()-started:.1f}s",
                    flush=True,
                )
        for candidate_id in candidate_ids:
            ce_paths = {}
            if receiver == "estimated":
                for label, values in (
                    ("error_energy", ce_error[candidate_id]),
                    ("true_energy", ce_signal[candidate_id]),
                    ("nmse", ce_nmse[candidate_id]),
                ):
                    path = _trial_array_path(
                        config,
                        scene,
                        "ce_trial_stats",
                        snr_db,
                        f"{candidate_id}_{label}",
                        trial_start,
                        trial_end,
                    )
                    path.parent.mkdir(parents=True, exist_ok=True)
                    np.save(path, values)
                    ce_paths[label] = _repo_relative(path)
            flag_path = _flag_path(
                config, scene, receiver, snr_db, candidate_id, trial_start, trial_end
            )
            flag_path.parent.mkdir(parents=True, exist_ok=True)
            np.save(flag_path, flags[candidate_id])
            errors = int(np.sum(flags[candidate_id]))
            low, high = wilson(errors, trial_count)
            row = {
                    "scenario_id": scenario_id,
                    "candidate_id": candidate_id,
                    "source_candidate_id": by_id[candidate_id].get("source_candidate_id", ""),
                    "family": by_id[candidate_id]["family"],
                    "receiver": receiver,
                    "n_rx": n_rx,
                    "snr_db": snr_db,
                    "trial_start": trial_start,
                    "trial_end": trial_end,
                    "trials": trial_count,
                    "tb_errors": errors,
                    "bler": errors / trial_count,
                    "bler_wilson95_lo": low,
                    "bler_wilson95_hi": high,
                    "seed": int(config["seed"]),
                    "error_flags": _repo_relative(flag_path),
                }
            if receiver == "estimated":
                row.update(
                    {
                        "ce_nmse_sum": float(np.sum(ce_nmse[candidate_id])),
                        "ce_nmse_sumsq": float(np.sum(ce_nmse[candidate_id] ** 2)),
                        "ce_error_energy_path": ce_paths["error_energy"],
                        "ce_true_energy_path": ce_paths["true_energy"],
                        "ce_nmse_trial_path": ce_paths["nmse"],
                        "analytic_trace_nmse": diagnostic["candidate_analytic_trace_nmse"][candidate_id],
                    }
                )
            stored.append(row)
        stored.sort(
            key=lambda row: (
                str(row["candidate_id"]),
                float(row["snr_db"]),
                int(row["trial_start"]),
            )
        )
        _write_rows_atomic(stored, _supplement_path(config, scene, receiver))
        progress_path = _scene_output(config, scene) / "run_progress.jsonl"
        with progress_path.open("a", encoding="utf-8") as handle:
            handle.write(
                json.dumps(
                    {
                        "snr_db": snr_db,
                        "trial_start": trial_start,
                        "trial_end": trial_end,
                        "trials": trial_count,
                        "elapsed_seconds": time.time() - started,
                        "receiver": receiver,
                        "n_rx": n_rx,
                        "errors": {
                            candidate_id: int(np.sum(flags[candidate_id]))
                            for candidate_id in candidate_ids
                        },
                    },
                    ensure_ascii=False,
                )
                + "\n"
            )


def _run_task_groups(
    config: dict,
    scene: dict,
    manifest: dict,
    candidates: Sequence[dict],
    receiver: str,
    tasks: Sequence[dict],
) -> None:
    if not tasks:
        return
    scenario_id = str(scene["scenario_id"])
    n_rx = _scene_n_rx(scene)
    by_id = {str(row["candidate_id"]): row for row in candidates}
    channel, grid, precoders, covariances = bler027.build_scene(
        manifest, candidates, scenario_id
    )
    mcs = get_mcs("nr_256qam", 8, None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(tb.cb_k_values, tb.cb_e_values, num_iter=8, llr_clip=50.0)
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    stored = _read_supplements(config, scene, receiver)
    groups: dict[float, list[dict]] = defaultdict(list)
    for task in tasks:
        groups[float(task["snr_db"])].append(dict(task))
    for snr_db, snr_tasks in sorted(groups.items()):
        task_by_id = {str(task["candidate_id"]): task for task in snr_tasks}
        candidate_ids = sorted(task_by_id)
        flags = {
            candidate_id: np.zeros(int(task_by_id[candidate_id]["trials"]), dtype=bool)
            for candidate_id in candidate_ids
        }
        nmse_sum = {candidate_id: 0.0 for candidate_id in candidate_ids}
        nmse_sumsq = {candidate_id: 0.0 for candidate_id in candidate_ids}
        first_trial = min(int(task["trial_start"]) for task in snr_tasks)
        last_trial = max(int(task["trial_end"]) for task in snr_tasks)
        snr_linear = 10.0 ** (snr_db / 10.0)
        noise_variance = 8.0 / snr_linear
        ls_noise_variance = noise_variance / 2.0
        filters = (
            {
                candidate_id: build_frequency_rmmse_filter(
                    covariances[candidate_id], pilot_local, ls_noise_variance, diagonal_loading=1e-10
                )
                for candidate_id in candidate_ids
            }
            if receiver == "estimated"
            else {}
        )
        started = time.time()
        for absolute_start in range(first_trial, last_trial + 1, int(config["batch_size"])):
            active = [
                candidate_id
                for candidate_id in candidate_ids
                if int(task_by_id[candidate_id]["trial_start"]) <= absolute_start
                <= int(task_by_id[candidate_id]["trial_end"])
            ]
            if not active:
                continue
            current = min(int(config["batch_size"]), last_trial - absolute_start + 1)
            payload_rng = np.random.default_rng(
                stable_seed(config["seed"], scenario_id, snr_db, absolute_start, "payload")
            )
            payload = [payload_rng.integers(0, 2, size=int(k), dtype=np.int8) for k in tb.cb_k_values]
            symbols = qam_modulate(np.concatenate(adapter.encode(payload)), int(mcs.qm))
            realization = generate_sionna_tdl_channel(
                grid,
                channel,
                n_tx=8,
                n_rx=n_rx,
                batch_size=current,
                seed=stable_seed(config["seed"], scenario_id, snr_db, absolute_start, "channel"),
            )
            noise_rng = np.random.default_rng(
                stable_seed(config["seed"], scenario_id, snr_db, absolute_start, "noise")
            )
            averaged_ls_noise = math.sqrt(ls_noise_variance / 2.0) * (
                noise_rng.normal(size=(current, n_rx, grid.pilot_count))
                + 1j * noise_rng.normal(size=(current, n_rx, grid.pilot_count))
            )
            data_noise = math.sqrt(noise_variance / 2.0) * (
                noise_rng.normal(size=(current, n_rx, grid.n_data_re))
                + 1j * noise_rng.normal(size=(current, n_rx, grid.n_data_re))
            )
            decode_batches = []
            for candidate_id in active:
                task = task_by_id[candidate_id]
                used = min(current, int(task["trial_end"]) - absolute_start + 1)
                offset = absolute_start - int(task["trial_start"])
                effective = equivalent_channel(realization.H, precoders[candidate_id].C)[:used]
                true_data = effective[:, :, grid.data_symbol_indices, data_local]
                received = true_data * symbols[None, None, :] + data_noise[:used]
                if receiver == "estimated":
                    pilot0 = effective[:, :, int(grid.pilot_symbol_indices[0]), pilot_local]
                    pilot1 = effective[:, :, int(grid.pilot_symbol_indices[-1]), pilot_local]
                    estimate_full = filters[candidate_id].estimate_full_band(
                        0.5 * (pilot0 + pilot1) + averaged_ls_noise[:used]
                    )
                    estimate_data = estimate_full[:, :, data_local]
                    _, _, ratios = _ce_trial_statistics(estimate_data, true_data)
                    nmse_sum[candidate_id] += float(np.sum(ratios))
                    nmse_sumsq[candidate_id] += float(np.sum(ratios**2))
                    equalized, effective_noise = _estimated_csi_mrc_equalize(
                        estimate_data, received, noise_variance
                    )
                else:
                    equalized, effective_noise = curves028.ideal_csi_equalize(
                        true_data, received, noise_variance
                    )
                decode_batches.append(
                    (candidate_id, offset, used, equalized, effective_noise)
                )
            combined_equalized = np.concatenate(
                [item[3] for item in decode_batches], axis=0
            )
            combined_noise = np.concatenate(
                [item[4] for item in decode_batches], axis=0
            )
            combined_flags = _decode_flags(
                adapter,
                combined_equalized,
                combined_noise,
                int(mcs.qm),
                payload,
            )
            cursor = 0
            for candidate_id, offset, used, _, _ in decode_batches:
                flags[candidate_id][offset : offset + used] = combined_flags[
                    cursor : cursor + used
                ]
                cursor += used
            completed = absolute_start + current - 1
            if absolute_start == first_trial or completed % 500 < current or completed == last_trial:
                print(
                    f"[BLER {scenario_id} {receiver}] snr={snr_db:g} "
                    f"absolute_trials={completed}/{last_trial} active={len(active)} "
                    f"elapsed={time.time()-started:.1f}s",
                    flush=True,
                )
        for candidate_id in candidate_ids:
            task = task_by_id[candidate_id]
            trial_start = int(task["trial_start"])
            trial_end = int(task["trial_end"])
            trial_count = int(task["trials"])
            flag_path = _flag_path(
                config, scene, receiver, snr_db, candidate_id, trial_start, trial_end
            )
            flag_path.parent.mkdir(parents=True, exist_ok=True)
            np.save(flag_path, flags[candidate_id])
            errors = int(np.sum(flags[candidate_id]))
            low, high = wilson(errors, trial_count)
            stored.append(
                {
                    "scenario_id": scenario_id,
                    "candidate_id": candidate_id,
                    "family": by_id[candidate_id]["family"],
                    "receiver": receiver,
                    "n_rx": n_rx,
                    "snr_db": snr_db,
                    "trial_start": trial_start,
                    "trial_end": trial_end,
                    "trials": trial_count,
                    "tb_errors": errors,
                    "bler": errors / trial_count,
                    "bler_wilson95_lo": low,
                    "bler_wilson95_hi": high,
                    "ce_nmse_sum": nmse_sum[candidate_id] if receiver == "estimated" else "",
                    "ce_nmse_sumsq": nmse_sumsq[candidate_id] if receiver == "estimated" else "",
                    "analytic_trace_nmse": (
                        mismatched_frequency_lmmse_trace_nmse(
                            covariances[candidate_id],
                            filters[candidate_id].weights,
                            pilot_local,
                            data_local,
                            ls_noise_variance,
                        )
                        if receiver == "estimated"
                        else ""
                    ),
                    "seed": int(config["seed"]),
                    "error_flags": _repo_relative(flag_path),
                }
            )
        stored.sort(
            key=lambda row: (
                str(row["candidate_id"]),
                float(row["snr_db"]),
                int(row["trial_start"]),
            )
        )
        _write_rows_atomic(stored, _supplement_path(config, scene, receiver))


def _write_styles_and_mapping(output: Path, candidates: Sequence[dict]) -> tuple[dict, dict, dict]:
    ordered = sorted(candidates, key=curves028._candidate_sort_key)
    styles = {}
    candidate_labels = {}
    abbreviations = {}
    for index, candidate in enumerate(ordered, start=1):
        candidate_id = str(candidate["candidate_id"])
        style = dict(candidate.get("style") or bler027.curve_style(candidate_id, str(candidate["family"])))
        styles[candidate_id] = {
            "color": style["color"],
            "linestyle": style["linestyle"],
            "marker": style["marker"],
        }
        candidate_labels[candidate_id] = f"candidate {index:02d}"
        abbreviations[candidate_id] = str(
            candidate.get("label") or bler027.curve_label(candidate_id)
        )
    output.mkdir(parents=True, exist_ok=True)
    (output / "curve_styles.json").write_text(
        json.dumps(styles, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return styles, candidate_labels, abbreviations


def merge_and_plot_scene(config: dict, scene: dict, manifest: dict, candidates: Sequence[dict]) -> None:
    output = _scene_output(config, scene) / "final"
    styles, candidate_labels, abbreviations = _write_styles_and_mapping(output, candidates)
    scenario_id = str(scene["scenario_id"])
    n_rx = _scene_n_rx(scene)
    merged_by_receiver = {}
    for receiver in RECEIVERS:
        if receiver not in scene.get("receivers", {}):
            continue
        base = _base_rows(scene, receiver, candidates)
        merged = merge_rows(base, _read_supplements(config, scene, receiver))
        policy = scene["receivers"][receiver]
        included = _policy_plot_keys(merged, policy)
        for row in merged:
            row["plot_included"] = (
                str(row["candidate_id"]), float(row["snr_db"])
            ) in included
        merged_by_receiver[receiver] = merged
        _write_rows_atomic(merged, output / f"{receiver}_csi_bler_points.csv")
        plotted = [
            row for row in merged if bool(row["plot_included"]) and int(row["trials"]) > 0
        ]
        for mode, mapping in (("candidate", candidate_labels), ("abbrev", abbreviations)):
            curves028._plot_metric(
                plotted,
                mapping,
                styles,
                "bler",
                f"{receiver.capitalize()}-CSI BLER",
                f"{scenario_id} {n_rx}Rx comb-{manifest['dmrs_spacing_subcarriers']} {receiver}-CSI BLER",
                output / f"{scenario_id.lower()}_comb{manifest['dmrs_spacing_subcarriers']}_{mode}_{receiver}_csi_bler.png",
                True,
            )
    if "estimated" in merged_by_receiver:
        plotted = [
            row
            for row in merged_by_receiver["estimated"]
            if bool(row["plot_included"]) and int(row["trials"]) > 0
        ]
        for mode, mapping in (("candidate", candidate_labels), ("abbrev", abbreviations)):
            curves028._plot_metric(
                plotted,
                mapping,
                styles,
                "ce_nmse_mean_db",
                "CE NMSE (dB)",
                f"{scenario_id} {n_rx}Rx comb-{manifest['dmrs_spacing_subcarriers']} CE NMSE",
                output / f"{scenario_id.lower()}_comb{manifest['dmrs_spacing_subcarriers']}_{mode}_ce_nmse.png",
                False,
            )


def validate_config(config: dict) -> None:
    seen = set()
    for scene in config["scenes"]:
        scenario_id = str(scene["scenario_id"])
        n_rx = _scene_n_rx(scene)
        if scenario_id in seen:
            raise ValueError(f"Duplicate scenario_id: {scenario_id}")
        seen.add(scenario_id)
        manifest, _ = _load_manifest(scene)
        candidates = _selected_candidates(scene, manifest)
        if _is_transparent_scene(scene):
            expected = {
                "A100_PRG_DFT8_4RB": (4, "cycle8_random_tail4"),
                "A100_PRG_DFT8_6RB": (6, "cycle_all"),
            }
            actual = {
                str(row["candidate_id"]): (int(row["prg_size_rb"]), str(row["mapping"]))
                for row in candidates
            }
            plan_section = int(scene.get("plan_section", 10))
            if plan_section in (14, 15):
                if not actual or any(
                    candidate_id not in expected or definition != expected[candidate_id]
                    for candidate_id, definition in actual.items()
                ):
                    raise ValueError(f"Transparent PRG candidate definitions changed: {actual}")
            elif actual != expected:
                raise ValueError(f"Transparent PRG candidate definitions changed: {actual}")
            receiver_set = set(scene.get("receivers", {}))
            if plan_section == 15:
                if actual != {"A100_PRG_DFT8_6RB": expected["A100_PRG_DFT8_6RB"]}:
                    raise ValueError("Plan-028 section 15 requires only transparent PRG 6 RB.")
                if receiver_set != {"estimated"}:
                    raise ValueError("Plan-028 section 15 requires only the estimated receiver.")
                if n_rx != 1:
                    raise ValueError("Plan-028 section 15 requires n_rx=1.")
            elif receiver_set != set(RECEIVERS):
                raise ValueError("Transparent PRG mode requires estimated and ideal receivers.")
            paired = scene.get("paired_policy") or {}
            for field in ("target_errors", "min_total_trials", "max_total_trials"):
                if int(paired.get(field, 0)) <= 0:
                    raise ValueError(f"paired_policy.{field} must be positive")
            if int(paired["min_total_trials"]) > int(paired["max_total_trials"]):
                raise ValueError("paired_policy min_total_trials exceeds max_total_trials")
            grids = set()
            for receiver in receiver_set:
                receiver_policy = scene["receivers"][receiver]
                if receiver_policy.get("snr_db"):
                    grid_values = tuple(
                        float(value) for value in receiver_policy.get("snr_db", [])
                    )
                elif receiver_policy.get("base_csv"):
                    grid_values = tuple(
                        sorted(
                            {
                                float(row["snr_db"])
                                for row in _base_rows(scene, receiver, candidates)
                            }
                        )
                    )
                else:
                    grid_values = ()
                grids.add(grid_values)
            if plan_section == 14:
                if len(grids) != 1:
                    raise ValueError("Plan-028 section 14 PRG receivers require one common SNR grid.")
                grid_values = next(iter(grids))
                if not grid_values or any(
                    right <= left for left, right in zip(grid_values, grid_values[1:])
                ):
                    raise ValueError("Plan-028 section 14 PRG SNR grid must increase.")
            elif plan_section == 15:
                expected_grid = (
                    (16.5, 17.0, 17.5, 18.0)
                    if str(scene.get("run_kind", "formal")) == "prescan"
                    else (16.25, 16.5, 16.75, 17.0, 17.25, 17.5, 17.75)
                )
                if grids != {expected_grid}:
                    raise ValueError(
                        "Plan-028 section 15 PRG SNR grid differs from the frozen grid."
                    )
                run_grid = tuple(float(value) for value in scene.get("run_snr_db", expected_grid))
                if not run_grid or any(value not in expected_grid for value in run_grid):
                    raise ValueError("Plan-028 section 15 run_snr_db must be a non-empty grid subset.")
                if len(set(run_grid)) != len(run_grid):
                    raise ValueError("Plan-028 section 15 run_snr_db contains duplicates.")
            elif grids != {(14.0, 14.25, 14.5, 14.75, 15.0, 15.5, 15.75, 16.0)}:
                raise ValueError("Transparent PRG SNR grid must match result-028 section 3.6.")
            _, grid, _ = _build_transparent_scene(manifest, scenario_id)
            pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
            for row in candidates:
                prg_sc = 12 * int(row["prg_size_rb"])
                counts = [
                    int(np.sum((pilot_local >= start) & (pilot_local < start + prg_sc)))
                    for start in range(0, grid.n_sc, prg_sc)
                ]
                expected_pilots = 8 if int(row["prg_size_rb"]) == 4 else 12
                if set(counts) != {expected_pilots}:
                    raise RuntimeError(f"Unexpected PRG pilot counts for {row['candidate_id']}: {counts}")
        if _is_transparent_cdd_scene(scene):
            expected = {
                "A100_B0_QC_TRANSPARENT_CDD": "A100_B0_QC",
                "A100_AP_RMS_T1_TRANSPARENT_CDD": "A100_AP_RMS_T1",
                "A100_AP_TU_NT_TRANSPARENT_CDD": "A100_AP_TU_NT",
                "A100_S0_SIDON_TRANSPARENT_CDD": "A100_S0_SIDON",
                "A100_MEFF_T2_06_TRANSPARENT_CDD": "A100_MEFF_T2_06",
            }
            actual_sources = {
                str(row["candidate_id"]): str(row["source_candidate_id"])
                for row in candidates
            }
            for candidate_id, source_id in actual_sources.items():
                if candidate_id in expected and source_id != expected[candidate_id]:
                    raise ValueError(
                        f"Transparent CDD candidate mapping changed for {candidate_id}: {source_id}"
                    )
            receiver_modes = sorted(str(value) for value in scene.get("receivers", {}))
            if len(receiver_modes) != 1 or receiver_modes[0] not in RECEIVERS:
                raise ValueError(
                    "Transparent CDD mode requires exactly one estimated or ideal receiver."
                )
            receiver = receiver_modes[0]
            covariance_mode = str(scene.get("covariance_mode", "physical")).strip().lower()
            if covariance_mode not in {"physical", "matched_effective"}:
                raise ValueError(
                    "Transparent CDD covariance_mode must be physical or matched_effective."
                )
            plan_section = int(scene.get("plan_section", 12))
            if receiver == "ideal" and covariance_mode != "physical":
                raise ValueError("Ideal CDD receiver does not accept covariance_mode overrides.")
            if plan_section == 16:
                expected_coordinates = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75]
                if n_rx != 1 or len(candidates) != 1:
                    raise ValueError("Plan-028 section 16 requires one 1Rx small-delay CDD candidate.")
                if [float(value) for value in candidates[0]["delay_grid_coordinates"]] != expected_coordinates:
                    raise ValueError("Plan-028 section 16 small-delay CDD coordinates changed.")
                if receiver == "estimated" and covariance_mode != "matched_effective":
                    raise ValueError(
                        "Plan-028 section 16 estimated receiver requires matched_effective covariance."
                    )
            paired = scene.get("paired_policy") or {}
            for field in ("target_errors", "min_total_trials", "max_total_trials"):
                if int(paired.get(field, 0)) <= 0:
                    raise ValueError(f"paired_policy.{field} must be positive")
            if int(paired["min_total_trials"]) > int(paired["max_total_trials"]):
                raise ValueError("paired_policy min_total_trials exceeds max_total_trials")
            receiver_policy = scene["receivers"][receiver]
            if receiver_policy.get("snr_db"):
                grid_values = tuple(
                    float(value) for value in receiver_policy.get("snr_db", [])
                )
            elif receiver_policy.get("base_csv"):
                grid_values = tuple(
                    sorted(
                        {
                            float(row["snr_db"])
                            for row in _base_rows(scene, receiver, candidates)
                        }
                    )
                )
            else:
                grid_values = ()
            if not grid_values or any(
                right <= left for left, right in zip(grid_values, grid_values[1:])
            ):
                raise ValueError("Transparent CDD SNR grid must be non-empty and increasing.")
            if plan_section == 16 and receiver == "ideal" and grid_values != (
                14.0,
                14.5,
                15.0,
                15.5,
                16.0,
                16.5,
                17.0,
                17.5,
                18.0,
                18.5,
                19.0,
            ):
                raise ValueError("Plan-028 section 16 ideal SNR grid differs from the frozen grid.")
            if plan_section == 16 and receiver == "estimated" and grid_values not in {
                (14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0),
                (
                    14.0,
                    14.5,
                    15.0,
                    15.5,
                    16.0,
                    16.5,
                    17.0,
                    17.5,
                    18.0,
                    18.5,
                    19.0,
                    19.5,
                    20.0,
                ),
            }:
                raise ValueError(
                    "Plan-028 section 16 matched SNR grid differs from prescan/formal grids."
                )
            if all(str(row.get("source_candidate_id", "")) for row in candidates):
                if receiver != "estimated":
                    raise ValueError(
                        "Legacy source transparent CDD scenes require estimated receiver."
                    )
                if grid_values != (
                    14.0,
                    14.25,
                    14.5,
                    14.75,
                    15.0,
                    15.5,
                    15.75,
                    16.0,
                ):
                    raise ValueError(
                        "Transparent CDD SNR grid must match result-028 section 3.6."
                    )
            channel, grid, precoders, true_covariances = bler027.build_scene(
                manifest, candidates, scenario_id
            )
            base_covariance = tdl_active_frequency_covariance(grid, channel)
            if not np.allclose(np.real(np.diag(base_covariance)), 1.0, atol=1e-12):
                raise RuntimeError("Physical covariance power normalization changed.")
            for row in candidates:
                candidate_id = str(row["candidate_id"])
                source_id = str(row.get("source_candidate_id", ""))
                if source_id:
                    source = next(
                        item
                        for item in bler027.scenario_candidates(manifest, scenario_id)
                        if str(item["candidate_id"]) == source_id
                    )
                    if [float(value) for value in row["delay_grid_coordinates"]] != [
                        float(value) for value in source["delay_grid_coordinates"]
                    ]:
                        raise RuntimeError(f"CDD transmitter changed for {candidate_id}.")
                elif str(row.get("transmitter_definition")) != "explicit_delay_grid_coordinates":
                    raise RuntimeError(
                        f"Explicit transparent CDD transmitter is not frozen for {candidate_id}."
                    )
                expected_covariance = base_covariance * (
                    precoders[candidate_id].C @ precoders[candidate_id].C.conj().T
                )
                if not np.allclose(true_covariances[candidate_id], expected_covariance):
                    raise RuntimeError(f"True CDD covariance changed for {candidate_id}.")
        for receiver in scene.get("receivers", {}):
            if receiver not in RECEIVERS:
                raise ValueError(f"Unsupported receiver: {receiver}")
            base = _base_rows(scene, receiver, candidates)
            _validate_intervals(base, _read_supplements(config, scene, receiver))
            policy = scene["receivers"][receiver]
            if int(policy["max_total_trials"]) <= 0:
                raise ValueError("max_total_trials must be positive")
            overrides = policy.get("candidate_overrides", {})
            if not isinstance(overrides, dict):
                raise ValueError("candidate_overrides must be a mapping")
            selected_ids = {str(row["candidate_id"]) for row in candidates}
            unknown = set(overrides) - selected_ids
            if unknown:
                raise ValueError(f"Unknown candidate_overrides: {sorted(unknown)}")
            for candidate_id in overrides:
                effective = _candidate_policy(policy, candidate_id)
                if int(effective["max_total_trials"]) <= 0:
                    raise ValueError("max_total_trials must be positive")
        if _is_transparent_scene(scene):
            receiver_ids = tuple(
                receiver for receiver in RECEIVERS if receiver in scene.get("receivers", {})
            )
            rows_by_receiver = {
                receiver: merge_rows(
                    _base_rows(scene, receiver, candidates),
                    _read_supplements(config, scene, receiver),
                )
                for receiver in receiver_ids
            }
            build_paired_tasks(rows_by_receiver, scene["paired_policy"], int(config["batch_size"]))
        if _is_transparent_cdd_scene(scene):
            receiver = next(iter(scene["receivers"]))
            rows = merge_rows(
                _base_rows(scene, receiver, candidates),
                _read_supplements(config, scene, receiver),
            )
            build_paired_tasks(
                {receiver: rows}, scene["paired_policy"], int(config["batch_size"])
            )
        print(
            f"[validate] {scenario_id}: candidates={len(candidates)} "
            f"n_rx={n_rx} receivers={sorted(scene.get('receivers', {}))}",
            flush=True,
        )


def run_config(config: dict) -> None:
    for scene in config["scenes"]:
        manifest, digest = _load_manifest(scene)
        candidates = _selected_candidates(scene, manifest)
        scene_output = _scene_output(config, scene)
        scene_output.mkdir(parents=True, exist_ok=True)
        transparent_manifest_path = None
        transparent_manifest_sha256 = None
        if _is_transparent_scene(scene):
            transparent_manifest_path, transparent_manifest_sha256 = _write_transparent_manifest(
                config, scene, digest, candidates
            )
        if _is_transparent_cdd_scene(scene):
            transparent_manifest_path, transparent_manifest_sha256 = (
                _write_transparent_cdd_manifest(config, scene, digest, candidates)
            )
        resolved = {
            "schema": SCHEMA,
            "scenario_id": scene["scenario_id"],
            "n_rx": _scene_n_rx(scene),
            "manifest_sha256": digest,
            "seed": config["seed"],
            "batch_size": config["batch_size"],
            "config_path": config["config_path"],
            "config_sha256": _sha256(Path(config["config_path"])),
            "receivers": scene.get("receivers", {}),
            "absolute_trial_seeded": True,
            "mode": scene.get("mode", "cdd_manifest"),
            "transparent_manifest_path": _repo_relative(transparent_manifest_path)
            if transparent_manifest_path
            else None,
            "transparent_manifest_sha256": transparent_manifest_sha256,
        }
        (scene_output / "resolved_run.json").write_text(
            json.dumps(resolved, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        if _is_transparent_scene(scene):
            receiver_ids = tuple(
                receiver for receiver in RECEIVERS if receiver in scene.get("receivers", {})
            )
            for _ in range(int(config["max_passes"])):
                rows_by_receiver = {
                    receiver: merge_rows(
                        _base_rows(scene, receiver, candidates),
                        _read_supplements(config, scene, receiver),
                    )
                    for receiver in receiver_ids
                }
                tasks = build_paired_tasks(
                    rows_by_receiver, scene["paired_policy"], int(config["batch_size"])
                )
                tasks = _filter_run_snr_tasks(scene, tasks)
                if not tasks:
                    break
                task_limit = int(scene.get("max_snr_tasks_per_run", 0))
                if task_limit > 0:
                    tasks = tasks[:task_limit]
                _run_transparent_paired_tasks(config, scene, manifest, candidates, tasks)
            merge_and_plot_scene(config, scene, manifest, candidates)
            continue
        if _is_transparent_cdd_scene(scene):
            receiver = next(iter(scene["receivers"]))
            for _ in range(int(config["max_passes"])):
                rows = merge_rows(
                    _base_rows(scene, receiver, candidates),
                    _read_supplements(config, scene, receiver),
                )
                tasks = build_paired_tasks(
                    {receiver: rows}, scene["paired_policy"], int(config["batch_size"])
                )
                tasks = _filter_run_snr_tasks(scene, tasks)
                if not tasks:
                    break
                task_limit = int(scene.get("max_snr_tasks_per_run", 0))
                if task_limit > 0:
                    tasks = tasks[:task_limit]
                _run_transparent_cdd_paired_tasks(
                    config, scene, manifest, candidates, tasks
                )
            merge_and_plot_scene(config, scene, manifest, candidates)
            continue
        for receiver in RECEIVERS:
            if receiver not in scene.get("receivers", {}):
                continue
            for _ in range(int(config["max_passes"])):
                base = _base_rows(scene, receiver, candidates)
                supplements = _read_supplements(config, scene, receiver)
                merged = merge_rows(base, supplements)
                tasks = build_tasks(merged, scene["receivers"][receiver], int(config["batch_size"]))
                if not tasks:
                    break
                _run_task_groups(config, scene, manifest, candidates, receiver, tasks)
        merge_and_plot_scene(config, scene, manifest, candidates)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--stage", choices=("validate", "run", "merge", "all"), default="all")
    args = parser.parse_args()
    config = load_runner_config(args.config.resolve())
    debug_seconds = int(config.get("debug_traceback_seconds", 0))
    debug_handle = None
    if debug_seconds > 0:
        debug_path = Path(config["output_dir"]) / "debug_stack.log"
        debug_path.parent.mkdir(parents=True, exist_ok=True)
        debug_handle = debug_path.open("a", encoding="utf-8")
        faulthandler.dump_traceback_later(debug_seconds, repeat=True, file=debug_handle)
    validate_config(config)
    if args.stage in ("run", "all"):
        run_config(config)
    elif args.stage == "merge":
        for scene in config["scenes"]:
            manifest, _ = _load_manifest(scene)
            merge_and_plot_scene(config, scene, manifest, _selected_candidates(scene, manifest))
    if debug_seconds > 0:
        faulthandler.cancel_dump_traceback_later()
        if debug_handle is not None:
            debug_handle.close()


if __name__ == "__main__":
    main()
