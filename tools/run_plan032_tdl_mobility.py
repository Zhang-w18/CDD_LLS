"""Run plan-032 A100 TDL-A mobility curves with two-DMRS 2D RMMSE."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import yaml

from cdd_lls.core.config import ChannelConfig
from cdd_lls.core.mcs import build_tb_layout, get_mcs
from cdd_lls.phy.channel_tdl import generate_sionna_tdl_channel_active
from cdd_lls.phy.estimators import (
    TDLTimeFrequencyCovariance,
    build_prg_time_frequency_rmmse_filter,
    build_time_frequency_rmmse_filter,
    tdl_active_time_frequency_covariance,
)
from cdd_lls.phy.ldpc import SionnaLDPCAdapter
from cdd_lls.phy.precoding import (
    build_active_dft_grid_precoder,
    build_prg_dft_precoder_batch,
    equivalent_channel,
)
from cdd_lls.phy.qam import qam_modulate
from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
from tools import run_plan027_bler as bler027
from tools.run_bler_curves import (
    _ce_trial_statistics,
    _decode_flags,
    _estimated_csi_mrc_equalize,
)
from tools.run_plan025_delay_matched_tdl import (
    read_csv_rows,
    stable_seed,
    wilson,
    write_csv_rows,
)


SCHEMA = "plan032-tdl-mobility-v1"
SCENARIO = "A100_V60"
PRG_CANDIDATE = "A100_PRG_DFT8_6RB"
MRT_CANDIDATE = "A100_AGED_CSI_MRT_PRG6_SLOTS10"
ORIGINAL_CURVE_SET = "original"
TRANSPARENT_CDD_SUPPLEMENT = "transparent_cdd_supplement"
TRANSPARENT_CDD_SOURCES = ("A100_S0_SIDON", "A100_B0_QC")
TRANSPARENT_CDD_IDS = {
    source_id: f"{source_id}_TRANSPARENT_CDD" for source_id in TRANSPARENT_CDD_SOURCES
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _resolve(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def _repo_relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return str(path.resolve())


def _write_rows_atomic(rows: Sequence[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    write_csv_rows(rows, temporary)
    temporary.replace(path)


def _snr_key(snr_db: float) -> str:
    return f"{float(snr_db):+.2f}".replace("+", "p").replace("-", "m").replace(".", "p")


def _safe(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def _array_path(output: Path, category: str, snr_db: float, candidate_id: str,
                trial_start: int, trial_end: int) -> Path:
    return (
        output
        / category
        / _snr_key(snr_db)
        / f"{_safe(candidate_id)}_t{trial_start:06d}_{trial_end:06d}.npy"
    )


def load_config(path: Path) -> dict:
    config = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if config.get("schema") != SCHEMA:
        raise ValueError(f"config.schema must be {SCHEMA!r}")
    config["config_path"] = str(path.resolve())
    config["output_dir"] = str(_resolve(config["output_dir"]))
    config["seed"] = int(config.get("seed", 20260727))
    config["batch_size"] = int(config.get("batch_size", 10))
    config["snr_db"] = [float(value) for value in config["snr_db"]]
    scalar_target = config.get("target_total_trials")
    per_snr_targets = config.get("target_total_trials_by_snr")
    if (scalar_target is None) == (per_snr_targets is None):
        raise ValueError(
            "Exactly one of target_total_trials and target_total_trials_by_snr is required"
        )
    if scalar_target is not None:
        config["target_total_trials"] = int(scalar_target)
    else:
        if not isinstance(per_snr_targets, dict):
            raise ValueError("target_total_trials_by_snr must be a mapping")
        config["target_total_trials_by_snr"] = {
            float(snr_db): int(total) for snr_db, total in per_snr_targets.items()
        }
    config["speed_kmh"] = float(config.get("speed_kmh", 60.0))
    config["feedback_period_slots"] = int(config.get("feedback_period_slots", 10))
    config["history_symbol_index"] = int(config.get("history_symbol_index", 140))
    config["run_kind"] = str(config["run_kind"])
    config["curve_set"] = str(config.get("curve_set", ORIGINAL_CURVE_SET))
    return config


def _trial_target(config: dict, snr_db: float) -> int:
    if "target_total_trials" in config:
        return int(config["target_total_trials"])
    return int(config["target_total_trials_by_snr"][float(snr_db)])


def _load_source_manifest(config: dict) -> tuple[dict, str]:
    bundle = config["manifest"]
    path = _resolve(bundle["path"])
    digest = _sha256(path)
    sha_path = _resolve(bundle["sha256_path"])
    recorded = sha_path.read_text(encoding="utf-8").split()[0]
    if recorded.lower() != digest.lower():
        raise RuntimeError("Source manifest hash mismatch.")
    approval = json.loads(_resolve(bundle["approval_path"]).read_text(encoding="utf-8"))
    if approval.get("approved") is not True:
        raise RuntimeError("Source manifest approval is not affirmative.")
    if str(approval.get("manifest_sha256", "")).lower() != digest.lower():
        raise RuntimeError("Source approval hash mismatch.")
    return json.loads(path.read_text(encoding="utf-8")), digest


def _candidate_definitions(
    manifest: dict, curve_set: str = ORIGINAL_CURVE_SET
) -> list[dict]:
    rows = bler027.scenario_candidates(manifest, "A100")
    if len(rows) != 10:
        raise RuntimeError("Plan-032 requires the ten A100 source CDD candidates.")
    rows = [dict(row) for row in rows]
    if curve_set == TRANSPARENT_CDD_SUPPLEMENT:
        by_id = {str(row["candidate_id"]): row for row in rows}
        missing = [
            source_id for source_id in TRANSPARENT_CDD_SOURCES if source_id not in by_id
        ]
        if missing:
            raise RuntimeError(f"Missing transparent CDD source candidates: {missing}")
        return [
            {
                "candidate_id": TRANSPARENT_CDD_IDS[source_id],
                "source_candidate_id": source_id,
                "family": "TRANSPARENT_CDD_PHYSCOV",
                "label": f"{source_id.removeprefix('A100_')}, transparent receiver",
                "receiver_knowledge": "speed and physical TDL covariance; no CDD delays",
            }
            for source_id in TRANSPARENT_CDD_SOURCES
        ]
    if curve_set != ORIGINAL_CURVE_SET:
        raise ValueError(
            f"curve_set must be {ORIGINAL_CURVE_SET!r} or "
            f"{TRANSPARENT_CDD_SUPPLEMENT!r}"
        )
    rows.append(
        {
            "candidate_id": PRG_CANDIDATE,
            "family": "TRANSPARENT_PRG_DFT",
            "label": "transparent PRG 6 RB",
            "receiver_knowledge": "speed and physical TDL covariance; no precoder",
        }
    )
    rows.append(
        {
            "candidate_id": MRT_CANDIDATE,
            "family": "AGED_CSI_MRT",
            "label": "aged-CSI MRT, 6 RB, slots10",
            "receiver_knowledge": "speed and physical TDL covariance; no precoder",
        }
    )
    return rows


def _build_scene(manifest: dict, speed_kmh: float):
    resource = bler027.resource_config(int(manifest["dmrs_spacing_subcarriers"]))
    grid = build_resource_grid(resource)
    channel = ChannelConfig(
        backend="sionna_tdl",
        model="3gpp_tr38901_tdl",
        tdl_profile="A",
        delay_spread_ns=100.0,
        carrier_frequency_hz=3.5e9,
        ue_speed_kmh=float(speed_kmh),
        num_sinusoids=20,
    )
    source = bler027.scenario_candidates(manifest, "A100")
    precoders = {
        str(row["candidate_id"]): build_active_dft_grid_precoder(
            grid,
            [float(value) for value in row["delay_grid_coordinates"]],
            n_tx=8,
            normalize=False,
        )
        for row in source
    }
    base = tdl_active_time_frequency_covariance(grid, channel)
    return channel, grid, precoders, base


def build_aged_mrt_prg_precoder(old_channel: np.ndarray, prg_size_rb: int = 6) -> np.ndarray:
    """Return exact unquantized PRG MRT weights from stale [batch,rx,tx,sc] CSI."""
    old = np.asarray(old_channel, dtype=np.complex128)
    if old.ndim != 4 or old.shape[1] != 1:
        raise ValueError("Aged MRT requires old_channel shape [batch,1,n_tx,n_sc].")
    batch, _, n_tx, n_sc = old.shape
    prg_size = 12 * int(prg_size_rb)
    if n_sc % prg_size:
        raise ValueError("PRG size must divide the active bandwidth.")
    weights = np.empty((batch, n_sc, n_tx), dtype=np.complex128)
    for start in range(0, n_sc, prg_size):
        stop = start + prg_size
        h = old[:, 0, :, start:stop]
        gram = np.einsum("bik,bjk->bij", h.conj(), h, optimize=True)
        _, eigenvectors = np.linalg.eigh(gram)
        vector = math.sqrt(float(n_tx)) * eigenvectors[:, :, -1]
        weights[:, start:stop, :] = vector[:, None, :]
    return weights


def _effective_channels(
    realization_h: np.ndarray,
    grid,
    precoders: dict[str, object],
    candidates: Sequence[dict] | None = None,
):
    old = realization_h[:, :, :, 0, :]
    current = realization_h[:, :, :, 1:, :]
    output = {
        candidate_id: equivalent_channel(current, precoder.C)
        for candidate_id, precoder in precoders.items()
    }
    orders = np.tile(np.arange(8, dtype=np.int64), (len(current), 1))
    cycling = build_prg_dft_precoder_batch(
        grid, n_tx=8, prg_size_rb=6, prg_vector_indices=orders, normalize=False
    )
    output[PRG_CANDIDATE] = equivalent_channel(current, cycling)
    aged_mrt = build_aged_mrt_prg_precoder(old, prg_size_rb=6)
    output[MRT_CANDIDATE] = equivalent_channel(current, aged_mrt)
    if candidates is not None:
        for row in candidates:
            source_id = row.get("source_candidate_id")
            if source_id is not None:
                output[str(row["candidate_id"])] = output[str(source_id)]
    return old, current, output, aged_mrt


def _build_filters(
    grid,
    base,
    precoders: dict[str, object],
    noise_variance: float,
    candidates: Sequence[dict] | None = None,
):
    filters = {}
    requested = (
        [str(row["candidate_id"]) for row in candidates]
        if candidates is not None
        else [*precoders, PRG_CANDIDATE, MRT_CANDIDATE]
    )
    for candidate_id in requested:
        if candidate_id not in precoders:
            continue
        precoder = precoders[candidate_id]
        covariance = TDLTimeFrequencyCovariance(
            time=base.time,
            frequency=base.frequency * (precoder.C @ precoder.C.conj().T),
            covariance_type=f"matched_cdd:{candidate_id}",
        )
        filters[candidate_id] = build_time_frequency_rmmse_filter(
            grid, covariance, noise_variance, diagonal_loading=1e-10
        )
    transparent_covariance = TDLTimeFrequencyCovariance(
        time=base.time,
        frequency=8.0 * base.frequency,
        covariance_type="transparent_8x_physical",
    )
    transparent_cdd_ids = [
        str(row["candidate_id"])
        for row in (candidates or [])
        if str(row.get("family")) == "TRANSPARENT_CDD_PHYSCOV"
    ]
    if transparent_cdd_ids:
        shared_cdd = build_time_frequency_rmmse_filter(
            grid,
            transparent_covariance,
            noise_variance,
            diagonal_loading=1e-10,
        )
        filters.update(
            {candidate_id: shared_cdd for candidate_id in transparent_cdd_ids}
        )
    prg_ids = [
        candidate_id
        for candidate_id in (PRG_CANDIDATE, MRT_CANDIDATE)
        if candidate_id in requested
    ]
    if prg_ids:
        shared_prg = build_prg_time_frequency_rmmse_filter(
            grid,
            transparent_covariance,
            prg_size_subcarriers=72,
            noise_variance=noise_variance,
            diagonal_loading=1e-10,
        )
        filters.update({candidate_id: shared_prg for candidate_id in prg_ids})
    if set(filters) != set(requested):
        missing = sorted(set(requested) - set(filters))
        raise RuntimeError(f"Missing receiver filters: {missing}")
    return filters


def validate_config(config: dict) -> tuple[dict, str, list[dict]]:
    if config["run_kind"] not in {"smoke", "prescan", "formal"}:
        raise ValueError("run_kind must be smoke, prescan, or formal")
    if "target_total_trials_by_snr" in config and set(
        config["target_total_trials_by_snr"]
    ) != set(config["snr_db"]):
        raise ValueError("target_total_trials_by_snr keys must exactly match snr_db")
    targets = [_trial_target(config, snr_db) for snr_db in config["snr_db"]]
    if config["batch_size"] <= 0 or any(target <= 0 for target in targets):
        raise ValueError("batch_size and trial targets must be positive")
    if any(target % config["batch_size"] for target in targets):
        raise ValueError("trial targets must align to batch_size")
    if config["speed_kmh"] != 60.0:
        raise ValueError("Plan-032 is frozen to 60 km/h")
    if config["feedback_period_slots"] != 10 or config["history_symbol_index"] != 140:
        raise ValueError("Plan-032 aged-CSI baseline is frozen to slots10/index140")
    if not config["snr_db"] or any(
        right <= left for left, right in zip(config["snr_db"], config["snr_db"][1:])
    ):
        raise ValueError("snr_db must be non-empty and strictly increasing")
    manifest, digest = _load_source_manifest(config)
    if int(manifest["dmrs_spacing_subcarriers"]) != 6:
        raise RuntimeError("Plan-032 requires comb-6 source candidates")
    candidates = _candidate_definitions(manifest, config["curve_set"])
    channel, grid, precoders, base = _build_scene(manifest, config["speed_kmh"])
    if grid.n_dmrs_re != 192 or list(np.unique(grid.pilot_symbol_indices)) != [2, 7]:
        raise RuntimeError("Plan-032 requires two DMRS symbols [2,7] and 192 pilot REs")
    if base.time.shape != (10, 10) or base.frequency.shape != (576, 576):
        raise RuntimeError("Unexpected TDL covariance dimensions")
    for candidate_id, precoder in precoders.items():
        power = np.sum(np.abs(precoder.C) ** 2, axis=1)
        if not np.allclose(power, 8.0, atol=1e-12):
            raise RuntimeError(f"Precoder power changed for {candidate_id}")
    print(
        f"[validate] {SCENARIO}: curve_set={config['curve_set']} "
        f"curves={len(candidates)} speed=60 km/h "
        f"dmrs_symbols=[2,7] pilot_re={grid.n_dmrs_re}",
        flush=True,
    )
    return manifest, digest, candidates


def _write_resolved(config: dict, manifest_sha256: str, candidates: Sequence[dict], grid) -> None:
    output = Path(config["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    actual_age_s = config["history_symbol_index"] * float(grid.ofdm_symbol_duration_s)
    payload = {
        "schema": SCHEMA,
        "plan": "research/plan-032-PDSCH-A100-60kmh.md",
        "run_kind": config["run_kind"],
        "curve_set": config["curve_set"],
        "config_path": _repo_relative(Path(config["config_path"])),
        "config_sha256": _sha256(Path(config["config_path"])),
        "source_manifest_sha256": manifest_sha256,
        "scenario_id": SCENARIO,
        "system": {
            "profile": "TDL-A",
            "delay_spread_ns": 100.0,
            "carrier_frequency_hz": 3.5e9,
            "ue_speed_kmh": config["speed_kmh"],
            "n_tx": 8,
            "n_rx": 1,
            "n_symbols": 10,
            "dmrs_symbol_indices": [2, 7],
            "pilot_re": grid.n_dmrs_re,
            "data_re": grid.n_data_re,
            "ici": False,
            "carrier_synchronized": True,
        },
        "receiver": {
            "type": "two-DMRS joint 2D time-frequency RMMSE",
            "knows_speed": True,
            "cdd_knows_precoder": config["curve_set"] == ORIGINAL_CURVE_SET,
            "transparent_cdd_knows_precoder": False,
            "transparent_prg_knows_precoder": False,
            "aged_mrt_knows_precoder": False,
            "per_ls_observation_noise_variance": "8/SNR_linear",
        },
        "aged_csi": {
            "nominal_feedback_period_slots": config["feedback_period_slots"],
            "nominal_feedback_period_ms": 5.0,
            "history_symbol_index": config["history_symbol_index"],
            "actual_age_ms": actual_age_s * 1e3,
            "precoder": "unquantized dominant eigenvector per 6-RB PRG",
            "codebook": None,
        },
        "snr_db": config["snr_db"],
        "target_total_trials": config.get("target_total_trials"),
        "target_total_trials_by_snr": {
            f"{snr_db:g}": _trial_target(config, snr_db) for snr_db in config["snr_db"]
        },
        "batch_size": config["batch_size"],
        "seed": config["seed"],
        "candidates": [dict(row) for row in candidates],
    }
    path = output / "resolved_run.json"
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    digest = _sha256(path)
    (output / "resolved_run.sha256").write_text(f"{digest}  {path.name}\n", encoding="utf-8")


def _merge_intervals(output: Path, candidates: Sequence[dict]) -> list[dict]:
    path = output / "intervals.csv"
    rows = [dict(row) for row in read_csv_rows(path)] if path.exists() else []
    grouped: dict[tuple[str, float], list[dict]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["candidate_id"]), float(row["snr_db"]))].append(row)
    by_id = {str(row["candidate_id"]): row for row in candidates}
    merged = []
    for (candidate_id, snr_db), values in sorted(grouped.items()):
        values.sort(key=lambda row: int(row["trial_start"]))
        expected = 1
        for row in values:
            if int(row["trial_start"]) != expected:
                raise RuntimeError(f"Trial gap/overlap for {candidate_id} at {snr_db:g} dB")
            expected = int(row["trial_end"]) + 1
        trials = sum(int(row["trials"]) for row in values)
        errors = sum(int(row["tb_errors"]) for row in values)
        nmse_sum = sum(float(row["ce_nmse_sum"]) for row in values)
        nmse_sumsq = sum(float(row["ce_nmse_sumsq"]) for row in values)
        low, high = wilson(errors, trials)
        mean = nmse_sum / trials
        variance = max((nmse_sumsq - trials * mean * mean) / max(trials - 1, 1), 0.0)
        se = math.sqrt(variance / trials)
        merged.append(
            {
                "scenario_id": SCENARIO,
                "candidate_id": candidate_id,
                "family": by_id[candidate_id]["family"],
                "receiver": "estimated",
                "snr_db": snr_db,
                "trials": trials,
                "tb_errors": errors,
                "bler": errors / trials,
                "bler_wilson95_lo": low,
                "bler_wilson95_hi": high,
                "ce_nmse_mean": mean,
                "ce_nmse_mean_db": 10.0 * math.log10(mean),
                "ce_nmse_standard_error": se,
                "ce_nmse_ci95_lo": max(mean - 1.96 * se, 0.0),
                "ce_nmse_ci95_hi": mean + 1.96 * se,
            }
        )
    final = output / "final"
    _write_rows_atomic(merged, final / "estimated_csi_bler_points.csv")
    return merged


def run(config: dict, manifest: dict, digest: str, candidates: Sequence[dict]) -> None:
    output = Path(config["output_dir"])
    channel, grid, precoders, base = _build_scene(manifest, config["speed_kmh"])
    _write_resolved(config, digest, candidates, grid)
    by_id = {str(row["candidate_id"]): row for row in candidates}
    candidate_ids = [str(row["candidate_id"]) for row in candidates]
    intervals_path = output / "intervals.csv"
    stored = [dict(row) for row in read_csv_rows(intervals_path)] if intervals_path.exists() else []
    mcs = get_mcs("nr_256qam", 8, None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(tb.cb_k_values, tb.cb_e_values, num_iter=8, llr_clip=50.0)
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarrier_indices)
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    requested_times = [0, *range(config["history_symbol_index"], config["history_symbol_index"] + 10)]
    for snr_db in config["snr_db"]:
        target_total_trials = _trial_target(config, snr_db)
        counts = {
            candidate_id: sum(
                int(row["trials"])
                for row in stored
                if str(row["candidate_id"]) == candidate_id
                and math.isclose(float(row["snr_db"]), snr_db, abs_tol=1e-12)
            )
            for candidate_id in candidate_ids
        }
        if len(set(counts.values())) != 1:
            raise RuntimeError(f"Paired trial counts differ at {snr_db:g} dB")
        current_total = next(iter(counts.values()))
        if current_total >= target_total_trials:
            print(f"[resume] snr={snr_db:g} trials={current_total}", flush=True)
            continue
        trial_start = current_total + 1
        trial_end = target_total_trials
        trial_count = trial_end - trial_start + 1
        if trial_count % config["batch_size"]:
            raise RuntimeError("Pending trial interval must align to batch_size")
        noise_variance = 8.0 / (10.0 ** (snr_db / 10.0))
        filters = _build_filters(grid, base, precoders, noise_variance, candidates)
        diagnostic = {
            candidate_id: {
                "condition_number": float(filters[candidate_id].condition_number),
                "numerical_jitter": float(filters[candidate_id].numerical_jitter),
                "covariance_type": filters[candidate_id].covariance_type,
            }
            for candidate_id in candidate_ids
        }
        diag_path = output / "filter_diagnostics" / f"{_snr_key(snr_db)}.json"
        diag_path.parent.mkdir(parents=True, exist_ok=True)
        diag_path.write_text(json.dumps(diagnostic, indent=2) + "\n", encoding="utf-8")
        flags = {candidate_id: np.zeros(trial_count, dtype=bool) for candidate_id in candidate_ids}
        ce_nmse = {
            candidate_id: np.zeros(trial_count, dtype=np.float64) for candidate_id in candidate_ids
        }
        old_current_cross = 0.0 + 0.0j
        old_power = 0.0
        current_power = 0.0
        current_edge_cross = 0.0 + 0.0j
        current_first_power = 0.0
        current_last_power = 0.0
        started = time.time()
        for absolute_start in range(trial_start, trial_end + 1, config["batch_size"]):
            batch = min(config["batch_size"], trial_end - absolute_start + 1)
            offset = absolute_start - trial_start
            payload_rng = np.random.default_rng(
                stable_seed(config["seed"], SCENARIO, snr_db, absolute_start, "payload")
            )
            payload = [payload_rng.integers(0, 2, size=int(k), dtype=np.int8) for k in tb.cb_k_values]
            symbols = qam_modulate(np.concatenate(adapter.encode(payload)), int(mcs.qm))
            realization = generate_sionna_tdl_channel_active(
                grid,
                channel,
                n_tx=8,
                n_rx=1,
                batch_size=batch,
                seed=stable_seed(config["seed"], SCENARIO, snr_db, absolute_start, "channel"),
                time_sample_indices=requested_times,
            )
            old, current, effective, aged_mrt = _effective_channels(
                realization.H, grid, precoders, candidates
            )
            if not np.allclose(np.sum(np.abs(aged_mrt) ** 2, axis=2), 8.0, atol=1e-10):
                raise RuntimeError("Aged MRT power normalization changed")
            old_current_cross += np.sum(old.conj() * current[:, :, :, 0, :])
            old_power += float(np.sum(np.abs(old) ** 2))
            current_power += float(np.sum(np.abs(current[:, :, :, 0, :]) ** 2))
            current_edge_cross += np.sum(
                current[:, :, :, 0, :].conj() * current[:, :, :, -1, :]
            )
            current_first_power += float(np.sum(np.abs(current[:, :, :, 0, :]) ** 2))
            current_last_power += float(np.sum(np.abs(current[:, :, :, -1, :]) ** 2))
            noise_rng = np.random.default_rng(
                stable_seed(config["seed"], SCENARIO, snr_db, absolute_start, "noise")
            )
            ls_noise = math.sqrt(noise_variance / 2.0) * (
                noise_rng.normal(size=(batch, 1, grid.n_dmrs_re))
                + 1j * noise_rng.normal(size=(batch, 1, grid.n_dmrs_re))
            )
            data_noise = math.sqrt(noise_variance / 2.0) * (
                noise_rng.normal(size=(batch, 1, grid.n_data_re))
                + 1j * noise_rng.normal(size=(batch, 1, grid.n_data_re))
            )
            decode_batches = []
            for candidate_id in candidate_ids:
                true = effective[candidate_id]
                true_pilot = true[:, :, grid.pilot_symbol_indices, pilot_local]
                true_data = true[:, :, grid.data_symbol_indices, data_local]
                estimate_data = filters[candidate_id].estimate_data(true_pilot + ls_noise)
                _, _, ratios = _ce_trial_statistics(estimate_data, true_data)
                ce_nmse[candidate_id][offset : offset + batch] = ratios
                received = true_data * symbols[None, None, :] + data_noise
                equalized, effective_noise = _estimated_csi_mrc_equalize(
                    estimate_data, received, noise_variance
                )
                decode_batches.append((candidate_id, equalized, effective_noise))
            combined_flags = _decode_flags(
                adapter,
                np.concatenate([item[1] for item in decode_batches], axis=0),
                np.concatenate([item[2] for item in decode_batches], axis=0),
                int(mcs.qm),
                payload,
            )
            cursor = 0
            for candidate_id, _, _ in decode_batches:
                flags[candidate_id][offset : offset + batch] = combined_flags[cursor : cursor + batch]
                cursor += batch
            completed = absolute_start + batch - 1
            if absolute_start == trial_start or completed % 100 < batch or completed == trial_end:
                print(
                    f"[{config['run_kind']} {SCENARIO}] snr={snr_db:g} "
                    f"absolute_trials={completed}/{trial_end} elapsed={time.time()-started:.1f}s",
                    flush=True,
                )
        interval_rows = []
        for candidate_id in candidate_ids:
            flag_path = _array_path(
                output, "error_flags", snr_db, candidate_id, trial_start, trial_end
            )
            nmse_path = _array_path(
                output, "ce_nmse_trials", snr_db, candidate_id, trial_start, trial_end
            )
            flag_path.parent.mkdir(parents=True, exist_ok=True)
            nmse_path.parent.mkdir(parents=True, exist_ok=True)
            np.save(flag_path, flags[candidate_id])
            np.save(nmse_path, ce_nmse[candidate_id])
            errors = int(np.sum(flags[candidate_id]))
            low, high = wilson(errors, trial_count)
            interval_rows.append(
                {
                    "scenario_id": SCENARIO,
                    "candidate_id": candidate_id,
                    "family": by_id[candidate_id]["family"],
                    "receiver": "estimated",
                    "snr_db": snr_db,
                    "trial_start": trial_start,
                    "trial_end": trial_end,
                    "trials": trial_count,
                    "tb_errors": errors,
                    "bler": errors / trial_count,
                    "bler_wilson95_lo": low,
                    "bler_wilson95_hi": high,
                    "ce_nmse_sum": float(np.sum(ce_nmse[candidate_id])),
                    "ce_nmse_sumsq": float(np.sum(ce_nmse[candidate_id] ** 2)),
                    "error_flags": _repo_relative(flag_path),
                    "ce_nmse_trial_path": _repo_relative(nmse_path),
                    "seed": config["seed"],
                }
            )
        stored.extend(interval_rows)
        stored.sort(
            key=lambda row: (
                str(row["candidate_id"]),
                float(row["snr_db"]),
                int(row["trial_start"]),
            )
        )
        _write_rows_atomic(stored, intervals_path)
        correlation = {
            "snr_db": snr_db,
            "trial_start": trial_start,
            "trial_end": trial_end,
            "old_to_current_first_complex_correlation": [
                float(np.real(old_current_cross / math.sqrt(old_power * current_power))),
                float(np.imag(old_current_cross / math.sqrt(old_power * current_power))),
            ],
            "current_symbol0_to9_complex_correlation": [
                float(np.real(current_edge_cross / math.sqrt(current_first_power * current_last_power))),
                float(np.imag(current_edge_cross / math.sqrt(current_first_power * current_last_power))),
            ],
        }
        progress = output / "channel_correlation_intervals.jsonl"
        with progress.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(correlation, ensure_ascii=False) + "\n")
        _merge_intervals(output, candidates)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--stage", choices=("validate", "run", "merge"), default="run")
    args = parser.parse_args()
    config = load_config(args.config.resolve())
    manifest, digest, candidates = validate_config(config)
    if args.stage == "run":
        run(config, manifest, digest, candidates)
    elif args.stage == "merge":
        _merge_intervals(Path(config["output_dir"]), candidates)


if __name__ == "__main__":
    main()
