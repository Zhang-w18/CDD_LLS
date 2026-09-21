from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
from typing import Any

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ChannelConfig  # noqa: E402
from cdd_lls.phy.channel_tdl import generate_sionna_tdl_channel_active  # noqa: E402
from cdd_lls.phy.pdcch import PDCCHResourceConfig, build_pdcch_grid  # noqa: E402
from cdd_lls.sim.pdcch import _stable_seed  # noqa: E402
from cdd_lls.sim.rsrp import (  # noqa: E402
    noiseless_block_average_power,
    noiseless_block_average_power_per_rx,
)


SCHEMA = "plan034-rsrp-cdf-v1"
PER_RX_SCHEMA = "plan034-rsrp-per-rx-4waveform-v1"
CANDIDATES = ("FIXED_DFT0", "FREQ_SIDON_0137")
PER_RX_CANDIDATES = ("FIXED_DFT0", "FREQ_SIDON_0137", "CDD911", "CDD130")
PER_RX_METADATA = {
    "FIXED_DFT0": {
        "delay_grid_coordinates": [0.0, 0.0, 0.0, 0.0],
        "delay_ns": [0.0, 0.0, 0.0, 0.0],
        "phase_denominator": None,
        "allow_duplicate_delays": False,
    },
    "FREQ_SIDON_0137": {
        "delay_grid_coordinates": [0.0, 1.0, 3.0, 7.0],
        "delay_ns": [0.0, 925.925926, 2777.777778, 6481.481481],
        "phase_denominator": 36,
        "allow_duplicate_delays": False,
    },
    "CDD911": {
        "delay_grid_coordinates": [0.0, 0.0, 0.98388, 0.98388],
        "delay_ns": [0.0, 0.0, 911.0, 911.0],
        "phase_denominator": 36,
        "allow_duplicate_delays": True,
    },
    "CDD130": {
        "delay_grid_coordinates": [0.0, 0.0, 0.1404, 0.1404],
        "delay_ns": [0.0, 0.0, 130.0, 130.0],
        "phase_denominator": 36,
        "allow_duplicate_delays": True,
    },
}
FIELDNAMES = (
    "absolute_trial",
    "trial_key",
    "fixed_block_power_linear",
    "fixed_block_power_db",
    "fixed_single_re_power_linear",
    "fixed_single_re_power_db",
    "freq_block_power_linear",
    "freq_block_power_db",
    "freq_single_re_power_linear",
    "freq_single_re_power_db",
)
PER_RX_FIELDNAMES = (
    "absolute_trial",
    "trial_key",
    "channel_key",
    "rx_index",
    "candidate_id",
    "block_power_linear",
    "block_power_db",
    "delay_grid_coordinates",
    "delay_ns",
    "phase_denominator",
    "allow_duplicate_delays",
    "total_transmit_power",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _resource(config: dict[str, Any]) -> PDCCHResourceConfig:
    fields = PDCCHResourceConfig.__dataclass_fields__
    return PDCCHResourceConfig(
        **{key: value for key, value in config["resource"].items() if key in fields}
    )


def _channel(config: dict[str, Any]) -> ChannelConfig:
    fields = ChannelConfig.__dataclass_fields__
    return ChannelConfig(
        **{key: value for key, value in config["channel"].items() if key in fields}
    )


def load_config(path: str | Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    validate_config(config)
    return config


def build_precoders(n_sc: int, n_tx: int, phase_denominator: int) -> dict[str, np.ndarray]:
    if int(n_tx) != 4 or int(phase_denominator) != 36:
        raise ValueError("Plan-034 freezes n_tx=4 and phase_denominator=36.")
    fixed = np.ones((int(n_sc), 4), dtype=np.complex128) / 2.0
    local = np.arange(int(n_sc), dtype=np.float64)
    coordinates = np.asarray([0.0, 1.0, 3.0, 7.0], dtype=np.float64)
    freq = np.exp(
        -1j * 2.0 * np.pi * local[:, None] * coordinates[None, :] / 36.0
    ) / 2.0
    return {"FIXED_DFT0": fixed, "FREQ_SIDON_0137": freq}


def build_per_rx_precoders(
    n_sc: int, n_tx: int, phase_denominator: int
) -> dict[str, np.ndarray]:
    """Build the four frozen waveforms for the Plan-034 per-Rx supplement."""

    if int(n_tx) != 4 or int(phase_denominator) != 36 or int(n_sc) != 36:
        raise ValueError("Plan-034 per-Rx supplement freezes n_sc=36, n_tx=4, denominator=36.")
    local = np.arange(36, dtype=np.float64)
    out: dict[str, np.ndarray] = {}
    for candidate_id in PER_RX_CANDIDATES:
        coordinates = np.asarray(
            PER_RX_METADATA[candidate_id]["delay_grid_coordinates"], dtype=np.float64
        )
        out[candidate_id] = np.exp(
            -1j * 2.0 * np.pi * local[:, None] * coordinates[None, :] / 36.0
        ) / 2.0
    return out


def _validate_common_config(config: dict[str, Any]) -> tuple[Any, np.ndarray]:
    """Validate the frozen physical geometry shared by both Plan-034 RSRP schemas."""

    grid = build_pdcch_grid(_resource(config))
    channel = _channel(config)
    if config.get("antenna") != {"n_tx": 4, "n_rx": 4}:
        raise ValueError("Plan-034 RSRP requires exactly 4Tx/4Rx.")
    if (
        int(grid.config.n_rb) != 48
        or int(grid.config.duration_symbols) != 2
        or int(grid.config.aggregation_level) != 1
        or int(grid.config.first_cce) != 0
        or str(grid.config.cce_reg_mapping).lower() != "noninterleaved"
        or int(grid.config.reg_bundle_size) != 6
    ):
        raise ValueError("Plan-034 RSRP requires 48-RB, 2-symbol, AL1, first-CCE-0, L=6.")
    if grid.resource_grid.n_data_re != 54 or grid.resource_grid.n_dmrs_re != 18:
        raise ValueError("Unexpected AL1 PDCCH resource counts.")
    occupied = np.unique(
        np.concatenate(
            [grid.data_local_subcarrier_indices, grid.dmrs_local_subcarrier_indices]
        )
    )
    if occupied.tolist() != list(range(36)):
        raise ValueError("Plan-034 requires the first contiguous 36-subcarrier candidate band.")
    if not (
        str(channel.backend).lower() == "sionna_tdl"
        and str(channel.tdl_profile).upper() == "C"
        and abs(float(channel.delay_spread_ns) - 300.0) <= 1e-12
        and abs(float(channel.carrier_frequency_hz) - 4.0e9) <= 1e-3
        and abs(float(channel.ue_speed_kmh) - 3.0) <= 1e-12
        and int(channel.num_sinusoids) == 20
        and bool(channel.normalize)
    ):
        raise ValueError("Plan-034 requires normalized Sionna TDL-C 300 ns, 4 GHz, 3 km/h.")
    simulation = config["simulation"]
    trials = int(simulation["trials"])
    batch_size = int(simulation["batch_size"])
    smoke_trials = int(simulation.get("smoke_trials", 100))
    if trials <= 0 or batch_size <= 0 or smoke_trials <= 0:
        raise ValueError("trials, batch_size, and smoke_trials must be positive.")
    if int(config["seed"]) != 20260916:
        raise ValueError("Plan-034 freezes seed=20260916.")
    if str(config["random_stream_namespace"]) != "plan034_rsrp_cdf_v1":
        raise ValueError("Unexpected Plan-034 random stream namespace.")
    return grid, occupied


def _validate_original_config(config: dict[str, Any]) -> dict[str, Any]:
    if str(config.get("schema")) != SCHEMA:
        raise ValueError(f"schema must be {SCHEMA}.")
    grid, occupied = _validate_common_config(config)
    precoders = build_precoders(36, 4, int(config["precoding"]["phase_denominator"]))
    row_power = {
        name: np.sum(np.abs(matrix) ** 2, axis=1) for name, matrix in precoders.items()
    }
    if any(not np.allclose(value, 1.0, rtol=0.0, atol=1e-12) for value in row_power.values()):
        raise RuntimeError("Precoder row-power validation failed.")
    return {
        "candidate_ids": list(CANDIDATES),
        "candidate_cces": grid.candidate_cces.tolist(),
        "candidate_bundles": grid.candidate_bundles.tolist(),
        "candidate_re": 72,
        "n_data_re": 54,
        "n_dmrs_re": 18,
        "occupied_subcarriers": occupied.tolist(),
        "single_re": {"symbol": 0, "local_subcarrier": 18},
        "phase_denominator": 36,
        "delay_grid_coordinates": [0, 1, 3, 7],
        "q_seconds": 1.0 / (36.0 * 30000.0),
        "row_power_min_max": {
            name: [float(np.min(value)), float(np.max(value))]
            for name, value in row_power.items()
        },
    }


def _validate_per_rx_config(config: dict[str, Any]) -> dict[str, Any]:
    if str(config.get("schema")) != PER_RX_SCHEMA:
        raise ValueError(f"schema must be {PER_RX_SCHEMA}.")
    grid, occupied = _validate_common_config(config)
    precoding = config.get("precoding", {})
    if int(precoding.get("phase_denominator", -1)) != 36:
        raise ValueError("Per-Rx supplement freezes phase_denominator=36.")
    if float(precoding.get("total_transmit_power", -1.0)) != 1.0:
        raise ValueError("Per-Rx supplement freezes total_transmit_power=1.0.")
    configured = precoding.get("candidates")
    expected = [
        {"candidate_id": candidate_id, **PER_RX_METADATA[candidate_id]}
        for candidate_id in PER_RX_CANDIDATES
    ]
    if configured != expected:
        raise ValueError("Per-Rx candidate metadata differs from the frozen four-candidate table.")
    precoders = build_per_rx_precoders(36, 4, 36)
    row_power = {
        name: np.sum(np.abs(matrix) ** 2, axis=1) for name, matrix in precoders.items()
    }
    if any(not np.allclose(value, 1.0, rtol=0.0, atol=1e-12) for value in row_power.values()):
        raise RuntimeError("Per-Rx precoder row-power validation failed.")
    return {
        "schema": PER_RX_SCHEMA,
        "candidate_ids": list(PER_RX_CANDIDATES),
        "candidate_metadata": expected,
        "candidate_cces": grid.candidate_cces.tolist(),
        "candidate_bundles": grid.candidate_bundles.tolist(),
        "candidate_re": 72,
        "n_data_re": 54,
        "n_dmrs_re": 18,
        "occupied_subcarriers": occupied.tolist(),
        "rx_indices": [0, 1, 2, 3],
        "per_rx_statistic": True,
        "rx_average_applied": False,
        "q_seconds": 1.0 / (36.0 * 30000.0),
        "row_power_min_max": {
            name: [float(np.min(value)), float(np.max(value))]
            for name, value in row_power.items()
        },
    }


def validate_config(config: dict[str, Any]) -> dict[str, Any]:
    schema = str(config.get("schema"))
    if schema == SCHEMA:
        return _validate_original_config(config)
    if schema == PER_RX_SCHEMA:
        return _validate_per_rx_config(config)
    raise ValueError(f"schema must be one of: {SCHEMA}, {PER_RX_SCHEMA}.")


def _completed_trials(path: Path) -> int:
    if not path.exists():
        return 0
    last = 0
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != FIELDNAMES:
            raise RuntimeError("Existing RSRP CSV fields do not match Plan-034 schema.")
        for expected, row in enumerate(reader, start=1):
            actual = int(row["absolute_trial"])
            if actual != expected:
                raise RuntimeError(
                    f"RSRP resume requires contiguous trials; expected {expected}, found {actual}."
                )
            last = actual
    return last


def _completed_per_rx_trials(path: Path) -> int:
    if not path.exists():
        return 0
    rows_per_trial = len(PER_RX_CANDIDATES) * 4
    row_count = 0
    last_trial = 0
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != PER_RX_FIELDNAMES:
            raise RuntimeError("Existing per-Rx RSRP CSV fields do not match the supplement schema.")
        for row_index, row in enumerate(reader):
            trial = row_index // rows_per_trial + 1
            position = row_index % rows_per_trial
            candidate_index = position // 4
            rx_index = position % 4
            if int(row["absolute_trial"]) != trial:
                raise RuntimeError(f"Per-Rx resume expected trial {trial} at row {row_index + 2}.")
            if row["candidate_id"] != PER_RX_CANDIDATES[candidate_index]:
                raise RuntimeError(f"Per-Rx resume candidate order mismatch at row {row_index + 2}.")
            if int(row["rx_index"]) != rx_index:
                raise RuntimeError(f"Per-Rx resume Rx order mismatch at row {row_index + 2}.")
            row_count += 1
            last_trial = trial
    if row_count % rows_per_trial:
        raise RuntimeError("Existing per-Rx RSRP CSV ends with an incomplete trial group.")
    return last_trial


def _trial_key(config: dict[str, Any], absolute_trial: int) -> str:
    return f"{_stable_seed(config['seed'], config['random_stream_namespace'], absolute_trial):08x}"


def _append_rows(path: Path, rows: list[dict[str, object]]) -> None:
    exists = path.exists()
    with open(path, "a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        if not exists:
            writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())


def _append_per_rx_rows(path: Path, rows: list[dict[str, object]]) -> None:
    exists = path.exists()
    with open(path, "a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PER_RX_FIELDNAMES)
        if not exists:
            writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())


def _run(config: dict[str, Any], config_path: Path, *, target_trials: int, output: Path) -> None:
    receipt = validate_config(config)
    output.mkdir(parents=True, exist_ok=True)
    configs_dir = output.parent / "configs"
    validation_dir = output.parent / "validation"
    configs_dir.mkdir(parents=True, exist_ok=True)
    validation_dir.mkdir(parents=True, exist_ok=True)
    destination = configs_dir / config_path.name
    shutil.copyfile(config_path, destination)
    (configs_dir / f"{config_path.stem}.sha256").write_text(
        _sha256(config_path) + "\n", encoding="utf-8"
    )
    with open(configs_dir / "pdcch_plan034_rsrp_cdf.resolved.yaml", "w", encoding="utf-8") as handle:
        yaml.safe_dump(config, handle, allow_unicode=True, sort_keys=False)
    (validation_dir / "rsrp_validation.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    resource = _resource(config)
    grid = build_pdcch_grid(resource)
    occupied = np.arange(36, dtype=np.int64)
    band_grid = grid.resource_grid
    from dataclasses import replace

    band_grid = replace(
        band_grid,
        n_sc=36,
        subcarrier_indices=band_grid.subcarrier_indices[occupied],
        active_fft_indices=band_grid.active_fft_indices[occupied],
        pilot_subcarriers=np.unique(band_grid.pilot_subcarrier_indices),
    )
    symbols = np.repeat(np.arange(2, dtype=np.int64), 36)
    subcarriers = np.tile(np.arange(36, dtype=np.int64), 2)
    single_symbol = np.asarray([0], dtype=np.int64)
    single_subcarrier = np.asarray([18], dtype=np.int64)
    precoders = build_precoders(36, 4, 36)
    csv_path = output / "paired_rsrp_trials.csv"
    completed = _completed_trials(csv_path) if bool(config["simulation"].get("resume", True)) else 0
    if completed > target_trials:
        raise RuntimeError(f"Existing RSRP output has {completed} trials, above target {target_trials}.")
    batch_size = int(config["simulation"]["batch_size"])
    started = time.time()
    while completed < target_trials:
        current = min(batch_size, target_trials - completed)
        absolute_start = completed + 1
        channel_seed = _stable_seed(
            config["seed"],
            config["random_stream_namespace"],
            absolute_start,
            current,
            "paired_channel_batch",
        )
        channel = generate_sionna_tdl_channel_active(
            band_grid,
            _channel(config),
            n_tx=4,
            n_rx=4,
            batch_size=current,
            seed=channel_seed,
        ).H
        values: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for candidate_id, precoder in precoders.items():
            block = noiseless_block_average_power(
                channel,
                precoder,
                symbols,
                subcarriers,
                n_rx=4,
                transmit_power=1.0,
            )
            single = noiseless_block_average_power(
                channel,
                precoder,
                single_symbol,
                single_subcarrier,
                n_rx=4,
                transmit_power=1.0,
            )
            if np.any(~np.isfinite(block)) or np.any(block <= 0.0):
                raise RuntimeError(f"{candidate_id}: non-finite or non-positive block power.")
            if np.any(~np.isfinite(single)) or np.any(single <= 0.0):
                raise RuntimeError(f"{candidate_id}: non-finite or non-positive single-RE power.")
            values[candidate_id] = (block, single)
        rows = []
        for offset in range(current):
            absolute_trial = absolute_start + offset
            fixed_block, fixed_single = values["FIXED_DFT0"]
            freq_block, freq_single = values["FREQ_SIDON_0137"]
            rows.append(
                {
                    "absolute_trial": absolute_trial,
                    "trial_key": _trial_key(config, absolute_trial),
                    "fixed_block_power_linear": float(fixed_block[offset]),
                    "fixed_block_power_db": float(10.0 * np.log10(fixed_block[offset])),
                    "fixed_single_re_power_linear": float(fixed_single[offset]),
                    "fixed_single_re_power_db": float(10.0 * np.log10(fixed_single[offset])),
                    "freq_block_power_linear": float(freq_block[offset]),
                    "freq_block_power_db": float(10.0 * np.log10(freq_block[offset])),
                    "freq_single_re_power_linear": float(freq_single[offset]),
                    "freq_single_re_power_db": float(10.0 * np.log10(freq_single[offset])),
                }
            )
        _append_rows(csv_path, rows)
        completed += current
        print(f"RSRP progress {completed}/{target_trials}", flush=True)

    receipt = {
        "schema": "plan034-rsrp-run-receipt-v1",
        "trials": completed,
        "candidate_pairing": "same channel tensor per absolute trial row",
        "noise": "none",
        "per_realization_normalization": False,
        "elapsed_seconds_this_invocation": time.time() - started,
        "csv_sha256": _sha256(csv_path),
    }
    (output / "run_receipt.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def _run_per_rx(
    config: dict[str, Any], config_path: Path, *, target_trials: int, output: Path
) -> None:
    receipt = validate_config(config)
    output.mkdir(parents=True, exist_ok=True)
    configs_dir = output.parent / "configs"
    validation_dir = output.parent / "validation"
    configs_dir.mkdir(parents=True, exist_ok=True)
    validation_dir.mkdir(parents=True, exist_ok=True)
    destination = configs_dir / config_path.name
    shutil.copyfile(config_path, destination)
    (configs_dir / f"{config_path.stem}.sha256").write_text(
        _sha256(config_path) + "\n", encoding="utf-8"
    )
    with open(
        configs_dir / "pdcch_plan034_rsrp_per_rx_4waveform.resolved.yaml",
        "w",
        encoding="utf-8",
    ) as handle:
        yaml.safe_dump(config, handle, allow_unicode=True, sort_keys=False)
    (validation_dir / "per_rx_rsrp_validation.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    grid = build_pdcch_grid(_resource(config))
    occupied = np.arange(36, dtype=np.int64)
    from dataclasses import replace

    band_grid = replace(
        grid.resource_grid,
        n_sc=36,
        subcarrier_indices=grid.resource_grid.subcarrier_indices[occupied],
        active_fft_indices=grid.resource_grid.active_fft_indices[occupied],
        pilot_subcarriers=np.unique(grid.resource_grid.pilot_subcarrier_indices),
    )
    symbols = np.repeat(np.arange(2, dtype=np.int64), 36)
    subcarriers = np.tile(np.arange(36, dtype=np.int64), 2)
    precoders = build_per_rx_precoders(36, 4, 36)
    csv_path = output / "per_rx_rsrp_trials.csv"
    completed = (
        _completed_per_rx_trials(csv_path)
        if bool(config["simulation"].get("resume", True))
        else 0
    )
    if completed > target_trials:
        raise RuntimeError(
            f"Existing per-Rx RSRP output has {completed} trials, above target {target_trials}."
        )
    batch_size = int(config["simulation"]["batch_size"])
    started = time.time()
    maximum_linear_reconstruction_difference = 0.0
    while completed < target_trials:
        current = min(batch_size, target_trials - completed)
        absolute_start = completed + 1
        channel_seed = _stable_seed(
            config["seed"],
            config["random_stream_namespace"],
            absolute_start,
            current,
            "paired_channel_batch",
        )
        channel = generate_sionna_tdl_channel_active(
            band_grid,
            _channel(config),
            n_tx=4,
            n_rx=4,
            batch_size=current,
            seed=channel_seed,
        ).H
        values: dict[str, np.ndarray] = {}
        for candidate_id, precoder in precoders.items():
            power = noiseless_block_average_power_per_rx(
                channel,
                precoder,
                symbols,
                subcarriers,
                n_rx=4,
                transmit_power=1.0,
            )
            if power.shape != (current, 4):
                raise RuntimeError(f"{candidate_id}: unexpected per-Rx power shape {power.shape}.")
            if np.any(~np.isfinite(power)) or np.any(power <= 0.0):
                raise RuntimeError(f"{candidate_id}: non-finite or non-positive per-Rx power.")
            aggregate = noiseless_block_average_power(
                channel,
                precoder,
                symbols,
                subcarriers,
                n_rx=4,
                transmit_power=1.0,
            )
            difference = float(np.max(np.abs(np.mean(power, axis=1) - aggregate)))
            maximum_linear_reconstruction_difference = max(
                maximum_linear_reconstruction_difference, difference
            )
            if not np.allclose(np.mean(power, axis=1), aggregate, rtol=0.0, atol=1e-12):
                raise RuntimeError(
                    f"{candidate_id}: per-Rx powers do not reconstruct the Rx-average statistic."
                )
            values[candidate_id] = power

        rows: list[dict[str, object]] = []
        for offset in range(current):
            absolute_trial = absolute_start + offset
            trial_key = _trial_key(config, absolute_trial)
            for candidate_id in PER_RX_CANDIDATES:
                metadata = PER_RX_METADATA[candidate_id]
                delay_grid = json.dumps(
                    metadata["delay_grid_coordinates"], separators=(",", ":")
                )
                delay_ns = json.dumps(metadata["delay_ns"], separators=(",", ":"))
                phase_denominator = metadata["phase_denominator"]
                for rx_index in range(4):
                    linear = float(values[candidate_id][offset, rx_index])
                    rows.append(
                        {
                            "absolute_trial": absolute_trial,
                            "trial_key": trial_key,
                            "channel_key": trial_key,
                            "rx_index": rx_index,
                            "candidate_id": candidate_id,
                            "block_power_linear": linear,
                            "block_power_db": float(10.0 * np.log10(linear)),
                            "delay_grid_coordinates": delay_grid,
                            "delay_ns": delay_ns,
                            "phase_denominator": (
                                "" if phase_denominator is None else int(phase_denominator)
                            ),
                            "allow_duplicate_delays": str(
                                bool(metadata["allow_duplicate_delays"])
                            ).lower(),
                            "total_transmit_power": 1.0,
                        }
                    )
        _append_per_rx_rows(csv_path, rows)
        completed += current
        print(f"Per-Rx RSRP progress {completed}/{target_trials}", flush=True)

    run_receipt = {
        "schema": "plan034-per-rx-rsrp-run-receipt-v1",
        "trials": completed,
        "rows": completed * len(PER_RX_CANDIDATES) * 4,
        "candidate_ids": list(PER_RX_CANDIDATES),
        "rx_indices": [0, 1, 2, 3],
        "candidate_pairing": "same channel tensor per absolute trial and Rx",
        "rx_averaging": False,
        "per_rx_mean_reconstructs_legacy_average": True,
        "maximum_linear_reconstruction_difference": maximum_linear_reconstruction_difference,
        "noise": "none",
        "per_realization_normalization": False,
        "elapsed_seconds_this_invocation": time.time() - started,
        "csv_sha256": _sha256(csv_path),
    }
    (output / "run_receipt.json").write_text(
        json.dumps(run_receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def _environment_receipt(root: Path) -> None:
    try:
        import scipy
        import tensorflow as tf
        import sionna

        versions = {
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "tensorflow": tf.__version__,
            "sionna": getattr(sionna, "__version__", "unknown"),
        }
    except Exception as exc:  # pragma: no cover - diagnostic fallback
        versions = {"numpy": np.__version__, "dependency_probe_error": repr(exc)}
    peak = None
    if os.name != "nt":
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    try:
        git_head = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        git_dirty = bool(
            subprocess.run(
                ["git", "status", "--porcelain"],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
        )
    except Exception as exc:  # pragma: no cover - diagnostic fallback
        git_head = "unavailable"
        git_dirty = f"probe_failed:{exc!r}"
    payload = {
        "python": sys.version,
        "platform": platform.platform(),
        "versions": versions,
        "peak_rss_platform_units": peak,
        "git_head": git_head,
        "git_worktree_dirty": git_dirty,
    }
    root.mkdir(parents=True, exist_ok=True)
    (root / "environment_receipt.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Plan-034 paired noiseless RSRP statistics.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--stage", choices=("validate", "smoke", "run"), default="run")
    args = parser.parse_args()
    config_path = Path(args.config).resolve()
    config = load_config(config_path)
    summary = validate_config(config)
    if args.stage == "validate":
        print(json.dumps(summary, indent=2, ensure_ascii=False))
        return
    root = Path(config["output_dir"])
    target = (
        int(config["simulation"].get("smoke_trials", 100))
        if args.stage == "smoke"
        else int(config["simulation"]["trials"])
    )
    output = root / ("smoke" if args.stage == "smoke" else "rsrp")
    if str(config["schema"]) == PER_RX_SCHEMA:
        _run_per_rx(config, config_path, target_trials=target, output=output)
    else:
        _run(config, config_path, target_trials=target, output=output)
    _environment_receipt(root)


if __name__ == "__main__":
    main()
