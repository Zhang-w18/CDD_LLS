"""Run plan-033 4/8Tx, 4Rx TDL-A mobility curves with 2D RMMSE."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Sequence

import numpy as np
import psutil

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import yaml

from cdd_lls.core.config import ChannelConfig
from cdd_lls.core.mcs import build_tb_layout, get_mcs
from cdd_lls.phy.channel_tdl import generate_sionna_channel_active
from cdd_lls.phy.estimators import (
    TDLTimeFrequencyCovariance,
    build_prg_time_frequency_rmmse_filter,
    build_time_frequency_rmmse_filter,
    cdl_spatial_unaware_covariance,
    tdl_active_time_frequency_covariance,
)
from cdd_lls.phy.ldpc import SionnaLDPCAdapter
from cdd_lls.phy.precoding import (
    build_active_dft_grid_precoder,
    build_aged_mrt_prg_precoder,
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
from tools.run_plan025_delay_matched_tdl import read_csv_rows, stable_seed, wilson, write_csv_rows
from tools.run_plan032_tdl_mobility import _repo_relative, _safe, _sha256, _snr_key


SCHEMA = "plan033-tdl-mobility-mimo-v1"
N_RX = 4
PRG_SIZE_RB = 6
HISTORY_SYMBOL_INDEX = 140
FEEDBACK_PERIOD_SLOTS = 10
PRECODER_NORMALIZE = True
TRANSPARENT_COVARIANCE_SCALE = 1.0


def scenario_id(n_tx: int, speed_kmh: float, n_rx: int = N_RX) -> str:
    speed = int(speed_kmh) if float(speed_kmh).is_integer() else float(speed_kmh)
    return f"A100_NT{int(n_tx)}_NR{int(n_rx)}_V{speed}"


def configured_scenario_id(config: dict) -> str:
    base = scenario_id(config["n_tx"], config["speed_kmh"], config["n_rx"])
    return base.replace("A100_", "D100_", 1) if config["channel_model"] == "cdl_d" else base


def awgn_variance(snr_db: float) -> float:
    """Return per-Rx AWGN/LS-observation variance for unit transmit power."""
    return 1.0 / (10.0 ** (float(snr_db) / 10.0))


def delay_sets(n_tx: int) -> dict[str, list[float]]:
    if int(n_tx) == 4:
        return {
            "B0_QC": [0.0, 24.0, 48.0, 72.0],
            "S0_SIDON": [0.0, 1.0, 3.0, 7.0],
            "SMALL_CDD_QSTEP0P25": [0.0, 0.25, 0.5, 0.75],
        }
    if int(n_tx) == 8:
        return {
            "B0_QC": [0.0, 9.0, 18.0, 27.0, 36.0, 45.0, 54.0, 63.0],
            "S0_SIDON": [0.0, 1.0, 3.0, 7.0, 12.0, 20.0, 30.0, 65.0],
            "SMALL_CDD_QSTEP0P25": [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75],
        }
    raise ValueError("Plan-033 n_tx must be 4 or 8.")


def candidate_ids(n_tx: int) -> dict[str, str]:
    prefix = f"A100_NT{int(n_tx)}"
    return {
        "b0": f"{prefix}_B0_QC",
        "sidon": f"{prefix}_S0_SIDON",
        "prg": f"{prefix}_TRANSPARENT_PRG6",
        "mrt": f"{prefix}_AGED_MRT_PRG6",
        "small_transparent": f"{prefix}_SMALL_CDD_QSTEP0P25_TRANSPARENT",
        "small_matched": f"{prefix}_SMALL_CDD_QSTEP0P25_MATCHED",
    }


def candidate_definitions(n_tx: int) -> list[dict]:
    ids = candidate_ids(n_tx)
    delays = delay_sets(n_tx)
    return [
        {"candidate_id": ids["b0"], "family": "CDD_MATCHED", "label": "B0_QC", "delay_grid_coordinates": delays["B0_QC"], "receiver_knowledge": "matched CDD covariance"},
        {"candidate_id": ids["sidon"], "family": "CDD_MATCHED", "label": "S0_SIDON", "delay_grid_coordinates": delays["S0_SIDON"], "receiver_knowledge": "matched CDD covariance"},
        {"candidate_id": ids["prg"], "family": "TRANSPARENT_PRG_DFT", "label": "transparent PRG6", "receiver_knowledge": "physical covariance within each PRG"},
        {"candidate_id": ids["mrt"], "family": "AGED_CSI_MRT", "label": "aged-CSI MRT PRG6", "receiver_knowledge": "physical covariance within each PRG"},
        {"candidate_id": ids["small_transparent"], "family": "SMALL_CDD_TRANSPARENT", "label": "small-delay CDD, transparent", "delay_grid_coordinates": delays["SMALL_CDD_QSTEP0P25"], "receiver_knowledge": "physical covariance; CDD delay unknown"},
        {"candidate_id": ids["small_matched"], "family": "SMALL_CDD_MATCHED", "label": "small-delay CDD, non-transparent", "delay_grid_coordinates": delays["SMALL_CDD_QSTEP0P25"], "receiver_knowledge": "matched CDD covariance"},
    ]


def _resolve(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def _write_rows_atomic(rows: Sequence[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    write_csv_rows(rows, temporary)
    temporary.replace(path)


def _array_path(output: Path, category: str, snr_db: float, candidate_id: str, trial_start: int, trial_end: int, receiver: str = "estimated") -> Path:
    return output / category / receiver / _snr_key(snr_db) / f"{_safe(candidate_id)}_t{trial_start:06d}_{trial_end:06d}.npy"


def _interval_identity(
    config: dict,
    candidate_id: str,
    snr_db: float,
    trial_start: int,
    trial_end: int,
    receiver: str = "estimated",
) -> dict:
    return {
        "scenario_id": config["scenario_id"],
        "n_tx": config["n_tx"],
        "n_rx": config["n_rx"],
        "speed_kmh": config["speed_kmh"],
        "candidate_id": candidate_id,
        "receiver": receiver,
        "trial_key": (
            f"{config['scenario_id']}|snr={float(snr_db):g}|"
            f"trials={int(trial_start)}:{int(trial_end)}"
        ),
    }


def load_config(path: Path, *, expected_schema: str = SCHEMA) -> dict:
    config = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if config.get("schema") != expected_schema:
        raise ValueError(f"config.schema must be {expected_schema!r}")
    config["config_path"] = str(path.resolve())
    config["output_dir"] = str(_resolve(config["output_dir"]))
    config["seed"] = int(config.get("seed", 20260727))
    config["batch_size"] = int(config.get("batch_size", 25))
    config["n_tx"] = int(config["n_tx"])
    if expected_schema != SCHEMA and "n_rx" not in config:
        raise ValueError("n_rx must be explicit in plan-035 configs.")
    config["n_rx"] = int(config.get("n_rx", N_RX))
    config["receiver_modes"] = [str(value) for value in config.get("receiver_modes", ["estimated"])]
    config["speed_kmh"] = float(config["speed_kmh"])
    config["channel_model"] = str(config.get("channel_model", "tdl_a")).lower()
    config["feedback_period_slots"] = int(config.get("feedback_period_slots", FEEDBACK_PERIOD_SLOTS))
    config["history_symbol_index"] = int(config.get("history_symbol_index", HISTORY_SYMBOL_INDEX))
    config["run_kind"] = str(config["run_kind"])
    if "precoder_normalize" not in config:
        raise ValueError("precoder_normalize: true must be explicit in every plan-033 config.")
    if config["precoder_normalize"] is not True:
        raise ValueError("precoder_normalize must be the YAML boolean true.")
    config["precoder_normalize"] = True
    config["scenario_id"] = str(config.get("scenario_id", scenario_id(config["n_tx"], config["speed_kmh"], config["n_rx"])))
    config["snr_db"] = [float(value) for value in config["snr_db"]]
    scalar_target = config.get("target_total_trials")
    per_snr_targets = config.get("target_total_trials_by_snr")
    if (scalar_target is None) == (per_snr_targets is None):
        raise ValueError("Exactly one trial target form is required.")
    if scalar_target is not None:
        config["target_total_trials"] = int(scalar_target)
    else:
        config["target_total_trials_by_snr"] = {float(key): int(value) for key, value in per_snr_targets.items()}
    policy = config.get("adaptive_policy")
    if policy is not None:
        config["adaptive_policy"] = {
            "minimum_trials": int(policy["minimum_trials"]),
            "check_interval_trials": int(policy["check_interval_trials"]),
            "maximum_trials": int(policy["maximum_trials"]),
            "target_errors_10pct": int(policy["target_errors_10pct"]),
            "target_errors_1pct": int(policy["target_errors_1pct"]),
            "ten_pct_error_band": [float(value) for value in policy["ten_pct_error_band"]],
            "one_pct_error_band": [float(value) for value in policy["one_pct_error_band"]],
        }
    return config


def _trial_target(config: dict, snr_db: float) -> int:
    if "target_total_trials" in config:
        return int(config["target_total_trials"])
    return int(config["target_total_trials_by_snr"][float(snr_db)])


def _next_interval_end(current_total: int, target: int, run_kind: str) -> int:
    """Return the inclusive end of the next recoverable trial interval."""
    if current_total >= target:
        return current_total
    interval_limit = 1000 if run_kind == "formal" else target
    return min(target, current_total + interval_limit)


def _make_error_flags(
    candidate_ids_: Sequence[str], receiver_modes: Sequence[str], trial_count: int
) -> dict[tuple[str, str], np.ndarray]:
    """Allocate one independent error-flag vector per candidate/receiver pair."""
    return {
        (candidate_id, receiver): np.zeros(int(trial_count), dtype=bool)
        for candidate_id in candidate_ids_
        for receiver in receiver_modes
    }


def _adaptive_schedule(config: dict, candidates: Sequence[dict]) -> tuple[dict[float, int], dict]:
    """Return the next common per-SNR trial targets and an auditable decision payload."""
    policy = config["adaptive_policy"]
    output = Path(config["output_dir"])
    path = output / "intervals.csv"
    rows = [dict(row) for row in read_csv_rows(path)] if path.exists() else []
    ids = [str(row["candidate_id"]) for row in candidates]
    curves = [(candidate_id, receiver) for candidate_id in ids for receiver in config["receiver_modes"]]
    grouped: dict[tuple[str, str, float], dict[str, int]] = {}
    for row in rows:
        key = (str(row["candidate_id"]), str(row["receiver"]), float(row["snr_db"]))
        state = grouped.setdefault(key, {"trials": 0, "errors": 0})
        state["trials"] += int(row["trials"])
        state["errors"] += int(row["tb_errors"])
    current_by_snr = {}
    for snr_db in config["snr_db"]:
        counts = {grouped.get((candidate_id, receiver, snr_db), {"trials": 0})["trials"] for candidate_id, receiver in curves}
        if len(counts) != 1:
            raise RuntimeError(f"Adaptive paired trial counts differ at {snr_db:g} dB.")
        current_by_snr[snr_db] = counts.pop()
    minimum = policy["minimum_trials"]
    interval = policy["check_interval_trials"]
    maximum = policy["maximum_trials"]
    if any(value < minimum for value in current_by_snr.values()):
        targets = {snr: min(minimum, current + interval) if current < minimum else current for snr, current in current_by_snr.items()}
        return targets, {"status": "collecting_minimum", "targets": targets, "current_trials": current_by_snr}
    endpoint_reasons: dict[float, list[dict]] = defaultdict(list)
    missing = []
    for candidate_id, receiver in curves:
        points = []
        for snr_db in config["snr_db"]:
            state = grouped[(candidate_id, receiver, snr_db)]
            points.append((snr_db, state["errors"] / state["trials"], state["errors"], state["trials"]))
        for level, error_key in ((0.10, "target_errors_10pct"), (0.01, "target_errors_1pct")):
            brackets = []
            for left, right in zip(points, points[1:]):
                if left[1] >= level and right[1] <= level and (left[1] > level or right[1] < level):
                    brackets.append((left, right))
            if not brackets:
                missing.append({"candidate_id": candidate_id, "receiver": receiver, "target_bler": level})
                continue
            for left, right in brackets:
                for snr_db, bler, errors, trials in (left, right):
                    band = policy["ten_pct_error_band"] if level == 0.10 else policy["one_pct_error_band"]
                    require = band[0] <= bler <= band[1]
                    if require:
                        endpoint_reasons[snr_db].append({"candidate_id": candidate_id, "receiver": receiver, "target_bler": level, "observed_bler": bler, "errors": errors, "required_errors": policy[error_key]})
    if missing:
        return current_by_snr, {"status": "missing_bracket", "missing": missing, "current_trials": current_by_snr}
    targets = dict(current_by_snr)
    capped = []
    for snr_db, reasons in endpoint_reasons.items():
        current = current_by_snr[snr_db]
        unmet = [reason for reason in reasons if reason["errors"] < reason["required_errors"]]
        if unmet and current < maximum:
            targets[snr_db] = min(current + interval, maximum)
        elif unmet:
            capped.append({"snr_db": snr_db, "unmet": unmet})
    status = ("complete_with_caps" if capped else "complete") if targets == current_by_snr else "append_required"
    return targets, {"status": status, "targets": targets, "current_trials": current_by_snr, "endpoint_reasons": dict(endpoint_reasons), "maximum_trials": maximum, "capped_endpoints": capped}


def run_adaptive(config: dict, candidates: Sequence[dict], channel, grid, precoders, base) -> None:
    """Run recoverable staged intervals until the predeclared adaptive rule stops."""
    output = Path(config["output_dir"])
    while True:
        targets, audit = _adaptive_schedule(config, candidates)
        output.mkdir(parents=True, exist_ok=True)
        (output / "adaptive_status.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
        if audit["status"] in {"complete", "complete_with_caps"}:
            print(f"[adaptive] {audit['status']}", flush=True)
            return
        if audit["status"] == "missing_bracket":
            raise RuntimeError("Adaptive run stopped because at least one 10%/1% bracket is missing; inspect adaptive_status.json.")
        staged = dict(config)
        staged.pop("target_total_trials", None)
        staged["target_total_trials_by_snr"] = targets
        run(staged, candidates, channel, grid, precoders, base)


def _build_scene(config: dict):
    grid = build_resource_grid(bler027.resource_config(6))
    if config["channel_model"] == "tdl_a":
        channel = ChannelConfig(
            backend="sionna_tdl", model="3gpp_tr38901_tdl", tdl_profile="A",
            delay_spread_ns=100.0, carrier_frequency_hz=3.5e9,
            ue_speed_kmh=config["speed_kmh"], num_sinusoids=20,
        )
        base = tdl_active_time_frequency_covariance(grid, channel)
    elif config["channel_model"] == "cdl_d":
        channel = ChannelConfig(
            backend="sionna_cdl", model="3gpp_tr38901_cdl", cdl_profile="D",
            delay_spread_ns=100.0, carrier_frequency_hz=3.5e9,
            ue_speed_kmh=config["speed_kmh"], cdl_direction="downlink",
            cdl_tx_array_rows=1, cdl_tx_array_cols=config["n_tx"],
            cdl_rx_array_rows=1, cdl_rx_array_cols=config["n_rx"],
            cdl_polarization="single", cdl_polarization_type="V",
            cdl_antenna_pattern="omni", cdl_element_vertical_spacing=0.5,
            cdl_element_horizontal_spacing=0.5,
        )
        base = cdl_spatial_unaware_covariance(grid, channel)
    else:
        raise ValueError("channel_model must be tdl_a or cdl_d")
    n_tx = config["n_tx"]
    runtime_candidates = config.get("_candidate_manifest")
    delays = (
        {
            str(row["candidate_id"]): [float(value) for value in row["delay_grid_coordinates"]]
            for row in runtime_candidates
            if row["family"] == "CDD_MATCHED"
        }
        if runtime_candidates is not None
        else delay_sets(n_tx)
    )
    precoders = {
        name: build_active_dft_grid_precoder(
            grid, values, n_tx=n_tx, normalize=config["precoder_normalize"]
        )
        for name, values in delays.items()
    }
    return channel, grid, precoders, base


def _delay_audit(grid, n_tx: int, name: str, values: Sequence[float], *, modulo_pair_sums: bool = False) -> dict:
    delay = np.asarray(values, dtype=np.float64)
    pilot_local = np.unique(local_indices_for_subcarriers(grid, grid.pilot_subcarrier_indices))
    matrix = np.exp(-1j * 2.0 * np.pi * pilot_local[:, None] * delay[None, :] / float(grid.n_sc))
    pair_sums = sorted(
        float((delay[i] + delay[j]) % grid.n_sc) if modulo_pair_sums else float(delay[i] + delay[j])
        for i in range(n_tx) for j in range(i, n_tx)
    )
    unique_pair_sums = sorted(set(pair_sums))
    pair_gaps = [right - left for left, right in zip(unique_pair_sums, unique_pair_sums[1:])]
    if modulo_pair_sums and len(unique_pair_sums) > 1:
        pair_gaps.append(unique_pair_sums[0] + grid.n_sc - unique_pair_sums[-1])
    residues = [float(value % 96.0) for value in delay]
    unique_residues = sorted(set(residues))
    fold_gaps = [((unique_residues[(i + 1) % len(unique_residues)] - unique_residues[i]) % 96.0) for i in range(len(unique_residues))]
    return {
        "name": name,
        "delay_grid_coordinates": [float(value) for value in delay],
        "delay_ns": [float(value / (grid.n_sc * grid.scs_khz * 1e3) * 1e9) for value in delay],
        "delay_fft_samples": [float(value * grid.n_fft / grid.n_sc) for value in delay],
        "comb6_residues_mod_96": residues,
        "unordered_pair_sums": pair_sums,
        "unordered_pair_sums_unique": len(set(pair_sums)) == len(pair_sums),
        "minimum_pair_sum_gap": float(min(pair_gaps)) if pair_gaps else 0.0,
        "minimum_fold_gap": float(min(fold_gaps)) if fold_gaps else 0.0,
        "pilot_rank": int(np.linalg.matrix_rank(matrix, tol=1e-10)),
        "pilot_condition_number": float(np.linalg.cond(matrix)),
        "selection_rule": "plan-033 section 4 frozen analytic geometry",
    }


def validate_config(config: dict, *, plan: str = "033"):
    if config["run_kind"] not in {"smoke", "prescan", "formal"}:
        raise ValueError("run_kind must be smoke, prescan, or formal")
    if config["precoder_normalize"] is not PRECODER_NORMALIZE:
        raise ValueError("Plan-033 requires precoder_normalize: true.")
    if plan == "033":
        valid_scene = config["n_tx"] in {4, 8} and config["n_rx"] == N_RX and config["speed_kmh"] in {3.0, 60.0}
        valid_receivers = config["receiver_modes"] == ["estimated"]
    else:
        valid_scene = config["n_tx"] in {4, 8} and config["n_rx"] == 2 and config["speed_kmh"] == 60.0
        valid_receivers = config["receiver_modes"] == ["estimated", "ideal"]
    expected_channel = "cdl_d" if plan == "035-cdl" else "tdl_a"
    if config["channel_model"] != expected_channel:
        raise ValueError(f"Plan-{plan} requires channel_model: {expected_channel}.")
    if not valid_scene or not valid_receivers:
        raise ValueError(f"Plan-{plan} antenna, speed, or receiver_modes are outside the frozen scope.")
    expected_scenario = configured_scenario_id(config)
    if config["scenario_id"] != expected_scenario:
        raise ValueError(f"scenario_id must be {expected_scenario}")
    if config["feedback_period_slots"] != FEEDBACK_PERIOD_SLOTS or config["history_symbol_index"] != HISTORY_SYMBOL_INDEX:
        raise ValueError("Plan-033 aged CSI is frozen to slots10/index140.")
    if not config["snr_db"] or any(right <= left for left, right in zip(config["snr_db"], config["snr_db"][1:])):
        raise ValueError("snr_db must be non-empty and strictly increasing.")
    lower, upper = (-6.0, 20.0) if plan == "033" else (-2.0, 22.0)
    if min(config["snr_db"]) < lower or max(config["snr_db"]) > upper:
        raise ValueError(f"Plan-{plan} SNR must remain within [{lower:g},{upper:g}] dB.")
    if "target_total_trials_by_snr" in config and set(config["target_total_trials_by_snr"]) != set(config["snr_db"]):
        raise ValueError("Per-SNR target keys must exactly match snr_db.")
    targets = [_trial_target(config, value) for value in config["snr_db"]]
    if config["batch_size"] <= 0 or any(value <= 0 for value in targets):
        raise ValueError("Batch size and trial targets must be positive.")
    if config["run_kind"] != "smoke" and any(value % config["batch_size"] for value in targets):
        raise ValueError("Non-smoke trial targets must align to batch_size.")
    if config["run_kind"] == "smoke" and targets != [20] * len(targets):
        raise ValueError("Plan-033 smoke requires 20 trials per point.")
    if config["run_kind"] == "prescan" and targets != [400] * len(targets):
        raise ValueError("Plan-033 prescan requires 400 trials per point.")
    channel, grid, precoders, base = _build_scene(config)
    if grid.n_dmrs_re != 192 or list(np.unique(grid.pilot_symbol_indices)) != [2, 7] or grid.n_data_re != 5568:
        raise RuntimeError("Plan-033 resource grid changed.")
    if base.time.shape != (10, 10) or base.frequency.shape != (576, 576):
        raise RuntimeError("Unexpected covariance dimensions.")
    for name, precoder in precoders.items():
        if not np.allclose(np.sum(np.abs(precoder.C) ** 2, axis=1), 1.0, atol=1e-12):
            raise RuntimeError(f"Precoder power changed for {name}.")
    candidates = candidate_definitions(config["n_tx"])
    print(f"[validate] {config['scenario_id']}: curves=6 snr_points={len(config['snr_db'])} trials={targets[0]} batch={config['batch_size']}", flush=True)
    return candidates, channel, grid, precoders, base


def _effective_channels(realization_h: np.ndarray, grid, n_tx: int, precoders: dict, candidates: Sequence[dict] | None = None):
    ids = candidate_ids(n_tx)
    old = realization_h[:, :, :, 0, :]
    current = realization_h[:, :, :, 1:, :]
    if candidates is not None:
        effective = {}
        diagnostic_precoder = None
        for row in candidates:
            candidate_id = str(row["candidate_id"])
            if row["family"] == "CDD_MATCHED":
                weights = precoders[candidate_id].C
                effective[candidate_id] = equivalent_channel(current, weights)
                diagnostic_precoder = weights[None, :, :]
            elif row["family"] == "TRANSPARENT_PRG_DFT":
                order = np.tile(np.asarray(row["prg_vector_indices"], dtype=np.int64), (len(current), 1))
                weights = build_prg_dft_precoder_batch(grid, n_tx, PRG_SIZE_RB, order, normalize=True)
                effective[candidate_id] = equivalent_channel(current, weights)
                diagnostic_precoder = weights
            else:
                raise ValueError(f"Unsupported manifest candidate family: {row['family']}")
        if diagnostic_precoder is None:
            raise ValueError("Candidate manifest is empty.")
        return old, current, effective, diagnostic_precoder
    b0 = equivalent_channel(current, precoders["B0_QC"].C)
    sidon = equivalent_channel(current, precoders["S0_SIDON"].C)
    small = equivalent_channel(current, precoders["SMALL_CDD_QSTEP0P25"].C)
    order = np.tile(np.arange(8, dtype=np.int64) % n_tx, (len(current), 1))
    prg = build_prg_dft_precoder_batch(grid, n_tx, PRG_SIZE_RB, order, normalize=True)
    mrt = build_aged_mrt_prg_precoder(old, PRG_SIZE_RB, normalize=True)
    return old, current, {
        ids["b0"]: b0, ids["sidon"]: sidon,
        ids["prg"]: equivalent_channel(current, prg),
        ids["mrt"]: equivalent_channel(current, mrt),
        ids["small_transparent"]: small, ids["small_matched"]: small,
    }, mrt


def _build_filters(grid, base, n_tx: int, precoders: dict, noise_variance: float, candidates: Sequence[dict] | None = None):
    if candidates is not None:
        filters = {}
        transparent = TDLTimeFrequencyCovariance(
            base.time, TRANSPARENT_COVARIANCE_SCALE * base.frequency,
            "transparent_unit_power_physical",
        )
        for row in candidates:
            candidate_id = str(row["candidate_id"])
            if row["family"] == "CDD_MATCHED":
                c = precoders[candidate_id].C
                covariance = TDLTimeFrequencyCovariance(
                    base.time, base.frequency * (c @ c.conj().T),
                    f"matched_cdd:{candidate_id}",
                )
                filters[candidate_id] = build_time_frequency_rmmse_filter(
                    grid, covariance, noise_variance, diagonal_loading=1e-10
                )
            elif row["family"] == "TRANSPARENT_PRG_DFT":
                filters[candidate_id] = build_prg_time_frequency_rmmse_filter(
                    grid, transparent, 72, noise_variance, diagonal_loading=1e-10
                )
            else:
                raise ValueError(f"Unsupported manifest candidate family: {row['family']}")
        return filters
    ids = candidate_ids(n_tx)
    matched = {}
    for key, candidate_key in (("B0_QC", "b0"), ("S0_SIDON", "sidon"), ("SMALL_CDD_QSTEP0P25", "small_matched")):
        c = precoders[key].C
        covariance = TDLTimeFrequencyCovariance(base.time, base.frequency * (c @ c.conj().T), f"matched_cdd:{key}")
        matched[ids[candidate_key]] = build_time_frequency_rmmse_filter(grid, covariance, noise_variance, diagonal_loading=1e-10)
    transparent = TDLTimeFrequencyCovariance(
        base.time,
        TRANSPARENT_COVARIANCE_SCALE * base.frequency,
        "transparent_unit_power_physical",
    )
    matched[ids["small_transparent"]] = build_time_frequency_rmmse_filter(grid, transparent, noise_variance, diagonal_loading=1e-10)
    prg_filter = build_prg_time_frequency_rmmse_filter(grid, transparent, 72, noise_variance, diagonal_loading=1e-10)
    matched[ids["prg"]] = prg_filter
    matched[ids["mrt"]] = prg_filter
    return matched


def _write_manifests(config: dict, candidates: Sequence[dict], grid) -> None:
    output = Path(config["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    if config.get("_candidate_manifest") is None:
        audits = [_delay_audit(grid, config["n_tx"], name, values) for name, values in delay_sets(config["n_tx"]).items()]
    else:
        audits = [
            _delay_audit(grid, config["n_tx"], str(row["candidate_id"]), row["delay_grid_coordinates"], modulo_pair_sums=True)
            for row in candidates if row["family"] == "CDD_MATCHED"
        ]
        for row in candidates:
            if row["family"] == "TRANSPARENT_PRG_DFT":
                audits.append({
                    "name": str(row["candidate_id"]), "prg_size_rb": PRG_SIZE_RB,
                    "spatial_dft_size": config["n_tx"],
                    "prg_vector_indices": list(row["prg_vector_indices"]),
                    "phase_sign": "negative", "normalized": True,
                })
    candidate_path = output / "candidate_receiver_manifest.json"
    expanded_candidates = [
        {
            **dict(row),
            "n_tx": config["n_tx"],
            "n_rx": config["n_rx"],
            "speed_kmh": config["speed_kmh"],
            "receiver_modes": config["receiver_modes"],
            "precoder_normalize": config["precoder_normalize"],
        }
        for row in candidates
    ]
    candidate_path.write_text(json.dumps({"schema": config["schema"], "scenario_id": config["scenario_id"], "candidates": expanded_candidates}, indent=2) + "\n", encoding="utf-8")
    audit_path = output / "delay_audit.json"
    audit_path.write_text(json.dumps(audits, indent=2) + "\n", encoding="utf-8")
    actual_age_ms = config["history_symbol_index"] * grid.ofdm_symbol_duration_s * 1e3
    resolved = {
        "schema": config["schema"], "plan": config.get("plan_path", "research/plan-033-PDSCH-4Rx-MIMO移动性.md"), "run_kind": config["run_kind"],
        "config_path": _repo_relative(Path(config["config_path"])), "config_sha256": _sha256(Path(config["config_path"])),
        "scenario_id": config["scenario_id"], "candidate_manifest_sha256": _sha256(candidate_path), "delay_audit_sha256": _sha256(audit_path),
        "system": {
            "profile": "CDL-D" if config["channel_model"] == "cdl_d" else "TDL-A",
            "channel_backend": "sionna_cdl" if config["channel_model"] == "cdl_d" else "sionna_tdl",
            "delay_spread_ns": 100.0, "carrier_frequency_hz": 3.5e9,
            "ue_speed_kmh": config["speed_kmh"], "n_tx": config["n_tx"], "n_rx": config["n_rx"],
            "num_sinusoids": 20 if config["channel_model"] == "tdl_a" else None,
            "cdl_array": ({"tx": [1, config["n_tx"]], "rx": [1, config["n_rx"]], "polarization": "single-V", "element_spacing_wavelengths": 0.5} if config["channel_model"] == "cdl_d" else None),
            "spatial_covariance_used_by_ce": False if config["channel_model"] == "cdl_d" else None,
            "dmrs_symbol_indices": [2, 7], "pilot_re": grid.n_dmrs_re, "data_re": grid.n_data_re, "ici": False,
        },
        "transmitter": {"precoder_normalize": config["precoder_normalize"], "precoder_normalization": "unit squared norm on every RE", "total_power": 1.0, "transparent_covariance_scale": TRANSPARENT_COVARIANCE_SCALE},
        "receiver": {"modes": config["receiver_modes"], "estimated_type": f"two-DMRS 2D time-frequency RMMSE plus {config['n_rx']}Rx MRC", "ideal_type": f"true data-RE channel plus {config['n_rx']}Rx MRC", "noise_variance": "1/SNR_linear", "ce_error_aware_llr": False},
        "aged_csi": {"nominal_feedback_period_slots": FEEDBACK_PERIOD_SLOTS, "history_symbol_index": HISTORY_SYMBOL_INDEX, "current_symbol_sample_indices": list(range(HISTORY_SYMBOL_INDEX, HISTORY_SYMBOL_INDEX + 10)), "actual_age_ms": actual_age_ms, "old_current_same_realization": True},
        "mobility": {"speed_mps": config["speed_kmh"] / 3.6, "maximum_doppler_hz": (config["speed_kmh"] / 3.6) * 3.5e9 / 299792458.0},
        "snr_db": config["snr_db"], "target_total_trials_by_snr": {f"{snr:g}": _trial_target(config, snr) for snr in config["snr_db"]}, "batch_size": config["batch_size"], "seed": config["seed"],
    }
    path = output / "resolved_run.json"
    path.write_text(json.dumps(resolved, indent=2) + "\n", encoding="utf-8")
    (output / "resolved_run.sha256").write_text(f"{_sha256(path)}  {path.name}\n", encoding="utf-8")
    source_config = Path(config["config_path"])
    (output / "original_config.yaml").write_text(
        source_config.read_text(encoding="utf-8"), encoding="utf-8"
    )
    (output / "seed_pairing_audit.json").write_text(
        json.dumps(
            {
                "base_seed": config["seed"],
                "trial_key_fields": ["scenario_id", "snr_db", "trial_start", "trial_end"],
                "shared_across_candidates": ["channel", "payload", "data_awgn", "dmrs_awgn"],
                "seed_components": ["base_seed", "scenario_id", "snr_db", "absolute_batch_start", "stream"],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
            capture_output=True, text=True, encoding="utf-8",
        ).stdout.strip()
        worktree = subprocess.run(
            ["git", "status", "--short"], cwd=ROOT, check=True,
            capture_output=True, text=True, encoding="utf-8",
        ).stdout.splitlines()
    except (OSError, subprocess.CalledProcessError):
        commit, worktree = "unavailable", ["unavailable"]
    (output / "code_state.json").write_text(
        json.dumps({"git_commit": commit, "git_status_short": worktree}, indent=2) + "\n",
        encoding="utf-8",
    )


def _write_runtime_diagnostics(
    output: Path,
    elapsed_s: float,
    process_peak_rss_bytes: int,
) -> None:
    import tensorflow as tf

    gpus = [device.name for device in tf.config.list_logical_devices("GPU")]
    memory = None
    if gpus:
        try:
            memory = {
                key: int(value)
                for key, value in tf.config.experimental.get_memory_info("GPU:0").items()
            }
        except (ValueError, RuntimeError):
            memory = None
    (output / "runtime_diagnostics.json").write_text(
        json.dumps(
            {
                "schema": "plan033-runtime-diagnostics-v2",
                "tensorflow_version": tf.__version__,
                "execution_device": "GPU" if gpus else "CPU",
                "logical_gpus": gpus,
                "logical_cpus": [
                    device.name for device in tf.config.list_logical_devices("CPU")
                ],
                "gpu_memory_bytes": memory,
                "process_peak_rss_bytes_sampled": int(process_peak_rss_bytes),
                "memory_measurement_scope": "sampled_during_latest_trial_interval",
                "elapsed_seconds_latest_interval": float(elapsed_s),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def _audit_smoke_output(output: Path, candidates: Sequence[dict], snr_db: Sequence[float]) -> None:
    rows = read_csv_rows(output / "final" / "estimated_csi_bler_points.csv")
    expected_ids = {str(row["candidate_id"]) for row in candidates}
    if len(rows) != len(expected_ids) * len(snr_db):
        raise RuntimeError("Smoke output does not contain six curves at both SNR points.")
    if {str(row["candidate_id"]) for row in rows} != expected_ids:
        raise RuntimeError("Smoke candidate set changed.")
    intervals = read_csv_rows(output / "intervals.csv")
    checked_arrays = 0
    for row in intervals:
        keys = ["error_flags"] + (["ce_nmse_trial_path"] if row["receiver"] == "estimated" else [])
        for key in keys:
            values = np.load(_resolve(row[key]))
            if not np.all(np.isfinite(values)):
                raise RuntimeError(f"Non-finite smoke array: {row[key]}")
            checked_arrays += 1
    runtime = json.loads((output / "runtime_diagnostics.json").read_text(encoding="utf-8"))
    if "execution_device" not in runtime:
        runtime.update(
            {
                "schema": "plan033-runtime-diagnostics-v2-upgraded",
                "execution_device": "GPU" if runtime.get("logical_gpus") else "CPU",
                "process_peak_rss_bytes_sampled": None,
                "process_rss_bytes_at_resume_audit": int(psutil.Process().memory_info().rss),
                "memory_measurement_scope": "legacy_completed_interval; peak RSS unavailable",
            }
        )
        (output / "runtime_diagnostics.json").write_text(
            json.dumps(runtime, indent=2) + "\n", encoding="utf-8"
        )
    peak_rss = runtime.get("process_peak_rss_bytes_sampled")
    payload = {
        "status": "passed" if peak_rss is not None else "passed_with_legacy_runtime_memory_unavailable",
        "candidate_count": len(expected_ids),
        "receiver_modes": sorted({str(row["receiver"]) for row in intervals}),
        "snr_point_count": len(snr_db),
        "checked_array_count": checked_arrays,
        "power_mrc_diagnostic_count": len(list((output / "power_mrc_diagnostics").glob("*.json"))),
        "old_current_same_realization": True,
        "logical_gpus": runtime["logical_gpus"],
        "gpu_memory_bytes": runtime["gpu_memory_bytes"],
        "execution_device": runtime["execution_device"],
        "process_peak_rss_bytes_sampled": peak_rss,
        "memory_measurement_scope": runtime.get("memory_measurement_scope"),
    }
    (output / "smoke_audit.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )


def _merge_intervals(output: Path, candidates: Sequence[dict], expected_scenario: str) -> list[dict]:
    path = output / "intervals.csv"
    rows = [dict(row) for row in read_csv_rows(path)] if path.exists() else []
    grouped = defaultdict(list)
    for row in rows:
        if str(row["scenario_id"]) != expected_scenario:
            raise RuntimeError("Cross-scenario interval contamination detected.")
        grouped[(str(row["candidate_id"]), str(row["receiver"]), float(row["snr_db"]))].append(row)
    by_id = {str(row["candidate_id"]): row for row in candidates}
    merged = []
    for (candidate_id, receiver, snr_db), values in sorted(grouped.items()):
        values.sort(key=lambda row: int(row["trial_start"]))
        expected = 1
        for row in values:
            if int(row["trial_start"]) != expected:
                raise RuntimeError(f"Trial gap/overlap for {candidate_id} at {snr_db:g} dB.")
            expected = int(row["trial_end"]) + 1
        trials = sum(int(row["trials"]) for row in values)
        errors = sum(int(row["tb_errors"]) for row in values)
        low, high = wilson(errors, trials)
        result = {"scenario_id": expected_scenario, "n_tx": int(values[0]["n_tx"]), "n_rx": int(values[0]["n_rx"]), "speed_kmh": float(values[0]["speed_kmh"]), "candidate_id": candidate_id, "family": by_id[candidate_id]["family"], "receiver": receiver, "snr_db": snr_db, "trials": trials, "tb_errors": errors, "bler": errors / trials, "bler_wilson95_lo": low, "bler_wilson95_hi": high}
        if receiver == "estimated":
            nmse_sum = sum(float(row["ce_nmse_sum"]) for row in values)
            nmse_sumsq = sum(float(row["ce_nmse_sumsq"]) for row in values)
            mean = nmse_sum / trials
            variance = max((nmse_sumsq - trials * mean * mean) / max(trials - 1, 1), 0.0)
            se = math.sqrt(variance / trials)
            result.update({"ce_nmse_mean": mean, "ce_nmse_mean_db": 10.0 * math.log10(mean), "ce_nmse_standard_error": se, "ce_nmse_ci95_lo": max(mean - 1.96 * se, 0.0), "ce_nmse_ci95_hi": mean + 1.96 * se})
        merged.append(result)
    for receiver in sorted({str(row["receiver"]) for row in merged}):
        _write_rows_atomic([row for row in merged if row["receiver"] == receiver], output / "final" / f"{receiver}_csi_bler_points.csv")
    return merged


def run(config: dict, candidates: Sequence[dict], channel, grid, precoders, base) -> None:
    output = Path(config["output_dir"])
    _write_manifests(config, candidates, grid)
    by_id = {str(row["candidate_id"]): row for row in candidates}
    ids = list(by_id)
    intervals_path = output / "intervals.csv"
    stored = [dict(row) for row in read_csv_rows(intervals_path)] if intervals_path.exists() else []
    mcs = get_mcs("nr_256qam", 8, None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(tb.cb_k_values, tb.cb_e_values, num_iter=8, llr_clip=50.0)
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarrier_indices)
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    requested_times = [0, *range(HISTORY_SYMBOL_INDEX, HISTORY_SYMBOL_INDEX + 10)]
    for snr_db in config["snr_db"]:
        target = _trial_target(config, snr_db)
        counts = {(candidate_id, receiver): sum(int(row["trials"]) for row in stored if str(row["candidate_id"]) == candidate_id and str(row["receiver"]) == receiver and math.isclose(float(row["snr_db"]), snr_db, abs_tol=1e-12)) for candidate_id in ids for receiver in config["receiver_modes"]}
        if len(set(counts.values())) != 1:
            raise RuntimeError(f"Paired trial counts differ at {snr_db:g} dB.")
        current_total = next(iter(counts.values()))
        if current_total >= target:
            print(f"[resume] {config['scenario_id']} snr={snr_db:g} trials={current_total}", flush=True)
            continue
        trial_start = current_total + 1
        trial_end = _next_interval_end(current_total, target, config["run_kind"])
        trial_count = trial_end - trial_start + 1
        if config["run_kind"] != "smoke" and trial_count % config["batch_size"]:
            raise RuntimeError("Pending trial interval must align to batch_size.")
        noise_variance = awgn_variance(snr_db)
        filters = _build_filters(grid, base, config["n_tx"], precoders, noise_variance, candidates if config.get("_candidate_manifest") is not None else None)
        diagnostic = {candidate_id: {"condition_number": float(filters[candidate_id].condition_number), "numerical_jitter": float(filters[candidate_id].numerical_jitter), "covariance_type": filters[candidate_id].covariance_type, "subfilters": list(filters[candidate_id].subfilter_diagnostics)} for candidate_id in ids}
        diag_path = output / "filter_diagnostics" / f"{_snr_key(snr_db)}.json"
        diag_path.parent.mkdir(parents=True, exist_ok=True)
        diag_path.write_text(json.dumps(diagnostic, indent=2) + "\n", encoding="utf-8")
        flags = _make_error_flags(ids, config["receiver_modes"], trial_count)
        ce_nmse = {candidate_id: np.zeros(trial_count, dtype=np.float64) for candidate_id in ids}
        mrc_denominator_bounds = {
            candidate_id: {"minimum": float("inf"), "maximum": 0.0}
            for candidate_id in ids
        }
        mrt_power_bounds = {"minimum": float("inf"), "maximum": 0.0}
        correlation = {"old_current_cross": 0j, "old_power": 0.0, "current_power": 0.0, "edge_cross": 0j, "first_power": 0.0, "last_power": 0.0}
        process = psutil.Process()
        process_peak_rss_bytes = int(process.memory_info().rss)
        started = time.time()
        for absolute_start in range(trial_start, trial_end + 1, config["batch_size"]):
            batch = min(config["batch_size"], trial_end - absolute_start + 1)
            offset = absolute_start - trial_start
            payload_rng = np.random.default_rng(stable_seed(config["seed"], config["scenario_id"], snr_db, absolute_start, "payload"))
            payload = [payload_rng.integers(0, 2, size=int(k), dtype=np.int8) for k in tb.cb_k_values]
            symbols = qam_modulate(np.concatenate(adapter.encode(payload)), int(mcs.qm))
            realization = generate_sionna_channel_active(grid, channel, n_tx=config["n_tx"], n_rx=config["n_rx"], batch_size=batch, seed=stable_seed(config["seed"], config["scenario_id"], snr_db, absolute_start, "channel"), time_sample_indices=requested_times)
            old, current, effective, mrt = _effective_channels(realization.H, grid, config["n_tx"], precoders, candidates if config.get("_candidate_manifest") is not None else None)
            if not np.all(np.isfinite(current)) or any(not np.all(np.isfinite(value)) for value in effective.values()):
                raise RuntimeError("Non-finite channel sample detected.")
            if not np.allclose(np.sum(np.abs(mrt) ** 2, axis=2), 1.0, atol=1e-10):
                raise RuntimeError("Aged MRT power normalization changed.")
            mrt_power = np.sum(np.abs(mrt) ** 2, axis=2)
            mrt_power_bounds["minimum"] = min(mrt_power_bounds["minimum"], float(np.min(mrt_power)))
            mrt_power_bounds["maximum"] = max(mrt_power_bounds["maximum"], float(np.max(mrt_power)))
            correlation["old_current_cross"] += np.sum(old.conj() * current[:, :, :, 0, :])
            correlation["old_power"] += float(np.sum(np.abs(old) ** 2))
            correlation["current_power"] += float(np.sum(np.abs(current[:, :, :, 0, :]) ** 2))
            correlation["edge_cross"] += np.sum(current[:, :, :, 0, :].conj() * current[:, :, :, -1, :])
            correlation["first_power"] += float(np.sum(np.abs(current[:, :, :, 0, :]) ** 2))
            correlation["last_power"] += float(np.sum(np.abs(current[:, :, :, -1, :]) ** 2))
            noise_rng = np.random.default_rng(stable_seed(config["seed"], config["scenario_id"], snr_db, absolute_start, "noise"))
            ls_noise = math.sqrt(noise_variance / 2.0) * (noise_rng.normal(size=(batch, config["n_rx"], grid.n_dmrs_re)) + 1j * noise_rng.normal(size=(batch, config["n_rx"], grid.n_dmrs_re)))
            data_noise = math.sqrt(noise_variance / 2.0) * (noise_rng.normal(size=(batch, config["n_rx"], grid.n_data_re)) + 1j * noise_rng.normal(size=(batch, config["n_rx"], grid.n_data_re)))
            decode_batches = []
            for candidate_id in ids:
                true = effective[candidate_id]
                true_pilot = true[:, :, grid.pilot_symbol_indices, pilot_local]
                true_data = true[:, :, grid.data_symbol_indices, data_local]
                estimate_data = filters[candidate_id].estimate_data(true_pilot + ls_noise)
                _, _, ratios = _ce_trial_statistics(estimate_data, true_data)
                if not np.all(np.isfinite(estimate_data)) or not np.all(np.isfinite(ratios)):
                    raise RuntimeError(f"Non-finite CE result for {candidate_id}.")
                ce_nmse[candidate_id][offset:offset + batch] = ratios
                received = true_data * symbols[None, None, :] + data_noise
                denominator = np.sum(np.abs(estimate_data) ** 2, axis=1)
                mrc_denominator_bounds[candidate_id]["minimum"] = min(
                    mrc_denominator_bounds[candidate_id]["minimum"], float(np.min(denominator))
                )
                mrc_denominator_bounds[candidate_id]["maximum"] = max(
                    mrc_denominator_bounds[candidate_id]["maximum"], float(np.max(denominator))
                )
                equalized, effective_noise = _estimated_csi_mrc_equalize(estimate_data, received, noise_variance)
                decode_batches.append((candidate_id, "estimated", equalized, effective_noise))
                if "ideal" in config["receiver_modes"]:
                    ideal_equalized, ideal_noise = _estimated_csi_mrc_equalize(true_data, received, noise_variance)
                    decode_batches.append((candidate_id, "ideal", ideal_equalized, ideal_noise))
            combined_flags = _decode_flags(adapter, np.concatenate([item[2] for item in decode_batches]), np.concatenate([item[3] for item in decode_batches]), int(mcs.qm), payload)
            cursor = 0
            for candidate_id, receiver, _, _ in decode_batches:
                flags[(candidate_id, receiver)][offset:offset + batch] = combined_flags[cursor:cursor + batch]
                cursor += batch
            completed = absolute_start + batch - 1
            process_peak_rss_bytes = max(
                process_peak_rss_bytes,
                int(process.memory_info().rss),
            )
            print(f"[{config['run_kind']} {config['scenario_id']}] snr={snr_db:g} absolute_trials={completed}/{trial_end} elapsed={time.time()-started:.1f}s", flush=True)
        interval_rows = []
        for candidate_id in ids:
          for receiver in config["receiver_modes"]:
            receiver_flags = flags[(candidate_id, receiver)]
            flag_path = _array_path(output, "error_flags", snr_db, candidate_id, trial_start, trial_end, receiver)
            nmse_path = _array_path(output, "ce_nmse_trials", snr_db, candidate_id, trial_start, trial_end)
            flag_path.parent.mkdir(parents=True, exist_ok=True)
            np.save(flag_path, receiver_flags)
            errors = int(np.sum(receiver_flags)); low, high = wilson(errors, trial_count)
            row = {**_interval_identity(config, candidate_id, snr_db, trial_start, trial_end, receiver), "family": by_id[candidate_id]["family"], "snr_db": snr_db, "trial_start": trial_start, "trial_end": trial_end, "trials": trial_count, "tb_errors": errors, "bler": errors / trial_count, "bler_wilson95_lo": low, "bler_wilson95_hi": high, "error_flags": _repo_relative(flag_path), "seed": config["seed"]}
            if receiver == "estimated":
                nmse_path.parent.mkdir(parents=True, exist_ok=True)
                np.save(nmse_path, ce_nmse[candidate_id])
                row.update({"ce_nmse_sum": float(np.sum(ce_nmse[candidate_id])), "ce_nmse_sumsq": float(np.sum(ce_nmse[candidate_id] ** 2)), "ce_nmse_trial_path": _repo_relative(nmse_path)})
            interval_rows.append(row)
        stored.extend(interval_rows)
        stored.sort(key=lambda row: (str(row["candidate_id"]), str(row["receiver"]), float(row["snr_db"]), int(row["trial_start"])))
        _write_rows_atomic(stored, intervals_path)
        old_corr = correlation["old_current_cross"] / math.sqrt(correlation["old_power"] * correlation["current_power"])
        edge_corr = correlation["edge_cross"] / math.sqrt(correlation["first_power"] * correlation["last_power"])
        with (output / "channel_correlation_intervals.jsonl").open("a", encoding="utf-8") as handle:
            handle.write(json.dumps({"scenario_id": config["scenario_id"], "snr_db": snr_db, "trial_start": trial_start, "trial_end": trial_end, "old_to_current_first_complex_correlation": [float(old_corr.real), float(old_corr.imag)], "current_symbol0_to9_complex_correlation": [float(edge_corr.real), float(edge_corr.imag)]}) + "\n")
        power_mrc_path = output / "power_mrc_diagnostics" / f"{_snr_key(snr_db)}_t{trial_start:06d}_{trial_end:06d}.json"
        power_mrc_path.parent.mkdir(parents=True, exist_ok=True)
        power_mrc_path.write_text(
            json.dumps(
                {
                    "scenario_id": config["scenario_id"],
                    "snr_db": snr_db,
                    "noise_variance": noise_variance,
                    "expected_precoder_power": 1.0,
                    "mrt_power_bounds": mrt_power_bounds,
                    "mrc_denominator_bounds": mrc_denominator_bounds,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        with (output / "run_log.jsonl").open("a", encoding="utf-8") as handle:
            handle.write(
                json.dumps(
                    {
                        "event": "interval_complete",
                        "scenario_id": config["scenario_id"],
                        "snr_db": snr_db,
                        "trial_start": trial_start,
                        "trial_end": trial_end,
                        "elapsed_seconds": time.time() - started,
                    }
                )
                + "\n"
            )
        _write_runtime_diagnostics(
            output,
            time.time() - started,
            process_peak_rss_bytes,
        )
        _merge_intervals(output, candidates, config["scenario_id"])
    if config["run_kind"] == "smoke":
        _audit_smoke_output(output, candidates, config["snr_db"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--stage", choices=("validate", "run", "merge"), default="run")
    args = parser.parse_args()
    config = load_config(args.config.resolve())
    candidates, channel, grid, precoders, base = validate_config(config)
    if args.stage == "run":
        run(config, candidates, channel, grid, precoders, base)
    elif args.stage == "merge":
        _merge_intervals(Path(config["output_dir"]), candidates, config["scenario_id"])


if __name__ == "__main__":
    main()
