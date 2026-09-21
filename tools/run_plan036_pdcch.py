from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import shutil
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.sim.pdcch import wilson_interval  # noqa: E402
from cdd_lls.sim.pdcch_cdd import (  # noqa: E402
    PDCCHCDDBLERSimulator,
    validate_pdcch_cdd_config,
)
from cdd_lls.sim.stats import save_csv, save_json  # noqa: E402


SCHEMA = "plan036-pdcch-paired-v1"
LIGHT_SPEED_MPS = 299_792_458.0
ESTIMATED_IDS = (
    "SIDON",
    "QC_DELAY_NT",
    "QC_DELAY_TRANSPARENT",
    "DFT4_CYCLING_BASELINE",
)
IDEAL_IDS = ("SIDON", "QC_DELAY_NT", "DFT4_CYCLING_BASELINE")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _tag(value: float) -> str:
    return str(float(value)).replace("-", "m").replace(".", "p")


def _read_rows(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def load_config(path: str | Path) -> dict[str, Any]:
    source = Path(path).resolve()
    config = yaml.safe_load(source.read_text(encoding="utf-8")) or {}
    config["config_path"] = str(source)
    config["output_dir"] = str((ROOT / str(config["output_dir"])).resolve())
    validate_config(config)
    return config


def _candidate_specs(config: dict[str, Any]) -> list[dict[str, Any]]:
    al = int(config["aggregation_level"])
    if al == 1:
        sidon = [0, 1, 3, 7]
        qc = [0.0, 0.24192, 0.48384, 0.72576]
        qc_family = "QC-delay-112"
        cycling = [0]
    else:
        sidon = [0, 11, 19, 64]
        qc = [0.0, 0.41472, 0.82944, 1.24416]
        qc_family = "QC-delay-96"
        cycling = [0, 1]
    return [
        {
            "candidate_id": "SIDON",
            "family": "SIDON",
            "scheme": "cdd",
            "label": f"Sidon AL{al} non-transparent",
            "delay_grid_coordinates": sidon,
            "strict_sidon": True,
            "receiver_covariance_mode": "matched_effective",
            "style": {"color": "#1f77b4", "marker": "o"},
        },
        {
            "candidate_id": "QC_DELAY_NT",
            "family": qc_family,
            "scheme": "cdd",
            "label": f"{qc_family} non-transparent",
            "delay_grid_coordinates": qc,
            "receiver_covariance_mode": "matched_effective",
            "style": {"color": "#ff7f0e", "marker": "s"},
        },
        {
            "candidate_id": "QC_DELAY_TRANSPARENT",
            "family": qc_family,
            "scheme": "cdd",
            "label": f"{qc_family} transparent",
            "delay_grid_coordinates": qc,
            "receiver_covariance_mode": "physical_fullband",
            "style": {"color": "#ff7f0e", "linestyle": "--", "marker": "^"},
        },
        {
            "candidate_id": "DFT4_CYCLING_BASELINE",
            "family": "DFTcodebook",
            "scheme": "reg_bundle_dft_cycling",
            "label": "DFT4 cycling transparent",
            "cycling_order": cycling,
            "receiver_covariance_mode": "physical_prg",
            "style": {"color": "#2ca02c", "marker": "D"},
        },
    ]


def _base_config(config: dict[str, Any], csi_mode: str) -> dict[str, Any]:
    al = int(config["aggregation_level"])
    candidates = _candidate_specs(config)
    if csi_mode == "ideal":
        candidates = [item for item in candidates if item["candidate_id"] in IDEAL_IDS]
    return {
        "schema": "pdcch-cdd-bler-v2",
        "output_dir": str(Path(config["output_dir"]) / csi_mode),
        "seed": int(config["seed"]),
        "random_stream_namespace": str(config["random_stream_namespace"]),
        "batch_size": int(config["batch_size"]),
        "antenna": {"n_tx": 4, "n_rx": int(config["n_rx"])},
        "resource": {
            "n_rb": 48,
            "duration_symbols": 1,
            "n_fft": 4096,
            "scs_khz": 30,
            "cyclic_prefix_length": 288,
            "cce_reg_mapping": "noninterleaved",
            "reg_bundle_size": 6,
            "interleaver_size": 2,
            "shift_index": 0,
            "aggregation_level": al,
            "first_cce": 0,
        },
        "pdcch": {
            "payload_bits": 40,
            "coded_bits": 108 * al,
            "crc_rnti": 65535,
            "scrambling_id": 208,
            "data_scrambling_rnti": 0,
            "polar_list_size": 8,
            "polar_decoder_type": "hybSCL",
            "polar_cpu_only": True,
            "coreset_start_symbol": 0,
            "slot_number": 0,
        },
        "channel": {
            "backend": "sionna_tdl",
            "model": "tdl",
            "tdl_profile": "C",
            "delay_spread_ns": 300.0,
            "carrier_frequency_hz": float(config["carrier_frequency_hz"]),
            "ue_speed_kmh": float(config["ue_speed_kmh"]),
            "num_sinusoids": 20,
            "normalize": True,
        },
        "receiver": {"channel_estimation": "ideal" if csi_mode == "ideal" else "frequency_lmmse", "llr_clip": 50.0},
        "simulation": {
            "snr_points_db": [float(value) for value in config["snr_points_db"]],
            "min_trials_per_snr": int(config["minimum_trials"]),
            "target_errors": int(config["target_errors"]),
            "max_trials_per_snr": int(config["maximum_trials"]),
            "save_trial_error_flags": True,
            "resume": True,
        },
        "candidates": candidates,
    }


def validate_config(config: dict[str, Any]) -> dict[str, Any]:
    if str(config.get("schema")) != SCHEMA:
        raise ValueError(f"schema must be {SCHEMA}.")
    al = int(config["aggregation_level"])
    if al not in (1, 2):
        raise ValueError("Plan-036 aggregation_level must be 1 or 2.")
    n_rx = int(config["n_rx"])
    if n_rx not in (2, 4):
        raise ValueError("This runner supports n_rx=2 or n_rx=4.")
    requested = float(config["requested_doppler_hz"])
    if requested not in (5.0, 1100.0):
        raise ValueError("Plan-036 requested_doppler_hz must be 5 or 1100 Hz.")
    carrier = float(config["carrier_frequency_hz"])
    if float(config["light_speed_mps"]) != LIGHT_SPEED_MPS:
        raise ValueError(f"light_speed_mps must be exactly {LIGHT_SPEED_MPS:.1f}.")
    speed = float(config["ue_speed_kmh"])
    expected_speed = requested * LIGHT_SPEED_MPS / carrier * 3.6
    if abs(speed - expected_speed) > 1e-9:
        raise ValueError(f"ue_speed_kmh must be the unrounded Doppler conversion {expected_speed:.12f}.")
    recovered = speed / 3.6 * carrier / LIGHT_SPEED_MPS
    snrs = [float(value) for value in config["snr_points_db"]]
    allowed_snr_grids = (
        [float(value) for value in range(-10, 6)],
        [float(value) for value in range(-10, 9)],
    )
    if snrs not in allowed_snr_grids:
        raise ValueError("The SNR grid must be exactly -10:1:5 dB or the approved -10:1:8 dB extension.")
    minimum = int(config["minimum_trials"])
    interval = int(config["check_interval_trials"])
    maximum = int(config["maximum_trials"])
    if (minimum, interval, maximum, int(config["target_errors"])) != (2000, 1000, 50000, 200):
        raise ValueError("Plan-036 adaptive budget must be 2000/1000/50000 trials and 200 errors.")
    if int(config["batch_size"]) != interval:
        raise ValueError("batch_size must equal the 1000-trial stopping interval.")
    estimated = _base_config(config, "estimated")
    ideal = _base_config(config, "ideal")
    estimated_receipt = validate_pdcch_cdd_config(estimated)
    ideal_receipt = validate_pdcch_cdd_config(ideal)
    return {
        "schema": SCHEMA,
        "aggregation_level": al,
        "n_rx": n_rx,
        "requested_doppler_hz": requested,
        "carrier_frequency_hz": carrier,
        "light_speed_mps": LIGHT_SPEED_MPS,
        "ue_speed_kmh": speed,
        "recovered_doppler_hz": recovered,
        "doppler_roundtrip_error_hz": recovered - requested,
        "estimated": estimated_receipt,
        "ideal": ideal_receipt,
    }


def _curve_keys() -> list[tuple[str, str]]:
    return [(value, "estimated") for value in ESTIMATED_IDS] + [
        (value, "ideal") for value in IDEAL_IDS
    ]


def _stop_reason(errors: int, trials: int, config: dict[str, Any]) -> str | None:
    if trials < int(config["minimum_trials"]):
        return None
    if errors >= int(config["target_errors"]):
        return "target_errors"
    if wilson_interval(errors, trials)[1] < 0.01:
        return "wilson95_upper_below_0p01"
    if trials >= int(config["maximum_trials"]):
        return "maximum_trials"
    return None


def _save_checkpoint(path: Path, trials: int, flags: dict[tuple[str, str], list[np.ndarray]], nmse: dict[str, list[np.ndarray]]) -> None:
    payload: dict[str, np.ndarray] = {"trials": np.asarray([trials], dtype=np.int64)}
    for candidate_id, mode in _curve_keys():
        payload[f"flags__{mode}__{candidate_id}"] = np.concatenate(flags[(candidate_id, mode)]).astype(bool)
    for candidate_id in ESTIMATED_IDS:
        payload[f"nmse__{candidate_id}"] = np.concatenate(nmse[candidate_id]).astype(np.float64)
    temporary = path.with_suffix(".tmp.npz")
    np.savez_compressed(temporary, **payload)
    temporary.replace(path)


def _load_checkpoint(path: Path) -> tuple[int, dict[tuple[str, str], list[np.ndarray]], dict[str, list[np.ndarray]]]:
    flags = {key: [] for key in _curve_keys()}
    nmse = {value: [] for value in ESTIMATED_IDS}
    if not path.exists():
        return 0, flags, nmse
    with np.load(path) as saved:
        trials = int(saved["trials"][0])
        for candidate_id, mode in _curve_keys():
            flags[(candidate_id, mode)].append(np.asarray(saved[f"flags__{mode}__{candidate_id}"], dtype=bool))
        for candidate_id in ESTIMATED_IDS:
            nmse[candidate_id].append(np.asarray(saved[f"nmse__{candidate_id}"], dtype=np.float64))
    return trials, flags, nmse


def _row(candidate, mode: str, snr_db: float, flags: np.ndarray, nmse: np.ndarray | None, reason: str, simulator: PDCCHCDDBLERSimulator) -> dict[str, Any]:
    errors = int(np.sum(flags))
    trials = int(len(flags))
    low, high = wilson_interval(errors, trials)
    return {
        "candidate_id": candidate.candidate_id,
        "family": candidate.family,
        "csi_mode": mode,
        "receiver_covariance_mode": candidate.receiver_covariance_mode,
        "snr_db": float(snr_db),
        "trials": trials,
        "errors": errors,
        "bler": errors / trials,
        "bler_wilson95_lo": low,
        "bler_wilson95_hi": high,
        "ce_nmse_db": "" if nmse is None else 10.0 * math.log10(max(float(np.mean(nmse)), 1e-30)),
        "stopping_reason": reason,
        "aggregation_level": int(simulator.resource_config.aggregation_level),
        "occupied_rb": int(simulator.k_active // 12),
        "k_active": simulator.k_active,
        "data_re": int(len(simulator.data_local)),
        "dmrs_re": int(len(simulator.pilot_local)),
        "coded_bits": int(simulator.codec.coded_bits),
        "pilot_rank": int(candidate.pilot_rank),
        "pilot_condition_number": float(candidate.pilot_condition_number),
        "ce_floor_nmse_db": 10.0 * math.log10(max(candidate.ce_floor_nmse, 1e-30)),
    }


def run(config: dict[str, Any], *, smoke: bool = False) -> None:
    working = copy.deepcopy(config)
    output = Path(working["output_dir"])
    if smoke:
        output = output / "smoke"
        working["snr_points_db"] = [-10.0, 5.0]
        working["minimum_trials"] = 20
        working["check_interval_trials"] = 20
        working["maximum_trials"] = 20
        working["target_errors"] = 1_000_000
        working["batch_size"] = 20
        working["random_stream_namespace"] += "_smoke"
    output.mkdir(parents=True, exist_ok=True)
    estimated_config = _base_config(working, "estimated")
    ideal_config = _base_config(working, "ideal")
    estimated = PDCCHCDDBLERSimulator(estimated_config)
    ideal = PDCCHCDDBLERSimulator(ideal_config)
    source = Path(config["config_path"])
    shutil.copyfile(source, output / source.name)
    (output / f"{source.stem}.sha256").write_text(_sha256(source) + "\n", encoding="utf-8")
    with (output / "resolved_config.yaml").open("w", encoding="utf-8") as handle:
        yaml.safe_dump(working, handle, allow_unicode=True, sort_keys=False)
    rows_path = output / "bler_points.csv"
    rows = _read_rows(rows_path)
    completed = {float(row["snr_db"]) for row in rows}
    flags_dir = output / "trial_data"
    flags_dir.mkdir(exist_ok=True)
    adaptive_path = output / "adaptive_status.json"
    adaptive_history = []
    if adaptive_path.exists():
        adaptive_history = list(json.loads(adaptive_path.read_text(encoding="utf-8")).get("history", []))
    started = time.time()
    for snr_db in working["snr_points_db"]:
        snr_db = float(snr_db)
        if snr_db in completed:
            continue
        checkpoint = flags_dir / f"snr_{_tag(snr_db)}_partial.npz"
        trials, flags, nmse = _load_checkpoint(checkpoint)
        while trials < int(working["maximum_trials"]):
            current = min(int(working["check_interval_trials"]), int(working["maximum_trials"]) - trials)
            estimated_result = estimated.run_batch(snr_db, trials + 1, current)
            ideal_result = ideal.run_batch(snr_db, trials + 1, current)
            for candidate_id in ESTIMATED_IDS:
                flags[(candidate_id, "estimated")].append(estimated_result[candidate_id]["error_flags"])
                nmse[candidate_id].append(estimated_result[candidate_id]["ce_nmse"])
            for candidate_id in IDEAL_IDS:
                flags[(candidate_id, "ideal")].append(ideal_result[candidate_id]["error_flags"])
            trials += current
            _save_checkpoint(checkpoint, trials, flags, nmse)
            reasons = {
                key: _stop_reason(int(np.sum(np.concatenate(value))), trials, working)
                for key, value in flags.items()
            }
            decision = {
                "snr_db": snr_db,
                "trials": trials,
                "curves": {
                    f"{mode}:{candidate_id}": {
                        "errors": int(np.sum(np.concatenate(flags[(candidate_id, mode)]))),
                        "wilson95": list(
                            wilson_interval(
                                int(np.sum(np.concatenate(flags[(candidate_id, mode)]))), trials
                            )
                        ),
                        "stop_reason": reasons[(candidate_id, mode)],
                    }
                    for candidate_id, mode in _curve_keys()
                },
            }
            adaptive_history.append(decision)
            save_json(
                {
                    "schema": "plan036-adaptive-status-v1",
                    "check_interval_trials": int(working["check_interval_trials"]),
                    "history": adaptive_history,
                    "latest": decision,
                },
                adaptive_path,
            )
            print(f"AL{working['aggregation_level']} fd={working['requested_doppler_hz']:g}Hz SNR={snr_db:g}dB trials={trials} " + " ".join(f"{mode}:{candidate}={int(np.sum(np.concatenate(flags[(candidate, mode)])))}" for candidate, mode in _curve_keys()), flush=True)
            if all(value is not None for value in reasons.values()):
                break
        estimated_by_id = {item.candidate_id: item for item in estimated.candidates}
        ideal_by_id = {item.candidate_id: item for item in ideal.candidates}
        reasons = {key: _stop_reason(int(np.sum(np.concatenate(value))), trials, working) or "maximum_trials" for key, value in flags.items()}
        for candidate_id, mode in _curve_keys():
            all_flags = np.concatenate(flags[(candidate_id, mode)]).astype(bool)
            all_nmse = np.concatenate(nmse[candidate_id]) if mode == "estimated" else None
            candidate = estimated_by_id[candidate_id] if mode == "estimated" else ideal_by_id[candidate_id]
            rows.append(_row(candidate, mode, snr_db, all_flags, all_nmse, reasons[(candidate_id, mode)], estimated if mode == "estimated" else ideal))
            np.save(flags_dir / f"{mode}_{candidate_id}_snr_{_tag(snr_db)}_error_flags.npy", all_flags)
            if all_nmse is not None:
                np.save(flags_dir / f"estimated_{candidate_id}_snr_{_tag(snr_db)}_ce_nmse.npy", all_nmse)
        save_csv(rows, rows_path)
        checkpoint.unlink(missing_ok=True)
    qc_nt = next(item for item in estimated.candidates if item.candidate_id == "QC_DELAY_NT")
    qc_tr = next(item for item in estimated.candidates if item.candidate_id == "QC_DELAY_TRANSPARENT")
    metadata = {
        "schema": "plan036-pdcch-paired-result-v1",
        "validation": validate_config(config),
        "elapsed_seconds_this_invocation": time.time() - started,
        "pairing": "estimated and ideal calls use identical seed namespace, SNR, absolute_start, payload, channel, and data-noise keys",
        "qc_transmit_matrices_elementwise_equal": bool(np.array_equal(qc_nt.precoder.C, qc_tr.precoder.C)),
        "ideal_ce_nmse_elementwise_zero_by_construction": True,
        "rows": rows,
    }
    save_json(metadata, output / "run_metadata.json")


def analyze(config: dict[str, Any], *, smoke: bool = False) -> None:
    import matplotlib.pyplot as plt

    output = Path(config["output_dir"]) / ("smoke" if smoke else "")
    rows = _read_rows(output / "bler_points.csv")
    if not rows:
        raise RuntimeError(f"No BLER rows found in {output}.")
    figures = output / "figures"
    figures.mkdir(exist_ok=True)
    published = ROOT / "docs" / "figures" / "result-036"
    if not smoke:
        published.mkdir(parents=True, exist_ok=True)
    styles = {
        "SIDON": {"color": "#1f77b4", "linestyle": "-", "marker": "o"},
        "QC_DELAY_NT": {"color": "#ff7f0e", "linestyle": "-", "marker": "s"},
        "QC_DELAY_TRANSPARENT": {"color": "#ff7f0e", "linestyle": "--", "marker": "^"},
        "DFT4_CYCLING_BASELINE": {"color": "#2ca02c", "linestyle": "-.", "marker": "D"},
    }
    rx_suffix = "" if int(config["n_rx"]) == 4 else f"_nr{int(config['n_rx'])}"
    scene = f"al{config['aggregation_level']}_fd{int(config['requested_doppler_hz'])}{rx_suffix}"
    for mode, ids in (("estimated", ESTIMATED_IDS), ("ideal", IDEAL_IDS)):
        fig, axis = plt.subplots(figsize=(7.6, 5.1))
        for candidate_id in ids:
            points = sorted((row for row in rows if row["csi_mode"] == mode and row["candidate_id"] == candidate_id), key=lambda row: float(row["snr_db"]))
            snr = np.asarray([float(row["snr_db"]) for row in points])
            trials = np.asarray([int(row["trials"]) for row in points])
            bler = np.asarray([float(row["bler"]) for row in points])
            axis.semilogy(snr, np.where(bler > 0, bler, 0.5 / trials), label=candidate_id, linewidth=2.5, markersize=8, **styles[candidate_id])
        axis.axhline(0.1, color="black", linestyle=":", linewidth=1.8)
        axis.axhline(0.01, color="black", linestyle="--", linewidth=1.8)
        axis.set(xlabel="SNR (dB)", ylabel="DCI BLER", title=f"4T{config['n_rx']}R AL{config['aggregation_level']} {config['requested_doppler_hz']:g} Hz {mode} CSI")
        axis.xaxis.label.set_size(16)
        axis.yaxis.label.set_size(16)
        axis.title.set_size(17)
        axis.tick_params(labelsize=14)
        axis.grid(True, which="both", alpha=0.3)
        axis.legend(fontsize=12, loc="best")
        fig.tight_layout()
        destination = figures / f"{mode}_bler.png"
        fig.savefig(destination, dpi=180)
        if not smoke:
            shutil.copyfile(destination, published / f"{scene}_{mode}_bler.png")
        plt.close(fig)
    fig, axis = plt.subplots(figsize=(7.6, 5.1))
    for candidate_id in ESTIMATED_IDS:
        points = sorted((row for row in rows if row["csi_mode"] == "estimated" and row["candidate_id"] == candidate_id), key=lambda row: float(row["snr_db"]))
        axis.plot([float(row["snr_db"]) for row in points], [float(row["ce_nmse_db"]) for row in points], label=candidate_id, linewidth=2.5, markersize=8, **styles[candidate_id])
    axis.set(xlabel="SNR (dB)", ylabel="Data-RE CE NMSE (dB)", title=f"4T{config['n_rx']}R AL{config['aggregation_level']} {config['requested_doppler_hz']:g} Hz estimated CSI")
    axis.xaxis.label.set_size(16)
    axis.yaxis.label.set_size(16)
    axis.title.set_size(17)
    axis.tick_params(labelsize=14)
    axis.grid(True, alpha=0.3)
    axis.legend(fontsize=12, loc="best")
    fig.tight_layout()
    destination = figures / "estimated_ce_nmse.png"
    fig.savefig(destination, dpi=180)
    if not smoke:
        shutil.copyfile(destination, published / f"{scene}_estimated_ce_nmse.png")
    plt.close(fig)
    save_json({"schema": "plan036-style-map-v1", "styles": styles}, output / "figure_styles.json")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run plan-036 paired estimated/ideal PDCCH simulations.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--stage", choices=("validate", "smoke", "run", "analyze", "analyze-smoke"), required=True)
    parser.add_argument("--validation-output")
    args = parser.parse_args()
    config = load_config(args.config)
    if args.stage == "validate":
        receipt = validate_config(config)
        rendered = json.dumps(receipt, indent=2, ensure_ascii=False)
        if args.validation_output:
            destination = Path(args.validation_output)
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text(rendered + "\n", encoding="utf-8")
        print(rendered)
    elif args.stage == "smoke":
        run(config, smoke=True)
    elif args.stage == "run":
        run(config)
    elif args.stage == "analyze-smoke":
        analyze(config, smoke=True)
    else:
        analyze(config)


if __name__ == "__main__":
    main()
