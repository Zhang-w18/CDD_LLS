from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.sim.pdcch import (  # noqa: E402
    PDCCHBLERSimulator,
    _stable_seed,
    load_pdcch_bler_config,
    validate_pdcch_bler_config,
)
from cdd_lls.sim.pdcch_cdd import (  # noqa: E402
    PDCCHCDDBLERSimulator,
    SCHEMA as PDCCH_CDD_SCHEMA,
    load_pdcch_cdd_config,
    validate_pdcch_cdd_config,
)
from cdd_lls.phy.channel_tdl import generate_sionna_tdl_channel_active  # noqa: E402


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _archive_config(source: Path, output: Path) -> None:
    if output.parent.name == "bler":
        experiment_root = output.parent.parent
    elif output.parent.name == "smoke":
        experiment_root = output.parent.parent
    else:
        experiment_root = output.parent
    destination_dir = experiment_root / "configs"
    destination_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination_dir / source.name)
    (destination_dir / f"{source.stem}.sha256").write_text(
        _sha256(source) + "\n", encoding="utf-8"
    )


def _candidate_stream_seed(config: dict, candidate_id: str) -> int:
    return _stable_seed(
        int(config["seed"]),
        str(config.get("random_stream_namespace", "default")),
        int(config["resource"]["aggregation_level"]),
        str(candidate_id),
        "independent_candidate_stream",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Sionna PDCCH DCI BLER curves.")
    parser.add_argument("--config", required=True, help="PDCCH YAML configuration path.")
    parser.add_argument("--stage", choices=("validate", "smoke", "run"), default="run")
    parser.add_argument(
        "--candidate",
        help="For pdcch-cdd-bler-v2, run one candidate in an independent derived random stream.",
    )
    parser.add_argument(
        "--validation-output",
        help="Optional JSON path for the complete validation receipt.",
    )
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}
    if str(raw.get("schema")) == PDCCH_CDD_SCHEMA:
        config = load_pdcch_cdd_config(args.config)
        if args.candidate:
            selected = [
                candidate
                for candidate in config["candidates"]
                if str(candidate["candidate_id"]) == str(args.candidate)
            ]
            if len(selected) != 1:
                raise ValueError(f"Unknown or duplicate candidate_id {args.candidate!r}.")
            config = copy.deepcopy(config)
            config["candidates"] = selected
            config["seed"] = _candidate_stream_seed(config, str(args.candidate))
            config["output_dir"] = str(Path(config["output_dir"]) / str(args.candidate))
            config["sample_mode"] = "independent_candidate_stream"
        else:
            config["sample_mode"] = "shared_candidates"
        if args.stage == "smoke":
            config = copy.deepcopy(config)
            configured_output = Path(config["output_dir"])
            if configured_output.parent.name == "bler":
                experiment_root = configured_output.parent.parent
            else:
                experiment_root = configured_output.parent
            config["output_dir"] = str(experiment_root / "smoke" / "bler")
            config["random_stream_namespace"] = (
                str(config.get("random_stream_namespace", "default")) + "_smoke"
            )
            config["simulation"].update(
                {
                    "snr_points_db": [2.0],
                    "min_trials_per_snr": 20,
                    "target_errors": 1000000,
                    "max_trials_per_snr": 20,
                    "resume": False,
                    "save_batch_checkpoints": False,
                }
            )
            config["batch_size"] = 20
        summary = validate_pdcch_cdd_config(config)
        simulator = PDCCHCDDBLERSimulator(config)
        if args.stage == "validate":
            rng = np.random.default_rng(20260908)
            payload = rng.integers(
                0, 2, size=(2, int(simulator.codec.payload_bits)), dtype="int8"
            )
            transmitted = simulator.codec.encode(payload)
            decoded = simulator.codec.decode((2.0 * transmitted - 1.0) * 30.0)
            noiseless_ok = bool(
                np.all(decoded.crc_status)
                and np.array_equal(decoded.payload_bits, payload)
            )
            if not noiseless_ok:
                raise RuntimeError("Noiseless PDCCH encode/decode validation failed.")
            rx_channel_receipt = None
            if int(config["antenna"]["n_rx"]) > 1:
                n_rx = int(config["antenna"]["n_rx"])
                channel = generate_sionna_tdl_channel_active(
                    simulator.band_grid,
                    simulator.channel_config,
                    n_tx=int(config["antenna"]["n_tx"]),
                    n_rx=n_rx,
                    batch_size=16,
                    seed=_stable_seed(
                        int(config["seed"]),
                        str(config.get("random_stream_namespace", "default")),
                        int(config["resource"]["aggregation_level"]),
                        "two_rx_validation_channel",
                    ),
                ).H
                correlations = []
                branches_identical = False
                for left in range(n_rx):
                    for right in range(left + 1, n_rx):
                        first = channel[:, left, ...].reshape(-1)
                        second = channel[:, right, ...].reshape(-1)
                        denominator = np.sqrt(
                            np.vdot(first, first).real * np.vdot(second, second).real
                        )
                        correlation = np.vdot(first, second) / max(float(denominator), 1e-30)
                        correlations.append(
                            {
                                "left": left,
                                "right": right,
                                "real": float(correlation.real),
                                "imag": float(correlation.imag),
                                "magnitude": float(abs(correlation)),
                            }
                        )
                        branches_identical = branches_identical or bool(
                            np.array_equal(first, second)
                        )
                rx_channel_receipt = {
                    "shape": list(channel.shape),
                    "batch_size": 16,
                    "n_rx": n_rx,
                    "any_branch_pair_elementwise_identical": branches_identical,
                    "empirical_complex_cross_correlations": correlations,
                    "two_symbols_elementwise_identical": bool(
                        np.array_equal(channel[..., 0, :], channel[..., 1, :])
                    ),
                    "interpretation": "finite-sample smoke receipt, not proof of independence",
                }
            ideal_link_receipt = None
            if str(config["receiver"]["channel_estimation"]).lower() == "ideal":
                ideal_result = simulator.run_batch(
                    snr_db=float("inf"), absolute_start=1, batch_size=2
                )
                ideal_ce_exact = all(
                    np.array_equal(value["ce_nmse"], np.zeros(2))
                    for value in ideal_result.values()
                )
                ideal_decode_ok = all(
                    not np.any(value["error_flags"])
                    for value in ideal_result.values()
                )
                if not ideal_ce_exact or not ideal_decode_ok or simulator._filters:
                    raise RuntimeError("Ideal-CSI noiseless link validation failed.")
                ideal_link_receipt = {
                    "payload_count": 2,
                    "candidate_count": len(ideal_result),
                    "all_crc_and_payload_checks_passed": True,
                    "ce_nmse_elementwise_zero": True,
                    "lmmse_filter_count": 0,
                }
            summary.update(
                {
                    "polar_mother_code_length": int(simulator.codec.encoder.n_polar),
                    "candidate_cces": simulator.pdcch_grid.candidate_cces.tolist(),
                    "candidate_bundles": simulator.pdcch_grid.candidate_bundles.tolist(),
                    "pilot_diagnostics": {
                        candidate.candidate_id: {
                            "receiver_covariance_mode": candidate.receiver_covariance_mode,
                            "rank": candidate.pilot_rank,
                            "condition_number": candidate.pilot_condition_number,
                            "ce_floor_nmse": candidate.ce_floor_nmse,
                        }
                        for candidate in simulator.candidates
                    },
                    "candidate_geometry": {
                        candidate.candidate_id: {
                            "delay_grid_coordinates": candidate.precoder.metadata.get(
                                "delay_grid_coordinates"
                            ),
                            "delay_seconds": candidate.precoder.metadata.get("delay_seconds"),
                            "q_seconds": candidate.precoder.metadata.get("q_seconds"),
                            "phase_denominator": candidate.precoder.metadata.get(
                                "phase_denominator"
                            ),
                        }
                        for candidate in simulator.candidates
                    },
                    "noiseless_decode": {
                        "payload_count": 2,
                        "all_crc_passed": True,
                        "all_payloads_equal": True,
                    },
                    "ideal_noiseless_link": ideal_link_receipt,
                    "multi_rx_channel_smoke": rx_channel_receipt,
                    "noise_variance_per_rx_branch": "10**(-snr_db/10), no n_rx scaling",
                }
            )
            rendered = json.dumps(summary, indent=2, ensure_ascii=False)
            if args.validation_output:
                destination = Path(args.validation_output)
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_text(rendered + "\n", encoding="utf-8")
            print(rendered)
            return
        output = Path(config["output_dir"])
        output.mkdir(parents=True, exist_ok=True)
        _archive_config(Path(args.config).resolve(), output)
        with open(output / "resolved_config.yaml", "w", encoding="utf-8") as handle:
            yaml.safe_dump(config, handle, allow_unicode=True, sort_keys=False)
        simulator.run()
        return

    if args.candidate:
        raise ValueError("--candidate is only supported by pdcch-cdd-bler-v2.")
    config = load_pdcch_bler_config(args.config)
    summary = validate_pdcch_bler_config(config)
    if args.stage == "validate":
        simulator = PDCCHBLERSimulator(config)
        summary.update(
            {
                "polar_mother_code_length": int(simulator.codec.encoder.n_polar),
                "candidate_cces": simulator.pdcch_grid.candidate_cces.tolist(),
                "candidate_bundles": simulator.pdcch_grid.candidate_bundles.tolist(),
            }
        )
        print(json.dumps(summary, indent=2, ensure_ascii=False))
        return

    output = Path(config["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    with open(output / "resolved_config.yaml", "w", encoding="utf-8") as handle:
        yaml.safe_dump(config, handle, allow_unicode=True, sort_keys=False)
    PDCCHBLERSimulator(config).run()


if __name__ == "__main__":
    main()
