from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any

import copy
import json
import math
import time

import numpy as np

from cdd_lls.core.config import ChannelConfig
from cdd_lls.phy.channel_tdl import generate_sionna_tdl_channel
from cdd_lls.phy.pdcch import (
    PDCCHGrid,
    PDCCHResourceConfig,
    build_pdcch_grid,
    build_reg_bundle_dft_precoder,
    pdcch_dmrs_symbols,
)
from cdd_lls.phy.pdcch_codec import SionnaPDCCHPolarCodec
from cdd_lls.phy.precoding import equivalent_channel
from cdd_lls.sim.stats import save_csv, save_json


DEFAULT_PDCCH_BLER_CONFIG: dict[str, Any] = {
    "schema": "pdcch-bler-v1",
    "output_dir": "outputs/pdcch",
    "seed": 42,
    "batch_size": 10,
    "antenna": {"n_tx": 8, "n_rx": 1},
    "resource": {
        "n_rb": 48,
        "duration_symbols": 2,
        "n_fft": 4096,
        "scs_khz": 30,
        "cyclic_prefix_length": 288,
        "cce_reg_mapping": "noninterleaved",
        "reg_bundle_size": 6,
        "interleaver_size": 2,
        "shift_index": 0,
        "aggregation_level": 4,
        "first_cce": 0,
    },
    "pdcch": {
        "payload_bits": 41,
        "coded_bits": 432,
        "crc_rnti": 0xFFFF,
        "scrambling_id": 0,
        "data_scrambling_rnti": 0,
        "polar_list_size": 8,
        "polar_decoder_type": "hybSCL",
        "polar_cpu_only": True,
        "coreset_start_symbol": 0,
        "slot_number": 0,
    },
    "precoding": {
        "scheme": "reg_bundle_dft_cycling",
        "cycling_order": list(range(8)),
    },
    "channel": {
        "backend": "sionna_tdl",
        "model": "tdl",
        "tdl_profile": "A",
        "delay_spread_ns": 100.0,
        "carrier_frequency_hz": 3.5e9,
        "ue_speed_kmh": 0.0,
        "num_sinusoids": 20,
        "normalize": True,
    },
    "receiver": {
        "channel_estimation": "bundle_ls_linear",
        "joint_time_estimation": True,
        "llr_clip": 50.0,
    },
    "simulation": {
        "snr_points_db": [-8.0, -6.0, -4.0, -2.0, 0.0],
        "min_trials_per_snr": 200,
        "target_errors": 50,
        "max_trials_per_snr": 1000,
        "target_bler": 0.01,
        "save_trial_error_flags": True,
    },
}


def _deep_merge(base: dict[str, Any], patch: dict[str, Any]) -> dict[str, Any]:
    out = copy.deepcopy(base)
    for key, value in (patch or {}).items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def load_pdcch_bler_config(path: str | Path) -> dict[str, Any]:
    import yaml

    with open(path, "r", encoding="utf-8") as handle:
        user = yaml.safe_load(handle) or {}
    config = _deep_merge(DEFAULT_PDCCH_BLER_CONFIG, user)
    validate_pdcch_bler_config(config)
    return config


def _resource_config(config: dict[str, Any]) -> PDCCHResourceConfig:
    fields = PDCCHResourceConfig.__dataclass_fields__
    return PDCCHResourceConfig(
        **{key: value for key, value in config["resource"].items() if key in fields}
    )


def _channel_config(config: dict[str, Any]) -> ChannelConfig:
    fields = ChannelConfig.__dataclass_fields__
    return ChannelConfig(
        **{key: value for key, value in config["channel"].items() if key in fields}
    )


def validate_pdcch_bler_config(config: dict[str, Any]) -> dict[str, Any]:
    if str(config.get("schema")) != "pdcch-bler-v1":
        raise ValueError("schema must be pdcch-bler-v1.")
    n_tx = int(config["antenna"]["n_tx"])
    n_rx = int(config["antenna"]["n_rx"])
    if n_tx not in (1, 2, 4, 8) or n_rx <= 0:
        raise ValueError("PDCCH antenna configuration requires n_tx in {1,2,4,8} and n_rx>0.")
    if int(config["batch_size"]) <= 0:
        raise ValueError("batch_size must be positive.")

    resource = _resource_config(config)
    grid = build_pdcch_grid(resource)
    pdcch = config["pdcch"]
    payload_bits = int(pdcch["payload_bits"])
    coded_bits = int(pdcch["coded_bits"])
    if coded_bits != grid.coded_bits:
        raise ValueError(
            f"pdcch.coded_bits={coded_bits} does not match 108*AL={grid.coded_bits}."
        )
    if payload_bits <= 0 or payload_bits > 140:
        raise ValueError("pdcch.payload_bits must lie in [1,140].")
    if payload_bits + 24 >= coded_bits:
        raise ValueError("payload_bits+24 must be smaller than coded_bits.")
    if coded_bits > 864:
        raise ValueError(
            "The validated PDCCH Polar adapter currently supports coded_bits<=864; "
            "use AL1/2/4/8."
        )
    if str(pdcch["polar_decoder_type"]) not in ("SCL", "hybSCL"):
        raise ValueError("pdcch.polar_decoder_type must be SCL or hybSCL.")

    if str(config["precoding"]["scheme"]).lower() != "reg_bundle_dft_cycling":
        raise ValueError("Only reg_bundle_dft_cycling is supported by the PDCCH branch.")
    order = [int(x) for x in config["precoding"].get("cycling_order", [])]
    if not order or any(x < 0 or x >= n_tx for x in order):
        raise ValueError("precoding.cycling_order must contain valid DFT vector indices.")

    channel = _channel_config(config)
    if str(channel.backend).lower() != "sionna_tdl":
        raise ValueError("The PDCCH branch currently requires channel.backend=sionna_tdl.")
    if str(channel.tdl_profile).upper() not in ("A", "B", "C", "D", "E"):
        raise ValueError("PDCCH Sionna TDL profile must be A, B, C, D, or E.")
    if float(channel.delay_spread_ns) <= 0 or float(channel.carrier_frequency_hz) <= 0:
        raise ValueError("Channel delay spread and carrier frequency must be positive.")
    if float(channel.ue_speed_kmh) < 0:
        raise ValueError("channel.ue_speed_kmh must be non-negative.")

    ce = str(config["receiver"]["channel_estimation"]).lower()
    if ce not in ("ideal", "bundle_ls_linear"):
        raise ValueError("receiver.channel_estimation must be ideal or bundle_ls_linear.")
    simulation = config["simulation"]
    snrs = [float(x) for x in simulation["snr_points_db"]]
    if not snrs:
        raise ValueError("simulation.snr_points_db must not be empty.")
    minimum = int(simulation["min_trials_per_snr"])
    maximum = int(simulation["max_trials_per_snr"])
    target_errors = int(simulation["target_errors"])
    if minimum <= 0 or maximum < minimum or target_errors <= 0:
        raise ValueError("Trial and error budgets are invalid.")
    return {
        "n_data_re": int(grid.resource_grid.n_data_re),
        "n_dmrs_re": int(grid.resource_grid.n_dmrs_re),
        "coded_bits": int(grid.coded_bits),
        "n_cce": int(grid.n_cce),
    }


def _stable_seed(*items: object) -> int:
    text = "|".join(str(item) for item in items)
    acc = 2166136261
    for char in text:
        acc ^= ord(char)
        acc = (acc * 16777619) % (2**32)
    return int(acc)


def _qpsk_modulate(bits: np.ndarray) -> np.ndarray:
    arr = np.asarray(bits, dtype=np.int8)
    if arr.shape[-1] % 2:
        raise ValueError("QPSK input width must be even.")
    pairs = arr.reshape(*arr.shape[:-1], -1, 2)
    return (
        (1.0 - 2.0 * pairs[..., 0]) + 1j * (1.0 - 2.0 * pairs[..., 1])
    ) / np.sqrt(2.0)


def _qpsk_llr(symbols: np.ndarray, noise_variance: np.ndarray) -> np.ndarray:
    z = np.asarray(symbols, dtype=np.complex128)
    no = np.maximum(np.asarray(noise_variance, dtype=np.float64), 1e-12)
    llr = np.empty(z.shape + (2,), dtype=np.float64)
    scale = -2.0 * np.sqrt(2.0) / no
    llr[..., 0] = scale * np.real(z)
    llr[..., 1] = scale * np.imag(z)
    return llr.reshape(z.shape[0], -1)


def _complex_noise(rng: np.random.Generator, shape: tuple[int, ...], variance: float) -> np.ndarray:
    return math.sqrt(float(variance) / 2.0) * (
        rng.normal(size=shape) + 1j * rng.normal(size=shape)
    )


def estimate_bundle_ls_linear(
    pilot_observations: np.ndarray,
    dmrs_symbols: np.ndarray,
    pdcch_grid: PDCCHGrid,
    joint_time: bool = True,
) -> np.ndarray:
    observations = np.asarray(pilot_observations, dtype=np.complex128)
    if observations.ndim != 3:
        raise ValueError("pilot_observations must have shape [batch,rx,n_dmrs_re].")
    ls = observations / np.asarray(dmrs_symbols, dtype=np.complex128)[None, None, :]
    data_count = int(pdcch_grid.resource_grid.n_data_re)
    estimate = np.empty((observations.shape[0], observations.shape[1], data_count), dtype=np.complex128)
    dmrs_sc = pdcch_grid.dmrs_local_subcarrier_indices
    data_sc = pdcch_grid.data_local_subcarrier_indices
    dmrs_symbols_idx = pdcch_grid.resource_grid.pilot_symbol_indices
    data_symbols_idx = pdcch_grid.resource_grid.data_symbol_indices

    for bundle in np.unique(pdcch_grid.data_bundle_indices):
        pilot_sel = np.flatnonzero(pdcch_grid.dmrs_bundle_indices == bundle)
        data_sel = np.flatnonzero(pdcch_grid.data_bundle_indices == bundle)
        if pilot_sel.size == 0 or data_sel.size == 0:
            raise RuntimeError(f"REG bundle {int(bundle)} lacks PDCCH pilot or data REs.")
        for batch in range(observations.shape[0]):
            for rx in range(observations.shape[1]):
                if bool(joint_time):
                    unique_sc = np.unique(dmrs_sc[pilot_sel])
                    values = np.asarray(
                        [
                            np.mean(ls[batch, rx, pilot_sel[dmrs_sc[pilot_sel] == sc]])
                            for sc in unique_sc
                        ],
                        dtype=np.complex128,
                    )
                    estimate[batch, rx, data_sel] = np.interp(
                        data_sc[data_sel], unique_sc, np.real(values)
                    ) + 1j * np.interp(data_sc[data_sel], unique_sc, np.imag(values))
                else:
                    for symbol in np.unique(data_symbols_idx[data_sel]):
                        local_data = data_sel[data_symbols_idx[data_sel] == symbol]
                        local_pilot = pilot_sel[dmrs_symbols_idx[pilot_sel] == symbol]
                        estimate[batch, rx, local_data] = np.interp(
                            data_sc[local_data], dmrs_sc[local_pilot], np.real(ls[batch, rx, local_pilot])
                        ) + 1j * np.interp(
                            data_sc[local_data], dmrs_sc[local_pilot], np.imag(ls[batch, rx, local_pilot])
                        )
    return estimate


def wilson_interval(errors: int, trials: int, z_value: float = 1.959963984540054) -> tuple[float, float]:
    if trials <= 0:
        return float("nan"), float("nan")
    p = float(errors) / float(trials)
    denom = 1.0 + z_value * z_value / trials
    center = (p + z_value * z_value / (2.0 * trials)) / denom
    half = z_value * math.sqrt(p * (1.0 - p) / trials + z_value * z_value / (4.0 * trials**2)) / denom
    return max(0.0, center - half), min(1.0, center + half)


class PDCCHBLERSimulator:
    def __init__(self, config: dict[str, Any]) -> None:
        validate_pdcch_bler_config(config)
        self.config = copy.deepcopy(config)
        self.resource_config = _resource_config(config)
        self.pdcch_grid = build_pdcch_grid(self.resource_config)
        self.channel_config = _channel_config(config)
        self.precoder = build_reg_bundle_dft_precoder(
            self.pdcch_grid,
            n_tx=int(config["antenna"]["n_tx"]),
            cycling_order=config["precoding"]["cycling_order"],
        )
        pdcch = config["pdcch"]
        self.codec = SionnaPDCCHPolarCodec(
            payload_bits=int(pdcch["payload_bits"]),
            coded_bits=int(pdcch["coded_bits"]),
            crc_rnti=int(pdcch["crc_rnti"]),
            scrambling_id=int(pdcch["scrambling_id"]),
            data_scrambling_rnti=int(pdcch["data_scrambling_rnti"]),
            list_size=int(pdcch["polar_list_size"]),
            decoder_type=str(pdcch["polar_decoder_type"]),
            cpu_only=bool(pdcch["polar_cpu_only"]),
        )
        self.dmrs_symbols = pdcch_dmrs_symbols(
            self.pdcch_grid,
            scrambling_id=int(pdcch["scrambling_id"]),
            slot_number=int(pdcch["slot_number"]),
            coreset_start_symbol=int(pdcch["coreset_start_symbol"]),
        )

    def run_batch(self, snr_db: float, absolute_start: int, batch_size: int) -> dict[str, np.ndarray]:
        base_seed = int(self.config["seed"])
        rng = np.random.default_rng(
            _stable_seed(base_seed, float(snr_db), int(absolute_start), "payload_noise")
        )
        payload = rng.integers(
            0,
            2,
            size=(int(batch_size), int(self.codec.payload_bits)),
            dtype=np.int8,
        )
        coded = self.codec.encode(payload)
        data_symbols = _qpsk_modulate(coded)

        channel = generate_sionna_tdl_channel(
            grid=self.pdcch_grid.resource_grid,
            channel=self.channel_config,
            n_tx=int(self.config["antenna"]["n_tx"]),
            n_rx=int(self.config["antenna"]["n_rx"]),
            batch_size=int(batch_size),
            seed=_stable_seed(base_seed, float(snr_db), int(absolute_start), "channel"),
        )
        effective = equivalent_channel(channel.H, self.precoder.C)
        grid = self.pdcch_grid.resource_grid
        data_true = effective[
            :,
            :,
            grid.data_symbol_indices,
            self.pdcch_grid.data_local_subcarrier_indices,
        ]
        pilot_true = effective[
            :,
            :,
            grid.pilot_symbol_indices,
            self.pdcch_grid.dmrs_local_subcarrier_indices,
        ]
        noise_variance = float(10.0 ** (-float(snr_db) / 10.0))
        data_received = data_true * data_symbols[:, None, :] + _complex_noise(
            rng, data_true.shape, noise_variance
        )
        pilot_received = pilot_true * self.dmrs_symbols[None, None, :] + _complex_noise(
            rng, pilot_true.shape, noise_variance
        )

        ce_method = str(self.config["receiver"]["channel_estimation"]).lower()
        if ce_method == "ideal":
            data_estimate = data_true
        else:
            data_estimate = estimate_bundle_ls_linear(
                pilot_received,
                self.dmrs_symbols,
                self.pdcch_grid,
                joint_time=bool(self.config["receiver"]["joint_time_estimation"]),
            )

        denominator = np.maximum(np.sum(np.abs(data_estimate) ** 2, axis=1), 1e-12)
        equalized = np.sum(np.conj(data_estimate) * data_received, axis=1) / denominator
        effective_noise = noise_variance / denominator
        llr = _qpsk_llr(equalized, effective_noise)
        llr = np.clip(
            llr,
            -float(self.config["receiver"]["llr_clip"]),
            float(self.config["receiver"]["llr_clip"]),
        )
        decoded = self.codec.decode(llr)
        payload_error = np.any(decoded.payload_bits != payload, axis=1)
        error_flags = payload_error | ~decoded.crc_status
        ce_nmse = np.sum(np.abs(data_estimate - data_true) ** 2, axis=(1, 2)) / np.maximum(
            np.sum(np.abs(data_true) ** 2, axis=(1, 2)), 1e-30
        )
        return {
            "error_flags": error_flags.astype(bool),
            "crc_status": decoded.crc_status.astype(bool),
            "payload_error": payload_error.astype(bool),
            "ce_nmse": ce_nmse.astype(np.float64),
        }

    def run(self) -> list[dict[str, Any]]:
        output = Path(self.config["output_dir"])
        output.mkdir(parents=True, exist_ok=True)
        flags_dir = output / "trial_error_flags"
        flags_dir.mkdir(parents=True, exist_ok=True)
        log_path = output / "run.log"
        log_lines: list[str] = []

        def log(message: str) -> None:
            line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {message}"
            print(line, flush=True)
            log_lines.append(line)
            log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")

        rows: list[dict[str, Any]] = []
        simulation = self.config["simulation"]
        minimum = int(simulation["min_trials_per_snr"])
        maximum = int(simulation["max_trials_per_snr"])
        target_errors = int(simulation["target_errors"])
        batch_size = int(self.config["batch_size"])
        started = time.time()
        for snr_db in [float(x) for x in simulation["snr_points_db"]]:
            flags: list[np.ndarray] = []
            nmse: list[np.ndarray] = []
            trials = 0
            errors = 0
            while trials < maximum:
                current = min(batch_size, maximum - trials)
                result = self.run_batch(snr_db, trials + 1, current)
                flags.append(result["error_flags"])
                nmse.append(result["ce_nmse"])
                trials += current
                errors += int(np.sum(result["error_flags"]))
                if trials >= minimum and errors >= target_errors:
                    break
            all_flags = np.concatenate(flags).astype(bool)
            all_nmse = np.concatenate(nmse).astype(np.float64)
            low, high = wilson_interval(errors, trials)
            row = {
                "snr_db": snr_db,
                "trials": trials,
                "errors": errors,
                "bler": float(errors) / float(trials),
                "bler_wilson95_lo": low,
                "bler_wilson95_hi": high,
                "ce_nmse_db": float(10.0 * np.log10(max(float(np.mean(all_nmse)), 1e-30))),
                "channel_estimation": str(self.config["receiver"]["channel_estimation"]),
                "payload_bits": int(self.codec.payload_bits),
                "crc_bits": 24,
                "polar_input_bits": int(self.codec.payload_bits + 24),
                "coded_bits": int(self.codec.coded_bits),
                "aggregation_level": int(self.resource_config.aggregation_level),
                "data_re": int(self.pdcch_grid.resource_grid.n_data_re),
                "dmrs_re": int(self.pdcch_grid.resource_grid.n_dmrs_re),
                "n_tx": int(self.config["antenna"]["n_tx"]),
                "n_rx": int(self.config["antenna"]["n_rx"]),
                "reg_bundle_size": int(self.resource_config.reg_bundle_size),
            }
            rows.append(row)
            if bool(simulation.get("save_trial_error_flags", True)):
                tag = str(snr_db).replace("-", "m").replace(".", "p")
                np.save(flags_dir / f"snr_{tag}_error_flags.npy", all_flags)
                np.save(flags_dir / f"snr_{tag}_ce_nmse.npy", all_nmse)
            log(
                f"SNR={snr_db:g}dB errors={errors} trials={trials} "
                f"BLER={row['bler']:.6g} CE_NMSE={row['ce_nmse_db']:.3f}dB"
            )

        save_csv(rows, output / "bler_points.csv")
        metadata = {
            "schema": "pdcch-bler-result-v1",
            "elapsed_seconds": float(time.time() - started),
            "resource_config": asdict(self.resource_config),
            "candidate_cces": self.pdcch_grid.candidate_cces.tolist(),
            "candidate_regs": self.pdcch_grid.candidate_regs.tolist(),
            "candidate_bundles": self.pdcch_grid.candidate_bundles.tolist(),
            "precoder": self.precoder.metadata,
            "polar": {
                "payload_bits": int(self.codec.payload_bits),
                "crc_bits": 24,
                "coded_bits": int(self.codec.coded_bits),
                "mother_code_length": int(self.codec.encoder.n_polar),
                "sionna_target_bits": int(self.codec.sionna_target_bits),
                "rate_matching_mode": self.codec.rate_matching_mode,
                "external_repetition": bool(self.codec.external_repetition),
                "list_size": int(self.codec.list_size),
                "decoder_type": self.codec.decoder_type,
                "crc_rnti": int(self.codec.crc_rnti),
                "dci_affine_crc_offset_applied": True,
            },
            "snr_definition": "Es/N0 per active unit-energy PDCCH RE at unit total transmit power",
            "rows": rows,
        }
        save_json(metadata, output / "run_metadata.json")
        self._plot(rows, output / "bler_curve.png")
        return rows

    @staticmethod
    def _plot(rows: list[dict[str, Any]], path: Path) -> None:
        import matplotlib.pyplot as plt

        snr = np.asarray([float(row["snr_db"]) for row in rows])
        bler = np.asarray([float(row["bler"]) for row in rows])
        trials = np.asarray([int(row["trials"]) for row in rows])
        display = np.where(bler > 0, bler, 0.5 / np.maximum(trials, 1))
        fig, axis = plt.subplots(figsize=(8.2, 5.8))
        axis.semilogy(snr, display, marker="o", linewidth=2.5, markersize=8)
        axis.axhline(0.01, color="black", linestyle="--", linewidth=1.8, label="BLER=0.01")
        axis.set_xlabel("SNR (dB)", fontsize=16)
        axis.set_ylabel("DCI BLER", fontsize=16)
        axis.set_title("PDCCH DCI BLER", fontsize=17)
        axis.grid(True, which="both", alpha=0.35)
        axis.tick_params(labelsize=14)
        axis.legend(fontsize=14)
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
