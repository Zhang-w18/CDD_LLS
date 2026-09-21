from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

import copy
import csv
import json
import math
import time

import numpy as np

from cdd_lls.core.config import ChannelConfig
from cdd_lls.phy.channel_tdl import generate_sionna_tdl_channel_active
from cdd_lls.phy.estimators import (
    TDLTimeFrequencyCovariance,
    build_frequency_rmmse_filter,
    build_prg_frequency_rmmse_filter,
    build_prg_time_frequency_rmmse_filter,
    build_time_frequency_rmmse_filter,
    linear_estimator_closed_form_nmse,
    tdl_active_frequency_covariance,
    tdl_active_time_frequency_covariance,
)
from cdd_lls.phy.pdcch import (
    PDCCHGrid,
    PDCCHResourceConfig,
    build_pdcch_grid,
    build_reg_bundle_dft_precoder,
    pdcch_dmrs_symbols,
)
from cdd_lls.phy.pdcch_codec import SionnaPDCCHPolarCodec
from cdd_lls.phy.precoding import PrecoderResult, equivalent_channel
from cdd_lls.sim.pdcch import (
    _channel_config,
    _complex_noise,
    _qpsk_llr,
    _qpsk_modulate,
    _resource_config,
    _stable_seed,
    wilson_interval,
)
from cdd_lls.sim.stats import save_csv, save_json


SCHEMA = "pdcch-cdd-bler-v2"
TOTAL_TRANSMIT_POWER = 1.0
RX_CORRELATION_MODEL = "identity_no_spatial_correlation"


@dataclass(frozen=True)
class CandidateRuntime:
    candidate_id: str
    family: str
    scheme: str
    label: str
    style: dict[str, Any]
    receiver_covariance_mode: str
    snr_points_db: tuple[float, ...] | None
    precoder: PrecoderResult
    covariance: np.ndarray
    time_frequency_covariance: TDLTimeFrequencyCovariance
    pilot_rank: int
    pilot_condition_number: float
    ce_floor_nmse: float


def _candidate_band(grid: PDCCHGrid) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    occupied = np.unique(
        np.concatenate(
            [grid.data_local_subcarrier_indices, grid.dmrs_local_subcarrier_indices]
        )
    )
    expected = (
        72 * int(grid.config.aggregation_level) // int(grid.config.duration_symbols)
    )
    if occupied.tolist() != list(range(expected)):
        raise ValueError(
            "Plan-031 requires first CCE 0 and a contiguous candidate band starting at RB 0."
        )
    pilot = np.asarray(grid.dmrs_local_subcarrier_indices, dtype=np.int64)
    data = np.asarray(grid.data_local_subcarrier_indices, dtype=np.int64)
    return occupied, pilot, data


def _validate_candidate(
    candidate: dict[str, Any], *, n_tx: int, k_active: int, n_p: int
) -> dict[str, Any]:
    candidate_id = str(candidate.get("candidate_id", "")).strip()
    if not candidate_id:
        raise ValueError("Every candidate requires candidate_id.")
    scheme = str(candidate.get("scheme", "")).lower()
    default_mode = "physical_prg" if scheme == "reg_bundle_dft_cycling" else "matched_effective"
    receiver_mode = str(candidate.get("receiver_covariance_mode", default_mode)).lower()
    if scheme == "cdd":
        if receiver_mode not in ("matched_effective", "physical_fullband"):
            raise ValueError(
                f"{candidate_id}: CDD receiver_covariance_mode must be "
                "matched_effective or physical_fullband."
            )
        delays = np.asarray(candidate.get("delay_grid_coordinates", []), dtype=np.float64)
        if delays.shape != (n_tx,) or not np.all(np.isfinite(delays)):
            raise ValueError(f"{candidate_id}: CDD delay_grid_coordinates must contain {n_tx} finite values.")
        if np.any(delays < 0.0) or np.any(delays >= float(k_active)):
            raise ValueError(f"{candidate_id}: CDD coordinates must satisfy 0<=j<K.")
        residues = np.mod(delays, float(n_p))
        pairwise = np.abs(residues[:, None] - residues[None, :])
        circular = np.minimum(pairwise, float(n_p) - pairwise)
        circular += np.eye(n_tx) * float(n_p)
        raw_allow_duplicates = candidate.get("allow_duplicate_delays", False)
        if not isinstance(raw_allow_duplicates, (bool, np.bool_)):
            raise ValueError(f"{candidate_id}: allow_duplicate_delays must be boolean.")
        allow_duplicates = bool(raw_allow_duplicates)
        if float(np.min(circular)) <= 1e-10 and not allow_duplicates:
            raise ValueError(f"{candidate_id}: folded CDD residues are not distinct.")
        strict_sidon = bool(candidate.get("strict_sidon", False))
        pair_sum_count = None
        if strict_sidon:
            rounded = np.rint(delays).astype(np.int64)
            if not np.allclose(delays, rounded, rtol=0.0, atol=1e-12):
                raise ValueError(f"{candidate_id}: strict Sidon coordinates must be integers.")
            pair_sums = [
                int((rounded[left] + rounded[right]) % k_active)
                for left in range(n_tx)
                for right in range(left, n_tx)
            ]
            pair_sum_count = len(pair_sums)
            if len(set(pair_sums)) != pair_sum_count:
                raise ValueError(
                    f"{candidate_id}: unordered pair sums are not unique modulo K."
                )
    elif scheme == "reg_bundle_dft_cycling":
        if receiver_mode != "physical_prg":
            raise ValueError(
                f"{candidate_id}: cycling receiver_covariance_mode must be physical_prg."
            )
        order = [int(value) for value in candidate.get("cycling_order", [])]
        if len(order) != int(candidate["aggregation_level"] if "aggregation_level" in candidate else 0):
            raise ValueError(
                f"{candidate_id}: cycling_order must contain one vector per occupied 6-REG bundle."
            )
        if len(set(order)) != len(order) or any(value < 0 or value >= n_tx for value in order):
            raise ValueError(f"{candidate_id}: cycling_order must contain distinct valid DFT indices.")
    else:
        raise ValueError(f"{candidate_id}: unsupported scheme {scheme!r}.")
    if "snr_points_db" in candidate:
        snrs = [float(value) for value in candidate["snr_points_db"]]
        if not snrs or any(not np.isfinite(value) for value in snrs) or len(set(snrs)) != len(snrs):
            raise ValueError(f"{candidate_id}: snr_points_db must contain distinct finite values.")
    return {
        "strict_sidon": bool(candidate.get("strict_sidon", False)),
        "pair_sum_count": pair_sum_count if scheme == "cdd" else None,
        "fold_residue_count": (
            len(np.unique(np.round(residues, decimals=12))) if scheme == "cdd" else None
        ),
        "allow_duplicate_delays": allow_duplicates if scheme == "cdd" else False,
    }


def validate_pdcch_cdd_config(config: dict[str, Any]) -> dict[str, Any]:
    if str(config.get("schema")) != SCHEMA:
        raise ValueError(f"schema must be {SCHEMA}.")
    resource = _resource_config(config)
    grid = build_pdcch_grid(resource)
    n_tx = int(config["antenna"]["n_tx"])
    raw_n_rx = config["antenna"]["n_rx"]
    if isinstance(raw_n_rx, (bool, np.bool_)) or not isinstance(
        raw_n_rx, (int, np.integer)
    ):
        raise ValueError("antenna.n_rx must be a positive integer.")
    n_rx = int(raw_n_rx)
    if n_rx <= 0:
        raise ValueError("antenna.n_rx must be a positive integer.")
    duration = int(resource.duration_symbols)
    aggregation = int(resource.aggregation_level)
    if (
        int(resource.n_rb) != 48
        or str(resource.cce_reg_mapping).lower() != "noninterleaved"
        or int(resource.reg_bundle_size) != 6
        or int(resource.first_cce) != 0
    ):
        raise ValueError(
            "The PDCCH CDD entry requires 48 RB, non-interleaved L=6, and first CCE 0."
        )
    legacy_scene = duration == 1 and n_tx == 8 and aggregation in (2, 4, 8)
    c300_scene = duration == 2 and n_tx == 4 and aggregation in (1, 2, 4)
    plan036_scene = duration == 1 and n_tx == 4 and aggregation in (1, 2)
    if not (legacy_scene or c300_scene or plan036_scene):
        raise ValueError(
            "The PDCCH CDD entry supports 8Tx/1-symbol/AL2,4,8, "
            "4Tx/2-symbol/AL1,2,4, or plan-036 4Tx/1-symbol/AL1,2."
        )
    pdcch = config["pdcch"]
    payload_bits = int(pdcch["payload_bits"])
    if payload_bits not in (40, 41) or int(pdcch["coded_bits"]) != grid.coded_bits:
        raise ValueError("The PDCCH CDD entry requires A in {40,41} and coded_bits=108*AL.")
    channel_estimation = str(config["receiver"]["channel_estimation"]).lower()
    if channel_estimation not in (
        "frequency_lmmse",
        "matched_lmmse",
        "ideal",
    ):
        raise ValueError(
            "Plan-031 requires receiver.channel_estimation=frequency_lmmse, ideal "
            "(legacy alias: matched_lmmse)."
        )
    channel = _channel_config(config)
    common_channel_ok = str(channel.backend).lower() == "sionna_tdl"
    legacy_channel_ok = (
        str(channel.tdl_profile).upper() == "A"
        and abs(float(channel.delay_spread_ns) - 100.0) <= 1e-12
        and float(channel.ue_speed_kmh) == 0.0
    )
    c300_channel_ok = (
        str(channel.tdl_profile).upper() == "C"
        and abs(float(channel.delay_spread_ns) - 300.0) <= 1e-12
        and abs(float(channel.ue_speed_kmh) - 3.0) <= 1e-12
        and abs(float(channel.carrier_frequency_hz) - 4.0e9) <= 1e-3
        and int(channel.num_sinusoids) == 20
    )
    plan036_channel_ok = (
        str(channel.tdl_profile).upper() == "C"
        and abs(float(channel.delay_spread_ns) - 300.0) <= 1e-12
        and float(channel.ue_speed_kmh) > 0.0
        and abs(float(channel.carrier_frequency_hz) - 4.0e9) <= 1e-3
        and int(channel.num_sinusoids) == 20
    )
    if not common_channel_ok or not (
        (legacy_scene and legacy_channel_ok)
        or (c300_scene and c300_channel_ok)
        or (plan036_scene and plan036_channel_ok)
    ):
        raise ValueError(
            "The channel must match its A100 legacy, C300 two-symbol, or plan-036 scene."
        )
    snrs = [float(value) for value in config["simulation"]["snr_points_db"]]
    if not snrs or any(not np.isfinite(value) for value in snrs) or len(set(snrs)) != len(snrs):
        raise ValueError("simulation.snr_points_db must contain distinct finite values.")
    minimum = int(config["simulation"]["min_trials_per_snr"])
    maximum = int(config["simulation"]["max_trials_per_snr"])
    target_errors = int(config["simulation"]["target_errors"])
    progress_every_batches = int(
        config["simulation"].get("progress_every_batches", 0)
    )
    if (
        minimum <= 0
        or maximum < minimum
        or target_errors <= 0
        or int(config["batch_size"]) <= 0
        or progress_every_batches < 0
    ):
        raise ValueError("Trial, error, and batch budgets are invalid.")
    candidates = list(config.get("candidates", []))
    if not candidates:
        raise ValueError("Plan-031 requires a non-empty candidates list.")
    ids = [str(candidate.get("candidate_id", "")) for candidate in candidates]
    if len(set(ids)) != len(ids):
        raise ValueError("candidate_id values must be unique.")
    occupied, pilot, data = _candidate_band(grid)
    k_active = int(len(occupied))
    n_p = int(k_active // 4)
    candidate_diagnostics = {}
    for candidate in candidates:
        candidate_for_validation = dict(candidate)
        candidate_for_validation["aggregation_level"] = aggregation
        candidate_diagnostics[str(candidate["candidate_id"])] = _validate_candidate(
            candidate_for_validation, n_tx=n_tx, k_active=k_active, n_p=n_p
        )
    return {
        "aggregation_level": int(resource.aggregation_level),
        "occupied_rb": int(k_active // 12),
        "k_active": k_active,
        "n_p": n_p,
        "n_data_re": int(len(data)),
        "n_dmrs_re": int(len(pilot)),
        "coded_bits": int(grid.coded_bits),
        "n_rx": n_rx,
        "layers": 1,
        "total_transmit_power": TOTAL_TRANSMIT_POWER,
        "rx_correlation_model": RX_CORRELATION_MODEL,
        "channel_estimation": channel_estimation,
        "candidate_ids": ids,
        "candidate_constraint_diagnostics": candidate_diagnostics,
    }


def load_pdcch_cdd_config(path: str | Path) -> dict[str, Any]:
    import yaml

    with open(path, "r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    validate_pdcch_cdd_config(config)
    return config


def _build_cdd_precoder(
    grid: PDCCHGrid,
    delays: list[float],
    n_tx: int,
) -> PrecoderResult:
    k_active = len(_candidate_band(grid)[0])
    j = np.asarray(delays, dtype=np.float64)
    local = np.arange(int(grid.resource_grid.n_sc), dtype=np.float64)
    matrix = np.exp(-1j * 2.0 * np.pi * local[:, None] * j[None, :] / float(k_active))
    matrix /= np.sqrt(float(n_tx))
    q_seconds = 1.0 / (float(k_active) * float(grid.resource_grid.scs_khz) * 1e3)
    return PrecoderResult(
        C=matrix.astype(np.complex128),
        label="PDCCH_CDD_ACTIVE_BAND",
        metadata={
            "normalized": True,
            "vector_power": 1.0,
            "phase_reference": "first_candidate_subcarrier",
            "phase_denominator": k_active,
            "delay_grid_coordinates": [float(value) for value in j],
            "delay_seconds": [float(value) * q_seconds for value in j],
            "q_seconds": q_seconds,
        },
    )


def _conditional_nmse_floor(
    true_covariance: np.ndarray,
    pilot: np.ndarray,
    data: np.ndarray,
    assumed_covariance: np.ndarray | None = None,
) -> float:
    assumed = true_covariance if assumed_covariance is None else assumed_covariance
    assumed_pp = assumed[np.ix_(pilot, pilot)]
    assumed_dp = assumed[np.ix_(data, pilot)]
    weights = assumed_dp @ np.linalg.pinv(assumed_pp, rcond=1e-12)
    true_pp = true_covariance[np.ix_(pilot, pilot)]
    true_dp = true_covariance[np.ix_(data, pilot)]
    true_dd = true_covariance[np.ix_(data, data)]
    error = (
        true_dd
        - weights @ true_dp.conj().T
        - true_dp @ weights.conj().T
        + weights @ true_pp @ weights.conj().T
    )
    return max(float(np.real(np.trace(error))) / float(np.real(np.trace(true_dd))), 0.0)


def _pilot_diagnostics(precoder: np.ndarray, pilot: np.ndarray) -> tuple[int, float]:
    matrix = np.asarray(precoder, dtype=np.complex128)[pilot, :]
    singular = np.linalg.svd(matrix, compute_uv=False)
    tolerance = max(matrix.shape) * np.finfo(np.float64).eps * singular[0]
    rank = int(np.sum(singular > tolerance))
    condition = float("inf") if singular[-1] <= tolerance else float(singular[0] / singular[-1])
    return rank, condition


class PDCCHCDDBLERSimulator:
    def __init__(self, config: dict[str, Any]) -> None:
        self.validation = validate_pdcch_cdd_config(config)
        self.config = copy.deepcopy(config)
        self.resource_config: PDCCHResourceConfig = _resource_config(config)
        self.pdcch_grid = build_pdcch_grid(self.resource_config)
        self.channel_config: ChannelConfig = _channel_config(config)
        self.channel_estimation = str(config["receiver"]["channel_estimation"]).lower()
        self.ideal_csi = self.channel_estimation == "ideal"
        self.occupied, self.pilot_local, self.data_local = _candidate_band(self.pdcch_grid)
        self.k_active = int(len(self.occupied))
        self.pilot_band = self.pilot_local.copy()
        self.data_band = self.data_local.copy()
        full_grid = self.pdcch_grid.resource_grid
        band_subcarriers = full_grid.subcarrier_indices[self.occupied]
        band_active_fft = full_grid.active_fft_indices[self.occupied]
        self.band_grid = replace(
            full_grid,
            n_sc=self.k_active,
            subcarrier_indices=band_subcarriers,
            active_fft_indices=band_active_fft,
            pilot_subcarriers=np.unique(full_grid.pilot_subcarrier_indices),
        )
        n_tx = int(config["antenna"]["n_tx"])
        self.physical_time_frequency_covariance = tdl_active_time_frequency_covariance(
            self.band_grid, self.channel_config
        )
        self.physical_covariance = self.physical_time_frequency_covariance.frequency
        self.candidates: list[CandidateRuntime] = []
        for spec in config["candidates"]:
            scheme = str(spec["scheme"]).lower()
            default_mode = "physical_prg" if scheme == "reg_bundle_dft_cycling" else "matched_effective"
            receiver_mode = str(spec.get("receiver_covariance_mode", default_mode)).lower()
            if scheme == "cdd":
                precoder = _build_cdd_precoder(
                    self.pdcch_grid,
                    [float(value) for value in spec["delay_grid_coordinates"]],
                    n_tx,
                )
            else:
                precoder = build_reg_bundle_dft_precoder(
                    self.pdcch_grid,
                    n_tx=n_tx,
                    cycling_order=spec["cycling_order"],
                )
            band_precoder = precoder.C[self.occupied, :]
            row_power = np.sum(np.abs(band_precoder) ** 2, axis=1)
            if not np.allclose(
                row_power,
                TOTAL_TRANSMIT_POWER,
                rtol=1e-12,
                atol=1e-12,
            ):
                raise RuntimeError(
                    f"{spec['candidate_id']}: precoder row power is not normalized to "
                    f"{TOTAL_TRANSMIT_POWER:g}."
                )
            precoder = PrecoderResult(
                C=band_precoder,
                label=precoder.label,
                metadata={
                    **precoder.metadata,
                    "runtime_band_subcarriers": self.k_active,
                    "total_transmit_power_per_re": TOTAL_TRANSMIT_POWER,
                },
            )
            factor = band_precoder @ band_precoder.conj().T
            covariance = self.physical_covariance * factor
            true_time_frequency = TDLTimeFrequencyCovariance(
                time=self.physical_time_frequency_covariance.time,
                frequency=covariance,
                covariance_type="tdl_effective_known_precoder",
            )
            rank, condition = _pilot_diagnostics(band_precoder, self.pilot_band)
            assumed_covariance = (
                self.physical_covariance
                if receiver_mode == "physical_fullband"
                else covariance
            )
            if self.ideal_csi:
                floor = 0.0
            elif int(self.resource_config.duration_symbols) == 1:
                floor = _conditional_nmse_floor(
                    covariance,
                    self.pilot_band,
                    self.data_band,
                    assumed_covariance=assumed_covariance,
                )
            else:
                assumed_tf = TDLTimeFrequencyCovariance(
                    time=self.physical_time_frequency_covariance.time,
                    frequency=assumed_covariance,
                    covariance_type=f"{receiver_mode}_zero_noise",
                )
                if receiver_mode == "physical_prg":
                    floor_filter = build_prg_time_frequency_rmmse_filter(
                        self.band_grid,
                        assumed_tf,
                        prg_size_subcarriers=(
                            12
                            * int(self.resource_config.reg_bundle_size)
                            // int(self.resource_config.duration_symbols)
                        ),
                        noise_variance=0.0,
                        diagonal_loading=1e-12,
                    )
                else:
                    floor_filter = build_time_frequency_rmmse_filter(
                        self.band_grid,
                        assumed_tf,
                        noise_variance=0.0,
                        diagonal_loading=1e-12,
                    )
                floor = linear_estimator_closed_form_nmse(
                    self.band_grid,
                    floor_filter,
                    true_time_frequency,
                    noise_variance=0.0,
                )
            self.candidates.append(
                CandidateRuntime(
                    candidate_id=str(spec["candidate_id"]),
                    family=str(spec.get("family", spec["candidate_id"])),
                    scheme=scheme,
                    label=str(spec.get("label", spec["candidate_id"])),
                    style=dict(spec.get("style", {})),
                    receiver_covariance_mode=receiver_mode,
                    snr_points_db=(
                        tuple(float(value) for value in spec["snr_points_db"])
                        if "snr_points_db" in spec
                        else None
                    ),
                    precoder=precoder,
                    covariance=covariance,
                    time_frequency_covariance=true_time_frequency,
                    pilot_rank=rank,
                    pilot_condition_number=condition,
                    ce_floor_nmse=floor,
                )
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
        self._filters: dict[tuple[str, float], Any] = {}

    def _filter(self, candidate: CandidateRuntime, noise_variance: float) -> Any:
        if self.ideal_csi:
            raise RuntimeError("The ideal-CSI branch must not construct or call an LMMSE filter.")
        key = (candidate.candidate_id, float(noise_variance))
        if key not in self._filters:
            if int(self.resource_config.duration_symbols) > 1:
                if candidate.receiver_covariance_mode == "physical_prg":
                    value = build_prg_time_frequency_rmmse_filter(
                        self.band_grid,
                        self.physical_time_frequency_covariance,
                        prg_size_subcarriers=(
                            12
                            * int(self.resource_config.reg_bundle_size)
                            // int(self.resource_config.duration_symbols)
                        ),
                        noise_variance=float(noise_variance),
                    )
                elif candidate.receiver_covariance_mode == "physical_fullband":
                    value = build_time_frequency_rmmse_filter(
                        self.band_grid,
                        self.physical_time_frequency_covariance,
                        noise_variance=float(noise_variance),
                    )
                else:
                    value = build_time_frequency_rmmse_filter(
                        self.band_grid,
                        candidate.time_frequency_covariance,
                        noise_variance=float(noise_variance),
                    )
            elif candidate.receiver_covariance_mode == "physical_prg":
                value = build_prg_frequency_rmmse_filter(
                    self.physical_covariance,
                    self.pilot_band,
                    prg_size_subcarriers=72,
                    noise_variance=float(noise_variance),
                )
            elif candidate.receiver_covariance_mode == "physical_fullband":
                value = build_frequency_rmmse_filter(
                    self.physical_covariance,
                    self.pilot_band,
                    noise_variance=float(noise_variance),
                )
            else:
                value = build_frequency_rmmse_filter(
                    candidate.covariance,
                    self.pilot_band,
                    noise_variance=float(noise_variance),
                )
            self._filters[key] = value
        return self._filters[key]

    @staticmethod
    def _coherent_mrc(
        data_received: np.ndarray,
        data_estimate: np.ndarray,
        noise_variance: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        denominator = np.maximum(np.sum(np.abs(data_estimate) ** 2, axis=1), 1e-12)
        equalized = np.sum(np.conj(data_estimate) * data_received, axis=1) / denominator
        effective_noise = float(noise_variance) / denominator
        return equalized, effective_noise, denominator

    def run_batch(
        self, snr_db: float, absolute_start: int, batch_size: int
    ) -> dict[str, dict[str, np.ndarray]]:
        base_seed = int(self.config["seed"])
        if int(self.codec.payload_bits) == 40:
            base_seed = _stable_seed(
                base_seed,
                str(self.config.get("random_stream_namespace", "default")),
            )
        rng = np.random.default_rng(
            _stable_seed(base_seed, float(snr_db), int(absolute_start), "shared_payload_noise")
        )
        payload = rng.integers(
            0, 2, size=(int(batch_size), int(self.codec.payload_bits)), dtype=np.int8
        )
        coded = self.codec.encode(payload)
        data_symbols = _qpsk_modulate(coded)
        channel = generate_sionna_tdl_channel_active(
            grid=self.band_grid,
            channel=self.channel_config,
            n_tx=int(self.config["antenna"]["n_tx"]),
            n_rx=int(self.config["antenna"]["n_rx"]),
            batch_size=int(batch_size),
            seed=_stable_seed(base_seed, float(snr_db), int(absolute_start), "shared_channel"),
        )
        noise_variance = float(10.0 ** (-float(snr_db) / 10.0))
        data_unit_noise = _complex_noise(
            rng,
            (int(batch_size), int(self.config["antenna"]["n_rx"]), len(self.data_local)),
            1.0,
        )
        pilot_unit_noise = None
        if not self.ideal_csi:
            pilot_unit_noise = _complex_noise(
                rng,
                (int(batch_size), int(self.config["antenna"]["n_rx"]), len(self.pilot_local)),
                1.0,
            )
        candidate_state: list[tuple[CandidateRuntime, np.ndarray, np.ndarray]] = []
        all_llr: list[np.ndarray] = []
        for candidate in self.candidates:
            effective = equivalent_channel(channel.H, candidate.precoder.C)
            data_true = effective[
                :, :, self.band_grid.data_symbol_indices, self.data_band
            ]
            data_received = (
                data_true * data_symbols[:, None, :] + math.sqrt(noise_variance) * data_unit_noise
            )
            if self.ideal_csi:
                data_estimate = data_true
            else:
                pilot_true = effective[
                    :, :, self.band_grid.pilot_symbol_indices, self.pilot_band
                ]
                pilot_received = (
                    pilot_true * self.dmrs_symbols[None, None, :]
                    + math.sqrt(noise_variance) * pilot_unit_noise
                )
                ls = pilot_received / self.dmrs_symbols[None, None, :]
                estimator = self._filter(candidate, noise_variance)
                if int(self.resource_config.duration_symbols) > 1:
                    data_estimate = estimator.estimate_data(ls)
                else:
                    full_estimate = estimator.estimate_full_band(ls)
                    data_estimate = full_estimate[:, :, self.data_band]
            equalized, effective_noise, denominator = self._coherent_mrc(
                data_received, data_estimate, noise_variance
            )
            llr = np.clip(
                _qpsk_llr(equalized, effective_noise),
                -float(self.config["receiver"]["llr_clip"]),
                float(self.config["receiver"]["llr_clip"]),
            )
            all_llr.append(llr)
            ce_nmse = np.sum(np.abs(data_estimate - data_true) ** 2, axis=(1, 2)) / np.maximum(
                np.sum(np.abs(data_true) ** 2, axis=(1, 2)), 1e-30
            )
            candidate_state.append(
                (candidate, data_true, ce_nmse, llr, denominator, effective_noise)
            )
        decoded = self.codec.decode(np.concatenate(all_llr, axis=0))
        out: dict[str, dict[str, np.ndarray]] = {}
        for index, (candidate, _, ce_nmse, llr, denominator, effective_noise) in enumerate(candidate_state):
            start = index * int(batch_size)
            stop = start + int(batch_size)
            decoded_payload = decoded.payload_bits[start:stop]
            crc_status = decoded.crc_status[start:stop]
            payload_error = np.any(decoded_payload != payload, axis=1)
            out[candidate.candidate_id] = {
                "error_flags": (payload_error | ~crc_status).astype(bool),
                "crc_status": crc_status.astype(bool),
                "payload_error": payload_error.astype(bool),
                "ce_nmse": ce_nmse.astype(np.float64),
                "llr": llr.astype(np.float64),
                "mrc_denominator": denominator.astype(np.float64),
                "mrc_denominator_trial_mean": np.mean(denominator, axis=1).astype(np.float64),
                "effective_noise_variance": effective_noise.astype(np.float64),
            }
        return out

    @staticmethod
    def _tag(value: float) -> str:
        return str(float(value)).replace("-", "m").replace(".", "p")

    @staticmethod
    def _safe_id(value: str) -> str:
        return "".join(char if char.isalnum() or char in "-_" else "_" for char in value)

    def _load_completed_rows(self, path: Path) -> list[dict[str, Any]]:
        if not path.exists() or not bool(self.config["simulation"].get("resume", True)):
            return []
        with open(path, "r", encoding="utf-8", newline="") as handle:
            return [dict(row) for row in csv.DictReader(handle)]

    def run(self) -> list[dict[str, Any]]:
        output = Path(self.config["output_dir"])
        output.mkdir(parents=True, exist_ok=True)
        flags_dir = output / "trial_error_flags"
        flags_dir.mkdir(parents=True, exist_ok=True)
        log_path = output / "run.log"
        log_lines = log_path.read_text(encoding="utf-8").splitlines() if log_path.exists() else []

        def log(message: str) -> None:
            line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {message}"
            print(line, flush=True)
            log_lines.append(line)
            log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")

        csv_path = output / "bler_points.csv"
        rows: list[dict[str, Any]] = self._load_completed_rows(csv_path)
        completed = {
            (str(row["candidate_id"]), float(row["snr_db"])) for row in rows
        }
        simulation = self.config["simulation"]
        minimum = int(simulation["min_trials_per_snr"])
        maximum = int(simulation["max_trials_per_snr"])
        target_errors = int(simulation["target_errors"])
        batch_size = int(self.config["batch_size"])
        progress_every_batches = int(simulation.get("progress_every_batches", 0))
        started = time.time()
        default_snrs = tuple(float(value) for value in simulation["snr_points_db"])
        snr_grid = sorted(
            {
                snr
                for candidate in self.candidates
                for snr in (candidate.snr_points_db or default_snrs)
            }
        )
        for snr_db in snr_grid:
            active_candidates = [
                candidate
                for candidate in self.candidates
                if snr_db in (candidate.snr_points_db or default_snrs)
            ]
            if all((candidate.candidate_id, snr_db) in completed for candidate in active_candidates):
                log(f"SNR={snr_db:g}dB already complete; skipped by resume.")
                continue
            flags = {candidate.candidate_id: [] for candidate in active_candidates}
            nmse = {candidate.candidate_id: [] for candidate in active_candidates}
            mrc_denominator = {candidate.candidate_id: [] for candidate in active_candidates}
            errors = {candidate.candidate_id: 0 for candidate in active_candidates}
            trials = 0
            checkpoint_enabled = bool(simulation.get("save_batch_checkpoints", False))
            if checkpoint_enabled and bool(simulation.get("resume", True)):
                loaded_trials: set[int] = set()
                for candidate in active_candidates:
                    candidate_id = candidate.candidate_id
                    stem = f"{self._safe_id(candidate_id)}_snr_{self._tag(snr_db)}"
                    checkpoint = flags_dir / f"{stem}_partial.npz"
                    if not checkpoint.exists():
                        loaded_trials.add(0)
                        continue
                    with np.load(checkpoint) as saved:
                        saved_flags = np.asarray(saved["error_flags"], dtype=bool)
                        saved_nmse = np.asarray(saved["ce_nmse"], dtype=np.float64)
                        saved_mrc = np.asarray(saved["mrc_denominator_trial_mean"], dtype=np.float64)
                    if not (len(saved_flags) == len(saved_nmse) == len(saved_mrc)):
                        raise RuntimeError(f"{checkpoint}: inconsistent partial checkpoint lengths.")
                    flags[candidate_id].append(saved_flags)
                    nmse[candidate_id].append(saved_nmse)
                    mrc_denominator[candidate_id].append(saved_mrc)
                    errors[candidate_id] = int(np.sum(saved_flags))
                    loaded_trials.add(len(saved_flags))
                if len(loaded_trials) != 1:
                    raise RuntimeError(
                        f"SNR={snr_db:g}: candidate partial checkpoints are not aligned: "
                        f"{sorted(loaded_trials)}"
                    )
                trials = loaded_trials.pop()
                if trials > maximum:
                    raise RuntimeError(
                        f"SNR={snr_db:g}: partial checkpoint has {trials} trials above maximum {maximum}."
                    )
            while trials < maximum and not (
                trials >= minimum and all(value >= target_errors for value in errors.values())
            ):
                current = min(batch_size, maximum - trials)
                result = self.run_batch(snr_db, trials + 1, current)
                for candidate in active_candidates:
                    candidate_id = candidate.candidate_id
                    flags[candidate_id].append(result[candidate_id]["error_flags"])
                    nmse[candidate_id].append(result[candidate_id]["ce_nmse"])
                    mrc_denominator[candidate_id].append(
                        result[candidate_id]["mrc_denominator_trial_mean"]
                    )
                    errors[candidate_id] += int(np.sum(result[candidate_id]["error_flags"]))
                trials += current
                if checkpoint_enabled:
                    for candidate in active_candidates:
                        candidate_id = candidate.candidate_id
                        stem = f"{self._safe_id(candidate_id)}_snr_{self._tag(snr_db)}"
                        destination = flags_dir / f"{stem}_partial.npz"
                        temporary = flags_dir / f"{stem}_partial.tmp.npz"
                        np.savez_compressed(
                            temporary,
                            error_flags=np.concatenate(flags[candidate_id]).astype(bool),
                            ce_nmse=np.concatenate(nmse[candidate_id]).astype(np.float64),
                            mrc_denominator_trial_mean=np.concatenate(
                                mrc_denominator[candidate_id]
                            ).astype(np.float64),
                        )
                        temporary.replace(destination)
                batch_index = int(math.ceil(trials / batch_size))
                if progress_every_batches > 0 and batch_index % progress_every_batches == 0:
                    log(
                        f"SNR={snr_db:g}dB progress={trials}/{maximum} "
                        + " ".join(
                            f"{candidate.candidate_id}:{errors[candidate.candidate_id]}"
                            for candidate in active_candidates
                        )
                    )
                if trials >= minimum and all(value >= target_errors for value in errors.values()):
                    break
            new_rows: list[dict[str, Any]] = []
            for candidate in active_candidates:
                candidate_id = candidate.candidate_id
                all_flags = np.concatenate(flags[candidate_id]).astype(bool)
                all_nmse = np.concatenate(nmse[candidate_id]).astype(np.float64)
                all_mrc = np.concatenate(mrc_denominator[candidate_id]).astype(np.float64)
                low, high = wilson_interval(errors[candidate_id], trials)
                stopping_reason = (
                    "minimum_trials_and_target_errors"
                    if trials >= minimum and errors[candidate_id] >= target_errors
                    else "maximum_trials"
                )
                row = {
                    "candidate_id": candidate_id,
                    "family": candidate.family,
                    "scheme": candidate.scheme,
                    "receiver_covariance_mode": candidate.receiver_covariance_mode,
                    "n_tx": int(self.config["antenna"]["n_tx"]),
                    "n_rx": int(self.config["antenna"]["n_rx"]),
                    "csi_mode": self.channel_estimation,
                    "snr_definition": (
                        "unit_total_tx_power_over_single_rx_branch_noise_power"
                    ),
                    "total_transmit_power_normalization": "unit_total_power_per_re",
                    "snr_db": snr_db,
                    "trials": trials,
                    "errors": errors[candidate_id],
                    "bler": float(errors[candidate_id]) / float(trials),
                    "bler_wilson95_lo": low,
                    "bler_wilson95_hi": high,
                    "ce_nmse_db": float(10.0 * np.log10(max(float(np.mean(all_nmse)), 1e-30))),
                    "mrc_denominator_mean": float(np.mean(all_mrc)),
                    "stopping_reason": stopping_reason,
                    "aggregation_level": int(self.resource_config.aggregation_level),
                    "occupied_rb": int(self.k_active // 12),
                    "k_active": self.k_active,
                    "data_re": int(len(self.data_local)),
                    "dmrs_re": int(len(self.pilot_local)),
                    "coded_bits": int(self.codec.coded_bits),
                    "pilot_rank": int(candidate.pilot_rank),
                    "pilot_condition_number": float(candidate.pilot_condition_number),
                    "ce_floor_nmse_db": float(
                        10.0 * np.log10(max(candidate.ce_floor_nmse, 1e-30))
                    ),
                }
                new_rows.append(row)
                if bool(simulation.get("save_trial_error_flags", True)):
                    stem = f"{self._safe_id(candidate_id)}_snr_{self._tag(snr_db)}"
                    np.save(flags_dir / f"{stem}_error_flags.npy", all_flags)
                    np.save(flags_dir / f"{stem}_ce_nmse.npy", all_nmse)
                    np.save(flags_dir / f"{stem}_mrc_denominator_trial_mean.npy", all_mrc)
                    partial = flags_dir / f"{stem}_partial.npz"
                    if partial.exists():
                        partial.unlink()
            rows.extend(new_rows)
            save_csv(rows, csv_path)
            log(
                f"SNR={snr_db:g}dB trials={trials} "
                + " ".join(
                    f"{candidate.candidate_id}:{errors[candidate.candidate_id]}/{trials}"
                    for candidate in active_candidates
                )
            )
        metadata = {
            "schema": "pdcch-cdd-bler-result-v2",
            "elapsed_seconds_this_invocation": float(time.time() - started),
            "resource_config": asdict(self.resource_config),
            "candidate_cces": self.pdcch_grid.candidate_cces.tolist(),
            "candidate_regs": self.pdcch_grid.candidate_regs.tolist(),
            "occupied_subcarriers": self.occupied.tolist(),
            "sample_mode": str(self.config.get("sample_mode", "shared_candidates")),
            "sample_scope": (
                "independent derived seed for this candidate"
                if str(self.config.get("sample_mode")) == "independent_candidate_stream"
                else "same payload, TDL realization, data noise, and DMRS noise per AL/SNR/absolute trial"
            ),
            "snr_definition": "Es/N0 per active unit-energy PDCCH RE at unit total transmit power",
            "system": {
                "n_tx": int(self.config["antenna"]["n_tx"]),
                "n_rx": int(self.config["antenna"]["n_rx"]),
                "layers": 1,
                "total_transmit_power": TOTAL_TRANSMIT_POWER,
                "rx_correlation_model": RX_CORRELATION_MODEL,
                "channel_estimation_across_rx": "independent_per_rx_branch",
                "combining": "coherent_mrc_across_rx",
            },
            "channel_estimation": self.channel_estimation,
            "candidates": [
                {
                    "candidate_id": candidate.candidate_id,
                    "family": candidate.family,
                    "scheme": candidate.scheme,
                    "receiver_covariance_mode": candidate.receiver_covariance_mode,
                    "label": candidate.label,
                    "style": candidate.style,
                    "precoder": candidate.precoder.metadata,
                    "pilot_rank": candidate.pilot_rank,
                    "pilot_condition_number": candidate.pilot_condition_number,
                    "ce_floor_nmse": candidate.ce_floor_nmse,
                }
                for candidate in self.candidates
            ],
            "polar": {
                "payload_bits": int(self.codec.payload_bits),
                "crc_bits": 24,
                "coded_bits": int(self.codec.coded_bits),
                "mother_code_length": int(self.codec.encoder.n_polar),
                "list_size": int(self.codec.list_size),
                "decoder_type": self.codec.decoder_type,
                "crc_rnti": int(self.codec.crc_rnti),
            },
            "rows": rows,
        }
        save_json(metadata, output / "run_metadata.json")
        self._plot(rows, output / "bler_curve.png")
        return rows

    def _plot(self, rows: list[dict[str, Any]], path: Path) -> None:
        import matplotlib.pyplot as plt

        fig, axis = plt.subplots(figsize=(10.2, 6.8))
        for candidate in self.candidates:
            points = sorted(
                [row for row in rows if str(row["candidate_id"]) == candidate.candidate_id],
                key=lambda row: float(row["snr_db"]),
            )
            snr = np.asarray([float(row["snr_db"]) for row in points])
            bler = np.asarray([float(row["bler"]) for row in points])
            trials = np.asarray([int(row["trials"]) for row in points])
            display = np.where(bler > 0, bler, 0.5 / np.maximum(trials, 1))
            kwargs = {
                "color": candidate.style.get("color"),
                "linestyle": candidate.style.get("linestyle", "-"),
                "marker": candidate.style.get("marker", "o"),
            }
            axis.semilogy(
                snr,
                display,
                label=candidate.label,
                linewidth=2.5,
                markersize=8,
                **{key: value for key, value in kwargs.items() if value is not None},
            )
        axis.axhline(0.10, color="black", linestyle=":", linewidth=1.8, label="BLER=0.10")
        axis.axhline(0.01, color="black", linestyle="--", linewidth=1.8, label="BLER=0.01")
        axis.set_xlabel("SNR (dB)", fontsize=16)
        axis.set_ylabel("Estimated-CSI DCI BLER", fontsize=16)
        axis.set_title(f"PDCCH AL{self.resource_config.aggregation_level} CDD comparison", fontsize=17)
        axis.grid(True, which="both", alpha=0.35)
        axis.tick_params(labelsize=14)
        axis.legend(fontsize=12, ncol=2, loc="best")
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
