from __future__ import annotations

import copy

import numpy as np
import pytest

from cdd_lls.core.config import ChannelConfig
from cdd_lls.phy.channel_tdl import (
    generate_sionna_tdl_channel,
    generate_sionna_tdl_channel_active,
)
from cdd_lls.phy.pdcch import (
    PDCCHResourceConfig,
    build_pdcch_grid,
    build_reg_bundle_dft_precoder,
)
from cdd_lls.phy.pdcch_codec import SionnaPDCCHPolarCodec, dci_mother_code_length
from cdd_lls.phy.precoding import spatial_dft_codebook
from cdd_lls.phy.estimators import build_frequency_rmmse_filter
from cdd_lls.sim.pdcch import (
    DEFAULT_PDCCH_BLER_CONFIG,
    PDCCHBLERSimulator,
    validate_pdcch_bler_config,
)
from cdd_lls.sim.pdcch_cdd import (
    PDCCHCDDBLERSimulator,
    _build_cdd_precoder,
    validate_pdcch_cdd_config,
)
from tools.run_pdcch_bler_curves import _candidate_stream_seed


def _plan031_config(aggregation_level: int = 2) -> dict:
    coded_bits = 108 * int(aggregation_level)
    k_active = 72 * int(aggregation_level)
    return {
        "schema": "pdcch-cdd-bler-v2",
        "output_dir": ".pytest_tmp_pdcch031",
        "seed": 20260903,
        "batch_size": 1,
        "antenna": {"n_tx": 8, "n_rx": 1},
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
            "aggregation_level": aggregation_level,
            "first_cce": 0,
        },
        "pdcch": {
            "payload_bits": 41,
            "coded_bits": coded_bits,
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
            "tdl_profile": "A",
            "delay_spread_ns": 100.0,
            "carrier_frequency_hz": 3.5e9,
            "ue_speed_kmh": 0.0,
            "num_sinusoids": 20,
            "normalize": True,
        },
        "receiver": {"channel_estimation": "frequency_lmmse", "llr_clip": 50.0},
        "simulation": {
            "snr_points_db": [20.0],
            "min_trials_per_snr": 1,
            "target_errors": 1,
            "max_trials_per_snr": 1,
            "save_trial_error_flags": True,
        },
        "candidates": [
            {
                "candidate_id": "CDD",
                "family": "TEST",
                "scheme": "cdd",
                "receiver_covariance_mode": "matched_effective",
                "delay_grid_coordinates": [
                    float(index) * float(k_active) / 64.0 for index in range(8)
                ],
            },
            {
                "candidate_id": "PRG",
                "family": "TEST",
                "scheme": "reg_bundle_dft_cycling",
                "receiver_covariance_mode": "physical_prg",
                "cycling_order": list(range(aggregation_level)),
            },
        ],
    }


def _plan031_c300_config(aggregation_level: int = 1) -> dict:
    config = _plan031_config(2)
    k_active = 36 * int(aggregation_level)
    config["seed"] = 20260907
    config["antenna"] = {"n_tx": 4, "n_rx": 1}
    config["resource"]["duration_symbols"] = 2
    config["resource"]["aggregation_level"] = int(aggregation_level)
    config["pdcch"]["coded_bits"] = 108 * int(aggregation_level)
    config["channel"].update(
        {
            "tdl_profile": "C",
            "delay_spread_ns": 300.0,
            "carrier_frequency_hz": 4.0e9,
            "ue_speed_kmh": 3.0,
            "num_sinusoids": 20,
        }
    )
    config["candidates"] = [
        {
            "candidate_id": "C300_S0_SIDON",
            "family": "S0_SIDON",
            "scheme": "cdd",
            "receiver_covariance_mode": "matched_effective",
            "delay_grid_coordinates": [0.0, 1.0, 3.0, 7.0],
        },
        {
            "candidate_id": "C300_SMALL_TRANSPARENT",
            "family": "SMALL_CDD",
            "scheme": "cdd",
            "receiver_covariance_mode": "physical_fullband",
            "delay_grid_coordinates": [
                0.0,
                k_active / 2304.0,
                2.0 * k_active / 2304.0,
                3.0 * k_active / 2304.0,
            ],
        },
        {
            "candidate_id": "C300_PRG_DFT4_6REG",
            "family": "PRG_DFT",
            "scheme": "reg_bundle_dft_cycling",
            "receiver_covariance_mode": "physical_prg",
            "cycling_order": list(range(aggregation_level)),
        },
    ]
    return config


def test_al4_resource_counts_and_noninterleaved_mapping() -> None:
    config = PDCCHResourceConfig(
        n_rb=48,
        duration_symbols=2,
        cce_reg_mapping="noninterleaved",
        reg_bundle_size=6,
        aggregation_level=4,
        first_cce=0,
    )
    grid = build_pdcch_grid(config)
    assert grid.n_cce == 16
    assert grid.candidate_cces.tolist() == [0, 1, 2, 3]
    assert grid.candidate_bundles.tolist() == [0, 1, 2, 3]
    assert grid.candidate_regs.tolist() == list(range(24))
    assert grid.resource_grid.n_data_re == 216
    assert grid.resource_grid.n_dmrs_re == 72
    assert grid.coded_bits == 432
    assert set((grid.dmrs_local_subcarrier_indices % 12).tolist()) == {1, 5, 9}
    assert set((grid.data_local_subcarrier_indices % 12).tolist()) == {
        0,
        2,
        3,
        4,
        6,
        7,
        8,
        10,
        11,
    }


def test_interleaved_mapping_and_bundle_constant_dft_precoder() -> None:
    config = PDCCHResourceConfig(
        n_rb=48,
        duration_symbols=2,
        cce_reg_mapping="interleaved",
        reg_bundle_size=2,
        interleaver_size=3,
        shift_index=1,
        aggregation_level=4,
        first_cce=0,
    )
    grid = build_pdcch_grid(config)
    assert len(grid.candidate_regs) == 24
    assert len(np.unique(grid.candidate_regs)) == 24
    assert len(grid.candidate_bundles) == 12

    precoder = build_reg_bundle_dft_precoder(grid, n_tx=8, cycling_order=range(8))
    assert precoder.C.shape == (576, 8)
    assert np.allclose(np.sum(np.abs(precoder.C) ** 2, axis=1), 1.0)
    for rb in range(48):
        assert np.allclose(precoder.C[rb * 12 : (rb + 1) * 12], precoder.C[rb * 12])
    assert not np.allclose(precoder.C[0], precoder.C[12])


def test_al8_resource_counts_and_full_dft8_bundle_cycle() -> None:
    config = PDCCHResourceConfig(
        n_rb=48,
        duration_symbols=2,
        cce_reg_mapping="noninterleaved",
        reg_bundle_size=6,
        aggregation_level=8,
        first_cce=0,
    )
    grid = build_pdcch_grid(config)
    assert grid.candidate_cces.tolist() == list(range(8))
    assert grid.candidate_bundles.tolist() == list(range(8))
    assert grid.candidate_regs.tolist() == list(range(48))
    assert grid.resource_grid.n_data_re == 432
    assert grid.resource_grid.n_dmrs_re == 144
    assert grid.coded_bits == 864

    precoder = build_reg_bundle_dft_precoder(grid, n_tx=8, cycling_order=range(8))
    used = [precoder.metadata["vector_by_physical_bundle"][bundle] for bundle in range(8)]
    assert used == list(range(8))
    for bundle in range(8):
        first_rb = 3 * bundle
        reference = precoder.C[first_rb * 12]
        assert np.allclose(
            precoder.C[first_rb * 12 : (first_rb + 3) * 12], reference
        )


@pytest.mark.parametrize(
    ("aggregation_level", "occupied_rb", "k_active", "n_p", "data_re", "dmrs_re"),
    [(2, 12, 144, 36, 108, 36), (4, 24, 288, 72, 216, 72), (8, 48, 576, 144, 432, 144)],
)
def test_plan031_resource_and_candidate_validation(
    aggregation_level: int,
    occupied_rb: int,
    k_active: int,
    n_p: int,
    data_re: int,
    dmrs_re: int,
) -> None:
    summary = validate_pdcch_cdd_config(_plan031_config(aggregation_level))
    assert summary["occupied_rb"] == occupied_rb
    assert summary["k_active"] == k_active
    assert summary["n_p"] == n_p
    assert summary["n_data_re"] == data_re
    assert summary["n_dmrs_re"] == dmrs_re


def test_pdcch_cdd_accepts_two_independent_rx_with_unit_total_power() -> None:
    config = _plan031_config(2)
    config["antenna"]["n_rx"] = 2
    summary = validate_pdcch_cdd_config(config)
    assert summary["n_rx"] == 2
    assert summary["layers"] == 1
    assert summary["total_transmit_power"] == 1.0
    assert summary["rx_correlation_model"] == "identity_no_spatial_correlation"

    simulator = PDCCHCDDBLERSimulator(config)
    for candidate in simulator.candidates:
        row_power = np.sum(np.abs(candidate.precoder.C) ** 2, axis=1)
        assert np.allclose(row_power, 1.0)
        assert candidate.precoder.metadata["total_transmit_power_per_re"] == 1.0

    estimator = simulator._filter(simulator.candidates[0], noise_variance=0.1)
    observations = np.zeros((3, 2, len(simulator.pilot_band)), dtype=np.complex128)
    observations[:, 0, :] = 1.0 + 0.5j
    estimate = estimator.estimate_full_band(observations)
    assert estimate.shape == (3, 2, simulator.k_active)
    assert np.any(np.abs(estimate[:, 0, :]) > 0.0)
    assert np.array_equal(estimate[:, 1, :], np.zeros_like(estimate[:, 1, :]))


def test_pdcch_cdd_two_rx_estimated_csi_batch_runs_through_mrc() -> None:
    config = _plan031_config(2)
    config["antenna"]["n_rx"] = 2
    config["candidates"] = config["candidates"][:1]
    simulator = PDCCHCDDBLERSimulator(config)

    state = simulator.run_batch(snr_db=30.0, absolute_start=1, batch_size=1)["CDD"]
    assert state["error_flags"].shape == (1,)
    assert state["ce_nmse"].shape == (1,)
    assert state["llr"].shape == (1, config["pdcch"]["coded_bits"])
    assert state["mrc_denominator"].shape == (
        1,
        simulator.validation["n_data_re"],
    )
    assert state["effective_noise_variance"].shape == state["mrc_denominator"].shape
    assert np.all(np.isfinite(state["ce_nmse"]))
    assert np.all(state["mrc_denominator"] > 0.0)


@pytest.mark.parametrize("n_rx", [0, -1, True, 1.5, "2"])
def test_pdcch_cdd_rejects_non_positive_or_non_integer_n_rx(n_rx: object) -> None:
    config = _plan031_config(2)
    config["antenna"]["n_rx"] = n_rx
    with pytest.raises(ValueError, match="antenna.n_rx must be a positive integer"):
        validate_pdcch_cdd_config(config)


@pytest.mark.parametrize(
    ("aggregation_level", "occupied_rb", "k_active", "n_p", "data_re", "dmrs_re"),
    [(1, 3, 36, 9, 54, 18), (2, 6, 72, 18, 108, 36), (4, 12, 144, 36, 216, 72)],
)
def test_plan031_c300_two_symbol_resource_counts(
    aggregation_level: int,
    occupied_rb: int,
    k_active: int,
    n_p: int,
    data_re: int,
    dmrs_re: int,
) -> None:
    summary = validate_pdcch_cdd_config(_plan031_c300_config(aggregation_level))
    assert summary["occupied_rb"] == occupied_rb
    assert summary["k_active"] == k_active
    assert summary["n_p"] == n_p
    assert summary["n_data_re"] == data_re
    assert summary["n_dmrs_re"] == dmrs_re
    assert summary["coded_bits"] == 108 * aggregation_level


def test_plan031_cdd_duplicate_delays_require_explicit_opt_in() -> None:
    config = _plan031_c300_config(1)
    candidate = config["candidates"][0]
    candidate["delay_grid_coordinates"] = [0.0, 0.0, 0.98388, 0.98388]

    with pytest.raises(ValueError, match="folded CDD residues are not distinct"):
        validate_pdcch_cdd_config(config)

    candidate["allow_duplicate_delays"] = True
    summary = validate_pdcch_cdd_config(config)
    diagnostics = summary["candidate_constraint_diagnostics"]["C300_S0_SIDON"]
    assert diagnostics["allow_duplicate_delays"] is True
    assert diagnostics["fold_residue_count"] == 2


@pytest.mark.parametrize("invalid", [1, "true", None])
def test_plan031_rejects_non_boolean_duplicate_delay_opt_in(invalid: object) -> None:
    config = _plan031_c300_config(1)
    config["candidates"][0]["allow_duplicate_delays"] = invalid
    with pytest.raises(ValueError, match="allow_duplicate_delays must be boolean"):
        validate_pdcch_cdd_config(config)


def test_plan031_c300_precoders_and_two_dimensional_filters() -> None:
    simulator = PDCCHCDDBLERSimulator(_plan031_c300_config(2))
    assert simulator.band_grid.n_symbols == 2
    assert simulator.band_grid.n_sc == 72
    assert simulator.candidates[0].precoder.C.shape == (72, 4)
    assert np.allclose(
        np.sum(np.abs(simulator.candidates[0].precoder.C) ** 2, axis=1), 1.0
    )
    matched_filter = simulator._filter(simulator.candidates[0], 0.1)
    transparent_filter = simulator._filter(simulator.candidates[1], 0.1)
    cycling_filter = simulator._filter(simulator.candidates[2], 0.1)
    assert matched_filter.weights.shape == (108, 36)
    assert transparent_filter.weights.shape == (108, 36)
    assert cycling_filter.weights.shape == (108, 36)
    assert len(cycling_filter.subfilter_diagnostics) == 2
    pilot_local = simulator.pilot_band
    data_local = simulator.data_band
    cross_bundle = (data_local[:, None] // 36) != (pilot_local[None, :] // 36)
    assert np.allclose(cycling_filter.weights[cross_bundle], 0.0)
    assert not np.allclose(matched_filter.weights, transparent_filter.weights)


def test_plan031_al1_dft_and_lte_codebook_vectors() -> None:
    codebook = spatial_dft_codebook(4, normalize=True)
    assert np.allclose(codebook[:, 0], np.asarray([1, 1, 1, 1]) / 2.0)
    assert np.allclose(codebook[:, 2], np.asarray([1, -1, 1, -1]) / 2.0)


def test_plan031_c300_channel_varies_across_two_symbols_and_batch_runs() -> None:
    config = _plan031_c300_config(1)
    simulator = PDCCHCDDBLERSimulator(config)
    channel = generate_sionna_tdl_channel_active(
        simulator.band_grid,
        ChannelConfig(**config["channel"]),
        n_tx=4,
        n_rx=1,
        batch_size=1,
        seed=17,
    )
    assert channel.H.shape == (1, 1, 4, 2, 36)
    assert not np.array_equal(channel.H[:, :, :, 0, :], channel.H[:, :, :, 1, :])
    result = simulator.run_batch(snr_db=30.0, absolute_start=1, batch_size=1)
    assert set(result) == {
        "C300_S0_SIDON",
        "C300_SMALL_TRANSPARENT",
        "C300_PRG_DFT4_6REG",
    }
    assert all(np.isfinite(value["ce_nmse"][0]) for value in result.values())


@pytest.mark.parametrize("aggregation_level", [1, 2, 4])
def test_plan031_c300_ideal_csi_is_exact_and_skips_lmmse(
    aggregation_level: int,
) -> None:
    config = _plan031_c300_config(aggregation_level)
    config["receiver"]["channel_estimation"] = "ideal"
    simulator = PDCCHCDDBLERSimulator(config)
    assert simulator._filters == {}
    assert all(candidate.ce_floor_nmse == 0.0 for candidate in simulator.candidates)
    result = simulator.run_batch(snr_db=float("inf"), absolute_start=1, batch_size=2)
    assert simulator._filters == {}
    assert all(np.array_equal(value["ce_nmse"], np.zeros(2)) for value in result.values())
    assert all(not np.any(value["error_flags"]) for value in result.values())


def test_plan031_c300_ideal_csi_mrc_denominator_and_noise_variance() -> None:
    received = np.asarray([[[1.0 + 2.0j, 3.0 - 1.0j], [2.0, -1.0j]]])
    estimate = np.asarray([[[1.0j, 2.0], [1.0, 1.0 - 1.0j]]])
    equalized, effective_noise, denominator = PDCCHCDDBLERSimulator._coherent_mrc(
        received, estimate, 0.25
    )
    expected_denominator = np.sum(np.abs(estimate) ** 2, axis=1)
    expected_equalized = np.sum(np.conj(estimate) * received, axis=1) / expected_denominator
    assert np.allclose(denominator, expected_denominator)
    assert np.allclose(equalized, expected_equalized)
    assert np.allclose(effective_noise, 0.25 / expected_denominator)


def test_plan031_c300_small_cdd_receiver_modes_are_identical_with_ideal_csi() -> None:
    config = _plan031_c300_config(2)
    config["receiver"]["channel_estimation"] = "ideal"
    matched = copy.deepcopy(config["candidates"][1])
    matched["candidate_id"] = "SMALL_MATCHED"
    matched["receiver_covariance_mode"] = "matched_effective"
    transparent = copy.deepcopy(matched)
    transparent["candidate_id"] = "SMALL_TRANSPARENT"
    transparent["receiver_covariance_mode"] = "physical_fullband"
    config["candidates"] = [matched, transparent]
    simulator = PDCCHCDDBLERSimulator(config)
    result = simulator.run_batch(snr_db=5.0, absolute_start=11, batch_size=2)
    assert np.array_equal(result["SMALL_MATCHED"]["llr"], result["SMALL_TRANSPARENT"]["llr"])
    assert np.array_equal(
        result["SMALL_MATCHED"]["error_flags"],
        result["SMALL_TRANSPARENT"]["error_flags"],
    )


def test_plan031_candidate_stream_seed_is_independent_by_stage_al_and_candidate() -> None:
    base = _plan031_c300_config(1)
    base["seed"] = 20260910
    base["random_stream_namespace"] = "prescan"
    values = {_candidate_stream_seed(base, "A")}
    formal = copy.deepcopy(base)
    formal["random_stream_namespace"] = "formal"
    values.add(_candidate_stream_seed(formal, "A"))
    al2 = copy.deepcopy(base)
    al2["resource"]["aggregation_level"] = 2
    values.add(_candidate_stream_seed(al2, "A"))
    values.add(_candidate_stream_seed(base, "B"))
    assert len(values) == 4


def test_plan031_cdd_precoder_uses_candidate_band_denominator() -> None:
    config = _plan031_config(2)
    grid = build_pdcch_grid(PDCCHResourceConfig(**config["resource"]))
    precoder = _build_cdd_precoder(grid, [0, 1, 3, 7, 12, 20, 30, 65], 8)
    assert precoder.metadata["phase_denominator"] == 144
    assert precoder.metadata["q_seconds"] == pytest.approx(1.0 / (144 * 30e3))
    assert np.allclose(np.sum(np.abs(precoder.C) ** 2, axis=1), 1.0)
    assert np.allclose(
        precoder.C[1, 1] / precoder.C[0, 1], np.exp(-2j * np.pi / 144)
    )


def test_plan031_matched_lmmse_candidates_run_high_snr_batch() -> None:
    simulator = PDCCHCDDBLERSimulator(_plan031_config(2))
    assert [candidate.pilot_rank for candidate in simulator.candidates] == [8, 2]
    result = simulator.run_batch(snr_db=30.0, absolute_start=1, batch_size=1)
    assert set(result) == {"CDD", "PRG"}
    assert all(value["error_flags"].shape == (1,) for value in result.values())
    assert all(np.isfinite(value["ce_nmse"][0]) for value in result.values())


def test_plan031_transparent_cdd_uses_only_physical_covariance() -> None:
    config = _plan031_config(2)
    transparent = copy.deepcopy(config["candidates"][0])
    transparent["candidate_id"] = "CDD_TRANSPARENT"
    transparent["receiver_covariance_mode"] = "physical_fullband"
    config["candidates"].insert(1, transparent)
    simulator = PDCCHCDDBLERSimulator(config)
    matched, transparent_runtime = simulator.candidates[:2]
    assert np.allclose(matched.covariance, transparent_runtime.covariance)
    matched_filter = simulator._filter(matched, 0.1)
    transparent_filter = simulator._filter(transparent_runtime, 0.1)
    assert not np.allclose(matched_filter.weights, transparent_filter.weights)
    expected = build_frequency_rmmse_filter(
        simulator.physical_covariance,
        simulator.pilot_band,
        noise_variance=0.1,
    )
    assert np.allclose(transparent_filter.weights, expected.weights)


def test_plan031_rejects_invalid_receiver_covariance_mode() -> None:
    config = _plan031_config(2)
    config["candidates"][0]["receiver_covariance_mode"] = "knows_cdd_but_claims_transparent"
    with pytest.raises(ValueError, match="receiver_covariance_mode"):
        validate_pdcch_cdd_config(config)


def test_active_only_sionna_tdl_matches_full_fft_slice() -> None:
    config = _plan031_config(2)
    grid = build_pdcch_grid(PDCCHResourceConfig(**config["resource"])).resource_grid
    channel = ChannelConfig(**config["channel"])
    full = generate_sionna_tdl_channel(grid, channel, n_tx=8, n_rx=1, batch_size=1, seed=91)
    active = generate_sionna_tdl_channel_active(
        grid, channel, n_tx=8, n_rx=1, batch_size=1, seed=91
    )
    assert active.H.shape == full.H.shape
    assert np.allclose(active.H, full.H, rtol=1e-12, atol=1e-12)


def test_pdcch_config_rejects_coded_length_mismatch() -> None:
    config = copy.deepcopy(DEFAULT_PDCCH_BLER_CONFIG)
    config["pdcch"]["coded_bits"] = 216
    with pytest.raises(ValueError, match="does not match"):
        validate_pdcch_bler_config(config)


def test_sionna_dci_affine_crc_and_noiseless_decode() -> None:
    from sionna.phy.fec.polar import PolarEncoder

    codec = SionnaPDCCHPolarCodec(
        payload_bits=41,
        coded_bits=432,
        crc_rnti=0xFFFF,
        scrambling_id=208,
        data_scrambling_rnti=0,
        list_size=8,
        decoder_type="hybSCL",
        cpu_only=True,
    )
    rng = np.random.default_rng(7)
    payload = rng.integers(0, 2, size=(2, 41), dtype=np.int8)
    exact = codec.encode_unscrambled(payload)
    transmitted = codec.encode(payload)
    initialized_crc = codec._crc_encode(
        np.concatenate([np.ones((2, 24), dtype=np.int8), payload], axis=1)
    )[:, -24:]
    standard_information = np.concatenate(
        [payload, initialized_crc ^ np.asarray([0] * 8 + [1] * 16, dtype=np.int8)],
        axis=1,
    )
    interleaved = standard_information[:, np.asarray(codec.encoder._ind_input_int)]
    mother = PolarEncoder(
        codec.encoder.frozen_pos,
        codec.encoder.n_polar,
        precision="double",
    )(codec.tf.constant(interleaved, dtype=codec.tf.float64)).numpy().astype(np.int8)
    manual_standard = mother[:, np.asarray(codec.encoder._ind_rate_matching)]
    assert exact.shape == (2, 432)
    assert np.array_equal(exact, manual_standard)
    assert np.array_equal(
        transmitted,
        exact ^ codec.scrambling_sequence[None, :],
    )
    llr = (2.0 * transmitted - 1.0) * 30.0
    decoded = codec.decode(llr)
    assert np.all(decoded.crc_status)
    assert np.array_equal(decoded.payload_bits, payload)


def test_al8_standard_repetition_and_noiseless_decode() -> None:
    codec = SionnaPDCCHPolarCodec(
        payload_bits=41,
        coded_bits=864,
        crc_rnti=0xFFFF,
        scrambling_id=208,
        data_scrambling_rnti=0,
        list_size=8,
        decoder_type="hybSCL",
        cpu_only=True,
    )
    assert dci_mother_code_length(41, 864) == 512
    assert codec.encoder.n_polar == 512
    assert codec.sionna_target_bits == 512
    assert codec.rate_matching_mode == "repetition"
    assert codec.external_repetition

    rng = np.random.default_rng(8)
    payload = rng.integers(0, 2, size=(2, 41), dtype=np.int8)
    exact = codec.encode_unscrambled(payload)
    assert exact.shape == (2, 864)
    assert np.array_equal(exact[:, 512:], exact[:, :352])

    probe = np.arange(864, dtype=np.float64)[None, :]
    recovered = codec._recover_external_repetition(probe)
    assert recovered.shape == (1, 512)
    assert np.array_equal(recovered[0, :352], probe[0, :352] + probe[0, 512:])
    assert np.array_equal(recovered[0, 352:], probe[0, 352:512])

    transmitted = codec.encode(payload)
    decoded = codec.decode((2.0 * transmitted - 1.0) * 30.0)
    assert np.all(decoded.crc_status)
    assert np.array_equal(decoded.payload_bits, payload)


def test_pdcch_sionna_tdl_one_batch_runs() -> None:
    config = copy.deepcopy(DEFAULT_PDCCH_BLER_CONFIG)
    config["batch_size"] = 1
    config["receiver"]["channel_estimation"] = "ideal"
    simulator = PDCCHBLERSimulator(config)
    result = simulator.run_batch(snr_db=30.0, absolute_start=1, batch_size=1)
    assert result["error_flags"].shape == (1,)
    assert not bool(result["error_flags"][0])
    assert float(result["ce_nmse"][0]) == pytest.approx(0.0)


def test_pdcch_al8_sionna_tdl_one_batch_runs() -> None:
    config = copy.deepcopy(DEFAULT_PDCCH_BLER_CONFIG)
    config["batch_size"] = 1
    config["resource"]["aggregation_level"] = 8
    config["pdcch"]["coded_bits"] = 864
    config["receiver"]["channel_estimation"] = "ideal"
    simulator = PDCCHBLERSimulator(config)
    result = simulator.run_batch(snr_db=30.0, absolute_start=1, batch_size=1)
    assert result["error_flags"].shape == (1,)
    assert not bool(result["error_flags"][0])
    assert float(result["ce_nmse"][0]) == pytest.approx(0.0)
