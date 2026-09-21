from __future__ import annotations

from pathlib import Path

import numpy as np
import yaml

from cdd_lls.sim.pdcch_cdd import PDCCHCDDBLERSimulator, validate_pdcch_cdd_config
from cdd_lls.sim.rsrp import (
    noiseless_block_average_power,
    noiseless_block_average_power_per_rx,
)
from tools.run_plan034_rsrp_cdf import (
    PER_RX_CANDIDATES,
    build_per_rx_precoders,
    build_precoders,
    load_config,
)


ROOT = Path(__file__).resolve().parents[1]


def test_noiseless_block_average_power_matches_elementwise_formula() -> None:
    channel = np.asarray(
        [
            [
                [[[1 + 1j, 2 - 1j]], [[0.5 - 0.5j, -1 + 2j]]],
                [[[2 + 0j, -0.5 + 1j]], [[1 - 1j, 0.25 + 0.75j]]],
            ]
        ],
        dtype=np.complex128,
    )
    precoder = np.asarray([[1.0, 1.0j], [1.0, -1.0]], dtype=np.complex128) / np.sqrt(2.0)
    symbols = np.asarray([0, 0])
    subcarriers = np.asarray([0, 1])
    actual = noiseless_block_average_power(
        channel,
        precoder,
        symbols,
        subcarriers,
        n_rx=2,
        transmit_power=1.0,
    )
    expected_sum = 0.0
    for receive in range(2):
        for re_index, (symbol, subcarrier) in enumerate(zip(symbols, subcarriers)):
            effective = np.sum(
                channel[0, receive, :, symbol, subcarrier] * precoder[re_index]
            )
            expected_sum += float(np.abs(effective) ** 2)
    expected = expected_sum / (2 * 2)
    np.testing.assert_allclose(actual, [expected])

    per_rx = noiseless_block_average_power_per_rx(
        channel,
        precoder,
        symbols,
        subcarriers,
        n_rx=2,
        transmit_power=1.0,
    )
    assert per_rx.shape == (1, 2)
    np.testing.assert_allclose(np.mean(per_rx, axis=1), actual)


def test_plan034_precoders_and_rsrp_config_geometry() -> None:
    config = load_config(ROOT / "configs" / "pdcch_plan034_rsrp_cdf.yaml")
    assert config["simulation"]["trials"] == 100000
    precoders = build_precoders(36, 4, 36)
    fixed = precoders["FIXED_DFT0"]
    freq = precoders["FREQ_SIDON_0137"]
    np.testing.assert_allclose(fixed, np.ones((36, 4)) / 2.0)
    np.testing.assert_allclose(np.sum(np.abs(freq) ** 2, axis=1), 1.0)
    np.testing.assert_allclose(freq[1], np.exp(-2j * np.pi * np.asarray([0, 1, 3, 7]) / 36) / 2)


def test_plan034_per_rx_four_waveform_config_and_precoders() -> None:
    config = load_config(
        ROOT / "configs" / "pdcch_plan034_rsrp_per_rx_4waveform.yaml"
    )
    assert config["simulation"]["trials"] == 100000
    assert tuple(item["candidate_id"] for item in config["precoding"]["candidates"]) == (
        PER_RX_CANDIDATES
    )
    precoders = build_per_rx_precoders(36, 4, 36)
    assert tuple(precoders) == PER_RX_CANDIDATES
    for matrix in precoders.values():
        np.testing.assert_allclose(np.sum(np.abs(matrix) ** 2, axis=1), 1.0)
    np.testing.assert_allclose(precoders["FIXED_DFT0"], np.ones((36, 4)) / 2.0)
    np.testing.assert_allclose(
        precoders["CDD911"][1],
        np.exp(-2j * np.pi * np.asarray([0.0, 0.0, 0.98388, 0.98388]) / 36.0)
        / 2.0,
    )
    np.testing.assert_allclose(
        precoders["CDD130"][1],
        np.exp(-2j * np.pi * np.asarray([0.0, 0.0, 0.1404, 0.1404]) / 36.0)
        / 2.0,
    )


def test_plan034_a40_four_rx_config_and_manual_mrc() -> None:
    path = ROOT / "configs" / "pdcch_plan034_fixed_dft0_4rx_bler_mid.yaml"
    with open(path, "r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    summary = validate_pdcch_cdd_config(config)
    assert summary["n_rx"] == 4
    assert summary["n_data_re"] == 54
    assert summary["n_dmrs_re"] == 18
    rng = np.random.default_rng(20260916)
    received = rng.normal(size=(2, 4, 54)) + 1j * rng.normal(size=(2, 4, 54))
    estimate = rng.normal(size=(2, 4, 54)) + 1j * rng.normal(size=(2, 4, 54))
    equalized, effective_noise, denominator = PDCCHCDDBLERSimulator._coherent_mrc(
        received, estimate, 0.25
    )
    expected_denominator = np.sum(np.abs(estimate) ** 2, axis=1)
    np.testing.assert_allclose(denominator, expected_denominator)
    np.testing.assert_allclose(
        equalized, np.sum(np.conj(estimate) * received, axis=1) / expected_denominator
    )
    np.testing.assert_allclose(effective_noise, 0.25 / expected_denominator)
