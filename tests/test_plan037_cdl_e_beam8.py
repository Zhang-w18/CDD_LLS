from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from cdd_lls.core.config import dataclass_to_dict, config_from_dict, load_config
from cdd_lls.phy.beam8 import (
    DELAYS,
    build_beam8_precoder,
    build_dft_2x8_same_pol_codebook,
    equivalent_channel_two_paths,
    circular_pdp_moments,
    frequency_covariance,
    load_frozen_manifest,
    merge_discrete_pdp,
    pdp_moments,
    shifted_reference_pdp,
    write_manifest,
)
from cdd_lls.phy.resource_grid import build_resource_grid
from cdd_lls.phy.estimators import TDLTimeFrequencyCovariance, build_prg_time_frequency_rmmse_filter
from tools.analyze_plan037_cdl_e_beam8 import crossing_bracket


ROOT = Path(__file__).resolve().parents[1]


def test_dft_2x8_same_pol_definition_and_ids() -> None:
    codebook = build_dft_2x8_same_pol_codebook()
    assert codebook.shape == (32, 16)
    np.testing.assert_allclose(codebook.conj().T @ codebook, np.eye(16), atol=1e-13)
    np.testing.assert_allclose(codebook[:16], codebook[16:], atol=0.0)
    for beam in range(16):
        qv, qh = divmod(beam, 8)
        expected = np.asarray([
            np.exp(-1j * 2 * np.pi * (qv * m / 2 + qh * n / 8)) / np.sqrt(32)
            for m in range(2) for n in range(8)
        ])
        np.testing.assert_allclose(codebook[:16, beam], expected, atol=1e-14)


def test_all_beam8_precoders_are_unit_power_and_directly_equivalent() -> None:
    cfg = load_config(ROOT / "configs" / "plan037_cdl_e_32tx_2rx_beam8_statistics.yaml")
    grid = build_resource_grid(cfg.resource)
    codebook = build_dft_2x8_same_pol_codebook()
    selected = [0, 2, 4, 6, 9, 11, 13, 15]
    rng = np.random.default_rng(37)
    h = rng.normal(size=(2, 576, 32)) + 1j * rng.normal(size=(2, 576, 32))
    for scheme in ("BEAM8_B0_QC", "BEAM8_S0_SIDON", "BEAM8_PRECODER_CYCLING"):
        result = build_beam8_precoder(grid, scheme, selected, codebook, 6)
        np.testing.assert_allclose(np.sum(np.abs(result.C) ** 2, axis=1), 1.0, atol=1e-12)
        direct, branch_first = equivalent_channel_two_paths(h, result.C, codebook[:, selected], scheme)
        np.testing.assert_allclose(direct, branch_first, atol=1e-11)
    assert DELAYS["BEAM8_B0_QC"] == (0, 9, 18, 27, 36, 45, 54, 63)
    assert DELAYS["BEAM8_S0_SIDON"] == (0, 1, 3, 7, 12, 20, 30, 65)


def test_manifest_hash_and_rejection(tmp_path: Path) -> None:
    manifest = {
        "status": "FROZEN", "codebook_type": "dft_2x8_same_pol",
        "selected_beam_indices": list(range(8)), "mean_rsrp": [float(index + 1) for index in range(16)],
    }
    path = tmp_path / "manifest.json"
    digest = write_manifest(path, manifest)
    assert load_frozen_manifest(path, digest)["selected_beam_indices"] == list(range(8))
    tampered = json.loads(path.read_text(encoding="utf-8"))
    tampered["selected_beam_indices"] = list(range(1, 9))
    path.write_text(json.dumps(tampered), encoding="utf-8")
    with pytest.raises(ValueError, match="hash mismatch"):
        load_frozen_manifest(path, digest)


def test_link_config_requires_profile_native_angles_and_manifest() -> None:
    cfg = load_config(ROOT / "configs" / "plan037_cdl_e_32tx_2rx_beam8_statistics.yaml")
    data = dataclass_to_dict(cfg)
    data["fixed_cdl_statistics"]["selection_only"] = False
    with pytest.raises(ValueError, match="frozen manifest"):
        config_from_dict(data)
    data["fixed_cdl_statistics"].update({
        "selection_only": True, "profile_native_angles": False,
    })
    with pytest.raises(ValueError, match="profile_native_angles"):
        config_from_dict(data)


def test_zero_error_prescan_point_still_forms_a_bracket() -> None:
    assert crossing_bracket([(16.0, 0.6925), (18.0, 0.0)], 0.10) == (16.0, 18.0)
    assert crossing_bracket([(16.0, 0.6475), (18.0, 0.0025)], 0.01) == (16.0, 18.0)


def test_pdp_merge_moments_and_circular_moments() -> None:
    delays, powers = merge_discrete_pdp([0.0, 1e-9, 1e-9], [1.0, 1.0, 2.0])
    np.testing.assert_allclose(delays, [0.0, 1e-9])
    np.testing.assert_allclose(powers, [0.25, 0.75])
    mean, rms = pdp_moments(delays, powers)
    assert mean == pytest.approx(0.75e-9)
    assert rms == pytest.approx(np.sqrt(0.1875) * 1e-9)
    circular_mean, circular_rms = circular_pdp_moments([0.99, 0.01], [0.5, 0.5], 1.0)
    assert circular_mean == pytest.approx(0.0, abs=1e-12)
    assert circular_rms == pytest.approx(0.01)


def test_shifted_reference_pdp_matches_hadamard_frequency_covariance() -> None:
    delays = np.asarray([0.0, 120e-9])
    powers = np.asarray([0.7, 0.3])
    indices = [0, 1, 3, 7, 12, 20, 30, 65]
    n_sc = 576
    scs_hz = 30e3
    shifted_d, shifted_p = shifted_reference_pdp(delays, powers, indices, n_sc, scs_hz)
    direct = frequency_covariance(shifted_d, shifted_p, n_sc, scs_hz)
    reference = frequency_covariance(delays, powers, n_sc, scs_hz)
    delta = np.arange(n_sc)[:, None] - np.arange(n_sc)[None, :]
    factor = np.mean(np.exp(-1j * 2 * np.pi * delta[:, :, None] * np.asarray(indices) / n_sc), axis=2)
    np.testing.assert_allclose(direct, reference * factor, atol=1e-11)


def test_prg_2d_lmmse_has_no_cross_prg_weights() -> None:
    cfg = load_config(ROOT / "configs" / "plan037_cdl_e_32tx_2rx_beam8_smoke.yaml")
    grid = build_resource_grid(cfg.resource)
    covariance = TDLTimeFrequencyCovariance(np.eye(grid.n_symbols), np.eye(grid.n_sc), "unit_test")
    estimator = build_prg_time_frequency_rmmse_filter(grid, covariance, 72, 0.1, 1e-10)
    pilot_local = estimator.pilot_coordinates[:, 1] - grid.subcarrier_indices[0]
    data_local = estimator.data_coordinates[:, 1] - grid.subcarrier_indices[0]
    for data_row, local in enumerate(data_local):
        outside = pilot_local // 72 != local // 72
        assert np.all(estimator.weights[data_row, outside] == 0.0)


def test_beam8_receiver_mode_validation() -> None:
    cfg = load_config(ROOT / "configs" / "plan037_cdl_e_32tx_2rx_beam8_smoke.yaml")
    data = dataclass_to_dict(cfg)
    data["transmission"]["tx_scheme"] = "BEAM8_PRECODER_CYCLING"
    data["channel_estimation"]["ce_method"] = "BEAM8_CDD_AWARE_LMMSE"
    with pytest.raises(ValueError, match="cycling requires"):
        config_from_dict(data)
    data["transmission"]["tx_scheme"] = "BEAM8_B0_QC"
    data["channel_estimation"]["ce_method"] = "BEAM8_PRG_LMMSE"
    with pytest.raises(ValueError, match="CDD requires"):
        config_from_dict(data)
