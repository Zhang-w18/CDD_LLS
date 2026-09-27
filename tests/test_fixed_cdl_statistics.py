from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from cdd_lls.core.config import config_from_dict, dataclass_to_dict, load_config
from cdd_lls.phy.channel_cdl_fixed import (
    FixedCDLStatisticsChannel,
    build_contiguous_txru_mapping,
    build_fixed_cdl_codebooks,
    build_fixed_cdl_precoder,
)
from cdd_lls.phy.resource_grid import build_resource_grid


ROOT = Path(__file__).resolve().parents[1]


def test_32tx_mapping_and_codebooks() -> None:
    cfg = load_config(ROOT / "configs" / "smoke_fixed_cdl_statistics_32tx_2rx.yaml")
    fixed = cfg.fixed_cdl_statistics
    mapping = build_contiguous_txru_mapping(8, 8, 2, 8)
    assert mapping.shape == (64, 16)
    np.testing.assert_allclose(mapping.conj().T @ mapping, np.eye(16), atol=1e-14)
    parent, secondary = build_fixed_cdl_codebooks(fixed)
    assert parent.shape == (32, 8)
    assert secondary.shape == (32, 16)
    np.testing.assert_allclose(np.sum(np.abs(parent) ** 2, axis=0), 1.0, atol=1e-14)
    np.testing.assert_allclose(np.sum(np.abs(secondary) ** 2, axis=0), 1.0, atol=1e-14)


def test_all_four_precoders_have_unit_power_on_every_subcarrier() -> None:
    cfg = load_config(ROOT / "configs" / "smoke_fixed_cdl_statistics_32tx_2rx.yaml")
    grid = build_resource_grid(cfg.resource)
    parent, secondary = build_fixed_cdl_codebooks(cfg.fixed_cdl_statistics)
    for scheme in ("BASELINE", "POLARIZATION_CYCLING", "BEAM_CYCLING", "BEAM_CDD"):
        result = build_fixed_cdl_precoder(
            grid,
            cfg.fixed_cdl_statistics,
            scheme,
            selected_ssb=3,
            parent_codebook=parent,
            secondary_codebook=secondary,
            prg_size_rb=cfg.resource.prg_size_rb,
        )
        assert result.C.shape == (grid.n_sc, 32)
        np.testing.assert_allclose(np.sum(np.abs(result.C) ** 2, axis=1), 1.0, atol=1e-12)


def test_fixed_cdl_config_checks_dimensions_and_seed_separation() -> None:
    cfg = load_config(ROOT / "configs" / "smoke_fixed_cdl_statistics_32tx_4rx.yaml")
    assert cfg.antenna.n_tx == 32
    assert cfg.antenna.n_rx == 4
    invalid = dataclass_to_dict(cfg)
    invalid["fixed_cdl_statistics"]["realization_seed"] = invalid["fixed_cdl_statistics"]["statistics_seed"]
    with pytest.raises(ValueError, match="seeds must differ"):
        config_from_dict(invalid)
    invalid = dataclass_to_dict(cfg)
    invalid["fixed_cdl_statistics"]["ue_horizontal_elements"] = 1
    with pytest.raises(ValueError, match="UE port count"):
        config_from_dict(invalid)

    single_4rx = dataclass_to_dict(cfg)
    single_4rx["fixed_cdl_statistics"].update({
        "ue_vertical_elements": 2,
        "ue_horizontal_elements": 2,
        "ue_polarizations": 1,
    })
    assert config_from_dict(single_4rx).antenna.n_rx == 4

    dual_2rx = load_config(ROOT / "configs" / "smoke_fixed_cdl_statistics_32tx_2rx.yaml")
    single_2rx = dataclass_to_dict(dual_2rx)
    single_2rx["fixed_cdl_statistics"].update({
        "ue_vertical_elements": 1,
        "ue_horizontal_elements": 2,
        "ue_polarizations": 1,
    })
    assert config_from_dict(single_2rx).antenna.n_rx == 2


def test_fixed_cdl_replays_realization_but_resamples_small_scale_phase() -> None:
    base = load_config(ROOT / "configs" / "smoke_fixed_cdl_statistics_32tx_2rx.yaml")
    data = dataclass_to_dict(base)
    data["antenna"]["n_tx"] = 4
    data["fixed_cdl_statistics"].update({
        "covariance_realizations": 1,
        "bs_vertical_aes": 2,
        "bs_horizontal_aes": 2,
        "bs_vertical_txrus_per_pol": 1,
        "bs_horizontal_txrus_per_pol": 2,
        "ssb_horizontal_beams": 2,
        "secondary_horizontal_beams": 4,
    })
    cfg = config_from_dict(data)
    grid = build_resource_grid(cfg.resource)
    channel = FixedCDLStatisticsChannel(cfg, grid)
    first = channel.generate(0)
    replay = channel.generate(0)
    second = channel.generate(1)
    assert first.H.shape == (1, 2, 4, grid.n_symbols, grid.n_sc)
    np.testing.assert_array_equal(first.H, replay.H)
    assert not np.array_equal(first.H, second.H)
    assert channel.transmit_covariance.shape == (4, 4)
    assert channel.reference_receive_power > 0.0
    analytic = channel.analytic_transmit_covariance()
    assert analytic.shape == (4, 4)
    np.testing.assert_allclose(analytic, analytic.conj().T, atol=1e-12)
    assert np.min(np.linalg.eigvalsh(analytic)) >= -1e-10
