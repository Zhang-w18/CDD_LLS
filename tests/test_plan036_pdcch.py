from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from cdd_lls.sim.pdcch_cdd import PDCCHCDDBLERSimulator
from tools.run_plan036_pdcch import (
    ESTIMATED_IDS,
    IDEAL_IDS,
    _base_config,
    _candidate_specs,
    load_config,
    validate_config,
)


ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(
    ("name", "al", "k_active", "data_re", "dmrs_re", "coded_bits"),
    [
        ("pdcch_plan036_al1_fd5_formal.yaml", 1, 72, 54, 18, 108),
        ("pdcch_plan036_al2_fd1100_formal.yaml", 2, 144, 108, 36, 216),
    ],
)
def test_plan036_validation_and_resource_counts(
    name: str, al: int, k_active: int, data_re: int, dmrs_re: int, coded_bits: int
) -> None:
    config = load_config(ROOT / "configs" / name)
    receipt = validate_config(config)
    assert receipt["aggregation_level"] == al
    assert receipt["estimated"]["k_active"] == k_active
    assert receipt["estimated"]["n_data_re"] == data_re
    assert receipt["estimated"]["n_dmrs_re"] == dmrs_re
    assert receipt["estimated"]["coded_bits"] == coded_bits
    assert receipt["estimated"]["n_rx"] == 4
    assert abs(receipt["doppler_roundtrip_error_hz"]) < 1e-12


def test_plan036_al2_geometry_and_qc_transmit_identity() -> None:
    config = load_config(ROOT / "configs" / "pdcch_plan036_al2_fd5_formal.yaml")
    specs = _candidate_specs(config)
    sidon = next(item for item in specs if item["candidate_id"] == "SIDON")
    assert sidon["delay_grid_coordinates"] == [0, 11, 19, 64]
    assert [value % 36 for value in sidon["delay_grid_coordinates"]] == [0, 11, 19, 28]
    pair_sums = [
        (sidon["delay_grid_coordinates"][left] + sidon["delay_grid_coordinates"][right]) % 144
        for left in range(4)
        for right in range(left, 4)
    ]
    assert len(set(pair_sums)) == 10
    simulator = PDCCHCDDBLERSimulator(_base_config(config, "estimated"))
    runtimes = {item.candidate_id: item for item in simulator.candidates}
    assert tuple(runtimes) == ESTIMATED_IDS
    assert np.array_equal(
        runtimes["QC_DELAY_NT"].precoder.C,
        runtimes["QC_DELAY_TRANSPARENT"].precoder.C,
    )
    assert np.allclose(
        np.sum(np.abs(runtimes["SIDON"].precoder.C) ** 2, axis=1), 1.0
    )


def test_plan036_ideal_has_three_physical_waveforms_and_no_filters() -> None:
    config = load_config(ROOT / "configs" / "pdcch_plan036_al1_fd5_formal.yaml")
    simulator = PDCCHCDDBLERSimulator(_base_config(config, "ideal"))
    assert tuple(item.candidate_id for item in simulator.candidates) == IDEAL_IDS
    assert simulator._filters == {}
    result = simulator.run_batch(snr_db=float("inf"), absolute_start=1, batch_size=2)
    assert simulator._filters == {}
    assert all(np.array_equal(item["ce_nmse"], np.zeros(2)) for item in result.values())
    assert all(not np.any(item["error_flags"]) for item in result.values())


@pytest.mark.parametrize(
    "name",
    ["pdcch_plan036_al1_fd5_4t2r_formal.yaml", "pdcch_plan036_al2_fd5_4t2r_formal.yaml"],
)
def test_plan036_4t2r_configs_reuse_the_same_waveforms(name: str) -> None:
    config = load_config(ROOT / "configs" / name)
    receipt = validate_config(config)
    assert receipt["n_rx"] == 2
    assert receipt["estimated"]["n_rx"] == 2
    simulator = PDCCHCDDBLERSimulator(_base_config(config, "estimated"))
    assert simulator.config["antenna"] == {"n_tx": 4, "n_rx": 2}
    expected = _candidate_specs(config)
    assert [item.candidate_id for item in simulator.candidates] == [
        item["candidate_id"] for item in expected
    ]


def test_plan036_4t2r_al1_extension_reuses_output_and_adds_6_7_8_db() -> None:
    original = load_config(ROOT / "configs" / "pdcch_plan036_al1_fd5_4t2r_formal.yaml")
    extension = load_config(ROOT / "configs" / "pdcch_plan036_al1_fd5_4t2r_snr8_extension.yaml")
    assert extension["output_dir"] == original["output_dir"]
    assert extension["random_stream_namespace"] == original["random_stream_namespace"]
    assert extension["snr_points_db"][:-3] == original["snr_points_db"]
    assert extension["snr_points_db"][-4:] == [5.0, 6.0, 7.0, 8.0]
    assert validate_config(extension)["estimated"]["n_rx"] == 2
