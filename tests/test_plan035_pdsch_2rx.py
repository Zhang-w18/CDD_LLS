from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import numpy as np

from cdd_lls.phy.precoding import build_aged_mrt_prg_precoder
from tools.run_plan033_tdl_mobility_mimo import (
    _adaptive_schedule,
    _interval_identity,
    _make_error_flags,
    candidate_definitions,
)
from tools.run_plan035_pdsch_2rx import load_config, validate_config
from tools.run_plan035_cdl_pdsch_2rx import load_config as load_cdl_config
from tools.run_plan035_cdl_pdsch_2rx import validate_config as validate_cdl_config


ROOT = Path(__file__).resolve().parents[1]


def test_plan035_configs_and_twelve_record_identities() -> None:
    for n_tx in (4, 8):
        config = load_config(ROOT / "configs" / f"plan035_nt{n_tx}_nr2_v60_smoke.yaml")
        candidates, *_ = validate_config(config)
        records = [_interval_identity(config, row["candidate_id"], 4.0, 1, 20, receiver) for row in candidates for receiver in config["receiver_modes"]]
        assert len(records) == 12
        assert {row["receiver"] for row in records} == {"estimated", "ideal"}
        assert len({row["trial_key"] for row in records}) == 1


def test_two_rx_aged_mrt_uses_both_receive_branches() -> None:
    old = np.zeros((1, 2, 4, 72), dtype=np.complex128)
    old[0, 0, 0, :] = 1.0
    old[0, 1, 1, :] = 2.0
    weights = build_aged_mrt_prg_precoder(old, prg_size_rb=6, normalize=True)
    np.testing.assert_allclose(np.abs(weights[0, 0]), [0.0, 1.0, 0.0, 0.0], atol=1e-12)
    np.testing.assert_allclose(np.sum(np.abs(weights) ** 2, axis=2), 1.0, atol=2e-12)


def test_candidate_sets_match_frozen_plan() -> None:
    assert [row["delay_grid_coordinates"] for row in candidate_definitions(4) if "delay_grid_coordinates" in row] == [[0.0, 24.0, 48.0, 72.0], [0.0, 1.0, 3.0, 7.0], [0.0, 0.25, 0.5, 0.75], [0.0, 0.25, 0.5, 0.75]]


def test_cdl_d_smoke_config_builds_spatially_correlated_scene() -> None:
    path = ROOT / "configs" / "plan035_cdl_d_nt4_nr2_v60_smoke.yaml"
    config = load_cdl_config(path)
    _, channel, _, _, covariance = validate_cdl_config(config)
    assert channel.backend == "sionna_cdl"
    assert channel.cdl_profile == "D"
    assert (channel.cdl_tx_array_rows, channel.cdl_tx_array_cols) == (1, 4)
    assert (channel.cdl_rx_array_rows, channel.cdl_rx_array_cols) == (1, 2)
    assert covariance.covariance_type == "cdl_pdp_doppler_spatial_unaware"


def test_error_flags_are_keyed_by_candidate_and_receiver() -> None:
    flags = _make_error_flags(["candidate"], ["estimated", "ideal"], 3)
    flags[("candidate", "estimated")][0] = True
    assert flags[("candidate", "estimated")].tolist() == [True, False, False]
    assert flags[("candidate", "ideal")].tolist() == [False, False, False]
    assert "candidate" not in flags


def test_adaptive_schedule_only_appends_unmet_bracket_endpoint() -> None:
    rows = []
    for receiver in ("estimated", "ideal"):
        for snr_db, errors in ((8.0, 500), (10.0, 100), (12.0, 0)):
            rows.append({"candidate_id": "C", "receiver": receiver, "snr_db": snr_db, "trials": 1000, "tb_errors": errors})
    config = {
        "output_dir": "unused",
        "snr_db": [8.0, 10.0, 12.0],
        "receiver_modes": ["estimated", "ideal"],
        "adaptive_policy": {"minimum_trials": 1000, "check_interval_trials": 1000, "maximum_trials": 50000, "target_errors_10pct": 200, "target_errors_1pct": 200, "ten_pct_error_band": [0.05, 0.2], "one_pct_error_band": [0.005, 0.02]},
    }
    with patch("tools.run_plan033_tdl_mobility_mimo.read_csv_rows", return_value=rows):
        with patch("tools.run_plan033_tdl_mobility_mimo.Path.exists", return_value=True):
            targets, audit = _adaptive_schedule(config, [{"candidate_id": "C"}])
    assert audit["status"] == "append_required"
    assert targets == {8.0: 1000, 10.0: 2000, 12.0: 1000}
