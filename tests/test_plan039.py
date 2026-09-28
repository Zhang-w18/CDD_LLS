from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from cdd_lls.core.config import PlatformConfig, config_from_dict, dataclass_to_dict, load_config
from cdd_lls.phy.cdl_beam_platform import beam_domain_cdd_precoder
from cdd_lls.phy.plan039 import (
    CE_METHODS,
    build_angular_full_coverage_codebook,
    frequency_covariance,
    search_strict_sidon_top8,
    sidon_candidate_audit,
)
from cdd_lls.sim.orchestrator import CDDLinkLevelOrchestrator
from tools.run_plan039_cdl_beam_bler import (
    _adaptive_points,
    _aged_mrt_grid,
    _completed_method1_sidon_intervals,
    _completed_batch_covers,
    _evaluation_variants,
    _stage1b_grid,
)


ROOT = Path(__file__).resolve().parents[1]


def _nonorthogonal_beams() -> np.ndarray:
    rng = np.random.default_rng(39)
    beams = rng.normal(size=(12, 8)) + 1j * rng.normal(size=(12, 8))
    return beams / np.linalg.norm(beams, axis=0, keepdims=True)


def test_nonorthogonal_cdd_uses_exact_per_subcarrier_normalization() -> None:
    precoder, alpha, orthogonal = beam_domain_cdd_precoder(
        _nonorthogonal_beams(), [0, 1, 3, 7, 12, 20, 30, 65], np.arange(576), 576
    )
    assert not orthogonal
    assert np.ptp(alpha) > 1e-5
    np.testing.assert_allclose(np.sum(np.abs(precoder) ** 2, axis=1), 1.0, atol=1e-12)


def test_two_declared_cdd_covariances_are_hermitian_psd_and_alpha_aware() -> None:
    rng = np.random.default_rng(3902)
    factors = rng.normal(size=(3, 8, 4)) + 1j * rng.normal(size=(3, 8, 4))
    joint = np.einsum("lmi,lni->lmn", factors, factors.conj())
    powers = np.real(np.diagonal(joint, axis1=1, axis2=2)).T
    arrays = {
        "path_delays_s": np.asarray([0.0, 50e-9, 120e-9]),
        "reference_path_powers": np.asarray([0.7, 0.2, 0.1]),
        "beam_path_powers": powers,
        "joint_path_covariances": joint,
    }
    frequencies = np.arange(18) * 30e3
    delays = [0, 1, 3, 7, 12, 20, 30, 65]
    alpha = np.linspace(0.2, 0.5, len(frequencies))
    assert "PLAN039_BEAM_JOINT_COVARIANCE" not in CE_METHODS
    for method in ("PLAN039_COMMON_REFERENCE_PDP", "PLAN039_BEAM_SPECIFIC_PDP_INDEPENDENT"):
        covariance = frequency_covariance(method, frequencies, arrays, delays, alpha)
        np.testing.assert_allclose(covariance, covariance.conj().T, atol=1e-10)
        assert np.min(np.linalg.eigvalsh(covariance)) >= -1e-8
        changed = frequency_covariance(method, frequencies, arrays, delays, np.full(len(frequencies), 0.3))
        assert not np.allclose(covariance, changed)


def test_plan039_covariances_match_direct_576_denominator_formula() -> None:
    frequencies = np.arange(5, dtype=float) * 30e3
    physical_delays = np.asarray([0.0, 85e-9])
    reference_powers = np.asarray([0.7, 0.2])
    beam_powers = np.asarray([[0.1 + 0.01*m, 0.02 + 0.005*m] for m in range(8)])
    delay_indices = np.asarray([0, 9, 18, 27, 36, 45, 54, 63])
    arrays = {"path_delays_s": physical_delays, "reference_path_powers": reference_powers,
              "beam_path_powers": beam_powers}
    alpha = np.full(frequencies.size, 1.0 / np.sqrt(8.0))
    delta_k = np.arange(frequencies.size)[:, None] - np.arange(frequencies.size)[None, :]
    physical = np.exp(-1j * 2*np.pi * delta_k[:, :, None] * 30e3 * physical_delays[None, None, :])
    artificial = np.exp(-1j * 2*np.pi * delta_k[:, :, None] * delay_indices[None, None, :] / 576.0)
    expected_common = np.sum(physical * reference_powers[None, None, :], axis=2)
    expected_common *= np.mean(artificial, axis=2)
    expected_specific = np.zeros_like(expected_common)
    for branch in range(8):
        per_branch = np.sum(physical * beam_powers[branch][None, None, :], axis=2)
        expected_specific += artificial[:, :, branch] * per_branch / 8.0
    actual_common = frequency_covariance(
        "PLAN039_COMMON_REFERENCE_PDP", frequencies, arrays, delay_indices, alpha)
    actual_specific = frequency_covariance(
        "PLAN039_BEAM_SPECIFIC_PDP_INDEPENDENT", frequencies, arrays, delay_indices, alpha)
    np.testing.assert_allclose(actual_common, expected_common, atol=1e-12, rtol=1e-12)
    np.testing.assert_allclose(actual_specific, expected_specific, atol=1e-12, rtol=1e-12)
    wrong_4096 = np.mean(np.exp(-1j * 2*np.pi * delta_k[:, :, None]
                                * delay_indices[None, None, :] / 4096.0), axis=2)
    assert not np.allclose(actual_common, np.sum(physical * reference_powers, axis=2) * wrong_4096)


def test_ce_only_run_spec_does_not_construct_ldpc(monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = PlatformConfig()
    cfg.simulation.ce_only = True
    cfg.simulation.snr_points_db = [0.0]
    cfg.channel.backend = "fixed_cdl_statistics"
    cfg.channel_estimation.ce_method = "PLAN039_COMMON_REFERENCE_PDP"
    runner = CDDLinkLevelOrchestrator(cfg)
    monkeypatch.setattr("cdd_lls.sim.orchestrator.SionnaLDPCAdapter",
                        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("LDPC constructed")))
    captured: dict[str, object] = {}
    def fake_run_snr(**kwargs):
        captured.update(kwargs)
        return {"tb_errors": 0, "n_trials": 1}
    monkeypatch.setattr(runner, "_run_snr", fake_run_snr)
    runner._run_spec(cfg)
    assert captured["adapter"] is None


def test_stage1b_adaptive_budget_requests_paired_bracket_endpoints() -> None:
    rows = []
    for variant in ("b0__ideal", "sidon_selected__ideal", "cycling__ideal"):
        rows.extend([
            {"variant_id": variant, "snr_db": 4.0, "bler": 0.15,
             "n_trials": 1000, "tb_errors": 150},
            {"variant_id": variant, "snr_db": 4.25, "bler": 0.08,
             "n_trials": 1000, "tb_errors": 80},
        ])
    variants = ("b0__ideal", "sidon_selected__ideal", "cycling__ideal")
    assert _adaptive_points(rows, variants) == {1000: [4.0, 4.25]}


def test_stage1b_uses_family_specific_refinement_windows() -> None:
    raw = {"plan039": {"stage1b": {
        "snr_min_db": 7.5,
        "snr_max_db": 20.0,
        "coarse_step_db": 0.5,
        "fine_step_db": 0.25,
        "cdd_refinement_db": [10.0, 15.0],
        "cycling_refinement_db": [15.0, 20.0],
    }}}
    cdd = _stage1b_grid(raw, "cdd")
    cycling = _stage1b_grid(raw, "cycling")
    assert cdd[0] == cycling[0] == 7.5
    assert cdd[-1] == cycling[-1] == 20.0
    assert 10.25 in cdd and 10.25 not in cycling
    assert 19.75 in cycling and 19.75 not in cdd


def test_aged_mrt_formal_grid_and_config_are_frozen() -> None:
    raw = {"plan039": {"aged_mrt": {
        "snr_min_db": 7.5,
        "snr_max_db": 20.0,
        "snr_step_db": 0.5,
    }}}
    grid = _aged_mrt_grid(raw)
    assert len(grid) == 26
    assert grid[0] == 7.5 and grid[-1] == 20.0

    base = load_config(ROOT / "configs" / "plan039_cdl_wide_beam_bler.yaml")
    data = dataclass_to_dict(base)
    data["channel"]["ue_speed_kmh"] = 60.0
    data["transmission"].update({
        "tx_scheme": "PLAN039_AGED_MRT_PRG6",
        "aged_csi_ms": 40.0,
    })
    data["channel_estimation"]["ce_method"] = "PLAN039_PRG_COMMON_REFERENCE_PDP"
    cfg = config_from_dict(data)
    assert cfg.transmission.aged_csi_ms == 40.0
    data["transmission"]["aged_csi_ms"] = 0.0
    with pytest.raises(ValueError, match="aged_csi_ms"):
        config_from_dict(data)


def test_stage1b_runs_only_method1_and_prg_estimated_curves() -> None:
    freeze = {"selected": {
        "manifest": "/tmp/selected.json",
        "manifest_sha256": "a" * 64,
    }}
    variants = _evaluation_variants(freeze)
    assert len(variants) == 3
    methods = [row["channel_estimation"]["ce_method"] for row in variants]
    assert "IDEAL" not in methods
    assert "PLAN039_BEAM_SPECIFIC_PDP_INDEPENDENT" not in methods
    assert methods.count("PLAN039_COMMON_REFERENCE_PDP") == 2
    assert "PLAN039_TRANSPARENT_COMMON_REFERENCE_PDP" not in methods
    assert methods.count("PLAN039_PRG_COMMON_REFERENCE_PDP") == 1


def test_stage1b_can_reuse_completed_legacy_superset_batch(tmp_path) -> None:
    directory = tmp_path / "legacy"
    directory.mkdir()
    (directory / "summary.csv").write_text(
        "variant_id,snr_db,absolute_trial_start,n_trials\n"
        "b0__plan039_common_reference_pdp,10.0,0,1000\n"
        "b0__plan039_transparent_common_reference_pdp,10.0,0,1000\n",
        encoding="utf-8",
    )
    (directory / "trial_metrics.csv").write_text("variant_id\n", encoding="utf-8")
    (directory / "batch_receipt.json").write_text("{}\n", encoding="utf-8")
    assert _completed_batch_covers(
        directory, ("b0__plan039_common_reference_pdp",), [10.0], 0, 1000
    )


def test_transparent_sidon_copies_only_completed_method1_intervals(tmp_path) -> None:
    batches = tmp_path / "estimated_confirm" / "batches"
    completed = batches / "completed"
    completed.mkdir(parents=True)
    (completed / "summary.csv").write_text(
        "variant_id,snr_db,absolute_trial_start,absolute_trial_stop,n_trials\n"
        "sidon_selected__plan039_common_reference_pdp,13.75,0,1000,1000\n"
        "sidon_selected__plan039_common_reference_pdp,14.0,0,1000,1000\n"
        "b0__plan039_common_reference_pdp,13.75,0,1000,1000\n",
        encoding="utf-8",
    )
    (completed / "trial_metrics.csv").write_text("variant_id\n", encoding="utf-8")
    (completed / "batch_receipt.json").write_text("{}\n", encoding="utf-8")
    incomplete = batches / "incomplete"
    incomplete.mkdir()
    (incomplete / "summary.csv").write_text(
        "variant_id,snr_db,absolute_trial_start,absolute_trial_stop,n_trials\n"
        "sidon_selected__plan039_common_reference_pdp,16.0,1000,2000,1000\n",
        encoding="utf-8",
    )
    assert _completed_method1_sidon_intervals(tmp_path) == {
        (0, 1000): [13.75, 14.0],
    }


def test_prg_common_reference_covariance_omits_cdd_shift() -> None:
    arrays = {
        "path_delays_s": np.asarray([0.0, 100e-9]),
        "reference_path_powers": np.asarray([2.0, 1.0]),
        "beam_path_powers": np.ones((8, 2)),
        "joint_path_covariances": np.zeros((2, 8, 8), complex),
    }
    covariance = frequency_covariance(
        "PLAN039_PRG_COMMON_REFERENCE_PDP", np.arange(12) * 30e3, arrays, [], np.ones(12)
    )
    np.testing.assert_allclose(np.real(np.diag(covariance)), 3.0, atol=1e-12)


def test_transparent_common_reference_covariance_omits_cdd_shift_fullband() -> None:
    arrays = {
        "path_delays_s": np.asarray([0.0, 100e-9]),
        "reference_path_powers": np.asarray([2.0, 1.0]),
        "beam_path_powers": np.ones((8, 2)),
        "joint_path_covariances": np.zeros((2, 8, 8), complex),
    }
    frequencies = np.arange(12) * 30e3
    transparent = frequency_covariance(
        "PLAN039_TRANSPARENT_COMMON_REFERENCE_PDP",
        frequencies,
        arrays,
        [0, 9, 18, 27, 36, 45, 54, 63],
        np.full(12, 1.0 / np.sqrt(8.0)),
    )
    prg = frequency_covariance(
        "PLAN039_PRG_COMMON_REFERENCE_PDP", frequencies, arrays, [], np.ones(12)
    )
    np.testing.assert_allclose(transparent, prg, atol=1e-12, rtol=1e-12)


def test_shifted_dft_codebook_is_orthogonal_and_ordered() -> None:
    reference, narrow, audit = build_angular_full_coverage_codebook()
    assert reference.shape == (32,)
    assert narrow.shape == (32, 8)
    np.testing.assert_allclose(narrow.conj().T @ narrow, np.eye(8), atol=1e-12)
    np.testing.assert_allclose(np.linalg.norm(reference), 1.0, atol=1e-12)
    centers = np.asarray(audit["horizontal_center_aod_deg"])
    assert np.all(np.diff(centers) > 0.0)
    assert centers[0] >= -61.5 and centers[-1] <= 61.5


def test_sidon_search_is_rule_only_max_fold_and_stable() -> None:
    result = search_strict_sidon_top8()
    assert result["uses_channel_or_performance_data"] is False
    assert result["theoretical_fold_upper_bound"] == 12
    assert result["attempted_fold_gaps"][12]["search_complete"] == 1
    assert result["attempted_fold_gaps"][12]["strict_sidon_candidates"] == 0
    assert result["maximum_fold_gap"] == 11
    assert len(result["top8"]) == 8
    assert result["top8"][0]["delay_indices"] == [0, 11, 28, 148, 170, 233, 277, 351]
    for row in result["top8"]:
        assert row["strict_sidon"] is True
        assert row["fold_gap"] == 11
        assert row["pilot_rank"] == 8
        assert len(row["unordered_pair_sums"]) == 36
        assert len(set(row["unordered_pair_sums"])) == 36
    assert result["top8"] == search_strict_sidon_top8()["top8"]


def test_historical_s0_is_strict_but_not_maximum_fold_gap() -> None:
    audit = sidon_candidate_audit([0, 1, 3, 7, 12, 20, 30, 65])
    assert audit["strict_sidon"] is True
    assert audit["fold_gap"] < 12
