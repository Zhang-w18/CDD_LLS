from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from cdd_lls.phy.cdl_beam_platform import (
    CoverageRectangle,
    beam_domain_cdd_precoder,
    covariance_beam_joint,
    covariance_beam_specific_independent,
    covariance_common_reference_pdp,
    dft_steering_beams,
    maximum_normalized_pattern_capture,
    minimum_power_rectangle,
    partition_angular_region,
    regular_dft_codebook,
    sha256_json,
    split_rectangle,
    validate_regular_txru_mapping,
    validate_cache,
    write_cache,
)


def test_adaptive_dft_shapes_norms_phase_and_polarization() -> None:
    for nv, nh, pol, ov, oh in ((1, 4, 1, 1, 2), (2, 3, 2, 2, 1), (3, 2, 2, 1, 1)):
        codebook, metadata = regular_dft_codebook(nv, nh, pol, ov, oh)
        assert codebook.shape == (nv * nh * pol, ov * nv * oh * nh)
        assert len(metadata) == codebook.shape[1]
        np.testing.assert_allclose(np.linalg.norm(codebook, axis=0), 1.0, atol=1e-14)
        if pol == 2:
            np.testing.assert_allclose(codebook[:nv*nh], codebook[nv*nh:], atol=0.0)
        qv, qh = 1 % (ov*nv), 1 % (oh*nh)
        beam = qv*(oh*nh) + qh
        expected = np.exp(-1j*2*np.pi*(qv*1/(ov*nv) + qh*1/(oh*nh))) / np.sqrt(nv*nh*pol)
        if nv > 1 and nh > 1:
            assert codebook[nh + 1, beam] == pytest.approx(expected)


def test_irregular_mapping_is_rejected() -> None:
    mapping = np.eye(8, 4, dtype=np.complex128)
    with pytest.raises(ValueError, match="REGULAR_DFT_UNSUPPORTED_GEOMETRY"):
        validate_regular_txru_mapping(mapping, 2, 2, 1, 2, 2)


def test_minimum_rectangle_tie_break_and_splits() -> None:
    x = [-1, -1, 1, 1, 0]
    y = [-1, 1, -1, 1, 0]
    p = [.2, .2, .2, .2, .2]
    rectangle = minimum_power_rectangle(x, y, p, .8)
    assert rectangle.coverage >= .8
    assert rectangle.horizontal_min == -1
    for count in (1, 4, 8):
        regions = split_rectangle(rectangle, count, "auto", 4.0)
        assert len(regions) == count
        assert sum(region.area for region in regions) == pytest.approx(rectangle.area)


def test_two_stage_ssb_grid_and_narrow_dft_centers() -> None:
    ssb = partition_angular_region([-60.0, 60.0], [90.0, 110.0], 2, 4)
    assert len(ssb) == 8
    assert (ssb[0].aod_min_deg, ssb[0].aod_max_deg, ssb[0].zod_min_deg, ssb[0].zod_max_deg) == (
        -60.0, -30.0, 90.0, 100.0
    )
    assert (ssb[-1].aod_min_deg, ssb[-1].aod_max_deg, ssb[-1].zod_min_deg, ssb[-1].zod_max_deg) == (
        30.0, 60.0, 100.0, 110.0
    )
    strongest = ssb[5]
    narrow = partition_angular_region(
        [strongest.aod_min_deg, strongest.aod_max_deg],
        [strongest.zod_min_deg, strongest.zod_max_deg], 2, 4,
    )
    assert len(narrow) == 8
    assert narrow[0].aod_max_deg - narrow[0].aod_min_deg == pytest.approx(7.5)
    assert narrow[0].zod_max_deg - narrow[0].zod_min_deg == pytest.approx(5.0)
    response = np.asarray([[1, 1j], [1, -1j]], dtype=np.complex128)
    beams = dft_steering_beams(narrow[:2], response)
    np.testing.assert_allclose(np.linalg.norm(beams, axis=0), 1.0, atol=1e-14)
    assert np.argmax(np.abs(response @ beams), axis=0).tolist() == [0, 1]


def test_general_cdd_normalization_and_direct_equivalence() -> None:
    rng = np.random.default_rng(37)
    raw = rng.normal(size=(6, 4)) + 1j*rng.normal(size=(6, 4))
    q, _ = np.linalg.qr(raw)
    delay = [0, 1, 3, 7]
    k = np.arange(24)
    for weights, expected_orthogonal in ((q, True), (np.column_stack((q[:, :3], (q[:, 0]+q[:, 1])/np.sqrt(2))), False)):
        precoder, alpha, orthogonal = beam_domain_cdd_precoder(weights, delay, k, 24)
        assert orthogonal is expected_orthogonal
        np.testing.assert_allclose(np.sum(np.abs(precoder)**2, axis=1), 1.0, atol=1e-12)
        h = rng.normal(size=(2, 24, 6)) + 1j*rng.normal(size=(2, 24, 6))
        direct = np.einsum("rkt,kt->rk", h, precoder)
        branches = np.einsum("rkt,tm->rkm", h, weights)
        phase = np.exp(-1j*2*np.pi*k[:, None]*np.asarray(delay)[None, :]/24)
        reconstructed = np.einsum("rkm,km->rk", branches, alpha[:, None]*phase)
        np.testing.assert_allclose(direct, reconstructed, atol=1e-12)


def test_maximum_normalized_pattern_capture_removes_absolute_gain_scale() -> None:
    weight = np.asarray([1.0, 0.0], dtype=np.complex128)
    ray_response = np.asarray([[2.0, 0.0], [1.0, 3.0]], dtype=np.complex128)
    grid_response = np.asarray([[2.0, 0.0], [0.0, 4.0]], dtype=np.complex128)
    ratio, peak = maximum_normalized_pattern_capture(weight, ray_response, [.75, .25], grid_response)
    assert peak == pytest.approx(4.0)
    assert ratio == pytest.approx(.75 * 1.0 + .25 * .25)
    scaled_ratio, scaled_peak = maximum_normalized_pattern_capture(
        weight, 7.0 * ray_response, [.75, .25], 7.0 * grid_response
    )
    assert scaled_ratio == pytest.approx(ratio)
    assert scaled_peak == pytest.approx(49.0 * peak)
    ratio_with_off_grid_peak, peak_with_off_grid_peak = maximum_normalized_pattern_capture(
        weight, 2.0 * ray_response, [.75, .25], grid_response
    )
    assert peak_with_off_grid_peak == pytest.approx(16.0)
    assert ratio_with_off_grid_peak == pytest.approx(.75 * 1.0 + .25 * .25)


def test_three_covariances_match_manual_small_examples() -> None:
    f = np.asarray([-1.0, 0.0, 1.0])
    tau = np.asarray([0.0, .125])
    delay = np.asarray([0.0, .25])
    common = covariance_common_reference_pdp(f, tau, [.7, .3], delay)
    specific_power = np.asarray([[.7, .3], [.2, .8]])
    specific = covariance_beam_specific_independent(f, tau, specific_power, delay)
    joint = np.asarray([np.diag(specific_power[:, path]) for path in range(2)])
    ideal = covariance_beam_joint(f, tau, joint, delay)
    np.testing.assert_allclose(specific, ideal, atol=1e-12)
    for matrix in (common, specific, ideal):
        np.testing.assert_allclose(matrix, matrix.conj().T, atol=1e-12)
        assert np.min(np.linalg.eigvalsh(matrix)) >= -1e-12


def test_cache_round_trip_and_corruption_rejection(tmp_path: Path) -> None:
    regular, _ = regular_dft_codebook(1, 2, 1)
    arrays = {"regular_weights": regular, "selected_regular_weights": regular,
              "ssb_weights": regular, "wide_weight": regular[:, 0], "split_weights": regular}
    key = "a" * 64
    write_cache(tmp_path, {"cache_key": key}, arrays)
    assert validate_cache(tmp_path, key) is not None
    with np.load(tmp_path / "weights.npz") as loaded:
        broken = {name: np.asarray(loaded[name]) for name in loaded.files}
    broken["wide_weight"] = np.zeros_like(broken["wide_weight"])
    np.savez_compressed(tmp_path / "weights.npz", **broken)
    assert validate_cache(tmp_path, key) is None


def test_cache_key_changes_with_two_stage_angular_grid() -> None:
    baseline = {"algorithm_version": "v5", "ssb_grid": {"aod_range_deg": [-60.0, 60.0],
                                                           "zod_range_deg": [90.0, 110.0],
                                                           "vertical_beams": 2, "horizontal_beams": 4}}
    changed = {**baseline, "ssb_grid": {**baseline["ssb_grid"], "aod_range_deg": [-55.0, 60.0]}}
    assert sha256_json(baseline) != sha256_json(changed)
