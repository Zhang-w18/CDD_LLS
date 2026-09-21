import numpy as np

from cdd_lls.phy.estimators import build_prg_frequency_rmmse_filter
from cdd_lls.phy.precoding import (
    build_prg_dft_precoder,
    build_prg_dft_precoder_batch,
    equivalent_channel,
    spatial_dft_codebook,
)
from cdd_lls.phy.resource_grid import build_resource_grid
from tools import run_bler_curves as runner
from tools import run_plan027_bler as bler027


def test_spatial_dft_codebook_is_orthogonal_with_power_eight() -> None:
    codebook = spatial_dft_codebook(8, normalize=False)
    np.testing.assert_allclose(codebook.conj().T @ codebook, 8.0 * np.eye(8), atol=1e-12)
    np.testing.assert_allclose(np.abs(codebook), 1.0, atol=1e-12)


def test_prg_boundaries_and_batched_precoder_orders() -> None:
    grid = build_resource_grid(bler027.resource_config(6))
    fixed = build_prg_dft_precoder(grid, 8, 6, list(range(8)), normalize=False)
    codebook = spatial_dft_codebook(8, normalize=False)
    for prg in range(8):
        start = 72 * prg
        np.testing.assert_allclose(
            fixed.C[start : start + 72], np.repeat(codebook[:, prg][None, :], 72, axis=0)
        )
    orders = np.asarray([list(range(8)), list(reversed(range(8)))], dtype=np.int64)
    batched = build_prg_dft_precoder_batch(grid, 8, 6, orders, normalize=False)
    assert batched.shape == (2, 576, 8)
    np.testing.assert_allclose(batched[0], fixed.C)
    np.testing.assert_allclose(
        batched[1, :72], np.repeat(codebook[:, 7][None, :], 72, axis=0)
    )


def test_batched_equivalent_channel_matches_trialwise_call() -> None:
    grid = build_resource_grid(bler027.resource_config(6))
    rng = np.random.default_rng(5)
    h = rng.normal(size=(2, 4, 8, 2, 576)) + 1j * rng.normal(size=(2, 4, 8, 2, 576))
    orders = np.asarray([list(range(8)), list(reversed(range(8)))], dtype=np.int64)
    precoders = build_prg_dft_precoder_batch(grid, 8, 6, orders, normalize=False)
    actual = equivalent_channel(h, precoders)
    assert actual.shape == (2, 4, 2, 576)
    for trial in range(2):
        expected = equivalent_channel(h[trial : trial + 1], precoders[trial])
        np.testing.assert_allclose(actual[trial : trial + 1], expected, atol=1e-12)


def test_random_tail_is_trial_reproducible_and_without_replacement() -> None:
    candidate = {
        "candidate_id": "A100_PRG_DFT8_4RB",
        "prg_size_rb": 4,
        "mapping": "cycle8_random_tail4",
    }
    orders = [runner.transparent_prg_order(candidate, 20260727, 14.0, trial) for trial in range(1, 20)]
    assert all(order[:8] == list(range(8)) for order in orders)
    assert all(len(set(order[8:])) == 4 for order in orders)
    assert len({tuple(order[8:]) for order in orders}) > 1
    assert orders[6] == runner.transparent_prg_order(candidate, 20260727, 14.0, 7)


def test_prg_lmmse_is_block_local_and_scaling_equivalent() -> None:
    n_sc = 24
    rho = 0.92
    covariance = rho ** np.abs(np.arange(n_sc)[:, None] - np.arange(n_sc)[None, :])
    pilots = np.asarray([0, 6, 12, 18], dtype=np.int64)
    scaled = build_prg_frequency_rmmse_filter(
        8.0 * covariance, pilots, 12, noise_variance=4.0, diagonal_loading=1e-10
    )
    normalized = build_prg_frequency_rmmse_filter(
        covariance, pilots, 12, noise_variance=0.5, diagonal_loading=1.25e-11
    )
    np.testing.assert_allclose(scaled.weights, normalized.weights, atol=1e-10, rtol=1e-10)
    np.testing.assert_allclose(scaled.weights[:12, 2:], 0.0, atol=0.0)
    np.testing.assert_allclose(scaled.weights[12:, :2], 0.0, atol=0.0)
    assert [row["pilot_count"] for row in scaled.prg_diagnostics] == [2, 2]


def test_two_dmrs_average_halves_complex_noise_variance() -> None:
    rng = np.random.default_rng(9)
    count = 200000
    one = (rng.normal(size=count) + 1j * rng.normal(size=count)) / np.sqrt(2.0)
    two = (rng.normal(size=count) + 1j * rng.normal(size=count)) / np.sqrt(2.0)
    averaged = 0.5 * (one + two)
    np.testing.assert_allclose(np.mean(np.abs(averaged) ** 2), 0.5, atol=0.01)
