import unittest

import numpy as np

from cdd_lls.design import (
    FrequencyCEMetric,
    delay_ns_to_grid_coordinates,
    effective_moments_direct,
    effective_moments_lag,
    equivalent_cdd_delay_spread,
    phase_matrix_from_grid_coordinates,
)
from tools.run_plan027_meff_design import paired_bootstrap_gain


class Plan027EffectiveMomentsTest(unittest.TestCase):
    def test_direct_and_lag_implementations_match(self):
        period = 32
        correlation = 0.92
        indices = np.arange(period)
        covariance = correlation ** np.abs(indices[:, None] - indices[None, :])
        coordinates = [0.0, 1.25, 3.0, 7.5]
        direct = effective_moments_direct(coordinates, covariance)
        lag = effective_moments_lag(coordinates, covariance)
        self.assertAlmostEqual(direct["m2_eff"], lag["m2_eff"], places=10)
        self.assertAlmostEqual(direct["m4_eff"], lag["m4_eff"], places=10)

    def test_continuous_physical_delay_is_not_rounded(self):
        delay_ns = np.arange(8, dtype=np.float64) * 5.0
        coordinates = delay_ns_to_grid_coordinates(delay_ns)
        self.assertFalse(np.allclose(coordinates, np.round(coordinates)))
        phase = phase_matrix_from_grid_coordinates(coordinates)
        expected = np.exp(-2j * np.pi * 17 * 30e3 * delay_ns * 1e-9)
        np.testing.assert_allclose(phase[17], expected, atol=1e-14)

    def test_ce_supports_fractional_grid_coordinates(self):
        period = 24
        indices = np.arange(period)
        covariance = 0.9 ** np.abs(indices[:, None] - indices[None, :])
        metric = FrequencyCEMetric(
            covariance,
            pilot_local_indices=np.arange(0, period, 3),
            data_local_indices=np.arange(period),
            n_tx=4,
        )
        result = metric.evaluate_grid_coordinates([0.0, 0.2, 0.4, 0.6], 16.0)
        self.assertTrue(np.isfinite(float(result["nmse_db"])))
        self.assertEqual(int(result["pilot_rank"]), 4)

    def test_paired_bootstrap_replays_and_identical_samples_give_zero(self):
        grid = np.arange(0.0, 5.0, 0.5)
        thresholds = np.linspace(0.1, 4.9, 2000, dtype=np.float32)
        curve = np.asarray([np.mean(thresholds > value) for value in grid])
        first = paired_bootstrap_gain(
            thresholds,
            thresholds,
            curve,
            curve,
            grid,
            target=0.10,
            repeats=100,
            seed=7,
        )
        second = paired_bootstrap_gain(
            thresholds,
            thresholds,
            curve,
            curve,
            grid,
            target=0.10,
            repeats=100,
            seed=7,
        )
        np.testing.assert_array_equal(first, second)
        np.testing.assert_allclose(first, 0.0, atol=1e-14)

    def test_equivalent_delay_spread_respects_circular_cdd_period(self):
        result = equivalent_cdd_delay_spread(
            physical_delays_s=[0.0],
            physical_powers=[1.0],
            delay_grid_coordinates=8 * np.arange(8),
            period=64,
            subcarrier_spacing_hz=1.0,
            contained_probability=0.99,
        )
        self.assertAlmostEqual(
            result["equivalent_circular_rms_delay_spread_s"],
            np.sqrt(21.0) / 16.0,
            places=12,
        )
        self.assertAlmostEqual(
            result["equivalent_circular_support_width_s"],
            7.0 / 8.0,
            places=12,
        )


if __name__ == "__main__":
    unittest.main()
