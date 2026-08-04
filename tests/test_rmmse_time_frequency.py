from __future__ import annotations

import unittest

import numpy as np

from cdd_lls.core.config import ChannelConfig, ResourceConfig
from cdd_lls.phy.estimators import (
    build_frequency_rmmse_filter,
    build_time_frequency_rmmse_filter,
    linear_estimator_closed_form_nmse,
    tdl_known_delay_covariance,
    tdl_unknown_delay_covariance,
)
from cdd_lls.phy.resource_grid import build_resource_grid
from cdd_lls.sim.orchestrator import construct_ls_observations


class TimeFrequencyRMMSETests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.grid = build_resource_grid(ResourceConfig(n_prbs=8, dmrs_spacing_sc=24))
        cls.channel3 = ChannelConfig(
            backend="sionna_tdl",
            tdl_profile="A",
            delay_spread_ns=100.0,
            carrier_frequency_hz=3.5e9,
            ue_speed_kmh=3.0,
        )

    def test_covariance_structure_and_zero_cdd_identity(self):
        unknown = tdl_unknown_delay_covariance(self.grid, self.channel3)
        zero = tdl_known_delay_covariance(self.grid, self.channel3, [0, 0, 0, 0])
        np.testing.assert_allclose(zero.frequency, unknown.frequency, atol=1e-14, rtol=1e-14)
        for matrix in (zero.time, zero.frequency):
            np.testing.assert_allclose(matrix, matrix.conj().T, atol=1e-13, rtol=1e-13)
            self.assertGreaterEqual(float(np.min(np.linalg.eigvalsh(matrix))), -1e-10)
            np.testing.assert_allclose(np.diag(matrix), 1.0, atol=1e-13, rtol=1e-13)
        known_filter = build_time_frequency_rmmse_filter(self.grid, zero, 0.1, 1e-12)
        unknown_filter = build_time_frequency_rmmse_filter(self.grid, unknown, 0.1, 1e-12)
        np.testing.assert_allclose(known_filter.weights, unknown_filter.weights, atol=1e-13, rtol=1e-13)
        observations = (
            np.random.default_rng(21).normal(size=(2, self.grid.pilot_re_count))
            + 1j * np.random.default_rng(22).normal(size=(2, self.grid.pilot_re_count))
        )
        np.testing.assert_allclose(
            known_filter.estimate_data(observations),
            unknown_filter.estimate_data(observations),
            atol=1e-13,
            rtol=1e-13,
        )

    def test_high_speed_long_lag_correlation_is_lower(self):
        channel60 = ChannelConfig(**{
            **self.channel3.__dict__,
            "ue_speed_kmh": 60.0,
        })
        slow = tdl_unknown_delay_covariance(self.grid, self.channel3).time
        fast = tdl_unknown_delay_covariance(self.grid, channel60).time
        self.assertLess(abs(fast[0, -1]), abs(slow[0, -1]))

    def test_known_expected_nmse_not_worse_than_unknown(self):
        true = tdl_known_delay_covariance(self.grid, self.channel3, [0, 9, 18, 27])
        unknown = tdl_unknown_delay_covariance(self.grid, self.channel3)
        known_filter = build_time_frequency_rmmse_filter(self.grid, true, 0.1, 1e-12)
        unknown_filter = build_time_frequency_rmmse_filter(self.grid, unknown, 0.1, 1e-12)
        known_nmse = linear_estimator_closed_form_nmse(self.grid, known_filter, true, 0.1)
        unknown_nmse = linear_estimator_closed_form_nmse(self.grid, unknown_filter, true, 0.1)
        self.assertLessEqual(known_nmse, unknown_nmse + 1e-12)

    def test_each_dmrs_observation_has_unhalved_noise_variance(self):
        noise_variance = 0.2
        true = np.zeros((20000, 2, 4), dtype=np.complex128)
        observed = construct_ls_observations(true, noise_variance, np.random.default_rng(1234))
        empirical_by_symbol = np.mean(np.abs(observed) ** 2, axis=(0, 2))
        np.testing.assert_allclose(empirical_by_symbol, noise_variance, atol=0.006, rtol=0.0)

    def test_frequency_rmmse_filter_matches_direct_solve(self):
        covariance = np.asarray(
            [[1.0, 0.4 - 0.1j, 0.2], [0.4 + 0.1j, 1.0, 0.3], [0.2, 0.3, 1.0]],
            dtype=np.complex128,
        )
        pilots = np.asarray([0, 2], dtype=np.int64)
        noise = 0.25
        filt = build_frequency_rmmse_filter(covariance, pilots, noise)
        expected = covariance[:, pilots] @ np.linalg.inv(
            covariance[np.ix_(pilots, pilots)] + noise * np.eye(2)
        )
        np.testing.assert_allclose(filt.weights, expected, atol=1e-14, rtol=1e-14)
        observations = np.asarray([[1.0 + 0.5j, -0.25j]], dtype=np.complex128)
        np.testing.assert_allclose(
            filt.estimate_full_band(observations), observations @ expected.T, atol=1e-14, rtol=1e-14
        )


if __name__ == "__main__":
    unittest.main()
