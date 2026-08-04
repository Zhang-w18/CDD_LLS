from __future__ import annotations

import unittest

import numpy as np

from cdd_lls.core.config import ResourceConfig
from cdd_lls.phy.precoding import build_active_dft_grid_precoder, normalize_delay_vector
from cdd_lls.phy.resource_grid import build_resource_grid


class PrecodingTests(unittest.TestCase):
    def setUp(self):
        self.grid = build_resource_grid(ResourceConfig(n_prbs=48, n_fft=4096, scs_khz=30))

    def test_fractional_fft_sample_delays_are_preserved(self):
        values = normalize_delay_vector([0.0, 64.0 / 9.0], n_tx=2)
        self.assertEqual(values[0], 0)
        self.assertAlmostEqual(float(values[1]), 64.0 / 9.0, places=14)

    def test_active_dft_grid_precoder_matches_plan024(self):
        for indices in (
            [0, 9, 18, 27, 36, 45, 54, 63],
            [0, 1, 3, 7, 12, 20, 30, 65],
        ):
            result = build_active_dft_grid_precoder(self.grid, indices, n_tx=8, normalize=False)
            local = np.arange(576, dtype=np.float64)
            expected = np.exp(
                -2j * np.pi * local[:, None] * np.asarray(indices, dtype=np.float64)[None, :] / 576.0
            )
            np.testing.assert_allclose(result.C, expected, rtol=0.0, atol=1e-12)
            self.assertEqual(result.metadata["phase_reference"], "first_active_subcarrier")
            self.assertEqual(result.metadata["phase_denominator"], 576)
            self.assertFalse(result.metadata["normalized"])

    def test_active_grid_to_fft_sample_conversion(self):
        indices = [0, 1, 3, 7, 12, 20, 30, 65]
        result = build_active_dft_grid_precoder(self.grid, indices, n_tx=8, normalize=False)
        expected_samples = np.asarray(indices, dtype=np.float64) * 4096.0 / 576.0
        expected_seconds = np.asarray(indices, dtype=np.float64) / (576.0 * 30e3)
        np.testing.assert_allclose(result.metadata["cdd_delay_fft_samples"], expected_samples, atol=1e-12)
        np.testing.assert_allclose(result.metadata["cdd_delay_seconds"], expected_seconds, atol=1e-18)


if __name__ == "__main__":
    unittest.main()
