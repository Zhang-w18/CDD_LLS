import unittest

import numpy as np

from cdd_lls.design.cdd_metrics import (
    FrequencyCEMetric,
    array_factor,
    build_cdd_kernel,
    jcdd_metrics,
)


class CDDMetricsTest(unittest.TestCase):
    def test_array_factor_matches_direct_sum(self):
        delays = [0, 1, 3]
        factor = array_factor(delays, period=16)
        direct = np.asarray([
            np.mean(np.exp(-2j * np.pi * lag * np.asarray(delays) / 16.0))
            for lag in range(16)
        ])
        np.testing.assert_allclose(factor, direct, atol=1e-14)

    def test_kernel_is_zero_at_independence_and_monotone(self):
        db_grid = np.arange(-20.0, 35.1, 0.1)
        values = 4.0 / (1.0 + np.exp(-(db_grid - 5.0) / 3.0))
        kernel = build_cdd_kernel(
            16.0,
            db_grid,
            values,
            rho_step=0.1,
            quadrature_order=5,
        )
        self.assertAlmostEqual(float(kernel.covariance_values[0]), 0.0)
        self.assertTrue(np.all(np.diff(kernel.covariance_values) >= -1e-12))

    def test_jcdd_and_ce_are_finite(self):
        period = 24
        lags = np.arange(period)
        first_row = np.exp(-0.05 * lags)
        covariance = np.empty((period, period), dtype=np.complex128)
        for first in range(period):
            for second in range(period):
                covariance[first, second] = np.exp(-0.05 * abs(first - second))
        db_grid = np.arange(-20.0, 35.1, 0.1)
        values = 4.0 / (1.0 + np.exp(-(db_grid - 5.0) / 3.0))
        kernel = build_cdd_kernel(
            16.0,
            db_grid,
            values,
            rho_step=0.1,
            quadrature_order=5,
        )
        metrics = jcdd_metrics([0, 1, 3, 7], covariance, kernel)
        self.assertGreaterEqual(metrics["jcdd"], 0.0)
        ce = FrequencyCEMetric(
            covariance,
            pilot_local_indices=[0, 4, 8, 12, 16, 20],
            data_local_indices=np.arange(period),
            n_tx=4,
        ).evaluate([0, 1, 2, 3], 16.0)
        self.assertTrue(np.isfinite(float(ce["nmse_db"])))
        self.assertGreaterEqual(float(ce["nmse"]), 0.0)


if __name__ == "__main__":
    unittest.main()
