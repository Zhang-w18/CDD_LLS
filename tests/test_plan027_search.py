import unittest

import numpy as np

from cdd_lls.design import (
    canonical_delay_set,
    continuous_arithmetic_delays_ns,
    delay_ns_to_grid_coordinates,
    thick_sidon_requirements,
    unique_random_delay_sets,
)
from tools.run_plan027_dense_dmrs import useful_symbol_nt_minus_one_delays_ns
from tools.run_plan027_meff_design import (
    AP_TALIAS_NT,
    AP_TU_NT,
    candidate_geometry,
)


class Plan027SearchTest(unittest.TestCase):
    def test_system_geometry_baseline_delays_and_pilot_ranks(self):
        self.assertEqual(
            AP_TU_NT,
            (0, 72, 144, 216, 288, 360, 432, 504),
        )
        self.assertEqual(
            AP_TALIAS_NT,
            (0, 3, 6, 9, 12, 15, 18, 21),
        )
        self.assertEqual(candidate_geometry(AP_TU_NT)["pilot_rank"], 1)
        self.assertEqual(candidate_geometry(AP_TALIAS_NT)["pilot_rank"], 8)
        self.assertEqual(candidate_geometry(AP_TU_NT, 12)["pilot_rank"], 2)
        self.assertEqual(candidate_geometry(AP_TU_NT, 6)["pilot_rank"], 4)
        self.assertEqual(
            candidate_geometry((0, 6, 12, 18, 24, 30, 36, 42), 12)[
                "pilot_rank"
            ],
            8,
        )
        self.assertEqual(
            candidate_geometry((0, 12, 24, 36, 48, 60, 72, 84), 6)[
                "pilot_rank"
            ],
            8,
        )

    def test_t1_arithmetic_interface_uses_only_spacing(self):
        self.assertEqual(
            continuous_arithmetic_delays_ns(5.0),
            (0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0),
        )

    def test_tu_nt_minus_one_baseline_preserves_period_endpoint(self):
        delays_ns = useful_symbol_nt_minus_one_delays_ns()
        coordinates = delay_ns_to_grid_coordinates(delays_ns, 576, 30_000.0)
        np.testing.assert_allclose(coordinates, np.linspace(0.0, 576.0, 8))
        pilot_subcarriers = np.arange(0, 576, 6, dtype=np.float64)
        phase = np.exp(
            -2j
            * np.pi
            * pilot_subcarriers[:, None]
            * coordinates[None, :]
            / 576.0
        )
        self.assertEqual(np.linalg.matrix_rank(phase, tol=1e-10), 7)

    def test_thick_sidon_packing_marks_a100_infeasible(self):
        requirements = thick_sidon_requirements(
            support_ns=479.66,
            grid_resolution_ns=57.87037037037037,
        )
        self.assertEqual(requirements["required_pair_gap_indices"], 17)
        self.assertEqual(requirements["required_fold_gap_indices"], 9)
        self.assertFalse(requirements["hard_feasible_by_packing_bound"])

    def test_unique_search_replays_and_canonicalizes(self):
        first = unique_random_delay_sets(20260726, 100)
        second = unique_random_delay_sets(20260726, 100)
        self.assertEqual(first, second)
        self.assertEqual(len(first), len(set(first)))
        self.assertTrue(all(value == canonical_delay_set(value) for value in first))


if __name__ == "__main__":
    unittest.main()
