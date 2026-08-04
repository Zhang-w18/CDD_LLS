import unittest

import numpy as np

from cdd_lls.design.cdd_metrics import (
    delay_indices_to_ns,
    fold_min_gap_indices,
    pair_sum_min_gap_indices,
    tdl_effective_support,
)
from cdd_lls.design.cdd_search import (
    arithmetic_delay_sets,
    canonical_delay_set,
    residues_and_lifts,
)


class CDDDesignTest(unittest.TestCase):
    def test_canonicalization_removes_permutation_and_common_shift(self):
        first = canonical_delay_set([0, 9, 18, 27, 36, 45, 54, 63])
        second = canonical_delay_set([117, 90, 108, 72, 99, 63, 81, 126])
        self.assertEqual(first, second)

    def test_sidon_and_uniform_fold_geometry(self):
        sidon = [0, 1, 3, 7, 12, 20, 30, 65]
        uniform = [0, 3, 6, 9, 12, 15, 18, 21]
        self.assertEqual(pair_sum_min_gap_indices(sidon), 1)
        self.assertEqual(fold_min_gap_indices(uniform), 3)
        self.assertEqual(fold_min_gap_indices([0, 24, 1, 2, 3, 4, 5, 6]), 0)

    def test_residue_lift_and_delay_conversion(self):
        residues, lifts = residues_and_lifts([0, 27, 54, 9])
        self.assertEqual(residues, (0, 3, 6, 9))
        self.assertEqual(lifts, (0, 1, 2, 0))
        delays = delay_indices_to_ns([0, 1, 9])
        np.testing.assert_allclose(delays, [0.0, 57.87037037037037, 520.8333333333334])

    def test_arithmetic_sets_are_unique_canonical_representatives(self):
        candidates = arithmetic_delay_sets()
        self.assertEqual(len(candidates), len(set(candidates)))
        self.assertTrue(all(candidate[0] == 0 for candidate in candidates))

    def test_shortest_effective_support_can_drop_weak_early_path(self):
        support = tdl_effective_support(
            delays_s=[0.0, 2e-9, 4e-9, 20e-9],
            powers=[0.01, 0.60, 0.35, 0.04],
            epsilon=0.05,
        )
        self.assertAlmostEqual(float(support["support_start_s"]), 2e-9)
        self.assertAlmostEqual(float(support["support_stop_s"]), 4e-9)
        self.assertGreaterEqual(float(support["contained_power"]), 0.95)


if __name__ == "__main__":
    unittest.main()
