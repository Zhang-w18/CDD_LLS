from __future__ import annotations

import unittest

import numpy as np

from tools.run_v_design_piecewise_tradeoff import (
    balanced_random_permutation_indices,
    centered_delay_alphabet,
)
from tools.run_v_design_new_front_link import apply_interleaver, cdd_reference, undo_interleaver


class VDesignTests(unittest.TestCase):
    def test_cdd_reference_depends_on_active_bandwidth(self):
        self.assertEqual(cdd_reference(24, "N3"), ("C6", 7))
        self.assertEqual(cdd_reference(24, "N23"), ("C4", 5))
        self.assertEqual(cdd_reference(36, "N3"), ("C5", 6))
        self.assertEqual(cdd_reference(48, "N3"), ("C0", 1))

    def test_coded_bit_interleaver_round_trip(self):
        bits = np.arange(32, dtype=np.float64)
        permutation = np.random.default_rng(17).permutation(bits.size)
        mapped = apply_interleaver(bits, permutation, "full")
        restored = undo_interleaver(mapped, permutation, "full")
        np.testing.assert_array_equal(restored, bits)

    def test_centered_delay_alphabet_has_positive_and_negative_slopes(self):
        delays = centered_delay_alphabet(n_tx=8, slope_step=4)
        expected = np.array([-14, -10, -6, -2, 2, 6, 10, 14], dtype=np.float64)
        np.testing.assert_array_equal(delays, expected)

    def test_balanced_random_permutations_are_unique_by_segment_and_tx(self):
        idx = balanced_random_permutation_indices(
            n_segments=8,
            n_tx=8,
            rng=np.random.default_rng(20260619),
        )

        self.assertEqual(idx.shape, (8, 8))
        expected = list(range(8))
        for row in idx:
            self.assertEqual(sorted(row.astype(int).tolist()), expected)
        for column in idx.T:
            self.assertEqual(sorted(column.astype(int).tolist()), expected)

    def test_balanced_random_permutations_reject_impossible_constraint(self):
        with self.assertRaises(ValueError):
            balanced_random_permutation_indices(
                n_segments=9,
                n_tx=8,
                rng=np.random.default_rng(20260619),
            )


if __name__ == "__main__":
    unittest.main()
