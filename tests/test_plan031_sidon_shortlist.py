import unittest

from tools.prepare_plan031_sidon_shortlist import (
    S0,
    canonical_delay_set,
    choose_shortlist,
    enumerate_four_tx,
    is_valid,
)


class Plan031SidonShortlistTest(unittest.TestCase):
    def test_known_sidon_and_non_sidon_sets(self):
        self.assertTrue(is_valid(S0[4], period=36, pilot_period=9))
        self.assertTrue(is_valid(S0[8], period=144, pilot_period=36))
        self.assertFalse(is_valid((0, 1, 2, 3), period=36, pilot_period=9))
        self.assertFalse(is_valid((0, 1, 3, 9), period=36, pilot_period=9))

    def test_four_tx_al1_complete_counts(self):
        candidates, valid_with_zero = enumerate_four_tx(period=36, pilot_period=9)
        self.assertEqual(valid_with_zero, 1920)
        self.assertEqual(len(candidates), 480)
        self.assertEqual(len(candidates), len(set(candidates)))
        self.assertTrue(all(is_valid(row, 36, 9) for row in candidates))

    def test_shortlist_is_deterministic_and_contains_s0(self):
        candidates, _ = enumerate_four_tx(period=36, pilot_period=9)
        first, _ = choose_shortlist(candidates, 36, 9, 4, 8)
        second, _ = choose_shortlist(candidates, 36, 9, 4, 8)
        self.assertEqual(first, second)
        self.assertEqual(first[0], canonical_delay_set(S0[4], 36))
        self.assertEqual(len(first), 8)
        self.assertEqual(len(first), len(set(first)))


if __name__ == "__main__":
    unittest.main()
