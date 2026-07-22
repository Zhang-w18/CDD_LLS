from __future__ import annotations

import unittest

import numpy as np

from cdd_lls.core.config import ResourceConfig
from cdd_lls.phy.resource_grid import (
    build_resource_grid,
    full_fft_indices_for_subcarriers,
)


class ResourceGridTimeFrequencyTests(unittest.TestCase):
    def test_8_and_48_prb_coordinate_counts(self):
        for n_prbs, n_sc in ((8, 96), (48, 576)):
            with self.subTest(n_prbs=n_prbs):
                grid = build_resource_grid(ResourceConfig(
                    n_prbs=n_prbs,
                    pdsch_n_symbols=10,
                    dmrs_symbol_indices=[2, 7],
                    dmrs_spacing_sc=6,
                    dmrs_offset_sc=0,
                ))
                pilots_per_symbol = (n_sc + 5) // 6
                self.assertEqual(grid.n_sc, n_sc)
                self.assertEqual(grid.n_symbols, 10)
                self.assertEqual(grid.pilot_re_count, 2 * pilots_per_symbol)
                self.assertEqual(grid.n_data_re, 10 * n_sc - 2 * pilots_per_symbol)
                self.assertEqual(len(np.unique(grid.pilot_coordinates, axis=0)), grid.n_dmrs_re)
                self.assertEqual(len(np.unique(grid.data_coordinates, axis=0)), grid.n_data_re)
                pilot_set = {tuple(x) for x in grid.pilot_coordinates.tolist()}
                data_set = {tuple(x) for x in grid.data_coordinates.tolist()}
                self.assertFalse(pilot_set & data_set)

    def test_full_fft_mapping_is_centered_and_contiguous(self):
        grid = build_resource_grid(ResourceConfig(n_prbs=8, n_fft=4096))
        mapped = full_fft_indices_for_subcarriers(grid, grid.subcarrier_indices)
        np.testing.assert_array_equal(mapped, grid.active_fft_indices)
        self.assertEqual(int(mapped[0]), 4096 // 2 - 48)
        self.assertEqual(int(mapped[-1]), 4096 // 2 + 47)

    def test_ofdm_symbol_duration_includes_cyclic_prefix(self):
        grid = build_resource_grid(ResourceConfig(
            n_fft=4096,
            scs_khz=30,
            cyclic_prefix_length=288,
        ))
        expected = (1.0 + 288.0 / 4096.0) / 30e3
        self.assertAlmostEqual(grid.ofdm_symbol_duration_s, expected, places=15)


if __name__ == "__main__":
    unittest.main()
