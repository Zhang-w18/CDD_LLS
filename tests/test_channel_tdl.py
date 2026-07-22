from __future__ import annotations

import unittest

import numpy as np

from cdd_lls.core.config import ChannelConfig, ResourceConfig
from cdd_lls.phy.channel_tdl import generate_sionna_tdl_channel, speed_kmh_to_mps
from cdd_lls.phy.precoding import build_precoder, equivalent_channel
from cdd_lls.phy.resource_grid import build_resource_grid
from cdd_lls.core.config import TransmissionConfig


class SionnaTDLChannelTests(unittest.TestCase):
    def setUp(self):
        self.grid = build_resource_grid(ResourceConfig(n_prbs=8))
        self.channel = ChannelConfig(
            backend="sionna_tdl",
            tdl_profile="A",
            delay_spread_ns=100.0,
            carrier_frequency_hz=3.5e9,
            ue_speed_kmh=3.0,
            num_sinusoids=20,
        )

    def test_speed_conversion(self):
        self.assertAlmostEqual(speed_kmh_to_mps(3.0), 3.0 / 3.6)
        self.assertAlmostEqual(speed_kmh_to_mps(60.0), 60.0 / 3.6)

    def test_shape_energy_and_seed_replay(self):
        first = generate_sionna_tdl_channel(self.grid, self.channel, 2, 1, batch_size=64, seed=1729)
        replay = generate_sionna_tdl_channel(self.grid, self.channel, 2, 1, batch_size=64, seed=1729)
        self.assertEqual(first.H.shape, (64, 1, 2, 10, 96))
        np.testing.assert_array_equal(first.H, replay.H)
        energy = float(np.mean(np.abs(first.H) ** 2))
        self.assertLess(abs(energy - 1.0), 0.20)

    def test_cdd_time_dimension_matches_each_symbol_slice(self):
        realization = generate_sionna_tdl_channel(self.grid, self.channel, 2, 1, batch_size=1, seed=9)
        tx = TransmissionConfig(cdd_delay_vector=[0, 9])
        precoder = build_precoder(self.grid, ResourceConfig(n_prbs=8), tx, 2)
        all_symbols = equivalent_channel(realization.H, precoder.C)
        for symbol in range(self.grid.n_symbols):
            static = equivalent_channel(realization.H[0, :, :, symbol, :], precoder.C)
            np.testing.assert_allclose(all_symbols[0, :, symbol, :], static, atol=1e-13, rtol=1e-13)

    def test_cdd_phase_uses_negative_frequency_slope(self):
        delays = [0, 9]
        tx = TransmissionConfig(cdd_delay_vector=delays)
        precoder = build_precoder(self.grid, ResourceConfig(n_prbs=8), tx, 2)
        local = 17
        k = float(self.grid.subcarrier_indices[local])
        expected = np.exp(-1j * 2.0 * np.pi * k * np.asarray(delays) / self.grid.n_fft) / np.sqrt(2.0)
        np.testing.assert_allclose(precoder.C[local], expected, atol=1e-15, rtol=1e-15)


if __name__ == "__main__":
    unittest.main()
