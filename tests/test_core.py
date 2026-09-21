from __future__ import annotations

import unittest
import numpy as np

from cdd_lls.core.config import config_from_dict, dataclass_to_dict, load_config
from cdd_lls.core.mcs import build_tb_layout, get_mcs
from cdd_lls.phy.estimators import shifted_pdp
from cdd_lls.phy.qam import qam_demapper_maxlog, qam_modulate
from cdd_lls.phy.resource_grid import build_resource_grid


class CoreTests(unittest.TestCase):
    def test_smoke_config_loads(self):
        cfg = load_config("configs/smoke.yaml")
        grid = build_resource_grid(cfg.resource)
        self.assertGreater(grid.n_data_re, 0)
        self.assertGreater(grid.pilot_count, 0)

    def test_cdl_smoke_config_loads_and_validates_array_size(self):
        cfg = load_config("configs/smoke_sionna_cdl.yaml")
        self.assertEqual(cfg.channel.backend, "sionna_cdl")
        self.assertEqual(cfg.channel.cdl_profile, "A")
        invalid = dataclass_to_dict(cfg)
        invalid["channel"]["cdl_tx_array_cols"] = 4
        with self.assertRaisesRegex(ValueError, "CDL Tx array"):
            config_from_dict(invalid)

    def test_tb_layout_matches_resource_bits(self):
        cfg = load_config("configs/smoke.yaml")
        grid = build_resource_grid(cfg.resource)
        mcs = get_mcs(cfg.mcs.table, cfg.mcs.index)
        tb = build_tb_layout(grid.n_data_re, mcs)
        self.assertEqual(sum(tb.cb_e_values), grid.n_data_re * tb.qm)
        self.assertEqual(sum(tb.cb_k_values), tb.tb_size)

    def test_qam_llr_sign_convention_round_trip(self):
        bits = np.array([0, 0, 0, 1, 1, 1, 1, 0], dtype=np.int8)
        syms = qam_modulate(bits, qm=4)
        llr = qam_demapper_maxlog(syms, noise_var=0.01, qm=4)
        hard = (llr > 0).astype(np.int8)
        self.assertTrue(np.array_equal(hard, bits))

    def test_shifted_pdp_normalizes_energy(self):
        pdp = np.array([0.7, 0.2, 0.1])
        sg = shifted_pdp(pdp, [0, 2])
        self.assertTrue(np.isclose(np.sum(sg), 1.0))
        self.assertEqual(len(sg), 5)


if __name__ == "__main__":
    unittest.main()
