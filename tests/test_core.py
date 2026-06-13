from __future__ import annotations

import unittest
import numpy as np

from cdd_lls.core.config import load_config
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
