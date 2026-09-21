from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np

from cdd_lls.phy.precoding import (
    build_active_dft_grid_precoder,
    build_aged_mrt_prg_precoder,
    build_prg_dft_precoder_batch,
    spatial_dft_codebook,
)
from cdd_lls.phy.resource_grid import build_resource_grid
from tools.run_plan027_bler import resource_config
from tools.run_plan033_tdl_mobility_mimo import (
    TRANSPARENT_COVARIANCE_SCALE,
    _delay_audit,
    _interval_identity,
    _next_interval_end,
    awgn_variance,
    candidate_definitions,
    delay_sets,
    load_config,
)


ROOT = Path(__file__).resolve().parents[1]


class Plan033MobilityMimoTests(unittest.TestCase):
    def test_four_configs_and_unit_power_model_have_no_tx_hardcoding(self):
        expected = {
            "plan033_nt8_nr4_v3_prescan.yaml": (8, 4, 3.0, "A100_NT8_NR4_V3"),
            "plan033_nt8_nr4_v60_prescan.yaml": (8, 4, 60.0, "A100_NT8_NR4_V60"),
            "plan033_nt4_nr4_v3_prescan.yaml": (4, 4, 3.0, "A100_NT4_NR4_V3"),
            "plan033_nt4_nr4_v60_prescan.yaml": (4, 4, 60.0, "A100_NT4_NR4_V60"),
        }
        grid = build_resource_grid(resource_config(6))
        for name, values in expected.items():
            config = load_config(ROOT / "configs" / name)
            self.assertEqual(
                (config["n_tx"], config["n_rx"], config["speed_kmh"], config["scenario_id"]),
                values,
            )
            self.assertTrue(config["precoder_normalize"])
            cdd = build_active_dft_grid_precoder(
                grid, delay_sets(config["n_tx"])["B0_QC"], config["n_tx"], normalize=True
            )
            np.testing.assert_allclose(np.sum(np.abs(cdd.C) ** 2, axis=1), 1.0, atol=2e-12)
        self.assertEqual(TRANSPARENT_COVARIANCE_SCALE, 1.0)
        self.assertAlmostEqual(awgn_variance(0.0), 1.0)
        self.assertAlmostEqual(awgn_variance(10.0), 0.1)

    def test_four_tx_delays_audit_and_dft_prg_mapping(self):
        grid = build_resource_grid(resource_config(6))
        frozen = delay_sets(4)
        self.assertEqual(frozen["B0_QC"], [0.0, 24.0, 48.0, 72.0])
        self.assertEqual(frozen["S0_SIDON"], [0.0, 1.0, 3.0, 7.0])
        self.assertEqual(frozen["SMALL_CDD_QSTEP0P25"], [0.0, 0.25, 0.5, 0.75])
        b0 = _delay_audit(grid, 4, "B0_QC", frozen["B0_QC"])
        sidon = _delay_audit(grid, 4, "S0_SIDON", frozen["S0_SIDON"])
        self.assertEqual(b0["comb6_residues_mod_96"], [0.0, 24.0, 48.0, 72.0])
        self.assertEqual((b0["minimum_fold_gap"], b0["pilot_rank"]), (24.0, 4))
        self.assertAlmostEqual(b0["pilot_condition_number"], 1.0, places=12)
        self.assertFalse(b0["unordered_pair_sums_unique"])
        self.assertEqual(sidon["unordered_pair_sums"], [0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 7.0, 8.0, 10.0, 14.0])
        self.assertTrue(sidon["unordered_pair_sums_unique"])
        order = np.asarray([[0, 1, 2, 3, 0, 1, 2, 3]])
        prg = build_prg_dft_precoder_batch(grid, 4, 6, order, normalize=True)
        codebook = spatial_dft_codebook(4, normalize=True)
        for prg_index, vector_index in enumerate(order[0]):
            np.testing.assert_allclose(prg[0, prg_index * 72], codebook[:, vector_index])

    def test_four_rx_aged_mrt_matches_manual_cross_rx_gram(self):
        old = np.zeros((1, 4, 4, 72), dtype=np.complex128)
        old[0, 0, 0, :] = 2.0
        old[0, 1:, 1, :] = 2.0
        weights = build_aged_mrt_prg_precoder(old, prg_size_rb=6, normalize=True)
        np.testing.assert_allclose(np.abs(weights[0, 0]), [0.0, 1.0, 0.0, 0.0], atol=1e-12)
        np.testing.assert_allclose(np.sum(np.abs(weights) ** 2, axis=2), 1.0, atol=2e-12)
        np.testing.assert_allclose(weights[:, :72], np.repeat(weights[:, :1], 72, axis=1))

    def test_six_candidates_share_trial_key_and_keep_output_identity(self):
        config = load_config(ROOT / "configs" / "plan033_nt4_nr4_v3_smoke.yaml")
        candidates = candidate_definitions(config["n_tx"])
        records = [
            _interval_identity(config, row["candidate_id"], -2.0, 1, 20)
            for row in candidates
        ]
        self.assertEqual(len(records), 6)
        self.assertIn("A100_NT4_AGED_MRT_PRG6", {record["candidate_id"] for record in records})
        self.assertIn(
            "A100_NT4_SMALL_CDD_QSTEP0P25_MATCHED",
            {record["candidate_id"] for record in records},
        )
        self.assertEqual(len({record["candidate_id"] for record in records}), 6)
        self.assertEqual(len({record["trial_key"] for record in records}), 1)
        self.assertEqual(
            {(record["n_tx"], record["n_rx"], record["speed_kmh"], record["receiver"]) for record in records},
            {(4, 4, 3.0, "estimated")},
        )
        self.assertEqual(_next_interval_end(0, 10_000, "formal"), 1_000)
        self.assertEqual(_next_interval_end(9_000, 10_000, "formal"), 10_000)
        self.assertEqual(_next_interval_end(0, 400, "prescan"), 400)


if __name__ == "__main__":
    unittest.main()
