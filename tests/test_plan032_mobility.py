from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from tools.run_plan032_tdl_mobility import (
    _trial_target,
    build_aged_mrt_prg_precoder,
    load_config,
)


class Plan032MobilityTests(unittest.TestCase):
    def test_per_snr_trial_targets_are_loaded(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "config.yaml"
            path.write_text(
                "\n".join(
                    [
                        "schema: plan032-tdl-mobility-v1",
                        "run_kind: formal",
                        "output_dir: outputs/test",
                        "batch_size: 20",
                        "snr_db: [16, 16.25]",
                        "target_total_trials_by_snr:",
                        "  16: 10000",
                        "  16.25: 17660",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            config = load_config(path)
        self.assertEqual(_trial_target(config, 16.0), 10000)
        self.assertEqual(_trial_target(config, 16.25), 17660)

    def test_aged_mrt_is_power_normalized_and_prg_constant(self):
        rng = np.random.default_rng(32)
        old = rng.normal(size=(4, 1, 8, 576)) + 1j * rng.normal(
            size=(4, 1, 8, 576)
        )
        weights = build_aged_mrt_prg_precoder(old, prg_size_rb=6)
        self.assertEqual(weights.shape, (4, 576, 8))
        np.testing.assert_allclose(
            np.sum(np.abs(weights) ** 2, axis=2), 8.0, atol=2e-12, rtol=0.0
        )
        for start in range(0, 576, 72):
            np.testing.assert_allclose(
                weights[:, start : start + 72],
                np.repeat(weights[:, start : start + 1], 72, axis=1),
                atol=0.0,
                rtol=0.0,
            )

    def test_aged_mrt_matches_prg_gram_dominant_eigenvalue(self):
        rng = np.random.default_rng(320)
        old = rng.normal(size=(2, 1, 8, 576)) + 1j * rng.normal(
            size=(2, 1, 8, 576)
        )
        weights = build_aged_mrt_prg_precoder(old, prg_size_rb=6)
        h = old[:, 0, :, :72]
        gram = np.einsum("bik,bjk->bij", h.conj(), h, optimize=True)
        eigenvalues = np.linalg.eigvalsh(gram)[:, -1]
        achieved = np.einsum(
            "bi,bij,bj->b",
            weights[:, 0].conj(),
            gram,
            weights[:, 0],
            optimize=True,
        ).real
        np.testing.assert_allclose(achieved, 8.0 * eigenvalues, atol=1e-9, rtol=1e-12)



if __name__ == "__main__":
    unittest.main()
