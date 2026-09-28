from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from tools.run_plan032_tdl_mobility import (
    TRANSPARENT_CDD_IDS,
    TRANSPARENT_CDD_SUPPLEMENT,
    _build_filters,
    _candidate_definitions,
    _trial_target,
    build_aged_mrt_prg_precoder,
    load_config,
)


class Plan032MobilityTests(unittest.TestCase):
    def test_transparent_cdd_supplement_defines_only_requested_sources(self):
        manifest = {
            "scene_counts": {"A100": 10},
            "dmrs_spacing_subcarriers": 6,
            "candidates": [
                {
                    "scenario_id": "A100",
                    "candidate_id": candidate_id,
                    "delay_grid_coordinates": [0] * 8,
                    "dmrs_spacing_subcarriers": 6,
                }
                for candidate_id in [
                    "A100_S0_SIDON",
                    "A100_B0_QC",
                    *[f"unused-{index}" for index in range(8)],
                ]
            ],
        }
        candidates = _candidate_definitions(manifest, TRANSPARENT_CDD_SUPPLEMENT)
        self.assertEqual(
            [row["candidate_id"] for row in candidates],
            [
                TRANSPARENT_CDD_IDS["A100_S0_SIDON"],
                TRANSPARENT_CDD_IDS["A100_B0_QC"],
            ],
        )
        self.assertTrue(
            all(row["family"] == "TRANSPARENT_CDD_PHYSCOV" for row in candidates)
        )

    def test_transparent_cdd_supplement_shares_fullband_physical_filter(self):
        candidates = [
            {
                "candidate_id": candidate_id,
                "family": "TRANSPARENT_CDD_PHYSCOV",
            }
            for candidate_id in TRANSPARENT_CDD_IDS.values()
        ]
        base = SimpleNamespace(
            time=np.eye(2, dtype=np.complex128),
            frequency=np.eye(3, dtype=np.complex128),
        )
        shared_filter = object()
        with patch(
            "tools.run_plan032_tdl_mobility.build_time_frequency_rmmse_filter",
            return_value=shared_filter,
        ) as build:
            filters = _build_filters(
                grid=object(),
                base=base,
                precoders={},
                noise_variance=0.25,
                candidates=candidates,
            )
        self.assertEqual(set(filters), set(TRANSPARENT_CDD_IDS.values()))
        self.assertTrue(all(value is shared_filter for value in filters.values()))
        covariance = build.call_args.args[1]
        np.testing.assert_allclose(covariance.frequency, 8.0 * base.frequency)
        self.assertEqual(covariance.covariance_type, "transparent_8x_physical")

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
