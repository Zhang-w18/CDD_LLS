from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.plot_result032_selected_bler import MRT, PRG, _plot, _plot_nmse


class Result032SelectedBlerPlotTests(unittest.TestCase):
    def test_plot_accepts_original_and_transparent_inputs(self):
        cdd_id = "A100_S0_SIDON"
        transparent_id = "A100_S0_SIDON_TRANSPARENT_CDD"
        original = []
        transparent = []
        for snr_db, bler in ((14.0, 0.2), (15.0, 0.05)):
            original.extend(
                {
                    "candidate_id": candidate_id,
                    "snr_db": str(snr_db),
                    "bler": str(bler),
                    "trials": "10000",
                }
                for candidate_id in (cdd_id, PRG, MRT)
            )
            transparent.append(
                {
                    "candidate_id": transparent_id,
                    "snr_db": str(snr_db),
                    "bler": str(bler * 1.2),
                    "trials": "10000",
                }
            )
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "selected.png"
            _plot(original, transparent, cdd_id, transparent_id, output)
            self.assertTrue(output.is_file())
            self.assertGreater(output.stat().st_size, 0)

    def test_nmse_plot_accepts_both_cdd_receiver_modes(self):
        original = []
        transparent = []
        for snr_db in (14.0, 20.0):
            for candidate_id, nmse_db in (
                ("A100_S0_SIDON", -18.0),
                ("A100_B0_QC", -17.0),
            ):
                original.append(
                    {
                        "candidate_id": candidate_id,
                        "snr_db": str(snr_db),
                        "ce_nmse_mean_db": str(nmse_db),
                    }
                )
            for candidate_id, nmse_db in (
                ("A100_S0_SIDON_TRANSPARENT_CDD", -2.8),
                ("A100_B0_QC_TRANSPARENT_CDD", -1.0),
            ):
                transparent.append(
                    {
                        "candidate_id": candidate_id,
                        "snr_db": str(snr_db),
                        "ce_nmse_mean_db": str(nmse_db),
                    }
                )
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "nmse.png"
            _plot_nmse(original, transparent, output, preview=True)
            self.assertTrue(output.is_file())
            self.assertGreater(output.stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
