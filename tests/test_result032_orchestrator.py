from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from tools.run_result032_result028_scale import read_completed_trials


class Result032OrchestratorTests(unittest.TestCase):
    @staticmethod
    def _write_intervals(path: Path, rows: list[dict]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["candidate_id", "snr_db", "trial_start", "trial_end"],
            )
            writer.writeheader()
            writer.writerows(rows)

    def test_completed_trials_sum_contiguous_intervals(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            rows = []
            for candidate in range(12):
                rows.extend(
                    [
                        {
                            "candidate_id": f"candidate-{candidate}",
                            "snr_db": 16,
                            "trial_start": 1,
                            "trial_end": 1000,
                        },
                        {
                            "candidate_id": f"candidate-{candidate}",
                            "snr_db": 16,
                            "trial_start": 1001,
                            "trial_end": 3000,
                        },
                    ]
                )
            self._write_intervals(output / "intervals.csv", rows)
            self.assertEqual(read_completed_trials(output), {16.0: 3000})

    def test_trial_gap_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            rows = []
            for candidate in range(12):
                rows.extend(
                    [
                        {
                            "candidate_id": f"candidate-{candidate}",
                            "snr_db": 16,
                            "trial_start": 1,
                            "trial_end": 1000,
                        },
                        {
                            "candidate_id": f"candidate-{candidate}",
                            "snr_db": 16,
                            "trial_start": 1002,
                            "trial_end": 3000,
                        },
                    ]
                )
            self._write_intervals(output / "intervals.csv", rows)
            with self.assertRaisesRegex(RuntimeError, "gap/overlap"):
                read_completed_trials(output)


if __name__ == "__main__":
    unittest.main()
