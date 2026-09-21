from __future__ import annotations

import unittest

from tools.run_plan033_append import plan_next_step


class Plan033AppendTests(unittest.TestCase):
    def test_one_step_groups_paired_snr_and_stops_completed_endpoints(self) -> None:
        requirements = [
            {"scenario_id": "scene", "candidate_id": "a", "snr_db": 6.0},
            {"scenario_id": "scene", "candidate_id": "b", "snr_db": 6.0},
            {"scenario_id": "scene", "candidate_id": "c", "snr_db": 6.25},
        ]
        rows = []
        for candidate, errors_6, errors_625 in (("a", 199, 10), ("b", 220, 10), ("c", 10, 200)):
            rows.extend(
                [
                    {"candidate_id": candidate, "snr_db": "6", "trials": "10000", "tb_errors": str(errors_6)},
                    {"candidate_id": candidate, "snr_db": "6.25", "trials": "10000", "tb_errors": str(errors_625)},
                ]
            )
        state = plan_next_step(requirements, rows)
        self.assertEqual(state["status"], "pending")
        self.assertEqual(state["next_target_total_trials_by_snr"], {6.0: 11000})

    def test_max_trials_is_a_valid_stop(self) -> None:
        requirements = [{"scenario_id": "scene", "candidate_id": "a", "snr_db": 7.0}]
        rows = [
            {"candidate_id": candidate, "snr_db": "7", "trials": "50000", "tb_errors": "150"}
            for candidate in ("a", "b", "c", "d", "e", "f")
        ]
        state = plan_next_step(requirements, rows)
        self.assertEqual(state["status"], "complete")
        self.assertEqual(state["next_target_total_trials_by_snr"], {})
        self.assertEqual(state["endpoints"][0]["completion_reason"], "max_total_trials")

    def test_rejects_unpaired_trial_counts(self) -> None:
        requirements = [{"scenario_id": "scene", "candidate_id": "a", "snr_db": 7.0}]
        rows = [
            {"candidate_id": "a", "snr_db": "7", "trials": "10000", "tb_errors": "100"},
            {"candidate_id": "b", "snr_db": "7", "trials": "11000", "tb_errors": "100"},
        ]
        with self.assertRaisesRegex(ValueError, "Paired trial counts differ"):
            plan_next_step(requirements, rows)


if __name__ == "__main__":
    unittest.main()
