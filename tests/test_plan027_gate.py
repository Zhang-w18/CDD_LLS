import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from tools.build_plan027_bler_gate import envelope_link_candidates
from tools.run_plan027_bler import (
    find_crossing_interval,
    validate_approved_manifest,
    validate_prior_scene_approval,
)
from tools.run_plan027_meff_design import read_csv


class Plan027GateTest(unittest.TestCase):
    def test_unapproved_manifest_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.json"
            manifest.write_text('{"candidates":[]}', encoding="utf-8")
            digest = hashlib.sha256(manifest.read_bytes()).hexdigest()
            with self.assertRaisesRegex(RuntimeError, "approval receipt is missing"):
                validate_approved_manifest(manifest, digest, root / "approval.json")

    def test_hash_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.json"
            manifest.write_text('{"candidates":[]}', encoding="utf-8")
            approval = root / "approval.json"
            approval.write_text(
                json.dumps({"approved": True, "manifest_sha256": "bad"}),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(RuntimeError, "does not match"):
                validate_approved_manifest(manifest, "bad", approval)

    def test_envelope_link_selection_has_53_physical_curves(self):
        targets = read_csv(
            Path("outputs")
            / "experiment027_meff_sidon"
            / "20260726_main"
            / "e1_outage"
            / "e1_outage_targets.csv"
        )
        rows = envelope_link_candidates(targets)
        self.assertEqual(len(rows), 53)
        counts = {
            scenario: sum(row["scenario_id"] == scenario for row in rows)
            for scenario in ("A1", "A5", "A10", "A30", "A100")
        }
        self.assertEqual(
            counts,
            {"A1": 10, "A5": 11, "A10": 12, "A30": 11, "A100": 9},
        )
        a5_rms = next(
            row
            for row in rows
            if row["scenario_id"] == "A5" and row["family"] == "AP_RMS_T1"
        )
        self.assertFalse(a5_rms["grid_aligned"])
        self.assertIsNone(a5_rms["delay_indices"])
        self.assertAlmostEqual(a5_rms["delay_grid_coordinates"][1], 0.0864)
        a5_tu = next(row for row in rows if row["candidate_id"] == "A5_AP_TU_NT")
        a5_alias = next(
            row for row in rows if row["candidate_id"] == "A5_AP_TALIAS_NT"
        )
        self.assertEqual(
            a5_tu["delay_indices"],
            [0, 72, 144, 216, 288, 360, 432, 504],
        )
        self.assertEqual(
            a5_alias["delay_indices"],
            [0, 3, 6, 9, 12, 15, 18, 21],
        )
        for scenario in ("A10", "A30"):
            families = {
                row["family"] for row in rows if row["scenario_id"] == scenario
            }
            self.assertIn("AP_TU_NT", families)
            self.assertIn("AP_TALIAS_NT", families)

    def test_candidate_specific_crossing_interval(self):
        rows = [
            {"candidate_id": "A", "snr_db": 14.0, "bler": 0.20},
            {"candidate_id": "A", "snr_db": 14.5, "bler": 0.09},
            {"candidate_id": "A", "snr_db": 15.0, "bler": 0.04},
        ]
        self.assertEqual(
            find_crossing_interval(rows, "A", 0.10),
            (14.0, 14.5),
        )

    def test_approved_47_curve_manifest_is_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.json"
            payload = {
                "schema": "plan027-e4-envelope-link-manifest-v2",
                "curve_count": 47,
                "dmrs_spacing_subcarriers": 24,
                "physical_definition": {
                    "phase_denominator": 576,
                    "delay_input_field": "delay_grid_coordinates",
                },
                "candidates": [],
            }
            manifest.write_text(json.dumps(payload), encoding="utf-8")
            digest = hashlib.sha256(manifest.read_bytes()).hexdigest()
            approval = root / "approval.json"
            approval.write_text(
                json.dumps(
                    {
                        "approved": True,
                        "manifest_sha256": digest,
                    }
                ),
                encoding="utf-8",
            )
            accepted = validate_approved_manifest(manifest, digest, approval)
            self.assertEqual(accepted["curve_count"], 47)

    def test_approved_49_curve_addendum_manifest_is_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.json"
            payload = {
                "schema": "plan027-e4-envelope-link-manifest-v3",
                "curve_count": 49,
                "dmrs_spacing_subcarriers": 24,
                "physical_definition": {
                    "phase_denominator": 576,
                    "delay_input_field": "delay_grid_coordinates",
                },
                "candidates": [],
            }
            manifest.write_text(json.dumps(payload), encoding="utf-8")
            digest = hashlib.sha256(manifest.read_bytes()).hexdigest()
            approval = root / "approval.json"
            approval.write_text(
                json.dumps({"approved": True, "manifest_sha256": digest}),
                encoding="utf-8",
            )
            accepted = validate_approved_manifest(manifest, digest, approval)
            self.assertEqual(accepted["curve_count"], 49)

    def test_approved_53_curve_addendum_manifest_is_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.json"
            payload = {
                "schema": "plan027-e4-envelope-link-manifest-v4",
                "curve_count": 53,
                "dmrs_spacing_subcarriers": 24,
                "physical_definition": {
                    "phase_denominator": 576,
                    "delay_input_field": "delay_grid_coordinates",
                },
                "candidates": [],
            }
            manifest.write_text(json.dumps(payload), encoding="utf-8")
            digest = hashlib.sha256(manifest.read_bytes()).hexdigest()
            approval = root / "approval.json"
            approval.write_text(
                json.dumps({"approved": True, "manifest_sha256": digest}),
                encoding="utf-8",
            )
            accepted = validate_approved_manifest(manifest, digest, approval)
            self.assertEqual(accepted["curve_count"], 53)

    def test_approved_a30_dense_dmrs_manifest_is_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.json"
            payload = {
                "schema": "plan027-e5-a30-dense-dmrs-manifest-v1",
                "curve_count": 9,
                "scene_counts": {"A30": 9},
                "dmrs_spacing_subcarriers": 12,
                "physical_definition": {
                    "phase_denominator": 576,
                    "delay_input_field": "delay_grid_coordinates",
                },
                "candidates": [],
            }
            manifest.write_text(json.dumps(payload), encoding="utf-8")
            digest = hashlib.sha256(manifest.read_bytes()).hexdigest()
            approval = root / "approval.json"
            approval.write_text(
                json.dumps({"approved": True, "manifest_sha256": digest}),
                encoding="utf-8",
            )
            accepted = validate_approved_manifest(manifest, digest, approval)
            self.assertEqual(accepted["dmrs_spacing_subcarriers"], 12)

    def test_approved_a100_comb6_manifest_is_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.json"
            payload = {
                "schema": "plan027-e6-a100-dense-dmrs-manifest-v1",
                "curve_count": 9,
                "scene_counts": {"A100": 9},
                "dmrs_spacing_subcarriers": 6,
                "physical_definition": {
                    "phase_denominator": 576,
                    "delay_input_field": "delay_grid_coordinates",
                },
                "candidates": [],
            }
            manifest.write_text(json.dumps(payload), encoding="utf-8")
            digest = hashlib.sha256(manifest.read_bytes()).hexdigest()
            approval = root / "approval.json"
            approval.write_text(
                json.dumps({"approved": True, "manifest_sha256": digest}),
                encoding="utf-8",
            )
            accepted = validate_approved_manifest(manifest, digest, approval)
            self.assertEqual(set(accepted["scene_counts"]), {"A100"})

    def test_a100_dense_manifest_rejects_comb12(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.json"
            payload = {
                "schema": "plan027-e6-a100-dense-dmrs-manifest-v1",
                "curve_count": 9,
                "scene_counts": {"A100": 9},
                "dmrs_spacing_subcarriers": 12,
                "physical_definition": {
                    "phase_denominator": 576,
                    "delay_input_field": "delay_grid_coordinates",
                },
                "candidates": [],
            }
            manifest.write_text(json.dumps(payload), encoding="utf-8")
            digest = hashlib.sha256(manifest.read_bytes()).hexdigest()
            approval = root / "approval.json"
            approval.write_text(
                json.dumps({"approved": True, "manifest_sha256": digest}),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(RuntimeError, "fixed to comb 6"):
                validate_approved_manifest(manifest, digest, approval)

    def test_a10_a30_execution_order_uses_completed_prior_requested_scene(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            a5_approval = root / "a5.json"
            a5_approval.write_text(
                json.dumps({"approved": True, "scenario_id": "A5"}),
                encoding="utf-8",
            )
            manifest = {"execution_order": ["A5", "A10", "A30"]}
            validate_prior_scene_approval(manifest, "A10", a5_approval)
            with self.assertRaisesRegex(RuntimeError, "completed scene A10"):
                validate_prior_scene_approval(manifest, "A30", a5_approval)


if __name__ == "__main__":
    unittest.main()
