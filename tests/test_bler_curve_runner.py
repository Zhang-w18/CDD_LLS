import json
from pathlib import Path

import pytest

from tools import run_bler_curves as runner


def _base(candidate: str = "A30_B0_QC", trials: int = 400, errors: int = 40) -> dict:
    return {
        "scenario_id": "A30",
        "candidate_id": candidate,
        "family": "B0_QC",
        "snr_db": 14.0,
        "trials": trials,
        "tb_errors": errors,
        "bler": errors / trials,
        "ce_nmse_mean": 0.01,
        "ce_nmse_mean_db": -20.0,
    }


def test_merge_adds_counts_and_linear_nmse() -> None:
    extra = {
        "candidate_id": "A30_B0_QC",
        "snr_db": 14.0,
        "trial_start": 401,
        "trial_end": 1000,
        "trials": 600,
        "tb_errors": 30,
        "ce_nmse_sum": 12.0,
    }
    row = runner.merge_rows([_base()], [extra])[0]
    assert int(row["trials"]) == 1000
    assert int(row["tb_errors"]) == 70
    assert float(row["bler"]) == pytest.approx(0.07)
    assert float(row["ce_nmse_mean"]) == pytest.approx((4.0 + 12.0) / 1000.0)


def test_merge_rejects_overlap_or_gap() -> None:
    extra = {
        "candidate_id": "A30_B0_QC",
        "snr_db": 14.0,
        "trial_start": 500,
        "trial_end": 900,
        "trials": 401,
        "tb_errors": 10,
    }
    with pytest.raises(RuntimeError, match="expected start 401"):
        runner.merge_rows([_base()], [extra])


def test_adaptive_task_targets_error_count_and_batch_alignment() -> None:
    tasks = runner.build_tasks(
        [_base()],
        {"target_errors": 100, "max_total_trials": 10000},
        batch_size=100,
    )
    assert tasks == [
        {
            "candidate_id": "A30_B0_QC",
            "snr_db": 14.0,
            "trial_start": 401,
            "trial_end": 1000,
            "trials": 600,
        }
    ]


def test_plot_cutoff_keeps_first_point_below_threshold() -> None:
    rows = []
    for snr, bler in [(14.0, 0.02), (14.5, 0.009), (15.0, 0.004), (15.5, 0.002)]:
        row = _base()
        row.update({"snr_db": snr, "bler": bler})
        rows.append(row)
    keys = runner._plot_keys(rows, 0.005)
    assert ("A30_B0_QC", 15.0) in keys
    assert ("A30_B0_QC", 15.5) not in keys


def test_candidate_override_extends_tail_and_trial_budget() -> None:
    rows = []
    for candidate, errors in (("A100_B0_QC", 2), ("A100_AP_RMS_T1", 2)):
        for snr, bler in ((16.0, 0.006), (16.5, 0.002), (17.0, 0.001)):
            row = _base(candidate=candidate, trials=1000, errors=errors)
            row.update({"snr_db": snr, "bler": bler})
            rows.append(row)
    policy = {
        "target_errors": 50,
        "max_total_trials": 5000,
        "stop_below_bler": 0.005,
        "candidate_overrides": {
            "A100_AP_RMS_T1": {
                "target_errors": 100,
                "max_total_trials": 20000,
                "stop_below_bler": None,
            }
        },
    }
    tasks = runner.build_tasks(rows, policy, batch_size=100)
    by_key = {(row["candidate_id"], row["snr_db"]): row for row in tasks}
    assert ("A100_B0_QC", 17.0) not in by_key
    assert by_key[("A100_B0_QC", 16.5)]["trial_end"] == 5000
    assert by_key[("A100_AP_RMS_T1", 17.0)]["trial_end"] == 20000


def test_load_config_requires_schema(tmp_path: Path) -> None:
    path = tmp_path / "bad.yaml"
    path.write_text("schema: wrong\nscenes: []\n", encoding="utf-8")
    with pytest.raises(ValueError, match="config.schema"):
        runner.load_runner_config(path)


def test_manifest_approval_is_schema_agnostic(tmp_path: Path) -> None:
    manifest = {
        "schema": "future-bler-manifest-v9",
        "scene_counts": {"A30": 1},
        "dmrs_spacing_subcarriers": 6,
        "physical_definition": {
            "phase_denominator": 576,
            "delay_input_field": "delay_grid_coordinates",
        },
        "candidates": [
            {
                "scenario_id": "A30",
                "candidate_id": "A30_TEST",
                "family": "TEST",
                "dmrs_spacing_subcarriers": 6,
                "delay_grid_coordinates": [0, 1, 2, 3, 4, 5, 6, 7],
            }
        ],
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    digest = runner._sha256(manifest_path)
    approval_path = tmp_path / "approval.json"
    approval_path.write_text(
        json.dumps({"approved": True, "manifest_sha256": digest}), encoding="utf-8"
    )
    scene = {
        "scenario_id": "A30",
        "manifest": {
            "path": str(manifest_path),
            "approval_path": str(approval_path),
        },
        "candidate_ids": "all",
    }
    loaded, actual = runner._load_manifest(scene)
    assert loaded["schema"] == "future-bler-manifest-v9"
    assert actual == digest
