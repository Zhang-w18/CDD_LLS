import csv
from pathlib import Path

import pytest
import yaml

from cdd_lls.sim.pdcch_cdd import _validate_candidate, validate_pdcch_cdd_config
from tools.run_plan031_strict_sidon_prescan import (
    COARSE_CONFIG_PATHS,
    CONFIG_PATHS,
    EXTENSION_CONFIG_PATHS,
    FINE_CONFIG_PATHS,
    FINE_SUPPLEMENT_CONFIG_PATHS,
    EXPECTED_TOTAL_CANDIDATES,
    _completed_points,
    load_jobs,
)


def test_strict_sidon_prescan_configs_match_frozen_shortlist() -> None:
    scenes, jobs = load_jobs()
    assert len(scenes) == 6
    assert len(jobs) == EXPECTED_TOTAL_CANDIDATES
    assert {scene["name"] for scene in scenes} == {
        "a100_al2",
        "a100_al4",
        "a100_al8",
        "c300_al1",
        "c300_al2",
        "c300_al4",
    }
    for config_path in CONFIG_PATHS:
        config = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
        assert config["seed"] == 20260908
        assert config["simulation"]["min_trials_per_snr"] == 300
        assert config["simulation"]["max_trials_per_snr"] == 300
        assert config["simulation"]["progress_every_batches"] == 1
        assert len(config["candidates"]) == 8
        assert all(
            row.get("receiver_covariance_mode", "matched_effective")
            == "matched_effective"
            for row in config["candidates"]
        )
        summary = validate_pdcch_cdd_config(config)
        assert len(summary["candidate_ids"]) == 8
        assert all(
            row["strict_sidon"]
            for row in summary["candidate_constraint_diagnostics"].values()
        )


def test_strict_sidon_validation_rejects_nonunique_pair_sums() -> None:
    candidate = {
        "candidate_id": "NOT_SIDON",
        "scheme": "cdd",
        "strict_sidon": True,
        "delay_grid_coordinates": [0, 1, 2, 4],
    }
    with pytest.raises(ValueError, match="pair sums"):
        _validate_candidate(candidate, n_tx=4, k_active=36, n_p=9)


def test_strict_sidon_validation_reports_complete_constraints() -> None:
    candidate = {
        "candidate_id": "SIDON",
        "scheme": "cdd",
        "strict_sidon": True,
        "delay_grid_coordinates": [0, 1, 3, 7],
    }
    diagnostics = _validate_candidate(candidate, n_tx=4, k_active=36, n_p=9)
    assert diagnostics == {
        "strict_sidon": True,
        "pair_sum_count": 10,
        "fold_residue_count": 4,
    }


def test_extension_contains_only_the_five_unbracketed_candidates() -> None:
    scenes, jobs = load_jobs(
        EXTENSION_CONFIG_PATHS,
        expected_total_candidates=5,
        require_full_shortlist=False,
    )
    assert len(scenes) == 3
    assert {job["candidate_id"] for job in jobs} == {
        "C300_SIDON4_AL1_05",
        "C300_SIDON4_AL1_06",
        "C300_SIDON4_AL1_08",
        "C300_SIDON4_AL2_05",
        "C300_SIDON4_AL4_05",
    }
    assert sum(job["snr_count"] for job in jobs) == 54


def test_coarse_confirmation_freezes_43_candidates_and_254_points() -> None:
    scenes, jobs = load_jobs(
        COARSE_CONFIG_PATHS,
        expected_total_candidates=43,
        require_full_shortlist=False,
    )
    assert len(scenes) == 6
    assert len(jobs) == 43
    assert sum(job["snr_count"] for job in jobs) == 254
    assert sum(
        job["snr_count"] * job["trials_per_point"] for job in jobs
    ) == 762_000
    assert all(job["trials_per_point"] == 3000 for job in jobs)
    for config_path in COARSE_CONFIG_PATHS:
        config = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
        summary = validate_pdcch_cdd_config(config)
        assert len(summary["candidate_ids"]) == len(config["candidates"])
        assert all(
            value["strict_sidon"]
            for value in summary["candidate_constraint_diagnostics"].values()
        )


def test_fine_scan_freezes_18_candidates_and_adaptive_budget() -> None:
    scenes, jobs = load_jobs(
        FINE_CONFIG_PATHS,
        expected_total_candidates=18,
        require_full_shortlist=False,
    )
    assert len(scenes) == 6
    assert len(jobs) == 18
    assert sum(job["snr_count"] for job in jobs) == 149
    assert all(job["minimum_trials_per_point"] == 10_000 for job in jobs)
    assert all(job["target_errors"] == 200 for job in jobs)
    assert all(job["maximum_trials_per_point"] == 50_000 for job in jobs)
    for config_path in FINE_CONFIG_PATHS:
        config = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
        assert config["random_stream_namespace"] == "plan031_strict_sidon_final_fine_v1"
        summary = validate_pdcch_cdd_config(config)
        assert len(summary["candidate_ids"]) == len(config["candidates"])


def test_adaptive_completion_accepts_target_or_maximum(tmp_path: Path) -> None:
    candidate_id = "SIDON"
    output = tmp_path / candidate_id
    output.mkdir()
    with (output / "bler_points.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["candidate_id", "snr_db", "trials", "errors"]
        )
        writer.writeheader()
        writer.writerow(
            {"candidate_id": candidate_id, "snr_db": 1.0, "trials": 10_000, "errors": 200}
        )
        writer.writerow(
            {"candidate_id": candidate_id, "snr_db": 2.0, "trials": 50_000, "errors": 100}
        )
    assert _completed_points(tmp_path, candidate_id, [1.0, 2.0], 10_000, 50_000, 200) == (
        2,
        60_000,
    )


def test_fine_supplement_contains_only_new_boundary_points() -> None:
    scenes, jobs = load_jobs(
        FINE_SUPPLEMENT_CONFIG_PATHS,
        expected_total_candidates=10,
        require_full_shortlist=False,
    )
    assert len(scenes) == 4
    assert len(jobs) == 10
    assert sum(job["snr_count"] for job in jobs) == 23
    assert all(job["minimum_trials_per_point"] == 10_000 for job in jobs)
    assert all(job["target_errors"] == 200 for job in jobs)
    assert all(job["maximum_trials_per_point"] == 50_000 for job in jobs)
