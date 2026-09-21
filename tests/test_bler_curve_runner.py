import json
from pathlib import Path

import numpy as np
import pytest

from tools import run_bler_curves as runner
from tools import search_plan028_cdd_nonuniform as nonuniform_search


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


def test_unpaired_tasks_cap_each_resumable_interval() -> None:
    row = _base()
    row.update({"trials": 0, "tb_errors": 0, "bler": 0.0})
    tasks = runner.build_tasks(
        [row],
        {
            "target_errors": 200,
            "min_total_trials": 10000,
            "max_total_trials": 50000,
            "max_interval_trials": 1000,
        },
        batch_size=25,
    )
    assert tasks[0]["trial_start"] == 1
    assert tasks[0]["trial_end"] == 1000


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


def test_load_config_defaults_to_one_rx_and_rejects_invalid_n_rx(tmp_path: Path) -> None:
    valid = tmp_path / "valid.yaml"
    valid.write_text(
        "schema: bler-curve-runner-v1\n"
        f"output_dir: {tmp_path.as_posix()}/out\n"
        "scenes:\n  - scenario_id: A100\n",
        encoding="utf-8",
    )
    loaded = runner.load_runner_config(valid)
    assert loaded["scenes"][0]["n_rx"] == 1

    for index, value in enumerate((0, -1, True, 1.5, "4")):
        invalid = tmp_path / f"invalid_n_rx_{index}.yaml"
        invalid.write_text(
            "schema: bler-curve-runner-v1\n"
            f"output_dir: {tmp_path.as_posix()}/out\n"
            f"scenes:\n  - scenario_id: A100\n    n_rx: {json.dumps(value)}\n",
            encoding="utf-8",
        )
        with pytest.raises(ValueError, match="scene.n_rx"):
            runner.load_runner_config(invalid)


def test_four_rx_ce_ratio_and_mrc_match_manual_sums() -> None:
    rng = np.random.default_rng(20260811)
    truth = rng.normal(size=(3, 4, 7)) + 1j * rng.normal(size=(3, 4, 7))
    estimate = truth + 0.1 * (
        rng.normal(size=truth.shape) + 1j * rng.normal(size=truth.shape)
    )
    received = truth * (0.7 - 0.2j) + 0.03 * (
        rng.normal(size=truth.shape) + 1j * rng.normal(size=truth.shape)
    )
    noise_variance = 0.4

    error, signal, ratio = runner._ce_trial_statistics(estimate, truth)
    np.testing.assert_allclose(error, np.sum(np.abs(estimate - truth) ** 2, axis=(1, 2)))
    np.testing.assert_allclose(signal, np.sum(np.abs(truth) ** 2, axis=(1, 2)))
    np.testing.assert_allclose(ratio, error / signal)

    equalized, effective_noise = runner._estimated_csi_mrc_equalize(
        estimate, received, noise_variance
    )
    denominator = np.sum(np.abs(estimate) ** 2, axis=1)
    np.testing.assert_allclose(
        equalized, np.sum(np.conj(estimate) * received, axis=1) / denominator
    )
    np.testing.assert_allclose(effective_noise, noise_variance / denominator)

    ideal_equalized, ideal_noise = runner.curves028.ideal_csi_equalize(
        truth, received, noise_variance
    )
    ideal_denominator = np.sum(np.abs(truth) ** 2, axis=1)
    np.testing.assert_allclose(
        ideal_equalized, np.sum(np.conj(truth) * received, axis=1) / ideal_denominator
    )
    np.testing.assert_allclose(ideal_noise, noise_variance / ideal_denominator)


def test_section14_transparent_prg_manifest_records_four_rx(tmp_path: Path) -> None:
    config = {"output_dir": str(tmp_path), "seed": 20260727}
    scene = {
        "scenario_id": "A100",
        "n_rx": 4,
        "plan_section": 14,
        "paired_policy": {
            "target_errors": 1,
            "min_total_trials": 300,
            "max_total_trials": 300,
        },
    }
    candidate = {
        "candidate_id": "A100_PRG_DFT8_6RB",
        "family": "TRANSPARENT_PRG_DFT",
        "prg_size_rb": 6,
        "mapping": "cycle_all",
    }
    path, _ = runner._write_transparent_manifest(
        config, scene, "source-sha256", [candidate]
    )
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["system"]["n_rx"] == 4
    assert payload["plan"].startswith("research/plan-028-comb6三类CSI-TDL-A300ns.md#14-")


def test_section15_transparent_prg_manifest_records_extension_provenance(tmp_path: Path) -> None:
    config = {"output_dir": str(tmp_path), "seed": 20260727}
    scene = {
        "scenario_id": "A100",
        "n_rx": 1,
        "plan_section": 15,
        "paired_policy": {
            "target_errors": 200,
            "min_total_trials": 10000,
            "max_total_trials": 50000,
        },
    }
    candidate = {
        "candidate_id": "A100_PRG_DFT8_6RB",
        "family": "TRANSPARENT_PRG_DFT",
        "prg_size_rb": 6,
        "mapping": "cycle_all",
    }
    path, _ = runner._write_transparent_manifest(
        config, scene, "source-sha256", [candidate]
    )
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["system"]["n_rx"] == 1
    assert payload["plan"].startswith("research/plan-028-comb6三类CSI-TDL-A300ns.md#15-")


def test_paired_tasks_accept_estimated_only_receiver() -> None:
    row = _base(candidate="A100_PRG_DFT8_6RB", trials=10000, errors=120)
    row["snr_db"] = 17.5
    tasks = runner.build_paired_tasks(
        {"estimated": [row]},
        {
            "target_errors": 200,
            "min_total_trials": 10000,
            "max_total_trials": 50000,
            "max_interval_trials": 1000,
        },
        batch_size=20,
    )
    assert tasks == [
        {
            "snr_db": 17.5,
            "trial_start": 10001,
            "trial_end": 11000,
            "trials": 1000,
        }
    ]


def test_run_snr_partition_filters_without_changing_trial_intervals() -> None:
    tasks = [
        {"snr_db": 16.25, "trial_start": 5001, "trial_end": 6000, "trials": 1000},
        {"snr_db": 17.0, "trial_start": 3001, "trial_end": 4000, "trials": 1000},
        {"snr_db": 17.75, "trial_start": 3001, "trial_end": 4000, "trials": 1000},
    ]
    selected = runner._filter_run_snr_tasks({"run_snr_db": [17.0, 17.75]}, tasks)
    assert selected == tasks[1:]


def test_paired_tasks_use_common_maximum_budget() -> None:
    rows_by_receiver = {}
    for receiver, error_pair in (("estimated", (250, 220)), ("ideal", (80, 50))):
        rows = []
        for candidate, errors in zip(("A100_PRG_DFT8_4RB", "A100_PRG_DFT8_6RB"), error_pair):
            row = _base(candidate=candidate, trials=10000, errors=errors)
            row["snr_db"] = 16.0
            rows.append(row)
        rows_by_receiver[receiver] = rows
    tasks = runner.build_paired_tasks(
        rows_by_receiver,
        {"target_errors": 200, "min_total_trials": 10000, "max_total_trials": 50000},
        batch_size=100,
    )
    assert tasks == [
        {
            "snr_db": 16.0,
            "trial_start": 10001,
            "trial_end": 40000,
            "trials": 30000,
        }
    ]


def test_paired_tasks_cap_each_resumable_interval() -> None:
    rows_by_receiver = {
        receiver: [
            {
                **_base(candidate=candidate),
                "snr_db": 14.0,
                "trials": 0,
                "tb_errors": 0,
                "bler": 0.0,
            }
            for candidate in ("A100_PRG_DFT8_4RB", "A100_PRG_DFT8_6RB")
        ]
        for receiver in ("estimated", "ideal")
    }
    tasks = runner.build_paired_tasks(
        rows_by_receiver,
        {
            "target_errors": 200,
            "min_total_trials": 10000,
            "max_total_trials": 50000,
            "max_interval_trials": 1000,
        },
        batch_size=20,
    )
    assert tasks == [
        {
            "snr_db": 14.0,
            "trial_start": 1,
            "trial_end": 1000,
            "trials": 1000,
        }
    ]


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


def test_mismatched_trace_nmse_matches_scalar_matched_closed_form() -> None:
    covariance = np.asarray([[1.0, 0.5], [0.5, 1.0]], dtype=np.complex128)
    noise_variance = 0.2
    estimator = runner.build_frequency_rmmse_filter(
        covariance,
        np.asarray([0]),
        noise_variance,
        diagonal_loading=0.0,
    )
    actual = runner.mismatched_frequency_lmmse_trace_nmse(
        covariance,
        estimator.weights,
        np.asarray([0]),
        np.asarray([1]),
        noise_variance,
    )
    expected = 1.0 - 0.5**2 / (1.0 + noise_variance)
    assert actual == pytest.approx(expected, abs=1e-12)


def test_lag_nmse_metric_matches_full_covariance_formula() -> None:
    base = np.asarray(
        [[1.0, 0.6 + 0.1j, 0.3], [0.6 - 0.1j, 1.0, 0.6 + 0.1j], [0.3, 0.6 - 0.1j, 1.0]],
        dtype=np.complex128,
    )
    pilots = np.asarray([0, 2], dtype=np.int64)
    data = np.asarray([1], dtype=np.int64)
    noise = 0.2
    assumed = 8.0 * base
    estimator = runner.build_frequency_rmmse_filter(
        assumed, pilots, noise, diagonal_loading=0.0
    )
    coordinates = np.asarray([-0.75, -0.05, 0.65, 1.35, 2.05, 2.75, 3.45, 4.15])
    indices = np.arange(3, dtype=np.float64)
    phase = np.exp(-2j * np.pi * indices[:, None] * coordinates[None, :] / 3.0)
    true_covariance = base * (phase @ phase.conj().T)
    expected = runner.mismatched_frequency_lmmse_trace_nmse(
        true_covariance, estimator.weights, pilots, data, noise
    )
    signed_lags = np.arange(-2, 3, dtype=np.float64)
    factor = np.sum(
        np.exp(-2j * np.pi * signed_lags[:, None] * coordinates[None, :] / 3.0),
        axis=1,
    )
    metric = nonuniform_search.build_lag_nmse_metric(
        base, estimator.weights, pilots, data, noise
    )
    actual = nonuniform_search.evaluate_lag_nmse(metric, factor)
    assert actual == pytest.approx(expected, abs=1e-12)


def test_transparent_cdd_candidate_preserves_source_transmitter() -> None:
    source = {
        "scenario_id": "A100",
        "candidate_id": "A100_SOURCE",
        "family": "SOURCE",
        "dmrs_spacing_subcarriers": 6,
        "delay_grid_coordinates": [0, 1, 2, 3, 4, 5, 6, 7],
    }
    manifest = {
        "scene_counts": {"A100": 1},
        "dmrs_spacing_subcarriers": 6,
        "candidates": [source],
    }
    scene = {
        "scenario_id": "A100",
        "transparent_cdd_baselines": [
            {
                "candidate_id": "A100_SOURCE_TRANSPARENT_CDD",
                "source_candidate_id": "A100_SOURCE",
                "label": "transparent",
            }
        ],
    }
    candidate = runner._transparent_cdd_candidates(scene, manifest)[0]
    assert candidate["candidate_id"] == "A100_SOURCE_TRANSPARENT_CDD"
    assert candidate["source_candidate_id"] == "A100_SOURCE"
    assert candidate["delay_grid_coordinates"] == source["delay_grid_coordinates"]
    assert candidate["source_family"] == "SOURCE"


def test_transparent_cdd_candidate_accepts_frozen_explicit_delays() -> None:
    manifest = {
        "scene_counts": {"A100": 0},
        "dmrs_spacing_subcarriers": 6,
        "candidates": [],
    }
    scene = {
        "scenario_id": "A100",
        "transparent_cdd_baselines": [
            {
                "candidate_id": "A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD",
                "delay_grid_coordinates": [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75],
            }
        ],
    }
    candidate = runner._transparent_cdd_candidates(scene, manifest)[0]
    assert candidate["source_candidate_id"] == ""
    assert candidate["transmitter_definition"] == "explicit_delay_grid_coordinates"
    assert candidate["delay_grid_coordinates"] == [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75]
    assert candidate["dmrs_spacing_subcarriers"] == 6


def test_transparent_cdd_candidate_rejects_nonmonotone_explicit_delays() -> None:
    manifest = {
        "scene_counts": {"A100": 0},
        "dmrs_spacing_subcarriers": 6,
        "candidates": [],
    }
    scene = {
        "scenario_id": "A100",
        "transparent_cdd_baselines": [
            {
                "candidate_id": "BAD",
                "delay_grid_coordinates": [0, 1, 0.5, 1.5, 2, 2.5, 3, 3.5],
            }
        ],
    }
    with pytest.raises(ValueError, match="nondecreasing"):
        runner._transparent_cdd_candidates(scene, manifest)


def test_transparent_cdd_receiver_filters_support_physical_and_matched() -> None:
    physical = np.eye(3, dtype=np.complex128)
    true = {
        "candidate": np.array(
            [[8.0, 1.0, 0.0], [1.0, 8.0, 2.0], [0.0, 2.0, 8.0]],
            dtype=np.complex128,
        )
    }
    pilots = np.array([0, 2], dtype=np.int64)
    transparent, transparent_labels = runner._transparent_cdd_receiver_filters(
        "physical", ["candidate"], true, physical, pilots, 0.25
    )
    matched, matched_labels = runner._transparent_cdd_receiver_filters(
        "matched_effective", ["candidate"], true, physical, pilots, 0.25
    )
    expected_transparent = runner.build_frequency_rmmse_filter(8.0 * physical, pilots, 0.25, 1e-10)
    expected_matched = runner.build_frequency_rmmse_filter(true["candidate"], pilots, 0.25, 1e-10)
    np.testing.assert_allclose(transparent["candidate"].weights, expected_transparent.weights)
    np.testing.assert_allclose(matched["candidate"].weights, expected_matched.weights)
    assert transparent_labels["candidate"].startswith("8 * physical")
    assert matched_labels["candidate"].startswith("candidate-specific")
    assert not np.allclose(transparent["candidate"].weights, matched["candidate"].weights)
