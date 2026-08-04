import csv
import json
from pathlib import Path

import numpy as np

from tools import run_plan027_dense_dmrs as dense027
from tools import run_plan028_csi_curves as plan028


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _write_point(root: Path, stage: str, trials: int, bler: float) -> None:
    row = {
        "scenario_id": "A30",
        "candidate_id": "A30_B0_QC",
        "family": "B0_QC",
        "snr_db": 14.0,
        "trials": trials,
        "tb_errors": int(round(trials * bler)),
        "bler": bler,
        "ce_nmse_mean_db": -20.0,
    }
    stage_root = root / stage
    _write_csv(stage_root / "bler_points.csv", [row])
    flag = plan028._flag_path(stage_root, 14.0, "A30_B0_QC")
    flag.parent.mkdir(parents=True, exist_ok=True)
    np.save(flag, np.zeros(trials, dtype=bool))


def test_a300_uses_six_baselines_and_dynamic_support() -> None:
    dense027.configure_scenario("A300", 6)
    assert dense027.SPREAD_NS == 300.0
    assert dense027.active_baseline_families() == (
        "B0_QC",
        "AP_RMS_T1",
        "AP_TEPS_T1",
        "AP_TU_NT",
        "AP_TU_NTM1",
        "AP_TALIAS_NT",
    )
    delays = dense027.useful_symbol_nt_minus_one_delays_ns()
    assert len(delays) == 8
    assert np.isclose(delays[-1], dense027.USEFUL_SYMBOL_NS)


def test_curve_merge_prefers_trials_then_refine_priority(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(plan028, "ROOT", tmp_path)
    _write_point(tmp_path, "prescan", 400, 0.2)
    _write_point(tmp_path, "refine_10pct", 3000, 0.1)
    _write_point(tmp_path, "refine_1pct", 3000, 0.01)
    rows = plan028._canonical_curve_rows(tmp_path, "A30")
    assert len(rows) == 1
    assert rows[0]["source_stage"] == "refine_1pct"
    assert int(rows[0]["trials"]) == 3000


def test_ideal_equalizer_uses_exact_channel() -> None:
    true_data = np.asarray([[[2.0 + 0.0j, 0.0 + 1.0j]]])
    symbols = np.asarray([[1.0 + 1.0j, -1.0 + 0.5j]])
    received = true_data * symbols[:, None, :]
    equalized, effective_noise = plan028.ideal_csi_equalize(
        true_data, received, 0.25
    )
    np.testing.assert_allclose(equalized, symbols)
    np.testing.assert_allclose(effective_noise, np.asarray([[0.0625, 0.25]]))


def test_plot_readability_constants_meet_result_spec() -> None:
    assert plan028.PLOT_FONT >= 16
    assert plan028.TICK_FONT >= 14
    assert plan028.LINE_WIDTH >= 2.5
    assert plan028.MARKER_SIZE >= 8


def test_candidate_numbering_is_bijective(tmp_path) -> None:
    manifest = {
        "scene_counts": {"A300": 2},
        "dmrs_spacing_subcarriers": 6,
        "candidates": [
            {
                "scenario_id": "A300",
                "candidate_id": "A300_S0_SIDON",
                "family": "S0_SIDON",
                "dmrs_spacing_subcarriers": 6,
                "delay_grid_coordinates": [0, 1, 3, 7, 12, 20, 30, 65],
                "delay_ns": [0, 1, 3, 7, 12, 20, 30, 65],
                "pilot_rank": 8,
                "pilot_condition_number": 1.0,
            },
            {
                "scenario_id": "A300",
                "candidate_id": "A300_B0_QC",
                "family": "B0_QC",
                "dmrs_spacing_subcarriers": 6,
                "delay_grid_coordinates": [0, 9, 18, 27, 36, 45, 54, 63],
                "delay_ns": [0, 9, 18, 27, 36, 45, 54, 63],
                "pilot_rank": 8,
                "pilot_condition_number": 1.0,
            },
        ],
    }
    plan028.write_candidate_assets(tmp_path, "A300", manifest)
    rows = list(
        csv.DictReader(
            (tmp_path / "a300_comb6/final/candidate_number_delay_table.csv").open(
                "r", encoding="utf-8", newline=""
            )
        )
    )
    assert [row["candidate_number"] for row in rows] == [
        "candidate 01",
        "candidate 02",
    ]
    assert [row["candidate_id"] for row in rows] == ["A300_B0_QC", "A300_S0_SIDON"]
    styles = json.loads(
        (tmp_path / "a300_comb6/final/curve_styles.json").read_text(encoding="utf-8")
    )
    assert set(styles) == {row["candidate_id"] for row in rows}
