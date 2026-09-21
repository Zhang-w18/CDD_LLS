from __future__ import annotations

import numpy as np

from tools.analyze_plan031 import analyze, interpolate_target, isotonic_decreasing


def test_isotonic_decreasing_pools_adjacent_violation() -> None:
    fitted = isotonic_decreasing(
        np.asarray([0.2, 0.08, 0.10, 0.01]),
        np.asarray([100, 100, 100, 100]),
    )
    assert np.allclose(fitted, [0.2, 0.09, 0.09, 0.01])


def test_log_bler_target_interpolation_is_local() -> None:
    value, index = interpolate_target(
        np.asarray([0.0, 0.25, 0.5]),
        np.asarray([0.2, 0.1, 0.01]),
        0.05,
    )
    assert index == 1
    assert 0.25 < value < 0.5


def test_analyze_marks_unbracketed_target_without_extrapolation() -> None:
    configs = {
        2: {
            "candidates": [
                {"candidate_id": "B0_QC"},
                {"candidate_id": "A100_PRG_DFT8_6RB"},
            ]
        }
    }
    points = []
    for candidate_id in ("B0_QC", "A100_PRG_DFT8_6RB"):
        points.extend(
            [
                {
                    "aggregation_level": 2,
                    "candidate_id": candidate_id,
                    "snr_db": 0.0,
                    "errors": 20,
                    "trials": 100,
                    "bler": 0.2,
                },
                {
                    "aggregation_level": 2,
                    "candidate_id": candidate_id,
                    "snr_db": 0.25,
                    "errors": 15,
                    "trials": 100,
                    "bler": 0.15,
                },
            ]
        )

    targets, gains = analyze(configs, points, bootstrap_repeats=20)

    assert all(row["status"] == "unbracketed" for row in targets)
    assert all(row["target_snr_db"] is None for row in targets)
    assert all(row["status"] == "unavailable_unbracketed_target" for row in gains)
    assert all(row["gain_db"] is None for row in gains)
