from __future__ import annotations

from pathlib import Path

import yaml

from tools.prepare_plan038_formal import TRANSPARENT_TO_MATCHED, _grid_for
from tools.run_plan038_refinement import WINDOWS, _quarter_db_grid


ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {
    "SIDON_MATCHED": "matched_effective",
    "SIDON_TRANSPARENT": "physical_fullband",
    "B0QC_MATCHED": "matched_effective",
    "B0QC_TRANSPARENT": "physical_fullband",
    "CDD911_MATCHED": "matched_effective",
    "CDD911_TRANSPARENT": "physical_fullband",
    "CDD130_MATCHED": "matched_effective",
    "CDD130_TRANSPARENT": "physical_fullband",
    "PRG_DFT4_TRANSPARENT": "physical_prg",
}


def test_plan038_prescan_configs_are_frozen_consistently() -> None:
    expected_order = {1: [0], 2: [0, 1], 4: [0, 1, 2, 3]}
    for al in (1, 2, 4):
        path = ROOT / f"configs/pdcch_result038_al{al}_prescan.yaml"
        config = yaml.safe_load(path.read_text(encoding="utf-8"))
        assert config["antenna"] == {"n_tx": 4, "n_rx": 2}
        assert config["resource"]["duration_symbols"] == 2
        assert config["resource"]["aggregation_level"] == al
        assert config["pdcch"]["payload_bits"] == 40
        assert config["pdcch"]["coded_bits"] == 108 * al
        candidates = {item["candidate_id"]: item for item in config["candidates"]}
        assert {key: item["receiver_covariance_mode"] for key, item in candidates.items()} == EXPECTED
        assert candidates["PRG_DFT4_TRANSPARENT"]["cycling_order"] == expected_order[al]
        for stem in ("SIDON", "B0QC", "CDD911", "CDD130"):
            matched = candidates[f"{stem}_MATCHED"]
            transparent = candidates[f"{stem}_TRANSPARENT"]
            assert matched["delay_grid_coordinates"] == transparent["delay_grid_coordinates"]


def test_plan038_formal_grid_uses_quarter_db_and_brackets_one_percent() -> None:
    grid, receipt = _grid_for([(0.0, 0.4), (2.0, 0.08), (4.0, 0.008)])
    assert all(abs(value * 4 - round(value * 4)) < 1e-12 for value in grid)
    assert receipt["bler_0.01"]["status"] == "bracketed"
    assert grid[0] == 2.0 and grid[-1] == 4.0


def test_plan038_no_bracket_uses_declared_diagnostic_window() -> None:
    grid, receipt = _grid_for([(0.0, 0.9), (2.0, 0.8), (4.0, 0.7)])
    assert receipt["bler_0.01"]["status"] == "nearest-window-no-prescan-bracket"
    assert grid[-1] == 4.5


def test_plan038_transparent_cdd_grid_sources_are_explicit() -> None:
    assert TRANSPARENT_TO_MATCHED == {
        "SIDON_TRANSPARENT": "SIDON_MATCHED",
        "B0QC_TRANSPARENT": "B0QC_MATCHED",
        "CDD911_TRANSPARENT": "CDD911_MATCHED",
        "CDD130_TRANSPARENT": "CDD130_MATCHED",
    }


def test_plan038_refinement_has_all_curves_and_uniform_quarter_db_grids() -> None:
    assert len(WINDOWS) == 27
    for (al, candidate_id), (start, stop) in WINDOWS.items():
        assert al in {1, 2, 4}
        assert candidate_id in EXPECTED
        grid = _quarter_db_grid(start, stop)
        assert grid[0] == start
        assert grid[-1] == stop
        assert all(abs((right - left) - 0.25) < 1e-12 for left, right in zip(grid, grid[1:]))
    assert WINDOWS[(1, "SIDON_TRANSPARENT")] == (4.0, 5.0)
    assert WINDOWS[(4, "CDD130_MATCHED")] == (-2.5, -1.0)
    assert WINDOWS[(4, "CDD130_TRANSPARENT")] == (-2.5, -1.0)
    assert WINDOWS[(2, "CDD911_TRANSPARENT")] == (2.0, 3.0)
    assert WINDOWS[(2, "SIDON_TRANSPARENT")] == (4.0, 5.0)
    assert WINDOWS[(4, "B0QC_TRANSPARENT")] == (-1.0, 0.0)
    assert WINDOWS[(4, "SIDON_TRANSPARENT")] == (-2.0, -1.0)
