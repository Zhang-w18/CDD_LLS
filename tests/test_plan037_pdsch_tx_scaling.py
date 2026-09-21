from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from cdd_lls.phy.precoding import build_prg_dft_precoder
from tools.run_plan037_pdsch_tx_scaling import (
    FROZEN_DELAYS,
    audit_delay_set,
    candidate_definitions,
    geometry_audit,
    load_config,
    prg_vector_indices,
    validate_config,
)


ROOT = Path(__file__).resolve().parents[1]


def test_frozen_delay_geometry() -> None:
    audit4 = audit_delay_set(4, "S0_SIDON", FROZEN_DELAYS[4]["S0_SIDON"])
    assert audit4["strict_sidon"] is True
    assert audit4["minimum_fold_gap"] == 24
    assert audit4["pilot_rank"] == 4
    assert np.isclose(audit4["pilot_condition_number"], 1.0)

    audit16 = audit_delay_set(16, "S0_SIDON", FROZEN_DELAYS[16]["S0_SIDON"])
    assert audit16["unique_pair_sum_count"] == 136
    assert len(set(audit16["comb6_residues_mod_96"])) == 16
    assert audit16["minimum_fold_gap"] == 4
    assert audit16["pilot_rank"] == 16
    assert np.isclose(audit16["pilot_condition_number"], 1.0)


def test_grouped_b0_rank_and_repetition() -> None:
    for n_tx, repeats in ((16, 2), (32, 4)):
        values = FROZEN_DELAYS[n_tx]["B0_QC_GROUPED"]
        assert values == [0, 9, 18, 27, 36, 45, 54, 63] * repeats
        audit = audit_delay_set(n_tx, "B0_QC_GROUPED", values)
        assert audit["pilot_rank"] == 8
        assert audit["pilot_condition_number"] == "inf"


def test_32tx_sidon_proof_and_run_rejection() -> None:
    proof = geometry_audit()["sidon_32tx"]
    assert proof["required_distinct_nonzero_directed_differences"] == 992
    assert proof["available_nonzero_residues_mod_576"] == 575
    assert proof["inequality_holds"] is False
    with pytest.raises(ValueError, match="strict Sidon is infeasible"):
        candidate_definitions(32)


def test_prg_mapping_and_power() -> None:
    assert prg_vector_indices(4) == [0, 1, 2, 3, 0, 1, 2, 3]
    assert prg_vector_indices(16) == list(range(8))
    for n_tx in (4, 16):
        config = load_config(ROOT / "configs" / f"plan037_nt{n_tx}_nr1_v60_smoke.yaml")
        candidates, _, grid, _, _ = validate_config(config)
        prg_row = next(row for row in candidates if row["family"] == "TRANSPARENT_PRG_DFT")
        precoder = build_prg_dft_precoder(grid, n_tx, 6, prg_row["prg_vector_indices"], normalize=True)
        np.testing.assert_allclose(np.sum(np.abs(precoder.C) ** 2, axis=1), 1.0, atol=1e-12)
        for prg, vector in enumerate(prg_row["prg_vector_indices"]):
            block = precoder.C[prg * 72:(prg + 1) * 72]
            np.testing.assert_allclose(block, np.repeat(block[0][None, :], 72, axis=0), atol=1e-12)
            assert vector < n_tx


def test_smoke_configs_have_three_candidates_and_six_curves() -> None:
    for n_tx in (4, 16):
        config = load_config(ROOT / "configs" / f"plan037_nt{n_tx}_nr1_v60_smoke.yaml")
        candidates, *_ = validate_config(config)
        assert len(candidates) == 3
        assert len(candidates) * len(config["receiver_modes"]) == 6
