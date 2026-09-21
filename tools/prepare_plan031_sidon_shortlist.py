"""Prepare a geometry-only strict-Sidon shortlist for plan-031.

This tool performs no channel or link simulation.  Four-transmit-branch spaces are
enumerated exactly.  Eight-transmit-branch pools are sampled deterministically.
The final representatives cover geometry extrema before deterministic farthest-
point filling; no PDP, NMSE, CE, or BLER value is used.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.design.cdd_search import canonical_delay_set


DEFAULT_OUTPUT = (
    ROOT
    / "outputs"
    / "experiment031_pdcch_cdd"
    / "20260908_strict_sidon_search"
    / "sidon_shortlist.json"
)
SCENARIOS = (
    ("C300", 4, 1, 36, 9),
    ("C300", 4, 2, 72, 18),
    ("C300", 4, 4, 144, 36),
    ("A100", 8, 2, 144, 36),
    ("A100", 8, 4, 288, 72),
    ("A100", 8, 8, 576, 144),
)
S0 = {
    4: (0, 1, 3, 7),
    8: (0, 1, 3, 7, 12, 20, 30, 65),
}


def circular_gaps(values: tuple[int, ...], period: int) -> tuple[int, ...]:
    ordered = tuple(sorted(int(value) % int(period) for value in values))
    return tuple(
        sorted(
            [
                ordered[index + 1] - ordered[index]
                for index in range(len(ordered) - 1)
            ]
            + [ordered[0] + int(period) - ordered[-1]]
        )
    )


def pair_sums(delays: tuple[int, ...], period: int) -> tuple[int, ...]:
    return tuple(
        sorted(
            (delays[left] + delays[right]) % int(period)
            for left in range(len(delays))
            for right in range(left, len(delays))
        )
    )


def is_valid(delays: tuple[int, ...], period: int, pilot_period: int) -> bool:
    if len(delays) != len(set(delays)):
        return False
    if len({value % int(pilot_period) for value in delays}) != len(delays):
        return False
    sums = pair_sums(delays, period)
    return len(sums) == len(set(sums))


def geometry(delays: tuple[int, ...], period: int, pilot_period: int) -> dict[str, object]:
    sums = pair_sums(delays, period)
    pair_gap_profile = circular_gaps(sums, period)
    fold_gap_profile = circular_gaps(
        tuple(value % int(pilot_period) for value in delays), pilot_period
    )
    delay_gap_profile = circular_gaps(delays, period)
    lags = np.arange(1, int(pilot_period), dtype=np.float64)
    values = np.asarray(delays, dtype=np.float64)
    factor = np.mean(
        np.exp(-2j * np.pi * lags[:, None] * values[None, :] / float(period)),
        axis=1,
    )
    power = np.abs(factor) ** 2
    pilot_rows = np.arange(int(pilot_period), dtype=np.float64)
    pilot_matrix = np.exp(
        -2j
        * np.pi
        * pilot_rows[:, None]
        * values[None, :]
        / float(pilot_period)
    )
    singular_values = np.linalg.svd(pilot_matrix, compute_uv=False)
    return {
        "pair_gap": int(pair_gap_profile[0]),
        "fold_gap": int(fold_gap_profile[0]),
        "circular_span": int(period - delay_gap_profile[-1]),
        "low_lag_energy": float(np.mean(power)),
        "low_lag_peak": float(np.max(power)),
        "pilot_rank": int(np.linalg.matrix_rank(pilot_matrix)),
        "pilot_condition": float(singular_values[0] / singular_values[-1]),
        "delay_gap_profile": [float(value) / float(period) for value in delay_gap_profile],
        "pair_gap_profile": [float(value) / float(period) for value in pair_gap_profile],
        "fold_gap_profile": [
            float(value) / float(pilot_period) for value in fold_gap_profile
        ],
    }


def enumerate_four_tx(period: int, pilot_period: int) -> tuple[list[tuple[int, ...]], int]:
    valid_with_zero = 0
    representatives: set[tuple[int, ...]] = set()
    for tail in itertools.combinations(range(1, int(period)), 3):
        candidate = (0, *tail)
        if not is_valid(candidate, period, pilot_period):
            continue
        valid_with_zero += 1
        representatives.add(canonical_delay_set(candidate, period))
    return sorted(representatives), valid_with_zero


def sample_eight_tx(
    period: int,
    pilot_period: int,
    seed: int,
    target: int,
    proposal_limit: int,
) -> tuple[list[tuple[int, ...]], int]:
    rng = np.random.default_rng(int(seed))
    representatives: set[tuple[int, ...]] = set()
    s0 = canonical_delay_set(S0[8], period)
    if is_valid(s0, period, pilot_period):
        representatives.add(s0)
    proposals = 0
    batch_size = 20_000
    pair_i, pair_j = np.triu_indices(8)
    while len(representatives) < int(target) and proposals < int(proposal_limit):
        count = min(batch_size, int(proposal_limit) - proposals)
        proposals += count
        batch = np.sort(
            rng.integers(0, int(period), size=(count, 8), dtype=np.int32), axis=1
        )
        unique_values = np.all(np.diff(batch, axis=1) != 0, axis=1)
        residues = np.sort(np.mod(batch, int(pilot_period)), axis=1)
        unique_residues = np.all(np.diff(residues, axis=1) != 0, axis=1)
        sums = np.sort(
            np.mod(batch[:, pair_i] + batch[:, pair_j], int(period)), axis=1
        )
        unique_sums = np.all(np.diff(sums, axis=1) != 0, axis=1)
        for row in batch[unique_values & unique_residues & unique_sums]:
            representatives.add(canonical_delay_set(row.tolist(), period))
            if len(representatives) >= int(target):
                break
    return sorted(representatives), proposals


def feature_vector(metric: dict[str, object], period: int, pilot_period: int) -> np.ndarray:
    scalars = [
        float(metric["pair_gap"]) / float(period),
        float(metric["fold_gap"]) / float(pilot_period),
        float(metric["circular_span"]) / float(period),
        float(metric["low_lag_energy"]),
        float(metric["low_lag_peak"]),
    ]
    return np.asarray(
        scalars
        + list(metric["delay_gap_profile"])
        + list(metric["pair_gap_profile"])
        + list(metric["fold_gap_profile"]),
        dtype=np.float64,
    )


def choose_shortlist(
    candidates: list[tuple[int, ...]],
    period: int,
    pilot_period: int,
    n_tx: int,
    count: int,
) -> tuple[list[tuple[int, ...]], dict[tuple[int, ...], dict[str, object]]]:
    metrics = {
        candidate: geometry(candidate, period, pilot_period) for candidate in candidates
    }
    selected: list[tuple[int, ...]] = []

    def add(candidate: tuple[int, ...]) -> None:
        if candidate not in selected and len(selected) < int(count):
            selected.append(candidate)

    add(canonical_delay_set(S0[n_tx], period))
    add(max(candidates, key=lambda row: (metrics[row]["pair_gap"], metrics[row]["fold_gap"], tuple(-x for x in row))))
    add(max(candidates, key=lambda row: (metrics[row]["fold_gap"], metrics[row]["pair_gap"], tuple(-x for x in row))))
    add(min(candidates, key=lambda row: (metrics[row]["circular_span"], row)))
    add(max(candidates, key=lambda row: (metrics[row]["circular_span"], tuple(-x for x in row))))
    add(min(candidates, key=lambda row: (metrics[row]["low_lag_energy"], row)))
    add(max(candidates, key=lambda row: (metrics[row]["low_lag_energy"], tuple(-x for x in row))))

    matrix = np.vstack(
        [feature_vector(metrics[candidate], period, pilot_period) for candidate in candidates]
    )
    minima = np.min(matrix, axis=0)
    ranges = np.max(matrix, axis=0) - minima
    ranges[ranges == 0.0] = 1.0
    matrix = (matrix - minima) / ranges
    index = {candidate: position for position, candidate in enumerate(candidates)}
    while len(selected) < min(int(count), len(candidates)):
        selected_matrix = matrix[[index[candidate] for candidate in selected]]
        distances = np.min(
            np.sum((matrix[:, None, :] - selected_matrix[None, :, :]) ** 2, axis=2),
            axis=1,
        )
        for candidate in selected:
            distances[index[candidate]] = -1.0
        best_distance = float(np.max(distances))
        tied = np.flatnonzero(np.isclose(distances, best_distance, rtol=0.0, atol=1e-15))
        add(min(candidates[int(position)] for position in tied))
    return selected, metrics


def prepare(args: argparse.Namespace) -> dict[str, object]:
    rows: list[dict[str, object]] = []
    receipts: list[dict[str, object]] = []
    for family, n_tx, al, period, pilot_period in SCENARIOS:
        if n_tx == 4:
            candidates, examined = enumerate_four_tx(period, pilot_period)
            generation = "complete_enumeration"
            proposals = math.comb(period - 1, 3)
            valid_with_zero = examined
        else:
            seed = int(args.seed) + 10_000 * n_tx + 100 * al + period
            candidates, proposals = sample_eight_tx(
                period,
                pilot_period,
                seed,
                int(args.eight_tx_pool),
                int(args.eight_tx_proposal_limit),
            )
            generation = "deterministic_random_pool"
            valid_with_zero = None
        shortlist, metrics = choose_shortlist(
            candidates, period, pilot_period, n_tx, int(args.shortlist_size)
        )
        if len(shortlist) != len(set(shortlist)):
            raise AssertionError("Shortlist contains duplicate canonical representatives.")
        for rank, candidate in enumerate(shortlist, start=1):
            if not is_valid(candidate, period, pilot_period):
                raise AssertionError(f"Invalid shortlist candidate: {candidate}")
            if int(metrics[candidate]["pilot_rank"]) != int(n_tx):
                raise AssertionError(f"Rank-deficient shortlist candidate: {candidate}")
            rows.append(
                {
                    "family": family,
                    "n_tx": n_tx,
                    "aggregation_level": al,
                    "k_active": period,
                    "pilot_period": pilot_period,
                    "candidate_index": rank,
                    "delay_indices": list(candidate),
                    "delay_ns": [
                        float(value) * 1e9 / (float(period) * 30e3)
                        for value in candidate
                    ],
                    **{
                        key: value
                        for key, value in metrics[candidate].items()
                        if not key.endswith("_profile")
                    },
                }
            )
        receipts.append(
            {
                "family": family,
                "n_tx": n_tx,
                "aggregation_level": al,
                "k_active": period,
                "pilot_period": pilot_period,
                "generation": generation,
                "proposals_or_zero_anchored_combinations": proposals,
                "valid_zero_anchored_before_canonicalization": valid_with_zero,
                "unique_valid_canonical_pool": len(candidates),
                "shortlist_size": len(shortlist),
            }
        )
    return {
        "seed": int(args.seed),
        "eight_tx_pool_target": int(args.eight_tx_pool),
        "eight_tx_proposal_limit": int(args.eight_tx_proposal_limit),
        "shortlist_size_per_scenario": int(args.shortlist_size),
        "selection_uses_channel_or_link_metrics": False,
        "exact_equivalence": "antenna permutation and common cyclic shift only",
        "receipts": receipts,
        "candidates": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20260908)
    parser.add_argument("--eight-tx-pool", type=int, default=20_000)
    parser.add_argument("--eight-tx-proposal-limit", type=int, default=20_000_000)
    parser.add_argument("--shortlist-size", type=int, default=8)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = prepare(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(result["receipts"], indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
