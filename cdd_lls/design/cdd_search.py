from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class CDDCandidate:
    candidate_id: str
    family: str
    delay_indices: tuple[int, ...]
    metadata: Dict[str, object] = field(default_factory=dict)


def canonical_delay_set(delays: Sequence[int], period: int = 576) -> tuple[int, ...]:
    """Canonicalize permutation and common circular-shift equivalent delay sets."""
    values = tuple(sorted(int(value) % int(period) for value in delays))
    if not values:
        raise ValueError("CDD delay set must not be empty.")
    representatives = []
    for anchor in values:
        representatives.append(tuple(sorted((value - anchor) % int(period) for value in values)))
    return min(representatives)


def residues_and_lifts(
    delays: Sequence[int],
    pilot_period: int = 24,
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    values = tuple(int(value) for value in delays)
    residues = tuple(value % int(pilot_period) for value in values)
    lifts = tuple((value - residue) // int(pilot_period) for value, residue in zip(values, residues))
    return residues, lifts


def arithmetic_delay_sets(
    n_tx: int = 8,
    period: int = 576,
) -> list[tuple[int, ...]]:
    unique = {
        canonical_delay_set([(branch * step) % int(period) for branch in range(int(n_tx))], period)
        for step in range(1, int(period))
    }
    return sorted(unique)


def continuous_arithmetic_delays_ns(
    spacing_ns: float,
    n_tx: int = 8,
) -> tuple[float, ...]:
    """T1 arithmetic baseline; deliberately accepts no PDP/covariance input."""
    if float(spacing_ns) < 0.0:
        raise ValueError("Arithmetic delay spacing must be non-negative.")
    return tuple(float(branch) * float(spacing_ns) for branch in range(int(n_tx)))


def thick_sidon_requirements(
    support_ns: float,
    grid_resolution_ns: float,
    period: int = 576,
    n_tx: int = 8,
    pilot_period: int = 24,
) -> dict[str, int | float | bool]:
    """Geometry-only T1 thick-Sidon thresholds and packing necessary condition."""
    if float(support_ns) < 0.0 or float(grid_resolution_ns) <= 0.0:
        raise ValueError("Support must be non-negative and grid resolution positive.")
    pair_count = int(n_tx) * (int(n_tx) + 1) // 2
    required_pair = int(np.floor(2.0 * float(support_ns) / float(grid_resolution_ns)) + 1)
    required_fold = int(np.floor(float(support_ns) / float(grid_resolution_ns)) + 1)
    upper_pair = int(period) // pair_count
    upper_fold = int(pilot_period) // int(n_tx)
    return {
        "support_ns": float(support_ns),
        "required_pair_gap_indices": required_pair,
        "required_fold_gap_indices": required_fold,
        "pair_gap_packing_upper_bound": upper_pair,
        "fold_gap_packing_upper_bound": upper_fold,
        "hard_feasible_by_packing_bound": bool(
            required_pair <= upper_pair and required_fold <= upper_fold
        ),
    }


def unique_random_delay_sets(
    seed: int,
    count: int,
    initial: Iterable[Sequence[int]] = (),
    n_tx: int = 8,
    period: int = 576,
    pilot_period: int = 24,
) -> list[tuple[int, ...]]:
    """Deterministic unique-residue candidate states with canonical replay."""
    requested = int(count)
    if requested <= 0:
        return []
    rng = np.random.default_rng(int(seed))
    values: dict[tuple[int, ...], None] = {}
    for delays in initial:
        key = canonical_delay_set(delays, period)
        residues, _ = residues_and_lifts(key, pilot_period)
        if len(key) == int(n_tx) and len(set(residues)) == int(n_tx):
            values[key] = None
    while len(values) < requested:
        candidate = random_unique_residue_delays(
            rng,
            n_tx=n_tx,
            period=period,
            pilot_period=pilot_period,
        )
        values[canonical_delay_set(candidate, period)] = None
    return list(values)[:requested]


def coordinate_exchange(
    initial: Sequence[int],
    neighbors: Callable[[tuple[int, ...], int], Iterable[Sequence[int]]],
    objective: Callable[[tuple[int, ...]], float],
    max_sweeps: int = 20,
    no_improvement_sweeps: int = 2,
) -> tuple[tuple[int, ...], list[dict[str, object]]]:
    """Deterministic coordinate exchange over a caller-defined neighborhood."""
    current = tuple(int(value) for value in initial)
    current_score = float(objective(current))
    trace: list[dict[str, object]] = [{
        "sweep": 0,
        "coordinate": -1,
        "score": current_score,
        "accepted": True,
        "delay_indices": list(current),
    }]
    stagnant = 0
    for sweep in range(1, int(max_sweeps) + 1):
        improved = False
        for coordinate in range(len(current)):
            best = current
            best_score = current_score
            for proposal_values in neighbors(current, coordinate):
                proposal = tuple(int(value) for value in proposal_values)
                score = float(objective(proposal))
                if score < best_score - 1e-12 or (
                    abs(score - best_score) <= 1e-12 and proposal < best
                ):
                    best = proposal
                    best_score = score
            accepted = best != current
            if accepted:
                current = best
                current_score = best_score
                improved = True
            trace.append({
                "sweep": sweep,
                "coordinate": coordinate,
                "score": current_score,
                "accepted": accepted,
                "delay_indices": list(current),
            })
        stagnant = 0 if improved else stagnant + 1
        if stagnant >= int(no_improvement_sweeps):
            break
    return current, trace


def random_unique_residue_delays(
    rng: np.random.Generator,
    n_tx: int = 8,
    period: int = 576,
    pilot_period: int = 24,
) -> tuple[int, ...]:
    residues = rng.choice(int(pilot_period), size=int(n_tx), replace=False)
    lifts = rng.integers(0, int(period) // int(pilot_period), size=int(n_tx))
    return tuple(int(residue + int(pilot_period) * lift) for residue, lift in zip(residues, lifts))
