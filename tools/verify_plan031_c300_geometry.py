"""Exhaustively verify the four-branch GEO candidates in plan-031 C300."""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs" / "experiment031_pdcch_cdd" / "20260907_c300_4tx_2sym"
EXPECTED = {
    1: {"balance": [0, 2, 6, 14], "pair": [0, 2, 6, 14], "fold": [0, 2, 5, 16]},
    2: {"balance": [0, 4, 12, 33], "pair": [0, 5, 16, 49], "fold": [0, 4, 10, 32]},
    4: {"balance": [0, 10, 30, 97], "pair": [0, 11, 44, 66], "fold": [0, 9, 18, 27]},
}
THRESHOLDS = {1: (5, 3), 2: (9, 5), 4: (17, 9)}
PAIR_I, PAIR_J = np.triu_indices(4)


def circular_minimum(values: np.ndarray, modulus: int) -> np.ndarray:
    ordered = np.sort(np.mod(values, modulus), axis=1)
    gaps = np.diff(ordered, axis=1)
    wrap = ordered[:, :1] + modulus - ordered[:, -1:]
    return np.min(np.concatenate((gaps, wrap), axis=1), axis=1)


def lexicographically_smaller(left: tuple[int, ...] | None, right: np.ndarray) -> bool:
    return left is None or tuple(int(value) for value in right) < left


def search(al: int, batch_size: int = 100_000) -> dict[str, object]:
    k_active = 36 * al
    n_p = k_active // 4
    g_pair, g_fold = THRESHOLDS[al]
    best: dict[str, dict[str, object]] = {
        name: {"coordinates": None, "pair_gap": -1, "fold_gap": -1, "score": -1.0}
        for name in ("balance", "pair", "fold")
    }
    enumerated = 0
    feasible = 0
    iterator = itertools.combinations(range(1, k_active), 3)
    while True:
        chunk = list(itertools.islice(iterator, batch_size))
        if not chunk:
            break
        tail = np.asarray(chunk, dtype=np.int16)
        coords = np.column_stack((np.zeros(len(tail), dtype=np.int16), tail))
        enumerated += len(coords)
        residues = np.mod(coords, n_p)
        valid = np.asarray([len(set(row.tolist())) == 4 for row in residues], dtype=bool)
        coords = coords[valid]
        residues = residues[valid]
        feasible += len(coords)
        pair_sums = coords[:, PAIR_I] + coords[:, PAIR_J]
        pair_gap = circular_minimum(pair_sums, k_active)
        fold_gap = circular_minimum(residues, n_p)
        score = np.minimum(pair_gap / g_pair, fold_gap / g_fold)

        # Pair/fold representatives: primary gap, other gap, then lexicographic order.
        for name, primary, secondary in (
            ("pair", pair_gap, fold_gap),
            ("fold", fold_gap, pair_gap),
        ):
            current = best[name]
            maximum = int(np.max(primary))
            indices = np.flatnonzero(primary == maximum)
            other_max = int(np.max(secondary[indices]))
            indices = indices[secondary[indices] == other_max]
            row = coords[int(indices[0])]
            candidate = (maximum, other_max)
            incumbent = (int(current[name + "_gap"]), int(current["fold_gap" if name == "pair" else "pair_gap"]))
            if candidate > incumbent or (candidate == incumbent and lexicographically_smaller(current["coordinates"], row)):
                current.update(
                    coordinates=tuple(int(value) for value in row),
                    pair_gap=int(pair_gap[int(indices[0])]),
                    fold_gap=int(fold_gap[int(indices[0])]),
                    score=float(score[int(indices[0])]),
                )

        # Balance representative: score, larger remaining normalized gap, then lexicographic.
        current = best["balance"]
        maximum = float(np.max(score))
        indices = np.flatnonzero(np.isclose(score, maximum, rtol=0.0, atol=1e-15))
        other = np.maximum(pair_gap[indices] / g_pair, fold_gap[indices] / g_fold)
        other_max = float(np.max(other))
        indices = indices[np.isclose(other, other_max, rtol=0.0, atol=1e-15)]
        row = coords[int(indices[0])]
        candidate_key = (maximum, other_max)
        incumbent_other = max(float(current["pair_gap"]) / g_pair, float(current["fold_gap"]) / g_fold)
        incumbent_key = (float(current["score"]), incumbent_other)
        if candidate_key > incumbent_key or (
            candidate_key == incumbent_key and lexicographically_smaller(current["coordinates"], row)
        ):
            current.update(
                coordinates=tuple(int(value) for value in row),
                pair_gap=int(pair_gap[int(indices[0])]),
                fold_gap=int(fold_gap[int(indices[0])]),
                score=float(score[int(indices[0])]),
            )

    serialized = {}
    for name, row in best.items():
        coordinates = list(row["coordinates"])
        if coordinates != EXPECTED[al][name]:
            raise AssertionError(f"AL{al} {name}: got {coordinates}, expected {EXPECTED[al][name]}")
        serialized[name] = {**row, "coordinates": coordinates}
    return {
        "aggregation_level": al,
        "k_active": k_active,
        "n_p": n_p,
        "g_pair": g_pair,
        "g_fold": g_fold,
        "enumerated_canonical_coordinate_sets": enumerated,
        "fold_distinct_coordinate_sets": feasible,
        "pair_sum_definition": "10 unordered sums with self, 0 <= i <= j < 4, modulo K",
        "gap_definition": "minimum circular spacing after modular reduction",
        "tie_break": "primary objective, larger other objective, lexicographically smaller coordinates",
        "representatives": serialized,
    }


def main() -> None:
    receipts = [search(al) for al in (1, 2, 4)]
    path = OUTPUT / "geo_exhaustive_receipt.json"
    path.write_text(json.dumps(receipts, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(receipts, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
