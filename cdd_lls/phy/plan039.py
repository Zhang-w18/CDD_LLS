"""Deterministic geometry and frozen-link helpers for Plan-039."""

from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from cdd_lls.phy.cdl_beam_platform import (
    beam_domain_cdd_precoder,
    covariance_beam_specific_independent,
    covariance_common_reference_pdp,
)

CODEBOOK_TYPE = "angular_full_coverage_ultrawide"
CE_METHODS = {
    "PLAN039_COMMON_REFERENCE_PDP",
    "PLAN039_BEAM_SPECIFIC_PDP_INDEPENDENT",
    "PLAN039_TRANSPARENT_COMMON_REFERENCE_PDP",
    "PLAN039_PRG_COMMON_REFERENCE_PDP",
}


def _wide_factor(
    element_count: int,
    spatial_frequencies: np.ndarray,
    regularizations: Sequence[float],
) -> tuple[np.ndarray, dict[str, Any]]:
    """Choose a deterministic least-squares wide factor by frozen lexicographic score."""
    n = np.arange(int(element_count), dtype=np.float64)
    frequencies = np.asarray(spatial_frequencies, dtype=np.float64).reshape(-1)
    response = np.exp(-1j * np.pi * frequencies[:, None] * n[None, :])
    records: list[tuple[tuple[float, ...], np.ndarray, dict[str, Any]]] = []
    for regularization in regularizations:
        normal = response.conj().T @ response + float(regularization) * np.eye(element_count)
        weight = np.linalg.solve(normal, response.conj().T @ np.ones(frequencies.size, dtype=np.complex128))
        weight /= np.linalg.norm(weight)
        gain = np.abs(response @ weight) ** 2
        p5, p95 = np.percentile(gain, [5.0, 95.0])
        metrics = {
            "regularization": float(regularization),
            "minimum_gain_linear": float(np.min(gain)),
            "mean_gain_linear": float(np.mean(gain)),
            "p5_gain_linear": float(p5),
            "p95_gain_linear": float(p95),
            "p95_p5_ripple_linear": float(p95 - p5),
        }
        score = (-metrics["minimum_gain_linear"], metrics["p95_p5_ripple_linear"],
                 -metrics["mean_gain_linear"], float(regularization))
        records.append((score, weight, metrics))
    _, selected, metrics = min(records, key=lambda item: item[0])
    metrics["candidate_count"] = len(records)
    metrics["selection_priority"] = ["max_min_gain", "min_p95_p5_ripple", "max_mean_gain", "min_regularization"]
    return selected, metrics


def build_angular_full_coverage_codebook(
    regularizations: Sequence[float] = (1e-8, 1e-6, 1e-4, 1e-2, 1.0),
    horizontal_samples: int = 241,
    vertical_samples: int = 81,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Build the frozen separable reference beam and eight shifted-DFT beams."""
    aod = np.deg2rad(np.linspace(-60.0, 60.0, int(horizontal_samples)))
    zod = np.deg2rad(np.linspace(90.0, 110.0, int(vertical_samples)))
    # Half-wavelength horizontal TXRU spacing gives pi*sin(AoD).  The two
    # vertical TXRUs are four 0.8-lambda AE spacings apart.
    horizontal_frequency = np.sin(aod)
    vertical_frequency = 2.0 * 3.2 * np.cos(zod)
    horizontal_wide, horizontal_metrics = _wide_factor(8, horizontal_frequency, regularizations)
    vertical_wide, vertical_metrics = _wide_factor(2, vertical_frequency, regularizations)
    polarization = np.asarray([1.0, 1.0], dtype=np.complex128) / np.sqrt(2.0)
    xi = (-7.0 + 2.0 * np.arange(8)) / 8.0
    horizontal_narrow = np.exp(1j * np.pi * np.arange(8)[:, None] * xi[None, :]) / np.sqrt(8.0)
    reference = np.kron(polarization, np.kron(vertical_wide, horizontal_wide))
    narrow = np.column_stack([
        np.kron(polarization, np.kron(vertical_wide, horizontal_narrow[:, index]))
        for index in range(8)
    ])
    gram_error = float(np.max(np.abs(narrow.conj().T @ narrow - np.eye(8))))
    if abs(np.linalg.norm(reference) - 1.0) > 1e-12 or gram_error > 1e-12:
        raise RuntimeError("Plan-039 codebook normalization or orthogonality failed.")
    centers = np.rad2deg(np.arcsin(np.clip(xi, -1.0, 1.0)))
    return reference, narrow, {
        "type": "ANGULAR_FULL_COVERAGE_ULTRAWIDE",
        "flatten_order": ["polarization", "vertical", "horizontal"],
        "horizontal_spatial_frequencies": xi.tolist(),
        "horizontal_center_aod_deg": centers.tolist(),
        "reference_norm_error": float(abs(np.linalg.norm(reference) - 1.0)),
        "narrow_gram_max_error": gram_error,
        "horizontal_wide": horizontal_metrics,
        "vertical_wide": vertical_metrics,
    }


def circular_min_gap(values: Sequence[int], modulus: int) -> int:
    """Return the minimum adjacent circular gap of distinct modular values."""
    ordered = np.unique(np.mod(np.asarray(values, dtype=np.int64), int(modulus)))
    if ordered.size < 2:
        return int(modulus)
    gaps = np.diff(np.r_[ordered, ordered[0] + int(modulus)])
    return int(np.min(gaps))


def unordered_pair_sums(indices: Sequence[int], modulus: int = 576) -> np.ndarray:
    values = tuple(int(value) % int(modulus) for value in indices)
    return np.asarray(
        [(values[left] + values[right]) % int(modulus)
         for left in range(len(values)) for right in range(left, len(values))],
        dtype=np.int64,
    )


def canonical_sidon(indices: Sequence[int], modulus: int = 576) -> tuple[int, ...]:
    """Canonicalize common shifts, reflection, and branch permutations."""
    values = np.mod(np.asarray(indices, dtype=np.int64), int(modulus))
    representatives: list[tuple[int, ...]] = []
    for sign in (1, -1):
        reflected = np.mod(sign * values, int(modulus))
        for origin in reflected:
            representatives.append(tuple(sorted(np.mod(reflected - origin, int(modulus)).tolist())))
    return min(representatives)


def sidon_candidate_audit(indices: Sequence[int], modulus: int = 576, pilot_period: int = 96) -> dict[str, Any]:
    values = canonical_sidon(indices, modulus)
    sums = unordered_pair_sums(values, modulus)
    residues = np.mod(values, int(pilot_period))
    pilot_index = np.arange(int(pilot_period), dtype=np.float64)
    pilot = np.exp(-1j * 2.0 * np.pi * pilot_index[:, None] * residues[None, :] / float(pilot_period))
    singular = np.linalg.svd(pilot, compute_uv=False)
    return {
        "delay_indices": list(values),
        "delay_grid_rationals": [f"{value}/{modulus}" for value in values],
        "delay_seconds": [value / (float(modulus) * 30e3) for value in values],
        "fft_sample_equivalents": [value * 4096.0 / float(modulus) for value in values],
        "comb6_residues": residues.tolist(),
        "unordered_pair_sums": sums.tolist(),
        "strict_sidon": bool(np.unique(sums).size == sums.size),
        "pair_gap": circular_min_gap(sums, modulus),
        "fold_gap": circular_min_gap(residues, pilot_period),
        "pilot_rank": int(np.linalg.matrix_rank(pilot, tol=1e-10)),
        "pilot_condition_number": float(singular[0] / singular[-1]),
    }


@lru_cache(maxsize=1)
def search_strict_sidon_top8(modulus: int = 576, pilot_period: int = 96, branches: int = 8) -> dict[str, Any]:
    """Infer the top strict-Sidon sets from the declared modular rules.

    For the maximum possible fold gap, eight residues on Z_96 must be exactly
    twelve apart.  The remaining finite search is therefore the six possible
    Z_576 lifts of each residue.  No channel realization or performance value
    enters this construction.
    """
    if (modulus, pilot_period, branches) != (576, 96, 8):
        raise ValueError("Plan-039 freezes K=576, Np=96, and eight branches.")
    theoretical_fold_upper_bound = pilot_period // branches

    def compositions(total: int, count: int, prefix: tuple[int, ...] = ()):
        if count == 1:
            yield (*prefix, total)
            return
        for value in range(total + 1):
            yield from compositions(total - value, count - 1, (*prefix, value))

    attempted: dict[int, dict[str, int]] = {}
    ranked: list[dict[str, Any]] = []
    for target_gap in range(theoretical_fold_upper_bound, 0, -1):
        extra = pilot_period - branches * target_gap
        unique: dict[tuple[int, ...], dict[str, Any]] = {}
        residue_sets = 0
        lift_nodes = 0
        pool_limit = 512
        stopped_at_pool_limit = False
        for additions in compositions(extra, branches):
            gaps = np.asarray(additions, dtype=np.int64) + target_gap
            residues = np.r_[0, np.cumsum(gaps[:-1])]
            # Common rotation/reflection equivalence: retain only the smallest
            # cyclic/reflected gap word.  This is a purely geometric reduction.
            gap_words = []
            for sign in (1, -1):
                word = gaps if sign == 1 else gaps[::-1]
                gap_words.extend(tuple(np.roll(word, shift).tolist()) for shift in range(branches))
            if tuple(gaps.tolist()) != min(gap_words):
                continue
            residue_sets += 1

            def extend(position: int, values: list[int], sums: set[int]) -> None:
                nonlocal lift_nodes, stopped_at_pool_limit
                if stopped_at_pool_limit:
                    return
                if position == branches:
                    audit = sidon_candidate_audit(values, modulus, pilot_period)
                    key = tuple(audit["delay_indices"])
                    unique.setdefault(key, audit)
                    if len(unique) >= pool_limit:
                        stopped_at_pool_limit = True
                    return
                residue = int(residues[position])
                for lift in range(modulus // pilot_period):
                    lift_nodes += 1
                    candidate = residue + pilot_period * lift
                    new_sums = [int((candidate + old) % modulus) for old in values]
                    new_sums.append(int((2 * candidate) % modulus))
                    if len(set(new_sums)) != len(new_sums) or any(value in sums for value in new_sums):
                        continue
                    extend(position + 1, [*values, candidate], sums.union(new_sums))

            extend(1, [0], {0})
            if stopped_at_pool_limit:
                break
        attempted[target_gap] = {"canonical_residue_sets": residue_sets, "lift_search_nodes": lift_nodes,
                                 "strict_sidon_candidates": len(unique),
                                 "search_complete": int(not stopped_at_pool_limit)}
        if unique:
            ranked = sorted(unique.values(), key=lambda row: (
                -int(row["fold_gap"]), -int(row["pair_gap"]),
                max(row["delay_indices"]), tuple(row["delay_indices"])))
            break
    if not ranked:
        raise RuntimeError("No strict Sidon set exists in the frozen Plan-039 search space.")
    return {
        "method": "rule_inference_residue_lift_enumeration",
        "uses_channel_or_performance_data": False,
        "modulus": modulus,
        "pilot_period": pilot_period,
        "branches": branches,
        "theoretical_fold_upper_bound": theoretical_fold_upper_bound,
        "attempted_fold_gaps": attempted,
        "canonical_feasible_count": len(ranked),
        "maximum_fold_gap": int(ranked[0]["fold_gap"]),
        "top8": ranked[:8],
    }


def _canonical(payload: dict[str, Any]) -> bytes:
    return (json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def load_manifest(path: str | Path, expected_sha256: str) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    source = Path(path)
    manifest = json.loads(source.read_text(encoding="utf-8"))
    embedded = str(manifest.pop("manifest_sha256", ""))
    actual = hashlib.sha256(_canonical(manifest)).hexdigest()
    manifest["manifest_sha256"] = embedded
    if actual != embedded or actual != str(expected_sha256).lower():
        raise ValueError(f"Plan-039 manifest hash mismatch: computed={actual}, embedded={embedded}, expected={expected_sha256}.")
    if manifest.get("schema") != "plan039-angular-full-coverage-manifest-v1" or manifest.get("status") != "FROZEN":
        raise ValueError("Plan-039 manifest is not a frozen wide-beam manifest.")
    arrays_path = (source.parent / str(manifest["arrays_file"])).resolve()
    if hashlib.sha256(arrays_path.read_bytes()).hexdigest() != manifest["arrays_sha256"]:
        raise ValueError("Plan-039 numeric artifact hash mismatch.")
    with np.load(arrays_path) as data:
        arrays = {name: data[name].copy() for name in data.files}
    if arrays["narrow_weights"].shape != (32, 8):
        raise ValueError("Plan-039 narrow codebook must have shape [32,8].")
    return manifest, arrays


def write_manifest(path: Path, payload: dict[str, Any]) -> str:
    body = dict(payload)
    body.pop("manifest_sha256", None)
    digest = hashlib.sha256(_canonical(body)).hexdigest()
    body["manifest_sha256"] = digest
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(_canonical(body))
    return digest


def frequency_covariance(
    method: str,
    frequencies_hz: Sequence[float],
    arrays: dict[str, np.ndarray],
    delay_indices: Sequence[int],
    alpha: Sequence[float],
) -> np.ndarray:
    method_u = str(method).upper()
    if method_u in {
        "PLAN039_TRANSPARENT_COMMON_REFERENCE_PDP",
        "PLAN039_PRG_COMMON_REFERENCE_PDP",
    }:
        return covariance_common_reference_pdp(
            frequencies_hz, arrays["path_delays_s"], arrays["reference_path_powers"], [0.0], np.ones(len(frequencies_hz))
        )
    artificial = np.asarray(delay_indices, dtype=float) / (576.0 * 30e3)
    if method_u == "PLAN039_COMMON_REFERENCE_PDP":
        return covariance_common_reference_pdp(frequencies_hz, arrays["path_delays_s"],
                                                arrays["reference_path_powers"], artificial, alpha)
    if method_u == "PLAN039_BEAM_SPECIFIC_PDP_INDEPENDENT":
        return covariance_beam_specific_independent(frequencies_hz, arrays["path_delays_s"],
                                                    arrays["beam_path_powers"], artificial, alpha)
    raise ValueError(f"Unsupported Plan-039 CDD covariance method {method!r}.")


def build_precoder(
    grid,
    scheme: str,
    weights: np.ndarray,
    sidon_delay_indices: Sequence[int] | None = None,
):
    from cdd_lls.phy.precoding import PrecoderResult

    scheme_u = str(scheme).upper()
    if scheme_u == "BEAM8_PRECODER_CYCLING":
        prg = np.arange(grid.n_sc, dtype=np.int64) // 72
        matrix = np.asarray(weights)[:, prg].T
        alpha = np.ones(grid.n_sc)
        delays: tuple[int, ...] = ()
        orthogonal = False
    else:
        if scheme_u == "BEAM8_B0_QC":
            delays = (0, 9, 18, 27, 36, 45, 54, 63)
        elif scheme_u == "BEAM8_SIDON_SELECTED":
            if sidon_delay_indices is None:
                raise ValueError("BEAM8_SIDON_SELECTED requires frozen Sidon delay indices.")
            delays = tuple(int(value) for value in sidon_delay_indices)
        else:
            raise ValueError(f"Unsupported Plan-039 transmission scheme {scheme!r}.")
        matrix, alpha, orthogonal = beam_domain_cdd_precoder(weights, delays, np.arange(grid.n_sc), 576)
    return PrecoderResult(C=matrix, label=scheme_u, metadata={
        "delay_grid_indices": list(delays), "alpha": alpha.tolist(), "orthogonal_branches": orthogonal,
        "raw_power_min": float(np.min(np.sum(np.abs(matrix) ** 2, axis=1))),
        "raw_power_max": float(np.max(np.sum(np.abs(matrix) ** 2, axis=1))),
    })
