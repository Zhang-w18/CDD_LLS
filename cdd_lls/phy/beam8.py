"""Plan-037 32-TXRU DFT2x8 beam reduction primitives.

This module is deliberately independent of the legacy parent/secondary
codebooks used by :mod:`channel_cdl_fixed`.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from cdd_lls.phy.cdl_beam_platform import beam_domain_cdd_precoder
from cdd_lls.phy.precoding import PrecoderResult
from cdd_lls.phy.resource_grid import ResourceGrid


CODEBOOK_TYPE = "dft_2x8_same_pol"
SCHEMES = ("BEAM8_B0_QC", "BEAM8_S0_SIDON", "BEAM8_PRECODER_CYCLING")
DELAYS = {
    "BEAM8_B0_QC": (0, 9, 18, 27, 36, 45, 54, 63),
    "BEAM8_S0_SIDON": (0, 1, 3, 7, 12, 20, 30, 65),
}


def merge_discrete_pdp(delays_s: Sequence[float], powers: Sequence[float], atol_s: float = 1e-15) -> tuple[np.ndarray, np.ndarray]:
    """Merge coincident delays, reject invalid mass, and return a normalized PDP."""
    delays = np.asarray(delays_s, dtype=np.float64).reshape(-1)
    mass = np.asarray(powers, dtype=np.float64).reshape(-1)
    if delays.size != mass.size or delays.size == 0:
        raise ValueError("PDP delays and powers must be non-empty arrays of equal length.")
    if not np.all(np.isfinite(delays)) or not np.all(np.isfinite(mass)) or np.any(mass < 0.0):
        raise ValueError("PDP delays/powers must be finite and powers non-negative.")
    order = np.argsort(delays, kind="stable")
    merged_d: list[float] = []
    merged_p: list[float] = []
    for delay, power in zip(delays[order], mass[order]):
        if merged_d and abs(float(delay) - merged_d[-1]) <= float(atol_s):
            merged_p[-1] += float(power)
        else:
            merged_d.append(float(delay)); merged_p.append(float(power))
    total = float(np.sum(merged_p))
    if not np.isfinite(total) or total <= 0.0:
        raise ValueError("PDP total power must be finite and positive.")
    return np.asarray(merged_d), np.asarray(merged_p) / total


def pdp_moments(delays_s: Sequence[float], powers: Sequence[float]) -> tuple[float, float]:
    delays, power = merge_discrete_pdp(delays_s, powers)
    mean = float(np.sum(delays * power))
    rms = float(np.sqrt(np.sum(power * (delays - mean) ** 2)))
    return mean, rms


def circular_pdp_moments(delays_s: Sequence[float], powers: Sequence[float], period_s: float) -> tuple[float, float]:
    """Return the minimum squared circular-distance mean and RMS."""
    delays, power = merge_discrete_pdp(delays_s, powers)
    period = float(period_s)
    if not np.isfinite(period) or period <= 0.0:
        raise ValueError("Circular PDP period must be finite and positive.")
    wrapped = np.mod(delays, period)
    candidates = []
    for anchor in wrapped:
        unwrapped = anchor + ((wrapped - anchor + period / 2.0) % period - period / 2.0)
        mean = float(np.sum(power * unwrapped))
        residual = (wrapped - mean + period / 2.0) % period - period / 2.0
        wrapped_mean = mean % period
        if np.isclose(wrapped_mean, period, atol=1e-12 * period, rtol=0.0):
            wrapped_mean = 0.0
        candidates.append((float(np.sum(power * residual**2)), wrapped_mean))
    objective, mean = min(candidates, key=lambda item: item[0])
    return float(mean), float(np.sqrt(max(objective, 0.0)))


def shifted_reference_pdp(delays_s: Sequence[float], powers: Sequence[float], delay_indices: Sequence[int], n_sc: int, scs_hz: float) -> tuple[np.ndarray, np.ndarray]:
    period = 1.0 / float(scs_hz)
    artificial = np.asarray(delay_indices, dtype=np.float64) / (float(n_sc) * float(scs_hz))
    delays = (np.asarray(delays_s, dtype=np.float64)[:, None] + artificial[None, :]).reshape(-1)
    mass = np.repeat(np.asarray(powers, dtype=np.float64) / len(artificial), len(artificial))
    return merge_discrete_pdp(np.mod(delays, period), mass)


def frequency_covariance(delays_s: Sequence[float], powers: Sequence[float], n_sc: int, scs_hz: float) -> np.ndarray:
    delays, power = merge_discrete_pdp(delays_s, powers)
    delta = np.arange(n_sc)[:, None] - np.arange(n_sc)[None, :]
    return np.sum(power[None, None, :] * np.exp(-1j * 2.0 * np.pi * delta[:, :, None] * float(scs_hz) * delays[None, None, :]), axis=2)


def canonical_json_bytes(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def manifest_payload(manifest: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in manifest.items() if key != "manifest_sha256"}


def manifest_sha256(manifest: dict[str, Any]) -> str:
    return sha256_bytes(canonical_json_bytes(manifest_payload(manifest)))


def write_manifest(path: Path, manifest: dict[str, Any]) -> str:
    payload = manifest_payload(manifest)
    digest = manifest_sha256(payload)
    payload["manifest_sha256"] = digest
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(canonical_json_bytes(payload))
    return digest


def load_frozen_manifest(path: str | Path, expected_sha256: str) -> dict[str, Any]:
    source = Path(path)
    manifest = json.loads(source.read_text(encoding="utf-8"))
    actual = manifest_sha256(manifest)
    embedded = str(manifest.get("manifest_sha256", ""))
    if actual != embedded or actual != str(expected_sha256).lower():
        raise ValueError(
            f"Frozen beam manifest hash mismatch: computed={actual}, embedded={embedded}, expected={expected_sha256}."
        )
    if manifest.get("status") != "FROZEN" or manifest.get("codebook_type") != CODEBOOK_TYPE:
        raise ValueError("Beam manifest is not a frozen dft_2x8_same_pol selection.")
    indices = [int(value) for value in manifest.get("selected_beam_indices", [])]
    if len(indices) != 8 or indices != sorted(set(indices)) or not all(0 <= value < 16 for value in indices):
        raise ValueError("Frozen beam manifest must contain eight unique ascending indices in [0,15].")
    powers = manifest.get("mean_rsrp", [])
    if len(powers) != 16 or not all(np.isfinite(float(value)) and float(value) > 0.0 for value in powers):
        raise ValueError("Frozen beam manifest must contain sixteen finite positive mean_rsrp values.")
    return manifest


def build_dft_2x8_same_pol_codebook() -> np.ndarray:
    """Return B[:, q_v*8+q_h] using vertical-major/horizontal-minor TXRU order."""
    rows = []
    for qv in range(2):
        for qh in range(8):
            spatial = np.asarray(
                [np.exp(-1j * 2.0 * np.pi * (qv * m / 2.0 + qh * n / 8.0)) / 4.0
                 for m in range(2) for n in range(8)],
                dtype=np.complex128,
            )
            rows.append(np.concatenate((spatial, spatial)) / np.sqrt(2.0))
    codebook = np.column_stack(rows)
    if not np.allclose(codebook.conj().T @ codebook, np.eye(16), atol=1e-13, rtol=0.0):
        raise RuntimeError("DFT2x8 codebook is not orthonormal.")
    return codebook


def build_beam8_precoder(
    grid: ResourceGrid,
    scheme: str,
    selected_beam_indices: Sequence[int],
    codebook: np.ndarray | None = None,
    prg_size_rb: int = 6,
) -> PrecoderResult:
    scheme = str(scheme).upper()
    if scheme not in SCHEMES:
        raise ValueError(f"Unsupported Beam8 scheme {scheme!r}.")
    indices = np.asarray(selected_beam_indices, dtype=np.int64)
    if indices.tolist() != sorted(set(indices.tolist())) or indices.size != 8:
        raise ValueError("Beam8 indices must be eight unique indices in ascending order.")
    full = build_dft_2x8_same_pol_codebook() if codebook is None else np.asarray(codebook)
    branches = full[:, indices]
    n_sc = int(grid.n_sc)
    if scheme == "BEAM8_PRECODER_CYCLING":
        if int(prg_size_rb) != 6 or n_sc != 576:
            raise ValueError("Plan-037 Beam8 cycling requires 576 active SC and 6-RB PRGs.")
        prg = np.arange(n_sc, dtype=np.int64) // 72
        if int(prg[-1]) != 7:
            raise RuntimeError("Beam8 cycling did not produce exactly eight PRGs.")
        weights = branches[:, prg].T
        delays = ()
    else:
        delays = DELAYS[scheme]
        weights, alpha, orthogonal = beam_domain_cdd_precoder(
            branches, delays, np.arange(n_sc, dtype=np.float64), active_subcarrier_count=576
        )
    power = np.sum(np.abs(weights) ** 2, axis=1).real
    if not np.allclose(power, 1.0, atol=1e-12, rtol=0.0):
        raise RuntimeError(f"Beam8 precoder violates unit power: max error={np.max(np.abs(power - 1.0))}.")
    return PrecoderResult(
        C=weights,
        label=scheme,
        metadata={
            "codebook_type": CODEBOOK_TYPE,
            "selected_beam_indices": indices.tolist(),
            "delay_grid_indices": list(delays),
            "phase_denominator": 576,
            "prg_size_rb": int(prg_size_rb),
            "raw_power_min": float(np.min(power)),
            "raw_power_max": float(np.max(power)),
            "normalization": "per_subcarrier_exact_unit_norm",
            "orthogonal_branches": bool(scheme == "BEAM8_PRECODER_CYCLING" or orthogonal),
            "alpha": [] if scheme == "BEAM8_PRECODER_CYCLING" else alpha.tolist(),
        },
    )


def equivalent_channel_two_paths(h: np.ndarray, precoder: np.ndarray, branches: np.ndarray, scheme: str) -> tuple[np.ndarray, np.ndarray]:
    """Return direct H@w and branch-first construction for a test/audit tensor [rx,sc,tx]."""
    direct = np.einsum("rkt,kt->rk", h, precoder, optimize=True)
    branch_h = np.einsum("rkt,tm->rkm", h, branches, optimize=True)
    if scheme == "BEAM8_PRECODER_CYCLING":
        prg = np.arange(h.shape[1], dtype=np.int64) // 72
        reconstructed = branch_h[np.arange(h.shape[0])[:, None], np.arange(h.shape[1])[None, :], prg[None, :]]
    else:
        delays = np.asarray(DELAYS[scheme], dtype=np.float64)
        phase = np.exp(-1j * 2.0 * np.pi * np.arange(h.shape[1])[:, None] * delays[None, :] / 576.0)
        reconstructed = np.einsum("rkm,km->rk", branch_h, phase, optimize=True) / np.sqrt(8.0)
    return direct, reconstructed
