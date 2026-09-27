"""Reusable CDL beam-design, PDP, covariance, and beam-domain CDD primitives.

The functions in this module implement the platform contract in Plan-037
section 13.  They intentionally do not depend on the legacy fixed DFT2x8
selection or on a link-level BLER loop.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np


ALGORITHM_VERSION = "plan037-section13-platform-v5"


def canonical_json(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json(value)).hexdigest()


def array_sha256(value: np.ndarray) -> str:
    array = np.ascontiguousarray(value)
    header = canonical_json({"dtype": array.dtype.str, "shape": list(array.shape)})
    return hashlib.sha256(header + array.view(np.uint8).tobytes()).hexdigest()


def regular_dft_codebook(
    n_vertical: int,
    n_horizontal: int,
    polarizations: int = 1,
    oversampling_vertical: int = 1,
    oversampling_horizontal: int = 1,
) -> tuple[np.ndarray, list[dict[str, float | int]]]:
    """Build the adaptive, dual-polarization-replicated regular DFT codebook."""
    nv, nh, npol = int(n_vertical), int(n_horizontal), int(polarizations)
    ov, oh = int(oversampling_vertical), int(oversampling_horizontal)
    if min(nv, nh, npol, ov, oh) <= 0:
        raise ValueError("DFT geometry and oversampling values must be positive.")
    beams: list[np.ndarray] = []
    metadata: list[dict[str, float | int]] = []
    for qv in range(ov * nv):
        for qh in range(oh * nh):
            spatial = np.asarray(
                [
                    np.exp(-1j * 2.0 * np.pi * (qv * row / (ov * nv) + qh * col / (oh * nh)))
                    for row in range(nv)
                    for col in range(nh)
                ],
                dtype=np.complex128,
            )
            vector = np.tile(spatial, npol)
            vector /= np.linalg.norm(vector)
            beams.append(vector)
            metadata.append(
                {
                    "beam_index": len(beams) - 1,
                    "q_vertical": qv,
                    "q_horizontal": qh,
                    "spatial_frequency_vertical_cycles": qv / (ov * nv),
                    "spatial_frequency_horizontal_cycles": qh / (oh * nh),
                }
            )
    return np.column_stack(beams), metadata


def select_beams(
    covariance: np.ndarray,
    codebook: np.ndarray,
    count: int,
    explicit_indices: Sequence[int] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    covariance = np.asarray(covariance, dtype=np.complex128)
    codebook = np.asarray(codebook, dtype=np.complex128)
    if covariance.shape != (codebook.shape[0], codebook.shape[0]):
        raise ValueError("Covariance and codebook dimensions do not match.")
    powers = np.einsum("tb,tu,ub->b", codebook.conj(), covariance, codebook, optimize=True).real
    if explicit_indices is None:
        order = np.lexsort((np.arange(powers.size), -powers))
        selected = order[: int(count)]
    else:
        selected = np.asarray(explicit_indices, dtype=np.int64)
    if selected.size != int(count) or len(set(selected.tolist())) != selected.size:
        raise ValueError("Selected beam indices must be unique and match num_branches.")
    if np.any(selected < 0) or np.any(selected >= codebook.shape[1]):
        raise ValueError("Selected beam index is outside the adaptive codebook.")
    return powers, selected


def direction_cosines(aod_rad: np.ndarray, zod_rad: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return horizontal and vertical direction cosines for a y-z panel."""
    aod = np.asarray(aod_rad, dtype=np.float64)
    zod = np.asarray(zod_rad, dtype=np.float64)
    return np.sin(zod) * np.sin(aod), np.cos(zod)


def txru_response_matrix(
    u_horizontal: np.ndarray,
    u_vertical: np.ndarray,
    vertical_aes: int,
    horizontal_aes: int,
    vertical_spacing_lambda: float,
    horizontal_spacing_lambda: float,
    mapping: np.ndarray,
    antenna_pattern: str = "omni",
) -> np.ndarray:
    """Return direct-AE TXRU response rows for a dual-polarized y-z panel."""
    uh = np.asarray(u_horizontal, dtype=np.float64).reshape(-1)
    uv = np.asarray(u_vertical, dtype=np.float64).reshape(-1)
    if uh.size != uv.size:
        raise ValueError("Horizontal and vertical coordinates must have equal size.")
    rows = np.arange(int(vertical_aes), dtype=np.float64)
    cols = np.arange(int(horizontal_aes), dtype=np.float64)
    position_v, position_h = np.meshgrid(rows, cols, indexing="ij")
    phase = 2.0 * np.pi * (
        uv[:, None] * float(vertical_spacing_lambda) * position_v.reshape(1, -1)
        + uh[:, None] * float(horizontal_spacing_lambda) * position_h.reshape(1, -1)
    )
    single = np.exp(1j * phase)
    pattern = str(antenna_pattern).lower()
    if pattern == "38.901":
        azimuth_deg = np.rad2deg(np.arcsin(np.clip(uh, -1.0, 1.0)))
        zenith_deg = np.rad2deg(np.arccos(np.clip(uv, -1.0, 1.0)))
        attenuation = np.minimum(
            12.0 * (azimuth_deg / 65.0) ** 2 + 12.0 * ((zenith_deg - 90.0) / 65.0) ** 2,
            30.0,
        )
        voltage_gain = np.sqrt(10.0 ** ((8.0 - attenuation) / 10.0))
        single *= voltage_gain[:, None]
    elif pattern != "omni":
        raise ValueError(f"Unsupported antenna pattern {antenna_pattern!r} for deterministic grid response.")
    ae = np.concatenate((single, single), axis=1)
    mapping = np.asarray(mapping, dtype=np.complex128)
    if ae.shape[1] != mapping.shape[0]:
        raise ValueError("AE response and AE-to-TXRU mapping dimensions do not match.")
    return ae @ mapping


def validate_regular_txru_mapping(
    mapping: np.ndarray,
    vertical_aes: int,
    horizontal_aes: int,
    vertical_txrus: int,
    horizontal_txrus: int,
    polarizations: int,
) -> None:
    """Reject mappings that are not polarization-isolated rectangular subarrays."""
    matrix = np.asarray(mapping, dtype=np.complex128)
    ma, na, mt, nt, pol = map(int, (vertical_aes, horizontal_aes, vertical_txrus, horizontal_txrus, polarizations))
    if ma % mt or na % nt or matrix.shape != (pol * ma * na, pol * mt * nt):
        raise ValueError("REGULAR_DFT_UNSUPPORTED_GEOMETRY: incompatible mapping shape or grid divisibility.")
    support = np.zeros(matrix.shape, dtype=bool)
    block_m, block_n = ma // mt, na // nt
    for polarization in range(pol):
        for row in range(mt):
            for col in range(nt):
                port = polarization * mt * nt + row * nt + col
                for ae_row in range(row * block_m, (row + 1) * block_m):
                    for ae_col in range(col * block_n, (col + 1) * block_n):
                        support[polarization * ma * na + ae_row * na + ae_col, port] = True
    if not np.array_equal(np.abs(matrix) > 1e-14, support):
        raise ValueError("REGULAR_DFT_UNSUPPORTED_GEOMETRY: mapping support is not a regular separable grid.")
    for port in range(matrix.shape[1]):
        values = matrix[support[:, port], port]
        if not np.allclose(values, values[0], atol=1e-12, rtol=0.0):
            raise ValueError("REGULAR_DFT_UNSUPPORTED_GEOMETRY: a TXRU subarray is not equal-phase.")


@dataclass(frozen=True)
class CoverageRectangle:
    horizontal_min: float
    horizontal_max: float
    vertical_min: float
    vertical_max: float
    covered_power: float
    total_power: float

    @property
    def area(self) -> float:
        return (self.horizontal_max - self.horizontal_min) * (self.vertical_max - self.vertical_min)

    @property
    def coverage(self) -> float:
        return self.covered_power / self.total_power


@dataclass(frozen=True)
class AngularRegion:
    """Closed AoD/ZoD region in degrees, ordered in vertical-major order."""

    aod_min_deg: float
    aod_max_deg: float
    zod_min_deg: float
    zod_max_deg: float

    @property
    def center_aod_deg(self) -> float:
        return 0.5 * (self.aod_min_deg + self.aod_max_deg)

    @property
    def center_zod_deg(self) -> float:
        return 0.5 * (self.zod_min_deg + self.zod_max_deg)


def partition_angular_region(
    aod_range_deg: Sequence[float],
    zod_range_deg: Sequence[float],
    vertical_beams: int,
    horizontal_beams: int,
) -> list[AngularRegion]:
    """Partition one angular rectangle into a vertical-major uniform grid."""
    aod = np.asarray(aod_range_deg, dtype=np.float64).reshape(-1)
    zod = np.asarray(zod_range_deg, dtype=np.float64).reshape(-1)
    kv, kh = int(vertical_beams), int(horizontal_beams)
    if aod.size != 2 or zod.size != 2 or kv <= 0 or kh <= 0:
        raise ValueError("Angular ranges must have two endpoints and beam-grid sizes must be positive.")
    if not np.all(np.isfinite(np.r_[aod, zod])) or aod[0] >= aod[1] or zod[0] >= zod[1]:
        raise ValueError("Angular ranges must be finite and strictly increasing.")
    aod_edges = np.linspace(aod[0], aod[1], kh + 1)
    zod_edges = np.linspace(zod[0], zod[1], kv + 1)
    return [
        AngularRegion(float(aod_edges[h]), float(aod_edges[h + 1]),
                      float(zod_edges[v]), float(zod_edges[v + 1]))
        for v in range(kv)
        for h in range(kh)
    ]


def dft_steering_beams(
    regions: Sequence[AngularRegion],
    response_at_centers: np.ndarray,
) -> np.ndarray:
    """Create unit-norm DFT/steering beams aimed at angular-region centers."""
    matrix = np.asarray(response_at_centers, dtype=np.complex128)
    if matrix.ndim != 2 or matrix.shape[0] != len(regions):
        raise ValueError("One TXRU response row is required for every angular region.")
    beams = matrix.conj().T
    norms = np.linalg.norm(beams, axis=0)
    if np.any(~np.isfinite(norms)) or np.any(norms <= 0.0):
        raise RuntimeError("DFT steering response contains a zero or non-finite vector.")
    return beams / norms[None, :]


def minimum_power_rectangle(
    u_horizontal: Sequence[float],
    u_vertical: Sequence[float],
    powers: Sequence[float],
    energy_coverage: float,
) -> CoverageRectangle:
    """Find the deterministic minimum-area axis-aligned weighted rectangle."""
    x = np.asarray(u_horizontal, dtype=np.float64).reshape(-1)
    y = np.asarray(u_vertical, dtype=np.float64).reshape(-1)
    p = np.asarray(powers, dtype=np.float64).reshape(-1)
    eta = float(energy_coverage)
    if x.size == 0 or x.size != y.size or x.size != p.size:
        raise ValueError("Coverage inputs must be non-empty arrays of equal length.")
    if not 0.0 < eta <= 1.0 or np.any(p < 0.0) or not np.all(np.isfinite(x + y + p)):
        raise ValueError("Coverage inputs or target are invalid.")
    total = float(np.sum(p))
    target = eta * total
    ux = np.unique(x)
    uy = np.unique(y)
    best: tuple[tuple[float, ...], CoverageRectangle] | None = None
    for left_index, left in enumerate(ux):
        for right in ux[left_index:]:
            horizontal = (x >= left) & (x <= right)
            if float(np.sum(p[horizontal])) + 1e-15 < target:
                continue
            vertical_mass = np.asarray([np.sum(p[horizontal & (y == value)]) for value in uy])
            start = 0
            running = 0.0
            for end, mass in enumerate(vertical_mass):
                running += float(mass)
                while start <= end and running - float(vertical_mass[start]) >= target - 1e-15:
                    running -= float(vertical_mass[start])
                    start += 1
                if running + 1e-15 < target:
                    continue
                bottom, top = float(uy[start]), float(uy[end])
                rectangle = CoverageRectangle(float(left), float(right), bottom, top, running, total)
                width_h = float(right - left)
                width_v = float(top - bottom)
                key = (rectangle.area, -running, width_h, width_v, float(left), bottom)
                if best is None or key < best[0]:
                    best = (key, rectangle)
    if best is None:
        raise RuntimeError("No power-coverage rectangle was found.")
    return best[1]


def split_rectangle(rectangle: CoverageRectangle, branches: int, split_axis: str, resolution_ratio: float) -> list[CoverageRectangle]:
    """Split a rectangle into equal-area, non-overlapping deterministic regions."""
    count = int(branches)
    if count <= 0:
        raise ValueError("Branch count must be positive.")
    axis = str(split_axis).lower()
    if axis not in {"auto", "horizontal", "vertical"}:
        raise ValueError("split_axis must be auto, horizontal, or vertical.")
    if axis == "horizontal":
        kv, kh = 1, count
    elif axis == "vertical":
        kv, kh = count, 1
    else:
        candidates = []
        width_h = max(rectangle.horizontal_max - rectangle.horizontal_min, np.finfo(float).eps)
        width_v = max(rectangle.vertical_max - rectangle.vertical_min, np.finfo(float).eps)
        for kv_candidate in range(1, count + 1):
            if count % kv_candidate:
                continue
            kh_candidate = count // kv_candidate
            cell_ratio = (width_h / kh_candidate) / (width_v / kv_candidate)
            score = abs(np.log(cell_ratio / max(float(resolution_ratio), np.finfo(float).eps)))
            candidates.append((score, -kh_candidate, kv_candidate, kh_candidate))
        _, _, kv, kh = min(candidates)
    h_edges = np.linspace(rectangle.horizontal_min, rectangle.horizontal_max, kh + 1)
    v_edges = np.linspace(rectangle.vertical_min, rectangle.vertical_max, kv + 1)
    return [
        CoverageRectangle(float(h_edges[h]), float(h_edges[h + 1]), float(v_edges[v]), float(v_edges[v + 1]), 0.0, rectangle.total_power)
        for v in range(kv)
        for h in range(kh)
    ]


def synthesize_region_beam(
    rectangle: CoverageRectangle,
    response: np.ndarray,
    grid_horizontal: np.ndarray,
    grid_vertical: np.ndarray,
    regularization: float,
    weight_constraint: str,
) -> np.ndarray:
    """Synthesize one deterministic unit-norm beam by weighted least squares."""
    matrix = np.asarray(response, dtype=np.complex128)
    gh = np.asarray(grid_horizontal).reshape(-1)
    gv = np.asarray(grid_vertical).reshape(-1)
    inside = (
        (gh >= rectangle.horizontal_min)
        & (gh <= rectangle.horizontal_max)
        & (gv >= rectangle.vertical_min)
        & (gv <= rectangle.vertical_max)
    )
    if not np.any(inside):
        center_h = 0.5 * (rectangle.horizontal_min + rectangle.horizontal_max)
        center_v = 0.5 * (rectangle.vertical_min + rectangle.vertical_max)
        distance = (gh - center_h) ** 2 + (gv - center_v) ** 2
        inside[int(np.argmin(distance))] = True
    target = inside.astype(np.complex128)
    sample_weight = np.where(inside, 4.0, 1.0)
    normal = matrix.conj().T @ (sample_weight[:, None] * matrix)
    rhs = matrix.conj().T @ (sample_weight * target)
    normal += float(regularization) * np.eye(matrix.shape[1])
    weight = np.linalg.solve(normal, rhs)
    constraint = str(weight_constraint).lower()
    if constraint == "constant_modulus":
        weight = np.exp(1j * np.angle(weight))
    elif constraint != "unit_norm":
        raise ValueError("weight_constraint must be unit_norm or constant_modulus.")
    norm = float(np.linalg.norm(weight))
    if not np.isfinite(norm) or norm <= 0.0:
        raise RuntimeError("Wide-beam synthesis produced a zero or non-finite vector.")
    return weight / norm


def beam_path_statistics(weights: np.ndarray, path_covariances: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return per-beam path powers and per-path beam joint covariance."""
    weights = np.asarray(weights, dtype=np.complex128)
    covariances = np.asarray(path_covariances, dtype=np.complex128)
    joint = np.einsum("tm,ltu,un->lmn", weights.conj(), covariances, weights, optimize=True)
    joint = 0.5 * (joint + joint.conj().transpose(0, 2, 1))
    powers = np.maximum(np.real(np.diagonal(joint, axis1=1, axis2=2)).T, 0.0)
    return powers, joint


def maximum_normalized_pattern_capture(
    weight: np.ndarray,
    ray_response: np.ndarray,
    ray_powers: Sequence[float],
    visible_grid_response: np.ndarray,
) -> tuple[float, float]:
    """Return raw-ray-power-weighted capture after peak-normalizing a beam pattern.

    The peak is evaluated on the supplied visible-direction grid.  The returned
    ratio is therefore insensitive to the absolute element/array gain scale,
    while retaining the angular shape of the configured antenna pattern.
    """
    vector = np.asarray(weight, dtype=np.complex128).reshape(-1)
    rays = np.asarray(ray_response, dtype=np.complex128)
    grid = np.asarray(visible_grid_response, dtype=np.complex128)
    powers = np.asarray(ray_powers, dtype=np.float64).reshape(-1)
    if rays.ndim != 2 or grid.ndim != 2 or rays.shape[1] != vector.size or grid.shape[1] != vector.size:
        raise ValueError("Beam weight and direction-response dimensions do not match.")
    if rays.shape[0] != powers.size or np.any(powers < 0.0) or not np.all(np.isfinite(powers)):
        raise ValueError("Ray powers must be finite, nonnegative, and aligned with ray responses.")
    total_power = float(np.sum(powers))
    if total_power <= 0.0:
        raise ValueError("Ray powers must have positive total power.")
    ray_gain = np.abs(rays @ vector) ** 2
    peak_gain = float(max(np.max(np.abs(grid @ vector) ** 2), np.max(ray_gain)))
    if not np.isfinite(peak_gain) or peak_gain <= 0.0:
        raise RuntimeError("Beam pattern has zero or non-finite peak gain.")
    normalized_ray_gain = ray_gain / peak_gain
    return float(np.sum(powers * normalized_ray_gain) / total_power), peak_gain


def beam_domain_cdd_precoder(
    weights: np.ndarray,
    delay_indices: Sequence[float],
    subcarrier_indices: Sequence[float],
    active_subcarrier_count: int | None = None,
) -> tuple[np.ndarray, np.ndarray, bool]:
    """Build unit-power CDD weights for orthogonal or non-orthogonal branches."""
    branches = np.asarray(weights, dtype=np.complex128)
    delays = np.asarray(delay_indices, dtype=np.float64).reshape(-1)
    subcarriers = np.asarray(subcarrier_indices, dtype=np.float64).reshape(-1)
    if branches.ndim != 2 or branches.shape[1] != delays.size:
        raise ValueError("CDD branch weights and delay indices do not match.")
    denominator = float(subcarriers.size if active_subcarrier_count is None else active_subcarrier_count)
    phase = np.exp(-1j * 2.0 * np.pi * subcarriers[:, None] * delays[None, :] / denominator)
    unnormalized = phase @ branches.T
    gram = branches.conj().T @ branches
    orthogonal = bool(np.allclose(gram, np.eye(delays.size), atol=1e-12, rtol=0.0))
    if orthogonal:
        alpha = np.full(subcarriers.size, 1.0 / np.sqrt(delays.size))
    else:
        raw_norm = np.linalg.norm(unnormalized, axis=1)
        if np.any(raw_norm <= 0.0):
            raise RuntimeError("Non-orthogonal CDD synthesis has a zero-power subcarrier.")
        alpha = 1.0 / raw_norm
    precoder = alpha[:, None] * unnormalized
    if not np.allclose(np.sum(np.abs(precoder) ** 2, axis=1), 1.0, atol=1e-11, rtol=0.0):
        raise RuntimeError("CDD precoder violates the per-subcarrier unit-power constraint.")
    return precoder, alpha, orthogonal


def _frequency_phase(frequencies_hz: np.ndarray, delays_s: np.ndarray) -> np.ndarray:
    return np.exp(-1j * 2.0 * np.pi * frequencies_hz[:, None] * delays_s[None, :])


def covariance_common_reference_pdp(
    frequencies_hz: Sequence[float], delays_s: Sequence[float], powers: Sequence[float], artificial_delays_s: Sequence[float],
    alpha: Sequence[float] | None = None,
) -> np.ndarray:
    f = np.asarray(frequencies_hz, dtype=np.float64)
    tau = np.asarray(delays_s, dtype=np.float64)
    p = np.asarray(powers, dtype=np.float64)
    artificial = np.asarray(artificial_delays_s, dtype=np.float64)
    joint = np.asarray([value * np.eye(artificial.size) for value in p], dtype=np.complex128)
    return covariance_beam_joint(f, tau, joint, artificial, alpha)


def covariance_beam_specific_independent(
    frequencies_hz: Sequence[float], delays_s: Sequence[float], beam_path_powers: np.ndarray, artificial_delays_s: Sequence[float],
    alpha: Sequence[float] | None = None,
) -> np.ndarray:
    f = np.asarray(frequencies_hz, dtype=np.float64)
    tau = np.asarray(delays_s, dtype=np.float64)
    power = np.asarray(beam_path_powers, dtype=np.float64)
    artificial = np.asarray(artificial_delays_s, dtype=np.float64)
    if power.shape != (artificial.size, tau.size):
        raise ValueError("Beam-specific powers must have shape [branch,path].")
    joint = np.asarray([np.diag(power[:, path]) for path in range(tau.size)], dtype=np.complex128)
    return covariance_beam_joint(f, tau, joint, artificial, alpha)


def covariance_beam_joint(
    frequencies_hz: Sequence[float], delays_s: Sequence[float], joint_path_covariances: np.ndarray,
    artificial_delays_s: Sequence[float], alpha: Sequence[float] | None = None,
) -> np.ndarray:
    f = np.asarray(frequencies_hz, dtype=np.float64)
    tau = np.asarray(delays_s, dtype=np.float64)
    joint = np.asarray(joint_path_covariances, dtype=np.complex128)
    artificial = np.asarray(artificial_delays_s, dtype=np.float64)
    if joint.shape != (tau.size, artificial.size, artificial.size):
        raise ValueError("Joint path covariance must have shape [path,branch,branch].")
    scale = np.full(f.size, 1.0 / np.sqrt(artificial.size)) if alpha is None else np.asarray(alpha, dtype=np.float64)
    coefficients = scale[:, None] * np.exp(-1j * 2.0 * np.pi * f[:, None] * artificial[None, :])
    result = np.zeros((f.size, f.size), dtype=np.complex128)
    for path, delay in enumerate(tau):
        spatial = coefficients @ joint[path] @ coefficients.conj().T
        physical = np.exp(-1j * 2.0 * np.pi * (f[:, None] - f[None, :]) * delay)
        result += spatial * physical
    return 0.5 * (result + result.conj().T)


def covariance_audit(covariance: np.ndarray) -> dict[str, float]:
    matrix = np.asarray(covariance, dtype=np.complex128)
    eig = np.linalg.eigvalsh(0.5 * (matrix + matrix.conj().T))
    return {
        "hermitian_max_error": float(np.max(np.abs(matrix - matrix.conj().T))),
        "minimum_eigenvalue": float(np.min(eig)),
        "diagonal_min": float(np.min(np.real(np.diag(matrix)))),
        "diagonal_max": float(np.max(np.real(np.diag(matrix)))),
    }


def validate_cache(cache_dir: Path, expected_key: str) -> tuple[dict[str, Any], dict[str, np.ndarray]] | None:
    manifest_path = cache_dir / "manifest.json"
    weights_path = cache_dir / "weights.npz"
    if not manifest_path.exists() or not weights_path.exists():
        return None
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    embedded = manifest.get("manifest_sha256", "")
    payload = {key: value for key, value in manifest.items() if key != "manifest_sha256"}
    if embedded != sha256_json(payload) or manifest.get("cache_key") != expected_key:
        return None
    with np.load(weights_path) as loaded:
        arrays = {name: np.asarray(loaded[name]) for name in loaded.files}
    for name in ("regular_weights", "selected_regular_weights", "ssb_weights", "wide_weight", "split_weights"):
        if name not in arrays or not np.all(np.isfinite(arrays[name])):
            return None
        norms = np.linalg.norm(arrays[name], axis=0) if arrays[name].ndim == 2 else np.asarray([np.linalg.norm(arrays[name])])
        if not np.allclose(norms, 1.0, atol=1e-10, rtol=0.0):
            return None
        if manifest.get("array_sha256", {}).get(name) != array_sha256(arrays[name]):
            return None
    for name, digest in manifest.get("artifact_sha256", {}).items():
        path = cache_dir / name
        if not path.exists() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            return None
    return manifest, arrays


def write_cache(cache_dir: Path, manifest: dict[str, Any], arrays: dict[str, np.ndarray]) -> str:
    cache_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache_dir / "weights.npz", **arrays)
    payload = dict(manifest)
    payload["array_sha256"] = {name: array_sha256(value) for name, value in arrays.items()}
    digest = sha256_json(payload)
    payload["manifest_sha256"] = digest
    (cache_dir / "manifest.json").write_bytes(canonical_json(payload))
    return digest


def rows_to_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    import csv

    materialized = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not materialized:
        raise ValueError("CSV output requires at least one row.")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(materialized[0]))
        writer.writeheader()
        writer.writerows(materialized)
