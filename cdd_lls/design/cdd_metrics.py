from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np


def delay_indices_to_ns(
    delay_indices: Sequence[int],
    n_active_subcarriers: int = 576,
    subcarrier_spacing_hz: float = 30e3,
) -> np.ndarray:
    scale = 1e9 / (float(n_active_subcarriers) * float(subcarrier_spacing_hz))
    return np.asarray(delay_indices, dtype=np.float64) * scale


def array_factor(delay_indices: Sequence[int], period: int = 576) -> np.ndarray:
    delays = np.asarray(delay_indices, dtype=np.float64)
    lags = np.arange(int(period), dtype=np.float64)
    return np.mean(
        np.exp(-2j * np.pi * lags[:, None] * delays[None, :] / float(period)),
        axis=1,
    )


def phase_matrix_from_grid_coordinates(
    delay_grid_coordinates: Sequence[float],
    period: int = 576,
) -> np.ndarray:
    """Return V[k,n] for integer or fractional active-band DFT coordinates."""
    coordinates = np.asarray(delay_grid_coordinates, dtype=np.float64).reshape(-1)
    subcarriers = np.arange(int(period), dtype=np.float64)
    return np.exp(
        -2j
        * np.pi
        * subcarriers[:, None]
        * coordinates[None, :]
        / float(period)
    )


def delay_ns_to_grid_coordinates(
    delay_ns: Sequence[float],
    n_active_subcarriers: int = 576,
    subcarrier_spacing_hz: float = 30e3,
) -> np.ndarray:
    return (
        np.asarray(delay_ns, dtype=np.float64)
        * 1e-9
        * float(n_active_subcarriers)
        * float(subcarrier_spacing_hz)
    )


def effective_moments_direct(
    delay_grid_coordinates: Sequence[float],
    physical_frequency_covariance: np.ndarray,
) -> dict[str, float]:
    """Compute full-matrix off-diagonal M2/M4 effective correlation moments."""
    covariance = np.asarray(physical_frequency_covariance, dtype=np.complex128)
    if covariance.ndim != 2 or covariance.shape[0] != covariance.shape[1]:
        raise ValueError("Physical frequency covariance must be square.")
    period = int(covariance.shape[0])
    diagonal = np.real(np.diag(covariance))
    if np.any(diagonal <= 0.0):
        raise ValueError("Physical frequency covariance must have positive diagonal.")
    normalized = covariance / np.sqrt(diagonal[:, None] * diagonal[None, :])
    phase = phase_matrix_from_grid_coordinates(delay_grid_coordinates, period)
    artificial = phase @ phase.conj().T / float(phase.shape[1])
    magnitude = np.abs(normalized * artificial)
    np.fill_diagonal(magnitude, 0.0)
    return {
        "m2_eff": float(np.sum(magnitude**2)),
        "m4_eff": float(np.sum(magnitude**4)),
        "max_effective_correlation": float(np.max(magnitude)),
        "mean_effective_correlation_offdiag": float(
            np.sum(magnitude) / float(period * (period - 1))
        ),
    }


def effective_moments_lag(
    delay_grid_coordinates: Sequence[float],
    physical_frequency_covariance: np.ndarray,
    stationarity_tolerance: float = 1e-10,
) -> dict[str, float]:
    """Compute M2/M4 by frequency lag for a stationary Toeplitz covariance."""
    covariance = np.asarray(physical_frequency_covariance, dtype=np.complex128)
    if covariance.ndim != 2 or covariance.shape[0] != covariance.shape[1]:
        raise ValueError("Physical frequency covariance must be square.")
    period = int(covariance.shape[0])
    diagonal = np.real(np.diag(covariance))
    if np.any(diagonal <= 0.0):
        raise ValueError("Physical frequency covariance must have positive diagonal.")
    diagonal_mean = float(np.mean(diagonal))
    if not np.allclose(diagonal, diagonal_mean, rtol=stationarity_tolerance, atol=stationarity_tolerance):
        raise ValueError("Lag implementation requires a constant covariance diagonal.")
    lags = np.arange(1, period, dtype=np.float64)
    reference = covariance[0, 1:] / diagonal_mean
    for offset in range(1, period):
        values = np.diag(covariance, k=offset) / diagonal_mean
        if not np.allclose(
            values,
            reference[offset - 1],
            rtol=stationarity_tolerance,
            atol=stationarity_tolerance,
        ):
            raise ValueError("Lag implementation requires a stationary Toeplitz covariance.")
    coordinates = np.asarray(delay_grid_coordinates, dtype=np.float64).reshape(-1)
    factor = np.mean(
        np.exp(
            -2j
            * np.pi
            * lags[:, None]
            * coordinates[None, :]
            / float(period)
        ),
        axis=1,
    )
    rho = np.abs(reference * factor)
    weights = 2.0 * (period - lags)
    return {
        "m2_eff": float(np.sum(weights * rho**2)),
        "m4_eff": float(np.sum(weights * rho**4)),
        "max_effective_correlation": float(np.max(rho)),
        "mean_effective_correlation_offdiag": float(
            np.sum(weights * rho) / float(period * (period - 1))
        ),
    }


def _circular_distance(a: int, b: int, period: int) -> int:
    delta = abs(int(a) - int(b)) % int(period)
    return int(min(delta, int(period) - delta))


def pair_sum_min_gap_indices(delay_indices: Sequence[int], period: int = 576) -> int:
    values = [int(value) % int(period) for value in delay_indices]
    pair_sums = np.sort(np.asarray([
        (values[first] + values[second]) % int(period)
        for first in range(len(values))
        for second in range(first, len(values))
    ], dtype=np.int64))
    adjacent = np.diff(pair_sums)
    wrap = int(pair_sums[0]) + int(period) - int(pair_sums[-1])
    return int(min(int(np.min(adjacent)), wrap))


def fold_min_gap_indices(delay_indices: Sequence[int], pilot_period: int = 24) -> int:
    residues = [int(value) % int(pilot_period) for value in delay_indices]
    if len(set(residues)) != len(residues):
        return 0
    minimum = int(pilot_period)
    for first in range(len(residues)):
        for second in range(first + 1, len(residues)):
            minimum = min(
                minimum,
                _circular_distance(residues[first], residues[second], int(pilot_period)),
            )
    return int(minimum)


def tdl_effective_support(
    delays_s: Sequence[float],
    powers: Sequence[float],
    epsilon: float = 0.01,
) -> dict[str, float | int]:
    """Shortest contiguous delay interval containing at least 1-epsilon power."""
    delays = np.asarray(delays_s, dtype=np.float64).reshape(-1)
    weights = np.asarray(powers, dtype=np.float64).reshape(-1)
    if delays.size != weights.size or delays.size == 0:
        raise ValueError("TDL delays and powers must have the same non-zero length.")
    if not 0.0 <= float(epsilon) < 1.0:
        raise ValueError("epsilon must lie in [0,1).")
    order = np.argsort(delays)
    delays = delays[order]
    weights = np.maximum(weights[order], 0.0)
    weights /= np.sum(weights)
    target = 1.0 - float(epsilon)
    best: tuple[float, float, float, float, int, int] | None = None
    for start in range(len(delays)):
        mass = 0.0
        for stop in range(start, len(delays)):
            mass += float(weights[stop])
            if mass >= target:
                candidate = (
                    float(delays[stop] - delays[start]),
                    float(delays[start]),
                    float(delays[stop]),
                    mass,
                    start,
                    stop,
                )
                if best is None or candidate < best:
                    best = candidate
                break
    if best is None:
        raise RuntimeError("Could not construct an effective TDL support interval.")
    return {
        "epsilon": float(epsilon),
        "target_power": target,
        "support_width_s": best[0],
        "support_start_s": best[1],
        "support_stop_s": best[2],
        "contained_power": best[3],
        "start_index": best[4],
        "stop_index": best[5],
    }


def equivalent_cdd_delay_spread(
    physical_delays_s: Sequence[float],
    physical_powers: Sequence[float],
    delay_grid_coordinates: Sequence[float],
    period: int = 576,
    subcarrier_spacing_hz: float = 30e3,
    contained_probability: float = 0.99,
) -> dict[str, float]:
    """Circular RMS and shortest-energy support of the CDD-composite PDP."""
    delays = np.asarray(physical_delays_s, dtype=np.float64).reshape(-1)
    powers = np.asarray(physical_powers, dtype=np.float64).reshape(-1)
    coordinates = np.asarray(delay_grid_coordinates, dtype=np.float64).reshape(-1)
    if delays.size == 0 or delays.size != powers.size:
        raise ValueError("Physical delays and powers must be non-empty and aligned.")
    if coordinates.size == 0:
        raise ValueError("At least one CDD delay coordinate is required.")
    if np.any(powers < 0.0) or not np.all(np.isfinite(powers)):
        raise ValueError("Physical powers must be finite and non-negative.")
    if not 0.0 < float(contained_probability) <= 1.0:
        raise ValueError("contained_probability must lie in (0, 1].")
    total_power = float(np.sum(powers))
    if total_power <= 0.0:
        raise ValueError("Physical powers must have positive total power.")

    symbol_period_s = 1.0 / float(subcarrier_spacing_hz)
    artificial_s = coordinates / (
        float(period) * float(subcarrier_spacing_hz)
    )
    composite = (
        delays[:, None] + artificial_s[None, :]
    ).reshape(-1) % symbol_period_s
    weights = np.repeat(
        powers / total_power / float(coordinates.size),
        coordinates.size,
    )

    # The squared geodesic objective is quadratic between antipodal
    # breakpoints. Enumerating those intervals gives the one-dimensional
    # Fréchet minimum up to floating-point arithmetic.
    breakpoints = np.unique(
        (composite + 0.5 * symbol_period_s) % symbol_period_s
    )
    boundaries = np.r_[breakpoints, breakpoints[0] + symbol_period_s]
    best_variance = math.inf
    best_center = 0.0
    for lower, upper in zip(boundaries[:-1], boundaries[1:]):
        midpoint = 0.5 * (lower + upper)
        unwrapped = midpoint + (
            (composite - midpoint + 0.5 * symbol_period_s) % symbol_period_s
            - 0.5 * symbol_period_s
        )
        center = float(np.sum(weights * unwrapped))
        center = min(max(center, float(lower)), float(upper))
        residual = (
            composite - center + 0.5 * symbol_period_s
        ) % symbol_period_s - 0.5 * symbol_period_s
        variance = float(np.sum(weights * residual**2))
        if variance < best_variance:
            best_variance = variance
            best_center = center % symbol_period_s

    order = np.argsort(composite, kind="mergesort")
    ordered = composite[order]
    ordered_weights = weights[order]
    doubled_delays = np.r_[ordered, ordered + symbol_period_s]
    doubled_weights = np.r_[ordered_weights, ordered_weights]
    target = float(contained_probability)
    best_width = symbol_period_s
    stop = 0
    running = 0.0
    for start in range(ordered.size):
        stop = max(stop, start)
        while stop < start + ordered.size and running + 1e-15 < target:
            running += float(doubled_weights[stop])
            stop += 1
        if running + 1e-15 >= target:
            best_width = min(
                best_width,
                float(doubled_delays[stop - 1] - doubled_delays[start]),
            )
        running -= float(doubled_weights[start])

    return {
        "equivalent_circular_rms_delay_spread_s": math.sqrt(
            max(best_variance, 0.0)
        ),
        "equivalent_circular_mean_delay_s": float(best_center),
        "equivalent_circular_support_width_s": float(best_width),
        "equivalent_support_contained_probability": target,
        "equivalent_delay_period_s": symbol_period_s,
    }


@dataclass(frozen=True)
class CDDKernel:
    snr_db: float
    rho_grid: np.ndarray
    covariance_values: np.ndarray

    def evaluate(self, rho: np.ndarray) -> np.ndarray:
        values = np.clip(np.asarray(rho, dtype=np.float64), 0.0, 1.0)
        return np.interp(values, self.rho_grid, self.covariance_values)


def _mi_from_power(
    power: np.ndarray,
    snr_linear: float,
    iqam_db_grid: np.ndarray,
    iqam_values: np.ndarray,
) -> np.ndarray:
    rho = np.maximum(float(snr_linear) * np.asarray(power, dtype=np.float64), 1e-300)
    rho_db = 10.0 * np.log10(rho)
    return np.interp(
        rho_db,
        np.asarray(iqam_db_grid, dtype=np.float64),
        np.asarray(iqam_values, dtype=np.float64),
        left=0.0,
        right=4.0,
    )


def build_cdd_kernel(
    snr_db: float,
    iqam_db_grid: np.ndarray,
    iqam_values: np.ndarray,
    rho_step: float = 0.005,
    quadrature_order: int = 10,
) -> CDDKernel:
    """Build the mutual-information covariance kernel for correlated CN(0,1) pairs."""
    nodes, weights = np.polynomial.hermite.hermgauss(int(quadrature_order))
    normal = np.sqrt(2.0) * nodes
    normal_weights = weights / math.sqrt(math.pi)
    x1, x2, x3, x4 = np.meshgrid(normal, normal, normal, normal, indexing="ij")
    w1, w2, w3, w4 = np.meshgrid(
        normal_weights,
        normal_weights,
        normal_weights,
        normal_weights,
        indexing="ij",
    )
    z1 = (x1.ravel() + 1j * x2.ravel()) / math.sqrt(2.0)
    independent = (x3.ravel() + 1j * x4.ravel()) / math.sqrt(2.0)
    product_weights = (w1 * w2 * w3 * w4).ravel()
    snr_linear = 10.0 ** (float(snr_db) / 10.0)
    first_mi = _mi_from_power(
        np.abs(z1) ** 2,
        snr_linear,
        iqam_db_grid,
        iqam_values,
    )
    mean_mi = float(np.sum(product_weights * first_mi))
    rho_grid = np.arange(0.0, 1.0 + 0.5 * float(rho_step), float(rho_step))
    covariance = np.empty_like(rho_grid)
    for index, rho in enumerate(rho_grid):
        z2 = float(rho) * z1 + math.sqrt(max(0.0, 1.0 - float(rho) ** 2)) * independent
        second_mi = _mi_from_power(
            np.abs(z2) ** 2,
            snr_linear,
            iqam_db_grid,
            iqam_values,
        )
        covariance[index] = (
            float(np.sum(product_weights * first_mi * second_mi)) - mean_mi * mean_mi
        )
    covariance[0] = 0.0
    covariance = np.maximum.accumulate(np.maximum(covariance, 0.0))
    return CDDKernel(float(snr_db), rho_grid, covariance)


def jcdd_metrics(
    delay_indices: Sequence[int],
    physical_frequency_covariance: np.ndarray,
    kernel: CDDKernel,
) -> dict[str, float]:
    covariance = np.asarray(physical_frequency_covariance, dtype=np.complex128)
    period = int(covariance.shape[0])
    if covariance.shape != (period, period):
        raise ValueError("Physical frequency covariance must be square.")
    diagonal = float(np.real(np.mean(np.diag(covariance))))
    if diagonal <= 0.0:
        raise ValueError("Physical frequency covariance must have positive diagonal.")
    factor = array_factor(delay_indices, period)
    physical = covariance[0, :] / diagonal
    rho = np.abs(physical[1:] * factor[1:])
    weights = 2.0 * (period - np.arange(1, period, dtype=np.float64))
    return {
        "jcdd": float(np.sum(weights * kernel.evaluate(rho))),
        "m2_eff": float(np.sum(weights * rho**2)),
        "m4_eff": float(np.sum(weights * rho**4)),
        "max_effective_correlation": float(np.max(rho)),
        "mean_effective_correlation": float(np.mean(rho)),
    }


class FrequencyCEMetric:
    """Matched full-band frequency-LMMSE NMSE for the static two-DMRS average."""

    def __init__(
        self,
        physical_frequency_covariance: np.ndarray,
        pilot_local_indices: Sequence[int],
        data_local_indices: Sequence[int],
        n_tx: int = 8,
    ) -> None:
        covariance = np.asarray(physical_frequency_covariance, dtype=np.complex128)
        if covariance.ndim != 2 or covariance.shape[0] != covariance.shape[1]:
            raise ValueError("Physical frequency covariance must be square.")
        self.covariance = covariance
        self.period = int(covariance.shape[0])
        self.pilots = np.asarray(pilot_local_indices, dtype=np.int64).reshape(-1)
        data = np.asarray(data_local_indices, dtype=np.int64).reshape(-1)
        self.data_counts = np.bincount(data, minlength=self.period).astype(np.float64)
        self.n_tx = int(n_tx)
        rows = np.arange(self.period, dtype=np.int64)
        self.full_pilot_difference = rows[:, None] - self.pilots[None, :]
        self.pilot_difference = self.pilots[:, None] - self.pilots[None, :]
        self.base_r_fp = covariance[:, self.pilots]
        self.base_r_pp = covariance[np.ix_(self.pilots, self.pilots)]
        self.signal_diagonal = (
            np.real(np.diag(covariance)) * float(self.n_tx)
        )

    def evaluate_grid_coordinates(
        self,
        delay_grid_coordinates: Sequence[float],
        snr_db: float,
    ) -> dict[str, float | int]:
        coordinates = np.asarray(delay_grid_coordinates, dtype=np.float64).reshape(-1)
        full_factor = np.mean(
            np.exp(
                -2j
                * np.pi
                * self.full_pilot_difference[:, :, None]
                * coordinates[None, None, :]
                / float(self.period)
            ),
            axis=2,
        ) * float(self.n_tx)
        pilot_factor = np.mean(
            np.exp(
                -2j
                * np.pi
                * self.pilot_difference[:, :, None]
                * coordinates[None, None, :]
                / float(self.period)
            ),
            axis=2,
        ) * float(self.n_tx)
        r_fp = self.base_r_fp * full_factor
        r_pp = self.base_r_pp * pilot_factor
        noise_variance = float(self.n_tx) / (10.0 ** (float(snr_db) / 10.0)) / 2.0
        system = r_pp + noise_variance * np.eye(len(self.pilots), dtype=np.complex128)
        jitter = 0.0
        for _ in range(8):
            try:
                chol = np.linalg.cholesky(system)
                break
            except np.linalg.LinAlgError:
                jitter = max(1e-14, 10.0 * jitter)
                system = r_pp + (noise_variance + jitter) * np.eye(
                    len(self.pilots), dtype=np.complex128
                )
        else:
            raise np.linalg.LinAlgError("Matched frequency-LMMSE system is not positive definite.")
        solved = np.linalg.solve(chol, r_fp.conj().T)
        solved = np.linalg.solve(chol.conj().T, solved)
        gain = np.real(np.sum(r_fp * solved.T, axis=1))
        mse = np.maximum(self.signal_diagonal - gain, 0.0)
        signal = float(np.sum(self.data_counts * self.signal_diagonal))
        nmse = float(np.sum(self.data_counts * mse) / signal)
        singular = np.linalg.svd(system, compute_uv=False)
        pilot_phase = np.exp(
            -2j
            * np.pi
            * self.pilots[:, None]
            * coordinates[None, :]
            / float(self.period)
        )
        pilot_singular = np.linalg.svd(pilot_phase, compute_uv=False)
        return {
            "snr_db": float(snr_db),
            "nmse": nmse,
            "nmse_db": 10.0 * math.log10(max(nmse, 1e-300)),
            "condition_number": float(singular[0] / singular[-1]),
            "minimum_singular_value": float(singular[-1]),
            "numerical_jitter": float(jitter),
            "pilot_rank": int(np.linalg.matrix_rank(pilot_phase, tol=1e-10)),
            "pilot_condition_number": float(pilot_singular[0] / pilot_singular[-1]),
            "pilot_minimum_singular_value": float(pilot_singular[-1]),
            "r_pp_trace": float(np.real(np.trace(r_pp))),
            "r_pp_frobenius_norm": float(np.linalg.norm(r_pp)),
            "r_fp_frobenius_norm": float(np.linalg.norm(r_fp)),
        }

    def evaluate(self, delay_indices: Sequence[int], snr_db: float) -> dict[str, float | int]:
        return self.evaluate_grid_coordinates(delay_indices, snr_db)
