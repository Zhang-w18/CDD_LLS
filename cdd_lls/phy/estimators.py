from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Dict, List, Optional
import numpy as np

from cdd_lls.core.config import ChannelConfig, ChannelEstimationConfig, ResourceConfig
from cdd_lls.phy.precoding import cdd_equivalent_from_branches, normalize_delay_vector
from cdd_lls.phy.resource_grid import ResourceGrid, local_indices_for_subcarriers


@dataclass
class EstimationResult:
    g_hat: np.ndarray
    h_hat: Optional[np.ndarray] = None
    ce_nmse_eff: float = 0.0
    ce_nmse_branch: float = float("nan")
    cond_number: float = float("nan")
    effective_rank: float = float("nan")
    metadata: Dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class TDLTimeFrequencyCovariance:
    time: np.ndarray
    frequency: np.ndarray
    covariance_type: str


@dataclass(frozen=True)
class TimeFrequencyRMMSEFilter:
    weights: np.ndarray
    pilot_coordinates: np.ndarray
    data_coordinates: np.ndarray
    noise_variance: float
    condition_number: float
    numerical_jitter: float
    covariance_type: str

    def estimate_data(self, ls_observations: np.ndarray) -> np.ndarray:
        obs = np.asarray(ls_observations, dtype=np.complex128)
        if obs.shape[-1] != self.weights.shape[1]:
            raise ValueError(
                f"LS observation has {obs.shape[-1]} pilot REs; expected {self.weights.shape[1]}."
            )
        original = obs.shape[:-1]
        flat = obs.reshape(-1, obs.shape[-1])
        estimated = flat @ self.weights.T
        return estimated.reshape(*original, self.weights.shape[0])


@lru_cache(maxsize=32)
def _cached_tdl_base_time_frequency_covariance(
    profile: str,
    subcarrier_spacing_hz: float,
    n_fft: int,
    delay_spread_s: float,
    speed_mps: float,
    carrier_frequency_hz: float,
    ofdm_symbol_duration_s: float,
    n_symbols: int,
    active_fft_indices: tuple[int, ...],
) -> tuple[np.ndarray, np.ndarray]:
    from sionna.phy.ofdm import tdl_freq_cov_mat, tdl_time_cov_mat

    if profile not in ("A", "B", "C", "D", "E"):
        raise ValueError("Sionna covariance helpers support TDL profiles A through E.")
    rf_full = np.asarray(
        tdl_freq_cov_mat(
            model=profile,
            subcarrier_spacing=float(subcarrier_spacing_hz),
            fft_size=int(n_fft),
            delay_spread=float(delay_spread_s),
            precision="double",
        ).numpy(),
        dtype=np.complex128,
    )
    active = np.asarray(active_fft_indices, dtype=np.int64)
    rf = rf_full[np.ix_(active, active)]
    rt = np.asarray(
        tdl_time_cov_mat(
            model=profile,
            speed=float(speed_mps),
            carrier_frequency=float(carrier_frequency_hz),
            ofdm_symbol_duration=float(ofdm_symbol_duration_s),
            num_ofdm_symbols=int(n_symbols),
            precision="double",
        ).numpy(),
        dtype=np.complex128,
    )
    return rt, rf


def _tdl_base_time_frequency_covariance(
    grid: ResourceGrid,
    channel: ChannelConfig,
) -> tuple[np.ndarray, np.ndarray]:
    return _cached_tdl_base_time_frequency_covariance(
        str(channel.tdl_profile).upper(),
        float(grid.scs_khz) * 1e3,
        int(grid.n_fft),
        float(channel.delay_spread_ns) * 1e-9,
        float(channel.ue_speed_kmh) / 3.6,
        float(channel.carrier_frequency_hz),
        float(grid.ofdm_symbol_duration_s),
        int(grid.n_symbols),
        tuple(int(x) for x in grid.active_fft_indices),
    )


def tdl_unknown_delay_covariance(
    grid: ResourceGrid,
    channel: ChannelConfig,
) -> TDLTimeFrequencyCovariance:
    """Return the baseline TDL covariance without accepting true CDD delays."""
    rt, rf = _tdl_base_time_frequency_covariance(grid, channel)
    return TDLTimeFrequencyCovariance(rt, rf, "tdl_unknown_delay")


def tdl_known_delay_covariance(
    grid: ResourceGrid,
    channel: ChannelConfig,
    delays: List[int],
) -> TDLTimeFrequencyCovariance:
    rt, rf = _tdl_base_time_frequency_covariance(grid, channel)
    normalized = normalize_delay_vector(delays, n_tx=len(delays))
    k = np.asarray(grid.subcarrier_indices, dtype=np.float64)
    delta = k[:, None] - k[None, :]
    d = np.asarray(normalized, dtype=np.float64)
    cdd_factor = np.mean(
        np.exp(-1j * 2.0 * np.pi * delta[:, :, None] * d[None, None, :] / float(grid.n_fft)),
        axis=2,
    )
    return TDLTimeFrequencyCovariance(rt, rf * cdd_factor, "tdl_cdd_known_delay")


def coordinate_covariance(
    grid: ResourceGrid,
    time_covariance: np.ndarray,
    frequency_covariance: np.ndarray,
    coordinates_a: np.ndarray,
    coordinates_b: np.ndarray,
) -> np.ndarray:
    a = np.asarray(coordinates_a, dtype=np.int64)
    b = np.asarray(coordinates_b, dtype=np.int64)
    if a.ndim != 2 or b.ndim != 2 or a.shape[1] != 2 or b.shape[1] != 2:
        raise ValueError("Coordinates must have shape [n_re,2] as (symbol,subcarrier).")
    a_local = local_indices_for_subcarriers(grid, a[:, 1])
    b_local = local_indices_for_subcarriers(grid, b[:, 1])
    return (
        np.asarray(time_covariance)[np.ix_(a[:, 0], b[:, 0])]
        * np.asarray(frequency_covariance)[np.ix_(a_local, b_local)]
    )


def build_time_frequency_rmmse_filter(
    grid: ResourceGrid,
    covariance: TDLTimeFrequencyCovariance,
    noise_variance: float,
    diagonal_loading: float = 0.0,
) -> TimeFrequencyRMMSEFilter:
    pilot = np.asarray(grid.pilot_coordinates, dtype=np.int64)
    data = np.asarray(grid.data_coordinates, dtype=np.int64)
    r_pp = coordinate_covariance(grid, covariance.time, covariance.frequency, pilot, pilot)
    r_dp = coordinate_covariance(grid, covariance.time, covariance.frequency, data, pilot)
    jitter = max(float(diagonal_loading), 0.0)
    identity = np.eye(len(pilot), dtype=np.complex128)
    system = r_pp + (float(noise_variance) + jitter) * identity
    for _ in range(8):
        try:
            chol = np.linalg.cholesky(system)
            break
        except np.linalg.LinAlgError:
            jitter = max(1e-14, 10.0 * jitter)
            system = r_pp + (float(noise_variance) + jitter) * identity
    else:
        raise np.linalg.LinAlgError("RMMSE pilot covariance is not numerically positive definite.")
    solved = np.linalg.solve(chol, r_dp.conj().T)
    solved = np.linalg.solve(chol.conj().T, solved)
    weights = solved.conj().T
    return TimeFrequencyRMMSEFilter(
        weights=weights,
        pilot_coordinates=pilot,
        data_coordinates=data,
        noise_variance=float(noise_variance),
        condition_number=float(np.linalg.cond(system)),
        numerical_jitter=float(jitter),
        covariance_type=str(covariance.covariance_type),
    )


def linear_estimator_closed_form_nmse(
    grid: ResourceGrid,
    estimator: TimeFrequencyRMMSEFilter,
    true_covariance: TDLTimeFrequencyCovariance,
    noise_variance: float,
) -> float:
    pilot = np.asarray(grid.pilot_coordinates, dtype=np.int64)
    data = np.asarray(grid.data_coordinates, dtype=np.int64)
    r_pp = coordinate_covariance(grid, true_covariance.time, true_covariance.frequency, pilot, pilot)
    r_dp = coordinate_covariance(grid, true_covariance.time, true_covariance.frequency, data, pilot)
    w = np.asarray(estimator.weights, dtype=np.complex128)
    data_t_diag = np.real(np.diag(true_covariance.time)[data[:, 0]])
    data_k_local = local_indices_for_subcarriers(grid, data[:, 1])
    data_f_diag = np.real(np.diag(true_covariance.frequency)[data_k_local])
    signal_power = float(np.sum(data_t_diag * data_f_diag))
    cross = np.sum(w * r_dp.conj())
    observation_covariance = r_pp + float(noise_variance) * np.eye(len(pilot), dtype=np.complex128)
    estimate_power = np.sum((w @ observation_covariance) * w.conj())
    mse = signal_power - 2.0 * float(np.real(cross)) + float(np.real(estimate_power))
    return mse / signal_power


def shifted_pdp(pdp: np.ndarray, delays: List[int]) -> np.ndarray:
    p = np.asarray(pdp, dtype=np.float64).reshape(-1)
    d = [int(x) for x in delays]
    if any(x < 0 for x in d):
        raise ValueError("Only non-negative CDD delays are supported in shifted_pdp.")
    out = np.zeros(len(p) + (max(d) if d else 0), dtype=np.float64)
    for delay in d:
        out[delay:delay + len(p)] += p / max(len(d), 1)
    s = float(np.sum(out))
    return out / s if s > 0 else out


def covariance_matrix(k_a: np.ndarray, k_b: np.ndarray, pdp: np.ndarray, n_fft: int) -> np.ndarray:
    ka = np.asarray(k_a, dtype=np.float64).reshape(-1)
    kb = np.asarray(k_b, dtype=np.float64).reshape(-1)
    p = np.asarray(pdp, dtype=np.float64).reshape(-1)
    taps = np.arange(len(p), dtype=np.float64)
    delta = ka[:, None] - kb[None, :]
    phase = np.exp(-1j * 2.0 * np.pi * delta[:, :, None] * taps[None, None, :] / float(n_fft))
    return np.tensordot(phase, p, axes=([2], [0]))


def _solve_lmmse(
    obs: np.ndarray,
    target_k: np.ndarray,
    pilot_k: np.ndarray,
    pdp: np.ndarray,
    n_fft: int,
    noise_var: float,
    loading: float,
) -> np.ndarray:
    if len(pilot_k) == 0:
        raise ValueError("At least one pilot is required for LMMSE estimation.")
    R_pp = covariance_matrix(pilot_k, pilot_k, pdp, n_fft)
    R_tp = covariance_matrix(target_k, pilot_k, pdp, n_fft)
    A = R_pp + (float(noise_var) + float(loading)) * np.eye(len(pilot_k), dtype=np.complex128)
    weights_t = np.linalg.solve(A, np.asarray(obs, dtype=np.complex128).T)
    return (R_tp @ weights_t).T


def linear_interpolate(ls_obs: np.ndarray, pilot_k: np.ndarray, target_k: np.ndarray) -> np.ndarray:
    obs = np.asarray(ls_obs, dtype=np.complex128)
    pk = np.asarray(pilot_k, dtype=np.float64)
    tk = np.asarray(target_k, dtype=np.float64)
    out = np.empty((obs.shape[0], len(tk)), dtype=np.complex128)
    for r in range(obs.shape[0]):
        real = np.interp(tk, pk, np.real(obs[r]))
        imag = np.interp(tk, pk, np.imag(obs[r]))
        out[r] = real + 1j * imag
    return out


def direct_rmmse(
    ls_obs: np.ndarray,
    grid: ResourceGrid,
    resource: ResourceConfig,
    pdp_assumed: np.ndarray,
    noise_var_ls: float,
    bundle_rb: Optional[int],
    loading: float,
) -> np.ndarray:
    target_k = grid.subcarrier_indices
    pilot_k = grid.pilot_subcarriers
    if bundle_rb is None:
        return _solve_lmmse(
            obs=ls_obs,
            target_k=target_k,
            pilot_k=pilot_k,
            pdp=pdp_assumed,
            n_fft=grid.n_fft,
            noise_var=noise_var_ls,
            loading=loading,
        )

    out = np.empty((ls_obs.shape[0], grid.n_sc), dtype=np.complex128)
    pilot_local = local_indices_for_subcarriers(grid, pilot_k)
    bundle_sc = int(bundle_rb) * 12
    for start in range(0, grid.n_sc, bundle_sc):
        stop = min(grid.n_sc, start + bundle_sc)
        target_sel = np.arange(start, stop, dtype=np.int64)
        pilot_mask = (pilot_local >= start) & (pilot_local < stop)
        if not np.any(pilot_mask):
            out[:, target_sel] = linear_interpolate(ls_obs, pilot_k, target_k[target_sel])
            continue
        out[:, target_sel] = _solve_lmmse(
            obs=ls_obs[:, pilot_mask],
            target_k=target_k[target_sel],
            pilot_k=pilot_k[pilot_mask],
            pdp=pdp_assumed,
            n_fft=grid.n_fft,
            noise_var=noise_var_ls,
            loading=loading,
        )
    return out


def reconstruction_pairwise(
    ls_obs: np.ndarray,
    grid: ResourceGrid,
    delays: List[int],
    reg: float,
    spacing_pilots: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    n_rx = int(ls_obs.shape[0])
    n_tx = len(delays)
    pilot_k = grid.pilot_subcarriers
    spacing = max(int(spacing_pilots), 1)
    group_len = n_tx
    anchors = []
    h_estimates = []
    conds = []
    d = np.asarray(delays, dtype=np.float64)

    max_start = len(pilot_k) - (group_len - 1) * spacing
    for start in range(max_start):
        idx = start + spacing * np.arange(group_len)
        pk = pilot_k[idx]
        A = np.exp(
            -1j * 2.0 * np.pi * pk[:, None] * d[None, :] / float(grid.n_fft)
        ) / np.sqrt(float(n_tx))
        conds.append(float(np.linalg.cond(A)))
        lhs = A.conj().T @ A + float(reg) * np.eye(n_tx, dtype=np.complex128)
        rhs = A.conj().T
        local = np.empty((n_rx, n_tx), dtype=np.complex128)
        for r in range(n_rx):
            local[r] = np.linalg.solve(lhs, rhs @ ls_obs[r, idx])
        anchors.append(float(np.mean(pk)))
        h_estimates.append(local)

    if not anchors:
        raise ValueError("Not enough pilots for pairwise reconstruction.")

    anchors_arr = np.asarray(anchors, dtype=np.float64)
    h_arr = np.asarray(h_estimates, dtype=np.complex128)  # [n_anchor,n_rx,n_tx]
    order = np.argsort(anchors_arr)
    anchors_arr = anchors_arr[order]
    h_arr = h_arr[order]

    h_hat = np.empty((n_rx, n_tx, grid.n_sc), dtype=np.complex128)
    tk = grid.subcarrier_indices.astype(np.float64)
    for r in range(n_rx):
        for m in range(n_tx):
            vals = h_arr[:, r, m]
            real = np.interp(tk, anchors_arr, np.real(vals))
            imag = np.interp(tk, anchors_arr, np.imag(vals))
            h_hat[r, m] = real + 1j * imag
    g_hat = cdd_equivalent_from_branches(h_hat, grid, delays)
    return g_hat, h_hat, float(np.nanmax(conds))


def _basis_support_indices(pdp: np.ndarray, mode: str, threshold: float) -> np.ndarray:
    p = np.asarray(pdp, dtype=np.float64).reshape(-1)
    if str(mode).lower() == "ideal":
        return np.arange(len(p), dtype=np.int64)
    c = np.cumsum(p)
    cutoff = int(np.searchsorted(c, float(threshold), side="left")) + 1
    cutoff = min(max(cutoff, 1), len(p))
    return np.arange(cutoff, dtype=np.int64)


def reconstruction_basis_lmmse(
    ls_obs: np.ndarray,
    grid: ResourceGrid,
    delays: List[int],
    pdp: np.ndarray,
    noise_var_ls: float,
    support_mode: str,
    energy_threshold: float,
    loading: float,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    n_rx = int(ls_obs.shape[0])
    n_tx = len(delays)
    support = _basis_support_indices(pdp, support_mode, energy_threshold)
    p_support = np.asarray(pdp, dtype=np.float64)[support]
    d = np.asarray(delays, dtype=np.float64)
    pilot_k = grid.pilot_subcarriers.astype(np.float64)
    n_cols = n_tx * len(support)
    Phi = np.empty((len(pilot_k), n_cols), dtype=np.complex128)
    col = 0
    for m in range(n_tx):
        shifted = support.astype(np.float64) + d[m]
        Phi[:, col:col + len(support)] = np.exp(
            -1j * 2.0 * np.pi * pilot_k[:, None] * shifted[None, :] / float(grid.n_fft)
        ) / np.sqrt(float(n_tx))
        col += len(support)

    rdiag = np.tile(p_support, n_tx)
    PhiR = Phi * rdiag[None, :]
    S = PhiR @ Phi.conj().T
    S += (float(noise_var_ls) + float(loading)) * np.eye(S.shape[0], dtype=np.complex128)
    B = (rdiag[:, None] * Phi.conj().T)

    h_taps = np.empty((n_rx, n_cols), dtype=np.complex128)
    for r in range(n_rx):
        h_taps[r] = B @ np.linalg.solve(S, ls_obs[r])

    tk = grid.subcarrier_indices.astype(np.float64)
    F = np.exp(-1j * 2.0 * np.pi * tk[:, None] * support[None, :] / float(grid.n_fft))
    h_hat = np.empty((n_rx, n_tx, grid.n_sc), dtype=np.complex128)
    col = 0
    for m in range(n_tx):
        taps_m = h_taps[:, col:col + len(support)]
        h_hat[:, m, :] = taps_m @ F.T
        col += len(support)

    g_hat = cdd_equivalent_from_branches(h_hat, grid, delays)
    weighted = Phi * np.sqrt(np.maximum(rdiag, 0.0))[None, :]
    svals = np.linalg.svd(weighted, compute_uv=False)
    eps = 1e-12
    cond = float(svals[0] / max(svals[-1], eps)) if len(svals) else float("nan")
    eff_rank = float((np.sum(svals) ** 2) / max(np.sum(svals ** 2), eps)) if len(svals) else float("nan")
    return g_hat, h_hat, cond, eff_rank


def nmse(true: np.ndarray, est: np.ndarray) -> float:
    denom = float(np.sum(np.abs(true) ** 2))
    if denom <= 0:
        return float("nan")
    return float(np.sum(np.abs(true - est) ** 2) / denom)


def estimate_channel(
    method: str,
    tx_scheme: str,
    ls_obs: np.ndarray,
    true_g: np.ndarray,
    true_H: np.ndarray,
    grid: ResourceGrid,
    resource: ResourceConfig,
    ce_cfg: ChannelEstimationConfig,
    pdp: np.ndarray,
    delays: List[int],
    noise_var_ls: float,
) -> EstimationResult:
    method_u = str(method).upper()
    tx_u = str(tx_scheme).upper()
    n_tx = int(true_H.shape[1])
    delays = normalize_delay_vector(delays, n_tx=n_tx)

    if method_u in ("IDEAL", "IDEAL_CSI"):
        return EstimationResult(
            g_hat=true_g.copy(),
            ce_nmse_eff=0.0,
            ce_nmse_branch=0.0,
            metadata={"ce_method": "IDEAL_CSI"},
        )

    if method_u == "LS_LINEAR":
        g_hat = linear_interpolate(ls_obs, grid.pilot_subcarriers, grid.subcarrier_indices)
        return EstimationResult(
            g_hat=g_hat,
            ce_nmse_eff=nmse(true_g, g_hat),
            metadata={"ce_method": method_u},
        )

    if method_u in ("PRG_RMMSE_4RB", "RMMSE_4RB_KNOWN", "RMMSE_4RB_UNKNOWN", "RMMSE_WB_KNOWN", "RMMSE_WB_UNKNOWN"):
        known = method_u.endswith("KNOWN") and not method_u.endswith("UNKNOWN")
        if tx_u.startswith("PRG"):
            assumed = np.asarray(pdp, dtype=np.float64)
            bundle = int(resource.prg_size_rb)
        else:
            assumed = shifted_pdp(pdp, delays) if known else np.asarray(pdp, dtype=np.float64)
            bundle = int(ce_cfg.rmmse_bundle_rb) if "_4RB_" in method_u or method_u == "PRG_RMMSE_4RB" else None
        g_hat = direct_rmmse(
            ls_obs=ls_obs,
            grid=grid,
            resource=resource,
            pdp_assumed=assumed,
            noise_var_ls=float(noise_var_ls),
            bundle_rb=bundle,
            loading=float(ce_cfg.diagonal_loading),
        )
        return EstimationResult(
            g_hat=g_hat,
            ce_nmse_eff=nmse(true_g, g_hat),
            metadata={
                "ce_method": method_u,
                "rmmse_processing": "4rb" if bundle else "wideband",
                "covariance": "cdd_shifted" if (known and not tx_u.startswith("PRG")) else "non_cdd",
            },
        )

    if method_u == "RECON_PAIRWISE":
        g_hat, h_hat, cond = reconstruction_pairwise(
            ls_obs=ls_obs,
            grid=grid,
            delays=delays,
            reg=float(ce_cfg.recon_regularization),
            spacing_pilots=int(ce_cfg.recon_pair_spacing_pilots),
        )
        return EstimationResult(
            g_hat=g_hat,
            h_hat=h_hat,
            ce_nmse_eff=nmse(true_g, g_hat),
            ce_nmse_branch=nmse(true_H, h_hat),
            cond_number=cond,
            metadata={"ce_method": method_u},
        )

    if method_u == "RECON_BASIS_LMMSE":
        g_hat, h_hat, cond, eff_rank = reconstruction_basis_lmmse(
            ls_obs=ls_obs,
            grid=grid,
            delays=delays,
            pdp=pdp,
            noise_var_ls=float(noise_var_ls),
            support_mode=str(ce_cfg.basis_support),
            energy_threshold=float(ce_cfg.basis_energy_threshold),
            loading=float(ce_cfg.diagonal_loading),
        )
        return EstimationResult(
            g_hat=g_hat,
            h_hat=h_hat,
            ce_nmse_eff=nmse(true_g, g_hat),
            ce_nmse_branch=nmse(true_H, h_hat),
            cond_number=cond,
            effective_rank=eff_rank,
            metadata={"ce_method": method_u, "basis_support": str(ce_cfg.basis_support)},
        )

    raise ValueError(f"Unsupported CE method={method}.")
