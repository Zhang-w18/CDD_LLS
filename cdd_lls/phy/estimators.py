from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional
import numpy as np

from cdd_lls.core.config import ChannelEstimationConfig, ResourceConfig
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
