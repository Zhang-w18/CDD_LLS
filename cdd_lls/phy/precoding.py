from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence
import numpy as np

from cdd_lls.core.config import ResourceConfig, TransmissionConfig
from cdd_lls.phy.resource_grid import ResourceGrid


@dataclass(frozen=True)
class PrecoderResult:
    C: np.ndarray
    label: str
    metadata: Dict[str, object]


def _canonical_delay(value: float) -> int | float:
    number = float(value)
    if not np.isfinite(number):
        raise ValueError("CDD delays must be finite.")
    return int(number) if number.is_integer() else number


def normalize_delay_vector(
    delays: Sequence[float] | None,
    n_tx: int,
    base_delay: float = 0,
) -> List[int | float]:
    if delays is None:
        return [_canonical_delay(float(i) * float(base_delay)) for i in range(int(n_tx))]
    out = [_canonical_delay(x) for x in delays]
    if len(out) == int(n_tx):
        return out
    if len(out) == 1:
        return [_canonical_delay(float(i) * float(out[0])) for i in range(int(n_tx))]
    raise ValueError(f"CDD delay vector length {len(out)} does not match n_tx={n_tx}.")


def build_active_dft_grid_precoder(
    grid: ResourceGrid,
    delay_indices: Sequence[float],
    n_tx: int,
    normalize: bool = False,
) -> PrecoderResult:
    """Build a CDD precoder on the active-band DFT grid.

    The phase reference is the first active subcarrier, matching the plan-024
    definition exp(-j*2*pi*m*j/K), where m is the local active-band index.
    """
    indices = normalize_delay_vector(delay_indices, n_tx=int(n_tx))
    j = np.asarray(indices, dtype=np.float64)
    local = np.arange(int(grid.n_sc), dtype=np.float64)
    C = np.exp(-1j * 2.0 * np.pi * local[:, None] * j[None, :] / float(grid.n_sc))
    if bool(normalize):
        C = C / np.sqrt(float(n_tx))
    fft_samples = j * float(grid.n_fft) / float(grid.n_sc)
    delays_s = j / (float(grid.n_sc) * float(grid.scs_khz) * 1e3)
    return PrecoderResult(
        C=C.astype(np.complex128),
        label="CDD_ACTIVE_DFT_GRID",
        metadata={
            "cdd_delay_unit": "active_dft_grid_index",
            "cdd_delay_indices": [_canonical_delay(x) for x in j],
            "cdd_delay_fft_samples": [float(x) for x in fft_samples],
            "cdd_delay_seconds": [float(x) for x in delays_s],
            "phase_reference": "first_active_subcarrier",
            "phase_denominator": int(grid.n_sc),
            "normalized": bool(normalize),
        },
    )


def qpsk_codebook(n_tx: int) -> np.ndarray:
    n_tx = int(n_tx)
    if n_tx == 1:
        return np.ones((1, 1), dtype=np.complex128)
    if n_tx == 2:
        return np.asarray(
            [
                [1, 1],
                [1, 1j],
                [1, -1],
                [1, -1j],
            ],
            dtype=np.complex128,
        ) / np.sqrt(2.0)
    if n_tx == 4:
        return np.asarray(
            [
                [1, 1, 1, 1],
                [1, 1j, -1, -1j],
                [1, -1, 1, -1],
                [1, -1j, -1, 1j],
            ],
            dtype=np.complex128,
        ) / 2.0
    raise ValueError(f"Unsupported PRG codebook n_tx={n_tx}.")


def build_precoder(
    grid: ResourceGrid,
    resource: ResourceConfig,
    tx: TransmissionConfig,
    n_tx: int,
) -> PrecoderResult:
    scheme = str(tx.tx_scheme).upper()
    n_tx = int(n_tx)
    if scheme in ("CDD", "NO_CDD"):
        delays = [0] * n_tx if scheme == "NO_CDD" else normalize_delay_vector(
            tx.cdd_delay_vector,
            n_tx=n_tx,
            base_delay=float(tx.cdd_base_delay),
        )
        d = np.asarray(delays, dtype=np.float64)
        C = np.exp(
            -1j
            * 2.0
            * np.pi
            * grid.subcarrier_indices[:, None]
            * d[None, :]
            / float(grid.n_fft)
        ) / np.sqrt(float(n_tx))
        return PrecoderResult(
            C=C.astype(np.complex128),
            label=scheme,
            metadata={
                "cdd_delay_vector": delays,
                "cdd_delay_seconds": [
                    float(x) / (float(grid.n_fft) * float(grid.scs_khz) * 1e3)
                    for x in delays
                ],
            },
        )

    if scheme in ("PRG_CYCLING_4RB", "PRG_CYCLING"):
        codebook = qpsk_codebook(n_tx)
        order = tx.prg_cycling_order
        if order is None:
            order = list(range(codebook.shape[0]))
        order = [int(x) for x in order]
        prg_size_rb = int(resource.prg_size_rb)
        C = np.empty((grid.n_sc, n_tx), dtype=np.complex128)
        for local_sc in range(grid.n_sc):
            rb = int(local_sc // 12)
            prg = int(rb // prg_size_rb)
            cb_idx = order[prg % len(order)] % codebook.shape[0]
            C[local_sc, :] = codebook[cb_idx, :]
        return PrecoderResult(
            C=C,
            label="PRG_CYCLING_4RB",
            metadata={
                "prg_size_rb": prg_size_rb,
                "prg_codebook": str(tx.prg_codebook),
                "prg_cycling_order": order,
            },
        )

    raise ValueError(f"Unsupported tx_scheme={tx.tx_scheme}.")


def equivalent_channel(H: np.ndarray, C: np.ndarray) -> np.ndarray:
    H_arr = np.asarray(H, dtype=np.complex128)
    C_arr = np.asarray(C, dtype=np.complex128)
    if H_arr.ndim == 3:
        if C_arr.shape != (H_arr.shape[2], H_arr.shape[1]):
            raise ValueError("C must have shape [n_sc,n_tx].")
        return np.einsum("rmk,km->rk", H_arr, C_arr, optimize=True)
    if H_arr.ndim == 4:
        if C_arr.shape != (H_arr.shape[3], H_arr.shape[1]):
            raise ValueError("C must have shape [n_sc,n_tx].")
        return np.einsum("rmsk,km->rsk", H_arr, C_arr, optimize=True)
    if H_arr.ndim == 5:
        if C_arr.shape != (H_arr.shape[4], H_arr.shape[2]):
            raise ValueError("C must have shape [n_sc,n_tx].")
        return np.einsum("brmsk,km->brsk", H_arr, C_arr, optimize=True)
    raise ValueError("H must have shape [rx,tx,sc], [rx,tx,symbol,sc], or [batch,rx,tx,symbol,sc].")


def cdd_equivalent_from_branches(H: np.ndarray, grid: ResourceGrid, delays: Sequence[float]) -> np.ndarray:
    h_arr = np.asarray(H)
    n_tx = int(h_arr.shape[2] if h_arr.ndim == 5 else h_arr.shape[1])
    d = np.asarray(normalize_delay_vector(delays, n_tx=n_tx), dtype=np.float64)
    C = np.exp(
        -1j
        * 2.0
        * np.pi
        * grid.subcarrier_indices[:, None]
        * d[None, :]
        / float(grid.n_fft)
    ) / np.sqrt(float(n_tx))
    return equivalent_channel(H, C)
