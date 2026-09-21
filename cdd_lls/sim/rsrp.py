from __future__ import annotations

import numpy as np


def _selected_effective_channel(
    channel: np.ndarray,
    precoder: np.ndarray,
    symbol_indices: np.ndarray,
    subcarrier_indices: np.ndarray,
    *,
    n_rx: int,
) -> tuple[np.ndarray, int]:
    """Return selected noiseless effective channels and the selected RE count."""

    h = np.asarray(channel, dtype=np.complex128)
    w = np.asarray(precoder, dtype=np.complex128)
    symbols = np.asarray(symbol_indices, dtype=np.int64).reshape(-1)
    subcarriers = np.asarray(subcarrier_indices, dtype=np.int64).reshape(-1)
    if h.ndim != 5:
        raise ValueError("channel must have shape [batch,rx,tx,symbol,subcarrier].")
    if symbols.size == 0 or symbols.shape != subcarriers.shape:
        raise ValueError("symbol_indices and subcarrier_indices must be non-empty and aligned.")
    if int(n_rx) != h.shape[1] or int(n_rx) <= 0:
        raise ValueError("n_rx must equal the channel receive-antenna dimension.")
    if np.any(symbols < 0) or np.any(symbols >= h.shape[3]):
        raise ValueError("symbol_indices lie outside the channel array.")
    if np.any(subcarriers < 0) or np.any(subcarriers >= h.shape[4]):
        raise ValueError("subcarrier_indices lie outside the channel array.")

    if w.ndim == 2 and w.shape == (h.shape[4], h.shape[2]):
        selected_w = w[subcarriers, :]
    elif w.ndim == 2 and w.shape == (symbols.size, h.shape[2]):
        selected_w = w
    elif w.ndim == 3 and w.shape == (h.shape[3], h.shape[4], h.shape[2]):
        selected_w = w[symbols, subcarriers, :]
    else:
        raise ValueError(
            "precoder must have shape [subcarrier,tx], [re,tx], or "
            "[symbol,subcarrier,tx]."
        )

    selected_h = h[:, :, :, symbols, subcarriers]
    effective = np.einsum("brmq,qm->brq", selected_h, selected_w, optimize=True)
    return effective, int(symbols.size)


def noiseless_block_average_power_per_rx(
    channel: np.ndarray,
    precoder: np.ndarray,
    symbol_indices: np.ndarray,
    subcarrier_indices: np.ndarray,
    *,
    n_rx: int,
    transmit_power: float = 1.0,
) -> np.ndarray:
    """Return normalized noiseless block-average power per trial and Rx.

    The result has shape ``[batch, rx]``. No sum or average is performed over
    receive antennas, and no per-trial normalization is applied.
    """

    if not np.isfinite(transmit_power) or float(transmit_power) <= 0.0:
        raise ValueError("transmit_power must be finite and positive.")
    effective, n_re = _selected_effective_channel(
        channel,
        precoder,
        symbol_indices,
        subcarrier_indices,
        n_rx=n_rx,
    )
    return (
        np.sum(np.abs(effective) ** 2, axis=2)
        / (float(n_re) * float(transmit_power))
    ).astype(np.float64)


def noiseless_block_average_power(
    channel: np.ndarray,
    precoder: np.ndarray,
    symbol_indices: np.ndarray,
    subcarrier_indices: np.ndarray,
    *,
    n_rx: int,
    transmit_power: float = 1.0,
) -> np.ndarray:
    """Return the normalized noiseless block-average receive power per trial.

    ``channel`` has shape ``[batch, rx, tx, symbol, subcarrier]``. ``precoder``
    may be indexed by subcarrier (``[subcarrier, tx]``), by selected RE
    (``[re, tx]``), or by symbol and subcarrier
    (``[symbol, subcarrier, tx]``). No per-trial normalization is applied.
    """

    if not np.isfinite(transmit_power) or float(transmit_power) <= 0.0:
        raise ValueError("transmit_power must be finite and positive.")
    effective, n_re = _selected_effective_channel(
        channel,
        precoder,
        symbol_indices,
        subcarrier_indices,
        n_rx=n_rx,
    )
    return (
        np.sum(np.abs(effective) ** 2, axis=(1, 2))
        / (float(n_rx) * float(n_re) * float(transmit_power))
    ).astype(np.float64)
