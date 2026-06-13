from __future__ import annotations

from functools import lru_cache
import math
import numpy as np


SUPPORTED_QM = (2, 4, 6, 8)


def _bits_to_int(bits: np.ndarray) -> int:
    out = 0
    for bit in bits:
        out = (out << 1) | int(bit)
    return int(out)


def _int_to_bits(x: int, width: int) -> np.ndarray:
    return np.asarray([(int(x) >> (width - 1 - i)) & 1 for i in range(width)], dtype=np.int8)


def _gray_to_binary(g: int) -> int:
    b = int(g)
    while g > 0:
        g >>= 1
        b ^= g
    return int(b)


@lru_cache(maxsize=None)
def qam_constellation(qm: int) -> tuple[np.ndarray, np.ndarray]:
    qm = int(qm)
    if qm not in SUPPORTED_QM:
        raise ValueError(f"Unsupported Qm={qm}. Supported: {SUPPORTED_QM}")
    axis_bits = qm // 2
    m_axis = 2 ** axis_bits
    levels = np.arange(-(m_axis - 1), m_axis, 2, dtype=np.float64)
    norm = math.sqrt((2.0 / 3.0) * ((2 ** qm) - 1))

    bits = []
    symbols = []
    for i_gray in range(m_axis):
        for q_gray in range(m_axis):
            i_idx = _gray_to_binary(i_gray)
            q_idx = _gray_to_binary(q_gray)
            b = np.concatenate([_int_to_bits(i_gray, axis_bits), _int_to_bits(q_gray, axis_bits)])
            bits.append(b)
            symbols.append((levels[i_idx] + 1j * levels[q_idx]) / norm)
    return np.asarray(bits, dtype=np.int8), np.asarray(symbols, dtype=np.complex128)


def qam_modulate(bits: np.ndarray, qm: int) -> np.ndarray:
    bits = np.asarray(bits, dtype=np.int8).reshape(-1)
    qm = int(qm)
    if bits.size % qm:
        raise ValueError(f"Bit length {bits.size} is not divisible by Qm={qm}.")
    const_bits, const_symbols = qam_constellation(qm)
    axis_bits = qm // 2
    out = np.empty(bits.size // qm, dtype=np.complex128)
    for i, chunk in enumerate(bits.reshape(-1, qm)):
        i_gray = _bits_to_int(chunk[:axis_bits])
        q_gray = _bits_to_int(chunk[axis_bits:])
        out[i] = const_symbols[(i_gray << axis_bits) | q_gray]
    return out


def qam_demapper_maxlog(symbols: np.ndarray, noise_var: np.ndarray | float, qm: int) -> np.ndarray:
    """Max-log LLRs with Sionna sign convention: positive means bit=1."""
    z = np.asarray(symbols, dtype=np.complex128).reshape(-1)
    no = np.asarray(noise_var, dtype=np.float64)
    if no.ndim == 0:
        no = np.full(z.shape, float(no), dtype=np.float64)
    else:
        no = no.reshape(-1)
    if no.size != z.size:
        raise ValueError("noise_var must be scalar or have one value per symbol.")
    no = np.maximum(no, 1e-12)
    const_bits, const_symbols = qam_constellation(int(qm))
    distances = np.abs(z[:, None] - const_symbols[None, :]) ** 2 / no[:, None]
    llr = np.empty((z.size, int(qm)), dtype=np.float64)
    for bit_idx in range(int(qm)):
        mask0 = const_bits[:, bit_idx] == 0
        mask1 = ~mask0
        m0 = np.min(distances[:, mask0], axis=1)
        m1 = np.min(distances[:, mask1], axis=1)
        llr[:, bit_idx] = m0 - m1
    return llr.reshape(-1)
