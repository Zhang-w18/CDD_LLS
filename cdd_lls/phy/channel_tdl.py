from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from cdd_lls.core.config import ChannelConfig
from cdd_lls.phy.resource_grid import ResourceGrid


@dataclass(frozen=True)
class TDLRealization:
    taps: np.ndarray
    H: np.ndarray
    pdp: np.ndarray
    tap_delays: np.ndarray
    sample_period_ns: float


def make_exponential_pdp(
    delay_spread_ns: float,
    sample_period_ns: float,
    max_delay_factor: float = 8.0,
) -> np.ndarray:
    tau = max(float(delay_spread_ns), 1e-6)
    ts = float(sample_period_ns)
    max_delay_ns = max(float(max_delay_factor) * tau, ts)
    n_taps = max(1, int(math.ceil(max_delay_ns / ts)) + 1)
    delays = np.arange(n_taps, dtype=np.float64) * ts
    pdp = np.exp(-delays / tau)
    pdp /= np.sum(pdp)
    return pdp.astype(np.float64)


def generate_tdl_channel(
    rng: np.random.Generator,
    grid: ResourceGrid,
    channel: ChannelConfig,
    n_tx: int,
    n_rx: int,
) -> TDLRealization:
    sample_rate_hz = float(grid.n_fft) * float(grid.scs_khz) * 1e3
    sample_period_ns = 1e9 / sample_rate_hz
    pdp = make_exponential_pdp(
        delay_spread_ns=float(channel.delay_spread_ns),
        sample_period_ns=sample_period_ns,
        max_delay_factor=float(channel.max_delay_factor),
    )
    if not bool(channel.normalize):
        pdp = pdp * len(pdp)

    sigma = np.sqrt(pdp / 2.0)
    real = rng.normal(size=(int(n_rx), int(n_tx), len(pdp)))
    imag = rng.normal(size=(int(n_rx), int(n_tx), len(pdp)))
    taps = (real + 1j * imag) * sigma[None, None, :]

    tap_delays = np.arange(len(pdp), dtype=np.float64)
    phase = np.exp(
        -1j
        * 2.0
        * np.pi
        * grid.subcarrier_indices[:, None]
        * tap_delays[None, :]
        / float(grid.n_fft)
    )
    H = np.einsum("rml,kl->rmk", taps, phase, optimize=True)
    return TDLRealization(
        taps=taps,
        H=H,
        pdp=pdp,
        tap_delays=tap_delays,
        sample_period_ns=float(sample_period_ns),
    )
