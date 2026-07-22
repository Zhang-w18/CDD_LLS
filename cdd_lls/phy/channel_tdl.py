from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math
from typing import Dict
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
    backend: str = "legacy_exponential"
    metadata: Dict[str, object] | None = None


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


def _generate_legacy_exponential_channel(
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
        backend="legacy_exponential",
        metadata={"normalization": "unit_ensemble_pdp" if bool(channel.normalize) else "legacy_unnormalized"},
    )


def speed_kmh_to_mps(speed_kmh: float) -> float:
    return float(speed_kmh) / 3.6


@lru_cache(maxsize=32)
def _sionna_generator(
    profile: str,
    delay_spread_s: float,
    carrier_frequency_hz: float,
    speed_mps: float,
    num_sinusoids: int,
    n_rx: int,
    n_tx: int,
    n_symbols: int,
    n_fft: int,
    subcarrier_spacing_hz: float,
    cyclic_prefix_length: int,
):
    from sionna.phy.channel import GenerateOFDMChannel
    from sionna.phy.channel.tr38901 import TDL
    from sionna.phy.ofdm import ResourceGrid as SionnaResourceGrid

    sionna_grid = SionnaResourceGrid(
        num_ofdm_symbols=int(n_symbols),
        fft_size=int(n_fft),
        subcarrier_spacing=float(subcarrier_spacing_hz),
        num_tx=1,
        num_streams_per_tx=1,
        cyclic_prefix_length=int(cyclic_prefix_length),
        precision="double",
    )
    tdl = TDL(
        model=str(profile).upper(),
        delay_spread=float(delay_spread_s),
        carrier_frequency=float(carrier_frequency_hz),
        num_sinusoids=int(num_sinusoids),
        min_speed=float(speed_mps),
        max_speed=float(speed_mps),
        num_rx_ant=int(n_rx),
        num_tx_ant=int(n_tx),
        precision="double",
    )
    generator = GenerateOFDMChannel(
        tdl,
        sionna_grid,
        normalize_channel=False,
        precision="double",
    )
    return generator, tdl


def generate_sionna_tdl_channel(
    grid: ResourceGrid,
    channel: ChannelConfig,
    n_tx: int,
    n_rx: int,
    batch_size: int = 1,
    seed: int = 0,
) -> TDLRealization:
    """Generate Sionna TDL responses as [batch,rx,tx,symbol,active_sc]."""
    import tensorflow as tf
    from sionna.phy import config as sionna_config

    seed = int(seed) % (2**31 - 1)
    tf.random.set_seed(seed)
    sionna_config.seed = seed
    speed_mps = speed_kmh_to_mps(float(channel.ue_speed_kmh))
    generator, tdl = _sionna_generator(
        str(channel.tdl_profile).upper(),
        float(channel.delay_spread_ns) * 1e-9,
        float(channel.carrier_frequency_hz),
        speed_mps,
        int(channel.num_sinusoids),
        int(n_rx),
        int(n_tx),
        int(grid.n_symbols),
        int(grid.n_fft),
        float(grid.scs_khz) * 1e3,
        int(grid.cyclic_prefix_length),
    )
    full = np.asarray(generator(batch_size=int(batch_size)).numpy(), dtype=np.complex128)
    expected_prefix = (int(batch_size), 1, int(n_rx), 1, int(n_tx), int(grid.n_symbols), int(grid.n_fft))
    if full.shape != expected_prefix:
        raise RuntimeError(f"Unexpected Sionna OFDM channel shape {full.shape}; expected {expected_prefix}.")
    full = full[:, 0, :, 0, :, :, :]
    active = full[..., grid.active_fft_indices]
    powers = np.asarray(tdl.mean_powers.numpy(), dtype=np.float64)
    delays_s = np.asarray(tdl.delays.numpy(), dtype=np.float64)
    return TDLRealization(
        taps=np.empty((0,), dtype=np.complex128),
        H=active,
        pdp=powers,
        tap_delays=delays_s,
        sample_period_ns=float("nan"),
        backend="sionna_tdl",
        metadata={
            "tdl_profile": str(channel.tdl_profile).upper(),
            "delay_spread_ns": float(channel.delay_spread_ns),
            "carrier_frequency_hz": float(channel.carrier_frequency_hz),
            "ue_speed_kmh": float(channel.ue_speed_kmh),
            "speed_mps": speed_mps,
            "num_sinusoids": int(channel.num_sinusoids),
            "normalization": "sionna_ensemble_no_realization_normalization",
            "full_fft_shape": list(full.shape),
            "active_fft_start": int(grid.active_fft_indices[0]),
            "active_fft_stop": int(grid.active_fft_indices[-1]),
        },
    )


def generate_tdl_channel(
    rng: np.random.Generator,
    grid: ResourceGrid,
    channel: ChannelConfig,
    n_tx: int,
    n_rx: int,
) -> TDLRealization:
    """Compatibility dispatcher for the legacy and Sionna TDL backends."""
    if str(getattr(channel, "backend", "legacy_exponential")).lower() == "sionna_tdl":
        seed = int(rng.integers(0, 2**31 - 1))
        return generate_sionna_tdl_channel(
            grid=grid,
            channel=channel,
            n_tx=n_tx,
            n_rx=n_rx,
            batch_size=1,
            seed=seed,
        )
    return _generate_legacy_exponential_channel(rng, grid, channel, n_tx, n_rx)
