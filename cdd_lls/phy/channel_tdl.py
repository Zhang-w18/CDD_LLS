from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math
from typing import Dict, Sequence
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
        metadata={
            "normalization": "unit_ensemble_pdp" if bool(channel.normalize) else "legacy_unnormalized",
            "rx_correlation_model": "identity_no_spatial_correlation",
        },
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


def _resolved_cdl_array_shape(channel: ChannelConfig, side: str, n_ant: int) -> tuple[int, int]:
    rows = int(getattr(channel, f"cdl_{side}_array_rows"))
    configured_cols = int(getattr(channel, f"cdl_{side}_array_cols"))
    cols = configured_cols or int(n_ant)
    if rows <= 0 or cols <= 0 or rows * cols != int(n_ant):
        raise ValueError(
            f"CDL {side} array rows*cols must equal {side.replace('tx', 'n_tx').replace('rx', 'n_rx')}="
            f"{int(n_ant)}; got {rows}x{cols}."
        )
    return rows, cols


@lru_cache(maxsize=32)
def _sionna_cdl_generator(
    profile: str,
    delay_spread_s: float,
    carrier_frequency_hz: float,
    speed_mps: float,
    direction: str,
    tx_rows: int,
    tx_cols: int,
    rx_rows: int,
    rx_cols: int,
    polarization: str,
    polarization_type: str,
    antenna_pattern: str,
    element_vertical_spacing: float,
    element_horizontal_spacing: float,
    n_symbols: int,
    n_fft: int,
    subcarrier_spacing_hz: float,
    cyclic_prefix_length: int,
):
    from sionna.phy.channel import GenerateOFDMChannel
    from sionna.phy.channel.tr38901 import CDL, PanelArray
    from sionna.phy.ofdm import ResourceGrid as SionnaResourceGrid

    def build_array(rows: int, cols: int) -> PanelArray:
        return PanelArray(
            num_rows_per_panel=int(rows),
            num_cols_per_panel=int(cols),
            polarization=str(polarization).lower(),
            polarization_type=str(polarization_type).upper(),
            antenna_pattern=str(antenna_pattern).lower(),
            carrier_frequency=float(carrier_frequency_hz),
            element_vertical_spacing=float(element_vertical_spacing),
            element_horizontal_spacing=float(element_horizontal_spacing),
            precision="double",
        )

    tx_array = build_array(tx_rows, tx_cols)
    rx_array = build_array(rx_rows, rx_cols)
    if str(direction).lower() == "downlink":
        ut_array, bs_array = rx_array, tx_array
    else:
        ut_array, bs_array = tx_array, rx_array
    cdl = CDL(
        model=str(profile).upper(),
        delay_spread=float(delay_spread_s),
        carrier_frequency=float(carrier_frequency_hz),
        ut_array=ut_array,
        bs_array=bs_array,
        direction=str(direction).lower(),
        min_speed=float(speed_mps),
        max_speed=float(speed_mps),
        precision="double",
    )
    sionna_grid = SionnaResourceGrid(
        num_ofdm_symbols=int(n_symbols),
        fft_size=int(n_fft),
        subcarrier_spacing=float(subcarrier_spacing_hz),
        num_tx=1,
        num_streams_per_tx=1,
        cyclic_prefix_length=int(cyclic_prefix_length),
        precision="double",
    )
    generator = GenerateOFDMChannel(
        cdl,
        sionna_grid,
        normalize_channel=False,
        precision="double",
    )
    return generator, cdl


def _cdl_generator_for_config(
    grid: ResourceGrid,
    channel: ChannelConfig,
    n_tx: int,
    n_rx: int,
    n_symbols: int,
):
    tx_rows, tx_cols = _resolved_cdl_array_shape(channel, "tx", n_tx)
    rx_rows, rx_cols = _resolved_cdl_array_shape(channel, "rx", n_rx)
    return _sionna_cdl_generator(
        str(channel.cdl_profile).upper(),
        float(channel.delay_spread_ns) * 1e-9,
        float(channel.carrier_frequency_hz),
        speed_kmh_to_mps(float(channel.ue_speed_kmh)),
        str(channel.cdl_direction).lower(),
        tx_rows,
        tx_cols,
        rx_rows,
        rx_cols,
        str(channel.cdl_polarization).lower(),
        str(channel.cdl_polarization_type).upper(),
        str(channel.cdl_antenna_pattern).lower(),
        float(channel.cdl_element_vertical_spacing),
        float(channel.cdl_element_horizontal_spacing),
        int(n_symbols),
        int(grid.n_fft),
        float(grid.scs_khz) * 1e3,
        int(grid.cyclic_prefix_length),
    )


def _cdl_metadata(channel: ChannelConfig, n_tx: int, n_rx: int) -> Dict[str, object]:
    tx_rows, tx_cols = _resolved_cdl_array_shape(channel, "tx", n_tx)
    rx_rows, rx_cols = _resolved_cdl_array_shape(channel, "rx", n_rx)
    return {
        "cdl_profile": str(channel.cdl_profile).upper(),
        "delay_spread_ns": float(channel.delay_spread_ns),
        "carrier_frequency_hz": float(channel.carrier_frequency_hz),
        "ue_speed_kmh": float(channel.ue_speed_kmh),
        "speed_mps": speed_kmh_to_mps(float(channel.ue_speed_kmh)),
        "direction": str(channel.cdl_direction).lower(),
        "tx_array_shape": [tx_rows, tx_cols],
        "rx_array_shape": [rx_rows, rx_cols],
        "polarization": str(channel.cdl_polarization).lower(),
        "polarization_type": str(channel.cdl_polarization_type).upper(),
        "antenna_pattern": str(channel.cdl_antenna_pattern).lower(),
        "element_vertical_spacing_wavelengths": float(channel.cdl_element_vertical_spacing),
        "element_horizontal_spacing_wavelengths": float(channel.cdl_element_horizontal_spacing),
        "normalization": "sionna_cdl_ensemble_no_realization_normalization",
        "spatial_correlation_model": "sionna_cdl_array_geometry",
    }


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
            "rx_correlation_model": "identity_no_spatial_correlation",
            "full_fft_shape": list(full.shape),
            "active_fft_start": int(grid.active_fft_indices[0]),
            "active_fft_stop": int(grid.active_fft_indices[-1]),
        },
    )


def generate_sionna_tdl_channel_active(
    grid: ResourceGrid,
    channel: ChannelConfig,
    n_tx: int,
    n_rx: int,
    batch_size: int = 1,
    seed: int = 0,
    time_sample_indices: Sequence[int] | None = None,
) -> TDLRealization:
    """Generate the same Sionna TDL CIR directly on active subcarriers.

    This avoids constructing an intermediate ``[..., n_fft]`` response when
    the simulation only consumes the active CORESET band.
    """
    import tensorflow as tf
    from sionna.phy import config as sionna_config
    from sionna.phy.channel.utils import cir_to_ofdm_channel

    seed = int(seed) % (2**31 - 1)
    tf.random.set_seed(seed)
    sionna_config.seed = seed
    speed_mps = speed_kmh_to_mps(float(channel.ue_speed_kmh))
    if time_sample_indices is None:
        selected_times = np.arange(int(grid.n_symbols), dtype=np.int64)
    else:
        selected_times = np.asarray(time_sample_indices, dtype=np.int64).reshape(-1)
        if selected_times.size == 0:
            raise ValueError("time_sample_indices must not be empty")
        if np.any(selected_times < 0):
            raise ValueError("time_sample_indices must be non-negative")
        if len(np.unique(selected_times)) != len(selected_times):
            raise ValueError("time_sample_indices must not contain duplicates")
    generated_time_steps = int(np.max(selected_times)) + 1
    generator, tdl = _sionna_generator(
        str(channel.tdl_profile).upper(),
        float(channel.delay_spread_ns) * 1e-9,
        float(channel.carrier_frequency_hz),
        speed_mps,
        int(channel.num_sinusoids),
        int(n_rx),
        int(n_tx),
        generated_time_steps,
        int(grid.n_fft),
        float(grid.scs_khz) * 1e3,
        int(grid.cyclic_prefix_length),
    )
    h, tau = tdl(
        int(batch_size),
        generated_time_steps,
        1.0 / float(grid.ofdm_symbol_duration_s),
    )
    h = tf.gather(h, tf.convert_to_tensor(selected_times, dtype=tf.int32), axis=-1)
    frequencies = tf.convert_to_tensor(
        np.asarray(grid.subcarrier_indices, dtype=np.float64)
        * float(grid.scs_khz)
        * 1e3,
        dtype=tf.float64,
    )
    active = np.asarray(
        cir_to_ofdm_channel(frequencies, h, tau, normalize=False).numpy(),
        dtype=np.complex128,
    )
    expected = (
        int(batch_size),
        1,
        int(n_rx),
        1,
        int(n_tx),
        int(len(selected_times)),
        int(grid.n_sc),
    )
    if active.shape != expected:
        raise RuntimeError(
            f"Unexpected active-only Sionna channel shape {active.shape}; expected {expected}."
        )
    active = active[:, 0, :, 0, :, :, :]
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
            "rx_correlation_model": "identity_no_spatial_correlation",
            "frequency_evaluation": "direct_active_subcarriers",
            "generated_time_steps": generated_time_steps,
            "selected_time_sample_indices": selected_times.astype(int).tolist(),
            "selected_time_seconds": (
                selected_times.astype(np.float64) * float(grid.ofdm_symbol_duration_s)
            ).tolist(),
            "active_shape": list(active.shape),
            "active_fft_start": int(grid.active_fft_indices[0]),
            "active_fft_stop": int(grid.active_fft_indices[-1]),
        },
    )


def generate_sionna_cdl_channel(
    grid: ResourceGrid,
    channel: ChannelConfig,
    n_tx: int,
    n_rx: int,
    batch_size: int = 1,
    seed: int = 0,
) -> TDLRealization:
    """Generate Sionna CDL responses as [batch,rx,tx,symbol,active_sc]."""
    import tensorflow as tf
    from sionna.phy import config as sionna_config

    seed = int(seed) % (2**31 - 1)
    tf.random.set_seed(seed)
    sionna_config.seed = seed
    generator, cdl = _cdl_generator_for_config(
        grid, channel, int(n_tx), int(n_rx), int(grid.n_symbols)
    )
    full = np.asarray(generator(batch_size=int(batch_size)).numpy(), dtype=np.complex128)
    expected = (
        int(batch_size),
        1,
        int(n_rx),
        1,
        int(n_tx),
        int(grid.n_symbols),
        int(grid.n_fft),
    )
    if full.shape != expected:
        raise RuntimeError(f"Unexpected Sionna CDL channel shape {full.shape}; expected {expected}.")
    full = full[:, 0, :, 0, :, :, :]
    active = full[..., grid.active_fft_indices]
    metadata = _cdl_metadata(channel, n_tx, n_rx)
    metadata.update(
        {
            "full_fft_shape": list(full.shape),
            "active_fft_start": int(grid.active_fft_indices[0]),
            "active_fft_stop": int(grid.active_fft_indices[-1]),
        }
    )
    return TDLRealization(
        taps=np.empty((0,), dtype=np.complex128),
        H=active,
        pdp=np.asarray(cdl.powers.numpy(), dtype=np.float64).reshape(-1),
        tap_delays=np.asarray(cdl.delays.numpy(), dtype=np.float64).reshape(-1),
        sample_period_ns=float("nan"),
        backend="sionna_cdl",
        metadata=metadata,
    )


def generate_sionna_cdl_channel_active(
    grid: ResourceGrid,
    channel: ChannelConfig,
    n_tx: int,
    n_rx: int,
    batch_size: int = 1,
    seed: int = 0,
    time_sample_indices: Sequence[int] | None = None,
) -> TDLRealization:
    """Generate a Sionna CDL CIR directly on the active subcarriers."""
    import tensorflow as tf
    from sionna.phy import config as sionna_config
    from sionna.phy.channel.utils import cir_to_ofdm_channel

    seed = int(seed) % (2**31 - 1)
    tf.random.set_seed(seed)
    sionna_config.seed = seed
    if time_sample_indices is None:
        selected_times = np.arange(int(grid.n_symbols), dtype=np.int64)
    else:
        selected_times = np.asarray(time_sample_indices, dtype=np.int64).reshape(-1)
        if selected_times.size == 0:
            raise ValueError("time_sample_indices must not be empty")
        if np.any(selected_times < 0):
            raise ValueError("time_sample_indices must be non-negative")
        if len(np.unique(selected_times)) != len(selected_times):
            raise ValueError("time_sample_indices must not contain duplicates")
    generated_time_steps = int(np.max(selected_times)) + 1
    _, cdl = _cdl_generator_for_config(
        grid, channel, int(n_tx), int(n_rx), generated_time_steps
    )
    h, tau = cdl(
        int(batch_size),
        generated_time_steps,
        1.0 / float(grid.ofdm_symbol_duration_s),
    )
    h = tf.gather(h, tf.convert_to_tensor(selected_times, dtype=tf.int32), axis=-1)
    frequencies = tf.convert_to_tensor(
        np.asarray(grid.subcarrier_indices, dtype=np.float64) * float(grid.scs_khz) * 1e3,
        dtype=tf.float64,
    )
    active = np.asarray(
        cir_to_ofdm_channel(frequencies, h, tau, normalize=False).numpy(),
        dtype=np.complex128,
    )
    expected = (
        int(batch_size),
        1,
        int(n_rx),
        1,
        int(n_tx),
        int(len(selected_times)),
        int(grid.n_sc),
    )
    if active.shape != expected:
        raise RuntimeError(
            f"Unexpected active-only Sionna CDL channel shape {active.shape}; expected {expected}."
        )
    active = active[:, 0, :, 0, :, :, :]
    metadata = _cdl_metadata(channel, n_tx, n_rx)
    metadata.update(
        {
            "frequency_evaluation": "direct_active_subcarriers",
            "generated_time_steps": generated_time_steps,
            "selected_time_sample_indices": selected_times.astype(int).tolist(),
            "selected_time_seconds": (
                selected_times.astype(np.float64) * float(grid.ofdm_symbol_duration_s)
            ).tolist(),
            "active_shape": list(active.shape),
            "active_fft_start": int(grid.active_fft_indices[0]),
            "active_fft_stop": int(grid.active_fft_indices[-1]),
        }
    )
    return TDLRealization(
        taps=np.empty((0,), dtype=np.complex128),
        H=active,
        pdp=np.asarray(cdl.powers.numpy(), dtype=np.float64).reshape(-1),
        tap_delays=np.asarray(cdl.delays.numpy(), dtype=np.float64).reshape(-1),
        sample_period_ns=float("nan"),
        backend="sionna_cdl",
        metadata=metadata,
    )


def generate_sionna_channel(
    grid: ResourceGrid,
    channel: ChannelConfig,
    n_tx: int,
    n_rx: int,
    batch_size: int = 1,
    seed: int = 0,
) -> TDLRealization:
    """Dispatch to the configured Sionna TDL or CDL backend."""
    backend = str(channel.backend).lower()
    if backend == "sionna_tdl":
        return generate_sionna_tdl_channel(grid, channel, n_tx, n_rx, batch_size, seed)
    if backend == "sionna_cdl":
        return generate_sionna_cdl_channel(grid, channel, n_tx, n_rx, batch_size, seed)
    raise ValueError(f"Sionna channel generation requires sionna_tdl or sionna_cdl, got {backend}.")


def generate_sionna_channel_active(
    grid: ResourceGrid,
    channel: ChannelConfig,
    n_tx: int,
    n_rx: int,
    batch_size: int = 1,
    seed: int = 0,
    time_sample_indices: Sequence[int] | None = None,
) -> TDLRealization:
    """Dispatch active-subcarrier generation to Sionna TDL or CDL."""
    backend = str(channel.backend).lower()
    if backend == "sionna_tdl":
        return generate_sionna_tdl_channel_active(
            grid, channel, n_tx, n_rx, batch_size, seed, time_sample_indices
        )
    if backend == "sionna_cdl":
        return generate_sionna_cdl_channel_active(
            grid, channel, n_tx, n_rx, batch_size, seed, time_sample_indices
        )
    raise ValueError(f"Sionna channel generation requires sionna_tdl or sionna_cdl, got {backend}.")


def generate_tdl_channel(
    rng: np.random.Generator,
    grid: ResourceGrid,
    channel: ChannelConfig,
    n_tx: int,
    n_rx: int,
) -> TDLRealization:
    """Compatibility dispatcher for the legacy and Sionna TDL backends."""
    if str(getattr(channel, "backend", "legacy_exponential")).lower() in (
        "sionna_tdl",
        "sionna_cdl",
    ):
        seed = int(rng.integers(0, 2**31 - 1))
        return generate_sionna_channel(
            grid=grid,
            channel=channel,
            n_tx=n_tx,
            n_rx=n_rx,
            batch_size=1,
            seed=seed,
        )
    return _generate_legacy_exponential_channel(rng, grid, channel, n_tx, n_rx)
