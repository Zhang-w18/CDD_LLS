from __future__ import annotations

from dataclasses import asdict, dataclass, field, is_dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional
import copy

try:
    import yaml
except ModuleNotFoundError:  # pragma: no cover - only used by minimal runtime environments.
    yaml = None


@dataclass
class AntennaConfig:
    n_tx: int = 2
    n_rx: int = 4


@dataclass
class ResourceConfig:
    carrier_bandwidth_mhz: float = 100.0
    scs_khz: int = 30
    n_fft: int = 4096
    n_prbs: int = 8
    pdsch_n_symbols: int = 10
    dmrs_symbol_indices: List[int] = field(default_factory=lambda: [2, 7])
    dmrs_spacing_sc: int = 6
    dmrs_offset_sc: int = 0
    cyclic_prefix_length: int = 288
    prg_size_rb: int = 4


@dataclass
class ChannelConfig:
    backend: str = "legacy_exponential"
    model: str = "tdl"
    tdl_profile: str = "A"
    cdl_profile: str = "A"
    delay_spread_ns: float = 30.0
    carrier_frequency_hz: float = 3.5e9
    ue_speed_kmh: float = 0.0
    num_sinusoids: int = 20
    cdl_direction: str = "downlink"
    cdl_tx_array_rows: int = 1
    cdl_tx_array_cols: int = 0
    cdl_rx_array_rows: int = 1
    cdl_rx_array_cols: int = 0
    cdl_polarization: str = "single"
    cdl_polarization_type: str = "V"
    cdl_antenna_pattern: str = "omni"
    cdl_element_vertical_spacing: float = 0.5
    cdl_element_horizontal_spacing: float = 0.5
    pdp: str = "exponential"
    max_delay_factor: float = 8.0
    normalize: bool = True


@dataclass
class FixedCDLStatisticsConfig:
    covariance_realizations: int = 1000
    statistics_seed: int = 20261001
    realization_seed: int = 20262001
    velocity_azimuth_deg: float = 0.0
    velocity_elevation_deg: float = 0.0
    mean_aod_deg: float = 0.0
    aod_scale: float = 1.0
    target_aod_asd_deg: Optional[float] = None
    aoa_scale: float = 1.0
    zod_scale: float = 1.0
    zoa_scale: float = 1.0
    bs_vertical_aes: int = 8
    bs_horizontal_aes: int = 8
    bs_polarizations: int = 2
    bs_vertical_txrus_per_pol: int = 2
    bs_horizontal_txrus_per_pol: int = 8
    bs_vertical_spacing_lambda: float = 0.8
    bs_horizontal_spacing_lambda: float = 0.5
    bs_antenna_pattern: str = "38.901"
    bs_polarization_type: str = "cross"
    bs_orientation_deg: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    ue_vertical_elements: int = 1
    ue_horizontal_elements: int = 1
    ue_polarizations: int = 2
    ue_vertical_spacing_lambda: float = 0.5
    ue_horizontal_spacing_lambda: float = 0.5
    ue_antenna_pattern: str = "omni"
    ue_polarization_type: str = "cross"
    ue_orientation_deg: List[float] = field(default_factory=lambda: [180.0, 0.0, 0.0])
    ssb_horizontal_beams: int = 8
    secondary_horizontal_beams: int = 16
    pol_cycling_phase_deg: List[float] = field(default_factory=lambda: [0.0, 90.0])
    beam_cdd_delay_grid_indices: List[float] = field(default_factory=lambda: [0.0, 1.0])
    codebook_type: str = "legacy_parent_secondary"
    profile_native_angles: bool = False
    frozen_beam_manifest: str = ""
    frozen_beam_manifest_sha256: str = ""
    selection_only: bool = False


@dataclass
class TransmissionConfig:
    tx_scheme: str = "CDD"
    cdd_delay_vector: Optional[List[float]] = field(default_factory=lambda: [0, 8])
    cdd_base_delay: float = 8
    prg_codebook: str = "qpsk_dft"
    prg_cycling_order: Optional[List[int]] = None


@dataclass
class ChannelEstimationConfig:
    ce_method: str = "RMMSE_WB_KNOWN"
    rmmse_bundle_rb: int = 4
    diagonal_loading: float = 1e-8
    recon_pair_spacing_pilots: int = 1
    recon_regularization: float = 1e-3
    basis_support: str = "truncated"
    basis_energy_threshold: float = 0.99
    cond_warning_threshold: float = 1e3


@dataclass
class ReceiverConfig:
    equalizer: str = "zf_mrc"
    max_ldpc_iterations: int = 20
    llr_clip: float = 50.0


@dataclass
class MCSConfig:
    table: str = "nr_256qam"
    index: int = 8
    qm: Optional[int] = None
    code_rate: Optional[float] = None


@dataclass
class SimulationConfig:
    snr_range_db: List[float] = field(default_factory=lambda: [-4, 12, 2])
    snr_points_db: List[float] = field(default_factory=list)
    n_trials_per_snr: int = 20
    min_block_errors: int = 0
    max_trials_per_snr: int = 200
    seed: int = 42
    output_dir: str = "outputs"
    run_id: Optional[str] = None
    bler_target: float = 0.10
    save_trial_metrics: bool = False
    common_random_numbers: bool = True
    absolute_trial_start: int = 0
    ce_only: bool = False


@dataclass
class SweepConfig:
    cdd_base_delays: List[int] = field(default_factory=list)
    dmrs_spacing_sc: List[int] = field(default_factory=list)


@dataclass
class PlotConfig:
    enabled: bool = True


@dataclass
class PlatformConfig:
    antenna: AntennaConfig = field(default_factory=AntennaConfig)
    resource: ResourceConfig = field(default_factory=ResourceConfig)
    channel: ChannelConfig = field(default_factory=ChannelConfig)
    fixed_cdl_statistics: FixedCDLStatisticsConfig = field(default_factory=FixedCDLStatisticsConfig)
    transmission: TransmissionConfig = field(default_factory=TransmissionConfig)
    channel_estimation: ChannelEstimationConfig = field(default_factory=ChannelEstimationConfig)
    receiver: ReceiverConfig = field(default_factory=ReceiverConfig)
    mcs: MCSConfig = field(default_factory=MCSConfig)
    simulation: SimulationConfig = field(default_factory=SimulationConfig)
    sweeps: SweepConfig = field(default_factory=SweepConfig)
    plots: PlotConfig = field(default_factory=PlotConfig)
    variants: List[Dict[str, Any]] = field(default_factory=list)
    scenarios: List[Dict[str, Any]] = field(default_factory=list)


def _deep_update(base: Dict[str, Any], patch: Dict[str, Any]) -> Dict[str, Any]:
    out = copy.deepcopy(base)
    for key, value in (patch or {}).items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_update(out[key], value)
        else:
            out[key] = value
    return out


def dataclass_to_dict(obj: Any) -> Any:
    if is_dataclass(obj):
        return {k: dataclass_to_dict(v) for k, v in asdict(obj).items()}
    if isinstance(obj, list):
        return [dataclass_to_dict(v) for v in obj]
    if isinstance(obj, dict):
        return {k: dataclass_to_dict(v) for k, v in obj.items()}
    return obj


def _construct_dataclass(cls, data: Dict[str, Any]):
    child_map = {
        "antenna": AntennaConfig,
        "resource": ResourceConfig,
        "channel": ChannelConfig,
        "fixed_cdl_statistics": FixedCDLStatisticsConfig,
        "transmission": TransmissionConfig,
        "channel_estimation": ChannelEstimationConfig,
        "receiver": ReceiverConfig,
        "mcs": MCSConfig,
        "simulation": SimulationConfig,
        "sweeps": SweepConfig,
        "plots": PlotConfig,
    }
    kwargs = {}
    for field_name in cls.__dataclass_fields__:
        if field_name not in data:
            continue
        value = data[field_name]
        if field_name in child_map and isinstance(value, dict):
            kwargs[field_name] = _construct_dataclass(child_map[field_name], value)
        else:
            kwargs[field_name] = value
    return cls(**kwargs)


def _validate_config(config: PlatformConfig, source_path: str = "") -> None:
    if int(config.antenna.n_tx) not in (1, 2, 4, 8, 16, 32):
        raise ValueError(f"antenna.n_tx must be 1, 2, 4, 8, 16, or 32. config={source_path}")
    if int(config.antenna.n_rx) <= 0:
        raise ValueError(f"antenna.n_rx must be positive. config={source_path}")
    if int(config.resource.n_prbs) <= 0:
        raise ValueError(f"resource.n_prbs must be positive. config={source_path}")
    if int(config.resource.dmrs_spacing_sc) <= 0:
        raise ValueError(f"resource.dmrs_spacing_sc must be positive. config={source_path}")
    if not config.resource.dmrs_symbol_indices:
        raise ValueError(f"resource.dmrs_symbol_indices must not be empty. config={source_path}")
    if int(config.resource.n_fft) <= int(config.resource.n_prbs) * 12:
        raise ValueError(f"resource.n_fft must exceed active subcarriers. config={source_path}")
    if int(config.resource.cyclic_prefix_length) < 0:
        raise ValueError(f"resource.cyclic_prefix_length must be non-negative. config={source_path}")
    if any(int(s) < 0 or int(s) >= int(config.resource.pdsch_n_symbols)
           for s in config.resource.dmrs_symbol_indices):
        raise ValueError(f"resource.dmrs_symbol_indices are outside the PDSCH grid. config={source_path}")
    backend = str(config.channel.backend).lower()
    if backend not in ("legacy_exponential", "sionna_tdl", "sionna_cdl", "fixed_cdl_statistics"):
        raise ValueError(
            "channel.backend must be legacy_exponential, sionna_tdl, sionna_cdl, or "
            "fixed_cdl_statistics. "
            f"config={source_path}"
        )
    if str(config.channel.tdl_profile).upper() not in ("A", "B", "C", "D", "E", "A30", "B100", "C300"):
        raise ValueError(f"channel.tdl_profile is invalid. config={source_path}")
    if str(config.channel.cdl_profile).upper() not in ("A", "B", "C", "D", "E"):
        raise ValueError(f"channel.cdl_profile must be A, B, C, D, or E. config={source_path}")
    if float(config.channel.delay_spread_ns) <= 0:
        raise ValueError(f"channel.delay_spread_ns must be positive. config={source_path}")
    if float(config.channel.carrier_frequency_hz) <= 0:
        raise ValueError(f"channel.carrier_frequency_hz must be positive. config={source_path}")
    if float(config.channel.ue_speed_kmh) < 0:
        raise ValueError(f"channel.ue_speed_kmh must be non-negative. config={source_path}")
    if int(config.channel.num_sinusoids) <= 0:
        raise ValueError(f"channel.num_sinusoids must be positive. config={source_path}")
    if str(config.channel.cdl_direction).lower() not in ("downlink", "uplink"):
        raise ValueError(f"channel.cdl_direction must be downlink or uplink. config={source_path}")
    if str(config.channel.cdl_polarization).lower() != "single":
        raise ValueError(
            "The platform currently supports channel.cdl_polarization=single only. "
            f"config={source_path}"
        )
    if str(config.channel.cdl_polarization_type).upper() not in ("V", "H"):
        raise ValueError(
            "Single-polarized CDL arrays require cdl_polarization_type V or H. "
            f"config={source_path}"
        )
    if str(config.channel.cdl_antenna_pattern).lower() not in ("omni", "38.901"):
        raise ValueError(
            "channel.cdl_antenna_pattern must be omni or 38.901. "
            f"config={source_path}"
        )
    for name in (
        "cdl_tx_array_rows",
        "cdl_tx_array_cols",
        "cdl_rx_array_rows",
        "cdl_rx_array_cols",
    ):
        if int(getattr(config.channel, name)) < 0:
            raise ValueError(f"channel.{name} must be non-negative. config={source_path}")
    for name in (
        "cdl_element_vertical_spacing",
        "cdl_element_horizontal_spacing",
    ):
        if float(getattr(config.channel, name)) <= 0:
            raise ValueError(f"channel.{name} must be positive. config={source_path}")
    if backend == "sionna_cdl":
        tx_rows = int(config.channel.cdl_tx_array_rows)
        tx_cols = int(config.channel.cdl_tx_array_cols) or int(config.antenna.n_tx)
        rx_rows = int(config.channel.cdl_rx_array_rows)
        rx_cols = int(config.channel.cdl_rx_array_cols) or int(config.antenna.n_rx)
        if tx_rows * tx_cols != int(config.antenna.n_tx):
            raise ValueError(
                "CDL Tx array rows*cols must equal antenna.n_tx; zero cols selects an automatic ULA. "
                f"config={source_path}"
            )
        if rx_rows * rx_cols != int(config.antenna.n_rx):
            raise ValueError(
                "CDL Rx array rows*cols must equal antenna.n_rx; zero cols selects an automatic ULA. "
                f"config={source_path}"
            )
    if backend == "fixed_cdl_statistics":
        fixed = config.fixed_cdl_statistics
        if str(config.channel.cdl_direction).lower() != "downlink":
            raise ValueError(f"fixed_cdl_statistics currently requires downlink. config={source_path}")
        if int(fixed.covariance_realizations) <= 0:
            raise ValueError(f"fixed CDL covariance_realizations must be positive. config={source_path}")
        if int(fixed.statistics_seed) == int(fixed.realization_seed):
            raise ValueError(f"fixed CDL statistics and realization seeds must differ. config={source_path}")
        if int(fixed.bs_polarizations) != 2:
            raise ValueError(f"fixed CDL BS array must use two polarizations. config={source_path}")
        if int(fixed.ue_polarizations) not in (1, 2):
            raise ValueError(f"fixed CDL UE polarizations must be 1 or 2. config={source_path}")
        if int(fixed.bs_vertical_aes) % int(fixed.bs_vertical_txrus_per_pol) or int(
            fixed.bs_horizontal_aes
        ) % int(fixed.bs_horizontal_txrus_per_pol):
            raise ValueError(f"fixed CDL BS AE dimensions must be divisible by TXRU dimensions. config={source_path}")
        expected_tx = (
            int(fixed.bs_vertical_txrus_per_pol)
            * int(fixed.bs_horizontal_txrus_per_pol)
            * int(fixed.bs_polarizations)
        )
        expected_rx = (
            int(fixed.ue_vertical_elements)
            * int(fixed.ue_horizontal_elements)
            * int(fixed.ue_polarizations)
        )
        if expected_tx != int(config.antenna.n_tx):
            raise ValueError(
                f"fixed CDL BS TXRU count {expected_tx} must equal antenna.n_tx={config.antenna.n_tx}. "
                f"config={source_path}"
            )
        if expected_rx != int(config.antenna.n_rx):
            raise ValueError(
                f"fixed CDL UE port count {expected_rx} must equal antenna.n_rx={config.antenna.n_rx}. "
                f"config={source_path}"
            )
        codebook_type = str(fixed.codebook_type).lower()
        if codebook_type not in {
            "legacy_parent_secondary", "dft_2x8_same_pol", "wide_beam_split",
            "angular_full_coverage_ultrawide",
        }:
            raise ValueError(f"Unsupported fixed CDL codebook_type={fixed.codebook_type!r}. config={source_path}")
        if codebook_type == "legacy_parent_secondary":
            if int(fixed.secondary_horizontal_beams) != 2 * int(fixed.ssb_horizontal_beams):
                raise ValueError(f"fixed CDL requires two secondary beams per SSB. config={source_path}")
            if len(fixed.pol_cycling_phase_deg) < 1 or len(fixed.beam_cdd_delay_grid_indices) != 2:
                raise ValueError(f"fixed CDL polarization/CDD patterns are invalid. config={source_path}")
            supported = {"BASELINE", "POLARIZATION_CYCLING", "BEAM_CYCLING", "BEAM_CDD"}
        else:
            supported = (
                {"BEAM8_B0_QC", "BEAM8_SIDON_SELECTED", "BEAM8_PRECODER_CYCLING"}
                if codebook_type == "angular_full_coverage_ultrawide"
                else {"BEAM8_B0_QC", "BEAM8_S0_SIDON", "BEAM8_PRECODER_CYCLING"}
            )
            if codebook_type == "angular_full_coverage_ultrawide":
                if bool(fixed.profile_native_angles):
                    raise ValueError(f"Plan-039 requires transformed AoD angles. config={source_path}")
                if fixed.target_aod_asd_deg is None or float(fixed.target_aod_asd_deg) <= 0.0:
                    raise ValueError(f"Plan-039 requires positive target_aod_asd_deg. config={source_path}")
                if any(float(value) != 1.0 for value in (fixed.aoa_scale, fixed.zod_scale, fixed.zoa_scale)):
                    raise ValueError(f"Plan-039 only permits AoD transformation. config={source_path}")
            else:
                if not bool(fixed.profile_native_angles):
                    raise ValueError(f"dft_2x8_same_pol requires profile_native_angles=true. config={source_path}")
                if any(float(value) != expected for value, expected in zip(
                    (fixed.mean_aod_deg, fixed.aod_scale, fixed.aoa_scale, fixed.zod_scale, fixed.zoa_scale),
                    (0.0, 1.0, 1.0, 1.0, 1.0),
                )):
                    raise ValueError(f"dft_2x8_same_pol forbids angle shifts/scales. config={source_path}")
            if (not bool(fixed.selection_only)) and (
                not fixed.frozen_beam_manifest or len(str(fixed.frozen_beam_manifest_sha256)) != 64
            ):
                raise ValueError(f"dft_2x8_same_pol requires a frozen manifest path and SHA-256. config={source_path}")
            if int(config.resource.n_prbs) != 48 or int(config.resource.n_fft) != 4096:
                raise ValueError(f"dft_2x8_same_pol requires 48 PRB and FFT 4096. config={source_path}")
        if str(config.transmission.tx_scheme).upper() not in supported:
            raise ValueError(
                f"fixed CDL transmission.tx_scheme must be one of {sorted(supported)}. config={source_path}"
            )
        ce_method = str(config.channel_estimation.ce_method).upper()
        if codebook_type in {"wide_beam_split", "angular_full_coverage_ultrawide"}:
            from cdd_lls.phy.plan039 import CE_METHODS
            allowed_ce = {"IDEAL", *CE_METHODS}
            if ce_method not in allowed_ce:
                raise ValueError(f"Plan-039 CE method must be one of {sorted(allowed_ce)}. config={source_path}")
            scheme = str(config.transmission.tx_scheme).upper()
            if scheme == "BEAM8_PRECODER_CYCLING" and ce_method not in {"IDEAL", "PLAN039_PRG_COMMON_REFERENCE_PDP"}:
                raise ValueError("Plan-039 cycling requires IDEAL or PRG common-reference PDP estimation.")
            if scheme != "BEAM8_PRECODER_CYCLING" and ce_method == "PLAN039_PRG_COMMON_REFERENCE_PDP":
                raise ValueError("Plan-039 PRG estimator is only valid for precoder cycling.")
            if scheme == "BEAM8_PRECODER_CYCLING" and ce_method == "PLAN039_TRANSPARENT_COMMON_REFERENCE_PDP":
                raise ValueError("Plan-039 full-band transparent estimator is only valid for CDD schemes.")
        elif codebook_type == "dft_2x8_same_pol":
            allowed_ce = {"IDEAL", "BEAM8_PRG_LMMSE", "BEAM8_CDD_AWARE_LMMSE"}
            if ce_method not in allowed_ce:
                raise ValueError(f"Beam8 fixed CDL CE method must be one of {sorted(allowed_ce)}. config={source_path}")
            scheme = str(config.transmission.tx_scheme).upper()
            if scheme == "BEAM8_PRECODER_CYCLING" and ce_method not in {"IDEAL", "BEAM8_PRG_LMMSE"}:
                raise ValueError(f"Beam8 cycling requires IDEAL or BEAM8_PRG_LMMSE. config={source_path}")
            if scheme in {"BEAM8_B0_QC", "BEAM8_S0_SIDON"} and ce_method == "BEAM8_PRG_LMMSE":
                raise ValueError(f"Beam8 CDD requires IDEAL or BEAM8_CDD_AWARE_LMMSE. config={source_path}")
        elif ce_method != "IDEAL":
            raise ValueError("Legacy fixed_cdl_statistics currently supports IDEAL only. " f"config={source_path}")
        if not bool(config.simulation.common_random_numbers):
            raise ValueError(
                "fixed_cdl_statistics requires simulation.common_random_numbers=true. "
                f"config={source_path}"
            )
    if int(config.simulation.n_trials_per_snr) <= 0:
        raise ValueError(f"simulation.n_trials_per_snr must be positive. config={source_path}")
    if int(config.simulation.absolute_trial_start) < 0:
        raise ValueError(f"simulation.absolute_trial_start must be nonnegative. config={source_path}")
    if len(config.simulation.snr_range_db) != 3:
        raise ValueError(f"simulation.snr_range_db must be [start, stop, step]. config={source_path}")


def load_config(path: str | Path) -> PlatformConfig:
    if yaml is None:
        raise ModuleNotFoundError("PyYAML is required to load YAML config files.")
    default = dataclass_to_dict(PlatformConfig())
    with open(path, "r", encoding="utf-8") as f:
        user = yaml.safe_load(f) or {}
    merged = _deep_update(default, user)
    config = _construct_dataclass(PlatformConfig, merged)
    _validate_config(config, source_path=str(path))
    return config


def save_resolved_config(config: PlatformConfig, path: str | Path) -> None:
    if yaml is None:
        raise ModuleNotFoundError("PyYAML is required to save YAML config files.")
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(dataclass_to_dict(config), f, allow_unicode=True, sort_keys=False)


def merged_config_dict(config: PlatformConfig, patch: Dict[str, Any]) -> Dict[str, Any]:
    return _deep_update(dataclass_to_dict(config), patch or {})


def config_from_dict(data: Dict[str, Any]) -> PlatformConfig:
    cfg = _construct_dataclass(PlatformConfig, data)
    _validate_config(cfg)
    return cfg
