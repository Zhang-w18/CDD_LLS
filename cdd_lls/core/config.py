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
    prg_size_rb: int = 4


@dataclass
class ChannelConfig:
    model: str = "tdl"
    delay_spread_ns: float = 30.0
    pdp: str = "exponential"
    max_delay_factor: float = 8.0
    normalize: bool = True


@dataclass
class TransmissionConfig:
    tx_scheme: str = "CDD"
    cdd_delay_vector: Optional[List[int]] = field(default_factory=lambda: [0, 8])
    cdd_base_delay: int = 8
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
    n_trials_per_snr: int = 20
    min_block_errors: int = 0
    max_trials_per_snr: int = 200
    seed: int = 42
    output_dir: str = "outputs"
    bler_target: float = 0.10
    save_trial_metrics: bool = False
    common_random_numbers: bool = True


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
    if int(config.antenna.n_tx) not in (1, 2, 4):
        raise ValueError(f"antenna.n_tx must be 1, 2, or 4. config={source_path}")
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
    if int(config.simulation.n_trials_per_snr) <= 0:
        raise ValueError(f"simulation.n_trials_per_snr must be positive. config={source_path}")
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
