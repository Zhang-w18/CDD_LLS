from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from cdd_lls.core.config import ResourceConfig


@dataclass(frozen=True)
class ResourceGrid:
    n_sc: int
    n_symbols: int
    n_fft: int
    scs_khz: int
    subcarrier_indices: np.ndarray
    active_fft_indices: np.ndarray
    pilot_subcarriers: np.ndarray
    pilot_symbol_indices: np.ndarray
    pilot_subcarrier_indices: np.ndarray
    data_symbol_indices: np.ndarray
    data_subcarrier_indices: np.ndarray
    cyclic_prefix_length: int
    ofdm_symbol_duration_s: float
    dmrs_overhead: float
    n_dmrs_re: int
    n_data_re: int

    @property
    def pilot_count(self) -> int:
        """Number of pilot subcarriers in one DMRS symbol (legacy meaning)."""
        return int(len(self.pilot_subcarriers))

    @property
    def pilot_re_count(self) -> int:
        return int(len(self.pilot_symbol_indices))

    @property
    def pilot_coordinates(self) -> np.ndarray:
        return np.column_stack((self.pilot_symbol_indices, self.pilot_subcarrier_indices))

    @property
    def data_coordinates(self) -> np.ndarray:
        return np.column_stack((self.data_symbol_indices, self.data_subcarrier_indices))


def build_resource_grid(resource: ResourceConfig) -> ResourceGrid:
    n_sc = int(resource.n_prbs) * 12
    n_symbols = int(resource.pdsch_n_symbols)
    n_fft = int(resource.n_fft)
    half = n_sc // 2
    if n_sc % 2 == 0:
        subcarrier_indices = np.arange(-half, half, dtype=np.int64)
    else:
        subcarrier_indices = np.arange(-half, half + 1, dtype=np.int64)

    pilot_local = np.arange(
        int(resource.dmrs_offset_sc),
        n_sc,
        int(resource.dmrs_spacing_sc),
        dtype=np.int64,
    )
    pilot_local = pilot_local[(pilot_local >= 0) & (pilot_local < n_sc)]
    if pilot_local.size == 0:
        raise ValueError("DMRS pattern produced no pilot subcarriers.")
    pilot_subcarriers = subcarrier_indices[pilot_local]

    dmrs_symbol_list = sorted(set(int(x) for x in resource.dmrs_symbol_indices))
    dmrs_symbols = set(dmrs_symbol_list)
    pilot_symbols = np.repeat(np.asarray(dmrs_symbol_list, dtype=np.int64), len(pilot_subcarriers))
    pilot_subcarrier_indices = np.tile(pilot_subcarriers, len(dmrs_symbol_list))
    data_symbols = []
    data_subcarriers = []
    for sym in range(n_symbols):
        if sym in dmrs_symbols:
            mask = np.ones(n_sc, dtype=bool)
            mask[pilot_local] = False
            ks = subcarrier_indices[mask]
        else:
            ks = subcarrier_indices
        data_symbols.extend([sym] * len(ks))
        data_subcarriers.extend(ks.tolist())

    n_dmrs_re = int(len(pilot_symbols))
    n_data_re = int(len(data_subcarriers))
    total_re = int(n_sc * n_symbols)
    return ResourceGrid(
        n_sc=n_sc,
        n_symbols=n_symbols,
        n_fft=n_fft,
        scs_khz=int(resource.scs_khz),
        subcarrier_indices=subcarrier_indices,
        active_fft_indices=subcarrier_indices + n_fft // 2,
        pilot_subcarriers=pilot_subcarriers,
        pilot_symbol_indices=pilot_symbols,
        pilot_subcarrier_indices=pilot_subcarrier_indices,
        data_symbol_indices=np.asarray(data_symbols, dtype=np.int64),
        data_subcarrier_indices=np.asarray(data_subcarriers, dtype=np.int64),
        cyclic_prefix_length=int(resource.cyclic_prefix_length),
        ofdm_symbol_duration_s=(1.0 + float(resource.cyclic_prefix_length) / float(n_fft))
        / (float(resource.scs_khz) * 1e3),
        dmrs_overhead=float(n_dmrs_re) / float(total_re),
        n_dmrs_re=n_dmrs_re,
        n_data_re=n_data_re,
    )


def local_indices_for_subcarriers(grid: ResourceGrid, k_values: np.ndarray) -> np.ndarray:
    first = int(grid.subcarrier_indices[0])
    return np.asarray(k_values, dtype=np.int64) - first


def full_fft_indices_for_subcarriers(grid: ResourceGrid, k_values: np.ndarray) -> np.ndarray:
    """Map centered project subcarrier numbers to Sionna's full-FFT axis."""
    indices = np.asarray(k_values, dtype=np.int64) + int(grid.n_fft) // 2
    if np.any(indices < 0) or np.any(indices >= int(grid.n_fft)):
        raise ValueError("Subcarrier lies outside the configured full FFT.")
    return indices
