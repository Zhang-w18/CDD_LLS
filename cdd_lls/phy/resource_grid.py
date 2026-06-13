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
    pilot_subcarriers: np.ndarray
    data_symbol_indices: np.ndarray
    data_subcarrier_indices: np.ndarray
    dmrs_overhead: float
    n_dmrs_re: int
    n_data_re: int

    @property
    def pilot_count(self) -> int:
        return int(len(self.pilot_subcarriers))


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

    dmrs_symbols = set(int(x) for x in resource.dmrs_symbol_indices)
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

    n_dmrs_re = int(len(dmrs_symbols) * len(pilot_subcarriers))
    n_data_re = int(len(data_subcarriers))
    total_re = int(n_sc * n_symbols)
    return ResourceGrid(
        n_sc=n_sc,
        n_symbols=n_symbols,
        n_fft=n_fft,
        scs_khz=int(resource.scs_khz),
        subcarrier_indices=subcarrier_indices,
        pilot_subcarriers=pilot_subcarriers,
        data_symbol_indices=np.asarray(data_symbols, dtype=np.int64),
        data_subcarrier_indices=np.asarray(data_subcarriers, dtype=np.int64),
        dmrs_overhead=float(n_dmrs_re) / float(total_re),
        n_dmrs_re=n_dmrs_re,
        n_data_re=n_data_re,
    )


def local_indices_for_subcarriers(grid: ResourceGrid, k_values: np.ndarray) -> np.ndarray:
    first = int(grid.subcarrier_indices[0])
    return np.asarray(k_values, dtype=np.int64) - first
