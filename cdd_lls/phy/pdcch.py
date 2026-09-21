from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from cdd_lls.phy.precoding import PrecoderResult, spatial_dft_codebook
from cdd_lls.phy.resource_grid import ResourceGrid


PDCCH_DATA_OFFSETS = np.asarray([0, 2, 3, 4, 6, 7, 8, 10, 11], dtype=np.int64)
PDCCH_DMRS_OFFSETS = np.asarray([1, 5, 9], dtype=np.int64)
SUPPORTED_AGGREGATION_LEVELS = (1, 2, 4, 8, 16)


@dataclass(frozen=True)
class PDCCHResourceConfig:
    n_rb: int = 48
    duration_symbols: int = 2
    n_fft: int = 4096
    scs_khz: int = 30
    cyclic_prefix_length: int = 288
    cce_reg_mapping: str = "noninterleaved"
    reg_bundle_size: int = 6
    interleaver_size: int = 2
    shift_index: int = 0
    aggregation_level: int = 4
    first_cce: int = 0


@dataclass(frozen=True)
class PDCCHGrid:
    resource_grid: ResourceGrid
    config: PDCCHResourceConfig
    candidate_cces: np.ndarray
    candidate_regs: np.ndarray
    candidate_bundles: np.ndarray
    bundle_by_rb: np.ndarray
    data_bundle_indices: np.ndarray
    dmrs_bundle_indices: np.ndarray
    data_local_subcarrier_indices: np.ndarray
    dmrs_local_subcarrier_indices: np.ndarray

    @property
    def n_cce(self) -> int:
        return int(self.config.n_rb * self.config.duration_symbols // 6)

    @property
    def n_reg(self) -> int:
        return int(self.config.n_rb * self.config.duration_symbols)

    @property
    def coded_bits(self) -> int:
        return int(2 * self.resource_grid.n_data_re)


def _validate_resource_config(config: PDCCHResourceConfig) -> None:
    n_rb = int(config.n_rb)
    duration = int(config.duration_symbols)
    n_reg = n_rb * duration
    mapping = str(config.cce_reg_mapping).lower()
    bundle = int(config.reg_bundle_size)
    aggregation = int(config.aggregation_level)
    first_cce = int(config.first_cce)

    if n_rb <= 0 or n_rb % 6:
        raise ValueError("PDCCH CORESET n_rb must be a positive multiple of 6.")
    if duration not in (1, 2, 3):
        raise ValueError("PDCCH CORESET duration_symbols must be 1, 2, or 3.")
    if int(config.n_fft) <= 12 * n_rb:
        raise ValueError("n_fft must exceed the active CORESET subcarrier count.")
    if int(config.cyclic_prefix_length) < 0:
        raise ValueError("cyclic_prefix_length must be non-negative.")
    if n_reg % 6:
        raise ValueError("The CORESET must contain an integer number of CCEs.")
    if aggregation not in SUPPORTED_AGGREGATION_LEVELS:
        raise ValueError(f"aggregation_level must be one of {SUPPORTED_AGGREGATION_LEVELS}.")
    if first_cce < 0 or first_cce + aggregation > n_reg // 6:
        raise ValueError("The configured PDCCH candidate lies outside the CORESET.")

    if mapping == "noninterleaved":
        if bundle != 6:
            raise ValueError("Non-interleaved PDCCH mapping requires reg_bundle_size=6.")
    elif mapping == "interleaved":
        allowed = (2, 6) if duration == 1 else (duration, 6)
        if bundle not in allowed:
            raise ValueError(
                f"Interleaved duration-{duration} mapping requires reg_bundle_size in {allowed}."
            )
        interleaver = int(config.interleaver_size)
        if interleaver not in (2, 3, 6):
            raise ValueError("interleaver_size must be 2, 3, or 6.")
        if n_reg % (bundle * interleaver):
            raise ValueError("N_REG/(L*R) must be an integer for interleaved mapping.")
    else:
        raise ValueError("cce_reg_mapping must be noninterleaved or interleaved.")

    if bundle % duration:
        raise ValueError(
            "This implementation requires a REG bundle to contain complete RBs across all CORESET symbols."
        )


def _interleaved_bundle_index(x: int, n_bundles: int, interleaver_size: int, shift_index: int) -> int:
    rows = int(interleaver_size)
    columns = int(n_bundles // rows)
    r = int(x % rows)
    c = int(x // rows)
    return int((r * columns + c + int(shift_index)) % n_bundles)


def cce_bundle_indices(config: PDCCHResourceConfig, cce_index: int) -> np.ndarray:
    _validate_resource_config(config)
    cce_index = int(cce_index)
    n_reg = int(config.n_rb) * int(config.duration_symbols)
    n_cce = n_reg // 6
    if cce_index < 0 or cce_index >= n_cce:
        raise ValueError("cce_index lies outside the CORESET.")
    bundle = int(config.reg_bundle_size)
    per_cce = 6 // bundle
    x_values = np.arange(cce_index * per_cce, (cce_index + 1) * per_cce, dtype=np.int64)
    if str(config.cce_reg_mapping).lower() == "noninterleaved":
        return x_values
    n_bundles = n_reg // bundle
    return np.asarray(
        [
            _interleaved_bundle_index(
                int(x), n_bundles, int(config.interleaver_size), int(config.shift_index)
            )
            for x in x_values
        ],
        dtype=np.int64,
    )


def build_pdcch_grid(config: PDCCHResourceConfig) -> PDCCHGrid:
    _validate_resource_config(config)
    n_rb = int(config.n_rb)
    duration = int(config.duration_symbols)
    n_sc = 12 * n_rb
    n_reg = n_rb * duration
    bundle_size = int(config.reg_bundle_size)

    candidate_cces = np.arange(
        int(config.first_cce),
        int(config.first_cce) + int(config.aggregation_level),
        dtype=np.int64,
    )
    cce_bundles = [cce_bundle_indices(config, int(cce)) for cce in candidate_cces]
    candidate_bundles = np.concatenate(cce_bundles).astype(np.int64)
    candidate_regs = np.concatenate(
        [
            np.arange(bundle * bundle_size, (bundle + 1) * bundle_size, dtype=np.int64)
            for bundle in candidate_bundles
        ]
    )
    if len(np.unique(candidate_regs)) != 6 * int(config.aggregation_level):
        raise RuntimeError("PDCCH CCE-to-REG mapping produced overlapping REGs.")

    data_rows: list[tuple[int, int, int]] = []
    dmrs_rows: list[tuple[int, int, int]] = []
    for reg in candidate_regs:
        rb = int(reg // duration)
        symbol = int(reg % duration)
        physical_bundle = int(reg // bundle_size)
        for offset in PDCCH_DATA_OFFSETS:
            data_rows.append((symbol, rb * 12 + int(offset), physical_bundle))
        for offset in PDCCH_DMRS_OFFSETS:
            dmrs_rows.append((symbol, rb * 12 + int(offset), physical_bundle))

    data_rows.sort(key=lambda row: (row[0], row[1]))
    dmrs_rows.sort(key=lambda row: (row[0], row[1]))
    data_array = np.asarray(data_rows, dtype=np.int64)
    dmrs_array = np.asarray(dmrs_rows, dtype=np.int64)

    half = n_sc // 2
    subcarrier_indices = np.arange(-half, half, dtype=np.int64)
    data_symbols = data_array[:, 0]
    data_local = data_array[:, 1]
    data_subcarriers = subcarrier_indices[data_local]
    dmrs_symbols = dmrs_array[:, 0]
    dmrs_local = dmrs_array[:, 1]
    dmrs_subcarriers = subcarrier_indices[dmrs_local]

    grid = ResourceGrid(
        n_sc=n_sc,
        n_symbols=duration,
        n_fft=int(config.n_fft),
        scs_khz=int(config.scs_khz),
        subcarrier_indices=subcarrier_indices,
        active_fft_indices=subcarrier_indices + int(config.n_fft) // 2,
        pilot_subcarriers=np.unique(dmrs_subcarriers),
        pilot_symbol_indices=dmrs_symbols,
        pilot_subcarrier_indices=dmrs_subcarriers,
        data_symbol_indices=data_symbols,
        data_subcarrier_indices=data_subcarriers,
        cyclic_prefix_length=int(config.cyclic_prefix_length),
        ofdm_symbol_duration_s=(
            1.0 + float(config.cyclic_prefix_length) / float(config.n_fft)
        )
        / (float(config.scs_khz) * 1e3),
        dmrs_overhead=float(len(dmrs_rows)) / float(len(dmrs_rows) + len(data_rows)),
        n_dmrs_re=int(len(dmrs_rows)),
        n_data_re=int(len(data_rows)),
    )

    bundle_by_rb = np.empty(n_rb, dtype=np.int64)
    for rb in range(n_rb):
        reg_ids = rb * duration + np.arange(duration, dtype=np.int64)
        physical = np.unique(reg_ids // bundle_size)
        if len(physical) != 1:
            raise RuntimeError("A physical RB crosses PDCCH REG-bundle boundaries.")
        bundle_by_rb[rb] = int(physical[0])

    return PDCCHGrid(
        resource_grid=grid,
        config=config,
        candidate_cces=candidate_cces,
        candidate_regs=candidate_regs,
        candidate_bundles=candidate_bundles,
        bundle_by_rb=bundle_by_rb,
        data_bundle_indices=data_array[:, 2],
        dmrs_bundle_indices=dmrs_array[:, 2],
        data_local_subcarrier_indices=data_local,
        dmrs_local_subcarrier_indices=dmrs_local,
    )


def build_reg_bundle_dft_precoder(
    pdcch_grid: PDCCHGrid,
    n_tx: int,
    cycling_order: Sequence[int] | None = None,
) -> PrecoderResult:
    n_tx = int(n_tx)
    if n_tx <= 0:
        raise ValueError("n_tx must be positive.")
    order = list(range(n_tx)) if cycling_order is None else [int(x) for x in cycling_order]
    if not order:
        raise ValueError("cycling_order must not be empty.")
    if any(x < 0 or x >= n_tx for x in order):
        raise ValueError("cycling_order entries must lie in [0,n_tx).")

    codebook = spatial_dft_codebook(n_tx, normalize=True)
    grid = pdcch_grid.resource_grid
    matrix = np.empty((grid.n_sc, n_tx), dtype=np.complex128)
    vector_by_bundle: dict[int, int] = {}
    for rb, physical_bundle in enumerate(pdcch_grid.bundle_by_rb):
        vector_index = order[int(physical_bundle) % len(order)]
        vector_by_bundle[int(physical_bundle)] = int(vector_index)
        matrix[rb * 12 : (rb + 1) * 12, :] = codebook[:, vector_index][None, :]

    return PrecoderResult(
        C=matrix,
        label="PDCCH_REG_BUNDLE_DFT_CYCLING",
        metadata={
            "reg_bundle_size": int(pdcch_grid.config.reg_bundle_size),
            "cycling_order": order,
            "vector_by_physical_bundle": vector_by_bundle,
            "normalized": True,
            "vector_power": 1.0,
        },
    )


def pdcch_dmrs_symbols(
    pdcch_grid: PDCCHGrid,
    scrambling_id: int,
    slot_number: int = 0,
    coreset_start_symbol: int = 0,
    symbols_per_slot: int = 14,
) -> np.ndarray:
    from sionna.phy.nr import generate_prng_seq

    config = pdcch_grid.config
    n_rb = int(config.n_rb)
    by_coordinate: dict[tuple[int, int], complex] = {}
    for local_symbol in range(int(config.duration_symbols)):
        symbol_in_slot = int(coreset_start_symbol) + local_symbol
        c_init = (
            (1 << 17)
            * (
                (int(symbols_per_slot) * int(slot_number) + symbol_in_slot + 1)
                * (2 * int(scrambling_id) + 1)
                + 2 * int(scrambling_id)
            )
        ) % (1 << 31)
        bits = generate_prng_seq(6 * n_rb, int(c_init)).astype(np.int8)
        values = (
            (1.0 - 2.0 * bits[0::2]) + 1j * (1.0 - 2.0 * bits[1::2])
        ) / np.sqrt(2.0)
        for rb in range(n_rb):
            for pos, offset in enumerate(PDCCH_DMRS_OFFSETS):
                by_coordinate[(local_symbol, rb * 12 + int(offset))] = values[3 * rb + pos]

    return np.asarray(
        [
            by_coordinate[(int(symbol), int(local_sc))]
            for symbol, local_sc in zip(
                pdcch_grid.resource_grid.pilot_symbol_indices,
                pdcch_grid.dmrs_local_subcarrier_indices,
            )
        ],
        dtype=np.complex128,
    )
