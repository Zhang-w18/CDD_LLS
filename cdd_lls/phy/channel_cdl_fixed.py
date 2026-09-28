"""Fixed-long-term-statistics 3GPP CDL channel and 32-TXRU beam codebooks."""

from __future__ import annotations

from dataclasses import asdict
from typing import Dict

import numpy as np

from cdd_lls.core.config import FixedCDLStatisticsConfig, PlatformConfig
from cdd_lls.phy.beam8 import build_dft_2x8_same_pol_codebook, load_frozen_manifest
from cdd_lls.phy.plan039 import load_manifest as load_plan039_manifest
from cdd_lls.phy.channel_tdl import TDLRealization, speed_kmh_to_mps
from cdd_lls.phy.precoding import PrecoderResult
from cdd_lls.phy.resource_grid import ResourceGrid


SCHEMES = ("BASELINE", "POLARIZATION_CYCLING", "BEAM_CYCLING", "BEAM_CDD")


def build_contiguous_txru_mapping(
    vertical_aes: int,
    horizontal_aes: int,
    vertical_txrus: int,
    horizontal_txrus: int,
) -> np.ndarray:
    """Map contiguous rectangular AE subarrays to unit-norm TXRU columns."""
    m, n = int(vertical_aes), int(horizontal_aes)
    mt, nt = int(vertical_txrus), int(horizontal_txrus)
    if min(m, n, mt, nt) <= 0 or m % mt or n % nt:
        raise ValueError("AE dimensions must be positive and divisible by TXRU dimensions.")
    block_m, block_n = m // mt, n // nt
    mapping = np.zeros((m * n, mt * nt), dtype=np.complex128)
    amplitude = 1.0 / np.sqrt(float(block_m * block_n))
    for tm in range(mt):
        for tn in range(nt):
            txru = tm * nt + tn
            for row in range(tm * block_m, (tm + 1) * block_m):
                for col in range(tn * block_n, (tn + 1) * block_n):
                    mapping[row * n + col, txru] = amplitude
    return mapping


def horizontal_dft_codebook(
    vertical_txrus: int,
    horizontal_txrus: int,
    beam_count: int,
) -> np.ndarray:
    """Return unit-norm horizontal steering beams spanning visible spatial frequency."""
    mt, nt, beam_count = int(vertical_txrus), int(horizontal_txrus), int(beam_count)
    if min(mt, nt, beam_count) <= 0:
        raise ValueError("TXRU dimensions and beam_count must be positive.")
    spatial_frequency = -1.0 + (np.arange(beam_count, dtype=np.float64) + 0.5) * 2.0 / beam_count
    horizontal = np.exp(1j * np.pi * np.arange(nt)[:, None] * spatial_frequency[None, :])
    weights = np.tile(horizontal[None, :, :], (mt, 1, 1)).reshape(mt * nt, beam_count)
    weights /= np.linalg.norm(weights, axis=0, keepdims=True)
    return weights.astype(np.complex128)


def dual_polarized(spatial: np.ndarray, phase_deg: float = 0.0) -> np.ndarray:
    spatial = np.asarray(spatial, dtype=np.complex128).reshape(-1)
    vector = np.concatenate((spatial, np.exp(1j * np.deg2rad(phase_deg)) * spatial))
    return vector / np.linalg.norm(vector)


def build_fixed_cdl_codebooks(settings: FixedCDLStatisticsConfig) -> tuple[np.ndarray, np.ndarray]:
    """Build full dual-polarized parent and secondary TXRU codebooks."""
    parent_single = horizontal_dft_codebook(
        settings.bs_vertical_txrus_per_pol,
        settings.bs_horizontal_txrus_per_pol,
        settings.ssb_horizontal_beams,
    )
    secondary_single = horizontal_dft_codebook(
        settings.bs_vertical_txrus_per_pol,
        settings.bs_horizontal_txrus_per_pol,
        settings.secondary_horizontal_beams,
    )
    parent = np.column_stack([dual_polarized(parent_single[:, i]) for i in range(parent_single.shape[1])])
    secondary = np.column_stack(
        [dual_polarized(secondary_single[:, i]) for i in range(secondary_single.shape[1])]
    )
    return parent, secondary


def build_fixed_cdl_precoder(
    grid: ResourceGrid,
    settings: FixedCDLStatisticsConfig,
    scheme: str,
    selected_ssb: int,
    parent_codebook: np.ndarray,
    secondary_codebook: np.ndarray,
    prg_size_rb: int,
) -> PrecoderResult:
    """Build one of the four fixed-CDL rank-1, unit-power precoders."""
    scheme = str(scheme).upper()
    if scheme not in SCHEMES:
        raise ValueError(f"Unsupported fixed-CDL scheme {scheme!r}.")
    n_sc = int(grid.n_sc)
    prg_sc = 12 * int(prg_size_rb)
    if prg_sc <= 0 or n_sc % prg_sc:
        raise ValueError("PRG size must be positive and divide the active allocation.")
    parent = int(selected_ssb)
    if parent < 0 or parent >= parent_codebook.shape[1]:
        raise ValueError("selected_ssb is outside the parent codebook.")
    single_count = parent_codebook.shape[0] // 2
    parent_single = np.asarray(parent_codebook[:single_count, parent], dtype=np.complex128) * np.sqrt(2.0)
    prg = np.arange(n_sc, dtype=np.int64) // prg_sc
    weights = np.empty((n_sc, parent_codebook.shape[0]), dtype=np.complex128)
    raw_power = np.ones(n_sc, dtype=np.float64)
    children = (2 * parent, 2 * parent + 1)

    if scheme == "BASELINE":
        weights[:] = parent_codebook[:, parent]
    elif scheme == "POLARIZATION_CYCLING":
        phases = np.asarray(settings.pol_cycling_phase_deg, dtype=np.float64)
        for k in range(n_sc):
            weights[k] = dual_polarized(parent_single, phases[prg[k] % phases.size])
    elif scheme == "BEAM_CYCLING":
        for k in range(n_sc):
            weights[k] = secondary_codebook[:, children[prg[k] % 2]]
    else:
        branches = secondary_codebook[:, np.asarray(children, dtype=np.int64)]
        indices = np.asarray(settings.beam_cdd_delay_grid_indices, dtype=np.float64)
        phase = np.exp(-1j * 2.0 * np.pi * np.arange(n_sc)[:, None] * indices[None, :] / n_sc)
        unnormalized = (phase @ branches.T) / np.sqrt(float(indices.size))
        raw_power = np.sum(np.abs(unnormalized) ** 2, axis=1).real
        if np.any(raw_power <= 0.0):
            raise RuntimeError("Beam-CDD produced a zero-power subcarrier.")
        weights = unnormalized / np.sqrt(raw_power)[:, None]

    power = np.sum(np.abs(weights) ** 2, axis=1).real
    if not np.allclose(power, 1.0, atol=1e-12, rtol=0.0):
        raise RuntimeError(f"Fixed-CDL precoder violates unit-power constraint: max error={np.max(np.abs(power-1))}.")
    return PrecoderResult(
        C=weights,
        label=scheme,
        metadata={
            "selected_ssb": parent,
            "secondary_children": list(children),
            "prg_size_rb": int(prg_size_rb),
            "normalization": "per_active_subcarrier_deterministic",
            "raw_power_min": float(np.min(raw_power)),
            "raw_power_max": float(np.max(raw_power)),
        },
    )


class FixedCDLStatisticsChannel:
    """Freeze CDL rays/topology and resample only initial small-scale phases."""

    def __init__(self, config: PlatformConfig, grid: ResourceGrid, prepare_reference: bool = True):
        import tensorflow as tf
        from sionna.phy import config as sionna_config
        from sionna.phy.channel.tr38901 import CDL, PanelArray, Rays, Topology

        self.config = config
        self.grid = grid
        self.settings = config.fixed_cdl_statistics
        self._tf = tf
        self._sionna_config = sionna_config
        self._topology_type = Topology
        carrier = float(config.channel.carrier_frequency_hz)
        precision = "double"

        self._bs_array = PanelArray(
            num_rows_per_panel=int(self.settings.bs_vertical_aes),
            num_cols_per_panel=int(self.settings.bs_horizontal_aes),
            polarization="dual",
            polarization_type=str(self.settings.bs_polarization_type),
            antenna_pattern=str(self.settings.bs_antenna_pattern),
            carrier_frequency=carrier,
            element_vertical_spacing=float(self.settings.bs_vertical_spacing_lambda),
            element_horizontal_spacing=float(self.settings.bs_horizontal_spacing_lambda),
            precision=precision,
        )
        ue_polarization = "dual" if int(self.settings.ue_polarizations) == 2 else "single"
        ue_pol_type = str(self.settings.ue_polarization_type) if ue_polarization == "dual" else "V"
        self._ut_array = PanelArray(
            num_rows_per_panel=int(self.settings.ue_vertical_elements),
            num_cols_per_panel=int(self.settings.ue_horizontal_elements),
            polarization=ue_polarization,
            polarization_type=ue_pol_type,
            antenna_pattern=str(self.settings.ue_antenna_pattern),
            carrier_frequency=carrier,
            element_vertical_spacing=float(self.settings.ue_vertical_spacing_lambda),
            element_horizontal_spacing=float(self.settings.ue_horizontal_spacing_lambda),
            precision=precision,
        )
        self._model = CDL(
            model=str(config.channel.cdl_profile).upper(),
            delay_spread=float(config.channel.delay_spread_ns) * 1e-9,
            carrier_frequency=carrier,
            ut_array=self._ut_array,
            bs_array=self._bs_array,
            direction="downlink",
            ut_orientation=tf.constant(np.deg2rad(self.settings.ue_orientation_deg), dtype=tf.float64),
            bs_orientation=tf.constant(np.deg2rad(self.settings.bs_orientation_deg), dtype=tf.float64),
            min_speed=speed_kmh_to_mps(float(config.channel.ue_speed_kmh)),
            max_speed=speed_kmh_to_mps(float(config.channel.ue_speed_kmh)),
            precision=precision,
        )
        self._profile_raw_arrays = {
            name: np.asarray(getattr(self._model, f"_{name}").numpy()).copy()
            for name in ("delays", "powers", "aoa", "aod", "zoa", "zod", "xpr", "k_factor")
        }
        self.angle_statistics = self._apply_angle_transform()

        single = build_contiguous_txru_mapping(
            self.settings.bs_vertical_aes,
            self.settings.bs_horizontal_aes,
            self.settings.bs_vertical_txrus_per_pol,
            self.settings.bs_horizontal_txrus_per_pol,
        )
        m, n = int(self.settings.bs_vertical_aes), int(self.settings.bs_horizontal_aes)
        self._sionna_to_row_major = np.asarray([col * m + row for row in range(m) for col in range(n)])
        self._mapping = np.zeros((2 * single.shape[0], 2 * single.shape[1]), dtype=np.complex128)
        self._mapping[: single.shape[0], : single.shape[1]] = single
        self._mapping[single.shape[0] :, single.shape[1] :] = single

        statistics_seed = int(self.settings.statistics_seed)
        tf.random.set_seed(statistics_seed % (2**31 - 1))
        sionna_config.seed = statistics_seed
        speed = speed_kmh_to_mps(float(config.channel.ue_speed_kmh))
        az = np.deg2rad(float(self.settings.velocity_azimuth_deg))
        el = np.deg2rad(float(self.settings.velocity_elevation_deg))
        velocity = [speed * np.cos(el) * np.cos(az), speed * np.cos(el) * np.sin(az), speed * np.sin(el)]
        tile3 = lambda value: tf.tile(value, [1, 1, 1])
        self._topology = Topology(
            velocities=tf.constant([[velocity]], self._model.rdtype),
            moving_end="rx",
            los_aoa=tile3(self._model._los_aoa),
            los_aod=tile3(self._model._los_aod),
            los_zoa=tile3(self._model._los_zoa),
            los_zod=tile3(self._model._los_zod),
            los=tf.fill([1, 1, 1], self._model._los),
            distance_3d=tf.zeros([1, 1, 1], self._model.rdtype),
            tx_orientations=tf.reshape(self._model._tx_orientation, [1, 1, 3]),
            rx_orientations=tf.reshape(self._model._rx_orientation, [1, 1, 3]),
        )
        tiled = lambda value, multiples: tf.tile(value, multiples)
        aoa, aod, zoa, zod = self._model._random_coupling(
            tiled(self._model._aoa, [1, 1, 1, 1, 1]),
            tiled(self._model._aod, [1, 1, 1, 1, 1]),
            tiled(self._model._zoa, [1, 1, 1, 1, 1]),
            tiled(self._model._zod, [1, 1, 1, 1, 1]),
        )
        self._rays = Rays(
            delays=tiled(self._model._delays * self._model._delay_spread, [1, 1, 1, 1]),
            powers=tiled(self._model._powers, [1, 1, 1, 1]),
            aoa=aoa,
            aod=aod,
            zoa=zoa,
            zod=zod,
            xpr=tiled(self._model._xpr, [1, 1, 1, 1, 1]),
        )
        self._k_factor = tiled(self._model._k_factor, [1, 1, 1])
        self.delays_s = np.asarray(self._rays.delays.numpy()[0, 0, 0], dtype=np.float64)
        self.powers = np.asarray(self._rays.powers.numpy()[0, 0, 0], dtype=np.float64)
        codebook_type = str(self.settings.codebook_type).lower()
        if codebook_type in {"wide_beam_split", "angular_full_coverage_ultrawide"}:
            if prepare_reference:
                manifest, arrays = load_plan039_manifest(
                    self.settings.frozen_beam_manifest, self.settings.frozen_beam_manifest_sha256
                )
                self._validate_beam8_manifest_context(manifest)
                self.frozen_beam_manifest = manifest
                self.plan039_arrays = arrays
                self.parent_codebook = arrays.get("narrow_weights", arrays.get("split_weights"))
                self.secondary_codebook = np.empty((32, 0), dtype=np.complex128)
                self.selected_beam_indices = list(range(8))
                self.ssb_long_term_powers = np.asarray(manifest["ssb_long_term_powers"], dtype=np.float64)
                self.reference_receive_power = float(manifest["reference_receive_power"])
                self.selected_ssb = int(manifest["selected_ssb_index"])
                self.transmit_covariance = np.empty((0, 0), dtype=np.complex128)
            else:
                self.parent_codebook = np.empty((32, 0), dtype=np.complex128)
                self.secondary_codebook = np.empty((32, 0), dtype=np.complex128)
        elif codebook_type == "dft_2x8_same_pol":
            self.parent_codebook = build_dft_2x8_same_pol_codebook()
            self.secondary_codebook = np.empty((self.parent_codebook.shape[0], 0), dtype=np.complex128)
            if prepare_reference:
                manifest = load_frozen_manifest(
                    self.settings.frozen_beam_manifest,
                    self.settings.frozen_beam_manifest_sha256,
                )
                self._validate_beam8_manifest_context(manifest)
                self.frozen_beam_manifest = manifest
                self.selected_beam_indices = [int(value) for value in manifest["selected_beam_indices"]]
                self.ssb_long_term_powers = np.asarray(manifest["mean_rsrp"], dtype=np.float64)
                self.reference_receive_power = float(np.max(self.ssb_long_term_powers))
                self.selected_ssb = -1
                self.transmit_covariance = np.empty((0, 0), dtype=np.complex128)
        else:
            self.parent_codebook, self.secondary_codebook = build_fixed_cdl_codebooks(self.settings)
            if prepare_reference:
                self._prepare_long_term_reference(int(self.settings.covariance_realizations))

    def _apply_angle_transform(self) -> Dict[str, Dict[str, float]]:
        tf = self._tf
        cluster_powers = np.asarray(self._model._powers.numpy(), dtype=np.float64).reshape(-1)
        statistics: Dict[str, Dict[str, float]] = {}
        for name in ("aod", "aoa", "zod", "zoa"):
            original_tensor = getattr(self._model, f"_{name}")
            values = np.rad2deg(np.asarray(original_tensor.numpy(), dtype=np.float64)).reshape(cluster_powers.size, -1)
            weights = np.broadcast_to(cluster_powers[:, None] / values.shape[1], values.shape)
            if name in {"aod", "aoa"}:
                old_mean = float(np.rad2deg(np.angle(np.sum(weights * np.exp(1j * np.deg2rad(values))))))
                offsets = (values - old_mean + 180.0) % 360.0 - 180.0
            else:
                old_mean = float(np.sum(weights * values) / np.sum(weights))
                offsets = values - old_mean
            old_spread = float(np.sqrt(np.sum(weights * offsets**2) / np.sum(weights)))
            if bool(self.settings.profile_native_angles):
                statistics[name] = {
                    "original_mean_deg": old_mean,
                    "original_rms_spread_deg": old_spread,
                    "applied_scale": 1.0,
                    "transformed_mean_deg": old_mean,
                    "transformed_rms_spread_deg": old_spread,
                }
                continue
            target = float(self.settings.mean_aod_deg) if name == "aod" else old_mean
            if name == "aod" and self.settings.target_aod_asd_deg is not None:
                scale = float(self.settings.target_aod_asd_deg) / old_spread
            else:
                scale = float(getattr(self.settings, f"{name}_scale"))
            shifted = target + scale * offsets
            if name in {"aod", "aoa"}:
                shifted = (shifted + 180.0) % 360.0 - 180.0
                provisional = float(np.rad2deg(np.angle(np.sum(weights * np.exp(1j * np.deg2rad(shifted))))))
                shifted = (shifted + target - provisional + 180.0) % 360.0 - 180.0
                new_mean = float(np.rad2deg(np.angle(np.sum(weights * np.exp(1j * np.deg2rad(shifted))))))
                new_offsets = (shifted - new_mean + 180.0) % 360.0 - 180.0
            else:
                shifted = np.clip(shifted, 0.0, 180.0)
                new_mean = float(np.sum(weights * shifted) / np.sum(weights))
                new_offsets = shifted - new_mean
            setattr(
                self._model,
                f"_{name}",
                tf.constant(np.deg2rad(shifted.reshape(original_tensor.shape)), dtype=self._model.rdtype),
            )
            statistics[name] = {
                "original_mean_deg": old_mean,
                "original_rms_spread_deg": old_spread,
                "applied_scale": scale,
                "transformed_mean_deg": new_mean,
                "transformed_rms_spread_deg": float(
                    np.sqrt(np.sum(weights * new_offsets**2) / np.sum(weights))
                ),
            }
        return statistics

    def _validate_beam8_manifest_context(self, manifest: Dict[str, object]) -> None:
        expected = {
            "cdl_profile": str(self.config.channel.cdl_profile).upper(),
            "delay_spread_ns": float(self.config.channel.delay_spread_ns),
            "carrier_frequency_hz": float(self.config.channel.carrier_frequency_hz),
            "ue_speed_kmh": float(self.config.channel.ue_speed_kmh),
            "n_sc": int(self.grid.n_sc),
            "n_symbols": int(self.grid.n_symbols),
            "statistics_seed": int(self.settings.statistics_seed),
        }
        actual = manifest.get("frozen_context")
        if actual != expected:
            raise ValueError(f"Frozen beam manifest context mismatch: expected={expected}, actual={actual}.")

    @property
    def txru_mapping(self) -> np.ndarray:
        return self._mapping.copy()

    def realization_beam_statistics(self, realization_index: int, codebook: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
        """Return realization-level R_t and all codebook powers using the statistics seed namespace."""
        response = self._response(self._coefficients(self._seed(0x434F56, int(realization_index))))
        covariance = np.einsum("nrkt,nrku->tu", response.conj(), response, optimize=True) / (
            int(self.grid.n_symbols) * int(self.grid.n_sc)
        )
        covariance = 0.5 * (covariance + covariance.conj().T)
        beams = self.parent_codebook if codebook is None else np.asarray(codebook, dtype=np.complex128)
        powers = np.einsum("tb,tu,ub->b", beams.conj(), covariance, beams, optimize=True).real
        return covariance, powers

    def profile_arrays(self) -> Dict[str, np.ndarray]:
        """Return the fixed expanded CDL arrays used by the CIR sampler."""
        expanded = {
            "delays_s": np.asarray(self._rays.delays.numpy()),
            "powers": np.asarray(self._rays.powers.numpy()),
            "aoa_rad": np.asarray(self._rays.aoa.numpy()),
            "aod_rad": np.asarray(self._rays.aod.numpy()),
            "zoa_rad": np.asarray(self._rays.zoa.numpy()),
            "zod_rad": np.asarray(self._rays.zod.numpy()),
            "xpr": np.asarray(self._rays.xpr.numpy()),
            "k_factor": np.asarray(self._k_factor.numpy()),
            "los_aoa_rad": np.asarray(self._model._los_aoa.numpy()),
            "los_aod_rad": np.asarray(self._model._los_aod.numpy()),
            "los_zoa_rad": np.asarray(self._model._los_zoa.numpy()),
            "los_zod_rad": np.asarray(self._model._los_zod.numpy()),
            "los_indicator": np.asarray(self._topology.los.numpy()),
        }
        expanded.update({f"profile_raw_{name}": value for name, value in self._profile_raw_arrays.items()})
        return expanded

    def analytic_transmit_covariance(self) -> np.ndarray:
        """Evaluate E[H^H H] by summing independent CDL ray/polarization phasors.

        Frequency- and time-dependent unit-magnitude phases cancel in each
        self outer product. Cross terms vanish because Sionna draws the four
        polarization phases independently for every cluster/ray.
        """
        tf = self._tf
        generator = self._model._cir_sampler
        rays = self._rays
        array = generator._step_11_array_offsets(
            self._topology, rays.aoa, rays.aod, rays.zoa, rays.zod
        )
        xpr_scale = tf.sqrt(1.0 / rays.xpr)
        shape = tuple(int(value) for value in rays.xpr.shape) + (2, 2)
        delays, covariances, _ = self.analytic_delay_covariances()
        covariance = np.sum(covariances, axis=0)
        covariance = 0.5 * (covariance + covariance.conj().T)
        if not np.all(np.isfinite(covariance)):
            raise RuntimeError("Analytic CDL covariance contains non-finite values.")
        return covariance

    def analytic_delay_covariances(self) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
        """Return analytic TXRU covariance contributions before equal-delay merging."""
        tf = self._tf
        generator = self._model._cir_sampler
        rays = self._rays
        array = generator._step_11_array_offsets(
            self._topology, rays.aoa, rays.aod, rays.zoa, rays.zod
        )
        xpr_scale = tf.sqrt(1.0 / rays.xpr)
        shape = tuple(int(value) for value in rays.xpr.shape) + (2, 2)
        cluster_count = int(np.asarray(rays.powers).reshape(-1).size)
        covariances = np.zeros((cluster_count, self._mapping.shape[1], self._mapping.shape[1]), dtype=np.complex128)
        per_pol = int(self._bs_array.num_ant) // 2
        order = np.r_[self._sionna_to_row_major, per_pol + self._sionna_to_row_major]
        los = bool(np.asarray(self._topology.los.numpy()).reshape(-1)[0])
        k_factor = float(np.asarray(self._k_factor.numpy()).reshape(-1)[0])
        nlos_scale = 1.0 / (k_factor + 1.0) if los else 1.0
        ray_scale = tf.sqrt(rays.powers / tf.cast(tf.shape(rays.xpr)[4], self._model.rdtype))
        for component in range(4):
            matrix = np.zeros(shape, dtype=np.complex128)
            row, col = divmod(component, 2)
            amplitude = np.asarray(xpr_scale.numpy()) if component in (1, 2) else 1.0
            matrix[..., row, col] = amplitude
            field = generator._step_11_field_matrix(
                self._topology, rays.aoa, rays.aod, rays.zoa, rays.zod,
                tf.constant(matrix, dtype=self._model.cdtype),
            )
            contribution = np.asarray((field * array).numpy()[0, 0, 0], dtype=np.complex128)
            contribution *= np.asarray(
                ray_scale.numpy()[0, 0, 0], dtype=np.float64
            )[:, None, None, None]
            contribution = contribution[..., order]
            txru = np.einsum("cmra,at->cmrt", contribution, self._mapping, optimize=True)
            covariances += nlos_scale * np.einsum("cmrt,cmru->ctu", txru.conj(), txru, optimize=True)
        delays = np.asarray(rays.delays.numpy()[0, 0, 0], dtype=np.float64).reshape(-1)
        component_kind = ["diffuse_cluster"] * cluster_count
        if los:
            t = tf.zeros([1], self._model.rdtype)
            direct = np.asarray(generator._step_11_los(self._topology, t).numpy()[0, 0, 0, 0, :, :, 0])
            direct = direct[:, order]
            txru = direct @ self._mapping
            specular = (k_factor / (k_factor + 1.0)) * txru.conj().T @ txru
            covariances = np.concatenate((covariances, specular[None, ...]), axis=0)
            delays = np.r_[delays, delays[0]]
            component_kind.append("los_specular")
        covariances = 0.5 * (covariances + covariances.conj().transpose(0, 2, 1))
        if not np.all(np.isfinite(covariances)):
            raise RuntimeError("Analytic CDL delay covariance contains non-finite values.")
        return delays, covariances, {
            "component_kind": component_kind,
            "los": los,
            "k_factor_linear": k_factor,
            "diffuse_scale": nlos_scale,
        }

    def analytic_ray_covariances(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return delay, Doppler, and TXRU covariance for each diffuse ray/LoS component."""
        tf = self._tf
        generator = self._model._cir_sampler
        rays = self._rays
        array = generator._step_11_array_offsets(self._topology, rays.aoa, rays.aod, rays.zoa, rays.zod)
        xpr_scale = tf.sqrt(1.0 / rays.xpr)
        shape = tuple(int(value) for value in rays.xpr.shape) + (2, 2)
        cluster_count = int(np.asarray(rays.powers).reshape(-1).size)
        ray_count = int(rays.xpr.shape[-1])
        covariance = np.zeros((cluster_count, ray_count, self._mapping.shape[1], self._mapping.shape[1]), dtype=np.complex128)
        per_pol = int(self._bs_array.num_ant) // 2
        order = np.r_[self._sionna_to_row_major, per_pol + self._sionna_to_row_major]
        los = bool(np.asarray(self._topology.los.numpy()).reshape(-1)[0])
        k_factor = float(np.asarray(self._k_factor.numpy()).reshape(-1)[0])
        nlos_scale = 1.0 / (k_factor + 1.0) if los else 1.0
        ray_scale = tf.sqrt(rays.powers / tf.cast(tf.shape(rays.xpr)[4], self._model.rdtype))
        for component in range(4):
            matrix = np.zeros(shape, dtype=np.complex128)
            row, col = divmod(component, 2)
            amplitude = np.asarray(xpr_scale.numpy()) if component in (1, 2) else 1.0
            matrix[..., row, col] = amplitude
            field = generator._step_11_field_matrix(
                self._topology, rays.aoa, rays.aod, rays.zoa, rays.zod,
                tf.constant(matrix, dtype=self._model.cdtype),
            )
            contribution = np.asarray((field * array).numpy()[0, 0, 0], dtype=np.complex128)
            contribution *= np.asarray(ray_scale.numpy()[0, 0, 0], dtype=np.float64)[:, None, None, None]
            txru = np.einsum("cmra,at->cmrt", contribution[..., order], self._mapping, optimize=True)
            covariance += nlos_scale * np.einsum("cmrt,cmru->cmtu", txru.conj(), txru, optimize=True)
        delays = np.repeat(np.asarray(rays.delays.numpy()[0, 0, 0], dtype=np.float64), ray_count)
        aoa = np.asarray(rays.aoa.numpy()[0, 0, 0], dtype=np.float64)
        zoa = np.asarray(rays.zoa.numpy()[0, 0, 0], dtype=np.float64)
        velocity = np.asarray(self._topology.velocities.numpy()[0, 0], dtype=np.float64)
        direction = np.stack((np.sin(zoa) * np.cos(aoa), np.sin(zoa) * np.sin(aoa), np.cos(zoa)), axis=-1)
        wavelength = 299792458.0 / float(self.config.channel.carrier_frequency_hz)
        dopplers = np.einsum("cmv,v->cm", direction, velocity).reshape(-1) / wavelength
        covariance = covariance.reshape(-1, covariance.shape[-2], covariance.shape[-1])
        if los:
            t = tf.zeros([1], self._model.rdtype)
            direct = np.asarray(generator._step_11_los(self._topology, t).numpy()[0, 0, 0, 0, :, :, 0])
            txru = direct[:, order] @ self._mapping
            covariance = np.concatenate((covariance, ((k_factor / (k_factor + 1.0)) * txru.conj().T @ txru)[None]), axis=0)
            los_aoa = float(np.asarray(self._model._los_aoa.numpy()).reshape(-1)[0])
            los_zoa = float(np.asarray(self._model._los_zoa.numpy()).reshape(-1)[0])
            los_direction = np.asarray([np.sin(los_zoa) * np.cos(los_aoa), np.sin(los_zoa) * np.sin(los_aoa), np.cos(los_zoa)])
            dopplers = np.r_[dopplers, float(los_direction @ velocity / wavelength)]
            delays = np.r_[delays, delays[0]]
        return delays, dopplers, covariance

    def _coefficients(
        self,
        seed: int,
        num_time_samples: int | None = None,
        sampling_frequency_hz: float | None = None,
        topology=None,
    ) -> np.ndarray:
        seed = int(seed)
        self._tf.random.set_seed(seed % (2**31 - 1))
        self._sionna_config.seed = seed
        time_samples = int(self.grid.n_symbols if num_time_samples is None else num_time_samples)
        sampling_frequency = float(
            1.0 / float(self.grid.ofdm_symbol_duration_s)
            if sampling_frequency_hz is None
            else sampling_frequency_hz
        )
        if time_samples <= 0 or not np.isfinite(sampling_frequency) or sampling_frequency <= 0.0:
            raise ValueError("Fixed-CDL time sampling must be finite and positive.")
        h, _ = self._model._cir_sampler(
            time_samples,
            sampling_frequency,
            self._k_factor,
            self._rays,
            self._topology if topology is None else topology,
        )
        h = self._tf.transpose(h, [0, 2, 4, 1, 5, 3, 6])
        ae = np.asarray(h.numpy()[0, 0, :, 0, :, :, :], dtype=np.complex128)
        per_pol = ae.shape[1] // 2
        order = np.r_[self._sionna_to_row_major, per_pol + self._sionna_to_row_major]
        return np.einsum("raln,at->nrtl", ae[:, order, :, :], self._mapping, optimize=True)

    def _response(self, coefficients: np.ndarray) -> np.ndarray:
        offsets_hz = np.asarray(self.grid.subcarrier_indices, dtype=np.float64) * float(self.grid.scs_khz) * 1e3
        phase = np.exp(-1j * 2.0 * np.pi * offsets_hz[:, None] * self.delays_s[None, :])
        return np.einsum("nrtl,kl->nrkt", coefficients, phase, optimize=True)

    def _seed(self, namespace: int, index: int) -> int:
        root = int(self.settings.statistics_seed if namespace == 0x434F56 else self.settings.realization_seed)
        return int(np.random.SeedSequence([root, namespace, int(index)]).generate_state(1)[0])

    def _prepare_long_term_reference(self, count: int) -> None:
        n_tx = self._mapping.shape[1]
        covariance = np.zeros((n_tx, n_tx), dtype=np.complex128)
        for index in range(int(count)):
            response = self._response(self._coefficients(self._seed(0x434F56, index)))
            covariance += np.einsum("nrkt,nrku->tu", response.conj(), response, optimize=True) / (
                int(self.grid.n_symbols) * int(self.grid.n_sc)
            )
        covariance /= float(count)
        self.transmit_covariance = 0.5 * (covariance + covariance.conj().T)
        self.ssb_long_term_powers = np.einsum(
            "tb,tu,ub->b",
            self.parent_codebook.conj(),
            self.transmit_covariance,
            self.parent_codebook,
            optimize=True,
        ).real
        self.selected_ssb = int(np.argmax(self.ssb_long_term_powers))
        self.reference_receive_power = float(self.ssb_long_term_powers[self.selected_ssb])
        if not np.isfinite(self.reference_receive_power) or self.reference_receive_power <= 0.0:
            raise RuntimeError("Fixed-CDL reference receive power must be finite and positive.")

    def generate(self, realization_index: int) -> TDLRealization:
        coefficients = self._coefficients(self._seed(0x43444C, int(realization_index)))
        response = self._response(coefficients)
        h = np.transpose(response, (1, 3, 0, 2))[None, ...]
        return TDLRealization(
            taps=np.empty((0,), dtype=np.complex128),
            H=h,
            pdp=self.powers.copy(),
            tap_delays=self.delays_s.copy(),
            sample_period_ns=float("nan"),
            backend="fixed_cdl_statistics",
            metadata=self.metadata(),
        )

    def generate_with_aged_csi(
        self,
        realization_index: int,
        age_s: float,
    ) -> tuple[TDLRealization, np.ndarray, float]:
        """Return the current slot and exact same-realization CSI from ``age_s`` earlier.

        The current slot is sampled at the normal OFDM-symbol spacing.  The stale
        sample is evaluated at negative time by replaying the same CDL random
        phases with reversed velocity and taking the sample at ``+age_s``.  This
        avoids generating every intermediate OFDM symbol over a long CSI age.
        ``old_channel`` has shape ``[1,n_rx,n_tx,n_sc]``.
        """
        age = float(age_s)
        if not np.isfinite(age) or age <= 0.0:
            raise ValueError("age_s must be finite and positive.")
        seed = self._seed(0x43444C, int(realization_index))
        current_coefficients = self._coefficients(seed)
        reverse_topology = self._topology_type(
            velocities=-self._topology.velocities,
            moving_end="rx",
            los_aoa=self._topology.los_aoa,
            los_aod=self._topology.los_aod,
            los_zoa=self._topology.los_zoa,
            los_zod=self._topology.los_zod,
            los=self._topology.los,
            distance_3d=self._topology.distance_3d,
            tx_orientations=self._topology.tx_orientations,
            rx_orientations=self._topology.rx_orientations,
        )
        reverse_pair = self._coefficients(
            seed,
            num_time_samples=2,
            sampling_frequency_hz=1.0 / age,
            topology=reverse_topology,
        )
        replay_error = float(np.max(np.abs(reverse_pair[0] - current_coefficients[0])))
        if replay_error > 1e-10:
            raise RuntimeError(
                "Fixed-CDL aged-CSI replay changed the t=0 channel: "
                f"max error={replay_error}."
            )
        current_response = self._response(current_coefficients)
        old_response = self._response(reverse_pair[-1:])
        current_h = np.transpose(current_response, (1, 3, 0, 2))[None, ...]
        old_h = np.transpose(old_response, (1, 3, 0, 2))[None, ..., 0, :]
        return (
            TDLRealization(
                taps=np.empty((0,), dtype=np.complex128),
                H=current_h,
                pdp=self.powers.copy(),
                tap_delays=self.delays_s.copy(),
                sample_period_ns=float("nan"),
                backend="fixed_cdl_statistics",
                metadata={**self.metadata(), "aged_csi_age_s": age},
            ),
            old_h,
            replay_error,
        )

    def metadata(self) -> Dict[str, object]:
        return {
            "backend": "fixed_cdl_statistics",
            "cdl_profile": str(self.config.channel.cdl_profile).upper(),
            "delay_spread_ns": float(self.config.channel.delay_spread_ns),
            "carrier_frequency_hz": float(self.config.channel.carrier_frequency_hz),
            "ue_speed_kmh": float(self.config.channel.ue_speed_kmh),
            "long_term_parameters_frozen": True,
            "statistics_and_bler_seeds_disjoint": True,
            "selected_ssb": int(getattr(self, "selected_ssb", -1)),
            "selected_beam_indices": list(getattr(self, "selected_beam_indices", [])),
            "reference_receive_power": float(getattr(self, "reference_receive_power", float("nan"))),
            "ssb_long_term_powers": np.asarray(getattr(self, "ssb_long_term_powers", [])).tolist(),
            "angle_statistics": self.angle_statistics,
            "settings": asdict(self.settings),
            "txru_mapping": "contiguous_equal_phase_column_normalized",
            "parent_codebook": "8-beam horizontal visible-region DFT grid",
            "secondary_codebook": "16-beam horizontal visible-region oversampled DFT grid",
        }
