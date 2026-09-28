from __future__ import annotations

from pathlib import Path
from typing import Dict, List
import copy
import datetime as dt
import math
import numpy as np

from cdd_lls.core.config import (
    PlatformConfig,
    config_from_dict,
    dataclass_to_dict,
    merged_config_dict,
    save_resolved_config,
)
from cdd_lls.core.mcs import build_tb_layout, get_mcs
from cdd_lls.phy.channel_tdl import generate_sionna_channel, generate_tdl_channel
from cdd_lls.phy.channel_cdl_fixed import (
    FixedCDLStatisticsChannel,
    build_fixed_cdl_precoder,
)
from cdd_lls.phy.beam8 import DELAYS as BEAM8_DELAYS, build_beam8_precoder, frequency_covariance
from cdd_lls.phy.plan039 import build_precoder as build_plan039_precoder, frequency_covariance as plan039_frequency_covariance
from cdd_lls.phy.estimators import (
    build_time_frequency_rmmse_filter,
    build_prg_time_frequency_rmmse_filter,
    cdl_spatial_unaware_covariance,
    estimate_channel,
    tdl_known_delay_covariance,
    tdl_unknown_delay_covariance,
    TDLTimeFrequencyCovariance,
)
from cdd_lls.phy.ldpc import SionnaLDPCAdapter
from cdd_lls.phy.precoding import (
    build_aged_mrt_prg_precoder,
    build_precoder,
    equivalent_channel,
    normalize_delay_vector,
)
from cdd_lls.phy.qam import qam_demapper_maxlog, qam_modulate
from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
from cdd_lls.sim.stats import interpolate_target_snr, save_csv, save_json, snr_values
from cdd_lls.utils.plotting import plot_bler, plot_nmse


def construct_ls_observations(
    true_pilot_channel: np.ndarray,
    noise_variance: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Construct z=g+w after division by unit-power known DMRS symbols."""
    true_pilot = np.asarray(true_pilot_channel, dtype=np.complex128)
    noise = math.sqrt(float(noise_variance) / 2.0) * (
        rng.normal(size=true_pilot.shape) + 1j * rng.normal(size=true_pilot.shape)
    )
    return true_pilot + noise


class CDDLinkLevelOrchestrator:
    def __init__(self, config: PlatformConfig):
        self.cfg = config
        self.summary_rows: List[Dict[str, object]] = []
        self.trial_rows: List[Dict[str, object]] = []
        self._fixed_cdl_contexts: Dict[str, FixedCDLStatisticsChannel] = {}

    def run(self) -> Path:
        out = self._make_output_dir()
        save_resolved_config(self.cfg, out / "resolved_config.yaml")
        run_specs = self._expanded_run_specs()
        for spec in run_specs:
            self._run_spec(spec)
        for index, context in enumerate(self._fixed_cdl_contexts.values()):
            prefix = f"fixed_cdl_statistics_{index}"
            save_json(context.metadata(), out / f"{prefix}.json")
            np.savez_compressed(
                out / f"{prefix}.npz",
                transmit_covariance=context.transmit_covariance,
                ssb_long_term_powers=context.ssb_long_term_powers,
                parent_codebook=context.parent_codebook,
                secondary_codebook=context.secondary_codebook,
            )
        save_csv(self.summary_rows, out / "summary.csv")
        save_json(self.summary_rows, out / "summary.json")
        save_csv(self._target_rows(), out / "summary_10pct_bler_snr.csv")
        if self.trial_rows:
            save_csv(self.trial_rows, out / "trial_metrics.csv")
        if bool(self.cfg.plots.enabled):
            plot_bler(self.summary_rows, out / "bler_curves.png")
            plot_nmse(self.summary_rows, out / "ce_nmse_curves.png")
        return out

    def _make_output_dir(self) -> Path:
        root = Path(self.cfg.simulation.output_dir)
        stamp = (
            str(self.cfg.simulation.run_id)
            if self.cfg.simulation.run_id
            else dt.datetime.now().strftime("sim_%Y%m%d_%H%M%S")
        )
        out = root / stamp
        out.mkdir(parents=True, exist_ok=True)
        return out

    def _expanded_run_specs(self) -> List[PlatformConfig]:
        base = dataclass_to_dict(self.cfg)
        scenarios = self.cfg.scenarios or [{"scenario_id": "default"}]
        variants = self.cfg.variants or [{"variant_id": "single"}]
        delay_sweep = list(self.cfg.sweeps.cdd_base_delays or [None])
        dmrs_sweep = list(self.cfg.sweeps.dmrs_spacing_sc or [None])

        specs = []
        for scenario in scenarios:
            for variant in variants:
                for delay in delay_sweep:
                    for dmrs in dmrs_sweep:
                        merged = copy.deepcopy(base)
                        merged = merged_config_dict(config_from_dict(merged), scenario)
                        merged = merged_config_dict(config_from_dict(merged), variant)
                        if delay is not None:
                            merged["transmission"]["cdd_base_delay"] = int(delay)
                            merged["transmission"]["cdd_delay_vector"] = None
                        if dmrs is not None:
                            merged["resource"]["dmrs_spacing_sc"] = int(dmrs)
                        cfg = config_from_dict(merged)
                        sid = str(scenario.get("scenario_id", scenario.get("id", "scenario")))
                        vid = str(variant.get("variant_id", variant.get("id", "variant")))
                        if delay is not None:
                            vid += f"_d{int(delay)}"
                        if dmrs is not None:
                            vid += f"_sf{int(dmrs)}"
                        cfg._scenario_id = sid  # type: ignore[attr-defined]
                        cfg._variant_id = vid  # type: ignore[attr-defined]
                        specs.append(cfg)
        return specs

    def _run_spec(self, cfg: PlatformConfig) -> None:
        scenario_id = str(getattr(cfg, "_scenario_id", "default"))
        variant_id = str(getattr(cfg, "_variant_id", "single"))
        grid = build_resource_grid(cfg.resource)
        mcs = get_mcs(cfg.mcs.table, cfg.mcs.index, cfg.mcs.qm, cfg.mcs.code_rate)
        tb = build_tb_layout(grid.n_data_re, mcs)
        ce_only = bool(getattr(cfg.simulation, "ce_only", False))
        if ce_only and str(cfg.channel_estimation.ce_method).upper() == "IDEAL":
            raise ValueError("simulation.ce_only requires a non-IDEAL channel estimator.")
        if ce_only and str(cfg.channel.backend).lower() != "fixed_cdl_statistics":
            raise ValueError("simulation.ce_only is currently supported only by fixed_cdl_statistics.")
        adapter = None if ce_only else SionnaLDPCAdapter(
            tb.cb_k_values,
            tb.cb_e_values,
            num_iter=int(cfg.receiver.max_ldpc_iterations),
            llr_clip=float(cfg.receiver.llr_clip),
        )

        tx_scheme = str(cfg.transmission.tx_scheme).upper()
        if str(cfg.channel.backend).lower() == "fixed_cdl_statistics":
            if str(cfg.fixed_cdl_statistics.codebook_type).lower() == "dft_2x8_same_pol":
                delays = list(BEAM8_DELAYS.get(tx_scheme, ()))
            else:
                delays = list(cfg.fixed_cdl_statistics.beam_cdd_delay_grid_indices)
        else:
            delays = normalize_delay_vector(
                cfg.transmission.cdd_delay_vector,
                int(cfg.antenna.n_tx),
                int(cfg.transmission.cdd_base_delay),
            )

        configured_snr = list(getattr(cfg.simulation, "snr_points_db", []) or [])
        for snr_db in (configured_snr if configured_snr else snr_values(cfg.simulation.snr_range_db)):
            print(
                f"[progress] start scenario={scenario_id} variant={variant_id} "
                f"snr_db={float(snr_db):g} trials={int(cfg.simulation.n_trials_per_snr)}",
                flush=True,
            )
            row = self._run_snr(
                cfg=cfg,
                grid=grid,
                adapter=adapter,
                scenario_id=scenario_id,
                variant_id=variant_id,
                tx_scheme=tx_scheme,
                delays=delays,
                snr_db=float(snr_db),
                tb=tb,
                qm=int(mcs.qm),
            )
            self.summary_rows.append(row)
            print(
                f"[progress] done scenario={scenario_id} variant={variant_id} "
                f"snr_db={float(snr_db):g} errors={int(row['tb_errors'])}/{int(row['n_trials'])}",
                flush=True,
            )

    def _run_snr(
        self,
        cfg: PlatformConfig,
        grid,
        adapter: SionnaLDPCAdapter,
        scenario_id: str,
        variant_id: str,
        tx_scheme: str,
        delays: List[int],
        snr_db: float,
        tb,
        qm: int,
    ) -> Dict[str, object]:
        if str(getattr(cfg.channel, "backend", "legacy_exponential")).lower() == "fixed_cdl_statistics":
            return self._run_snr_fixed_cdl(
                cfg=cfg,
                grid=grid,
                adapter=adapter,
                scenario_id=scenario_id,
                variant_id=variant_id,
                tx_scheme=tx_scheme,
                delays=delays,
                snr_db=snr_db,
                tb=tb,
                qm=qm,
            )
        if str(getattr(cfg.channel, "backend", "legacy_exponential")).lower() in (
            "sionna_tdl",
            "sionna_cdl",
        ):
            return self._run_snr_sionna_tdl(
                cfg=cfg,
                grid=grid,
                adapter=adapter,
                scenario_id=scenario_id,
                variant_id=variant_id,
                tx_scheme=tx_scheme,
                delays=delays,
                snr_db=snr_db,
                tb=tb,
                qm=qm,
            )
        if bool(getattr(cfg.simulation, "common_random_numbers", True)):
            seed = self._stable_seed(cfg.simulation.seed, scenario_id, snr_db)
        else:
            seed = self._stable_seed(cfg.simulation.seed, scenario_id, variant_id, snr_db)
        noise_var = float(10.0 ** (-float(snr_db) / 10.0))
        noise_var_ls = noise_var / float(len(cfg.resource.dmrs_symbol_indices))
        max_trials = int(max(cfg.simulation.n_trials_per_snr, cfg.simulation.max_trials_per_snr))
        target_trials = int(cfg.simulation.n_trials_per_snr)
        min_errors = int(cfg.simulation.min_block_errors)

        tb_errors = 0
        cb_errors = 0
        goodput_bits = 0
        ce_nmse_eff_values = []
        ce_nmse_branch_values = []
        cond_values = []
        rank_values = []
        trials = 0

        data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
        pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)

        while trials < max_trials:
            trials += 1
            rng = np.random.default_rng(self._stable_seed(seed, trials))
            channel = generate_tdl_channel(
                rng,
                grid,
                cfg.channel,
                n_tx=int(cfg.antenna.n_tx),
                n_rx=int(cfg.antenna.n_rx),
            )
            precoder = build_precoder(grid, cfg.resource, cfg.transmission, n_tx=int(cfg.antenna.n_tx))
            true_g = equivalent_channel(channel.H, precoder.C)

            payload = [
                rng.integers(0, 2, size=int(k), dtype=np.int8)
                for k in tb.cb_k_values
            ]
            coded = adapter.encode(payload)
            coded_cw = np.concatenate(coded)
            symbols = qam_modulate(coded_cw, qm)
            if len(symbols) != int(grid.n_data_re):
                raise RuntimeError("QAM symbol count does not match data RE count.")

            ls_noise = math.sqrt(noise_var_ls / 2.0) * (
                rng.normal(size=(int(cfg.antenna.n_rx), grid.pilot_count))
                + 1j * rng.normal(size=(int(cfg.antenna.n_rx), grid.pilot_count))
            )
            ls_obs = true_g[:, pilot_local] + ls_noise

            est = estimate_channel(
                method=str(cfg.channel_estimation.ce_method),
                tx_scheme=tx_scheme,
                ls_obs=ls_obs,
                true_g=true_g,
                true_H=channel.H,
                grid=grid,
                resource=cfg.resource,
                ce_cfg=cfg.channel_estimation,
                pdp=channel.pdp,
                delays=delays,
                noise_var_ls=noise_var_ls,
            )
            ce_nmse_eff_values.append(float(est.ce_nmse_eff))
            if math.isfinite(float(est.ce_nmse_branch)):
                ce_nmse_branch_values.append(float(est.ce_nmse_branch))
            if math.isfinite(float(est.cond_number)):
                cond_values.append(float(est.cond_number))
            if math.isfinite(float(est.effective_rank)):
                rank_values.append(float(est.effective_rank))

            y = self._apply_channel(true_g[:, data_local], symbols, noise_var, rng)
            z, no_eff = self._equalize_mrc(y, est.g_hat[:, data_local], noise_var)
            llr_cw = qam_demapper_maxlog(z, no_eff, qm)
            llrs_by_cb = self._split_llrs(llr_cw, tb.cb_e_values)
            dec = adapter.decode(llrs_by_cb, payload)
            if not dec.tb_success:
                tb_errors += 1
            cb_errors += sum(1 for ok in dec.cb_success if not ok)
            goodput_bits += int(dec.goodput_bits)

            if bool(cfg.simulation.save_trial_metrics):
                self.trial_rows.append({
                    "scenario_id": scenario_id,
                    "variant_id": variant_id,
                    "snr_db": float(snr_db),
                    "trial": int(trials),
                    "seed": int(seed),
                    "trial_seed": int(self._stable_seed(seed, trials)),
                    "tb_error": int(not dec.tb_success),
                    "cb_errors": int(sum(1 for ok in dec.cb_success if not ok)),
                    "ce_nmse_eff": float(est.ce_nmse_eff),
                    "ce_nmse_branch": float(est.ce_nmse_branch),
                    "cond_number": float(est.cond_number),
                })

            if trials >= target_trials and (min_errors <= 0 or tb_errors >= min_errors):
                break

        n_cb_trials = int(trials * len(tb.cb_k_values))
        row = self._base_metadata(cfg, grid, scenario_id, variant_id, tx_scheme, delays)
        row.update({
            "snr_db": float(snr_db),
            "noise_var": float(noise_var),
            "n_trials": int(trials),
            "tb_errors": int(tb_errors),
            "cb_errors": int(cb_errors),
            "bler": float(tb_errors) / float(trials),
            "cb_bler": float(cb_errors) / float(n_cb_trials),
            "ce_nmse_eff": self._mean(ce_nmse_eff_values),
            "ce_nmse_branch": self._mean(ce_nmse_branch_values),
            "cond_number": self._mean(cond_values),
            "effective_rank": self._mean(rank_values),
            "tbs_bits": int(tb.tb_size),
            "coded_bits": int(tb.coded_bits),
            "n_cbs": int(len(tb.cb_k_values)),
            "goodput_bits_per_slot": float(goodput_bits) / float(trials),
            "goodput_se_per_re": float(goodput_bits) / float(max(trials * grid.n_data_re, 1)),
            "common_random_numbers": bool(getattr(cfg.simulation, "common_random_numbers", True)),
            "base_seed": int(seed),
        })
        return row

    def _fixed_cdl_context(self, cfg: PlatformConfig, grid) -> FixedCDLStatisticsChannel:
        key = repr((
            dataclass_to_dict(cfg.channel),
            dataclass_to_dict(cfg.fixed_cdl_statistics),
            dataclass_to_dict(cfg.antenna),
            dataclass_to_dict(cfg.resource),
        ))
        if key not in self._fixed_cdl_contexts:
            self._fixed_cdl_contexts[key] = FixedCDLStatisticsChannel(cfg, grid)
        return self._fixed_cdl_contexts[key]

    def _run_snr_fixed_cdl(
        self,
        cfg: PlatformConfig,
        grid,
        adapter: SionnaLDPCAdapter | None,
        scenario_id: str,
        variant_id: str,
        tx_scheme: str,
        delays: List[int],
        snr_db: float,
        tb,
        qm: int,
    ) -> Dict[str, object]:
        context = self._fixed_cdl_context(cfg, grid)
        seed = self._stable_seed(cfg.simulation.seed, scenario_id, snr_db)
        noise_var = float(context.reference_receive_power * 10.0 ** (-float(snr_db) / 10.0))
        target_trials = int(cfg.simulation.n_trials_per_snr)
        max_trials = int(max(target_trials, cfg.simulation.max_trials_per_snr))
        min_errors = int(cfg.simulation.min_block_errors)
        codebook_type = str(cfg.fixed_cdl_statistics.codebook_type).lower()
        beam8_mode = codebook_type in {
            "dft_2x8_same_pol", "wide_beam_split", "angular_full_coverage_ultrawide"
        }
        aged_mrt = tx_scheme == "PLAN039_AGED_MRT_PRG6"
        if aged_mrt:
            precoder = None
        elif codebook_type in {"wide_beam_split", "angular_full_coverage_ultrawide"}:
            precoder = build_plan039_precoder(
                grid,
                tx_scheme,
                context.parent_codebook,
                context.plan039_arrays.get("sidon_selected_delay_indices"),
            )
        elif beam8_mode:
            precoder = build_beam8_precoder(
                grid,
                tx_scheme,
                context.selected_beam_indices,
                context.parent_codebook,
                int(cfg.resource.prg_size_rb),
            )
        else:
            precoder = build_fixed_cdl_precoder(
                grid,
                cfg.fixed_cdl_statistics,
                tx_scheme,
                context.selected_ssb,
                context.parent_codebook,
                context.secondary_codebook,
                int(cfg.resource.prg_size_rb),
            )
        data_sc_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
        pilot_sc_local = local_indices_for_subcarriers(grid, grid.pilot_subcarrier_indices)
        method = str(cfg.channel_estimation.ce_method).upper()
        ce_only = bool(getattr(cfg.simulation, "ce_only", False))
        estimator = None
        if method != "IDEAL":
            reference = context.frozen_beam_manifest.get("reference_pdp", {})
            delays_s = np.asarray(reference.get("delays_s", []), dtype=np.float64)
            powers = np.asarray(reference.get("powers", []), dtype=np.float64)
            time_real = np.asarray(reference.get("time_covariance_real", []), dtype=np.float64)
            time_imag = np.asarray(reference.get("time_covariance_imag", []), dtype=np.float64)
            if delays_s.size == 0 or powers.size == 0 or time_real.shape != (grid.n_symbols, grid.n_symbols):
                raise ValueError("Beam8 estimated CSI requires reference PDP and time covariance in the frozen manifest.")
            rf = frequency_covariance(delays_s, powers, grid.n_sc, float(grid.scs_khz) * 1e3)
            if codebook_type in {"wide_beam_split", "angular_full_coverage_ultrawide"}:
                frequencies = np.arange(grid.n_sc, dtype=float) * float(grid.scs_khz) * 1e3
                delay_indices = [] if aged_mrt else precoder.metadata["delay_grid_indices"]
                alpha = np.ones(grid.n_sc) if aged_mrt else precoder.metadata["alpha"]
                rf = plan039_frequency_covariance(
                    method, frequencies, context.plan039_arrays, delay_indices, alpha,
                )
            elif method == "BEAM8_CDD_AWARE_LMMSE":
                delta = np.arange(grid.n_sc)[:, None] - np.arange(grid.n_sc)[None, :]
                factor = np.mean(np.exp(-1j * 2.0 * np.pi * delta[:, :, None]
                                     * np.asarray(delays)[None, None, :] / float(grid.n_sc)), axis=2)
                rf = rf * factor
            assumed = TDLTimeFrequencyCovariance(time_real + 1j * time_imag, rf,
                "beam8_reference_cdd_aware" if method == "BEAM8_CDD_AWARE_LMMSE" else "beam8_reference_prg")
            normalized_noise = noise_var if codebook_type in {
                "wide_beam_split", "angular_full_coverage_ultrawide"
            } else (
                noise_var / float(context.reference_receive_power)
            )
            if method in {"BEAM8_PRG_LMMSE", "PLAN039_PRG_COMMON_REFERENCE_PDP"}:
                estimator = build_prg_time_frequency_rmmse_filter(
                    grid, assumed, int(cfg.resource.prg_size_rb) * 12, normalized_noise,
                    float(cfg.channel_estimation.diagonal_loading))
            else:
                estimator = build_time_frequency_rmmse_filter(
                    grid, assumed, normalized_noise, float(cfg.channel_estimation.diagonal_loading))
        tb_errors = 0
        cb_errors = 0
        goodput_bits = 0
        ce_nmse_values: List[float] = []
        precoder_power_min = float("inf")
        precoder_power_max = 0.0
        aged_replay_error_max = 0.0
        trials = 0
        while trials < max_trials:
            trials += 1
            absolute_trial = int(cfg.simulation.absolute_trial_start) + trials
            payload_seed = self._stable_seed(seed, absolute_trial, "payload")
            data_noise_seed = self._stable_seed(seed, absolute_trial, "data_noise")
            pilot_noise_seed = self._stable_seed(seed, absolute_trial, "pilot_noise")
            payload_rng = np.random.default_rng(payload_seed)
            data_rng = np.random.default_rng(data_noise_seed)
            pilot_rng = np.random.default_rng(pilot_noise_seed)
            if aged_mrt:
                age_s = float(cfg.transmission.aged_csi_ms) * 1e-3
                channel, old_channel, replay_error = context.generate_with_aged_csi(
                    absolute_trial - 1, age_s
                )
                weights = build_aged_mrt_prg_precoder(
                    old_channel,
                    prg_size_rb=int(cfg.resource.prg_size_rb),
                    normalize=True,
                )[0]
                power = np.sum(np.abs(weights) ** 2, axis=1)
                precoder_power_min = min(precoder_power_min, float(np.min(power)))
                precoder_power_max = max(precoder_power_max, float(np.max(power)))
                aged_replay_error_max = max(aged_replay_error_max, float(replay_error))
                true_g = equivalent_channel(channel.H, weights)[0]
            else:
                channel = context.generate(absolute_trial - 1)
                assert precoder is not None
                weights = precoder.C
                power = np.sum(np.abs(weights) ** 2, axis=1)
                precoder_power_min = min(precoder_power_min, float(np.min(power)))
                precoder_power_max = max(precoder_power_max, float(np.max(power)))
                true_g = equivalent_channel(channel.H, weights)[0]
            true_pilot = true_g[:, grid.pilot_symbol_indices, pilot_sc_local]
            true_data = true_g[:, grid.data_symbol_indices, data_sc_local]
            payload = None
            symbols = None
            if not ce_only:
                if adapter is None:
                    raise RuntimeError("LDPC adapter is unavailable for a BLER run.")
                payload = [payload_rng.integers(0, 2, size=int(k), dtype=np.int8) for k in tb.cb_k_values]
                coded_cw = np.concatenate(adapter.encode(payload))
                symbols = qam_modulate(coded_cw, qm)
                if len(symbols) != int(grid.n_data_re):
                    raise RuntimeError("QAM symbol count does not match data RE count.")
            if method == "IDEAL":
                g_hat_data = true_data
                ce_nmse = 0.0
            else:
                ls_obs = construct_ls_observations(true_pilot, noise_var, pilot_rng)
                if estimator is None:
                    raise RuntimeError("Beam8 LMMSE estimator was not constructed.")
                g_hat_data = estimator.estimate_data(ls_obs)
                ce_nmse = float(np.sum(np.abs(g_hat_data - true_data) ** 2) /
                                max(float(np.sum(np.abs(true_data) ** 2)), 1e-30))
            ce_nmse_values.append(ce_nmse)
            dec = None
            if not ce_only:
                assert symbols is not None and payload is not None and adapter is not None
                y = self._apply_channel(true_data, symbols, noise_var, data_rng)
                z, no_eff = self._equalize_mrc(y, g_hat_data, noise_var)
                llr_cw = qam_demapper_maxlog(z, no_eff, qm)
                dec = adapter.decode(self._split_llrs(llr_cw, tb.cb_e_values), payload)
                if not dec.tb_success:
                    tb_errors += 1
                cb_errors += sum(1 for ok in dec.cb_success if not ok)
                goodput_bits += int(dec.goodput_bits)
            if bool(cfg.simulation.save_trial_metrics):
                self.trial_rows.append({
                    "scenario_id": scenario_id,
                    "variant_id": variant_id,
                    "snr_db": float(snr_db),
                    "trial": int(absolute_trial),
                    "tb_error": -1 if dec is None else int(not dec.tb_success),
                    "cb_errors": -1 if dec is None else int(sum(1 for ok in dec.cb_success if not ok)),
                    "selected_ssb": int(context.selected_ssb),
                    "selected_beam_indices": ",".join(str(value) for value in getattr(context, "selected_beam_indices", [])),
                    "frozen_beam_manifest_sha256": str(getattr(cfg.fixed_cdl_statistics, "frozen_beam_manifest_sha256", "")),
                    "reference_receive_power": float(context.reference_receive_power),
                    "noise_var": float(noise_var),
                    "channel_realization_index": int(absolute_trial - 1),
                    "channel_realization_seed_root": int(cfg.fixed_cdl_statistics.realization_seed),
                    "payload_noise_seed": int(payload_seed),
                    "payload_seed": int(payload_seed),
                    "data_noise_seed": int(data_noise_seed),
                    "pilot_noise_seed": int(pilot_noise_seed),
                    "receiver_mode": "ce_only" if ce_only else ("ideal" if method == "IDEAL" else "estimated"),
                    "ce_method": method,
                    "ce_nmse_eff": ce_nmse,
                    "ce_error_energy": float(np.sum(np.abs(g_hat_data - true_data) ** 2)),
                    "ce_true_energy": float(np.sum(np.abs(true_data) ** 2)),
                    "filter_condition_number": float(estimator.condition_number) if estimator else float("nan"),
                    "filter_numerical_loading": float(estimator.numerical_jitter) if estimator else 0.0,
                    "filter_minimum_singular_value": float(estimator.minimum_singular_value) if estimator else float("nan"),
                    "precoder_power_min": float(np.min(power)),
                    "precoder_power_max": float(np.max(power)),
                    "precoder_type": "aged_mrt_prg6" if aged_mrt else "static",
                    "precoder_csi_age_ms": float(cfg.transmission.aged_csi_ms),
                    "aged_csi_replay_error": float(replay_error) if aged_mrt else 0.0,
                })
            if trials >= target_trials and (min_errors <= 0 or tb_errors >= min_errors):
                break

        row = self._base_metadata(cfg, grid, scenario_id, variant_id, tx_scheme, delays)
        n_cb_trials = int(trials * len(tb.cb_k_values))
        row.update({
            "snr_db": float(snr_db),
            "noise_var": float(noise_var),
            "noise_var_ls": float(noise_var),
            "n_trials": int(trials),
            "absolute_trial_start": int(cfg.simulation.absolute_trial_start),
            "absolute_trial_stop": int(cfg.simulation.absolute_trial_start) + int(trials),
            "tb_errors": int(tb_errors),
            "cb_errors": int(cb_errors),
            "bler": float("nan") if ce_only else float(tb_errors) / float(trials),
            "cb_bler": float("nan") if ce_only else float(cb_errors) / float(n_cb_trials),
            "ce_nmse_eff": self._mean(ce_nmse_values),
            "ce_nmse_sum": float(np.sum(ce_nmse_values)),
            "ce_nmse_sum_squares": float(np.sum(np.square(ce_nmse_values))),
            "ce_nmse_branch": float("nan"),
            "cond_number": float(estimator.condition_number) if estimator else float("nan"),
            "numerical_jitter": float(estimator.numerical_jitter) if estimator else 0.0,
            "minimum_singular_value": float(estimator.minimum_singular_value) if estimator else float("nan"),
            "effective_rank": float("nan"),
            "covariance_type": estimator.covariance_type if estimator else "ideal_csi",
            "receiver_mode": "ce_only" if ce_only else ("ideal" if method == "IDEAL" else "estimated"),
            "ce_method": method,
            "true_covariance_type": "fixed_cdl_statistics_monte_carlo_rt",
            "tbs_bits": int(tb.tb_size),
            "coded_bits": int(tb.coded_bits),
            "n_cbs": int(len(tb.cb_k_values)),
            "goodput_bits_per_slot": float(goodput_bits) / float(trials),
            "goodput_se_per_re": float(goodput_bits) / float(max(trials * grid.n_data_re, 1)),
            "common_random_numbers": True,
            "base_seed": int(seed),
            "selected_ssb": int(context.selected_ssb),
            "selected_beam_indices": ",".join(str(value) for value in getattr(context, "selected_beam_indices", [])),
            "frozen_beam_manifest_sha256": str(getattr(cfg.fixed_cdl_statistics, "frozen_beam_manifest_sha256", "")),
            "reference_receive_power": float(context.reference_receive_power),
            "statistics_seed": int(cfg.fixed_cdl_statistics.statistics_seed),
            "realization_seed": int(cfg.fixed_cdl_statistics.realization_seed),
            "covariance_realizations": int(cfg.fixed_cdl_statistics.covariance_realizations),
            "mean_aod_deg": float(cfg.fixed_cdl_statistics.mean_aod_deg),
            "precoder_raw_power_min": float(
                precoder_power_min if aged_mrt else precoder.metadata["raw_power_min"]
            ),
            "precoder_raw_power_max": float(
                precoder_power_max if aged_mrt else precoder.metadata["raw_power_max"]
            ),
            "precoder_type": "aged_mrt_prg6" if aged_mrt else "static",
            "precoder_csi_age_ms": float(cfg.transmission.aged_csi_ms),
            "aged_csi_replay_error_max": float(aged_replay_error_max),
        })
        return row

    def _run_snr_sionna_tdl(
        self,
        cfg: PlatformConfig,
        grid,
        adapter: SionnaLDPCAdapter,
        scenario_id: str,
        variant_id: str,
        tx_scheme: str,
        delays: List[int],
        snr_db: float,
        tb,
        qm: int,
    ) -> Dict[str, object]:
        if bool(getattr(cfg.simulation, "common_random_numbers", True)):
            seed = self._stable_seed(cfg.simulation.seed, scenario_id, snr_db)
        else:
            seed = self._stable_seed(cfg.simulation.seed, scenario_id, variant_id, snr_db)
        noise_var = float(10.0 ** (-float(snr_db) / 10.0))
        noise_var_ls = noise_var
        max_trials = int(max(cfg.simulation.n_trials_per_snr, cfg.simulation.max_trials_per_snr))
        target_trials = int(cfg.simulation.n_trials_per_snr)
        min_errors = int(cfg.simulation.min_block_errors)

        method = str(cfg.channel_estimation.ce_method).upper()
        backend = str(cfg.channel.backend).lower()
        if backend == "sionna_cdl":
            true_covariance = cdl_spatial_unaware_covariance(grid, cfg.channel, delays)
            unknown_covariance = cdl_spatial_unaware_covariance(grid, cfg.channel)
        else:
            true_covariance = tdl_known_delay_covariance(grid, cfg.channel, delays)
            unknown_covariance = tdl_unknown_delay_covariance(grid, cfg.channel)
        estimator = None
        if method in ("TF_RMMSE_UNKNOWN", "RMMSE_TF_UNKNOWN", "RMMSE_2D_UNKNOWN"):
            assumed_covariance = unknown_covariance
        elif method in ("TF_RMMSE_KNOWN", "RMMSE_TF_KNOWN", "RMMSE_2D_KNOWN"):
            assumed_covariance = true_covariance
        elif method == "IDEAL":
            assumed_covariance = None
        else:
            raise ValueError(
                "Sionna TDL/CDL backends require IDEAL, TF_RMMSE_KNOWN, or TF_RMMSE_UNKNOWN "
                "channel estimation."
            )
        if assumed_covariance is not None:
            estimator = build_time_frequency_rmmse_filter(
                grid,
                assumed_covariance,
                noise_variance=noise_var_ls,
                diagonal_loading=float(cfg.channel_estimation.diagonal_loading),
            )
        precoder = build_precoder(grid, cfg.resource, cfg.transmission, n_tx=int(cfg.antenna.n_tx))
        pilot_sc_local = local_indices_for_subcarriers(grid, grid.pilot_subcarrier_indices)
        data_sc_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)

        tb_errors = 0
        cb_errors = 0
        goodput_bits = 0
        ce_nmse_values: List[float] = []
        trials = 0
        while trials < max_trials:
            trials += 1
            rng = np.random.default_rng(self._stable_seed(seed, trials))
            channel_seed = self._stable_seed(seed, trials, backend)
            channel = generate_sionna_channel(
                grid=grid,
                channel=cfg.channel,
                n_tx=int(cfg.antenna.n_tx),
                n_rx=int(cfg.antenna.n_rx),
                batch_size=1,
                seed=channel_seed,
            )
            true_g = equivalent_channel(channel.H, precoder.C)[0]
            true_pilot = true_g[:, grid.pilot_symbol_indices, pilot_sc_local]
            true_data = true_g[:, grid.data_symbol_indices, data_sc_local]

            payload = [
                rng.integers(0, 2, size=int(k), dtype=np.int8)
                for k in tb.cb_k_values
            ]
            coded_cw = np.concatenate(adapter.encode(payload))
            symbols = qam_modulate(coded_cw, qm)
            if len(symbols) != int(grid.n_data_re):
                raise RuntimeError("QAM symbol count does not match data RE count.")

            if method == "IDEAL":
                g_hat_data = true_data
            else:
                ls_obs = construct_ls_observations(true_pilot, noise_var_ls, rng)
                if estimator is None:
                    raise RuntimeError("RMMSE estimator was not constructed.")
                g_hat_data = estimator.estimate_data(ls_obs)
            ce_nmse = float(
                np.sum(np.abs(g_hat_data - true_data) ** 2)
                / max(float(np.sum(np.abs(true_data) ** 2)), 1e-30)
            )
            ce_nmse_values.append(ce_nmse)

            y = self._apply_channel(true_data, symbols, noise_var, rng)
            z, no_eff = self._equalize_mrc(y, g_hat_data, noise_var)
            llr_cw = qam_demapper_maxlog(z, no_eff, qm)
            llrs_by_cb = self._split_llrs(llr_cw, tb.cb_e_values)
            dec = adapter.decode(llrs_by_cb, payload)
            if not dec.tb_success:
                tb_errors += 1
            cb_errors += sum(1 for ok in dec.cb_success if not ok)
            goodput_bits += int(dec.goodput_bits)

            if bool(cfg.simulation.save_trial_metrics):
                self.trial_rows.append({
                    "scenario_id": scenario_id,
                    "variant_id": variant_id,
                    "snr_db": float(snr_db),
                    "trial": int(trials),
                    "seed": int(seed),
                    "channel_seed": int(channel_seed),
                    "tb_error": int(not dec.tb_success),
                    "cb_errors": int(sum(1 for ok in dec.cb_success if not ok)),
                    "ce_nmse_eff": ce_nmse,
                    "cond_number": float(estimator.condition_number) if estimator else float("nan"),
                    "numerical_jitter": float(estimator.numerical_jitter) if estimator else 0.0,
                    "covariance_type": estimator.covariance_type if estimator else "ideal_csi",
                })

            if trials >= target_trials and (min_errors <= 0 or tb_errors >= min_errors):
                break

        n_cb_trials = int(trials * len(tb.cb_k_values))
        row = self._base_metadata(cfg, grid, scenario_id, variant_id, tx_scheme, delays)
        row.update({
            "snr_db": float(snr_db),
            "noise_var": float(noise_var),
            "noise_var_ls": float(noise_var_ls),
            "n_trials": int(trials),
            "tb_errors": int(tb_errors),
            "cb_errors": int(cb_errors),
            "bler": float(tb_errors) / float(trials),
            "cb_bler": float(cb_errors) / float(n_cb_trials),
            "ce_nmse_eff": self._mean(ce_nmse_values),
            "ce_nmse_branch": float("nan"),
            "cond_number": float(estimator.condition_number) if estimator else float("nan"),
            "effective_rank": float("nan"),
            "numerical_jitter": float(estimator.numerical_jitter) if estimator else 0.0,
            "covariance_type": estimator.covariance_type if estimator else "ideal_csi",
            "true_covariance_type": true_covariance.covariance_type,
            "tbs_bits": int(tb.tb_size),
            "coded_bits": int(tb.coded_bits),
            "n_cbs": int(len(tb.cb_k_values)),
            "goodput_bits_per_slot": float(goodput_bits) / float(trials),
            "goodput_se_per_re": float(goodput_bits) / float(max(trials * grid.n_data_re, 1)),
            "common_random_numbers": bool(getattr(cfg.simulation, "common_random_numbers", True)),
            "base_seed": int(seed),
        })
        return row

    @staticmethod
    def _apply_channel(g_by_re: np.ndarray, symbols: np.ndarray, noise_var: float, rng: np.random.Generator) -> np.ndarray:
        g = np.asarray(g_by_re, dtype=np.complex128)
        x = np.asarray(symbols, dtype=np.complex128).reshape(-1)
        noise = math.sqrt(float(noise_var) / 2.0) * (
            rng.normal(size=g.shape) + 1j * rng.normal(size=g.shape)
        )
        return g * x[None, :] + noise

    @staticmethod
    def _equalize_mrc(y: np.ndarray, g_hat: np.ndarray, noise_var: float) -> tuple[np.ndarray, np.ndarray]:
        gh = np.asarray(g_hat, dtype=np.complex128)
        yy = np.asarray(y, dtype=np.complex128)
        denom = np.sum(np.abs(gh) ** 2, axis=0)
        denom = np.maximum(denom, 1e-10)
        z = np.sum(np.conj(gh) * yy, axis=0) / denom
        no_eff = float(noise_var) / denom
        return z, no_eff

    @staticmethod
    def _split_llrs(llr: np.ndarray, cb_e_values: List[int]) -> List[np.ndarray]:
        out = []
        cursor = 0
        arr = np.asarray(llr, dtype=np.float64).reshape(-1)
        for e in cb_e_values:
            e = int(e)
            out.append(arr[cursor:cursor + e].copy())
            cursor += e
        if cursor != arr.size:
            raise ValueError(f"CB E values sum to {cursor}, but LLR stream has {arr.size}.")
        return out

    @staticmethod
    def _stable_seed(*items: object) -> int:
        text = "|".join(str(x) for x in items)
        acc = 2166136261
        for ch in text:
            acc ^= ord(ch)
            acc = (acc * 16777619) % (2 ** 32)
        return int(acc)

    @staticmethod
    def _mean(values: List[float]) -> float:
        vals = [float(x) for x in values if math.isfinite(float(x))]
        return float(np.mean(vals)) if vals else float("nan")

    def _base_metadata(
        self,
        cfg: PlatformConfig,
        grid,
        scenario_id: str,
        variant_id: str,
        tx_scheme: str,
        delays: List[int],
    ) -> Dict[str, object]:
        fixed_mode = str(getattr(cfg.channel, "backend", "legacy_exponential")).lower() == "fixed_cdl_statistics"
        delay_denominator = float(grid.n_sc) if fixed_mode else float(grid.n_fft)
        return {
            "scenario_id": scenario_id,
            "variant_id": variant_id,
            "channel_model": str(cfg.channel.model),
            "channel_backend": str(getattr(cfg.channel, "backend", "legacy_exponential")),
            "tdl_profile": "" if fixed_mode else str(getattr(cfg.channel, "tdl_profile", "")),
            "cdl_profile": str(getattr(cfg.channel, "cdl_profile", "")),
            "delay_spread_ns": float(cfg.channel.delay_spread_ns),
            "carrier_frequency_hz": float(getattr(cfg.channel, "carrier_frequency_hz", float("nan"))),
            "speed_kmh": float(getattr(cfg.channel, "ue_speed_kmh", 0.0)),
            "speed_mps": float(getattr(cfg.channel, "ue_speed_kmh", 0.0)) / 3.6,
            "n_tx": int(cfg.antenna.n_tx),
            "n_rx": int(cfg.antenna.n_rx),
            "pdsch_rb": int(cfg.resource.n_prbs),
            "pdsch_symbols": int(cfg.resource.pdsch_n_symbols),
            "dmrs_symbols": int(len(cfg.resource.dmrs_symbol_indices)),
            "dmrs_symbol_indices": ",".join(str(int(x)) for x in cfg.resource.dmrs_symbol_indices),
            "dmrs_spacing_sc": int(cfg.resource.dmrs_spacing_sc),
            "dmrs_overhead": float(grid.dmrs_overhead),
            "n_dmrs_re": int(grid.n_dmrs_re),
            "ofdm_symbol_duration_s": float(grid.ofdm_symbol_duration_s),
            "data_re": int(grid.n_data_re),
            "tx_scheme": tx_scheme,
            "ce_method": str(cfg.channel_estimation.ce_method),
            "cdd_delay_vector": ",".join(str(int(x)) for x in delays),
            "cdd_delay_seconds": ",".join(
                str(float(x) / (delay_denominator * float(grid.scs_khz) * 1e3)) for x in delays
            ),
            "cdd_delay_definition": "active_band_dft_grid" if fixed_mode else "fft_sample",
            "prg_size_rb": int(cfg.resource.prg_size_rb),
            "prg_codebook": str(cfg.transmission.prg_codebook),
            "mcs_table": str(cfg.mcs.table),
            "mcs_index": int(cfg.mcs.index),
        }

    def _target_rows(self) -> List[Dict[str, object]]:
        groups: Dict[tuple[str, str], List[Dict[str, object]]] = {}
        for row in self.summary_rows:
            key = (str(row["scenario_id"]), str(row["variant_id"]))
            groups.setdefault(key, []).append(row)
        out = []
        target = float(self.cfg.simulation.bler_target)
        baseline_by_scenario = {}
        reference_by_scenario = {}
        for key, rows in groups.items():
            snr = interpolate_target_snr(rows, target_bler=target)
            first = rows[0]
            item = {
                "scenario_id": key[0],
                "variant_id": key[1],
                "target_bler": target,
                "snr_at_target_bler_db": snr,
                "tx_scheme": first.get("tx_scheme"),
                "ce_method": first.get("ce_method"),
                "pdsch_rb": first.get("pdsch_rb"),
                "delay_spread_ns": first.get("delay_spread_ns"),
                "dmrs_spacing_sc": first.get("dmrs_spacing_sc"),
                "cdd_delay_vector": first.get("cdd_delay_vector"),
            }
            out.append(item)
            if str(first.get("tx_scheme", "")).startswith("PRG"):
                baseline_by_scenario.setdefault(key[0], snr)
                reference_by_scenario.setdefault(key[0], snr)
            if str(first.get("tx_scheme", "")).upper() == "BASELINE":
                reference_by_scenario.setdefault(key[0], snr)
        for item in out:
            base = baseline_by_scenario.get(str(item["scenario_id"]), float("nan"))
            reference = reference_by_scenario.get(str(item["scenario_id"]), float("nan"))
            snr = float(item["snr_at_target_bler_db"])
            item["snr_gain_vs_prg_db"] = float(base - snr) if math.isfinite(base) and math.isfinite(snr) else float("nan")
            item["snr_gain_vs_baseline_db"] = (
                float(reference - snr) if math.isfinite(reference) and math.isfinite(snr) else float("nan")
            )
        return out
