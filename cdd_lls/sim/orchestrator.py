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
from cdd_lls.phy.channel_tdl import generate_tdl_channel
from cdd_lls.phy.estimators import estimate_channel
from cdd_lls.phy.ldpc import SionnaLDPCAdapter
from cdd_lls.phy.precoding import build_precoder, equivalent_channel, normalize_delay_vector
from cdd_lls.phy.qam import qam_demapper_maxlog, qam_modulate
from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
from cdd_lls.sim.stats import interpolate_target_snr, save_csv, save_json, snr_values
from cdd_lls.utils.plotting import plot_bler, plot_nmse


class CDDLinkLevelOrchestrator:
    def __init__(self, config: PlatformConfig):
        self.cfg = config
        self.summary_rows: List[Dict[str, object]] = []
        self.trial_rows: List[Dict[str, object]] = []

    def run(self) -> Path:
        out = self._make_output_dir()
        save_resolved_config(self.cfg, out / "resolved_config.yaml")
        run_specs = self._expanded_run_specs()
        for spec in run_specs:
            self._run_spec(spec)
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
        stamp = dt.datetime.now().strftime("sim_%Y%m%d_%H%M%S")
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
        adapter = SionnaLDPCAdapter(
            tb.cb_k_values,
            tb.cb_e_values,
            num_iter=int(cfg.receiver.max_ldpc_iterations),
            llr_clip=float(cfg.receiver.llr_clip),
        )

        tx_scheme = str(cfg.transmission.tx_scheme).upper()
        delays = normalize_delay_vector(
            cfg.transmission.cdd_delay_vector,
            int(cfg.antenna.n_tx),
            int(cfg.transmission.cdd_base_delay),
        )

        for snr_db in snr_values(cfg.simulation.snr_range_db):
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
        return {
            "scenario_id": scenario_id,
            "variant_id": variant_id,
            "channel_model": str(cfg.channel.model),
            "delay_spread_ns": float(cfg.channel.delay_spread_ns),
            "speed_kmh": 0.0,
            "n_tx": int(cfg.antenna.n_tx),
            "n_rx": int(cfg.antenna.n_rx),
            "pdsch_rb": int(cfg.resource.n_prbs),
            "pdsch_symbols": int(cfg.resource.pdsch_n_symbols),
            "dmrs_symbols": int(len(cfg.resource.dmrs_symbol_indices)),
            "dmrs_spacing_sc": int(cfg.resource.dmrs_spacing_sc),
            "dmrs_overhead": float(grid.dmrs_overhead),
            "n_dmrs_re": int(grid.n_dmrs_re),
            "data_re": int(grid.n_data_re),
            "tx_scheme": tx_scheme,
            "ce_method": str(cfg.channel_estimation.ce_method),
            "cdd_delay_vector": ",".join(str(int(x)) for x in delays),
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
        for item in out:
            base = baseline_by_scenario.get(str(item["scenario_id"]), float("nan"))
            snr = float(item["snr_at_target_bler_db"])
            item["snr_gain_vs_prg_db"] = float(base - snr) if math.isfinite(base) and math.isfinite(snr) else float("nan")
        return out
