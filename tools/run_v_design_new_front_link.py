from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime
from pathlib import Path
import sys
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ChannelConfig, ResourceConfig
from cdd_lls.core.mcs import build_tb_layout, get_mcs
from cdd_lls.phy.channel_tdl import generate_tdl_channel, make_exponential_pdp
from cdd_lls.phy.estimators import covariance_matrix, nmse
from cdd_lls.phy.ldpc import SionnaLDPCAdapter
from cdd_lls.phy.precoding import equivalent_channel
from cdd_lls.phy.qam import qam_demapper_maxlog, qam_modulate
from cdd_lls.phy.resource_grid import ResourceGrid, build_resource_grid, local_indices_for_subcarriers
from tools.run_v_design_pareto_link_curves import interp_target_snr, svg_all_pareto_plot
from tools.run_v_design_piecewise_tradeoff import (
    VDesign,
    cdd_design,
    decode_same_tb_batch,
    design_from_phase,
    equalize_mrc,
    make_alg1_estimator,
    mean_or_nan,
    matrix_metrics,
    covariance_metrics,
    phase_from_piecewise_delays,
    stable_seed,
    write_csv,
)


SCAN_DIR = Path(
    "outputs/v_design_balanced_slope_scan/"
    "v_design_balanced_slope_scan_20260619_162622"
)

CDD_REFERENCES = {
    "N3": ("C0", 1),
    "N6": ("C0", 1),
    "N8": ("C1", 2),
    "N23": ("C2", 3),
}

CDD_REFERENCES_36PRB = {
    "N3": ("C5", 6),
    "N6": ("C6", 7),
    "N8": ("C3", 4),
    "N23": ("C4", 5),
}

CDD_REFERENCES_24PRB = {
    "N3": ("C6", 7),
    "N6": ("C6", 7),
    "N8": ("C6", 7),
    "N23": ("C4", 5),
}

CDD_IDS_BY_STEP = {1: "C0", 2: "C1", 3: "C2", 4: "C3", 5: "C4", 6: "C5", 7: "C6", 64: "C7"}


def cdd_reference(n_prbs: int, pareto_id: str) -> Tuple[str, int]:
    references_by_bandwidth = {
        24: CDD_REFERENCES_24PRB,
        36: CDD_REFERENCES_36PRB,
        48: CDD_REFERENCES,
    }
    references = references_by_bandwidth.get(int(n_prbs), CDD_REFERENCES)
    return references[str(pareto_id)]


def parse_int_list(text: str) -> List[int]:
    return [int(item.strip()) for item in str(text).split(",") if item.strip()]


def parse_float_list(text: str) -> List[float]:
    return [float(item.strip()) for item in str(text).split(",") if item.strip()]


def parse_text_list(text: str) -> List[str]:
    return [item.strip() for item in str(text).split(",") if item.strip()]


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def build_grid(n_prbs: int, dmrs_spacing_sc: int) -> Tuple[ResourceConfig, ResourceGrid]:
    resource = ResourceConfig(
        carrier_bandwidth_mhz=100.0,
        scs_khz=30,
        n_fft=4096,
        n_prbs=int(n_prbs),
        pdsch_n_symbols=10,
        dmrs_symbol_indices=[2, 7],
        dmrs_spacing_sc=int(dmrs_spacing_sc),
        prg_size_rb=4,
    )
    grid = build_resource_grid(resource)
    if grid.n_sc % 8 != 0:
        raise ValueError("The active subcarrier count must be divisible by eight segments.")
    return resource, grid


def load_selected_designs(
    grid: ResourceGrid,
    scan_dir: Path,
    selected_ids: Sequence[str],
    additional_cdd_steps: Sequence[int] = (),
    include_reference_cdd: bool = True,
) -> List[Tuple[str, VDesign, Dict[str, object]]]:
    pareto_by_id = {
        str(row["pareto_id"]): row
        for row in read_csv(scan_dir / "new_signed_balanced_pareto.csv")
    }
    with (scan_dir / "new_signed_balanced_implementation_manifest.json").open("r", encoding="utf-8") as handle:
        manifest_rows = json.load(handle)
    manifest_by_name = {str(row["design_name"]): row for row in manifest_rows}

    selected: List[Tuple[str, VDesign, Dict[str, object]]] = []
    for pareto_id in selected_ids:
        prow = pareto_by_id[pareto_id]
        name = str(prow["precoder"])
        manifest = manifest_by_name[name]
        delays = np.asarray(json.loads(manifest["segment_tx_delay_samples"]), dtype=np.float64)
        phase = phase_from_piecewise_delays(grid, delays, continuous=True)
        design = design_from_phase(
            label=name,
            family="piecewise_linear_continuous",
            phase=phase,
            metadata={
                "pareto_id": pareto_id,
                "delay_alphabet_samples": manifest["delay_alphabet_samples"],
                "layout_id": manifest["layout_id"],
                "layout_seed": manifest["layout_seed"],
                "segment_tx_delay_samples": manifest["segment_tx_delay_samples"],
            },
        )
        selected.append((pareto_id, design, manifest))

    n_prbs = grid.n_sc // 12
    reference_steps = (
        {cdd_reference(n_prbs, pareto_id)[1] for pareto_id in selected_ids}
        if include_reference_cdd
        else set()
    )
    cdd_steps = sorted(
        reference_steps | {int(step) for step in additional_cdd_steps}
    )
    for step in cdd_steps:
        selected.append((CDD_IDS_BY_STEP[step], cdd_design(grid, 8, step), {}))
    return selected


def fixed_interleaver(coded_bits: int, n_prbs: int, base_seed: int) -> Tuple[np.ndarray, int]:
    seed = stable_seed("full-coded-bit-interleaver", n_prbs, coded_bits, base=int(base_seed))
    permutation = np.random.default_rng(seed).permutation(int(coded_bits)).astype(np.int64)
    return permutation, int(seed)


def apply_interleaver(bits: np.ndarray, permutation: np.ndarray, mode: str) -> np.ndarray:
    if mode == "full":
        return np.asarray(bits)[permutation]
    if mode == "none":
        return np.asarray(bits)
    raise ValueError(f"Unknown interleaver mode={mode}")


def undo_interleaver(llrs: np.ndarray, permutation: np.ndarray, mode: str) -> np.ndarray:
    if mode == "full":
        return np.asarray(llrs)[np.argsort(permutation)]
    if mode == "none":
        return np.asarray(llrs)
    raise ValueError(f"Unknown interleaver mode={mode}")


def curve_family(pareto_id: str) -> str:
    return "cdd_linear_phase" if str(pareto_id).startswith("C") else "piecewise_linear_continuous"


def run_ideal_bler_for_bandwidth(
    args: argparse.Namespace,
    n_prbs: int,
    resource: ResourceConfig,
    grid: ResourceGrid,
    selected: Sequence[Tuple[str, VDesign, Dict[str, object]]],
    output_dir: Path,
    figure_dir: Path,
) -> List[Dict[str, object]]:
    if int(args.ideal_trials) <= 0:
        return []
    modes = parse_text_list(args.interleaver_modes)
    mcs = get_mcs("nr_256qam", int(args.mcs_index), None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(
        tb.cb_k_values,
        tb.cb_e_values,
        num_iter=int(args.ldpc_iterations),
        llr_clip=float(args.llr_clip),
    )
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    channel = ChannelConfig(
        delay_spread_ns=float(args.delay_spread_ns),
        max_delay_factor=float(args.max_delay_factor),
    )
    permutation, interleaver_seed = fixed_interleaver(tb.coded_bits, n_prbs, int(args.seed))
    np.save(output_dir / f"coded_bit_interleaver_{int(n_prbs)}prb.npy", permutation)

    rows: List[Dict[str, object]] = []
    for snr_db in parse_float_list(args.ideal_snrs):
        noise_var = 10.0 ** (-float(snr_db) / 10.0)
        accum = {
            (mode, pareto_id): {"tb_errors": 0, "cb_errors": 0}
            for mode in modes
            for pareto_id, _, _ in selected
        }
        for trial in range(1, int(args.ideal_trials) + 1):
            rng = np.random.default_rng(
                stable_seed("new-front-ideal", n_prbs, trial, base=int(args.seed))
            )
            realization = generate_tdl_channel(rng, grid, channel, n_tx=8, n_rx=1)
            payload = [rng.integers(0, 2, size=int(k), dtype=np.int8) for k in tb.cb_k_values]
            coded_bits = np.concatenate(adapter.encode(payload))
            noise = (
                rng.normal(size=(1, len(data_local)))
                + 1j * rng.normal(size=(1, len(data_local)))
            ) * math.sqrt(noise_var / 2.0)

            llrs: List[np.ndarray] = []
            trial_keys: List[Tuple[str, str]] = []
            for mode in modes:
                mapped_bits = apply_interleaver(coded_bits, permutation, mode)
                symbols = qam_modulate(mapped_bits, int(mcs.qm))
                for pareto_id, design, _ in selected:
                    g = equivalent_channel(realization.H, design.C)
                    y = g[:, data_local] * symbols[None, :] + noise
                    z, no_eff = equalize_mrc(y, g[:, data_local], noise_var)
                    llr_mapped = qam_demapper_maxlog(z, no_eff, int(mcs.qm))
                    llrs.append(undo_interleaver(llr_mapped, permutation, mode))
                    trial_keys.append((mode, pareto_id))
            decoded = decode_same_tb_batch(adapter, llrs, payload)
            for key, result in zip(trial_keys, decoded):
                accum[key]["tb_errors"] += int(not result.tb_success)
                accum[key]["cb_errors"] += int(sum(1 for ok in result.cb_success if not ok))
            if int(args.progress_every) > 0 and trial % int(args.progress_every) == 0:
                print(
                    f"[ideal] PRB={n_prbs} SNR={snr_db:g} trial={trial}/{int(args.ideal_trials)}",
                    flush=True,
                )

        for mode in modes:
            for pareto_id, design, manifest in selected:
                result = accum[(mode, pareto_id)]
                rows.append({
                    "bandwidth_prb": int(n_prbs),
                    "active_subcarriers": int(grid.n_sc),
                    "occupied_bandwidth_mhz": float(grid.n_sc * resource.scs_khz / 1000.0),
                    "segment_length_sc": int(grid.n_sc // 8),
                    "segment_bandwidth_mhz": float((grid.n_sc // 8) * resource.scs_khz / 1000.0),
                    "interleaver_mode": mode,
                    "interleaver_seed": int(interleaver_seed) if mode == "full" else "",
                    "pareto_id": pareto_id,
                    "precoder": design.label,
                    "family": curve_family(pareto_id),
                    "snr_db": float(snr_db),
                    "trials": int(args.ideal_trials),
                    "tb_errors": int(result["tb_errors"]),
                    "cb_errors": int(result["cb_errors"]),
                    "bler": float(result["tb_errors"]) / float(args.ideal_trials),
                    "csi_mode": "ideal",
                    "mcs_index": int(args.mcs_index),
                    "qm": int(mcs.qm),
                    "code_rate": float(mcs.code_rate),
                    "tbs_bits": int(tb.tb_size),
                    "coded_bits": int(tb.coded_bits),
                    "data_re": int(grid.n_data_re),
                    "dmrs_spacing_sc": int(args.dmrs_spacing_sc),
                    "delay_spread_ns": float(args.delay_spread_ns),
                    "reference_cdd_id": cdd_reference(n_prbs, pareto_id)[0] if pareto_id.startswith("N") else "",
                })
        write_csv(rows, output_dir / f"ideal_csi_bler_{int(n_prbs)}prb_partial.csv")

    for mode in modes:
        plot_rows = []
        for row in rows:
            if row["interleaver_mode"] != mode:
                continue
            item = dict(row)
            item["bler_plot"] = max(float(row["bler"]), 0.5 / float(row["trials"]))
            plot_rows.append(item)
        svg_all_pareto_plot(
            plot_rows,
            metric="bler_plot",
            title=f"Ideal-CSI BLER, {n_prbs} PRB, interleaver={mode}",
            ylabel="BLER",
            path=figure_dir / f"experiment21_ideal_csi_bler_{n_prbs}prb_{mode}.svg",
            log_y=True,
            threshold_lines=(0.1, 0.01),
            label_endpoints=False,
            direct_labels=True,
        )
    return rows


def run_estimated_bler_for_bandwidth(
    args: argparse.Namespace,
    n_prbs: int,
    resource: ResourceConfig,
    grid: ResourceGrid,
    selected: Sequence[Tuple[str, VDesign, Dict[str, object]]],
    output_dir: Path,
    figure_dir: Path,
) -> List[Dict[str, object]]:
    if int(args.estimated_trials) <= 0:
        return []
    modes = parse_text_list(args.estimated_interleaver_modes)
    mcs = get_mcs("nr_256qam", int(args.mcs_index), None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(
        tb.cb_k_values,
        tb.cb_e_values,
        num_iter=int(args.ldpc_iterations),
        llr_clip=float(args.llr_clip),
    )
    sample_period_ns = 1e9 / (resource.n_fft * resource.scs_khz * 1e3)
    pdp = make_exponential_pdp(
        delay_spread_ns=float(args.delay_spread_ns),
        sample_period_ns=sample_period_ns,
        max_delay_factor=float(args.max_delay_factor),
    )
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    channel = ChannelConfig(
        delay_spread_ns=float(args.delay_spread_ns),
        max_delay_factor=float(args.max_delay_factor),
    )
    permutation, interleaver_seed = fixed_interleaver(tb.coded_bits, n_prbs, int(args.seed))
    np.save(output_dir / f"coded_bit_interleaver_{int(n_prbs)}prb.npy", permutation)

    rows: List[Dict[str, object]] = []
    for snr_db in parse_float_list(args.estimated_snrs):
        noise_var = 10.0 ** (-float(snr_db) / 10.0)
        noise_var_ls = noise_var / float(len(resource.dmrs_symbol_indices))
        estimators = {
            pareto_id: make_alg1_estimator(
                grid, pdp, design, noise_var_ls, float(args.loading)
            )
            for pareto_id, design, _ in selected
        }
        accum = {
            (mode, pareto_id): {
                "tb_errors": 0,
                "cb_errors": 0,
                "ce_nmse": [],
            }
            for mode in modes
            for pareto_id, _, _ in selected
        }
        for trial in range(1, int(args.estimated_trials) + 1):
            rng = np.random.default_rng(
                stable_seed("new-front-estimated", n_prbs, trial, base=int(args.seed))
            )
            realization = generate_tdl_channel(rng, grid, channel, n_tx=8, n_rx=1)
            payload = [rng.integers(0, 2, size=int(k), dtype=np.int8) for k in tb.cb_k_values]
            coded_bits = np.concatenate(adapter.encode(payload))
            pilot_noise = (
                rng.normal(size=(1, len(pilot_local)))
                + 1j * rng.normal(size=(1, len(pilot_local)))
            ) * math.sqrt(noise_var_ls / 2.0)
            data_noise = (
                rng.normal(size=(1, len(data_local)))
                + 1j * rng.normal(size=(1, len(data_local)))
            ) * math.sqrt(noise_var / 2.0)

            llrs: List[np.ndarray] = []
            trial_keys: List[Tuple[str, str]] = []
            trial_nmse: List[float] = []
            for mode in modes:
                mapped_bits = apply_interleaver(coded_bits, permutation, mode)
                symbols = qam_modulate(mapped_bits, int(mcs.qm))
                for pareto_id, design, _ in selected:
                    g = equivalent_channel(realization.H, design.C)
                    ls_observation = g[:, pilot_local] + pilot_noise
                    g_hat = ls_observation @ estimators[pareto_id].matrix.T
                    y = g[:, data_local] * symbols[None, :] + data_noise
                    z, no_eff = equalize_mrc(y, g_hat[:, data_local], noise_var)
                    llr_mapped = qam_demapper_maxlog(z, no_eff, int(mcs.qm))
                    llrs.append(undo_interleaver(llr_mapped, permutation, mode))
                    trial_keys.append((mode, pareto_id))
                    trial_nmse.append(nmse(g, g_hat))
            decoded = decode_same_tb_batch(adapter, llrs, payload)
            for key, result, ce_nmse in zip(trial_keys, decoded, trial_nmse):
                accum[key]["tb_errors"] += int(not result.tb_success)
                accum[key]["cb_errors"] += int(sum(1 for ok in result.cb_success if not ok))
                accum[key]["ce_nmse"].append(float(ce_nmse))
            if int(args.progress_every) > 0 and trial % int(args.progress_every) == 0:
                print(
                    f"[estimated] PRB={n_prbs} SNR={snr_db:g} "
                    f"trial={trial}/{int(args.estimated_trials)}",
                    flush=True,
                )

        for mode in modes:
            for pareto_id, design, _ in selected:
                result = accum[(mode, pareto_id)]
                rows.append({
                    "bandwidth_prb": int(n_prbs),
                    "active_subcarriers": int(grid.n_sc),
                    "occupied_bandwidth_mhz": float(grid.n_sc * resource.scs_khz / 1000.0),
                    "segment_length_sc": int(grid.n_sc // 8),
                    "segment_bandwidth_mhz": float((grid.n_sc // 8) * resource.scs_khz / 1000.0),
                    "interleaver_mode": mode,
                    "interleaver_seed": int(interleaver_seed) if mode == "full" else "",
                    "pareto_id": pareto_id,
                    "precoder": design.label,
                    "family": curve_family(pareto_id),
                    "snr_db": float(snr_db),
                    "trials": int(args.estimated_trials),
                    "tb_errors": int(result["tb_errors"]),
                    "cb_errors": int(result["cb_errors"]),
                    "bler": float(result["tb_errors"]) / float(args.estimated_trials),
                    "ce_nmse_mean": mean_or_nan(result["ce_nmse"]),
                    "ce_nmse_median": float(np.median(result["ce_nmse"])),
                    "rpp_loaded_cond": float(estimators[pareto_id].cond),
                    "estimator": "Alg1 matched full-covariance RMMSE",
                    "csi_mode": "estimated",
                    "reference_cdd_id": cdd_reference(n_prbs, pareto_id)[0] if pareto_id.startswith("N") else "",
                    "mcs_index": int(args.mcs_index),
                    "qm": int(mcs.qm),
                    "code_rate": float(mcs.code_rate),
                    "tbs_bits": int(tb.tb_size),
                    "coded_bits": int(tb.coded_bits),
                    "data_re": int(grid.n_data_re),
                    "dmrs_spacing_sc": int(args.dmrs_spacing_sc),
                    "delay_spread_ns": float(args.delay_spread_ns),
                })
        write_csv(rows, output_dir / f"estimated_csi_bler_{int(n_prbs)}prb_partial.csv")

    for mode in modes:
        plot_rows = [row for row in rows if row["interleaver_mode"] == mode]
        svg_all_pareto_plot(
            plot_rows,
            metric="bler",
            title=f"Algorithm 1 estimated-CSI BLER, {n_prbs} PRB, {mode} interleaver",
            ylabel="BLER",
            path=figure_dir / f"experiment22_estimated_csi_bler_{n_prbs}prb_{mode}.svg",
            log_y=True,
            threshold_lines=(0.1, 0.01),
            label_endpoints=False,
            direct_labels=True,
        )
    return rows


def ideal_target_rows(rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    output: List[Dict[str, object]] = []
    groups = sorted({
        (int(row["bandwidth_prb"]), str(row["interleaver_mode"]), str(row["pareto_id"]))
        for row in rows
    })
    lookup: Dict[Tuple[int, str, str], Dict[str, object]] = {}
    for n_prbs, mode, pareto_id in groups:
        points = [
            row for row in rows
            if int(row["bandwidth_prb"]) == n_prbs
            and str(row["interleaver_mode"]) == mode
            and str(row["pareto_id"]) == pareto_id
        ]
        item = {
            "bandwidth_prb": n_prbs,
            "active_subcarriers": points[0]["active_subcarriers"],
            "occupied_bandwidth_mhz": points[0]["occupied_bandwidth_mhz"],
            "interleaver_mode": mode,
            "pareto_id": pareto_id,
            "precoder": points[0]["precoder"],
            "snr_at_bler_0p1_db": interp_target_snr(points, 0.1),
            "snr_at_bler_0p01_db": interp_target_snr(points, 0.01),
            "min_observed_bler": min(float(row["bler"]) for row in points),
        }
        lookup[(n_prbs, mode, pareto_id)] = item
        output.append(item)
    for item in output:
        pareto_id = str(item["pareto_id"])
        if not pareto_id.startswith("N"):
            item["reference_cdd_id"] = ""
            item["gain_at_bler_0p1_db"] = ""
            item["gain_at_bler_0p01_db"] = ""
            continue
        reference_id = cdd_reference(int(item["bandwidth_prb"]), pareto_id)[0]
        item["reference_cdd_id"] = reference_id
        ref = lookup.get((int(item["bandwidth_prb"]), str(item["interleaver_mode"]), reference_id))
        if ref is None:
            item["gain_at_bler_0p1_db"] = ""
            item["gain_at_bler_0p01_db"] = ""
            continue
        for suffix in ("0p1", "0p01"):
            n_snr = float(item[f"snr_at_bler_{suffix}_db"])
            c_snr = float(ref[f"snr_at_bler_{suffix}_db"])
            item[f"gain_at_bler_{suffix}_db"] = (
                c_snr - n_snr if np.isfinite(n_snr) and np.isfinite(c_snr) else ""
            )
    return output


def run_nmse_for_bandwidth(
    args: argparse.Namespace,
    n_prbs: int,
    resource: ResourceConfig,
    grid: ResourceGrid,
    selected: Sequence[Tuple[str, VDesign, Dict[str, object]]],
    output_dir: Path,
    figure_dir: Path,
) -> List[Dict[str, object]]:
    if int(args.nmse_trials) <= 0:
        return []
    sample_period_ns = 1e9 / (resource.n_fft * resource.scs_khz * 1e3)
    pdp = make_exponential_pdp(
        delay_spread_ns=float(args.delay_spread_ns),
        sample_period_ns=sample_period_ns,
        max_delay_factor=float(args.max_delay_factor),
    )
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    channel = ChannelConfig(
        delay_spread_ns=float(args.delay_spread_ns),
        max_delay_factor=float(args.max_delay_factor),
    )
    rows: List[Dict[str, object]] = []
    for snr_db in parse_float_list(args.nmse_snrs):
        noise_var = 10.0 ** (-float(snr_db) / 10.0)
        noise_var_ls = noise_var / float(len(resource.dmrs_symbol_indices))
        estimators = {
            pareto_id: make_alg1_estimator(grid, pdp, design, noise_var_ls, float(args.loading))
            for pareto_id, design, _ in selected
        }
        values = {pareto_id: [] for pareto_id, _, _ in selected}
        for trial in range(1, int(args.nmse_trials) + 1):
            rng = np.random.default_rng(
                stable_seed("new-front-nmse", n_prbs, trial, base=int(args.seed))
            )
            realization = generate_tdl_channel(rng, grid, channel, n_tx=8, n_rx=1)
            pilot_noise = (
                rng.normal(size=(1, len(pilot_local)))
                + 1j * rng.normal(size=(1, len(pilot_local)))
            ) * math.sqrt(noise_var_ls / 2.0)
            for pareto_id, design, _ in selected:
                g = equivalent_channel(realization.H, design.C)
                ls_observation = g[:, pilot_local] + pilot_noise
                g_hat = ls_observation @ estimators[pareto_id].matrix.T
                values[pareto_id].append(nmse(g, g_hat))
        for pareto_id, design, _ in selected:
            rows.append({
                "bandwidth_prb": int(n_prbs),
                "active_subcarriers": int(grid.n_sc),
                "occupied_bandwidth_mhz": float(grid.n_sc * resource.scs_khz / 1000.0),
                "pareto_id": pareto_id,
                "precoder": design.label,
                "family": curve_family(pareto_id),
                "snr_db": float(snr_db),
                "trials": int(args.nmse_trials),
                "ce_nmse_mean": mean_or_nan(values[pareto_id]),
                "ce_nmse_median": float(np.median(values[pareto_id])),
                "rpp_loaded_cond": float(estimators[pareto_id].cond),
                "estimator": "Alg1 matched full-covariance RMMSE",
                "dmrs_spacing_sc": int(args.dmrs_spacing_sc),
                "delay_spread_ns": float(args.delay_spread_ns),
                "reference_cdd_id": cdd_reference(n_prbs, pareto_id)[0] if pareto_id.startswith("N") else "",
            })
        write_csv(rows, output_dir / f"mmse_nmse_{int(n_prbs)}prb_partial.csv")
    svg_all_pareto_plot(
        rows,
        metric="ce_nmse_mean",
        title=f"Matched RMMSE NMSE, {n_prbs} PRB",
        ylabel="NMSE (log scale)",
        path=figure_dir / f"experiment21_mmse_nmse_{n_prbs}prb.svg",
        log_y=True,
        label_endpoints=False,
    )
    return rows


def nmse_relative_rows(rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    lookup = {
        (int(row["bandwidth_prb"]), float(row["snr_db"]), str(row["pareto_id"])): row
        for row in rows
    }
    output = []
    for row in rows:
        pareto_id = str(row["pareto_id"])
        if not pareto_id.startswith("N"):
            continue
        reference_id = cdd_reference(int(row["bandwidth_prb"]), pareto_id)[0]
        ref = lookup[(int(row["bandwidth_prb"]), float(row["snr_db"]), reference_id)]
        n_nmse = float(row["ce_nmse_mean"])
        c_nmse = float(ref["ce_nmse_mean"])
        output.append({
            "bandwidth_prb": row["bandwidth_prb"],
            "occupied_bandwidth_mhz": row["occupied_bandwidth_mhz"],
            "snr_db": row["snr_db"],
            "pareto_id": pareto_id,
            "reference_cdd_id": reference_id,
            "n_nmse": n_nmse,
            "cdd_nmse": c_nmse,
            "n_nmse_penalty_db": 10.0 * math.log10(n_nmse / c_nmse),
        })
    return output


def run(args: argparse.Namespace) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path(args.out) / f"v_design_new_front_link_{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir = Path(args.fig_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)
    selected_ids = parse_text_list(args.selected_ids)
    additional_cdd_steps = parse_int_list(args.additional_cdd_steps)
    scan_dir = Path(args.scan_dir)

    all_ideal_rows: List[Dict[str, object]] = []
    all_estimated_rows: List[Dict[str, object]] = []
    all_nmse_rows: List[Dict[str, object]] = []
    all_design_metric_rows: List[Dict[str, object]] = []
    bandwidth_rows: List[Dict[str, object]] = []
    design_rows: List[Dict[str, object]] = []
    for n_prbs in parse_int_list(args.n_prbs):
        resource, grid = build_grid(n_prbs, int(args.dmrs_spacing_sc))
        selected = load_selected_designs(
            grid,
            scan_dir,
            selected_ids,
            additional_cdd_steps=additional_cdd_steps,
            include_reference_cdd=not bool(args.skip_reference_cdd),
        )
        sample_period_ns = 1e9 / (resource.n_fft * resource.scs_khz * 1e3)
        metric_pdp = make_exponential_pdp(
            delay_spread_ns=float(args.delay_spread_ns),
            sample_period_ns=sample_period_ns,
            max_delay_factor=float(args.max_delay_factor),
        )
        base_covariance = covariance_matrix(
            grid.subcarrier_indices,
            grid.subcarrier_indices,
            metric_pdp,
            grid.n_fft,
        )
        for pareto_id, design, _ in selected:
            matrix_row = matrix_metrics(design, grid, 8)
            covariance_row, _, _ = covariance_metrics(
                design,
                base_covariance,
                grid.n_sc - 1,
            )
            right_censored = int(covariance_row["coherence_bw_abs_0p5_sc"] == "")
            if right_censored:
                covariance_row["coherence_bw_abs_0p5_sc"] = int(grid.n_sc)
            item: Dict[str, object] = {
                "bandwidth_prb": int(n_prbs),
                "active_subcarriers": int(grid.n_sc),
                "occupied_bandwidth_mhz": float(grid.n_sc * resource.scs_khz / 1000.0),
                "pareto_id": pareto_id,
                "precoder": design.label,
                "family": curve_family(pareto_id),
                "coherence_bw_abs_0p5_right_censored": right_censored,
            }
            item.update(matrix_row)
            item.update(covariance_row)
            all_design_metric_rows.append(item)
        mcs = get_mcs("nr_256qam", int(args.mcs_index), None, None)
        tb = build_tb_layout(grid.n_data_re, mcs)
        bandwidth_rows.append({
            "bandwidth_prb": int(n_prbs),
            "active_subcarriers": int(grid.n_sc),
            "occupied_bandwidth_mhz": float(grid.n_sc * resource.scs_khz / 1000.0),
            "segment_length_sc": int(grid.n_sc // 8),
            "segment_bandwidth_mhz": float((grid.n_sc // 8) * resource.scs_khz / 1000.0),
            "dmrs_spacing_sc": int(args.dmrs_spacing_sc),
            "pilot_subcarriers_per_symbol": int(grid.pilot_count),
            "data_re": int(grid.n_data_re),
            "coded_bits": int(tb.coded_bits),
            "tbs_bits": int(tb.tb_size),
            "n_code_blocks": int(len(tb.cb_k_values)),
            "cb_k_values": json.dumps(tb.cb_k_values),
            "cb_e_values": json.dumps(tb.cb_e_values),
        })
        if not design_rows:
            for pareto_id, design, manifest in selected:
                design_rows.append({
                    "pareto_id": pareto_id,
                    "precoder": design.label,
                    "family": curve_family(pareto_id),
                    "reference_cdd_id": cdd_reference(n_prbs, pareto_id)[0] if pareto_id.startswith("N") else "",
                    "delay_alphabet_samples": manifest.get("delay_alphabet_samples", ""),
                    "layout_id": manifest.get("layout_id", ""),
                    "layout_seed": manifest.get("layout_seed", ""),
                    "segment_tx_delay_samples": manifest.get("segment_tx_delay_samples", ""),
                })
        all_ideal_rows.extend(
            run_ideal_bler_for_bandwidth(
                args, n_prbs, resource, grid, selected, output_dir, figure_dir
            )
        )
        all_estimated_rows.extend(
            run_estimated_bler_for_bandwidth(
                args, n_prbs, resource, grid, selected, output_dir, figure_dir
            )
        )
        all_nmse_rows.extend(
            run_nmse_for_bandwidth(
                args, n_prbs, resource, grid, selected, output_dir, figure_dir
            )
        )

    write_csv(all_ideal_rows, output_dir / "ideal_csi_bler.csv")
    target_rows = ideal_target_rows(all_ideal_rows) if all_ideal_rows else []
    write_csv(target_rows, output_dir / "ideal_csi_bler_targets.csv")
    write_csv(all_estimated_rows, output_dir / "estimated_csi_bler.csv")
    estimated_target_rows = ideal_target_rows(all_estimated_rows) if all_estimated_rows else []
    write_csv(estimated_target_rows, output_dir / "estimated_csi_bler_targets.csv")
    write_csv(all_nmse_rows, output_dir / "mmse_nmse.csv")
    write_csv(nmse_relative_rows(all_nmse_rows), output_dir / "mmse_nmse_relative.csv")
    write_csv(all_design_metric_rows, output_dir / "selected_design_metrics_by_bandwidth.csv")
    write_csv(bandwidth_rows, output_dir / "bandwidth_configs.csv")
    write_csv(design_rows, output_dir / "selected_designs.csv")
    with (output_dir / "run_config.json").open("w", encoding="utf-8") as handle:
        json.dump(vars(args), handle, indent=2, ensure_ascii=False)
    print(output_dir)
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan-dir", default=str(SCAN_DIR))
    parser.add_argument("--selected-ids", default="N3,N6,N8,N23")
    parser.add_argument("--additional-cdd-steps", default="")
    parser.add_argument("--skip-reference-cdd", action="store_true")
    parser.add_argument("--n-prbs", default="48,36")
    parser.add_argument("--dmrs-spacing-sc", type=int, default=24)
    parser.add_argument("--delay-spread-ns", type=float, default=5.0)
    parser.add_argument("--max-delay-factor", type=float, default=8.0)
    parser.add_argument("--mcs-index", type=int, default=8)
    parser.add_argument("--ideal-snrs", default="8,10,12,14,16,18")
    parser.add_argument("--ideal-trials", type=int, default=400)
    parser.add_argument("--interleaver-modes", default="full,none")
    parser.add_argument("--estimated-snrs", default="")
    parser.add_argument("--estimated-trials", type=int, default=0)
    parser.add_argument("--estimated-interleaver-modes", default="full")
    parser.add_argument("--nmse-snrs", default="0,4,8,12,16,20")
    parser.add_argument("--nmse-trials", type=int, default=300)
    parser.add_argument("--ldpc-iterations", type=int, default=8)
    parser.add_argument("--llr-clip", type=float, default=50.0)
    parser.add_argument("--loading", type=float, default=1e-8)
    parser.add_argument("--seed", type=int, default=20260619)
    parser.add_argument("--progress-every", type=int, default=50)
    parser.add_argument("--out", default="outputs/v_design_new_front_link")
    parser.add_argument("--fig-dir", default="docs/figures")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
