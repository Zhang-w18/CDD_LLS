"""Run the plan-024 QC/Sidon comparison on the plan-025 Sionna TDL chain."""

from __future__ import annotations

import argparse
import csv
import importlib.metadata
import json
import math
import platform
import sys
import time
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import (
    AntennaConfig,
    ChannelConfig,
    ChannelEstimationConfig,
    MCSConfig,
    PlatformConfig,
    PlotConfig,
    ReceiverConfig,
    ResourceConfig,
    SimulationConfig,
    TransmissionConfig,
)
from cdd_lls.sim.orchestrator import CDDLinkLevelOrchestrator
import numpy as np

from cdd_lls.core.mcs import build_tb_layout, get_mcs
from cdd_lls.phy.channel_tdl import generate_sionna_tdl_channel
from cdd_lls.phy.estimators import build_time_frequency_rmmse_filter, tdl_known_delay_covariance
from cdd_lls.phy.ldpc import SionnaLDPCAdapter
from cdd_lls.phy.precoding import build_precoder, equivalent_channel
from cdd_lls.phy.qam import qam_demapper_maxlog, qam_modulate
from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
from tools.run_v_design_piecewise_tradeoff import decode_same_tb_batch


DESIGNS = {
    "QC": [0, 9, 18, 27, 36, 45, 54, 63],
    "Sidon": [0, 1, 3, 7, 12, 20, 30, 65],
}
SEED = 20260722


def parse_float_list(text: str) -> List[float]:
    return [float(item.strip()) for item in text.split(",") if item.strip()]


def stable_seed(*items: object) -> int:
    text = "|".join(str(item) for item in items)
    acc = 2166136261
    for char in text:
        acc ^= ord(char)
        acc = (acc * 16777619) % (2**32)
    return int(acc % (2**31 - 1))


def write_rows(rows: List[Dict[str, object]], path: Path) -> None:
    if not rows:
        return
    fields = sorted(set().union(*(row.keys() for row in rows)))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run_known_comparison(
    snrs_db: Sequence[float],
    trials: int,
    batch_size: int,
    output: Path,
) -> Path:
    output.mkdir(parents=True, exist_ok=True)
    resource = ResourceConfig(
        carrier_bandwidth_mhz=100.0,
        scs_khz=30,
        n_fft=4096,
        n_prbs=48,
        pdsch_n_symbols=10,
        dmrs_symbol_indices=[2, 7],
        dmrs_spacing_sc=24,
        dmrs_offset_sc=0,
        cyclic_prefix_length=288,
    )
    channel = ChannelConfig(
        backend="sionna_tdl",
        model="3gpp_tr38901_tdl",
        tdl_profile="A",
        delay_spread_ns=5.0,
        carrier_frequency_hz=3.5e9,
        ue_speed_kmh=0.0,
        num_sinusoids=20,
    )
    grid = build_resource_grid(resource)
    mcs = get_mcs("nr_256qam", 8, None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(tb.cb_k_values, tb.cb_e_values, num_iter=8, llr_clip=50.0)
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarrier_indices)
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    precoders = {
        design: build_precoder(
            grid,
            resource,
            TransmissionConfig(tx_scheme="CDD", cdd_delay_vector=delays),
            n_tx=8,
        )
        for design, delays in DESIGNS.items()
    }
    result_path = output / "sidon_qc_bler.csv"
    paired_path = output / "paired_error_counts.csv"
    if result_path.exists():
        with result_path.open("r", encoding="utf-8", newline="") as handle:
            rows: List[Dict[str, object]] = list(csv.DictReader(handle))
    else:
        rows = []
    if paired_path.exists():
        with paired_path.open("r", encoding="utf-8", newline="") as handle:
            paired_rows: List[Dict[str, object]] = list(csv.DictReader(handle))
    else:
        paired_rows = []

    receipt = {
        "experiment": "plan-025 migration of plan-024 E2",
        "channel": {
            "backend": "sionna_tdl",
            "profile": "A",
            "delay_spread_ns": 5.0,
            "carrier_frequency_hz": 3.5e9,
            "ue_speed_kmh": 0.0,
            "normalize_per_realization": False,
        },
        "resource": {
            "n_prbs": 48,
            "n_active_subcarriers": 576,
            "n_symbols": 10,
            "dmrs_symbol_indices": [2, 7],
            "dmrs_spacing_sc": 24,
            "n_dmrs_re": grid.n_dmrs_re,
            "n_data_re": grid.n_data_re,
            "n_fft": 4096,
            "cyclic_prefix_length": 288,
            "ofdm_symbol_duration_s": grid.ofdm_symbol_duration_s,
        },
        "designs": DESIGNS,
        "receiver": "known-delay matched 2D RMMSE",
        "noise_variance_ls": "N0 for every DMRS RE; no pre-averaging",
        "mcs": {"modulation": "16QAM", "index": 8, "code_rate": 553 / 1024},
        "ldpc_iterations": 8,
        "snrs_db": [float(value) for value in snrs_db],
        "trials_per_snr": int(trials),
        "batch_size": int(batch_size),
        "seed": SEED,
        "common_random_numbers": ["TDL realization", "payload", "LS noise", "data noise"],
    }
    (output / "resolved_experiment.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    package_names = ("numpy", "scipy", "matplotlib", "PyYAML", "tensorflow", "sionna")
    environment_receipt = {
        "python": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "packages": {name: importlib.metadata.version(name) for name in package_names},
        "command": [sys.executable, *sys.argv],
    }
    (output / "environment.json").write_text(
        json.dumps(environment_receipt, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    for snr_db in snrs_db:
        completed = [
            row for row in rows
            if float(row["snr_db"]) == float(snr_db)
            and int(row.get("trials", row.get("n_trials", 0))) == int(trials)
            and str(row["id"]) in DESIGNS
        ]
        if {str(row["id"]) for row in completed} == set(DESIGNS):
            print(f"[TDL E2] resume: skip completed snr={snr_db:g} ({trials} trials)", flush=True)
            continue
        rows = [row for row in rows if float(row["snr_db"]) != float(snr_db)]
        paired_rows = [row for row in paired_rows if float(row["snr_db"]) != float(snr_db)]
        noise_variance = 10.0 ** (-float(snr_db) / 10.0)
        filters = {}
        for design, delays in DESIGNS.items():
            covariance = tdl_known_delay_covariance(grid, channel, delays)
            filters[design] = build_time_frequency_rmmse_filter(
                grid, covariance, noise_variance, diagonal_loading=1e-10,
            )
        errors = {design: 0 for design in DESIGNS}
        nmse_sum = {design: 0.0 for design in DESIGNS}
        paired = {"both_error": 0, "qc_only_error": 0, "sidon_only_error": 0, "both_success": 0}
        started = time.time()
        completed_trials = 0
        for batch_start in range(1, int(trials) + 1, int(batch_size)):
            current = min(int(batch_size), int(trials) - batch_start + 1)
            payload_rng = np.random.default_rng(stable_seed(SEED, snr_db, batch_start, "payload"))
            payload = [
                payload_rng.integers(0, 2, size=int(k), dtype=np.int8)
                for k in tb.cb_k_values
            ]
            coded = np.concatenate(adapter.encode(payload))
            symbols = qam_modulate(coded, int(mcs.qm))
            realization = generate_sionna_tdl_channel(
                grid,
                channel,
                n_tx=8,
                n_rx=1,
                batch_size=current,
                seed=stable_seed(SEED, snr_db, batch_start, "channel"),
            )
            noise_rng = np.random.default_rng(stable_seed(SEED, snr_db, batch_start, "noise"))
            pilot_noise = math.sqrt(noise_variance / 2.0) * (
                noise_rng.normal(size=(current, 1, grid.pilot_re_count))
                + 1j * noise_rng.normal(size=(current, 1, grid.pilot_re_count))
            )
            data_noise = math.sqrt(noise_variance / 2.0) * (
                noise_rng.normal(size=(current, 1, grid.n_data_re))
                + 1j * noise_rng.normal(size=(current, 1, grid.n_data_re))
            )
            decoded_by_design = {}
            for design, precoder in precoders.items():
                g = equivalent_channel(realization.H, precoder.C)[:, 0:1]
                true_pilot = g[:, :, grid.pilot_symbol_indices, pilot_local]
                true_data = g[:, :, grid.data_symbol_indices, data_local]
                estimate = filters[design].estimate_data(true_pilot + pilot_noise)
                per_trial_nmse = np.sum(np.abs(estimate - true_data) ** 2, axis=(1, 2)) / np.maximum(
                    np.sum(np.abs(true_data) ** 2, axis=(1, 2)), 1e-30
                )
                nmse_sum[design] += float(np.sum(per_trial_nmse))
                y = true_data * symbols[None, None, :] + data_noise
                denominator = np.maximum(np.sum(np.abs(estimate) ** 2, axis=1), 1e-10)
                equalized = np.sum(np.conj(estimate) * y, axis=1) / denominator
                effective_noise = noise_variance / denominator
                llrs = [
                    qam_demapper_maxlog(equalized[index], effective_noise[index], int(mcs.qm))
                    for index in range(current)
                ]
                decoded = decode_same_tb_batch(adapter, llrs, payload)
                decoded_by_design[design] = decoded
                errors[design] += sum(int(not result.tb_success) for result in decoded)
            for qc_result, sidon_result in zip(decoded_by_design["QC"], decoded_by_design["Sidon"]):
                qc_error = not qc_result.tb_success
                sidon_error = not sidon_result.tb_success
                if qc_error and sidon_error:
                    paired["both_error"] += 1
                elif qc_error:
                    paired["qc_only_error"] += 1
                elif sidon_error:
                    paired["sidon_only_error"] += 1
                else:
                    paired["both_success"] += 1
            completed_trials += current
            if completed_trials % 100 < current or completed_trials == int(trials):
                print(
                    f"[TDL E2] snr={snr_db:g} trial={completed_trials}/{trials} "
                    f"errors={errors} elapsed={time.time()-started:.1f}s",
                    flush=True,
                )

        for design in DESIGNS:
            lo, hi = wilson(errors[design], int(trials))
            estimator = filters[design]
            rows.append({
                "id": design,
                "variant_id": f"{design}_known",
                "branch": "known",
                "snr_db": float(snr_db),
                "trials": int(trials),
                "n_trials": int(trials),
                "tb_errors": int(errors[design]),
                "bler": float(errors[design]) / float(trials),
                "bler_wilson95_lo": lo,
                "bler_wilson95_hi": hi,
                "ce_nmse_mean": nmse_sum[design] / float(trials),
                "ce_nmse_mean_dB": 10.0 * math.log10(nmse_sum[design] / float(trials)),
                "condition_number": estimator.condition_number,
                "numerical_jitter": estimator.numerical_jitter,
                "covariance_type": estimator.covariance_type,
                "noise_variance": noise_variance,
                "noise_variance_ls": noise_variance,
            })
        paired_rows.append({"snr_db": float(snr_db), "trials": int(trials), **paired})
        rows.sort(key=lambda row: (float(row["snr_db"]), str(row["id"])))
        paired_rows.sort(key=lambda row: float(row["snr_db"]))
        write_rows(rows, result_path)
        write_rows(rows, output / "summary.csv")
        write_rows(paired_rows, paired_path)
    return output


def wilson(errors: int, trials: int, z: float = 1.96) -> Tuple[float, float]:
    if trials <= 0:
        return float("nan"), float("nan")
    p = errors / trials
    denominator = 1.0 + z * z / trials
    center = (p + z * z / (2.0 * trials)) / denominator
    half = z * math.sqrt(p * (1.0 - p) / trials + z * z / (4.0 * trials * trials)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def target_snr(rows: Sequence[Dict[str, object]], variant_id: str, target: float, field: str = "bler") -> float:
    points = sorted(
        (float(row["snr_db"]), max(float(row[field]), 1e-12))
        for row in rows if str(row["variant_id"]) == variant_id
    )
    for (s0, b0), (s1, b1) in zip(points[:-1], points[1:]):
        if (b0 - target) * (b1 - target) <= 0.0 and b0 != b1:
            alpha = (math.log10(target) - math.log10(b0)) / (math.log10(b1) - math.log10(b0))
            return s0 + alpha * (s1 - s0)
    return float("nan")


def analyze(output: Path) -> Dict[str, object]:
    with (output / "summary.csv").open("r", encoding="utf-8", newline="") as handle:
        rows: List[Dict[str, object]] = list(csv.DictReader(handle))
    augmented = []
    for row in rows:
        lo, hi = wilson(int(row["tb_errors"]), int(row["n_trials"]))
        item = dict(row)
        item["bler_wilson95_lo"] = lo
        item["bler_wilson95_hi"] = hi
        augmented.append(item)
    fields = sorted(set().union(*(row.keys() for row in augmented)))
    with (output / "bler_with_intervals.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(augmented)

    variants = sorted({str(row["variant_id"]) for row in augmented})
    targets: Dict[str, object] = {}
    for target in (0.10, 0.01):
        label = "bler10" if target == 0.10 else "bler1"
        targets[label] = {}
        for variant in variants:
            errors_near_target = sum(
                int(row["tb_errors"]) for row in augmented
                if str(row["variant_id"]) == variant
                and 0.5 * target <= float(row["bler"]) <= 2.0 * target
            )
            targets[label][variant] = {
                "snr_db": target_snr(augmented, variant, target),
                "errors_in_half_to_twice_target_band": errors_near_target,
                "sufficient_for_deterministic_claim": errors_near_target >= 30,
            }
    for branch in ("known", "unknown"):
        qc = f"QC_{branch}"
        sidon = f"Sidon_{branch}"
        if qc in variants and sidon in variants:
            for label in ("bler10", "bler1"):
                q = float(targets[label][qc]["snr_db"])
                s = float(targets[label][sidon]["snr_db"])
                targets[label][f"Sidon_gain_vs_QC_{branch}_dB"] = q - s
    summary = {
        "targets": targets,
        "confidence_interval_note": (
            "Per-point BLER uses Wilson 95% intervals. Target-SNR difference intervals require "
            "the formal analysis/refinement stage and are not inferred when crossings are absent."
        ),
    }
    (output / "tdl_comparison_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--snrs", default="13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18")
    parser.add_argument("--trials", type=int, default=400)
    parser.add_argument("--branches", choices=("known", "unknown", "both"), default="known")
    parser.add_argument("--batch-size", type=int, default=20)
    parser.add_argument("--run-id", default="20260722_main/prescan")
    parser.add_argument(
        "--output-root", type=Path,
        default=ROOT / "outputs" / "experiment025_sionna_tdl_rmmse",
    )
    parser.add_argument("--analyze-only", type=Path)
    args = parser.parse_args()
    if args.analyze_only is not None:
        print(json.dumps(analyze(args.analyze_only), indent=2, ensure_ascii=False))
        return

    if args.branches == "known":
        output = args.output_root / Path(str(args.run_id))
        run_known_comparison(
            parse_float_list(args.snrs),
            trials=int(args.trials),
            batch_size=int(args.batch_size),
            output=output,
        )
        summary = analyze(output)
        print(f"Results: {output}")
        print(json.dumps(summary, indent=2, ensure_ascii=False))
        return

    branches = ("known", "unknown") if args.branches == "both" else (args.branches,)
    variants: List[Dict[str, object]] = []
    for design, delays in DESIGNS.items():
        for branch in branches:
            variants.append({
                "variant_id": f"{design}_{branch}",
                "transmission": {"tx_scheme": "CDD", "cdd_delay_vector": delays},
                "channel_estimation": {
                    "ce_method": "TF_RMMSE_KNOWN" if branch == "known" else "TF_RMMSE_UNKNOWN"
                },
            })
    config = PlatformConfig(
        antenna=AntennaConfig(n_tx=8, n_rx=1),
        resource=ResourceConfig(
            carrier_bandwidth_mhz=100.0,
            scs_khz=30,
            n_fft=4096,
            n_prbs=48,
            pdsch_n_symbols=10,
            dmrs_symbol_indices=[2, 7],
            dmrs_spacing_sc=24,
            dmrs_offset_sc=0,
            cyclic_prefix_length=288,
        ),
        channel=ChannelConfig(
            backend="sionna_tdl",
            model="3gpp_tr38901_tdl",
            tdl_profile="A",
            delay_spread_ns=5.0,
            carrier_frequency_hz=3.5e9,
            ue_speed_kmh=0.0,
            num_sinusoids=20,
        ),
        transmission=TransmissionConfig(
            tx_scheme="CDD",
            cdd_delay_vector=DESIGNS["QC"],
        ),
        channel_estimation=ChannelEstimationConfig(
            ce_method="TF_RMMSE_KNOWN",
            diagonal_loading=1e-10,
        ),
        receiver=ReceiverConfig(max_ldpc_iterations=8, llr_clip=50.0),
        mcs=MCSConfig(table="nr_256qam", index=8),
        simulation=SimulationConfig(
            snr_points_db=parse_float_list(args.snrs),
            n_trials_per_snr=int(args.trials),
            max_trials_per_snr=int(args.trials),
            min_block_errors=0,
            seed=20260722,
            output_dir=str(args.output_root),
            run_id=str(args.run_id),
            common_random_numbers=True,
            save_trial_metrics=True,
        ),
        plots=PlotConfig(enabled=True),
        variants=variants,
    )
    output = CDDLinkLevelOrchestrator(config).run()
    summary = analyze(output)
    print(f"Results: {output}")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
