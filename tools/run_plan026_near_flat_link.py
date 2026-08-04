"""Run plan-026 E3 paired estimated-CSI BLER experiments."""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ChannelConfig, ResourceConfig
from cdd_lls.core.mcs import build_tb_layout, get_mcs
from cdd_lls.phy.channel_tdl import generate_sionna_tdl_channel
from cdd_lls.phy.estimators import (
    build_frequency_rmmse_filter,
    tdl_active_frequency_covariance,
)
from cdd_lls.phy.ldpc import SionnaLDPCAdapter
from cdd_lls.phy.precoding import build_active_dft_grid_precoder, equivalent_channel
from cdd_lls.phy.qam import qam_demapper_maxlog, qam_modulate
from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
from tools.analyze_plan025_sionna_e2 import (
    paired_error_analysis,
    target_summary,
)
from tools.run_plan025_delay_matched_tdl import (
    read_csv_rows,
    stable_seed,
    wilson,
    write_csv_rows,
)
from tools.run_v_design_piecewise_tradeoff import decode_same_tb_batch


DESIGNS = {
    "QC": [0, 9, 18, 27, 36, 45, 54, 63],
    "Sidon": [0, 1, 3, 7, 12, 20, 30, 65],
}
SCENARIOS: Dict[str, Tuple[str, float]] = {
    "A_0p1ns": ("A", 0.1),
    "A_1ns": ("A", 1.0),
    "A_5ns": ("A", 5.0),
    "D_0p1ns": ("D", 0.1),
    "D_1ns": ("D", 1.0),
    "E_0p1ns": ("E", 0.1),
    "E_1ns": ("E", 1.0),
}
SEED = 20260726
DEFAULT_OUTPUT_ROOT = ROOT / "outputs" / "experiment026_cdd_design"
DEFAULT_RUN_ID = "20260724_main"
PRESCAN_SNRS = tuple(np.arange(13.0, 18.5 + 0.01, 0.5))


def parse_list(text: str) -> List[str]:
    return [item.strip() for item in str(text).split(",") if item.strip()]


def parse_float_list(text: str) -> List[float]:
    return [float(item.strip()) for item in str(text).split(",") if item.strip()]


def resource_config() -> ResourceConfig:
    return ResourceConfig(
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


def build_scenario(scenario_id: str):
    profile, delay_spread_ns = SCENARIOS[scenario_id]
    resource = resource_config()
    channel = ChannelConfig(
        backend="sionna_tdl",
        model="3gpp_tr38901_tdl",
        tdl_profile=profile,
        delay_spread_ns=delay_spread_ns,
        carrier_frequency_hz=3.5e9,
        ue_speed_kmh=0.0,
        num_sinusoids=20,
    )
    grid = build_resource_grid(resource)
    precoders = {
        name: build_active_dft_grid_precoder(grid, delays, n_tx=8, normalize=False)
        for name, delays in DESIGNS.items()
    }
    base_covariance = tdl_active_frequency_covariance(grid, channel)
    covariances = {
        name: base_covariance * (precoder.C @ precoder.C.conj().T)
        for name, precoder in precoders.items()
    }
    return channel, grid, precoders, covariances


def _pilot_geometry(precoder, pilot_local: np.ndarray) -> Dict[str, float]:
    singular = np.linalg.svd(precoder.C[pilot_local], compute_uv=False)
    return {
        "rank": int(np.linalg.matrix_rank(precoder.C[pilot_local], tol=1e-10)),
        "condition_number": float(singular[0] / singular[-1]),
        "minimum_singular_value": float(singular[-1]),
    }


def run_scenario(
    scenario_id: str,
    snrs_db: Sequence[float],
    trials: int,
    batch_size: int,
    output: Path,
) -> None:
    scenario_output = output / scenario_id
    scenario_output.mkdir(parents=True, exist_ok=True)
    channel, grid, precoders, covariances = build_scenario(scenario_id)
    mcs = get_mcs("nr_256qam", 8, None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(tb.cb_k_values, tb.cb_e_values, num_iter=8, llr_clip=50.0)
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    geometry = {
        name: _pilot_geometry(precoder, pilot_local)
        for name, precoder in precoders.items()
    }
    result_path = scenario_output / "sidon_qc_bler.csv"
    paired_path = scenario_output / "paired_error_counts.csv"
    rows = read_csv_rows(result_path)
    paired_rows = read_csv_rows(paired_path)
    profile, delay_spread_ns = SCENARIOS[scenario_id]

    receipt = {
        "experiment": "plan-026 E3",
        "scenario_id": scenario_id,
        "tdl_profile": profile,
        "delay_spread_ns": delay_spread_ns,
        "snrs_db": [float(value) for value in snrs_db],
        "trials_per_snr": int(trials),
        "batch_size": int(batch_size),
        "seed": SEED,
        "designs": DESIGNS,
        "phase_denominator": 576,
        "delay_grid_ns": 1e9 / (576 * 30e3),
        "receiver": "two-DMRS averaged, V-aware matched frequency LMMSE",
    }
    (scenario_output / "resolved_experiment.json").write_text(
        json.dumps(receipt, indent=2), encoding="utf-8"
    )

    for snr_db in snrs_db:
        completed = [
            row for row in rows
            if float(row["snr_db"]) == float(snr_db)
            and int(row.get("trials", 0)) == int(trials)
            and str(row["id"]) in DESIGNS
        ]
        if {str(row["id"]) for row in completed} == set(DESIGNS):
            print(f"[E3 {scenario_id}] resume skip snr={snr_db:g}", flush=True)
            continue
        rows = [row for row in rows if float(row["snr_db"]) != float(snr_db)]
        paired_rows = [
            row for row in paired_rows if float(row["snr_db"]) != float(snr_db)
        ]
        snr_linear = 10.0 ** (float(snr_db) / 10.0)
        noise_variance = 8.0 / snr_linear
        ls_noise_variance = noise_variance / 2.0
        filters = {
            name: build_frequency_rmmse_filter(
                covariance,
                pilot_local,
                ls_noise_variance,
                diagonal_loading=1e-10,
            )
            for name, covariance in covariances.items()
        }
        errors = {name: 0 for name in DESIGNS}
        nmse_sum = {name: 0.0 for name in DESIGNS}
        paired = {
            "both_error": 0,
            "qc_only_error": 0,
            "sidon_only_error": 0,
            "both_success": 0,
        }
        started = time.time()
        completed_trials = 0
        for batch_start in range(1, int(trials) + 1, int(batch_size)):
            current = min(int(batch_size), int(trials) - batch_start + 1)
            payload_rng = np.random.default_rng(
                stable_seed(SEED, scenario_id, snr_db, batch_start, "payload")
            )
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
                seed=stable_seed(SEED, scenario_id, snr_db, batch_start, "channel"),
            )
            noise_rng = np.random.default_rng(
                stable_seed(SEED, scenario_id, snr_db, batch_start, "noise")
            )
            averaged_ls_noise = math.sqrt(ls_noise_variance / 2.0) * (
                noise_rng.normal(size=(current, 1, grid.pilot_count))
                + 1j * noise_rng.normal(size=(current, 1, grid.pilot_count))
            )
            data_noise = math.sqrt(noise_variance / 2.0) * (
                noise_rng.normal(size=(current, 1, grid.n_data_re))
                + 1j * noise_rng.normal(size=(current, 1, grid.n_data_re))
            )
            decoded_by_design = {}
            for name, precoder in precoders.items():
                effective = equivalent_channel(realization.H, precoder.C)[:, 0:1]
                pilot0 = effective[:, :, int(grid.pilot_symbol_indices[0]), pilot_local]
                pilot1 = effective[:, :, int(grid.pilot_symbol_indices[-1]), pilot_local]
                true_pilot_average = 0.5 * (pilot0 + pilot1)
                estimate_full = filters[name].estimate_full_band(
                    true_pilot_average + averaged_ls_noise
                )
                estimate_data = estimate_full[:, :, data_local]
                true_data = effective[:, :, grid.data_symbol_indices, data_local]
                per_trial_nmse = np.sum(
                    np.abs(estimate_data - true_data) ** 2, axis=(1, 2)
                ) / np.maximum(np.sum(np.abs(true_data) ** 2, axis=(1, 2)), 1e-30)
                nmse_sum[name] += float(np.sum(per_trial_nmse))
                received = true_data * symbols[None, None, :] + data_noise
                denominator = np.maximum(
                    np.sum(np.abs(estimate_data) ** 2, axis=1), 1e-10
                )
                equalized = np.sum(
                    np.conj(estimate_data) * received, axis=1
                ) / denominator
                effective_noise = noise_variance / denominator
                llrs = [
                    qam_demapper_maxlog(
                        equalized[index],
                        effective_noise[index],
                        int(mcs.qm),
                    )
                    for index in range(current)
                ]
                decoded = decode_same_tb_batch(adapter, llrs, payload)
                decoded_by_design[name] = decoded
                errors[name] += sum(int(not item.tb_success) for item in decoded)
            for qc_result, sidon_result in zip(
                decoded_by_design["QC"], decoded_by_design["Sidon"]
            ):
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
                    f"[E3 {scenario_id}] snr={snr_db:g} "
                    f"trial={completed_trials}/{trials} errors={errors} "
                    f"elapsed={time.time()-started:.1f}s",
                    flush=True,
                )

        for name in DESIGNS:
            low, high = wilson(errors[name], int(trials))
            mean_nmse = nmse_sum[name] / float(trials)
            filt = filters[name]
            rows.append({
                "scenario_id": scenario_id,
                "tdl_profile": profile,
                "delay_spread_ns": delay_spread_ns,
                "id": name,
                "design": name,
                "snr_db": float(snr_db),
                "trials": int(trials),
                "tb_errors": int(errors[name]),
                "bler": float(errors[name]) / float(trials),
                "bler_wilson95_lo": low,
                "bler_wilson95_hi": high,
                "ce_nmse_mean": mean_nmse,
                "ce_nmse_mean_dB": 10.0 * math.log10(mean_nmse),
                "estimator_condition_number": filt.condition_number,
                "estimator_min_singular_value": filt.minimum_singular_value,
                "numerical_jitter": filt.numerical_jitter,
                "pilot_matrix_rank": geometry[name]["rank"],
                "pilot_matrix_condition_number": geometry[name]["condition_number"],
                "pilot_matrix_min_singular_value": geometry[name]["minimum_singular_value"],
                "noise_variance": noise_variance,
                "noise_variance_ls": ls_noise_variance,
            })
        paired_rows.append({
            "scenario_id": scenario_id,
            "tdl_profile": profile,
            "delay_spread_ns": delay_spread_ns,
            "snr_db": float(snr_db),
            "trials": int(trials),
            **paired,
        })
        rows.sort(key=lambda row: (float(row["snr_db"]), str(row["id"])))
        paired_rows.sort(key=lambda row: float(row["snr_db"]))
        write_csv_rows(rows, result_path)
        write_csv_rows(paired_rows, paired_path)


def crossing_interval(
    rows: Sequence[Dict[str, object]],
    candidate: str,
    target: float,
) -> Tuple[float, float]:
    selected = sorted(
        (row for row in rows if str(row["id"]) == candidate),
        key=lambda row: float(row["snr_db"]),
    )
    snrs = np.asarray([float(row["snr_db"]) for row in selected])
    bler = np.asarray([float(row["bler"]) for row in selected])
    monotone = np.minimum.accumulate(bler)
    for index in range(len(snrs) - 1):
        if monotone[index] >= target and monotone[index + 1] <= target:
            return float(snrs[index]), float(snrs[index + 1])
    raise ValueError(f"No {target:g} crossing for {candidate}.")


def ensure_10pct_crossings(
    scenarios: Sequence[str],
    trials: int,
    batch_size: int,
    output: Path,
) -> None:
    for scenario_id in scenarios:
        rows = read_csv_rows(output / scenario_id / "sidon_qc_bler.csv")
        extra: set[float] = set()
        for candidate in DESIGNS:
            selected = [
                row for row in rows if str(row["id"]) == candidate
            ]
            try:
                crossing_interval(selected, candidate, 0.10)
                continue
            except ValueError:
                pass
            ordered = sorted(selected, key=lambda row: float(row["snr_db"]))
            if all(float(row["bler"]) > 0.10 for row in ordered):
                extra.update((19.0, 19.5))
            elif all(float(row["bler"]) < 0.10 for row in ordered):
                extra.update((12.0, 12.5))
            else:
                raise RuntimeError(
                    f"{scenario_id}/{candidate}: non-monotone 10% crossing cannot "
                    "be resolved by the predetermined extension rule."
                )
        if extra:
            run_scenario(scenario_id, sorted(extra), trials, batch_size, output)


def generate_refinement_grid(
    scenarios: Sequence[str],
    base_output: Path,
) -> Dict[str, object]:
    source = base_output / "e3_prescan"
    result: Dict[str, object] = {
        "method": "union of QC/Sidon adjacent crossing brackets, 0.25 dB grid",
        "scenarios": {},
    }
    for scenario_id in scenarios:
        rows = read_csv_rows(source / scenario_id / "sidon_qc_bler.csv")
        scenario_result: Dict[str, object] = {}
        for label, target in (("bler10", 0.10), ("bler1", 0.01)):
            try:
                intervals = {
                    name: crossing_interval(rows, name, target)
                    for name in DESIGNS
                }
            except ValueError as error:
                scenario_result[label] = {
                    "available": False,
                    "reason": str(error),
                }
                continue
            grid_values: set[float] = set()
            for lower, upper in intervals.values():
                count = int(round((upper - lower) / 0.25)) + 1
                grid_values.update(
                    float(round(lower + 0.25 * index, 10))
                    for index in range(count)
                )
            grid = sorted(grid_values)
            scenario_result[label] = {
                "available": True,
                "probability": target,
                "candidate_crossing_intervals_db": {
                    name: list(value) for name, value in intervals.items()
                },
                "refinement_snrs_db": grid,
                "trials_per_snr": 3000,
            }
        result["scenarios"][scenario_id] = scenario_result
    (base_output / "e3_refinement_grid.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    return result


def analyze(base_output: Path, scenarios: Sequence[str]) -> Dict[str, object]:
    final = base_output / "final"
    final.mkdir(parents=True, exist_ok=True)
    summary: Dict[str, object] = {
        "experiment": "plan-026 E3",
        "scenarios": {},
    }
    combined: List[Dict[str, object]] = []
    for scenario_id in scenarios:
        scenario_summary: Dict[str, object] = {}
        for label, directory, target in (
            ("bler10", "e3_refine_10pct", 0.10),
            ("bler1", "e3_refine_1pct", 0.01),
        ):
            rows = read_csv_rows(
                base_output / directory / scenario_id / "sidon_qc_bler.csv"
            )
            if not rows:
                scenario_summary[label] = {"available": False}
                continue
            value = target_summary(rows, label, target)
            value["available"] = True
            scenario_summary[label] = value
            combined.extend(rows)
        paired_directories = [
            base_output / directory / scenario_id
            for directory in ("e3_refine_10pct", "e3_refine_1pct")
            if (base_output / directory / scenario_id).exists()
        ]
        scenario_summary["paired_error_analysis"] = (
            paired_error_analysis(paired_directories)
            if paired_directories else {}
        )
        summary["scenarios"][scenario_id] = scenario_summary
    write_csv_rows(combined, final / "e3_refined_combined.csv")
    (final / "e3_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage",
        choices=(
            "smoke",
            "prescan",
            "refine-grid",
            "refine-10pct",
            "refine-1pct",
            "analyze",
        ),
        required=True,
    )
    parser.add_argument("--run-id", default=DEFAULT_RUN_ID)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--scenarios", default=",".join(SCENARIOS))
    parser.add_argument(
        "--snrs",
        default=",".join(f"{value:g}" for value in PRESCAN_SNRS),
    )
    parser.add_argument("--trials", type=int, default=400)
    parser.add_argument("--batch-size", type=int, default=20)
    args = parser.parse_args()
    scenarios = parse_list(args.scenarios)
    unknown = sorted(set(scenarios) - set(SCENARIOS))
    if unknown:
        raise ValueError(f"Unknown scenarios: {unknown}")
    base_output = args.output_root / Path(str(args.run_id))

    if args.stage == "smoke":
        for scenario_id in scenarios:
            run_scenario(
                scenario_id, [15.0], 20, int(args.batch_size),
                base_output / "e3_smoke",
            )
        return
    if args.stage == "prescan":
        output = base_output / "e3_prescan"
        snrs = parse_float_list(args.snrs)
        for scenario_id in scenarios:
            run_scenario(
                scenario_id, snrs, int(args.trials), int(args.batch_size), output
            )
        ensure_10pct_crossings(
            scenarios, int(args.trials), int(args.batch_size), output
        )
        return
    if args.stage == "refine-grid":
        print(json.dumps(
            generate_refinement_grid(scenarios, base_output),
            indent=2,
        ))
        return
    grid = json.loads(
        (base_output / "e3_refinement_grid.json").read_text(encoding="utf-8")
    )
    if args.stage in ("refine-10pct", "refine-1pct"):
        label = "bler10" if args.stage == "refine-10pct" else "bler1"
        directory = (
            "e3_refine_10pct" if label == "bler10" else "e3_refine_1pct"
        )
        for scenario_id in scenarios:
            spec = grid["scenarios"][scenario_id][label]
            if not spec["available"]:
                print(
                    f"[E3 {scenario_id}] {label} skipped: {spec['reason']}",
                    flush=True,
                )
                continue
            run_scenario(
                scenario_id,
                spec["refinement_snrs_db"],
                int(spec["trials_per_snr"]),
                int(args.batch_size),
                base_output / directory,
            )
        return
    print(json.dumps(analyze(base_output, scenarios), indent=2))


if __name__ == "__main__":
    main()
