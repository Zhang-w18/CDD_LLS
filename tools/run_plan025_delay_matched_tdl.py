"""Run the plan-025 delay-matched TDL-A supplement to plan-024 E2."""

from __future__ import annotations

import argparse
import csv
import importlib.metadata
import json
import math
import platform
import subprocess
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
    tdl_unknown_delay_covariance,
)
from cdd_lls.phy.ldpc import SionnaLDPCAdapter
from cdd_lls.phy.precoding import build_active_dft_grid_precoder, equivalent_channel
from cdd_lls.phy.qam import qam_demapper_maxlog, qam_modulate
from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers
from tools.analyze_plan025_sionna_e2 import (
    ce_fairness,
    paired_error_analysis,
    read_rows,
    target_summary,
    write_rows as write_analysis_rows,
)
from tools.run_v_design_piecewise_tradeoff import decode_same_tb_batch


DESIGNS = {
    "QC": [0, 9, 18, 27, 36, 45, 54, 63],
    "Sidon": [0, 1, 3, 7, 12, 20, 30, 65],
}
SEED = 20260716
DEFAULT_OUTPUT_ROOT = ROOT / "outputs" / "experiment025_sionna_tdl_rmmse"
DEFAULT_BASE_RUN_ID = "20260723_delay_matched"


def parse_float_list(text: str) -> List[float]:
    return [float(item.strip()) for item in str(text).split(",") if item.strip()]


def stable_seed(*items: object) -> int:
    text = "|".join(str(item) for item in items)
    value = 2166136261
    for char in text:
        value ^= ord(char)
        value = (value * 16777619) % (2**32)
    return int(value % (2**31 - 1))


def write_csv_rows(rows: Sequence[Dict[str, object]], path: Path) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted(set().union(*(row.keys() for row in rows)))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_csv_rows(path: Path) -> List[Dict[str, object]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def wilson(errors: int, trials: int, z: float = 1.96) -> Tuple[float, float]:
    if trials <= 0:
        return float("nan"), float("nan")
    probability = float(errors) / float(trials)
    denominator = 1.0 + z * z / float(trials)
    center = (probability + z * z / (2.0 * trials)) / denominator
    half = z * math.sqrt(
        probability * (1.0 - probability) / trials + z * z / (4.0 * trials * trials)
    ) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def scenario():
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
    precoders = {
        name: build_active_dft_grid_precoder(grid, indices, n_tx=8, normalize=False)
        for name, indices in DESIGNS.items()
    }
    base_frequency_covariance = tdl_unknown_delay_covariance(grid, channel).frequency
    effective_covariances = {
        name: base_frequency_covariance * (precoder.C @ precoder.C.conj().T)
        for name, precoder in precoders.items()
    }
    return resource, channel, grid, precoders, effective_covariances


def correlation_proxy_analysis() -> Dict[str, object]:
    """Compare flat and TDL-weighted off-diagonal correlation moments."""
    _, _, grid, precoders, effective_covariances = scenario()
    n_tx = len(DESIGNS["QC"])
    sidon = np.asarray(DESIGNS["Sidon"], dtype=np.int64)
    pair_sums_mod_k = [
        int((sidon[a] + sidon[b]) % grid.n_sc)
        for a in range(n_tx)
        for b in range(a, n_tx)
    ]
    ordered_quadruples_mod_k = sum(
        int((a + c - b - d) % grid.n_sc == 0)
        for a in sidon
        for c in sidon
        for b in sidon
        for d in sidon
    )
    trivial_quadruples = 2 * n_tx * n_tx - n_tx
    rows: List[Dict[str, object]] = []
    by_design: Dict[str, object] = {}
    for name in DESIGNS:
        flat_correlation = (
            precoders[name].C @ precoders[name].C.conj().T
        ) / float(n_tx)
        tdl_correlation = effective_covariances[name] / float(n_tx)
        by_power: Dict[str, object] = {}
        for correlation_power in (2, 4):
            flat_total = float(np.sum(np.abs(flat_correlation) ** correlation_power))
            tdl_total = float(np.sum(np.abs(tdl_correlation) ** correlation_power))
            flat_offdiagonal = flat_total - float(grid.n_sc)
            tdl_offdiagonal = tdl_total - float(grid.n_sc)
            reduction_fraction = 1.0 - tdl_offdiagonal / flat_offdiagonal
            row = {
                "design": name,
                "correlation_power": correlation_power,
                "flat_offdiagonal_sum": flat_offdiagonal,
                "tdl_a_5ns_offdiagonal_sum": tdl_offdiagonal,
                "tdl_reduction_fraction": reduction_fraction,
                "tdl_reduction_percent": 100.0 * reduction_fraction,
            }
            rows.append(row)
            by_power[str(correlation_power)] = {
                key: value for key, value in row.items()
                if key not in ("design", "correlation_power")
            }
        by_design[name] = by_power
    return {
        "definition": (
            "Sum over k != l of |rho(k,l)|^p. "
            "Flat rho=(C C^H)/8; TDL rho=R_TDL hadamard (C C^H)/8."
        ),
        "n_active_subcarriers": int(grid.n_sc),
        "sidon_check": {
            "indices": sidon.tolist(),
            "unordered_pair_count_with_repetition": len(pair_sums_mod_k),
            "unique_pair_sums_mod_k": len(set(pair_sums_mod_k)),
            "ordered_quadruples_mod_k": ordered_quadruples_mod_k,
            "trivial_ordered_quadruples": trivial_quadruples,
            "nontrivial_ordered_quadruples": ordered_quadruples_mod_k - trivial_quadruples,
            "passed": bool(
                len(set(pair_sums_mod_k)) == len(pair_sums_mod_k)
                and ordered_quadruples_mod_k == trivial_quadruples
            ),
        },
        "rows": rows,
        "by_design": by_design,
    }


def environment_receipt() -> Dict[str, object]:
    packages = ("numpy", "scipy", "matplotlib", "PyYAML", "tensorflow", "sionna")
    return {
        "python": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "packages": {name: importlib.metadata.version(name) for name in packages},
    }


def save_run_receipts(
    output: Path,
    snrs_db: Sequence[float],
    trials: int,
    batch_size: int,
    grid,
    precoders,
) -> None:
    output.mkdir(parents=True, exist_ok=True)
    receipt = {
        "experiment": "plan-025 section 5 delay-matched TDL-A supplement",
        "only_active_change_vs_plan024": "flat branch channel -> Sionna TDL-A 5 ns",
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
            "n_active_subcarriers": int(grid.n_sc),
            "n_symbols": int(grid.n_symbols),
            "dmrs_symbol_indices": [2, 7],
            "dmrs_spacing_sc": 24,
            "dmrs_offset_sc": 0,
            "n_dmrs_re": int(grid.n_dmrs_re),
            "n_unique_frequency_pilots": int(grid.pilot_count),
            "n_data_re": int(grid.n_data_re),
            "n_fft": int(grid.n_fft),
            "scs_hz": float(grid.scs_khz) * 1e3,
            "cyclic_prefix_length_samples": int(grid.cyclic_prefix_length),
        },
        "designs": {name: precoder.metadata for name, precoder in precoders.items()},
        "receiver": "plan-024-compatible two-DMRS average plus known-delay matched full-band frequency LMMSE",
        "precoder_normalization": "unnormalized abs(V)=1",
        "data_noise_variance": "N0=8/SNR",
        "averaged_ls_noise_variance": "N0/2=4/SNR",
        "mcs": {"modulation": "16QAM", "index": 8, "code_rate": 553 / 1024},
        "ldpc_iterations": 8,
        "llr_noise": "N0 only; no CE-error-aware term",
        "snrs_db": [float(value) for value in snrs_db],
        "trials_per_snr": int(trials),
        "batch_size": int(batch_size),
        "seed": SEED,
        "common_random_numbers": ["TDL realization", "payload", "averaged LS noise", "data noise"],
    }
    (output / "resolved_experiment.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    environment = environment_receipt()
    environment["command"] = [sys.executable, *sys.argv]
    (output / "environment.json").write_text(
        json.dumps(environment, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (output / "commands.json").write_text(
        json.dumps({"command": [sys.executable, *sys.argv]}, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def run_preflight(base_output: Path) -> Dict[str, object]:
    validation = base_output / "validation"
    validation.mkdir(parents=True, exist_ok=True)
    _, _, grid, precoders, covariances = scenario()
    expected_local = np.arange(576, dtype=np.float64)
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    checks: Dict[str, object] = {}
    passed = True
    for name, indices in DESIGNS.items():
        expected = np.exp(
            -2j * np.pi * expected_local[:, None]
            * np.asarray(indices, dtype=np.float64)[None, :] / 576.0
        )
        precoder = precoders[name]
        phase_error = float(np.max(np.abs(precoder.C - expected)))
        sample_expected = np.asarray(indices, dtype=np.float64) * 4096.0 / 576.0
        seconds_expected = np.asarray(indices, dtype=np.float64) / (576.0 * 30e3)
        sample_error = float(np.max(np.abs(
            np.asarray(precoder.metadata["cdd_delay_fft_samples"]) - sample_expected
        )))
        seconds_error = float(np.max(np.abs(
            np.asarray(precoder.metadata["cdd_delay_seconds"]) - seconds_expected
        )))
        covariance = covariances[name]
        hermitian_error = float(np.max(np.abs(covariance - covariance.conj().T)))
        minimum_eigenvalue = float(np.min(np.linalg.eigvalsh((covariance + covariance.conj().T) / 2.0)))
        diagonal_error = float(np.max(np.abs(np.diag(covariance) - 8.0)))
        pilot_matrix = precoder.C[pilot_local]
        pilot_singular = np.linalg.svd(pilot_matrix, compute_uv=False)
        design_passed = bool(
            phase_error <= 1e-12
            and sample_error <= 1e-12
            and seconds_error <= 1e-18
            and hermitian_error <= 1e-10
            and minimum_eigenvalue >= -1e-9
            and diagonal_error <= 1e-10
            and np.max(np.abs(np.abs(precoder.C) - 1.0)) <= 1e-14
            and np.all(np.isfinite(covariance))
        )
        checks[name] = {
            "phase_max_abs_error": phase_error,
            "fft_sample_conversion_max_abs_error": sample_error,
            "seconds_conversion_max_abs_error": seconds_error,
            "covariance_hermitian_max_abs_error": hermitian_error,
            "covariance_minimum_eigenvalue": minimum_eigenvalue,
            "covariance_diagonal_max_abs_error_vs_8": diagonal_error,
            "pilot_matrix_rank": int(np.linalg.matrix_rank(pilot_matrix, tol=1e-10)),
            "pilot_matrix_condition_number": float(pilot_singular[0] / pilot_singular[-1]),
            "pilot_matrix_minimum_singular_value": float(pilot_singular[-1]),
            "passed": design_passed,
        }
        passed = passed and design_passed
    resource_check = {
        "n_active_subcarriers": int(grid.n_sc),
        "pilot_re": int(grid.n_dmrs_re),
        "unique_frequency_pilots": int(grid.pilot_count),
        "data_re": int(grid.n_data_re),
        "pilot_local_indices": [int(x) for x in pilot_local],
        "passed": bool(
            grid.n_sc == 576 and grid.n_dmrs_re == 48 and grid.pilot_count == 24
            and grid.n_data_re == 5712
            and np.array_equal(pilot_local, np.arange(0, 576, 24, dtype=np.int64))
        ),
    }
    passed = passed and bool(resource_check["passed"])
    summary = {
        "definition": {
            "K": 576,
            "N_FFT": 4096,
            "subcarrier_spacing_hz": 30000.0,
            "data_noise_variance": "8/SNR",
            "averaged_ls_noise_variance": "4/SNR",
        },
        "design_checks": checks,
        "resource_check": resource_check,
        "passed": bool(passed),
    }
    (validation / "preflight.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (validation / "environment.json").write_text(
        json.dumps(environment_receipt(), indent=2, ensure_ascii=False), encoding="utf-8"
    )
    if not passed:
        raise RuntimeError(f"Preflight failed; see {validation / 'preflight.json'}")
    return summary


def run_tests(base_output: Path) -> None:
    validation = base_output / "validation"
    command = [sys.executable, "-m", "pytest", "tests", "-q"]
    completed = subprocess.run(
        command,
        cwd=ROOT,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    log = completed.stdout + ("\n" + completed.stderr if completed.stderr else "")
    (validation / "tests.log").write_text(log, encoding="utf-8")
    (validation / "test_command.json").write_text(
        json.dumps({"command": command, "exit_code": completed.returncode}, indent=2), encoding="utf-8"
    )
    if completed.returncode != 0:
        raise RuntimeError(f"Tests failed; see {validation / 'tests.log'}")


def run_comparison(
    snrs_db: Sequence[float],
    trials: int,
    batch_size: int,
    output: Path,
) -> Path:
    output.mkdir(parents=True, exist_ok=True)
    _, channel, grid, precoders, covariances = scenario()
    mcs = get_mcs("nr_256qam", 8, None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(tb.cb_k_values, tb.cb_e_values, num_iter=8, llr_clip=50.0)
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    pilot_geometry = {}
    for name, precoder in precoders.items():
        singular = np.linalg.svd(precoder.C[pilot_local], compute_uv=False)
        pilot_geometry[name] = {
            "rank": int(np.linalg.matrix_rank(precoder.C[pilot_local], tol=1e-10)),
            "condition_number": float(singular[0] / singular[-1]),
            "minimum_singular_value": float(singular[-1]),
        }
    save_run_receipts(output, snrs_db, trials, batch_size, grid, precoders)
    result_path = output / "sidon_qc_bler.csv"
    paired_path = output / "paired_error_counts.csv"
    rows = read_csv_rows(result_path)
    paired_rows = read_csv_rows(paired_path)
    if result_path.exists():
        header = result_path.read_text(encoding="utf-8").splitlines()[0].split(",")
        if len({field.lower() for field in header}) != len(header):
            rows = []

    for snr_db in snrs_db:
        completed = [
            row for row in rows
            if float(row["snr_db"]) == float(snr_db)
            and int(row.get("trials", 0)) == int(trials)
            and str(row["id"]) in DESIGNS
        ]
        if {str(row["id"]) for row in completed} == set(DESIGNS):
            print(f"[delay-matched] resume skip snr={snr_db:g}", flush=True)
            continue
        rows = [row for row in rows if float(row["snr_db"]) != float(snr_db)]
        paired_rows = [row for row in paired_rows if float(row["snr_db"]) != float(snr_db)]
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
                g = equivalent_channel(realization.H, precoder.C)[:, 0:1]
                pilot_symbol_0 = g[:, :, int(grid.pilot_symbol_indices[0]), pilot_local]
                pilot_symbol_1 = g[:, :, int(grid.pilot_symbol_indices[-1]), pilot_local]
                true_pilot_average = 0.5 * (pilot_symbol_0 + pilot_symbol_1)
                estimate_full = filters[name].estimate_full_band(true_pilot_average + averaged_ls_noise)
                estimate_data = estimate_full[:, :, data_local]
                true_data = g[:, :, grid.data_symbol_indices, data_local]
                per_trial_nmse = np.sum(np.abs(estimate_data - true_data) ** 2, axis=(1, 2)) / np.maximum(
                    np.sum(np.abs(true_data) ** 2, axis=(1, 2)), 1e-30
                )
                nmse_sum[name] += float(np.sum(per_trial_nmse))
                y = true_data * symbols[None, None, :] + data_noise
                denominator = np.maximum(np.sum(np.abs(estimate_data) ** 2, axis=1), 1e-10)
                equalized = np.sum(np.conj(estimate_data) * y, axis=1) / denominator
                effective_noise = noise_variance / denominator
                llrs = [
                    qam_demapper_maxlog(equalized[index], effective_noise[index], int(mcs.qm))
                    for index in range(current)
                ]
                decoded = decode_same_tb_batch(adapter, llrs, payload)
                decoded_by_design[name] = decoded
                errors[name] += sum(int(not result.tb_success) for result in decoded)
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
                    f"[delay-matched] snr={snr_db:g} trial={completed_trials}/{trials} "
                    f"errors={errors} elapsed={time.time()-started:.1f}s",
                    flush=True,
                )

        for name in DESIGNS:
            low, high = wilson(errors[name], int(trials))
            filt = filters[name]
            mean_nmse = nmse_sum[name] / float(trials)
            rows.append({
                "id": name,
                "design": name,
                "variant_id": f"{name}_known",
                "branch": "known",
                "snr_db": float(snr_db),
                "trials": int(trials),
                "n_trials": int(trials),
                "tb_errors": int(errors[name]),
                "bler": float(errors[name]) / float(trials),
                "bler_wilson95_lo": low,
                "bler_wilson95_hi": high,
                "ce_nmse_mean": mean_nmse,
                "ce_nmse_mean_dB": 10.0 * math.log10(mean_nmse),
                "condition_number": filt.condition_number,
                "estimator_condition_number": filt.condition_number,
                "estimator_min_singular_value": filt.minimum_singular_value,
                "numerical_jitter": filt.numerical_jitter,
                "pilot_matrix_rank": pilot_geometry[name]["rank"],
                "pilot_matrix_condition_number": pilot_geometry[name]["condition_number"],
                "pilot_matrix_min_singular_value": pilot_geometry[name]["minimum_singular_value"],
                "covariance_type": "tdl_a_5ns_times_unnormalized_active_dft_grid_cdd",
                "noise_variance": noise_variance,
                "noise_variance_ls": ls_noise_variance,
            })
        paired_rows.append({"snr_db": float(snr_db), "trials": int(trials), **paired})
        rows.sort(key=lambda row: (float(row["snr_db"]), str(row["id"])))
        paired_rows.sort(key=lambda row: float(row["snr_db"]))
        write_csv_rows(rows, result_path)
        write_csv_rows(rows, output / "summary.csv")
        write_csv_rows(paired_rows, paired_path)
    return output


def crossing_interval(rows: Sequence[Dict[str, object]], candidate: str, target: float) -> Tuple[float, float]:
    selected = sorted(
        (row for row in rows if str(row["id"]) == candidate),
        key=lambda row: float(row["snr_db"]),
    )
    snrs = np.asarray([float(row["snr_db"]) for row in selected], dtype=np.float64)
    bler = np.asarray([float(row["bler"]) for row in selected], dtype=np.float64)
    monotone = np.minimum.accumulate(bler)
    for index in range(len(snrs) - 1):
        if monotone[index] >= target and monotone[index + 1] <= target:
            return float(snrs[index]), float(snrs[index + 1])
    raise ValueError(f"No {target:g} crossing for {candidate} in prescan.")


def generate_refinement_grid(base_output: Path) -> Dict[str, object]:
    rows = read_csv_rows(base_output / "prescan" / "sidon_qc_bler.csv")
    if not rows:
        raise ValueError("Prescan output is missing.")
    result: Dict[str, object] = {
        "method": "plan-024 two-stage scan; union of monotone 0.5 dB crossing brackets at 0.25 dB spacing",
        "source": str(base_output / "prescan" / "sidon_qc_bler.csv"),
        "targets": {},
    }
    for label, target in (("bler10", 0.10), ("bler1", 0.01)):
        intervals = {name: crossing_interval(rows, name, target) for name in DESIGNS}
        lower = min(interval[0] for interval in intervals.values())
        upper = max(interval[1] for interval in intervals.values())
        count = int(round((upper - lower) / 0.25)) + 1
        grid = [float(round(lower + 0.25 * index, 10)) for index in range(count)]
        result["targets"][label] = {
            "probability": target,
            "candidate_crossing_intervals_db": {name: list(interval) for name, interval in intervals.items()},
            "refinement_snrs_db": grid,
            "trials_per_snr": 3000,
        }
    (base_output / "refinement_grid.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return result


def plot_results(refine10, refine1, summary, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"QC": "#f77f00", "Sidon": "#111111"}
    markers = {"QC": "o", "Sidon": "s"}
    figure, axes = plt.subplots(1, 2, figsize=(12.4, 5.2))
    for axis, rows, label, target in zip(
        axes, (refine10, refine1), ("bler10", "bler1"), (0.10, 0.01)
    ):
        for candidate in DESIGNS:
            selected = sorted(
                (row for row in rows if row["id"] == candidate),
                key=lambda row: float(row["snr_db"]),
            )
            x = np.asarray([float(row["snr_db"]) for row in selected])
            y = np.asarray([max(float(row["bler"]), 0.5 / int(row["trials"])) for row in selected])
            low = np.asarray([max(float(row["bler_wilson95_lo"]), 1e-5) for row in selected])
            high = np.asarray([float(row["bler_wilson95_hi"]) for row in selected])
            fit = summary["targets"][label]["candidates"][candidate]
            xx = np.linspace(float(np.min(x)), float(np.max(x)), 200)
            probability = 1.0 / (
                1.0 + np.exp(-(float(fit["logit_intercept"]) + float(fit["logit_slope_per_db"]) * xx))
            )
            axis.semilogy(x, y, marker=markers[candidate], ls="none", color=colors[candidate], label=candidate)
            axis.vlines(x, low, high, color=colors[candidate], alpha=0.45)
            axis.semilogy(xx, probability, color=colors[candidate])
        target_result = summary["targets"][label]
        axis.axhline(target, color="#555555", ls="--", lw=1)
        axis.set_title(
            f"{target*100:g}% BLER: Sidon gain {target_result['sidon_gain_vs_qc_db']:.2f} dB\n"
            f"95% conservative CI [{target_result['sidon_gain_ci95_conservative_lo_db']:.2f}, "
            f"{target_result['sidon_gain_ci95_conservative_hi_db']:.2f}] dB"
        )
        axis.set_xlabel("average receive SNR (dB)")
        axis.set_ylabel("estimated-CSI TB BLER")
        axis.grid(alpha=0.3, which="both")
        axis.legend()
    figure.suptitle("Plan-025 supplement: plan-024 delays on Sionna TDL-A 5 ns")
    figure.tight_layout()
    figure.savefig(path, dpi=170)
    plt.close(figure)


def analyze_final(base_output: Path) -> Dict[str, object]:
    prescan = read_rows(base_output / "prescan")
    refine10 = read_rows(base_output / "refine_10pct")
    refine1 = read_rows(base_output / "refine_1pct")
    final = base_output / "final"
    final.mkdir(parents=True, exist_ok=True)
    write_analysis_rows([*refine10, *refine1], final / "sidon_qc_refined_combined.csv")
    targets = {
        "bler10": target_summary(refine10, "bler10", 0.10),
        "bler1": target_summary(refine1, "bler1", 0.01),
    }
    for label in targets:
        targets[label]["first_execution_4096_gain_db"] = -0.161 if label == "bler10" else -0.582
        targets[label]["gain_difference_vs_first_execution_db"] = (
            float(targets[label]["sidon_gain_vs_qc_db"])
            - float(targets[label]["first_execution_4096_gain_db"])
        )
    pilot_geometry = {}
    for row in refine10:
        pilot_geometry[row["id"]] = {
            "rank": int(row["pilot_matrix_rank"]),
            "condition_number": float(row["pilot_matrix_condition_number"]),
            "minimum_singular_value": float(row["pilot_matrix_min_singular_value"]),
        }
    correlation_proxy = correlation_proxy_analysis()
    write_csv_rows(correlation_proxy["rows"], final / "correlation_proxy.csv")
    summary = {
        "model": "Sionna 1.0.2 TDL-A, 5 ns, 0 km/h, 3.5 GHz",
        "receiver": "plan-024-compatible averaged-DMRS known-delay matched frequency LMMSE",
        "precoder": "plan-024 active-band DFT-grid V, unnormalized, physical delays preserved",
        "targets": targets,
        "ce_fairness": ce_fairness([*refine10, *refine1]),
        "pilot_geometry": pilot_geometry,
        "correlation_proxy": correlation_proxy,
        "paired_error_analysis": paired_error_analysis([
            base_output / "refine_10pct", base_output / "refine_1pct"
        ]),
        "prescan": {
            "snrs_db": sorted({float(row["snr_db"]) for row in prescan}),
            "trials_per_point": sorted({int(row["trials"]) for row in prescan}),
        },
        "refinement_grid": json.loads((base_output / "refinement_grid.json").read_text(encoding="utf-8")),
        "evidence_paths": {
            "base": str(base_output),
            "prescan": str(base_output / "prescan"),
            "refine10": str(base_output / "refine_10pct"),
            "refine1": str(base_output / "refine_1pct"),
        },
    }
    (final / "final_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    plot_results(refine10, refine1, summary, final / "sidon_qc_delay_matched_tdl_bler.png")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("validate", "smoke", "link", "refine-grid", "refine", "analyze"), required=True)
    parser.add_argument("--snrs", default="13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18")
    parser.add_argument("--trials", type=int, default=400)
    parser.add_argument("--batch-size", type=int, default=20)
    parser.add_argument("--run-id", default=DEFAULT_BASE_RUN_ID)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    args = parser.parse_args()
    base_output = args.output_root / DEFAULT_BASE_RUN_ID

    if args.stage == "validate":
        base_output = args.output_root / Path(str(args.run_id))
        summary = run_preflight(base_output)
        run_tests(base_output)
        print(json.dumps(summary, indent=2, ensure_ascii=False))
        return
    if args.stage == "smoke":
        base_output = args.output_root / Path(str(args.run_id))
        output = run_comparison(
            parse_float_list(args.snrs), int(args.trials), int(args.batch_size), base_output / "smoke"
        )
        print(f"Results: {output}")
        return
    if args.stage == "link":
        output = run_comparison(
            parse_float_list(args.snrs), int(args.trials), int(args.batch_size),
            args.output_root / Path(str(args.run_id)),
        )
        print(f"Results: {output}")
        return
    if args.stage == "refine-grid":
        base_output = args.output_root / Path(str(args.run_id))
        print(json.dumps(generate_refinement_grid(base_output), indent=2, ensure_ascii=False))
        return
    if args.stage == "refine":
        base_output = args.output_root / Path(str(args.run_id))
        grid = json.loads((base_output / "refinement_grid.json").read_text(encoding="utf-8"))
        run_comparison(
            grid["targets"]["bler10"]["refinement_snrs_db"], 3000, int(args.batch_size),
            base_output / "refine_10pct",
        )
        run_comparison(
            grid["targets"]["bler1"]["refinement_snrs_db"], 3000, int(args.batch_size),
            base_output / "refine_1pct",
        )
        return
    base_output = args.output_root / Path(str(args.run_id))
    print(json.dumps(analyze_final(base_output), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
