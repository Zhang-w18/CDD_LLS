"""Validate Sionna TDL covariance and 2D RMMSE against Monte Carlo samples."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ChannelConfig, ResourceConfig, TransmissionConfig
from cdd_lls.phy.channel_tdl import generate_sionna_tdl_channel
from cdd_lls.phy.estimators import (
    build_time_frequency_rmmse_filter,
    linear_estimator_closed_form_nmse,
    tdl_known_delay_covariance,
    tdl_unknown_delay_covariance,
)
from cdd_lls.phy.precoding import build_precoder, equivalent_channel
from cdd_lls.phy.resource_grid import build_resource_grid, local_indices_for_subcarriers


DESIGNS = {
    "zero": [0] * 8,
    "QC": [0, 9, 18, 27, 36, 45, 54, 63],
    "Sidon": [0, 1, 3, 7, 12, 20, 30, 65],
}


def write_csv(rows: List[Dict[str, object]], path: Path) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted(set().union(*(row.keys() for row in rows)))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def stable_seed(*items: object) -> int:
    text = "|".join(str(x) for x in items)
    acc = 2166136261
    for char in text:
        acc ^= ord(char)
        acc = (acc * 16777619) % (2**32)
    return int(acc % (2**31 - 1))


def parse_float_list(text: str) -> List[float]:
    return [float(value.strip()) for value in text.split(",") if value.strip()]


def scenario(n_prbs: int, speed_kmh: float):
    resource = ResourceConfig(
        n_fft=4096,
        scs_khz=30,
        n_prbs=int(n_prbs),
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
        delay_spread_ns=100.0,
        carrier_frequency_hz=3.5e9,
        ue_speed_kmh=float(speed_kmh),
        num_sinusoids=20,
    )
    return resource, channel, build_resource_grid(resource)


def representative_pairs(n_sc: int, n_symbols: int) -> List[tuple[str, int, int, int, int]]:
    center = n_sc // 2
    rows: List[tuple[str, int, int, int, int]] = []
    for offset in sorted(set((1, 6, 24, min(72, n_sc // 2 - 1)))):
        rows.append((f"freq_lag_{offset}", 0, center, 0, center + offset))
    for lag in (1, 5, n_symbols - 1):
        rows.append((f"time_lag_{lag}", 0, center, lag, center))
    rows.append(("mixed_lag_t5_f24", 0, center, min(5, n_symbols - 1), center + 24))
    return rows


def validate_covariance(samples: int, batch_size: int, output: Path) -> Dict[str, object]:
    rows: List[Dict[str, object]] = []
    structural: List[Dict[str, object]] = []
    for n_prbs in (8, 48):
        for speed in (3.0, 60.0):
            resource, channel, grid = scenario(n_prbs, speed)
            pairs = representative_pairs(grid.n_sc, grid.n_symbols)
            sums = {name: np.zeros(len(pairs), dtype=np.complex128) for name in DESIGNS}
            count = 0
            for start in range(0, int(samples), int(batch_size)):
                current = min(int(batch_size), int(samples) - start)
                realization = generate_sionna_tdl_channel(
                    grid, channel, n_tx=8, n_rx=1, batch_size=current,
                    seed=stable_seed("covariance", n_prbs, speed, start),
                )
                for name, delays in DESIGNS.items():
                    precoder = build_precoder(
                        grid, resource, TransmissionConfig(cdd_delay_vector=delays), n_tx=8,
                    )
                    g = equivalent_channel(realization.H, precoder.C)[:, 0]
                    for index, (_, sa, ka, sb, kb) in enumerate(pairs):
                        sums[name][index] += np.sum(g[:, sa, ka] * np.conj(g[:, sb, kb]))
                count += current

            for name, delays in DESIGNS.items():
                covariance = tdl_known_delay_covariance(grid, channel, delays)
                hermitian_error = float(np.max(np.abs(covariance.frequency - covariance.frequency.conj().T)))
                min_eigenvalue = float(np.min(np.linalg.eigvalsh(covariance.frequency)))
                diagonal_error = float(np.max(np.abs(np.diag(covariance.frequency) - 1.0)))
                structural.append({
                    "n_prbs": n_prbs,
                    "speed_kmh": speed,
                    "design": name,
                    "hermitian_max_error": hermitian_error,
                    "min_eigenvalue": min_eigenvalue,
                    "diagonal_max_error": diagonal_error,
                })
                errors = []
                for index, (label, sa, ka, sb, kb) in enumerate(pairs):
                    empirical = sums[name][index] / float(count)
                    theoretical = covariance.time[sa, sb] * covariance.frequency[ka, kb]
                    error = abs(empirical - theoretical)
                    errors.append(float(error))
                    rows.append({
                        "n_prbs": n_prbs,
                        "speed_kmh": speed,
                        "design": name,
                        "pair": label,
                        "samples": count,
                        "theory_real": float(np.real(theoretical)),
                        "theory_imag": float(np.imag(theoretical)),
                        "empirical_real": float(np.real(empirical)),
                        "empirical_imag": float(np.imag(empirical)),
                        "absolute_error": float(error),
                    })

    write_csv(rows, output / "covariance_representative.csv")
    write_csv(structural, output / "covariance_structure.csv")
    all_errors = np.asarray([float(row["absolute_error"]) for row in rows])
    rmse = float(np.sqrt(np.mean(all_errors**2)))
    max_error = float(np.max(all_errors))
    enough = int(samples) >= 1000
    summary = {
        "samples": int(samples),
        "representative_element_rmse": rmse,
        "representative_element_max_absolute_error": max_error,
        "threshold_rmse": 0.03,
        "threshold_max_absolute_error": 0.06,
        "sample_sufficiency": "formal" if enough else "smoke_only",
        "passed": bool(rmse <= 0.03 and max_error <= 0.06) if enough else None,
        "max_hermitian_error": max(float(row["hermitian_max_error"]) for row in structural),
        "min_eigenvalue": min(float(row["min_eigenvalue"]) for row in structural),
        "max_diagonal_error": max(float(row["diagonal_max_error"]) for row in structural),
    }
    (output / "covariance_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return summary


def validate_rmmse(
    samples: int,
    batch_size: int,
    snrs_db: Sequence[float],
    output: Path,
) -> Dict[str, object]:
    rows: List[Dict[str, object]] = []
    for n_prbs in (8, 48):
        for speed in (3.0, 60.0):
            resource, channel, grid = scenario(n_prbs, speed)
            pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarrier_indices)
            data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
            design_state: Dict[str, Dict[str, object]] = {}
            for design in ("QC", "Sidon"):
                delays = DESIGNS[design]
                true_covariance = tdl_known_delay_covariance(grid, channel, delays)
                unknown_covariance = tdl_unknown_delay_covariance(grid, channel)
                precoder = build_precoder(
                    grid, resource, TransmissionConfig(cdd_delay_vector=delays), n_tx=8,
                )
                by_snr: Dict[float, Dict[str, object]] = {}
                for snr_db in snrs_db:
                    noise_variance = 10.0 ** (-float(snr_db) / 10.0)
                    filters = {
                        "known": build_time_frequency_rmmse_filter(
                            grid, true_covariance, noise_variance, 1e-12,
                        ),
                        "unknown": build_time_frequency_rmmse_filter(
                            grid, unknown_covariance, noise_variance, 1e-12,
                        ),
                    }
                    closed = {
                        branch: linear_estimator_closed_form_nmse(
                            grid, estimator, true_covariance, noise_variance,
                        )
                        for branch, estimator in filters.items()
                    }
                    by_snr[float(snr_db)] = {
                        "noise_variance": noise_variance,
                        "filters": filters,
                        "closed": closed,
                        "error_power": {branch: 0.0 for branch in filters},
                    }
                design_state[design] = {
                    "precoder": precoder,
                    "by_snr": by_snr,
                    "signal_power": 0.0,
                }

            count = 0
            for start in range(0, int(samples), int(batch_size)):
                current = min(int(batch_size), int(samples) - start)
                realization = generate_sionna_tdl_channel(
                    grid, channel, n_tx=8, n_rx=1, batch_size=current,
                    seed=stable_seed("rmmse", n_prbs, speed, start),
                )
                for design, state in design_state.items():
                    precoder = state["precoder"]
                    by_snr = state["by_snr"]
                    g = equivalent_channel(realization.H, precoder.C)[:, 0]
                    true_pilot = g[:, grid.pilot_symbol_indices, pilot_local]
                    true_data = g[:, grid.data_symbol_indices, data_local]
                    state["signal_power"] = float(state["signal_power"]) + float(
                        np.sum(np.abs(true_data) ** 2)
                    )
                    for snr_db, snr_state in by_snr.items():
                        noise_variance = float(snr_state["noise_variance"])
                        filters = snr_state["filters"]
                        error_power = snr_state["error_power"]
                        rng = np.random.default_rng(
                            stable_seed("rmmse_noise", n_prbs, speed, snr_db, start)
                        )
                        noise = math.sqrt(noise_variance / 2.0) * (
                            rng.normal(size=true_pilot.shape) + 1j * rng.normal(size=true_pilot.shape)
                        )
                        observations = true_pilot + noise
                        for branch, estimator in filters.items():
                            estimate = estimator.estimate_data(observations)
                            error_power[branch] += float(np.sum(np.abs(estimate - true_data) ** 2))
                count += current

            for design, state in design_state.items():
                signal_power = float(state["signal_power"])
                for snr_db, snr_state in state["by_snr"].items():
                    noise_variance = float(snr_state["noise_variance"])
                    filters = snr_state["filters"]
                    closed = snr_state["closed"]
                    error_power = snr_state["error_power"]
                    for branch, estimator in filters.items():
                        empirical = error_power[branch] / signal_power
                        difference_db = 10.0 * math.log10(empirical / closed[branch])
                        rows.append({
                            "n_prbs": n_prbs,
                            "speed_kmh": speed,
                            "design": design,
                            "snr_db": float(snr_db),
                            "branch": branch,
                            "samples": count,
                            "empirical_nmse": empirical,
                            "empirical_nmse_db": 10.0 * math.log10(empirical),
                            "closed_form_nmse": closed[branch],
                            "closed_form_nmse_db": 10.0 * math.log10(closed[branch]),
                            "empirical_minus_closed_db": difference_db,
                            "condition_number": estimator.condition_number,
                            "numerical_jitter": estimator.numerical_jitter,
                            "covariance_type": estimator.covariance_type,
                        })

    write_csv(rows, output / "rmmse_monte_carlo.csv")
    max_difference = max(abs(float(row["empirical_minus_closed_db"])) for row in rows)
    ordering_ok = True
    keys = sorted({
        (row["n_prbs"], row["speed_kmh"], row["design"], row["snr_db"])
        for row in rows
    })
    for key in keys:
        subset = {str(row["branch"]): float(row["closed_form_nmse"]) for row in rows
                  if (row["n_prbs"], row["speed_kmh"], row["design"], row["snr_db"]) == key}
        ordering_ok = ordering_ok and subset["known"] <= subset["unknown"] + 1e-12
    enough = int(samples) >= 1000
    summary = {
        "samples": int(samples),
        "snrs_db": [float(x) for x in snrs_db],
        "max_absolute_empirical_minus_closed_db": max_difference,
        "threshold_db": 0.2,
        "known_expected_nmse_not_worse": bool(ordering_ok),
        "sample_sufficiency": "formal" if enough else "smoke_only",
        "passed": bool(max_difference <= 0.2 and ordering_ok) if enough else None,
    }
    (output / "rmmse_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("covariance", "rmmse", "all"), default="all")
    parser.add_argument("--samples", type=int, default=4000)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--snrs", default="0,10,20")
    parser.add_argument(
        "--out", type=Path,
        default=ROOT / "outputs" / "experiment025_sionna_tdl_rmmse" / "20260722_main" / "validation",
    )
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    receipt = {
        "tdl_profile": "A",
        "delay_spread_ns": 100.0,
        "carrier_frequency_hz": 3.5e9,
        "speeds_kmh": [3.0, 60.0],
        "n_prbs": [8, 48],
        "n_fft": 4096,
        "scs_hz": 30000.0,
        "n_symbols": 10,
        "dmrs_symbol_indices": [2, 7],
        "dmrs_spacing_sc": 24,
        "samples": int(args.samples),
        "batch_size": int(args.batch_size),
        "snrs_db": parse_float_list(args.snrs),
        "seed_scheme": "stable FNV-1a labels",
    }
    (args.out / "validation_config.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    result: Dict[str, object] = {}
    if args.mode in ("covariance", "all"):
        result["covariance"] = validate_covariance(args.samples, args.batch_size, args.out)
        print(json.dumps(result["covariance"], indent=2), flush=True)
    if args.mode in ("rmmse", "all"):
        result["rmmse"] = validate_rmmse(
            args.samples, args.batch_size, parse_float_list(args.snrs), args.out,
        )
        print(json.dumps(result["rmmse"], indent=2), flush=True)


if __name__ == "__main__":
    main()
