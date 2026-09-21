"""Deterministic nonuniform small-delay CDD search for plan-028 section 13."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.design.cdd_metrics import array_factor
from cdd_lls.phy.estimators import (
    build_frequency_rmmse_filter,
    tdl_active_frequency_covariance,
)
from tools import run_bler_curves as runner
from tools import run_plan027_bler as bler027


def build_lag_nmse_metric(
    base_covariance: np.ndarray,
    weights: np.ndarray,
    pilots: np.ndarray,
    data: np.ndarray,
    noise_variance: float,
) -> dict[str, np.ndarray | float | int]:
    period = base_covariance.shape[0]
    weights = np.asarray(weights, dtype=np.complex128)[data]
    dp_difference = data[:, None] - pilots[None, :]
    pp_difference = pilots[:, None] - pilots[None, :]
    lag_count = 2 * period - 1
    cross_coefficients = np.zeros(lag_count, dtype=np.complex128)
    np.add.at(
        cross_coefficients,
        (dp_difference + period - 1).ravel(),
        (weights * base_covariance[np.ix_(data, pilots)].conj()).ravel(),
    )
    quadratic = weights.T @ weights.conj()
    filter_coefficients = np.zeros(lag_count, dtype=np.complex128)
    np.add.at(
        filter_coefficients,
        (pp_difference + period - 1).ravel(),
        (quadratic * base_covariance[np.ix_(pilots, pilots)]).ravel(),
    )
    signal_trace = float(np.real(np.sum(8.0 * base_covariance[data, data])))
    return {
        "period": period,
        "cross_coefficients": cross_coefficients,
        "filter_coefficients": filter_coefficients,
        "signal_trace": signal_trace,
        "filtered_noise_trace": float(noise_variance) * float(np.sum(np.abs(weights) ** 2)),
    }


def evaluate_lag_nmse(metric: dict, factor: np.ndarray) -> float:
    cross_trace = np.sum(metric["cross_coefficients"] * factor.conj())
    filtered_trace = (
        np.sum(metric["filter_coefficients"] * factor)
        + float(metric["filtered_noise_trace"])
    )
    signal_trace = float(metric["signal_trace"])
    numerator = max(
        float(
            np.real(
                signal_trace
                - cross_trace
                - cross_trace.conjugate()
                + filtered_trace
            )
        ),
        0.0,
    )
    return numerator / signal_trace


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--count", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260727)
    parser.add_argument("--threshold-db", type=float, default=-16.5)
    parser.add_argument("--span-min", type=float, default=4.9)
    parser.add_argument("--span-max", type=float, default=10.5)
    parser.add_argument("--output-name", default="nonuniform_search.json")
    args = parser.parse_args()

    config = runner.load_runner_config(args.config.resolve())
    runner.validate_config(config)
    scene = config["scenes"][0]
    manifest, source_digest = runner._load_manifest(scene)
    seed_candidate = runner._selected_candidates(scene, manifest)[0]
    channel, grid, _, _ = bler027.build_scene(
        manifest, [seed_candidate], str(scene["scenario_id"])
    )
    base = tdl_active_frequency_covariance(grid, channel)
    assumed = 8.0 * base
    pilots = runner.local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    data = runner.local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    snrs = (14.0, 16.0, 20.0)
    filters = {}
    noises = {}
    lag_metrics = {}
    for snr in snrs:
        noise = 4.0 / (10.0 ** (snr / 10.0))
        noises[snr] = noise
        filters[snr] = runner.build_frequency_rmmse_filter(
            assumed, pilots, noise, diagonal_loading=1e-10
        )
        lag_metrics[snr] = build_lag_nmse_metric(
            base, filters[snr].weights, pilots, data, noise
        )

    rng = np.random.default_rng(int(args.seed))
    offsets = np.arange(-3.0, 2.0 + 1e-12, 0.25)
    physical = base[0, :] / float(np.real(base[0, 0]))
    lag_weights = 2.0 * (grid.n_sc - np.arange(1, grid.n_sc, dtype=np.float64))
    signed_lags = np.arange(-grid.n_sc + 1, grid.n_sc, dtype=np.float64)
    offset_phases = {
        float(offset): np.exp(
            -2j * np.pi * signed_lags * float(offset) / float(grid.n_sc)
        )
        for offset in offsets
    }
    accepted = []
    for index in range(int(args.count)):
        increments = rng.uniform(0.2, 2.0, size=7)
        target_span = rng.uniform(float(args.span_min), float(args.span_max))
        relative = np.concatenate(([0.0], np.cumsum(increments)))
        relative *= target_span / relative[-1]
        relative_factor = np.sum(
            np.exp(
                -2j
                * np.pi
                * signed_lags[:, None]
                * relative[None, :]
                / float(grid.n_sc)
            ),
            axis=1,
        )
        best = None
        for offset in offsets:
            coordinates = relative + offset
            factor = relative_factor * offset_phases[float(offset)]
            value = evaluate_lag_nmse(
                lag_metrics[14.0],
                factor,
            )
            value_db = 10.0 * math.log10(max(value, 1e-300))
            if best is None or value_db < best[0]:
                best = (value_db, float(offset), coordinates.copy())
        if best is None or best[0] > float(args.threshold_db):
            continue
        values = {"14": best[0]}
        best_factor = relative_factor * offset_phases[best[1]]
        for snr in (16.0, 20.0):
            value = evaluate_lag_nmse(
                lag_metrics[snr],
                best_factor,
            )
            values[f"{snr:g}"] = 10.0 * math.log10(max(value, 1e-300))
        if any(value > float(args.threshold_db) for value in values.values()):
            continue
        factor = array_factor(best[2], grid.n_sc)
        rho = np.abs(physical[1:] * factor[1:])
        m4_eff = float(np.sum(lag_weights * rho**4))
        accepted.append(
            {
                "sample_index": index,
                "offset": best[1],
                "delay_grid_coordinates": best[2].tolist(),
                "span": float(best[2][-1] - best[2][0]),
                "analytic_trace_nmse_db": values,
                "m4_eff": m4_eff,
            }
        )
    accepted.sort(key=lambda row: (float(row["m4_eff"]), int(row["sample_index"])))
    payload = {
        "schema": "plan028-cdd-nonuniform-search-v1",
        "plan": "research/plan-028-comb6三类CSI-TDL-A300ns.md#13-a100-小时延透明-cdd分集与信道估计联合验收",
        "source_manifest_sha256": source_digest,
        "seed": int(args.seed),
        "sample_count": int(args.count),
        "accepted_count": len(accepted),
        "threshold_db": float(args.threshold_db),
        "span_min": float(args.span_min),
        "span_max": float(args.span_max),
        "top_candidates": accepted[:20],
    }
    output = runner._scene_output(config, scene) / str(args.output_name)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    digest = runner._sha256(output)
    output.with_suffix(".sha256").write_text(
        f"{digest}  {output.name}\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2))
    print(f"sha256={digest}")


if __name__ == "__main__":
    main()
