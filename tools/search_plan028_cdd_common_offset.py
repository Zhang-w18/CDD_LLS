"""Search common CDD offsets for plan-028 section 13 using analytic CE NMSE."""

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

from cdd_lls.phy.estimators import (
    build_frequency_rmmse_filter,
    tdl_active_frequency_covariance,
)
from cdd_lls.phy.precoding import build_active_dft_grid_precoder
from tools import run_bler_curves as runner
from tools import run_plan027_bler as bler027


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--steps", nargs="+", type=float, required=True)
    parser.add_argument("--offset-min", type=float, default=-4.0)
    parser.add_argument("--offset-max", type=float, default=2.0)
    parser.add_argument("--offset-step", type=float, default=0.25)
    parser.add_argument("--snr-db", nargs="+", type=float, default=[14.0, 16.0, 20.0])
    parser.add_argument("--threshold-db", type=float, default=-16.5)
    args = parser.parse_args()

    config = runner.load_runner_config(args.config.resolve())
    runner.validate_config(config)
    scene = config["scenes"][0]
    manifest, source_digest = runner._load_manifest(scene)
    seed_candidate = runner._selected_candidates(scene, manifest)[0]
    channel, grid, _, _ = bler027.build_scene(
        manifest, [seed_candidate], str(scene["scenario_id"])
    )
    physical_covariance = tdl_active_frequency_covariance(grid, channel)
    assumed_covariance = 8.0 * physical_covariance
    pilot_local = runner.local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    data_local = runner.local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    filters = {}
    noise_variances = {}
    for snr_db in args.snr_db:
        ls_noise_variance = 4.0 / (10.0 ** (float(snr_db) / 10.0))
        noise_variances[float(snr_db)] = ls_noise_variance
        filters[float(snr_db)] = build_frequency_rmmse_filter(
            assumed_covariance,
            pilot_local,
            ls_noise_variance,
            diagonal_loading=1e-10,
        )

    offsets = np.arange(
        float(args.offset_min),
        float(args.offset_max) + 0.5 * float(args.offset_step),
        float(args.offset_step),
    )
    best_rows = []
    for step in args.steps:
        trials = []
        for offset in offsets:
            coordinates = offset + float(step) * np.arange(8, dtype=np.float64)
            precoder = build_active_dft_grid_precoder(
                grid, coordinates.tolist(), n_tx=8, normalize=False
            )
            true_covariance = physical_covariance * (
                precoder.C @ precoder.C.conj().T
            )
            values = {}
            for snr_db in args.snr_db:
                value = runner.mismatched_frequency_lmmse_trace_nmse(
                    true_covariance,
                    filters[float(snr_db)].weights,
                    pilot_local,
                    data_local,
                    noise_variances[float(snr_db)],
                )
                values[f"{snr_db:g}"] = 10.0 * math.log10(max(value, 1e-300))
            trials.append(
                {
                    "offset": float(offset),
                    "delay_grid_coordinates": coordinates.tolist(),
                    "analytic_trace_nmse_db": values,
                    "worst_nmse_db": max(values.values()),
                }
            )
        best = min(trials, key=lambda row: float(row["worst_nmse_db"]))
        best_rows.append(
            {
                "step": float(step),
                **best,
                "passes_all_points": all(
                    value <= float(args.threshold_db)
                    for value in best["analytic_trace_nmse_db"].values()
                ),
            }
        )

    payload = {
        "schema": "plan028-cdd-common-offset-search-v1",
        "plan": "research/plan-028-comb6三类CSI-TDL-A300ns.md#13-a100-小时延透明-cdd分集与信道估计联合验收",
        "source_manifest_sha256": source_digest,
        "steps": [float(value) for value in args.steps],
        "offset_grid": {
            "minimum": float(args.offset_min),
            "maximum": float(args.offset_max),
            "step": float(args.offset_step),
        },
        "snr_db": [float(value) for value in args.snr_db],
        "threshold_db": float(args.threshold_db),
        "best_by_step": best_rows,
    }
    output = runner._scene_output(config, scene) / "common_offset_search.json"
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
