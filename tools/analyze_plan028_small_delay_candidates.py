"""Analytic CE prescreen for plan-028 section 13 transparent small-delay CDD."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.phy.estimators import (
    build_frequency_rmmse_filter,
    tdl_active_frequency_covariance,
)
from tools import run_bler_curves as runner
from tools import run_plan027_bler as bler027


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--snr-db", nargs="+", type=float, default=[14.0, 16.0, 20.0])
    parser.add_argument("--threshold-db", type=float, default=-16.5)
    args = parser.parse_args()

    config = runner.load_runner_config(args.config.resolve())
    runner.validate_config(config)
    if len(config["scenes"]) != 1:
        raise ValueError("Analytic prescreen requires exactly one scene.")
    scene = config["scenes"][0]
    manifest, source_digest = runner._load_manifest(scene)
    candidates = runner._selected_candidates(scene, manifest)
    channel, grid, _, true_covariances = bler027.build_scene(
        manifest, candidates, str(scene["scenario_id"])
    )
    physical_covariance = tdl_active_frequency_covariance(grid, channel)
    assumed_covariance = 8.0 * physical_covariance
    pilot_local = runner.local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    data_local = runner.local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)

    rows = []
    for candidate in candidates:
        candidate_id = str(candidate["candidate_id"])
        nmse_db = {}
        for snr_db in args.snr_db:
            noise_variance = 8.0 / (10.0 ** (float(snr_db) / 10.0))
            averaged_ls_noise_variance = noise_variance / 2.0
            estimator = build_frequency_rmmse_filter(
                assumed_covariance,
                pilot_local,
                averaged_ls_noise_variance,
                diagonal_loading=1e-10,
            )
            value = runner.mismatched_frequency_lmmse_trace_nmse(
                true_covariances[candidate_id],
                estimator.weights,
                pilot_local,
                data_local,
                averaged_ls_noise_variance,
            )
            nmse_db[f"{snr_db:g}"] = 10.0 * math.log10(max(value, 1e-300))
        rows.append(
            {
                "candidate_id": candidate_id,
                "delay_grid_coordinates": [
                    float(value) for value in candidate["delay_grid_coordinates"]
                ],
                "analytic_trace_nmse_db": nmse_db,
                "passes_all_points": all(
                    value <= float(args.threshold_db) for value in nmse_db.values()
                ),
            }
        )

    payload = {
        "schema": "plan028-small-delay-cdd-analytic-prescreen-v1",
        "plan": "research/plan-028-comb6三类CSI-TDL-A300ns.md#13-a100-小时延透明-cdd分集与信道估计联合验收",
        "source_manifest_sha256": source_digest,
        "snr_db": [float(value) for value in args.snr_db],
        "threshold_db": float(args.threshold_db),
        "assumed_covariance": "8 * R_phy(static TDL-A 100 ns)",
        "averaged_ls_noise_variance": "4/SNR_linear",
        "candidates": rows,
    }
    output = runner._scene_output(config, scene) / "analytic_prescreen.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    digest = runner._sha256(output)
    output.with_suffix(".sha256").write_text(
        f"{digest}  {output.name}\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    print(f"sha256={digest}")


if __name__ == "__main__":
    main()
