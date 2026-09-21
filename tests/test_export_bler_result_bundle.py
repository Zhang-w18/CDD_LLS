import csv
import json
from pathlib import Path

import yaml

from tools.export_bler_result_bundle import collect_config_results, render_text


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def test_collect_config_results_and_interval_audit(tmp_path: Path) -> None:
    config_path = tmp_path / "quick.yaml"
    output = tmp_path / "output"
    config = {
        "schema": "bler-curve-runner-v1",
        "output_dir": str(output),
        "seed": 7,
        "batch_size": 25,
        "max_passes": 1,
        "scenes": [
            {
                "scenario_id": "A100",
                "n_rx": 4,
                "mode": "transparent_cdd_physical_covariance",
                "receivers": {"estimated": {"snr_db": [4.0]}},
            }
        ],
    }
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    scene_output = output / "a100"
    scene_output.mkdir(parents=True)
    (scene_output / "resolved_run.json").write_text(
        json.dumps({"n_rx": 4, "seed": 7}), encoding="utf-8"
    )
    (scene_output / "run_progress.jsonl").write_text(
        json.dumps({"snr_db": 4.0, "elapsed_seconds": 12.5}) + "\n",
        encoding="utf-8",
    )
    common = {
        "candidate_id": "small",
        "snr_db": 4.0,
        "trials": 1000,
        "tb_errors": 100,
        "bler": 0.1,
        "bler_wilson95_lo": 0.0829,
        "bler_wilson95_hi": 0.1201,
        "base_trials": 0,
        "ce_nmse_mean_db": -12.0,
    }
    _write_csv(scene_output / "final" / "estimated_csi_bler_points.csv", [common])
    supplemental = {
        **common,
        "trial_start": 1,
        "trial_end": 1000,
    }
    _write_csv(scene_output / "estimated" / "supplemental_points.csv", [supplemental])

    result = collect_config_results(config_path, root=tmp_path)

    scene = result["scenes"][0]
    assert scene["run_elapsed_seconds"] == 12.5
    assert scene["receivers"]["estimated"]["interval_audit"]["ok"] is True
    assert scene["receivers"]["estimated"]["points"][0]["tb_errors"] == 100
    text = render_text(
        {
            "schema": "bler-copy-bundle-v1",
            "environment": {"platform": "test", "python": "3.x", "tensorflow": {}},
            "repository": {"git_head": "abc"},
            "runs": [result],
        }
    )
    assert "small | 4 | 1000 | 100 | 0.1" in text
