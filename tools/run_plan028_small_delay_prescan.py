"""Run plan-028 section 13 prescan candidates in isolated Python processes.

The decoded-link stack retains a large TensorFlow graph when many transparent
CDD candidates are decoded in one process.  This controller preserves the
master YAML, seed derivation, and per-candidate definitions while giving each
candidate its own output directory and fresh decoder process.
"""

from __future__ import annotations

import argparse
import copy
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "tools" / "run_bler_curves.py"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--stage", choices=("validate", "run", "merge"), default="run")
    args = parser.parse_args()

    master_path = args.config.resolve()
    master = yaml.safe_load(master_path.read_text(encoding="utf-8"))
    scenes = master.get("scenes", [])
    if len(scenes) != 1:
        raise ValueError("Prescan controller requires exactly one scene.")
    scene = scenes[0]
    baselines = scene.get("transparent_cdd_baselines", [])
    if not baselines:
        raise ValueError("No transparent CDD candidates are configured.")
    receiver_modes = sorted(scene.get("receivers", {}))
    if len(receiver_modes) != 1:
        raise ValueError("Prescan controller requires exactly one receiver mode.")
    receiver = receiver_modes[0]
    output_root = Path(master["output_dir"])
    if not output_root.is_absolute():
        output_root = ROOT / output_root
    config_root = output_root / "control" / "isolated_configs" / receiver
    config_root.mkdir(parents=True, exist_ok=True)

    for index, baseline in enumerate(baselines, start=1):
        candidate_id = str(baseline["candidate_id"])
        short_id = f"c{index:02d}"
        isolated = copy.deepcopy(master)
        isolated["output_dir"] = str(output_root / short_id)
        isolated["scenes"][0]["transparent_cdd_baselines"] = [copy.deepcopy(baseline)]
        isolated_path = config_root / f"{short_id}.yaml"
        isolated_path.write_text(
            yaml.safe_dump(isolated, sort_keys=False, allow_unicode=True),
            encoding="utf-8",
        )
        print(
            f"[prescan controller] receiver={receiver} slot={short_id} "
            f"candidate={candidate_id}",
            flush=True,
        )
        subprocess.run(
            [
                sys.executable,
                str(RUNNER),
                "--config",
                str(isolated_path),
                "--stage",
                args.stage,
            ],
            cwd=ROOT,
            check=True,
        )


if __name__ == "__main__":
    main()
