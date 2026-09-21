"""Run the plan-035 supplement with a CDL-D channel and unchanged link settings."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools import run_plan033_tdl_mobility_mimo as core


SCHEMA = "plan035-cdl-d-pdsch-2rx-v1"


def load_config(path: Path) -> dict:
    config = core.load_config(path, expected_schema=SCHEMA)
    config["plan_path"] = "research/plan-035-PDSCH 2Rx仿真.md"
    return config


def validate_config(config: dict):
    return core.validate_config(config, plan="035-cdl")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--stage", choices=("validate", "run", "merge", "adaptive"), default="run")
    args = parser.parse_args()
    config = load_config(args.config.resolve())
    candidates, channel, grid, precoders, base = validate_config(config)
    if args.stage == "run":
        core.run(config, candidates, channel, grid, precoders, base)
    elif args.stage == "adaptive":
        if "adaptive_policy" not in config:
            raise ValueError("--stage adaptive requires adaptive_policy in the config.")
        core.run_adaptive(config, candidates, channel, grid, precoders, base)
    elif args.stage == "merge":
        core._merge_intervals(Path(config["output_dir"]), candidates, config["scenario_id"])


if __name__ == "__main__":
    main()
