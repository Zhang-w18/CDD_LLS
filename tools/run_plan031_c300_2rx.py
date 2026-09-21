"""Validate or run all plan-031 C300 4Tx/2Rx candidate streams."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
ALS = (1, 2, 4)
MODES = ("estimated", "ideal")


def _config_path(al: int, mode: str, stage: str) -> Path:
    return (
        ROOT
        / "configs"
        / f"pdcch_result031_c300_2rx_{mode}_al{al}_{stage}.yaml"
    )


def _run(command: list[str]) -> None:
    print(" ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("smoke", "prescan", "formal"), required=True)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    runner = ROOT / "tools" / "run_pdcch_bler_curves.py"
    for al in ALS:
        for mode in MODES:
            config_path = _config_path(al, mode, args.stage)
            if not config_path.exists():
                raise FileNotFoundError(config_path)
            if args.validate_only:
                receipt = (
                    ROOT
                    / "outputs"
                    / "experiment031_pdcch_cdd"
                    / "20260911_c300_4tx_2rx_2sym"
                    / "validation"
                    / f"{args.stage}_{mode}_al{al}.json"
                )
                _run(
                    [
                        sys.executable,
                        str(runner),
                        "--config",
                        str(config_path),
                        "--stage",
                        "validate",
                        "--validation-output",
                        str(receipt),
                    ]
                )
                continue
            with open(config_path, "r", encoding="utf-8") as handle:
                config = yaml.safe_load(handle) or {}
            for candidate in config["candidates"]:
                _run(
                    [
                        sys.executable,
                        str(runner),
                        "--config",
                        str(config_path),
                        "--stage",
                        "run",
                        "--candidate",
                        str(candidate["candidate_id"]),
                    ]
                )


if __name__ == "__main__":
    main()
