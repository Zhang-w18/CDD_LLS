"""Run plan-031 AL/candidate shards with independent deterministic streams."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import subprocess
import sys

import yaml


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "tools" / "run_pdcch_bler_curves.py"


def _candidate_ids(config_path: Path) -> list[str]:
    with open(config_path, "r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle) or {}
    return [str(candidate["candidate_id"]) for candidate in payload["candidates"]]


def _run_one(config_path: Path, candidate_id: str, stage: str) -> tuple[str, str, int, str]:
    command = [
        sys.executable,
        str(RUNNER),
        "--config",
        str(config_path),
        "--candidate",
        candidate_id,
        "--stage",
        stage,
    ]
    completed = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        encoding="utf-8",
        errors="replace",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    return config_path.name, candidate_id, int(completed.returncode), completed.stdout


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", action="append", required=True)
    parser.add_argument("--stage", choices=("validate", "run"), default="run")
    parser.add_argument("--max-workers", type=int, default=3)
    args = parser.parse_args()
    configs = [(ROOT / value).resolve() for value in args.config]
    tasks = [
        (config_path, candidate_id)
        for config_path in configs
        for candidate_id in _candidate_ids(config_path)
    ]
    failures: list[str] = []
    with ThreadPoolExecutor(max_workers=int(args.max_workers)) as pool:
        futures = {
            pool.submit(_run_one, config_path, candidate_id, args.stage): (config_path, candidate_id)
            for config_path, candidate_id in tasks
        }
        for future in as_completed(futures):
            config_name, candidate_id, returncode, output = future.result()
            print(f"===== {config_name} :: {candidate_id} :: exit {returncode} =====", flush=True)
            print(output.rstrip(), flush=True)
            if returncode:
                failures.append(f"{config_name}:{candidate_id}")
    if failures:
        raise SystemExit("Failed shards: " + ", ".join(failures))


if __name__ == "__main__":
    main()
