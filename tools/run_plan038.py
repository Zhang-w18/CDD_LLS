"""Orchestrate plan-038 validation, prescan, formal runs, and analysis."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import subprocess
import sys

import yaml


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "tools" / "run_pdcch_bler_curves.py"
PREPARE = ROOT / "tools" / "prepare_plan038_formal.py"
ANALYZE = ROOT / "tools" / "analyze_plan038.py"
PRESCAN_CONFIGS = [ROOT / f"configs/pdcch_result038_al{al}_prescan.yaml" for al in (1, 2, 4)]
FORMAL_CONFIGS = [ROOT / f"configs/pdcch_result038_al{al}_formal.yaml" for al in (1, 2, 4)]
PRESCAN_SOURCES = {
    "SIDON_MATCHED",
    "B0QC_MATCHED",
    "CDD911_MATCHED",
    "CDD130_MATCHED",
    "PRG_DFT4_TRANSPARENT",
}


def _candidate_ids(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    return [str(item["candidate_id"]) for item in config["candidates"]]


def _run_shard(path: Path, candidate: str, stage: str) -> tuple[str, int, str]:
    command = [
        sys.executable,
        str(RUNNER),
        "--config",
        str(path),
        "--candidate",
        candidate,
        "--stage",
        stage,
    ]
    result = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        encoding="utf-8",
        errors="replace",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    return f"{path.name}:{candidate}", int(result.returncode), result.stdout


def _run_configs(paths: list[Path], stage: str, workers: int, candidate_filter: set[str] | None = None) -> None:
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing configuration(s): " + ", ".join(missing))
    tasks = [
        (path, candidate)
        for path in paths
        for candidate in _candidate_ids(path)
        if candidate_filter is None or candidate in candidate_filter
    ]
    failures: list[str] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_run_shard, path, candidate, stage): (path, candidate) for path, candidate in tasks}
        for future in as_completed(futures):
            name, returncode, output = future.result()
            print(f"===== {name} :: exit {returncode} =====", flush=True)
            print(output.rstrip(), flush=True)
            if returncode:
                failures.append(name)
    if failures:
        raise SystemExit("Failed shards: " + ", ".join(failures))


def _call(script: Path) -> None:
    subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        required=True,
        choices=("validate", "prescan", "prepare", "validate-formal", "formal", "analyze", "all"),
    )
    parser.add_argument("--max-workers", type=int, default=3)
    args = parser.parse_args()
    workers = int(args.max_workers)
    if workers < 1:
        raise ValueError("--max-workers must be positive")

    if args.stage in {"validate", "all"}:
        _run_configs(PRESCAN_CONFIGS, "validate", workers)
    if args.stage in {"prescan", "all"}:
        _run_configs(PRESCAN_CONFIGS, "run", workers, PRESCAN_SOURCES)
    if args.stage in {"prepare", "all"}:
        _call(PREPARE)
    if args.stage in {"validate-formal", "all"}:
        _run_configs(FORMAL_CONFIGS, "validate", workers)
    if args.stage in {"formal", "all"}:
        _run_configs(FORMAL_CONFIGS, "run", workers)
    if args.stage in {"analyze", "all"}:
        _call(ANALYZE)


if __name__ == "__main__":
    main()
