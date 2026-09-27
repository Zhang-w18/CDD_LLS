"""Generate and run plan-038 uniform 0.25 dB refinement grids."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import copy
from pathlib import Path
import subprocess
import sys

import yaml


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "tools" / "run_pdcch_bler_curves.py"
CONFIG_DIR = ROOT / "outputs/experiment038_pdcch_a40/20260922_4tx_2rx_2sym/refinement_configs"

# Inclusive plotting/refinement windows. The two user-unspecified cases retain their
# existing formal ranges so the final figure still contains all nine schemes per AL.
WINDOWS: dict[tuple[int, str], tuple[float, float]] = {
    (1, "SIDON_MATCHED"): (4.0, 6.0),
    (1, "SIDON_TRANSPARENT"): (4.0, 5.0),
    (1, "B0QC_MATCHED"): (5.0, 8.0),
    (1, "B0QC_TRANSPARENT"): (5.0, 8.0),
    (1, "CDD911_MATCHED"): (6.0, 8.0),
    (1, "CDD911_TRANSPARENT"): (6.0, 8.0),
    (1, "CDD130_MATCHED"): (7.0, 9.0),
    (1, "CDD130_TRANSPARENT"): (7.0, 9.0),
    (1, "PRG_DFT4_TRANSPARENT"): (7.0, 10.0),
    (2, "SIDON_MATCHED"): (0.0, 1.0),
    (2, "SIDON_TRANSPARENT"): (4.0, 5.0),
    (2, "B0QC_MATCHED"): (0.0, 1.5),
    (2, "B0QC_TRANSPARENT"): (1.5, 2.5),
    (2, "CDD911_MATCHED"): (0.5, 2.5),
    (2, "CDD911_TRANSPARENT"): (2.0, 3.0),
    (2, "CDD130_MATCHED"): (2.5, 4.0),
    (2, "CDD130_TRANSPARENT"): (2.5, 4.0),
    (2, "PRG_DFT4_TRANSPARENT"): (1.5, 3.0),
    (4, "SIDON_MATCHED"): (-4.0, -2.5),
    (4, "SIDON_TRANSPARENT"): (-2.0, -1.0),
    (4, "B0QC_MATCHED"): (-4.0, -2.5),
    (4, "B0QC_TRANSPARENT"): (-1.0, 0.0),
    (4, "CDD911_MATCHED"): (-3.5, -2.0),
    (4, "CDD911_TRANSPARENT"): (-0.5, 0.5),
    (4, "CDD130_MATCHED"): (-2.5, -1.0),
    (4, "CDD130_TRANSPARENT"): (-2.5, -1.0),
    (4, "PRG_DFT4_TRANSPARENT"): (-3.5, -2.0),
}


def _quarter_db_grid(start: float, stop: float) -> list[float]:
    first, last = round(start * 4), round(stop * 4)
    if abs(first / 4 - start) > 1e-12 or abs(last / 4 - stop) > 1e-12:
        raise ValueError(f"Window is not aligned to 0.25 dB: {start}, {stop}")
    return [value / 4 for value in range(first, last + 1)]


def _write_configs() -> list[Path]:
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    generated: list[Path] = []
    for al in (1, 2, 4):
        source = ROOT / f"configs/pdcch_result038_al{al}_formal.yaml"
        with source.open("r", encoding="utf-8") as handle:
            config = yaml.safe_load(handle)
        refined = copy.deepcopy(config)
        for candidate in refined["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            candidate["snr_points_db"] = _quarter_db_grid(*WINDOWS[(al, candidate_id)])
        destination = CONFIG_DIR / f"pdcch_result038_al{al}_refinement.yaml"
        destination.write_text(
            yaml.safe_dump(refined, sort_keys=False, allow_unicode=True), encoding="utf-8"
        )
        generated.append(destination)
    return generated


def _run_one(config: Path, candidate: str, stage: str) -> tuple[str, int, str]:
    result = subprocess.run(
        [sys.executable, str(RUNNER), "--config", str(config), "--candidate", candidate, "--stage", stage],
        cwd=ROOT,
        text=True,
        encoding="utf-8",
        errors="replace",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    return f"{config.name}:{candidate}", int(result.returncode), result.stdout


def _needs_run(config: dict[str, object], candidate: dict[str, object]) -> bool:
    output = Path(str(config["output_dir"])) / str(candidate["candidate_id"]) / "bler_points.csv"
    if not output.is_file():
        return True
    import csv

    with output.open("r", encoding="utf-8", newline="") as handle:
        completed = {float(row["snr_db"]) for row in csv.DictReader(handle)}
    requested = {float(value) for value in candidate["snr_points_db"]}
    return not requested.issubset(completed)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("prepare", "validate", "run"), required=True)
    parser.add_argument("--max-workers", type=int, default=3)
    args = parser.parse_args()
    configs = _write_configs()
    if args.stage == "prepare":
        for path in configs:
            print(path.relative_to(ROOT))
        return
    tasks: list[tuple[Path, str]] = []
    for path in configs:
        with path.open("r", encoding="utf-8") as handle:
            config = yaml.safe_load(handle)
        tasks.extend(
            (path, str(item["candidate_id"]))
            for item in config["candidates"]
            if args.stage == "validate" or _needs_run(config, item)
        )
    if not tasks:
        print("All requested refinement points are already complete.")
        return
    failures: list[str] = []
    with ThreadPoolExecutor(max_workers=int(args.max_workers)) as pool:
        futures = {pool.submit(_run_one, path, candidate, args.stage): (path, candidate) for path, candidate in tasks}
        for future in as_completed(futures):
            name, returncode, output = future.result()
            print(f"===== {name} :: exit {returncode} =====", flush=True)
            print(output.rstrip(), flush=True)
            if returncode:
                failures.append(name)
    if failures:
        raise SystemExit("Failed refinement shards: " + ", ".join(failures))


if __name__ == "__main__":
    main()
