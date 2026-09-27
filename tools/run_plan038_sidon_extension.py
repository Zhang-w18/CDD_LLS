"""Extend the three requested plan-038 Sidon refinement windows."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import copy
from pathlib import Path

import yaml

from run_plan038_refinement import ROOT, _needs_run, _quarter_db_grid, _run_one


CONFIG_DIR = (
    ROOT
    / "outputs/experiment038_pdcch_a40/20260922_4tx_2rx_2sym"
    / "sidon_extension_configs"
)
WINDOWS: dict[tuple[int, str], tuple[float, float]] = {
    (1, "SIDON_MATCHED"): (4.0, 6.5),
    (2, "SIDON_MATCHED"): (0.0, 2.0),
    (2, "SIDON_TRANSPARENT"): (4.0, 6.0),
}


def _write_configs() -> list[Path]:
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    generated: list[Path] = []
    for al in (1, 2):
        source = ROOT / f"configs/pdcch_result038_al{al}_formal.yaml"
        with source.open("r", encoding="utf-8") as handle:
            config = yaml.safe_load(handle)
        extended = copy.deepcopy(config)
        selected = []
        for candidate in extended["candidates"]:
            key = (al, str(candidate["candidate_id"]))
            if key in WINDOWS:
                candidate["snr_points_db"] = _quarter_db_grid(*WINDOWS[key])
                selected.append(candidate)
        expected = {candidate_id for candidate_al, candidate_id in WINDOWS if candidate_al == al}
        actual = {str(candidate["candidate_id"]) for candidate in selected}
        if actual != expected:
            raise ValueError(f"Unexpected Sidon candidate set for AL{al}: {actual}")
        extended["candidates"] = selected
        destination = CONFIG_DIR / f"pdcch_result038_al{al}_sidon_extension.yaml"
        destination.write_text(
            yaml.safe_dump(extended, sort_keys=False, allow_unicode=True),
            encoding="utf-8",
        )
        generated.append(destination)
    return generated


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
            (path, str(candidate["candidate_id"]))
            for candidate in config["candidates"]
            if args.stage == "validate" or _needs_run(config, candidate)
        )
    if not tasks:
        print("All requested Sidon extension points are already complete.")
        return

    failures: list[str] = []
    with ThreadPoolExecutor(max_workers=int(args.max_workers)) as pool:
        futures = {
            pool.submit(_run_one, path, candidate, args.stage): (path, candidate)
            for path, candidate in tasks
        }
        for future in as_completed(futures):
            name, returncode, output = future.result()
            print(f"===== {name} :: exit {returncode} =====", flush=True)
            print(output.rstrip(), flush=True)
            if returncode:
                failures.append(name)
    if failures:
        raise SystemExit("Failed Sidon extension shards: " + ", ".join(failures))


if __name__ == "__main__":
    main()
