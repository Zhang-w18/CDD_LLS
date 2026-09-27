"""Freeze plan-038 formal grids from the completed 300-trial prescan."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import yaml


ROOT = Path(__file__).resolve().parents[1]
TARGETS = (0.01,)
TRANSPARENT_TO_MATCHED = {
    "SIDON_TRANSPARENT": "SIDON_MATCHED",
    "B0QC_TRANSPARENT": "B0QC_MATCHED",
    "CDD911_TRANSPARENT": "CDD911_MATCHED",
    "CDD130_TRANSPARENT": "CDD130_MATCHED",
}


def _read_points(path: Path) -> list[tuple[float, float]]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing prescan result: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Empty prescan result: {path}")
    return sorted((float(row["snr_db"]), float(row["bler"])) for row in rows)


def _grid_for(points: list[tuple[float, float]]) -> tuple[list[float], dict[str, object]]:
    selected: set[float] = set()
    receipt: dict[str, object] = {}
    x = np.asarray([item[0] for item in points], dtype=float)
    y = np.asarray([item[1] for item in points], dtype=float)
    for target in TARGETS:
        brackets = [
            (float(x[i]), float(x[i + 1]))
            for i in range(len(x) - 1)
            if y[i] >= target and y[i + 1] <= target
        ]
        key = f"bler_{target:g}"
        if brackets:
            lo, hi = brackets[0]
            start, stop = lo, hi
            status = "bracketed"
        else:
            distance = np.abs(np.log10(np.maximum(y, 0.5 / 300.0)) - np.log10(target))
            center = float(x[int(np.argmin(distance))])
            start, stop = center - 0.5, center + 0.5
            status = "nearest-window-no-prescan-bracket"
        values = np.arange(start, stop + 0.125, 0.25)
        selected.update(round(float(value), 6) for value in values)
        receipt[key] = {"status": status, "prescan_brackets": brackets, "formal_window_db": [start, stop]}
    return sorted(selected), receipt


def main() -> None:
    full_receipt: dict[str, object] = {"targets": list(TARGETS), "selection": {}}
    for al in (1, 2, 4):
        source = ROOT / f"configs/pdcch_result038_al{al}_prescan.yaml"
        destination = ROOT / f"configs/pdcch_result038_al{al}_formal.yaml"
        with source.open("r", encoding="utf-8") as handle:
            config = yaml.safe_load(handle)
        config["output_dir"] = f"outputs/experiment038_pdcch_a40/20260922_4tx_2rx_2sym/formal/al{al}"
        config["random_stream_namespace"] = f"plan038-a40-al{al}-formal-v1"
        config["simulation"].update(
            {
                "min_trials_per_snr": 10000,
                "target_errors": 200,
                "max_trials_per_snr": 50000,
                "target_bler": 0.01,
                "progress_every_batches": 50,
            }
        )
        frozen: dict[str, tuple[list[float], dict[str, object]]] = {}
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            if candidate_id in TRANSPARENT_TO_MATCHED:
                continue
            points = _read_points(
                ROOT
                / f"outputs/experiment038_pdcch_a40/20260922_4tx_2rx_2sym/prescan/al{al}"
                / candidate_id
                / "bler_points.csv"
            )
            grid, receipt = _grid_for(points)
            frozen[candidate_id] = (grid, receipt)
            candidate["snr_points_db"] = grid
            full_receipt["selection"][f"al{al}:{candidate_id}"] = receipt | {"formal_grid_db": grid}
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            if candidate_id not in TRANSPARENT_TO_MATCHED:
                continue
            source_id = TRANSPARENT_TO_MATCHED[candidate_id]
            grid = list(frozen[source_id][0])
            candidate["snr_points_db"] = grid
            full_receipt["selection"][f"al{al}:{candidate_id}"] = {
                "status": "reused-non-transparent-grid",
                "grid_source_candidate": source_id,
                "formal_grid_db": grid,
                "transparent_target_crossing_does_not_expand_grid": True,
            }
        destination.write_text(
            yaml.safe_dump(config, sort_keys=False, allow_unicode=True), encoding="utf-8"
        )
    receipt_path = ROOT / "outputs/experiment038_pdcch_a40/20260922_4tx_2rx_2sym/formal_grid_receipt.json"
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(full_receipt, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote formal configurations and {receipt_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
