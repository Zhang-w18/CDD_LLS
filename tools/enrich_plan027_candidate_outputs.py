"""Mechanically add required delay/rank fields to existing Plan-027 raw CSVs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.design import (
    canonical_delay_set,
    delay_indices_to_ns,
    residues_and_lifts,
)
from tools.run_plan027_meff_design import (
    DEFAULT_OUTPUT_ROOT,
    DEFAULT_RUN_ID,
    K,
    N_TX,
    PILOT_PERIOD,
    SCS_HZ,
)


EXTRA_FIELDS = (
    "delay_grid_coordinates",
    "delay_ns",
    "residues",
    "lifts",
    "canonical_key",
    "grid_aligned",
    "pilot_rank",
    "pilot_condition_number",
)


def encode(values) -> str:
    return json.dumps(list(values), separators=(",", ":"))


def enrich(path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with path.open("r", encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source)
        if all(field in (reader.fieldnames or []) for field in EXTRA_FIELDS):
            return
        fields = list(reader.fieldnames or [])
        insert_at = fields.index("delay_indices") + 1
        fields[insert_at:insert_at] = list(EXTRA_FIELDS)
        with temporary.open("w", encoding="utf-8", newline="") as destination:
            writer = csv.DictWriter(destination, fieldnames=fields)
            writer.writeheader()
            for row in reader:
                indices = canonical_delay_set(
                    [int(value) for value in json.loads(row["delay_indices"])],
                    K,
                )
                residues, lifts = residues_and_lifts(indices, PILOT_PERIOD)
                row.update(
                    {
                        "delay_indices": encode(indices),
                        "delay_grid_coordinates": encode(indices),
                        "delay_ns": encode(delay_indices_to_ns(indices, K, SCS_HZ)),
                        "residues": encode(residues),
                        "lifts": encode(lifts),
                        "canonical_key": encode(indices),
                        "grid_aligned": True,
                        "pilot_rank": N_TX,
                        "pilot_condition_number": 1.0,
                    }
                )
                writer.writerow(row)
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", default=DEFAULT_RUN_ID)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    args = parser.parse_args()
    output = (args.output_root / str(args.run_id)).resolve()
    workspace = DEFAULT_OUTPUT_ROOT.resolve().parent
    if workspace not in output.parents:
        raise RuntimeError("Resolved output path is outside the Plan-027 output root.")
    enrich(output / "e1_search" / "e1_all_candidates.csv")
    enrich(output / "e2_sidon_map" / "e2_sidon_candidates.csv")


if __name__ == "__main__":
    main()
