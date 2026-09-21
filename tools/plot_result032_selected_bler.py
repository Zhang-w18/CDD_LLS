"""Plot the selected three-curve estimated-CSI BLER views for result-032."""

from __future__ import annotations

import argparse
import csv
import shutil
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    ROOT
    / "outputs/experiment032_tdl_mobility/20260906_main/final"
    / "estimated_csi_bler_points.csv"
)
DEFAULT_OUTPUT_DIR = ROOT / "docs/figures/result-032"

PRG = "A100_PRG_DFT8_6RB"
MRT = "A100_AGED_CSI_MRT_PRG6_SLOTS10"

FIGURES = (
    (
        "A100_S0_SIDON",
        "a100_v60_estimated_csi_bler_sidon_comparison.png",
    ),
    (
        "A100_B0_QC",
        "a100_v60_estimated_csi_bler_b0qc_comparison.png",
    ),
)

DISPLAY = {
    "selected_cdd": {
        "label": "Non-transparent CDD",
        "color": "#1f77b4",
        "linestyle": "-",
        "marker": "o",
    },
    PRG: {
        "label": "Precoder cycling",
        "color": "#ff7f0e",
        "linestyle": "-.",
        "marker": "s",
    },
    MRT: {
        "label": "Closed-loop MRT",
        "color": "#2ca02c",
        "linestyle": "--",
        "marker": "^",
    },
}


def _read_points(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _plot(points: list[dict[str, str]], cdd_id: str, output: Path) -> None:
    requested = (cdd_id, PRG, MRT)
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in points:
        candidate_id = row["candidate_id"]
        if candidate_id in requested and float(row["snr_db"]) <= 20.0:
            grouped[candidate_id].append(row)

    missing = [candidate_id for candidate_id in requested if not grouped[candidate_id]]
    if missing:
        raise RuntimeError(f"Missing result-032 curves: {missing}")

    figure, axis = plt.subplots(figsize=(15.0, 9.0))
    for candidate_id in requested:
        rows = sorted(grouped[candidate_id], key=lambda row: float(row["snr_db"]))
        snr = np.asarray([float(row["snr_db"]) for row in rows])
        bler = np.asarray([float(row["bler"]) for row in rows])
        trials = np.asarray([float(row["trials"]) for row in rows])
        style = DISPLAY["selected_cdd" if candidate_id == cdd_id else candidate_id]
        axis.plot(
            snr,
            np.maximum(bler, 0.5 / trials),
            label=style["label"],
            color=style["color"],
            linestyle=style["linestyle"],
            marker=style["marker"],
            linewidth=3.5,
            markersize=10,
            markeredgewidth=1.5,
            markevery=1,
            zorder=3,
        )

    axis.set_yscale("log")
    axis.set_ylim(5e-3, 1.1)
    axis.axhline(1e-1, color="0.25", linewidth=1.8, linestyle=":", zorder=0)
    axis.axhline(1e-2, color="0.25", linewidth=1.8, linestyle=":", zorder=0)
    axis.set_xlabel("SNR (dB)", fontsize=21)
    axis.set_ylabel("Estimated-CSI BLER", fontsize=21)
    axis.tick_params(labelsize=17)
    axis.grid(True, which="both", alpha=0.28)
    axis.legend(
        loc="upper right",
        fontsize=17,
        frameon=False,
    )
    figure.subplots_adjust(bottom=0.12, left=0.10, right=0.98, top=0.98)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--copy-to-final",
        action="store_true",
        help="Also copy the figures beside the final result-032 CSV.",
    )
    args = parser.parse_args()

    points = _read_points(args.input)
    for cdd_id, filename in FIGURES:
        output = args.output_dir / filename
        _plot(points, cdd_id, output)
        if args.copy_to_final:
            destination = args.input.parent / filename
            if output.resolve() != destination.resolve():
                shutil.copy2(output, destination)
        print(f"[plot] {output}", flush=True)


if __name__ == "__main__":
    main()
