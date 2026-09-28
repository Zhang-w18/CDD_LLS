"""Plot selected result-032 BLER views including transparent CDD supplements."""

from __future__ import annotations

import argparse
import csv
import shutil
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    ROOT
    / "outputs/experiment032_tdl_mobility/20260906_main/final"
    / "estimated_csi_bler_points.csv"
)
DEFAULT_OUTPUT_DIR = ROOT / "docs/figures/result-032"
DEFAULT_TRANSPARENT_INPUT = (
    ROOT
    / "outputs/experiment032_tdl_mobility/20260906_main/transparent_cdd_supplement/final"
    / "estimated_csi_bler_points.csv"
)

PRG = "A100_PRG_DFT8_6RB"
MRT = "A100_AGED_CSI_MRT_PRG6_SLOTS10"
SIDON = "A100_S0_SIDON"
B0QC = "A100_B0_QC"
SIDON_TRANSPARENT = "A100_S0_SIDON_TRANSPARENT_CDD"
B0QC_TRANSPARENT = "A100_B0_QC_TRANSPARENT_CDD"
NMSE_FILENAME = "a100_v60_ce_nmse_transparent_cdd_comparison.png"
NMSE_PREVIEW_FILENAME = "a100_v60_ce_nmse_transparent_cdd_comparison_preview_13cm.png"

FIGURES = (
    (
        SIDON,
        SIDON_TRANSPARENT,
        "a100_v60_estimated_csi_bler_sidon_comparison.png",
    ),
    (
        B0QC,
        B0QC_TRANSPARENT,
        "a100_v60_estimated_csi_bler_b0qc_comparison.png",
    ),
)

CDD_DISPLAY = {
    SIDON: {
        "color": "#1f77b4",
        "linestyle": "-",
        "marker": "o",
    },
    SIDON_TRANSPARENT: {
        "color": "#1f77b4",
        "linestyle": ":",
        "marker": "X",
    },
    B0QC: {
        "color": "#9467bd",
        "linestyle": "-",
        "marker": "s",
    },
    B0QC_TRANSPARENT: {
        "color": "#9467bd",
        "linestyle": ":",
        "marker": "P",
    },
}

DISPLAY = {
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


def _plot(
    points: list[dict[str, str]],
    transparent_points: list[dict[str, str]],
    cdd_id: str,
    transparent_cdd_id: str,
    output: Path,
) -> None:
    requested = (cdd_id, transparent_cdd_id, PRG, MRT)
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in [*points, *transparent_points]:
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
        if candidate_id in CDD_DISPLAY:
            style = CDD_DISPLAY[candidate_id]
            label = (
                "Transparent CDD"
                if candidate_id == transparent_cdd_id
                else "Non-transparent CDD"
            )
        else:
            style = DISPLAY[candidate_id]
            label = style["label"]
        axis.plot(
            snr,
            np.maximum(bler, 0.5 / trials),
            label=label,
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


def _plot_nmse(
    points: list[dict[str, str]],
    transparent_points: list[dict[str, str]],
    output: Path,
    preview: bool = False,
) -> None:
    requested = (SIDON, SIDON_TRANSPARENT, B0QC, B0QC_TRANSPARENT)
    labels = {
        SIDON: "S0 Sidon, non-transparent",
        SIDON_TRANSPARENT: "S0 Sidon, transparent",
        B0QC: "B0QC, non-transparent",
        B0QC_TRANSPARENT: "B0QC, transparent",
    }
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in [*points, *transparent_points]:
        candidate_id = row["candidate_id"]
        if candidate_id in requested and float(row["snr_db"]) <= 20.0:
            grouped[candidate_id].append(row)
    missing = [candidate_id for candidate_id in requested if not grouped[candidate_id]]
    if missing:
        raise RuntimeError(f"Missing result-032 NMSE curves: {missing}")

    figure, axis = plt.subplots(figsize=(5.12, 4.1) if preview else (15.0, 9.0))
    for candidate_id in requested:
        rows = sorted(grouped[candidate_id], key=lambda row: float(row["snr_db"]))
        snr = np.asarray([float(row["snr_db"]) for row in rows])
        nmse_db = np.asarray([float(row["ce_nmse_mean_db"]) for row in rows])
        style = CDD_DISPLAY[candidate_id]
        axis.plot(
            snr,
            nmse_db,
            label=labels[candidate_id],
            color=style["color"],
            linestyle=style["linestyle"],
            marker=style["marker"],
            linewidth=2.5 if preview else 3.5,
            markersize=7 if preview else 10,
            markeredgewidth=1.2 if preview else 1.5,
            markevery=2 if preview else 1,
        )
    axis.set_xlabel("SNR (dB)", fontsize=12 if preview else 21)
    axis.set_ylabel("CE NMSE (dB)", fontsize=12 if preview else 21)
    axis.tick_params(labelsize=10 if preview else 17)
    axis.grid(True, alpha=0.28)
    axis.legend(
        loc="best",
        fontsize=8.5 if preview else 15,
        frameon=False,
        ncol=1,
    )
    figure.subplots_adjust(
        bottom=0.16 if preview else 0.12,
        left=0.16 if preview else 0.10,
        right=0.98,
        top=0.98,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument(
        "--transparent-input", type=Path, default=DEFAULT_TRANSPARENT_INPUT
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--copy-to-final",
        action="store_true",
        help="Also copy the figures beside the final result-032 CSV.",
    )
    args = parser.parse_args()

    points = _read_points(args.input)
    transparent_points = _read_points(args.transparent_input)
    for cdd_id, transparent_cdd_id, filename in FIGURES:
        output = args.output_dir / filename
        _plot(points, transparent_points, cdd_id, transparent_cdd_id, output)
        if args.copy_to_final:
            destination = args.input.parent / filename
            if output.resolve() != destination.resolve():
                shutil.copy2(output, destination)
        print(f"[plot] {output}", flush=True)
    nmse_output = args.output_dir / NMSE_FILENAME
    _plot_nmse(points, transparent_points, nmse_output)
    preview_output = args.output_dir / NMSE_PREVIEW_FILENAME
    _plot_nmse(points, transparent_points, preview_output, preview=True)
    if args.copy_to_final:
        destination = args.transparent_input.parent / NMSE_FILENAME
        if nmse_output.resolve() != destination.resolve():
            shutil.copy2(nmse_output, destination)
    print(f"[plot] {nmse_output}", flush=True)
    print(f"[plot] {preview_output}", flush=True)


if __name__ == "__main__":
    main()
