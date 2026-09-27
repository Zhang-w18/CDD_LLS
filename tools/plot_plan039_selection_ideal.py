"""Plot the three Plan-039 stage-1A ideal-CSI BLER prescan curves."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCENARIO = "C_UE2R_DUAL_ASD25"
CURVES = {
    "b0__ideal": ("B0 QC — ideal CSI", "#1f77b4", "o"),
    "sidon_01__ideal": ("Sidon [0,11,28,148,170,233,277,351] — ideal CSI", "#d62728", "s"),
    "cycling__ideal": ("PRG6 cycling — ideal CSI", "#2ca02c", "^")
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path,
                        default=ROOT / "outputs" / "experiment039_cdl_beam_bler" / "stage0_20260925")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    run_root = args.run_root.resolve()
    merged = run_root / "selection" / SCENARIO / "summary.csv"
    batch = run_root / "selection" / "batches" / "batch_01" / "summary.csv"
    source = merged if merged.exists() else batch
    if not source.exists():
        raise FileNotFoundError(
            "Stage-1A summary is not on disk yet. The running batch writes summary.csv only after all variants finish."
        )
    rows = read_rows(source)
    available = {row["variant_id"] for row in rows}
    missing = sorted(set(CURVES) - available)
    if missing:
        raise RuntimeError(f"Stage-1A ideal curves are incomplete: missing {missing}.")
    output = (args.output or run_root / "selection" / SCENARIO / "ideal_bler_prescan.png").resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    plot_rows: list[dict[str, object]] = []
    fig, ax = plt.subplots(figsize=(8.2, 5.2), constrained_layout=True)
    for variant, (label, color, marker) in CURVES.items():
        curve = sorted((row for row in rows if row["variant_id"] == variant),
                       key=lambda row: float(row["snr_db"]))
        snr = np.asarray([float(row["snr_db"]) for row in curve])
        trials = np.asarray([int(row["n_trials"]) for row in curve])
        errors = np.asarray([int(row["tb_errors"]) for row in curve])
        bler = errors / trials
        display = np.where(errors == 0, 0.5 / trials, bler)
        ax.semilogy(snr, display, color=color, marker=marker, linewidth=1.8,
                    markersize=5.5, label=label)
        zero = errors == 0
        if np.any(zero):
            ax.scatter(snr[zero], display[zero], color=color, marker="v", s=38, facecolors="none")
        for index in range(len(curve)):
            plot_rows.append({"variant_id": variant, "snr_db": float(snr[index]),
                              "trials": int(trials[index]), "errors": int(errors[index]),
                              "bler_raw": float(bler[index]), "bler_display": float(display[index]),
                              "zero_error_display_rule": "0.5/trials" if errors[index] == 0 else "raw"})
    ax.axhline(0.10, color="0.35", linestyle="--", linewidth=1.0, label="10% BLER")
    ax.axhline(0.01, color="0.55", linestyle=":", linewidth=1.0, label="1% BLER")
    ax.set_xlabel("Reference SNR (dB)")
    ax.set_ylabel("BLER")
    ax.set_title("Plan-039 stage 1A ideal-CSI BLER prescan")
    ax.grid(True, which="both", linestyle=":", linewidth=0.6, alpha=0.7)
    ax.legend(fontsize=8)
    fig.savefig(output, dpi=180)
    plt.close(fig)
    csv_path = output.with_name(output.stem + "_data.csv")
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(plot_rows[0]))
        writer.writeheader()
        writer.writerows(plot_rows)
    print(output)
    print(csv_path)


if __name__ == "__main__":
    main()
