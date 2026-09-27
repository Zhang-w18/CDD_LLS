"""Plot Plan-039 stage-1A CE-only NMSE prescan curves."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


DEFAULT_SELECTION_DIR = Path(
    "outputs/experiment039_cdl_beam_bler/stage0_20260925/selection"
)

CURVES = {
    "b0__plan039_common_reference_pdp": (
        "Method 1 | B0 QC | common reference PDP",
        "tab:blue",
        "o-",
    ),
    "b0__plan039_beam_specific_pdp_independent": (
        "Method 2 | B0 QC | beam-specific PDP",
        "tab:blue",
        "X--",
    ),
    "sidon_01__plan039_common_reference_pdp": (
        "Method 1 | Sidon | common reference PDP",
        "tab:orange",
        "o-",
    ),
    "sidon_01__plan039_beam_specific_pdp_independent": (
        "Method 2 | Sidon | beam-specific PDP",
        "tab:orange",
        "X--",
    ),
    "cycling__plan039_prg_common_reference_pdp": (
        "Method 1 | PRG6 cycling | PRG-local common PDP",
        "tab:green",
        "o-",
    ),
}


def _find_summary(selection_dir: Path) -> Path:
    merged = sorted(selection_dir.glob("*/summary.csv"))
    if merged:
        return merged[0]
    batches = sorted(selection_dir.glob("batches/batch_*/summary.csv"))
    if batches:
        return batches[-1]
    raise FileNotFoundError(f"No selection summary.csv found under {selection_dir}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selection-dir", type=Path, default=DEFAULT_SELECTION_DIR)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    summary_path = _find_summary(args.selection_dir)
    output_dir = args.output_dir or args.selection_dir / "C_UE2R_DUAL_ASD25"
    output_dir.mkdir(parents=True, exist_ok=True)

    points: dict[str, list[tuple[float, float, int]]] = {key: [] for key in CURVES}
    with summary_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            variant_id = row.get("variant_id", "")
            if variant_id not in points:
                continue
            nmse = float(row["ce_nmse_eff"])
            if not math.isfinite(nmse) or nmse <= 0.0:
                raise ValueError(f"Invalid ce_nmse_eff for {variant_id}: {nmse}")
            trials = int(float(row.get("trials", row.get("n_trials", "0"))))
            points[variant_id].append((float(row["snr_db"]), 10.0 * math.log10(nmse), trials))

    missing = [variant_id for variant_id, values in points.items() if not values]
    if missing:
        raise ValueError(f"Missing expected CE-only variants: {missing}")

    plot_data_path = output_dir / "nmse_prescan_data.csv"
    with plot_data_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["variant_id", "snr_db", "ce_nmse_eff", "ce_nmse_db", "trials"])
        for variant_id, values in points.items():
            for snr_db, nmse_db, trials in sorted(values):
                writer.writerow([variant_id, snr_db, 10.0 ** (nmse_db / 10.0), nmse_db, trials])

    fig, ax = plt.subplots(figsize=(10.5, 6.5), constrained_layout=True)
    for variant_id, (label, color, style) in CURVES.items():
        values = sorted(points[variant_id])
        ax.plot(
            [value[0] for value in values],
            [value[1] for value in values],
            style,
            color=color,
            linewidth=1.8,
            markersize=7 if "Method 2" in label else 5,
            markeredgewidth=1.3,
            label=label,
        )
    ax.set_title("Plan-039 stage 1A CE-only NMSE prescan")
    ax.set_xlabel("SNR (dB)")
    ax.set_ylabel("Effective-channel NMSE (dB)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    figure_path = output_dir / "nmse_prescan.png"
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)
    print(figure_path.resolve())
    print(plot_data_path.resolve())


if __name__ == "__main__":
    main()
