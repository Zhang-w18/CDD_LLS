"""Plot plan-031 4Tx/2Rx estimated-CSI BLER for AL2 and AL4."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import yaml


ROOT = Path(__file__).resolve().parents[1]
BASE_ROOT = (
    ROOT
    / "outputs/experiment031_pdcch_cdd/20260911_c300_4tx_2rx_2sym/formal/estimated"
)
NEW_ROOT = (
    ROOT
    / "outputs/experiment031_pdcch_cdd/20260920_c300_4tx_2rx_al2_al4_cdd911_cdd130/estimated"
)
ANALYSIS = NEW_ROOT.parent / "analysis_al2_al4_five_curves"

CURVES = {
    "C300_B0_QC": {
        "label": "B0 QC (matched)",
        "color": "#111111",
        "linestyle": "-",
        "marker": "o",
        "root": BASE_ROOT,
    },
    "C300_S0_SIDON": {
        "label": "Sidon (non-transparent)",
        "color": "#1f77b4",
        "linestyle": "--",
        "marker": "s",
        "root": BASE_ROOT,
    },
    "C300_PRG_DFT4_6REG": {
        "label": "PRG DFT4 (transparent)",
        "color": "#2ca02c",
        "linestyle": "--",
        "marker": ">",
        "root": BASE_ROOT,
    },
    "CDD911": {
        "label": "CDD911 (non-transparent)",
        "color": "#d62728",
        "linestyle": "--",
        "marker": "D",
        "root": NEW_ROOT,
    },
    "CDD130_transparent": {
        "label": "CDD130 (transparent)",
        "color": "#ff7f0e",
        "linestyle": ":",
        "marker": "X",
        "root": NEW_ROOT,
    },
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_yaml(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _load_curve(al: int, candidate_id: str, spec: dict[str, object]) -> tuple[list[dict[str, str]], dict[str, object]]:
    curve_dir = Path(spec["root"]) / f"al{al}" / candidate_id
    csv_path = curve_dir / "bler_points.csv"
    config_path = curve_dir / "resolved_config.yaml"
    if not csv_path.is_file() or not config_path.is_file():
        raise FileNotFoundError(f"Missing formal source for AL{al} {candidate_id}: {curve_dir}")
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Empty formal data: {csv_path}")
    snr_values = [float(row["snr_db"]) for row in rows]
    if len(snr_values) != len(set(snr_values)):
        raise ValueError(f"Duplicate SNR rows: {csv_path}")
    for row in rows:
        expected = {
            "aggregation_level": str(al),
            "candidate_id": candidate_id,
            "csi_mode": "frequency_lmmse",
            "n_tx": "4",
            "n_rx": "2",
            "snr_definition": "unit_total_tx_power_over_single_rx_branch_noise_power",
            "total_transmit_power_normalization": "unit_total_power_per_re",
        }
        for field, value in expected.items():
            if row[field] != value:
                raise ValueError(f"{csv_path}: {field}={row[field]!r}, expected {value!r}")
        trials = int(row["trials"])
        errors = int(row["errors"])
        if trials < 10_000 or (trials < 50_000 and errors < 200):
            raise ValueError(f"Invalid stopping state in {csv_path}: {row}")
    config = _load_yaml(config_path)
    if config["antenna"] != {"n_tx": 4, "n_rx": 2}:
        raise ValueError(f"Unexpected antenna configuration: {config_path}")
    if config["resource"]["aggregation_level"] != al:
        raise ValueError(f"Unexpected AL: {config_path}")
    if config["receiver"]["channel_estimation"] != "frequency_lmmse":
        raise ValueError(f"Unexpected receiver mode: {config_path}")
    return sorted(rows, key=lambda row: float(row["snr_db"])), {
        "csv": str(csv_path.relative_to(ROOT)).replace("\\", "/"),
        "csv_sha256": _sha256(csv_path),
        "resolved_config": str(config_path.relative_to(ROOT)).replace("\\", "/"),
        "resolved_config_sha256": _sha256(config_path),
        "points": len(rows),
        "trials": sum(int(row["trials"]) for row in rows),
        "errors": sum(int(row["errors"]) for row in rows),
    }


def main() -> None:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    from PIL import Image

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exclude-b0", action="store_true", help="omit C300_B0_QC from the figure")
    args = parser.parse_args()
    active_curves = {
        candidate_id: spec
        for candidate_id, spec in CURVES.items()
        if not (args.exclude_b0 and candidate_id == "C300_B0_QC")
    }
    figure = ANALYSIS / (
        "c300_4tx_2rx_estimated_bler_al2_al4_no_b0qc.png"
        if args.exclude_b0
        else "c300_4tx_2rx_estimated_bler_al2_al4.png"
    )
    artifact_suffix = "_no_b0qc" if args.exclude_b0 else ""

    ANALYSIS.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict[str, str]] = []
    receipt: dict[str, object] = {"scope": "4Tx/2Rx, two-symbol, estimated CSI, AL2/AL4", "sources": {}}
    loaded: dict[tuple[int, str], list[dict[str, str]]] = {}
    for al in (2, 4):
        for candidate_id, spec in active_curves.items():
            rows, source = _load_curve(al, candidate_id, spec)
            loaded[(al, candidate_id)] = rows
            receipt["sources"][f"al{al}:{candidate_id}"] = source
            for row in rows:
                all_rows.append({**row, "plot_label": str(spec["label"]), "source_csv": source["csv"]})

    fieldnames = list(dict.fromkeys(key for row in all_rows for key in row))
    with (ANALYSIS / f"plot_points{artifact_suffix}.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)
    (ANALYSIS / f"source_receipt{artifact_suffix}.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (ANALYSIS / f"curve_styles{artifact_suffix}.json").write_text(
        json.dumps({key: {k: v for k, v in value.items() if k != "root"} for key, value in active_curves.items()}, indent=2) + "\n",
        encoding="utf-8",
    )

    fig, axes = plt.subplots(1, 2, figsize=(18.5, 8.5), sharey=True)
    handles = []
    labels = []
    for axis, al in zip(axes, (2, 4), strict=True):
        for candidate_id, spec in active_curves.items():
            rows = loaded[(al, candidate_id)]
            x = np.asarray([float(row["snr_db"]) for row in rows])
            trials = np.asarray([int(row["trials"]) for row in rows])
            y_raw = np.asarray([float(row["bler"]) for row in rows])
            y = np.maximum(y_raw, 0.5 / trials)
            (line,) = axis.plot(
                x,
                y,
                label=spec["label"],
                color=spec["color"],
                linestyle=spec["linestyle"],
                marker=spec["marker"],
                linewidth=2.5,
                markersize=8,
            )
            if al == 2:
                handles.append(line)
                labels.append(str(spec["label"]))
        if not args.exclude_b0:
            axis.axhline(0.10, color="#555555", linestyle="-.", linewidth=1.8)
        axis.axhline(0.01, color="#555555", linestyle=(0, (2, 2)), linewidth=1.8)
        axis.set_yscale("log")
        axis.set_ylim(1e-4, 1.0)
        axis.set_xlabel("SNR (dB)", fontsize=16)
        axis.set_title(f"AL{al}", fontsize=18)
        axis.xaxis.set_major_locator(MultipleLocator(0.5))
        axis.xaxis.set_minor_locator(MultipleLocator(0.25))
        axis.tick_params(labelsize=14)
        axis.grid(True, which="major", axis="x", alpha=0.42, linewidth=0.9)
        axis.grid(True, which="minor", axis="x", alpha=0.22, linewidth=0.55)
        axis.grid(True, which="both", axis="y", alpha=0.35)
        if not args.exclude_b0:
            axis.text(0.99, 0.105, "10%", transform=axis.get_yaxis_transform(), ha="right", va="bottom", fontsize=14)
        axis.text(0.99, 0.0105, "1%", transform=axis.get_yaxis_transform(), ha="right", va="bottom", fontsize=14)
    axes[0].set_ylabel("Estimated-CSI DCI BLER", fontsize=16)
    fig.suptitle("C300 4Tx/2Rx two-symbol PDCCH BLER", fontsize=19)
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=16, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.13, 1, 0.95))
    fig.savefig(figure, dpi=180, bbox_inches="tight")
    plt.close(fig)

    preview_dir = ANALYSIS / "previews_13cm"
    preview_dir.mkdir(parents=True, exist_ok=True)
    with Image.open(figure) as source:
        width = 1535
        height = round(source.height * width / source.width)
        source.resize((width, height), Image.Resampling.LANCZOS).save(preview_dir / figure.name)
    with Image.open(figure) as source:
        print(json.dumps({"figure": str(figure), "size_px": list(source.size), "points": len(all_rows)}, indent=2))


if __name__ == "__main__":
    main()
