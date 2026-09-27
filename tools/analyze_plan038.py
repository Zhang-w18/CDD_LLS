"""Audit and plot all plan-038 BLER and CE-NMSE curves in one figure."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import yaml

from run_plan038_refinement import WINDOWS, _quarter_db_grid


ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "outputs/experiment038_pdcch_a40/20260922_4tx_2rx_2sym"
ANALYSIS = RUN_ROOT / "analysis"
CURVES = {
    "SIDON_MATCHED": ("Sidon non-transparent", "#1f77b4", "-", "s"),
    "SIDON_TRANSPARENT": ("Sidon transparent", "#1f77b4", ":", "s"),
    "B0QC_MATCHED": ("B0 QC non-transparent", "#111111", "-", "o"),
    "B0QC_TRANSPARENT": ("B0 QC transparent", "#111111", ":", "o"),
    "CDD911_MATCHED": ("CDD911 non-transparent", "#d62728", "-", "D"),
    "CDD911_TRANSPARENT": ("CDD911 transparent", "#d62728", ":", "D"),
    "CDD130_MATCHED": ("CDD130 non-transparent", "#ff7f0e", "-", "X"),
    "CDD130_TRANSPARENT": ("CDD130 transparent", "#ff7f0e", ":", "X"),
    "PRG_DFT4_TRANSPARENT": ("Precoder cycling transparent", "#2ca02c", "--", ">"),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _audit_rows(rows: list[dict[str, str]], al: int, candidate_id: str, formal: bool) -> None:
    if not rows:
        raise ValueError(f"No rows for AL{al} {candidate_id}")
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
            if row.get(field) != value:
                raise ValueError(f"AL{al} {candidate_id}: {field}={row.get(field)!r}, expected {value!r}")
        if formal:
            trials, errors = int(row["trials"]), int(row["errors"])
            if trials < 10000 or (trials < 50000 and errors < 200):
                raise ValueError(f"Incomplete formal point: AL{al} {candidate_id} {row['snr_db']} dB")


def _load_curve(al: int, candidate_id: str) -> tuple[list[dict[str, str]], dict[str, object]]:
    path = RUN_ROOT / "formal" / f"al{al}" / candidate_id / "bler_points.csv"
    all_rows = _read_csv(path)
    _audit_rows(all_rows, al, candidate_id, formal=True)
    start, stop = WINDOWS[(al, candidate_id)]
    rows = [row for row in all_rows if start - 1e-12 <= float(row["snr_db"]) <= stop + 1e-12]
    actual = sorted(float(row["snr_db"]) for row in rows)
    expected = _quarter_db_grid(start, stop)
    if actual != expected:
        missing = sorted(set(expected) - set(actual))
        extra = sorted(set(actual) - set(expected))
        raise ValueError(f"Incomplete uniform plot grid AL{al} {candidate_id}: missing={missing}, extra={extra}")
    source = {
        "stage": "formal-refined",
        "path": str(path.relative_to(ROOT)).replace("\\", "/"),
        "sha256": _sha256(path),
        "plot_window_db": [start, stop],
        "spacing_db": 0.25,
        "points": len(rows),
    }
    return sorted(rows, key=lambda row: float(row["snr_db"])), source


def _audit_config(al: int) -> None:
    path = ROOT / f"configs/pdcch_result038_al{al}_formal.yaml"
    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if config["antenna"] != {"n_tx": 4, "n_rx": 2}:
        raise ValueError(f"Unexpected antenna configuration: {path}")
    if config["resource"]["duration_symbols"] != 2 or config["resource"]["aggregation_level"] != al:
        raise ValueError(f"Unexpected resource configuration: {path}")
    if config["pdcch"]["payload_bits"] != 40 or config["pdcch"]["coded_bits"] != 108 * al:
        raise ValueError(f"Unexpected payload/coded bits: {path}")
    if {item["candidate_id"] for item in config["candidates"]} != set(CURVES):
        raise ValueError(f"Unexpected candidate whitelist: {path}")


def _save_png(fig: object, path: Path) -> None:
    """Save through an open handle so matplotlib never reparses a Windows path string."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        fig.savefig(handle, format="png", dpi=180, bbox_inches="tight")


def main() -> None:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator

    ANALYSIS.mkdir(parents=True, exist_ok=True)
    loaded: dict[tuple[int, str], list[dict[str, str]]] = {}
    receipt: dict[str, object] = {"scope": "4Tx/2Rx, two-symbol PDCCH, A=40, estimated CSI", "sources": {}}
    flat_rows: list[dict[str, str]] = []
    for al in (1, 2, 4):
        _audit_config(al)
        for candidate_id in CURVES:
            rows, source = _load_curve(al, candidate_id)
            loaded[(al, candidate_id)] = rows
            receipt["sources"][f"al{al}:{candidate_id}"] = source
            for row in rows:
                flat_rows.append({**row, "plot_label": CURVES[candidate_id][0]})

    fields = list(dict.fromkeys(key for row in flat_rows for key in row))
    with (ANALYSIS / "plot_points.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(flat_rows)
    (ANALYSIS / "source_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")

    fig, axes = plt.subplots(2, 3, figsize=(21, 12), sharex="col")
    handles, labels = [], []
    for column, al in enumerate((1, 2, 4)):
        for candidate_id, (label, color, linestyle, marker) in CURVES.items():
            rows = loaded[(al, candidate_id)]
            x = np.asarray([float(row["snr_db"]) for row in rows])
            trials = np.asarray([int(row["trials"]) for row in rows])
            bler = np.maximum(np.asarray([float(row["bler"]) for row in rows]), 0.5 / trials)
            nmse = np.asarray([float(row["ce_nmse_db"]) for row in rows])
            line, = axes[0, column].plot(x, bler, color=color, linestyle=linestyle, marker=marker, linewidth=1.7, markersize=4, label=label)
            axes[1, column].plot(x, nmse, color=color, linestyle=linestyle, marker=marker, linewidth=1.7, markersize=4)
            if column == 0:
                handles.append(line)
                labels.append(label)
        axes[0, column].set_yscale("log")
        axes[0, column].set_ylim(1e-4, 1.0)
        axes[0, column].axhline(0.01, color="#777777", linewidth=0.9, linestyle=(0, (2, 2)))
        axes[0, column].set_title(f"AL{al}")
        axes[1, column].set_xlabel("SNR (dB)")
        for row in range(2):
            axes[row, column].grid(True, which="both", alpha=0.28)
    axes[0, 0].set_ylabel("DCI BLER")
    axes[1, 0].set_ylabel("Data-RE CE NMSE (dB)")
    fig.suptitle("Plan-038: 4Tx/2Rx, two-symbol PDCCH, A=40")
    fig.legend(handles, labels, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.01))
    fig.tight_layout(rect=(0, 0.10, 1, 0.96))
    figure = ANALYSIS / "plan038_all_al_bler_nmse_updated_v2.png"
    _save_png(fig, figure)
    plt.close(fig)

    focused_ids = ("SIDON_MATCHED", "SIDON_TRANSPARENT", "PRG_DFT4_TRANSPARENT")
    al_colors = {1: "#1f77b4", 2: "#d62728", 4: "#2ca02c"}
    scheme_markers = {
        "SIDON_MATCHED": "s",
        "SIDON_TRANSPARENT": "o",
        "PRG_DFT4_TRANSPARENT": ">",
    }
    focused_fig, focused_axis = plt.subplots(1, 1, figsize=(11.5, 7.0))
    for al in (1, 2, 4):
        for candidate_id in focused_ids:
            if al == 1 and candidate_id == "SIDON_TRANSPARENT":
                continue
            label = {
                "SIDON_MATCHED": "CDD non-transparent",
                "SIDON_TRANSPARENT": "CDD transparent",
                "PRG_DFT4_TRANSPARENT": "Precoder cycling transparent",
            }[candidate_id]
            if al == 1 and candidate_id == "PRG_DFT4_TRANSPARENT":
                label = "Fixed precoder"
            rows = loaded[(al, candidate_id)]
            x = np.asarray([float(row["snr_db"]) for row in rows])
            trials = np.asarray([int(row["trials"]) for row in rows])
            bler = np.maximum(np.asarray([float(row["bler"]) for row in rows]), 0.5 / trials)
            focused_axis.plot(
                x,
                bler,
                color=al_colors[al],
                linestyle="--" if candidate_id == "SIDON_TRANSPARENT" else "-",
                marker=scheme_markers[candidate_id],
                linewidth=2.0,
                markersize=8,
                label=f"AL{al} — {label}",
            )
    focused_axis.set_yscale("log")
    focused_axis.set_ylim(1e-4, 1.0)
    focused_axis.axhline(0.01, color="#777777", linewidth=0.9, linestyle=(0, (2, 2)))
    focused_axis.set_xlabel("SNR (dB)", fontsize=15)
    focused_axis.set_ylabel("DCI BLER", fontsize=15)
    focused_axis.set_title("Sidon and precoder cycling: AL1/AL2/AL4", fontsize=16, pad=16)
    focused_axis.tick_params(axis="both", which="both", labelsize=13)
    focused_axis.grid(True, which="both", alpha=0.28)
    focused_axis.xaxis.set_minor_locator(AutoMinorLocator(2))
    focused_axis.grid(True, axis="x", which="minor", alpha=0.16)
    focused_axis.legend(loc="best", ncol=1, fontsize=12)
    focused_fig.tight_layout()
    focused_figure = ANALYSIS / "plan038_sidon_cycling.png"
    _save_png(focused_fig, focused_figure)
    plt.close(focused_fig)

    cdd911_ids = ("CDD911_MATCHED", "CDD911_TRANSPARENT", "PRG_DFT4_TRANSPARENT")
    cdd911_markers = {
        "CDD911_MATCHED": "D",
        "CDD911_TRANSPARENT": "o",
        "PRG_DFT4_TRANSPARENT": ">",
    }
    cdd911_fig, cdd911_axis = plt.subplots(1, 1, figsize=(11.5, 7.0))
    for al in (1, 2, 4):
        for candidate_id in cdd911_ids:
            label = {
                "CDD911_MATCHED": "CDD non-transparent",
                "CDD911_TRANSPARENT": "CDD transparent",
                "PRG_DFT4_TRANSPARENT": "Precoder cycling transparent",
            }[candidate_id]
            if al == 1 and candidate_id == "PRG_DFT4_TRANSPARENT":
                label = "Fixed precoder"
            rows = loaded[(al, candidate_id)]
            x = np.asarray([float(row["snr_db"]) for row in rows])
            trials = np.asarray([int(row["trials"]) for row in rows])
            bler = np.maximum(np.asarray([float(row["bler"]) for row in rows]), 0.5 / trials)
            cdd911_axis.plot(
                x,
                bler,
                color=al_colors[al],
                linestyle="--" if candidate_id == "CDD911_TRANSPARENT" else "-",
                marker=cdd911_markers[candidate_id],
                linewidth=2.0,
                markersize=8,
                label=f"AL{al} — {label}",
            )
    cdd911_axis.set_yscale("log")
    cdd911_axis.set_ylim(1e-4, 1.0)
    cdd911_axis.axhline(0.01, color="#777777", linewidth=0.9, linestyle=(0, (2, 2)))
    cdd911_axis.set_xlabel("SNR (dB)", fontsize=15)
    cdd911_axis.set_ylabel("DCI BLER", fontsize=15)
    cdd911_axis.set_title("CDD911 and precoder cycling: AL1/AL2/AL4", fontsize=16, pad=16)
    cdd911_axis.tick_params(axis="both", which="both", labelsize=13)
    cdd911_axis.grid(True, which="both", alpha=0.28)
    cdd911_axis.xaxis.set_minor_locator(AutoMinorLocator(2))
    cdd911_axis.grid(True, axis="x", which="minor", alpha=0.16)
    cdd911_axis.legend(loc="upper right", ncol=2, fontsize=12)
    cdd911_fig.tight_layout()
    cdd911_figure = ANALYSIS / "plan038_cdd911_cycling_al1_al2_al4_single_axis.png"
    _save_png(cdd911_fig, cdd911_figure)
    plt.close(cdd911_fig)
    print(
        json.dumps(
            {
                "all_curves_figure": str(figure.relative_to(ROOT)),
                "sidon_cycling_figure": str(focused_figure.relative_to(ROOT)),
                "cdd911_cycling_figure": str(cdd911_figure.relative_to(ROOT)),
                "curves": len(loaded),
                "rows": len(flat_rows),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
