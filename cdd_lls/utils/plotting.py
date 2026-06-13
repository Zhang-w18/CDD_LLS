from __future__ import annotations

from pathlib import Path
from typing import Dict, List
import math


def _import_pyplot():
    import os
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/cdd_lls_matplotlib")
    os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_bler(rows: List[Dict[str, object]], path: str | Path) -> None:
    if not rows:
        return
    plt = _import_pyplot()
    groups = {}
    for row in rows:
        key = str(row.get("scenario_id", "")) + " | " + str(row.get("variant_id", ""))
        groups.setdefault(key, []).append(row)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for label, items in sorted(groups.items()):
        pts = sorted((float(r["snr_db"]), max(float(r["bler"]), 1e-4)) for r in items)
        ax.semilogy([x for x, _ in pts], [y for _, y in pts], marker="o", linewidth=1.6, label=label)
    ax.grid(True, which="both", alpha=0.35)
    ax.set_xlabel("SNR (dB)")
    ax.set_ylabel("TB BLER")
    ax.set_ylim(1e-4, 1.0)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_nmse(rows: List[Dict[str, object]], path: str | Path) -> None:
    usable = [r for r in rows if math.isfinite(float(r.get("ce_nmse_eff", float("nan"))))]
    if not usable:
        return
    plt = _import_pyplot()
    groups = {}
    for row in usable:
        key = str(row.get("scenario_id", "")) + " | " + str(row.get("variant_id", ""))
        groups.setdefault(key, []).append(row)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for label, items in sorted(groups.items()):
        pts = sorted((float(r["snr_db"]), max(float(r["ce_nmse_eff"]), 1e-8)) for r in items)
        ax.semilogy([x for x, _ in pts], [y for _, y in pts], marker="s", linewidth=1.6, label=label)
    ax.grid(True, which="both", alpha=0.35)
    ax.set_xlabel("SNR (dB)")
    ax.set_ylabel("Equivalent-channel NMSE")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
