"""Plot result-031 effective PDPs and their frequency correlation."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PDP_RECEIPT = (
    ROOT
    / "outputs"
    / "experiment031_pdcch_cdd"
    / "20260907_c300_4tx_2sym"
    / "pdp_t0p01_receipt.json"
)
OUTPUT = (
    ROOT
    / "outputs"
    / "experiment031_pdcch_cdd"
    / "20260914_c300_al1_cdd911"
    / "analysis"
    / "pdp_comparison"
)
SUBCARRIER_SPACING_HZ = 30e3
ACTIVE_SUBCARRIERS = 36
DMRS_SUBCARRIERS = np.asarray([1, 5, 9, 13, 17, 21, 25, 29, 33], dtype=int)
Q_NS = 1e9 / (ACTIVE_SUBCARRIERS * SUBCARRIER_SPACING_HZ)
CDD_DELAY_SETS_NS = (
    ("CDD [0, 0, 130, 130] ns", np.asarray([0.0, 0.0, 130.0, 130.0])),
    ("CDD [0, 0, 911, 911] ns", np.asarray([0.0, 0.0, 911.0, 911.0])),
    ("Sidon [0, 1, 3, 7]q", Q_NS * np.asarray([0.0, 1.0, 3.0, 7.0])),
)


def rms_delay_spread_ns(delays_ns: np.ndarray, powers: np.ndarray) -> tuple[float, float]:
    powers = powers / np.sum(powers)
    mean_ns = float(np.sum(delays_ns * powers))
    rms_ns = float(np.sqrt(np.sum(powers * (delays_ns - mean_ns) ** 2)))
    return mean_ns, rms_ns


def effective_pdp(
    physical_delays_ns: np.ndarray,
    physical_powers: np.ndarray,
    cdd_delays_ns: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    # Equal-power Tx branches produce one shifted physical PDP per CDD delay.
    delays_ns = np.concatenate(
        [physical_delays_ns + cdd_delay_ns for cdd_delay_ns in cdd_delays_ns]
    )
    powers = np.tile(physical_powers / len(cdd_delays_ns), len(cdd_delays_ns))
    order = np.argsort(delays_ns, kind="stable")
    return delays_ns[order], powers[order]


def frequency_correlation(
    delays_ns: np.ndarray,
    powers: np.ndarray,
    lags: np.ndarray,
) -> np.ndarray:
    delays_s = delays_ns * 1e-9
    phase = -2j * np.pi * lags[:, None] * SUBCARRIER_SPACING_HZ * delays_s[None, :]
    correlation = np.exp(phase) @ powers
    return correlation / correlation[0]


def main() -> None:
    receipt = json.loads(PDP_RECEIPT.read_text(encoding="utf-8"))
    physical_delays_ns = np.asarray(receipt["tap_delays_ns"], dtype=np.float64)
    physical_powers = np.asarray(receipt["normalized_tap_powers"], dtype=np.float64)
    physical_powers /= np.sum(physical_powers)

    scenarios: list[tuple[str, np.ndarray, np.ndarray]] = [
        ("Original TDL-C 300 ns", physical_delays_ns, physical_powers)
    ]
    cdd_delay_sets: dict[str, list[float]] = {}
    for name, cdd_delays_ns in CDD_DELAY_SETS_NS:
        delays_ns, powers = effective_pdp(
            physical_delays_ns,
            physical_powers,
            cdd_delays_ns,
        )
        scenarios.append((name, delays_ns, powers))
        cdd_delay_sets[name] = cdd_delays_ns.tolist()

    OUTPUT.mkdir(parents=True, exist_ok=True)

    summary: list[dict[str, float | str | int]] = []
    with (OUTPUT / "effective_pdp_taps.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "scenario",
                "tap_index",
                "delay_ns",
                "power_linear",
                "power_db_total_normalized",
                "power_db_peak_normalized",
            ),
        )
        writer.writeheader()
        for name, delays_ns, powers in scenarios:
            mean_ns, rms_ns = rms_delay_spread_ns(delays_ns, powers)
            summary.append(
                {
                    "scenario": name,
                    "tap_count": int(len(delays_ns)),
                    "mean_delay_ns": mean_ns,
                    "rms_delay_spread_ns": rms_ns,
                }
            )
            for tap_index, (delay_ns, power) in enumerate(zip(delays_ns, powers, strict=True)):
                writer.writerow(
                    {
                        "scenario": name,
                        "tap_index": tap_index,
                        "delay_ns": f"{delay_ns:.12g}",
                        "power_linear": f"{power:.17g}",
                        "power_db_total_normalized": f"{10.0 * np.log10(power):.12g}",
                        "power_db_peak_normalized": f"{10.0 * np.log10(power / np.max(powers)):.12g}",
                    }
                )

    (OUTPUT / "effective_pdp_summary.json").write_text(
        json.dumps(
            {
                "source_receipt": str(PDP_RECEIPT.relative_to(ROOT)).replace("\\", "/"),
                "subcarrier_spacing_hz": SUBCARRIER_SPACING_HZ,
                "active_subcarriers": ACTIVE_SUBCARRIERS,
                "aggregation_level": 1,
                "ofdm_symbols": 2,
                "dmrs_subcarrier_indices_per_symbol": DMRS_SUBCARRIERS.tolist(),
                "dmrs_re_count": int(2 * len(DMRS_SUBCARRIERS)),
                "q_ns": Q_NS,
                "cdd_model": "equal-power Tx branches; effective PDP is the average of shifted physical PDPs",
                "cdd_delay_sets_ns": cdd_delay_sets,
                "scenarios": summary,
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
        }
    )
    colors = ("#333333", "#0072B2", "#D55E00", "#009E73")
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.2), sharey=True)
    for ax, color, (name, delays_ns, powers), metrics in zip(
        axes.flat, colors, scenarios, summary, strict=True
    ):
        powers_db = 10.0 * np.log10(powers / np.max(powers))
        floor_db = -25.0
        ax.vlines(delays_ns * 1e-3, floor_db, powers_db, color=color, linewidth=0.9, alpha=0.85)
        ax.plot(delays_ns * 1e-3, powers_db, "o", color=color, markersize=3.2)
        ax.set_title(
            f"{name}   (RMS delay spread = {float(metrics['rms_delay_spread_ns']):.1f} ns)",
            loc="left",
        )
        ax.set_ylabel("Tap power relative to peak (dB)")
        ax.set_xlabel("Excess delay (µs)")
        ax.grid(True, which="major", linewidth=0.5, alpha=0.35)
        ax.set_ylim(floor_db, 1.0)
        ax.set_xlim(-0.03, 1.03 * np.max(delays_ns) * 1e-3)

    fig.suptitle("Effective PDP for TDL-C 300 ns (each panel peak-normalized)", fontsize=11)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.97))
    pdp_figure_path = OUTPUT / "tdlc300_effective_pdp_comparison.png"
    fig.savefig(pdp_figure_path, dpi=220, bbox_inches="tight")
    plt.close(fig)

    lags = np.arange(ACTIVE_SUBCARRIERS, dtype=np.float64)
    correlations = {
        name: frequency_correlation(delays_ns, powers, lags)
        for name, delays_ns, powers in scenarios
    }
    with (OUTPUT / "frequency_correlation.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("scenario", "subcarrier_lag", "frequency_separation_khz", "magnitude", "magnitude_squared"),
        )
        writer.writeheader()
        for name, correlation in correlations.items():
            for lag, value in zip(lags.astype(int), correlation, strict=True):
                writer.writerow(
                    {
                        "scenario": name,
                        "subcarrier_lag": lag,
                        "frequency_separation_khz": f"{lag * SUBCARRIER_SPACING_HZ / 1e3:.12g}",
                        "magnitude": f"{abs(value):.17g}",
                        "magnitude_squared": f"{abs(value) ** 2:.17g}",
                    }
                )

    fig, ax = plt.subplots(figsize=(7.4, 4.5))
    for position in DMRS_SUBCARRIERS:
        ax.axvline(position, color="#9467BD", linewidth=0.7, alpha=0.25, zorder=0)
    for color, (name, correlation) in zip(colors, correlations.items(), strict=True):
        ax.plot(lags, np.abs(correlation), marker="o", markersize=2.8, linewidth=1.5, color=color, label=name)
    ax.scatter(
        DMRS_SUBCARRIERS,
        np.full_like(DMRS_SUBCARRIERS, 1.015, dtype=float),
        marker="v",
        s=25,
        color="#9467BD",
        clip_on=False,
        label="DMRS frequency locations (same in both symbols)",
        zorder=5,
    )
    ax.set_xlabel("Local subcarrier index / separation from subcarrier 0 (30 kHz per step)")
    ax.set_ylabel("Normalized correlation magnitude |ρ(0,k)|")
    ax.set_xlim(0, ACTIVE_SUBCARRIERS - 1)
    ax.set_ylim(0, 1.05)
    ax.grid(True, linewidth=0.5, alpha=0.35)
    ax.legend(loc="best")
    ax.set_title("2-symbol AL1 frequency correlation (3 RB; 18 DMRS RE at 9 shared frequencies)")
    fig.tight_layout()
    correlation_figure_path = OUTPUT / "tdlc300_frequency_correlation.png"
    fig.savefig(correlation_figure_path, dpi=220, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(8.2, 7.0), sharex=True, sharey=True)
    image = None
    for ax, (name, correlation) in zip(axes.flat, correlations.items(), strict=True):
        indices = np.arange(ACTIVE_SUBCARRIERS)
        matrix = np.abs(correlation[np.abs(indices[:, None] - indices[None, :])])
        image = ax.imshow(matrix, origin="lower", vmin=0.0, vmax=1.0, cmap="viridis", interpolation="nearest")
        ax.set_title(name, loc="left", fontsize=9)
        ax.set_xlabel("Subcarrier index l")
        ax.set_ylabel("Subcarrier index k")
    assert image is not None
    fig.colorbar(image, ax=axes.ravel().tolist(), shrink=0.82, label="|ρ(k,l)|")
    fig.suptitle("Frequency-correlation magnitude matrices (Toeplitz / frequency-stationary)", fontsize=11)
    fig.subplots_adjust(left=0.08, right=0.88, bottom=0.08, top=0.92, wspace=0.22, hspace=0.25)
    heatmap_path = OUTPUT / "tdlc300_frequency_correlation_heatmaps.png"
    fig.savefig(heatmap_path, dpi=220, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {pdp_figure_path}")
    print(f"Wrote {correlation_figure_path}")
    print(f"Wrote {heatmap_path}")
    for item in summary:
        print(
            f"{item['scenario']}: mean={float(item['mean_delay_ns']):.3f} ns, "
            f"RMS={float(item['rms_delay_spread_ns']):.3f} ns"
        )


if __name__ == "__main__":
    main()
