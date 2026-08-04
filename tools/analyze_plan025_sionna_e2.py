"""Combine plan-025 QC/Sidon known-delay runs and fit 10%/1% BLER targets."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
from scipy.optimize import minimize
from scipy.stats import binomtest


IDS = ("QC", "Sidon")
FLAT_RESULT024_GAIN_DB = {"bler10": 0.33, "bler1": 0.85}


def read_rows(directory: Path) -> List[Dict[str, str]]:
    path = directory / "sidon_qc_bler.csv"
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No rows in {path}")
    if {row["id"] for row in rows} != set(IDS):
        raise ValueError(f"Expected only QC and Sidon rows in {path}")
    if any(str(row.get("branch", "known")) != "known" for row in rows):
        raise ValueError(f"Only known-delay rows are allowed in {path}")
    return rows


def write_rows(rows: Sequence[Dict[str, object]], path: Path) -> None:
    fields = sorted(set().union(*(row.keys() for row in rows)))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def logistic_target(rows: Sequence[Dict[str, str]], candidate: str, target: float) -> Dict[str, object]:
    selected = sorted(
        (row for row in rows if row["id"] == candidate),
        key=lambda row: float(row["snr_db"]),
    )
    x = np.asarray([float(row["snr_db"]) for row in selected], dtype=np.float64)
    n = np.asarray([int(row["trials"]) for row in selected], dtype=np.float64)
    y = np.asarray([int(row["tb_errors"]) for row in selected], dtype=np.float64)
    if len(x) < 2 or not np.any(y > 0) or not np.any(y < n):
        raise ValueError(f"Insufficient finite-error points for {candidate} target {target}")
    matrix = np.column_stack((np.ones_like(x), x))

    def negative_log_likelihood(theta):
        linear = matrix @ theta
        return np.sum(n * np.logaddexp(0.0, linear) - y * linear)

    fit = minimize(negative_log_likelihood, [15.0, -1.0], method="BFGS")
    if not fit.success and not np.all(np.isfinite(fit.x)):
        raise RuntimeError(f"Logistic fit failed for {candidate}: {fit.message}")
    intercept, slope = fit.x
    fitted_probability = 1.0 / (1.0 + np.exp(-(matrix @ fit.x)))
    weight = n * fitted_probability * (1.0 - fitted_probability)
    covariance = np.linalg.inv(matrix.T @ (weight[:, None] * matrix))
    logit_target = math.log(target / (1.0 - target))
    snr = (logit_target - intercept) / slope
    gradient = np.asarray([-1.0 / slope, -snr / slope])
    standard_error = float(np.sqrt(gradient @ covariance @ gradient))
    near_errors = sum(
        int(row["tb_errors"]) for row in selected
        if 0.5 * target <= float(row["bler"]) <= 2.0 * target
    )
    return {
        "snr_db": float(snr),
        "standard_error_db": standard_error,
        "ci95_lo_db": float(snr - 1.96 * standard_error),
        "ci95_hi_db": float(snr + 1.96 * standard_error),
        "logit_intercept": float(intercept),
        "logit_slope_per_db": float(slope),
        "points": len(selected),
        "total_trials": int(np.sum(n)),
        "total_errors": int(np.sum(y)),
        "errors_in_half_to_twice_target_band": int(near_errors),
        "sufficient_for_deterministic_claim": bool(near_errors >= 30),
    }


def target_summary(rows: Sequence[Dict[str, str]], label: str, target: float) -> Dict[str, object]:
    candidates = {candidate: logistic_target(rows, candidate, target) for candidate in IDS}
    gain = float(candidates["QC"]["snr_db"]) - float(candidates["Sidon"]["snr_db"])
    gain_standard_error = math.sqrt(
        float(candidates["QC"]["standard_error_db"]) ** 2
        + float(candidates["Sidon"]["standard_error_db"]) ** 2
    )
    return {
        "target_probability": target,
        "candidates": candidates,
        "sidon_gain_vs_qc_db": gain,
        "sidon_gain_standard_error_independent_db": gain_standard_error,
        "sidon_gain_ci95_conservative_lo_db": gain - 1.96 * gain_standard_error,
        "sidon_gain_ci95_conservative_hi_db": gain + 1.96 * gain_standard_error,
        "result024_flat_gain_db": FLAT_RESULT024_GAIN_DB[label],
        "gain_difference_vs_result024_db": gain - FLAT_RESULT024_GAIN_DB[label],
        "confidence_note": (
            "QC/Sidon trials are paired, but the conservative interval ignores positive pairing covariance."
        ),
    }


def ce_fairness(rows: Sequence[Dict[str, str]]) -> Dict[str, object]:
    differences = []
    for snr in sorted({float(row["snr_db"]) for row in rows}):
        values = {
            row["id"]: float(row["ce_nmse_mean_dB"])
            for row in rows if float(row["snr_db"]) == snr
        }
        if set(values) == set(IDS):
            differences.append({
                "snr_db": snr,
                "sidon_minus_qc_nmse_db": values["Sidon"] - values["QC"],
            })
    return {
        "max_absolute_nmse_difference_db": max(
            abs(float(row["sidon_minus_qc_nmse_db"])) for row in differences
        ),
        "by_snr": differences,
    }


def paired_error_analysis(directories: Sequence[Path]) -> Dict[str, object]:
    rows = []
    for directory in directories:
        path = directory / "paired_error_counts.csv"
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows.extend(csv.DictReader(handle))

    by_snr = []
    total_qc_only = 0
    total_sidon_only = 0
    for row in sorted(rows, key=lambda item: float(item["snr_db"])):
        qc_only = int(row["qc_only_error"])
        sidon_only = int(row["sidon_only_error"])
        discordant = qc_only + sidon_only
        p_value = float(
            binomtest(min(qc_only, sidon_only), discordant, p=0.5).pvalue
        ) if discordant else 1.0
        by_snr.append({
            "snr_db": float(row["snr_db"]),
            "qc_only_error": qc_only,
            "sidon_only_error": sidon_only,
            "discordant_pairs": discordant,
            "mcnemar_exact_two_sided_p": p_value,
        })
        total_qc_only += qc_only
        total_sidon_only += sidon_only

    total_discordant = total_qc_only + total_sidon_only
    total_p_value = float(
        binomtest(
            min(total_qc_only, total_sidon_only), total_discordant, p=0.5
        ).pvalue
    ) if total_discordant else 1.0
    return {
        "aggregate_qc_only_error": total_qc_only,
        "aggregate_sidon_only_error": total_sidon_only,
        "aggregate_discordant_pairs": total_discordant,
        "aggregate_mcnemar_exact_two_sided_p": total_p_value,
        "aggregation_note": (
            "The aggregate test pools different SNR points and is supporting evidence; "
            "target-SNR fits remain the primary comparison."
        ),
        "by_snr": by_snr,
    }


def plot(refine10, refine1, summary, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"QC": "#f77f00", "Sidon": "#111111"}
    markers = {"QC": "o", "Sidon": "s"}
    figure, axes = plt.subplots(1, 2, figsize=(12.4, 5.2))
    for axis, rows, label, target in zip(
        axes, (refine10, refine1), ("bler10", "bler1"), (0.10, 0.01)
    ):
        for candidate in IDS:
            selected = sorted(
                (row for row in rows if row["id"] == candidate),
                key=lambda row: float(row["snr_db"]),
            )
            x = np.asarray([float(row["snr_db"]) for row in selected])
            y = np.asarray([max(float(row["bler"]), 0.5 / int(row["trials"])) for row in selected])
            low = np.asarray([max(float(row["bler_wilson95_lo"]), 1e-5) for row in selected])
            high = np.asarray([float(row["bler_wilson95_hi"]) for row in selected])
            fit = summary["targets"][label]["candidates"][candidate]
            xx = np.linspace(float(np.min(x)), float(np.max(x)), 200)
            probability = 1.0 / (
                1.0 + np.exp(-(float(fit["logit_intercept"]) + float(fit["logit_slope_per_db"]) * xx))
            )
            axis.semilogy(x, y, marker=markers[candidate], ls="none", color=colors[candidate], label=candidate)
            axis.vlines(x, low, high, color=colors[candidate], alpha=0.45)
            axis.semilogy(xx, probability, color=colors[candidate])
        target_result = summary["targets"][label]
        axis.axhline(target, color="#555555", ls="--", lw=1)
        axis.set_title(
            f"{target*100:g}% BLER: Sidon gain {target_result['sidon_gain_vs_qc_db']:.2f} dB\n"
            f"95% conservative CI [{target_result['sidon_gain_ci95_conservative_lo_db']:.2f}, "
            f"{target_result['sidon_gain_ci95_conservative_hi_db']:.2f}] dB"
        )
        axis.set_xlabel("average receive SNR (dB)")
        axis.set_ylabel("estimated-CSI TB BLER")
        axis.grid(alpha=0.3, which="both")
        axis.legend()
    figure.suptitle("Plan-025: QC vs Sidon, Sionna TDL-A 5 ns, known-delay 2D RMMSE")
    figure.tight_layout()
    figure.savefig(path, dpi=170)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prescan", type=Path, required=True)
    parser.add_argument("--refine10", type=Path, required=True)
    parser.add_argument("--refine1", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    prescan = read_rows(args.prescan)
    refine10 = read_rows(args.refine10)
    refine1 = read_rows(args.refine1)
    write_rows([*refine10, *refine1], args.out / "sidon_qc_refined_combined.csv")
    summary = {
        "model": "Sionna 1.0.2 TDL-A, 5 ns, 0 km/h, 3.5 GHz",
        "receiver": "known-delay matched/oracle 2D RMMSE",
        "targets": {
            "bler10": target_summary(refine10, "bler10", 0.10),
            "bler1": target_summary(refine1, "bler1", 0.01),
        },
        "ce_fairness": ce_fairness([*refine10, *refine1]),
        "paired_error_analysis": paired_error_analysis([args.refine10, args.refine1]),
        "prescan": {
            "snrs_db": sorted({float(row["snr_db"]) for row in prescan}),
            "trials_per_point": sorted({int(row["trials"]) for row in prescan}),
        },
        "evidence_paths": {
            "prescan": str(args.prescan),
            "refine10": str(args.refine10),
            "refine1": str(args.refine1),
        },
    }
    (args.out / "final_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    plot(refine10, refine1, summary, args.out / "sidon_qc_tdl_bler.png")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
