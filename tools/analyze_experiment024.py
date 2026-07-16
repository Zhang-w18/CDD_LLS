"""Analyze refined Experiment 024 link results and produce final figures."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "experiment024_segment_sidon_qc" / "20260716_main"
IDS = ("QC_arith_s9", "Sidon")


def read(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write(rows, path: Path):
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def logistic_target(rows, cid: str, target: float):
    rr = sorted((r for r in rows if r["id"] == cid), key=lambda r: float(r["snr_db"]))
    x = np.array([float(r["snr_db"]) for r in rr])
    n = np.array([int(r["trials"]) for r in rr])
    y = np.array([int(r["tb_errors"]) for r in rr])
    X = np.column_stack([np.ones_like(x), x])

    def nll(theta):
        z = X @ theta
        return np.sum(n * np.logaddexp(0.0, z) - y * z)

    opt = minimize(nll, [15.0, -1.1], method="BFGS")
    a, b = opt.x
    prob = 1.0 / (1.0 + np.exp(-(X @ opt.x)))
    weight = n * prob * (1.0 - prob)
    cov = np.linalg.inv(X.T @ (weight[:, None] * X))
    level = math.log(target / (1.0 - target))
    snr = (level - a) / b
    grad = np.array([-1.0 / b, -snr / b])
    se = float(np.sqrt(grad @ cov @ grad))
    return {
        "snr_dB": float(snr), "se_dB": se,
        "ci95_lo_dB": float(snr - 1.96 * se), "ci95_hi_dB": float(snr + 1.96 * se),
        "logit_intercept": float(a), "logit_slope_per_dB": float(b),
        "points": len(rr), "total_trials": int(n.sum()), "total_errors": int(y.sum()),
    }


def analyze():
    ten = read(OUT / "refine_10pct" / "sidon_qc_bler.csv") + read(OUT / "refine_10pct_14p75" / "sidon_qc_bler.csv")
    one = read(OUT / "refine_1pct" / "sidon_qc_bler.csv")
    combined = sorted(ten + one, key=lambda r: (float(r["snr_db"]), r["id"]))
    write(combined, OUT / "sidon_qc_bler_refined.csv")
    result = {"method": "binomial-logit fit versus SNR; delta-method 95% CI", "targets": {}}
    for label, target, rows in (("bler10", 0.10, ten), ("bler1", 0.01, one)):
        fits = {cid: logistic_target(rows, cid, target) for cid in IDS}
        gain = fits["QC_arith_s9"]["snr_dB"] - fits["Sidon"]["snr_dB"]
        # Candidate errors are paired, but only aggregate counts were retained. Ignoring
        # positive pairing covariance is conservative for the difference uncertainty.
        gain_se = math.sqrt(fits["QC_arith_s9"]["se_dB"] ** 2 + fits["Sidon"]["se_dB"] ** 2)
        result["targets"][label] = {
            "target_probability": target,
            "candidates": fits,
            "sidon_gain_dB": gain,
            "sidon_gain_se_independent_dB": gain_se,
            "sidon_gain_ci95_conservative_lo_dB": gain - 1.96 * gain_se,
            "sidon_gain_ci95_conservative_hi_dB": gain + 1.96 * gain_se,
        }
    ce_diffs = []
    for snr in sorted({float(r["snr_db"]) for r in combined}):
        d = {r["id"]: float(r["ce_nmse_mean_dB"]) for r in combined if float(r["snr_db"]) == snr}
        if len(d) == 2:
            ce_diffs.append({"snr_db": snr, "sidon_minus_qc_nmse_dB": d["Sidon"] - d["QC_arith_s9"]})
    result["ce_fairness"] = {
        "max_abs_nmse_difference_dB": max(abs(r["sidon_minus_qc_nmse_dB"]) for r in ce_diffs),
        "by_snr": ce_diffs,
    }
    result["sidon_check"] = {
        "j_list": [0, 1, 3, 7, 12, 20, 30, 65],
        "unordered_pair_sums": 36,
        "unique_integer_pair_sums": 36,
        "unique_pair_sums_mod_576": 36,
    }
    (OUT / "final_link_summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    plot(ten, one, result, OUT / "figures" / "sidon_qc_bler_refined.png")
    return result


def plot(ten, one, result, path: Path):
    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.2), sharey=False)
    colors = {"QC_arith_s9": "#f77f00", "Sidon": "#111111"}
    markers = {"QC_arith_s9": "o", "Sidon": "s"}
    for ax, rows, label, target in zip(axes, (ten, one), ("bler10", "bler1"), (0.10, 0.01)):
        for cid in IDS:
            rr = sorted((r for r in rows if r["id"] == cid), key=lambda r: float(r["snr_db"]))
            x = np.array([float(r["snr_db"]) for r in rr])
            y = np.array([float(r["bler"]) for r in rr])
            lo = np.array([float(r["bler_wilson95_lo"]) for r in rr])
            hi = np.array([float(r["bler_wilson95_hi"]) for r in rr])
            ax.semilogy(x, y, marker=markers[cid], ls="none", color=colors[cid], label=cid)
            ax.vlines(x, np.maximum(lo, 1e-4), hi, color=colors[cid], alpha=0.45, lw=1.2)
            fit = result["targets"][label]["candidates"][cid]
            xx = np.linspace(x.min(), x.max(), 150)
            pp = 1.0 / (1.0 + np.exp(-(fit["logit_intercept"] + fit["logit_slope_per_dB"] * xx)))
            ax.semilogy(xx, pp, "-", color=colors[cid], lw=1.5)
            ax.axvline(fit["snr_dB"], color=colors[cid], ls=":", lw=1)
        gain = result["targets"][label]["sidon_gain_dB"]
        ci_lo = result["targets"][label]["sidon_gain_ci95_conservative_lo_dB"]
        ci_hi = result["targets"][label]["sidon_gain_ci95_conservative_hi_dB"]
        ax.axhline(target, color="#555555", ls="--", lw=1)
        ax.set_title(f"{target*100:g}% BLER: Sidon gain {gain:.2f} dB\nconservative 95% CI [{ci_lo:.2f}, {ci_hi:.2f}] dB")
        ax.set_xlabel("average receive SNR (dB)")
        ax.set_ylabel("estimated-CSI TB BLER")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
    fig.suptitle("Experiment 024 E2: Sidon vs QC, V-aware matched LMMSE, 48 PRB flat channel")
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2))
