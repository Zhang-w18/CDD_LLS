"""Aggregate plan-026 E1/E2/E3 evidence and generate final figures."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_ROOT = ROOT / "outputs" / "experiment026_cdd_design"
DEFAULT_RUN_ID = "20260724_main"


def read_csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(rows: Sequence[Dict[str, object]], path: Path) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted(set().union(*(row.keys() for row in rows)))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def as_bool(value: object) -> bool:
    return str(value).strip().lower() in ("1", "true", "yes")


def e1_analysis(base: Path) -> Dict[str, object]:
    metrics = read_csv(base / "e1_search" / "e1_metrics.csv")
    targets = read_csv(base / "e1_outage" / "outage_targets.csv")
    metric_by_id = {row["candidate_id"]: row for row in metrics}
    baseline = next(
        row for row in metrics
        if json.loads(row["delay_indices"]) == list(range(8))
    )
    baseline_target = next(
        row for row in targets if row["candidate_id"] == baseline["candidate_id"]
    )
    rows: List[Dict[str, object]] = []
    for target in targets:
        metric = metric_by_id[target["candidate_id"]]
        gain = (
            float(baseline_target["outage10_snr_db"])
            - float(target["outage10_snr_db"])
        )
        ce16_degradation = (
            float(metric["ce_nmse_db_16db"])
            - float(baseline["ce_nmse_db_16db"])
        )
        floor_degradation = (
            float(metric["ce_nmse_db_80db"])
            - float(baseline["ce_nmse_db_80db"])
        )
        non_arithmetic = json.loads(metric["delay_indices"]) != list(range(8))
        passed = bool(
            non_arithmetic
            and gain >= 0.15
            and ce16_degradation <= 1.0
            and floor_degradation <= 1.0
            and int(float(metric["ce_pilot_rank_16db"])) == 8
            and math.isfinite(float(metric["ce_condition_number_16db"]))
        )
        rows.append({
            "candidate_id": target["candidate_id"],
            "family": target["family"],
            "delay_indices": target["delay_indices"],
            "delay_ns": target["delay_ns"],
            "outage10_snr_db": float(target["outage10_snr_db"]),
            "outage10_gain_vs_B1_db": gain,
            "ce_nmse_db_16db": float(metric["ce_nmse_db_16db"]),
            "ce16_degradation_vs_B1_db": ce16_degradation,
            "ce_nmse_db_80db": float(metric["ce_nmse_db_80db"]),
            "ce80_floor_degradation_vs_B1_db": floor_degradation,
            "pilot_rank": int(float(metric["ce_pilot_rank_16db"])),
            "h1_pass": passed,
        })
    write_csv(rows, base / "final" / "e1_h1_candidates.csv")
    passing = [row for row in rows if row["h1_pass"]]
    return {
        "baseline_candidate_id": baseline["candidate_id"],
        "baseline_delay_indices": baseline["delay_indices"],
        "tested_outage_candidates": len(rows),
        "passing_candidates": passing,
        "h1_pass": bool(passing),
        "best_outage_candidate": min(rows, key=lambda row: float(row["outage10_snr_db"])),
        "closest_joint_candidate": min(
            (row for row in rows if row["candidate_id"] != baseline["candidate_id"]),
            key=lambda row: (
                max(0.0, 0.15 - float(row["outage10_gain_vs_B1_db"]))
                + max(0.0, float(row["ce16_degradation_vs_B1_db"]) - 1.0)
                + max(0.0, float(row["ce80_floor_degradation_vs_B1_db"]) - 1.0)
            ),
        ),
    }


def e2_analysis(base: Path) -> Dict[str, object]:
    metrics = read_csv(base / "e2_thick_sidon" / "e2_metrics.csv")
    targets = read_csv(base / "e2_outage" / "outage_targets.csv")
    metric_by_id = {row["candidate_id"]: row for row in metrics}
    scenario_rows: List[Dict[str, object]] = []
    passing_scenarios: List[str] = []
    for scenario_id in sorted({row["scenario_id"] for row in targets}):
        selected = [row for row in targets if row["scenario_id"] == scenario_id]
        baseline_target = next(row for row in selected if row["family"] == "B1_T1")
        baseline_metric = metric_by_id[baseline_target["candidate_id"]]
        candidate_rows: List[Dict[str, object]] = []
        for target in selected:
            metric = metric_by_id[target["candidate_id"]]
            gain = (
                float(baseline_target["outage10_snr_db"])
                - float(target["outage10_snr_db"])
            )
            degradation = (
                float(metric["ce_nmse_db_16db"])
                - float(baseline_metric["ce_nmse_db_16db"])
            )
            hard = as_bool(metric["hard_thick_sidon"])
            passed = bool(
                hard
                and gain >= 0.10
                and degradation <= 1.0
                and int(float(metric["ce_pilot_rank_16db"])) == 8
            )
            candidate_rows.append({
                "scenario_id": scenario_id,
                "candidate_id": target["candidate_id"],
                "family": target["family"],
                "delay_indices": target["delay_indices"],
                "delay_ns": target["delay_ns"],
                "hard_thick_sidon": hard,
                "outage10_snr_db": float(target["outage10_snr_db"]),
                "outage10_gain_vs_B1_T1_db": gain,
                "ce_nmse_db_16db": float(metric["ce_nmse_db_16db"]),
                "ce16_degradation_vs_B1_T1_db": degradation,
                "h2_scenario_pass": passed,
            })
        scenario_rows.extend(candidate_rows)
        if any(row["h2_scenario_pass"] for row in candidate_rows):
            passing_scenarios.append(scenario_id)
    write_csv(scenario_rows, base / "final" / "e2_h2_candidates.csv")
    best = {}
    for scenario_id in sorted({row["scenario_id"] for row in scenario_rows}):
        selected = [
            row for row in scenario_rows
            if row["scenario_id"] == scenario_id and row["hard_thick_sidon"]
        ]
        best[scenario_id] = (
            max(selected, key=lambda row: float(row["outage10_gain_vs_B1_T1_db"]))
            if selected else None
        )
    return {
        "passing_scenarios": passing_scenarios,
        "passing_scenario_count": len(passing_scenarios),
        "h2_pass": len(passing_scenarios) >= 2,
        "best_hard_candidate_by_scenario": best,
    }


def e3_analysis(base: Path) -> Dict[str, object]:
    summary_path = base / "final" / "e3_summary.json"
    raw = json.loads(summary_path.read_text(encoding="utf-8"))
    rows: List[Dict[str, object]] = []
    for scenario_id, scenario in raw["scenarios"].items():
        for label in ("bler10", "bler1"):
            value = scenario[label]
            if not value.get("available"):
                rows.append({
                    "scenario_id": scenario_id,
                    "target": label,
                    "available": False,
                })
                continue
            qc = value["candidates"]["QC"]
            sidon = value["candidates"]["Sidon"]
            rows.append({
                "scenario_id": scenario_id,
                "target": label,
                "available": True,
                "qc_snr_db": qc["snr_db"],
                "sidon_snr_db": sidon["snr_db"],
                "sidon_gain_db": value["sidon_gain_vs_qc_db"],
                "gain_ci95_lo_db": value["sidon_gain_ci95_conservative_lo_db"],
                "gain_ci95_hi_db": value["sidon_gain_ci95_conservative_hi_db"],
                "qc_target_band_errors": qc["errors_in_half_to_twice_target_band"],
                "sidon_target_band_errors": sidon["errors_in_half_to_twice_target_band"],
                "qc_deterministic": qc["sufficient_for_deterministic_claim"],
                "sidon_deterministic": sidon["sufficient_for_deterministic_claim"],
            })
    write_csv(rows, base / "final" / "e3_targets.csv")
    a01 = next(
        row for row in rows
        if row["scenario_id"] == "A_0p1ns" and row["target"] == "bler10"
    )
    return {
        "targets": rows,
        "h3a_pass": bool(
            a01["available"]
            and float(a01["sidon_gain_db"]) >= 0.15
            and float(a01["gain_ci95_lo_db"]) >= 0.0
        ),
        "h3b_interval_supported": [
            row["scenario_id"] for row in rows
            if row["target"] == "bler10"
            and row["scenario_id"][0] in ("D", "E")
            and row["available"]
            and float(row["gain_ci95_lo_db"]) > 0.0
        ],
        "h3b_unresolved_target": ["D_1ns"],
    }


def candidate_catalog(base: Path) -> None:
    rows: List[Dict[str, object]] = []
    for source, phase in (
        (base / "e1_search" / "e1_metrics.csv", "E1"),
        (base / "e2_thick_sidon" / "e2_metrics.csv", "E2"),
    ):
        for row in read_csv(source):
            rows.append({
                "phase": phase,
                "scenario_id": row.get("scenario_id", "E1-A5"),
                "candidate_id": row["candidate_id"],
                "family": row["family"],
                "delay_indices": row["delay_indices"],
                "delay_ns": row["delay_ns"],
                "residues": row["residues"],
                "lifts": row["lifts"],
                "pair_sum_min_gap_ns": row["pair_sum_min_gap_ns"],
                "fold_min_gap_ns": row["fold_min_gap_ns"],
                "knowledge_level": "T2/R2" if phase == "E1" else "T1/R2",
            })
    write_csv(rows, base / "final" / "candidate_catalog.csv")


def _plot_setup(output: Path):
    os.environ.setdefault("MPLCONFIGDIR", str(output / "mplcache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_e1_outage(base: Path, path: Path) -> None:
    plt = _plot_setup(base / "final")
    rows = read_csv(base / "e1_outage" / "outage_curves.csv")
    figure, axis = plt.subplots(figsize=(8.2, 5.2))
    for candidate_id in sorted({row["candidate_id"] for row in rows}):
        selected = sorted(
            (row for row in rows if row["candidate_id"] == candidate_id),
            key=lambda row: float(row["snr_db"]),
        )
        family = selected[0]["family"]
        highlight = family in ("S0", "B0", "B1+E1_O")
        axis.semilogy(
            [float(row["snr_db"]) for row in selected],
            [max(float(row["outage"]), 2.5e-6) for row in selected],
            lw=2.1 if highlight else 0.9,
            alpha=1.0 if highlight else 0.55,
            label=f"{candidate_id} ({family})" if highlight else None,
        )
    axis.axhline(0.1, color="#555555", ls="--", lw=1)
    axis.axhline(0.01, color="#777777", ls=":", lw=1)
    axis.set(xlabel="average receive SNR (dB)", ylabel="ideal-CSI outage")
    axis.set_title("Plan-026 E1: frozen candidates on TDL-A 5 ns")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def plot_e2_gain(base: Path, path: Path) -> None:
    plt = _plot_setup(base / "final")
    rows = read_csv(base / "final" / "e2_h2_candidates.csv")
    scenarios = sorted({row["scenario_id"] for row in rows})
    values = []
    labels = []
    for scenario in scenarios:
        hard = [
            row for row in rows
            if row["scenario_id"] == scenario and as_bool(row["hard_thick_sidon"])
        ]
        if hard:
            best = max(hard, key=lambda row: float(row["outage10_gain_vs_B1_T1_db"]))
            values.append(float(best["outage10_gain_vs_B1_T1_db"]))
            labels.append(best["candidate_id"])
        else:
            values.append(float("nan"))
            labels.append("no tested hard candidate")
    figure, axis = plt.subplots(figsize=(8.2, 4.8))
    bars = axis.bar(scenarios, values, color="#2a9d8f")
    axis.axhline(0.10, color="#d62828", ls="--", label="H2 threshold")
    for bar, label in zip(bars, labels):
        if np.isfinite(bar.get_height()):
            axis.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.006,
                label,
                ha="center",
                va="bottom",
                fontsize=7,
                rotation=25,
            )
    axis.set_ylabel("10% outage gain vs T1 arithmetic baseline (dB)")
    axis.set_title("Plan-026 E2: best tested hard thick-Sidon candidate")
    axis.grid(True, axis="y", alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def _scenario_rows(base: Path, directory: str, scenario: str) -> List[Dict[str, str]]:
    return read_csv(base / directory / scenario / "sidon_qc_bler.csv")


def plot_e3_bler(base: Path, path: Path) -> None:
    plt = _plot_setup(base / "final")
    scenarios = (
        "A_0p1ns", "A_1ns", "A_5ns", "D_0p1ns",
        "D_1ns", "E_0p1ns", "E_1ns",
    )
    figure, axes = plt.subplots(2, 4, figsize=(15.2, 8.0))
    colors = {"QC": "#e76f51", "Sidon": "#264653"}
    for axis, scenario in zip(axes.flat, scenarios):
        rows = _scenario_rows(base, "e3_refine_10pct", scenario)
        if not rows:
            rows = _scenario_rows(base, "e3_prescan", scenario)
        for design in ("QC", "Sidon"):
            selected = sorted(
                (row for row in rows if row["id"] == design),
                key=lambda row: float(row["snr_db"]),
            )
            axis.semilogy(
                [float(row["snr_db"]) for row in selected],
                [max(float(row["bler"]), 1e-4) for row in selected],
                marker="o",
                color=colors[design],
                label=design,
            )
        axis.axhline(0.1, color="#555555", ls="--", lw=1)
        axis.set_title(scenario)
        axis.grid(True, which="both", alpha=0.25)
        axis.set_xlabel("SNR (dB)")
        axis.set_ylabel("estimated-CSI BLER")
    axes.flat[-1].axis("off")
    axes.flat[0].legend()
    figure.suptitle("Plan-026 E3: 10% BLER refinement (D1 uses prescan)")
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def plot_e3_gain(base: Path, path: Path) -> None:
    plt = _plot_setup(base / "final")
    rows = [
        row for row in read_csv(base / "final" / "e3_targets.csv")
        if row["target"] == "bler10" and as_bool(row["available"])
    ]
    labels = [row["scenario_id"] for row in rows]
    gains = np.asarray([float(row["sidon_gain_db"]) for row in rows])
    lower = gains - np.asarray([float(row["gain_ci95_lo_db"]) for row in rows])
    upper = np.asarray([float(row["gain_ci95_hi_db"]) for row in rows]) - gains
    figure, axis = plt.subplots(figsize=(9.2, 4.8))
    axis.bar(labels, gains, color=["#457b9d" if label.startswith("A") else "#2a9d8f" for label in labels])
    axis.errorbar(
        np.arange(len(labels)), gains, yerr=np.vstack([lower, upper]),
        fmt="none", color="#111111", capsize=4,
    )
    axis.axhline(0.0, color="#333333", lw=1)
    axis.set_ylabel("Sidon 10% BLER SNR gain vs QC (dB)")
    axis.set_title("Plan-026 E3 target-SNR gains (conservative 95% intervals)")
    axis.tick_params(axis="x", rotation=30)
    axis.grid(True, axis="y", alpha=0.25)
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def plot_e3_one_percent(base: Path, path: Path) -> None:
    plt = _plot_setup(base / "final")
    scenarios = ("A_0p1ns", "A_1ns", "A_5ns")
    figure, axes = plt.subplots(1, 3, figsize=(13.6, 4.3))
    colors = {"QC": "#e76f51", "Sidon": "#264653"}
    for axis, scenario in zip(axes, scenarios):
        rows = _scenario_rows(base, "e3_refine_1pct", scenario)
        for design in ("QC", "Sidon"):
            selected = sorted(
                (row for row in rows if row["id"] == design),
                key=lambda row: float(row["snr_db"]),
            )
            axis.semilogy(
                [float(row["snr_db"]) for row in selected],
                [max(float(row["bler"]), 1e-4) for row in selected],
                marker="o",
                color=colors[design],
                label=design,
            )
        axis.axhline(0.01, color="#555555", ls="--", lw=1)
        axis.set_title(scenario)
        axis.set_xlabel("SNR (dB)")
        axis.set_ylabel("estimated-CSI BLER")
        axis.grid(True, which="both", alpha=0.25)
    axes[0].legend()
    figure.suptitle("Plan-026 E3: conditional 1% BLER refinement")
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", default=DEFAULT_RUN_ID)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    args = parser.parse_args()
    base = args.output_root / str(args.run_id)
    final = base / "final"
    final.mkdir(parents=True, exist_ok=True)
    candidate_catalog(base)
    e1 = e1_analysis(base)
    e2 = e2_analysis(base)
    e3 = e3_analysis(base)
    plot_e1_outage(base, final / "e1_outage_curves.png")
    plot_e2_gain(base, final / "e2_h2_gain_map.png")
    plot_e3_bler(base, final / "e3_bler_10pct.png")
    plot_e3_gain(base, final / "e3_gain_summary.png")
    plot_e3_one_percent(base, final / "e3_bler_1pct_A.png")
    summary = {
        "experiment": "plan-026",
        "run_id": str(args.run_id),
        "e1": e1,
        "e2": e2,
        "e3": e3,
        "verdicts": {
            "H1": bool(e1["h1_pass"]),
            "H2": bool(e2["h2_pass"]),
            "H3a": bool(e3["h3a_pass"]),
            "H3b_interval_supported_scenarios": e3["h3b_interval_supported"],
            "H3b_unresolved_target_scenarios": e3["h3b_unresolved_target"],
        },
    }
    (final / "final_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary["verdicts"], indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
