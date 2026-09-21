"""Analyze plan-033 prescan data or available formal scenarios.

Formal plots use saved Monte Carlo points directly. No smoothing, monotonic
correction, or extrapolation is applied.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


SCENARIOS = (
    "A100_NT8_NR4_V3",
    "A100_NT8_NR4_V60",
    "A100_NT4_NR4_V3",
    "A100_NT4_NR4_V60",
)
TARGETS = (0.1, 0.01)
EXPECTED_POINTS = {
    "A100_NT8_NR4_V3": 33,
    "A100_NT8_NR4_V60": 17,
    "A100_NT4_NR4_V3": 29,
    "A100_NT4_NR4_V60": 17,
}
LABELS = {
    "AGED_MRT_PRG6": "Aged MRT PRG6 (5 ms)",
    "B0_QC": "B0-QC CDD",
    "S0_SIDON": "Sidon CDD",
    "SMALL_CDD_QSTEP0P25_MATCHED": "Small CDD, matched Rx",
    "SMALL_CDD_QSTEP0P25_TRANSPARENT": "Small CDD, transparent Rx",
    "TRANSPARENT_PRG6": "Transparent PRG6",
}
STYLE_ORDER = tuple(LABELS)
STYLE_DEFS = {
    "AGED_MRT_PRG6": {"color": "#E69F00", "linestyle": "--", "marker": "D"},
    "B0_QC": {"color": "#0072B2", "linestyle": "-", "marker": "s"},
    "S0_SIDON": {"color": "#009E73", "linestyle": "-", "marker": "^"},
    "SMALL_CDD_QSTEP0P25_MATCHED": {"color": "#D55E00", "linestyle": "-.", "marker": "v"},
    "SMALL_CDD_QSTEP0P25_TRANSPARENT": {"color": "#CC79A7", "linestyle": ":", "marker": "X"},
    "TRANSPARENT_PRG6": {"color": "#000000", "linestyle": "-", "marker": "o"},
}


def _read_points(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _candidate_key(candidate_id: str) -> str:
    for key in STYLE_ORDER:
        if candidate_id.endswith(key):
            return key
    raise ValueError(f"Unknown plan-033 candidate: {candidate_id}")


def _bracket(points: list[dict], target: float) -> dict:
    ordered = sorted(points, key=lambda row: float(row["snr_db"]))
    for left, right in zip(ordered, ordered[1:]):
        left_bler, right_bler = float(left["bler"]), float(right["bler"])
        if left_bler >= target >= right_bler:
            crossing = None
            if left_bler > 0.0 and right_bler > 0.0 and left_bler != right_bler:
                fraction = (
                    math.log10(target) - math.log10(left_bler)
                ) / (math.log10(right_bler) - math.log10(left_bler))
                crossing = float(left["snr_db"]) + fraction * (
                    float(right["snr_db"]) - float(left["snr_db"])
                )
            return {
                "status": "bracketed" if crossing is not None else "bracketed_zero_endpoint",
                "crossing_snr_db": crossing,
                "lower_snr_db": float(left["snr_db"]),
                "lower_bler": left_bler,
                "lower_errors": int(left["tb_errors"]),
                "lower_trials": int(left["trials"]),
                "lower_wilson95_lo": float(left["bler_wilson95_lo"]),
                "lower_wilson95_hi": float(left["bler_wilson95_hi"]),
                "upper_snr_db": float(right["snr_db"]),
                "upper_bler": right_bler,
                "upper_errors": int(right["tb_errors"]),
                "upper_trials": int(right["trials"]),
                "upper_wilson95_lo": float(right["bler_wilson95_lo"]),
                "upper_wilson95_hi": float(right["bler_wilson95_hi"]),
            }
    return {"status": "unbracketed", "crossing_snr_db": None}


def analyze_prescan(root: Path) -> dict:
    report = {"schema": "plan033-prescan-bracket-audit-v1", "root": str(root), "scenarios": {}}
    for scenario in SCENARIOS:
        path = root / scenario / "final" / "estimated_csi_bler_points.csv"
        rows = _read_points(path)
        by_candidate: dict[str, list[dict]] = {}
        for row in rows:
            if row["scenario_id"] != scenario:
                raise RuntimeError(f"Scenario mismatch in {path}")
            by_candidate.setdefault(row["candidate_id"], []).append(row)
        if len(by_candidate) != 6:
            raise RuntimeError(f"{scenario} has {len(by_candidate)} candidates; expected 6.")
        report["scenarios"][scenario] = {
            candidate: {f"bler_{target:g}": _bracket(points, target) for target in TARGETS}
            for candidate, points in sorted(by_candidate.items())
        }
    return report


def _available_formal(root: Path) -> dict[str, list[dict]]:
    available = {}
    for scenario in SCENARIOS:
        path = root / scenario / "final" / "estimated_csi_bler_points.csv"
        if not path.exists():
            continue
        rows = _read_points(path)
        candidates = {row["candidate_id"] for row in rows}
        snrs = {float(row["snr_db"]) for row in rows}
        if len(candidates) != 6 or len(snrs) != EXPECTED_POINTS[scenario]:
            continue
        if any(int(row["trials"]) < 10_000 for row in rows):
            continue
        available[scenario] = rows
    return available


def _plot_scenario(scenario: str, rows: list[dict], figure_dir: Path) -> list[dict]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_dir.mkdir(parents=True, exist_ok=True)
    scenario_slug = scenario.lower()
    outputs = []

    def draw(kind: str, preview: bool) -> Path:
        size = (5.12, 4.1) if preview else (9.5, 6.6)
        fig, ax = plt.subplots(figsize=size)
        for key in STYLE_ORDER:
            points = sorted(
                (row for row in rows if _candidate_key(row["candidate_id"]) == key),
                key=lambda row: float(row["snr_db"]),
            )
            x = [float(row["snr_db"]) for row in points]
            if kind == "bler":
                y = [float(row["bler"]) if float(row["bler"]) >= 0.005 else math.nan for row in points]
            else:
                y = [float(row["ce_nmse_mean_db"]) for row in points]
            style = STYLE_DEFS[key]
            ax.plot(
                x,
                y,
                label=LABELS[key],
                color=style["color"],
                linestyle=style["linestyle"],
                marker=style["marker"],
                linewidth=2.5,
                markersize=8,
                markeredgewidth=1.0,
            )
        if kind == "bler":
            ax.set_yscale("log")
            ax.set_ylim(0.005, 1.0)
            ax.axhline(0.1, color="#666666", linewidth=1.2, linestyle="--")
            ax.axhline(0.01, color="#666666", linewidth=1.2, linestyle=":")
            ax.set_ylabel("Transport-block error rate", fontsize=16)
            title_suffix = "estimated-CSI BLER"
        else:
            ax.set_ylabel("Data-RE CE NMSE (dB)", fontsize=16)
            title_suffix = "channel-estimation NMSE"
        ax.set_xlabel("SNR (dB), unit total transmit power", fontsize=16)
        ax.set_title(f"{scenario}: {title_suffix}", fontsize=16)
        ax.grid(True, which="both", alpha=0.28)
        ax.tick_params(axis="both", labelsize=14)
        ax.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.20),
            ncol=2,
            fontsize=10 if preview else 13,
            frameon=False,
        )
        fig.tight_layout()
        suffix = "_preview_13cm" if preview else ""
        path = figure_dir / f"{scenario_slug}_{kind}{suffix}.png"
        fig.savefig(path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        return path

    for kind in ("bler", "ce_nmse"):
        for preview in (False, True):
            path = draw(kind, preview)
            outputs.append(
                {
                    "scenario_id": scenario,
                    "kind": kind,
                    "preview_13cm": preview,
                    "path": str(path),
                    "bytes": path.stat().st_size,
                }
            )
    return outputs


def analyze_formal(root: Path, analysis_dir: Path, figure_dir: Path) -> dict:
    available = _available_formal(root)
    if not available:
        raise RuntimeError("No complete 10,000-trial formal scenario is available.")
    combined = []
    targets = []
    append_rows = []
    figures = []
    for scenario, rows in available.items():
        combined.extend(rows)
        by_candidate: dict[str, list[dict]] = {}
        for row in rows:
            by_candidate.setdefault(row["candidate_id"], []).append(row)
        for candidate_id, points in sorted(by_candidate.items()):
            for target in TARGETS:
                bracket = _bracket(points, target)
                targets.append(
                    {
                        "scenario_id": scenario,
                        "candidate_id": candidate_id,
                        "label": LABELS[_candidate_key(candidate_id)],
                        "target_bler": target,
                        **bracket,
                    }
                )
                if target == 0.01 and bracket["status"].startswith("bracketed"):
                    for side in ("lower", "upper"):
                        bler = float(bracket[f"{side}_bler"])
                        errors = int(bracket[f"{side}_errors"])
                        if 0.005 <= bler <= 0.02 and errors < 200:
                            append_rows.append(
                                {
                                    "scenario_id": scenario,
                                    "candidate_id": candidate_id,
                                    "snr_db": bracket[f"{side}_snr_db"],
                                    "bler": bler,
                                    "errors": errors,
                                    "trials": bracket[f"{side}_trials"],
                                    "reason": "1pct_bracket_endpoint_in_0p5_to_2pct_with_lt_200_errors",
                                }
                            )
        figures.extend(_plot_scenario(scenario, rows, figure_dir))

    _write_csv(analysis_dir / "available_formal_points.csv", combined)
    _write_csv(analysis_dir / "target_summary_initial.csv", targets)
    if append_rows:
        _write_csv(analysis_dir / "append_requirements_initial.csv", append_rows)
    style_path = analysis_dir / "curve_styles.json"
    style_path.parent.mkdir(parents=True, exist_ok=True)
    style_path.write_text(
        json.dumps(
            {key: {"label": LABELS[key], **STYLE_DEFS[key]} for key in STYLE_ORDER},
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    report = {
        "schema": "plan033-partial-formal-analysis-v1",
        "status": "provisional_initial_10000_trials",
        "formal_root": str(root),
        "available_scenarios": list(available),
        "missing_scenarios": [scenario for scenario in SCENARIOS if scenario not in available],
        "target_rows": len(targets),
        "append_requirement_rows": len(append_rows),
        "figures": figures,
        "plot_policy": "raw points only; BLER display clipped below 0.005; no smoothing, monotonic correction, or extrapolation",
    }
    report_path = analysis_dir / "partial_analysis.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prescan-root", type=Path)
    parser.add_argument("--formal-root", type=Path)
    parser.add_argument("--analysis-dir", type=Path)
    parser.add_argument("--figure-dir", type=Path)
    args = parser.parse_args()
    if bool(args.prescan_root) == bool(args.formal_root):
        parser.error("Specify exactly one of --prescan-root or --formal-root.")
    if args.prescan_root:
        root = args.prescan_root.resolve()
        report = analyze_prescan(root)
        output = root / "prescan_bracket_audit.json"
        output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    else:
        if args.analysis_dir is None or args.figure_dir is None:
            parser.error("Formal analysis requires --analysis-dir and --figure-dir.")
        output = args.analysis_dir.resolve() / "partial_analysis.json"
        report = analyze_formal(
            args.formal_root.resolve(),
            args.analysis_dir.resolve(),
            args.figure_dir.resolve(),
        )
    print(json.dumps(report, indent=2))
    print(f"[written] {output}")


if __name__ == "__main__":
    main()
