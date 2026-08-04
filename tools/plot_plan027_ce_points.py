"""Plot all evaluated comb-24 CE NMSE points used by result-027.

Every candidate present in the E3 comb-24 results is included. Exact
five-baseline comparison points added for Section 5 are included where
available. Colors, line styles, and markers are shared with E4.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np

from tools.run_plan027_bler import curve_style


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RUN_ROOT = (
    ROOT / "outputs" / "experiment027_meff_sidon" / "20260726_main"
)
DEFAULT_FIGURE = ROOT / "docs" / "figures" / "result-027" / "ce_nmse_points.png"
SCENARIOS = ("A1", "A5", "A10", "A30", "A100")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def add_point(
    points: dict[tuple[str, str, float], dict[str, object]],
    *,
    scenario_id: str,
    candidate_id: str,
    family: str,
    snr_db: float,
    nmse_db: float,
    source: str,
) -> None:
    key = (scenario_id, candidate_id, round(float(snr_db), 12))
    if key in points:
        old = points[key]
        if not math.isclose(
            float(old["ce_nmse_db"]), float(nmse_db), abs_tol=2e-9, rel_tol=0.0
        ):
            raise RuntimeError(
                f"Conflicting CE values for {scenario_id}/{candidate_id} at "
                f"{snr_db:.12g} dB: {old['ce_nmse_db']} vs {nmse_db}"
            )
        old["source"] = ";".join(
            sorted(set(str(old["source"]).split(";")) | {source})
        )
        return
    points[key] = {
        "scenario_id": scenario_id,
        "candidate_id": candidate_id,
        "family": family,
        "snr_db": float(snr_db),
        "ce_nmse_db": float(nmse_db),
        "source": source,
    }


def collect_points(run_root: Path) -> list[dict[str, object]]:
    e3_rows = read_csv(run_root / "e3_ce_density" / "e3_ce_at_targets.csv")
    five_rows = read_csv(
        run_root / "final" / "five_baseline_candidate_comparison.csv"
    )
    selected = {
        scenario: {
            row["candidate_id"]
            for row in e3_rows
            if row["scenario_id"] == scenario
            and int(row["dmrs_spacing_subcarriers"]) == 24
        }
        for scenario in SCENARIOS
    }
    family_by_id = {
        row["candidate_id"]: row["family"] for row in e3_rows
    }
    points: dict[tuple[str, str, float], dict[str, object]] = {}

    for row in e3_rows:
        scenario_id = row["scenario_id"]
        candidate_id = row["candidate_id"]
        if int(row["dmrs_spacing_subcarriers"]) != 24:
            continue
        add_point(
            points,
            scenario_id=scenario_id,
            candidate_id=candidate_id,
            family=row["family"],
            snr_db=float(row["snr_db"]),
            nmse_db=float(row["ce_nmse_db"]),
            source=f"e3:{row['working_point']}",
        )

    for row in five_rows:
        scenario_id = row["scenario_id"]
        candidate_id = row["representative_candidate_id"]
        if candidate_id in selected[scenario_id]:
            add_point(
                points,
                scenario_id=scenario_id,
                candidate_id=candidate_id,
                family=family_by_id[candidate_id],
                snr_db=float(row["baseline_target_snr_db"]),
                nmse_db=float(row["candidate_ce_at_baseline_target_nmse_db"]),
                source=f"five_baseline:{row['baseline_family']}",
            )
        baseline_id = row["baseline_candidate_id"]
        if baseline_id in selected[scenario_id]:
            add_point(
                points,
                scenario_id=scenario_id,
                candidate_id=baseline_id,
                family=row["baseline_family"],
                snr_db=float(row["baseline_target_snr_db"]),
                nmse_db=float(row["baseline_ce_at_own_target_nmse_db"]),
                source=f"five_baseline:{row['baseline_family']}",
            )

    rows = sorted(
        points.values(),
        key=lambda row: (
            SCENARIOS.index(str(row["scenario_id"])),
            str(row["candidate_id"]),
            float(row["snr_db"]),
        ),
    )
    missing = {
        scenario: sorted(
            candidate
            for candidate in selected[scenario]
            if not any(
                row["scenario_id"] == scenario
                and row["candidate_id"] == candidate
                for row in rows
            )
        )
        for scenario in SCENARIOS
    }
    if any(missing.values()):
        raise RuntimeError(f"E3 comb-24 candidates without CE points: {missing}")
    return rows


def write_points(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_points(rows: list[dict[str, object]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(5, 1, figsize=(18.0, 20.0))
    for axis, scenario_id in zip(axes, SCENARIOS):
        scenario_rows = [
            row for row in rows if row["scenario_id"] == scenario_id
        ]
        candidate_ids = sorted(
            {str(row["candidate_id"]) for row in scenario_rows}
        )
        for candidate_id in candidate_ids:
            selected = sorted(
                (
                    row
                    for row in scenario_rows
                    if row["candidate_id"] == candidate_id
                ),
                key=lambda row: float(row["snr_db"]),
            )
            family = str(selected[0]["family"])
            axis.plot(
                np.asarray([float(row["snr_db"]) for row in selected]),
                np.asarray([float(row["ce_nmse_db"]) for row in selected]),
                linewidth=1.3,
                markersize=4,
                label=candidate_id,
                **curve_style(candidate_id, family),
            )
        axis.set_xlabel("SNR (dB)")
        axis.set_ylabel("Matched CE NMSE (dB)")
        axis.set_title(
            f"{scenario_id}: evaluated comb-24 CE points "
            f"({len(candidate_ids)} candidates)"
        )
        axis.grid(True, alpha=0.25)
        axis.legend(
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            fontsize=7,
            frameon=False,
        )
    figure.tight_layout(rect=(0.0, 0.0, 0.78, 1.0))
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument(
        "--points-csv",
        type=Path,
        default=DEFAULT_RUN_ROOT / "final" / "ce_nmse_plotted_points.csv",
    )
    args = parser.parse_args()
    rows = collect_points(args.run_root)
    write_points(rows, args.points_csv)
    plot_points(rows, args.figure)
    print(
        f"wrote {args.figure.relative_to(ROOT)} and "
        f"{args.points_csv.relative_to(ROOT)} ({len(rows)} points)"
    )


if __name__ == "__main__":
    main()
