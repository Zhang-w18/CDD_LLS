"""Plot ideal-outage and estimated-CSI link target coordinates for E4.

For every E4-evaluated physical curve in A5, A10, and A30, each target plot
contains:

* the ideal-outage target SNR and matched CE NMSE at that ideal work point;
* the estimated-CSI BLER target SNR and interpolated CE NMSE at that link
  target, when the formal E4 scan closed.

Candidate colors, line styles, and markers are shared with the E4 BLER plots.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from tools.run_plan027_bler import curve_style


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RUN_ROOT = (
    ROOT / "outputs" / "experiment027_meff_sidon" / "20260726_main"
)
DEFAULT_FIGURE_ROOT = ROOT / "docs" / "figures" / "result-027"
SCENARIOS = ("A5", "A10", "A30")
TARGETS = (
    ("10pct", 0.1, "outage10_snr_db", "10%"),
    ("1pct", 0.01, "outage1_snr_db", "1%"),
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def collect_points(run_root: Path) -> list[dict[str, object]]:
    manifest_path = run_root / "e4_link_gate" / "e4_link_manifest.json"
    with manifest_path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)

    manifest_candidates = [
        candidate
        for candidate in manifest["candidates"]
        if candidate["scenario_id"] in SCENARIOS
    ]
    e3_rows = read_csv(run_root / "e3_ce_density" / "e3_ce_at_targets.csv")
    e3_own = {
        (
            row["scenario_id"],
            row["candidate_id"],
            round(float(row["target_outage_probability"]), 8),
        ): row
        for row in e3_rows
        if int(row["dmrs_spacing_subcarriers"]) == 24
        and row["working_point"] == "candidate_own_target"
    }
    link_rows = []
    for scenario_id in SCENARIOS:
        link_rows.extend(
            read_csv(
                run_root
                / "e4_link"
                / scenario_id
                / "final"
                / "target_summary.csv"
            )
        )
    link_by_key = {
        (row["scenario_id"], row["candidate_id"], row["target"]): row
        for row in link_rows
    }

    points: list[dict[str, object]] = []
    for candidate in manifest_candidates:
        scenario_id = candidate["scenario_id"]
        candidate_id = candidate["candidate_id"]
        for target, probability, outage_field, _ in TARGETS:
            e3_key = (scenario_id, candidate_id, round(probability, 8))
            if e3_key not in e3_own:
                raise RuntimeError(f"Missing E3 own-workpoint CE row: {e3_key}")
            e3_row = e3_own[e3_key]
            ideal_snr_db = float(candidate[outage_field])
            if not math.isclose(
                ideal_snr_db,
                float(e3_row["snr_db"]),
                rel_tol=0.0,
                abs_tol=2e-9,
            ):
                raise RuntimeError(
                    f"Manifest/E3 outage SNR mismatch for {e3_key}: "
                    f"{ideal_snr_db} vs {e3_row['snr_db']}"
                )
            link_row = link_by_key.get((scenario_id, candidate_id, target))
            point = {
                "scenario_id": scenario_id,
                "target": target,
                "target_probability": probability,
                "candidate_id": candidate_id,
                "family": candidate["family"],
                "ideal_outage_snr_db": ideal_snr_db,
                "ideal_workpoint_ce_nmse_db": float(e3_row["ce_nmse_db"]),
                "link_target_closed": link_row is not None,
                "link_target_snr_db": "",
                "link_target_ce_nmse_db": "",
                "link_target_ci95_low_db": "",
                "link_target_ci95_high_db": "",
            }
            if link_row is not None:
                if not math.isclose(
                    ideal_snr_db,
                    float(link_row["outage_target_snr_db"]),
                    rel_tol=0.0,
                    abs_tol=2e-9,
                ):
                    raise RuntimeError(
                        f"Manifest/E4 outage SNR mismatch for "
                        f"{scenario_id}/{candidate_id}/{target}"
                    )
                point.update(
                    {
                        "link_target_snr_db": float(
                            link_row["target_snr_db"]
                        ),
                        "link_target_ce_nmse_db": float(
                            link_row["target_ce_nmse_db"]
                        ),
                        "link_target_ci95_low_db": float(
                            link_row["ci95_low_db"]
                        ),
                        "link_target_ci95_high_db": float(
                            link_row["ci95_high_db"]
                        ),
                    }
                )
            points.append(point)
    return points


def write_points(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_target(
    rows: list[dict[str, object]],
    *,
    scenario_id: str,
    target: str,
    target_label: str,
    path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MultipleLocator

    selected = [
        row
        for row in rows
        if row["scenario_id"] == scenario_id and row["target"] == target
    ]
    closed_count = sum(bool(row["link_target_closed"]) for row in selected)
    figure, axis = plt.subplots(figsize=(15.5, 9.0))
    candidate_handles = []

    for row in selected:
        candidate_id = str(row["candidate_id"])
        family = str(row["family"])
        style = curve_style(candidate_id, family)
        color = style["color"]
        marker = style["marker"]
        linestyle = style["linestyle"]
        ideal_x = float(row["ideal_outage_snr_db"])
        ideal_y = float(row["ideal_workpoint_ce_nmse_db"])

        if row["link_target_closed"]:
            link_x = float(row["link_target_snr_db"])
            link_y = float(row["link_target_ce_nmse_db"])
            low = float(row["link_target_ci95_low_db"])
            high = float(row["link_target_ci95_high_db"])
            axis.plot(
                [ideal_x, link_x],
                [ideal_y, link_y],
                color=color,
                linestyle=linestyle,
                linewidth=1.25,
                alpha=0.65,
                zorder=2,
            )
            axis.errorbar(
                link_x,
                link_y,
                xerr=[[link_x - low], [high - link_x]],
                fmt=marker,
                color=color,
                markerfacecolor=color,
                markeredgecolor=color,
                markersize=7,
                elinewidth=1.0,
                capsize=3,
                alpha=0.95,
                zorder=4,
            )

        axis.plot(
            ideal_x,
            ideal_y,
            linestyle="None",
            marker=marker,
            color=color,
            markerfacecolor="none",
            markeredgecolor=color,
            markeredgewidth=1.5,
            markersize=8,
            zorder=5,
        )
        candidate_handles.append(
            Line2D(
                [0],
                [0],
                color=color,
                linestyle=linestyle,
                linewidth=1.3,
                marker=marker,
                markerfacecolor=color,
                markeredgecolor=color,
                markersize=6,
                label=candidate_id.removeprefix(f"{scenario_id}_"),
            )
        )

    x_values = [float(row["ideal_outage_snr_db"]) for row in selected]
    y_values = [float(row["ideal_workpoint_ce_nmse_db"]) for row in selected]
    for row in selected:
        if row["link_target_closed"]:
            x_values.extend(
                [
                    float(row["link_target_ci95_low_db"]),
                    float(row["link_target_ci95_high_db"]),
                ]
            )
            y_values.append(float(row["link_target_ce_nmse_db"]))
    axis.set_xlim(math.floor(min(x_values)) - 0.25, math.ceil(max(x_values)) + 0.25)
    axis.set_ylim(math.floor(min(y_values)) - 0.5, math.ceil(max(y_values)) + 0.5)

    axis.xaxis.set_major_locator(MultipleLocator(1.0))
    axis.xaxis.set_minor_locator(MultipleLocator(0.25))
    axis.yaxis.set_major_locator(MultipleLocator(2.0))
    axis.yaxis.set_minor_locator(MultipleLocator(0.5))
    axis.grid(which="major", color="#7f7f7f", linewidth=0.65, alpha=0.45)
    axis.grid(
        which="minor",
        color="#a6a6a6",
        linewidth=0.35,
        linestyle=":",
        alpha=0.35,
    )
    axis.set_axisbelow(True)
    axis.set_xlabel("Target SNR (dB)")
    axis.set_ylabel("Matched CE NMSE at target SNR (dB)")
    axis.set_title(
        f"{scenario_id} — {target_label} target: ideal outage and "
        f"estimated-CSI link ({closed_count}/{len(selected)} E4 targets closed)"
    )

    semantic_handles = [
        Line2D(
            [0],
            [0],
            color="#222222",
            linestyle="None",
            marker="o",
            markerfacecolor="none",
            markeredgewidth=1.5,
            markersize=8,
            label="Ideal outage + CE at ideal work point",
        ),
        Line2D(
            [0],
            [0],
            color="#222222",
            linestyle="None",
            marker="o",
            markerfacecolor="#222222",
            markersize=7,
            label="Estimated-CSI BLER target + CE at link target",
        ),
        Line2D(
            [0],
            [0],
            color="#555555",
            linestyle="-",
            linewidth=1.1,
            marker="|",
            markersize=8,
            label="Horizontal bar: link target SNR 95% CI",
        ),
    ]
    semantic_legend = axis.legend(
        handles=semantic_handles,
        loc="lower left",
        fontsize=8,
        frameon=True,
        framealpha=0.88,
    )
    axis.add_artist(semantic_legend)
    axis.legend(
        handles=candidate_handles,
        title="E4 candidates",
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        fontsize=8,
        title_fontsize=9,
        frameon=False,
    )
    figure.tight_layout(rect=(0.0, 0.0, 0.80, 1.0))
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument(
        "--figure-root", type=Path, default=DEFAULT_FIGURE_ROOT
    )
    parser.add_argument(
        "--points-csv",
        type=Path,
        default=(
            DEFAULT_RUN_ROOT
            / "final"
            / "outage_ce_link_scatter_points.csv"
        ),
    )
    args = parser.parse_args()

    rows = collect_points(args.run_root)
    write_points(rows, args.points_csv)
    figure_paths = []
    for scenario_id in SCENARIOS:
        for target, _, _, target_label in TARGETS:
            path = (
                args.figure_root
                / f"{scenario_id.lower()}_{target}_outage_ce_link.png"
            )
            plot_target(
                rows,
                scenario_id=scenario_id,
                target=target,
                target_label=target_label,
                path=path,
            )
            figure_paths.append(path)

    print(
        f"wrote {len(figure_paths)} figures and "
        f"{args.points_csv.relative_to(ROOT)} ({len(rows)} candidate-target rows)"
    )


if __name__ == "__main__":
    main()
