"""Analyze Plan-037 section 13 Beam8 prescan/formal link outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from statistics import NormalDist
from typing import Any

import numpy as np
import yaml


SCHEMES = ("beam8_b0_qc", "beam8_s0_sidon", "beam8_precoder_cycling")
RECEIVERS = ("ideal", "estimated")
CURVES = tuple(f"{scheme}_{receiver}" for scheme in SCHEMES for receiver in RECEIVERS)
TARGETS = (0.10, 0.01)
SCENARIO = "E100_NT32_NR2_V60_BEAM8"


def calibration_plots(gate_a_dir: Path, manifest_path: Path, output: Path) -> None:
    """Plot the frozen CDL angular spectrum and DFT2x8 beam powers."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    data = np.load(gate_a_dir / "gate_a_numeric.npz")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    aod = np.rad2deg(np.asarray(data["aod_rad"]).reshape(-1))
    zod = np.rad2deg(np.asarray(data["zod_rad"]).reshape(-1))
    cluster_power = np.asarray(data["powers"]).reshape(-1)
    rays_per_cluster = int(np.asarray(data["aod_rad"]).shape[-1])
    ray_power = np.repeat(cluster_power / rays_per_cluster, rays_per_cluster)
    k_factor = float(np.asarray(data["k_factor"]).reshape(-1)[0])
    los = bool(np.asarray(data["los_indicator"]).reshape(-1)[0])
    if los:
        ray_power = ray_power / (k_factor + 1.0)
        aod = np.r_[aod, np.rad2deg(np.asarray(data["los_aod_rad"]).reshape(-1)[0])]
        zod = np.r_[zod, np.rad2deg(np.asarray(data["los_zod_rad"]).reshape(-1)[0])]
        ray_power = np.r_[ray_power, k_factor / (k_factor + 1.0)]
    aod = (aod + 180.0) % 360.0 - 180.0
    bins_aod = np.linspace(-180.0, 180.0, 145)
    bins_zod = np.linspace(0.0, 180.0, 73)
    spectrum, _, _ = np.histogram2d(aod, zod, bins=(bins_aod, bins_zod), weights=ray_power)
    positive = spectrum[spectrum > 0.0]
    floor = float(np.max(positive) * 1e-5) if positive.size else 1e-12
    shown = np.where(spectrum.T > 0.0, spectrum.T, np.nan)

    output.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(11, 6.5))
    mesh = ax.pcolormesh(
        bins_aod, bins_zod, shown, shading="auto",
        norm=LogNorm(vmin=floor, vmax=float(np.nanmax(shown))), cmap="magma",
    )
    sizes = 25.0 + 130.0 * np.sqrt(ray_power / np.max(ray_power))
    ax.scatter(aod, zod, s=sizes, facecolors="none", edgecolors="cyan", linewidths=0.7, alpha=0.8)
    ax.set_xlabel("AoD (deg)", fontsize=16); ax.set_ylabel("ZoD (deg)", fontsize=16)
    ax.set_title("CDL-E transmit angular power spectrum", fontsize=17)
    ax.tick_params(labelsize=14); ax.grid(True, alpha=0.2)
    colorbar = fig.colorbar(mesh, ax=ax); colorbar.set_label("Ray power per angular bin", fontsize=14)
    fig.tight_layout(); fig.savefig(output / "cdl_e_aod_zod_power.png", dpi=180); plt.close(fig)

    analytic = 10.0 * np.log10(np.asarray(data["analytic_beam_powers"]) / np.max(data["analytic_beam_powers"]))
    mc = 10.0 * np.log10(np.asarray(manifest["mean_rsrp"]) / np.max(manifest["mean_rsrp"]))
    selected = set(int(value) for value in manifest["selected_beam_indices"])
    fig, axes = plt.subplots(2, 1, figsize=(12, 5.8), sharex=True)
    for ax, values, label in zip(axes, (analytic, mc), ("Analytic", f"Monte Carlo, D={manifest['D']}")):
        matrix = values.reshape(2, 8)
        image = ax.imshow(matrix, aspect="auto", vmin=min(-20.0, float(np.min(values))), vmax=0.0, cmap="viridis")
        for qv in range(2):
            for qh in range(8):
                beam = qv * 8 + qh
                marker = "*" if beam in selected else ""
                ax.text(qh, qv, f"b{beam}{marker}\n{matrix[qv,qh]:.2f} dB", ha="center", va="center",
                        color="white" if matrix[qv,qh] < -4 else "black", fontsize=11, fontweight="bold")
        ax.set_yticks([0, 1], ["qv=0", "qv=1"], fontsize=13); ax.set_title(label, fontsize=16)
    axes[-1].set_xticks(range(8), [f"qh={value}" for value in range(8)], fontsize=12)
    colorbar = fig.colorbar(image, ax=axes, pad=0.02); colorbar.set_label("Normalized beam power (dB)", fontsize=14)
    fig.suptitle(f"DFT2x8 beam powers; * selected ({manifest['status']})", fontsize=17)
    fig.subplots_adjust(left=.09, right=.88, bottom=.12, top=.86, hspace=.38)
    fig.savefig(output / "beam_power_heatmap.png", dpi=180); plt.close(fig)

    rows = [{"beam_id": beam, "q_v": beam // 8, "q_h": beam % 8,
             "analytic_normalized_db": analytic[beam], "mc_normalized_db": mc[beam],
             "selected": int(beam in selected), "selection_status": manifest["status"]}
            for beam in range(16)]
    write_csv(output / "beam_power_heatmap_data.csv", rows)
    report = {
        "stage": "calibration-plots", "manifest_status": manifest["status"], "D": manifest["D"],
        "selected_beam_indices": sorted(selected), "los_included": los,
        "k_factor_linear": k_factor, "angular_samples": int(len(ray_power)),
        "outputs": ["cdl_e_aod_zod_power.png", "beam_power_heatmap.png", "beam_power_heatmap_data.csv"],
    }
    (output / "calibration_plot_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


def pdp_plots(report_path: Path, output: Path) -> None:
    """Render Plan-037 section 13.5 figures from the saved numeric report."""
    import matplotlib.pyplot as plt

    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("status") != "PASS":
        raise ValueError("PDP plotting requires a PASS pdp_report.json.")
    output.mkdir(parents=True, exist_ok=True)

    def stem(ax, table: dict[str, Any], label: str, **kwargs: Any) -> None:
        ax.plot(np.asarray(table["delays_s"]) * 1e9, table["powers"], marker="o", label=label, **kwargs)

    fig, ax = plt.subplots(figsize=(10, 5.5))
    stem(ax, report["standard_final"], "Standard CDL-E final PDP")
    ax.set(xlabel="Delay from common physical origin (ns)", ylabel="Normalized power", title="Standard CDL-E PDP")
    ax.grid(True, alpha=.3); ax.legend(); fig.tight_layout(); fig.savefig(output / "standard_cdl_e_pdp.png", dpi=180); plt.close(fig)

    reference = int(report["reference_beam_id"])
    fig, ax = plt.subplots(figsize=(11, 6))
    strongest = max(float(value) for value in report["beam_total_power"].values())
    for beam_text, table in report["beam_pdp"].items():
        beam = int(beam_text); total = float(report["beam_total_power"][beam_text])
        power = np.asarray(table["powers"]) * total / strongest
        ax.plot(np.asarray(table["delays_s"]) * 1e9, power, marker="o", linewidth=2.5 if beam == reference else 1.2,
                label=f"b{beam}" + (" (reference)" if beam == reference else ""))
    ax.set(xlabel="Delay from common physical origin (ns)", ylabel="Power / strongest beam total",
           title="Selected Beam8 absolute-shape PDPs")
    ax.grid(True, alpha=.3); ax.legend(ncol=2); fig.tight_layout(); fig.savefig(output / "beam8_absolute_pdp.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 6))
    for beam_text, table in report["beam_pdp"].items():
        beam = int(beam_text); stem(ax, table, f"b{beam}" + (" (reference)" if beam == reference else ""),
                                    linewidth=2.5 if beam == reference else 1.2)
    ax.set(xlabel="Delay from common physical origin (ns)", ylabel="Per-beam normalized power",
           title="Selected Beam8 normalized PDPs")
    ax.grid(True, alpha=.3); ax.legend(ncol=2); fig.tight_layout(); fig.savefig(output / "beam8_normalized_pdp.png", dpi=180); plt.close(fig)

    rows = []
    with (report_path.parent / "pdp_tables.csv").open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    fig, ax = plt.subplots(figsize=(11, 6))
    for kind, label in (("reference", "Reference"), ("cdd_effective", "CDD")):
        groups: dict[str, list[dict[str, str]]] = {}
        for row in rows:
            if row["pdp_type"] == kind:
                groups.setdefault(row["beam_id"], []).append(row)
        for name, group in groups.items():
            ax.plot([float(row["delay_ns"]) for row in group], [float(row["power"]) for row in group], marker="o",
                    label=label if kind == "reference" else name)
    ax.set(xlabel="Delay modulo 1/SCS (ns)", ylabel="Normalized assumed power",
           title="Reference and CDD-aware assumed effective PDPs")
    ax.grid(True, alpha=.3); ax.legend(); fig.tight_layout(); fig.savefig(output / "reference_and_cdd_effective_pdp.png", dpi=180); plt.close(fig)
    (output / "pdp_plot_report.json").write_text(json.dumps({"status": "PASS", "source": str(report_path),
        "outputs": ["standard_cdl_e_pdp.png", "beam8_absolute_pdp.png", "beam8_normalized_pdp.png",
                    "reference_and_cdd_effective_pdp.png"]}, indent=2) + "\n", encoding="utf-8")


def read_rows(paths: list[Path], name: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for directory in paths:
        with (directory / name).open("r", encoding="utf-8", newline="") as handle:
            rows.extend(csv.DictReader(handle))
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else [])
        if rows:
            writer.writeheader(); writer.writerows(rows)


def aggregate_trials(rows: list[dict[str, str]]) -> dict[tuple[str, float], dict[int, int]]:
    curves: dict[tuple[str, float], dict[int, int]] = {}
    identities: dict[tuple[float, int], tuple[str, str, str]] = {}
    for row in rows:
        if row.get("scenario_id") != SCENARIO:
            raise ValueError(f"Unexpected scenario_id={row.get('scenario_id')!r}; expected {SCENARIO}.")
        scheme = row["variant_id"]
        snr = float(row["snr_db"])
        trial = int(row["trial"])
        key = (scheme, snr)
        if trial in curves.setdefault(key, {}):
            raise ValueError(f"Duplicate absolute trial: {key}, trial={trial}.")
        curves[key][trial] = int(row["tb_error"])
        identity = (row["channel_realization_index"], row["payload_noise_seed"], row["noise_var"])
        common_key = (snr, trial)
        if common_key in identities and identities[common_key] != identity:
            raise ValueError(f"Paired identity mismatch at snr={snr}, trial={trial}.")
        identities[common_key] = identity
    for snr in sorted({key[1] for key in curves}):
        sets = [set(curves[(scheme, snr)]) for scheme in CURVES]
        if len({tuple(sorted(value)) for value in sets}) != 1:
            raise ValueError(f"Schemes do not share identical absolute trials at {snr} dB.")
        ordered = sorted(sets[0])
        if not ordered or ordered[0] != 1 or ordered != list(range(1, ordered[-1] + 1)):
            raise ValueError(f"Absolute trials are not contiguous from trial 1 at {snr} dB.")
    return curves


def crossing(points: list[tuple[float, float]], target: float) -> float:
    points = sorted(points)
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if y0 >= target >= y1 and y0 > 0 and y1 > 0 and y0 != y1:
            return x0 + (math.log(target) - math.log(y0)) * (x1 - x0) / (math.log(y1) - math.log(y0))
    return float("nan")


def crossing_bracket(points: list[tuple[float, float]], target: float) -> tuple[float, float] | None:
    """Return a genuine adjacent-point bracket, allowing a zero-error lower endpoint."""
    points = sorted(points)
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if y0 >= target >= y1:
            return float(x0), float(x1)
    return None


def wilson(errors: int, trials: int) -> tuple[float, float]:
    z = NormalDist().inv_cdf(0.975); p = errors / trials
    d = 1 + z*z/trials
    c = (p + z*z/(2*trials))/d
    h = z*math.sqrt(p*(1-p)/trials + z*z/(4*trials*trials))/d
    return c-h, c+h


def point_rows(curves: dict[tuple[str, float], dict[int, int]]) -> list[dict[str, Any]]:
    rows = []
    for (scheme, snr), values in sorted(curves.items()):
        errors = sum(values.values()); trials = len(values); low, high = wilson(errors, trials)
        rows.append({"scheme": scheme, "snr_db": snr, "trials": trials, "errors": errors,
                     "bler": errors/trials, "wilson_95_low": low, "wilson_95_high": high,
                     "absolute_trial_min": min(values), "absolute_trial_max": max(values)})
    return rows


def make_formal_grid(curves: dict[tuple[str, float], dict[int, int]]) -> list[float]:
    bracket_edges = []
    for scheme in CURVES:
        points = [(snr, np.mean(list(values.values()))) for (name, snr), values in curves.items() if name == scheme]
        for target in TARGETS:
            bracket = crossing_bracket(points, target)
            if bracket is None:
                direction = "lower" if all(y < target for _, y in points) else "higher"
                if (direction == "higher" and max(x for x, _ in points) >= 22.0) or (
                    direction == "lower" and min(x for x, _ in points) <= -10.0
                ):
                    raise ValueError(
                        f"Prescan lacks a two-sided {target:g} bracket for {scheme} at the frozen "
                        f"[-10,22] dB hard boundary; stop without extrapolation."
                    )
                raise ValueError(f"Prescan lacks a two-sided {target:g} bracket for {scheme}; extend {direction} SNR.")
            bracket_edges.extend(bracket)
    low = math.floor((min(bracket_edges)-0.5)*4)/4
    high = math.ceil((max(bracket_edges)+0.5)*4)/4
    if low < -10 or high > 22:
        raise ValueError("Required formal grid exceeds the frozen [-10,22] dB boundary.")
    return np.round(np.arange(low, high + 0.001, 0.25), 2).tolist()


def bootstrap(curves: dict[tuple[str, float], dict[int, int]], repeats: int, seed: int) -> list[dict[str, Any]]:
    rng = np.random.default_rng(seed)
    snrs = sorted({key[1] for key in curves})
    observed = {scheme: {snr: np.mean(list(curves[(scheme, snr)].values())) for snr in snrs} for scheme in CURVES}
    estimates = {(scheme, target): crossing(list(observed[scheme].items()), target) for scheme in CURVES for target in TARGETS}
    samples = {(scheme, target): [] for scheme in CURVES for target in TARGETS}
    gains = {(scheme, receiver, target): [] for scheme in SCHEMES[:2] for receiver in RECEIVERS for target in TARGETS}
    for _ in range(repeats):
        replicate = {scheme: {} for scheme in CURVES}
        for snr in snrs:
            ids = sorted(curves[(CURVES[0], snr)])
            draw = rng.integers(0, len(ids), len(ids))
            for scheme in CURVES:
                values = np.asarray([curves[(scheme, snr)][trial] for trial in ids])
                replicate[scheme][snr] = float(np.mean(values[draw]))
        for target in TARGETS:
            values = {scheme: crossing(list(replicate[scheme].items()), target) for scheme in CURVES}
            if all(math.isfinite(value) for value in values.values()):
                for scheme in CURVES:
                    samples[(scheme, target)].append(values[scheme])
                for receiver in RECEIVERS:
                    cycling = f"{SCHEMES[2]}_{receiver}"
                    for scheme in SCHEMES[:2]:
                        key = f"{scheme}_{receiver}"
                        gains[(scheme, receiver, target)].append(values[cycling] - values[key])
    rows = []
    for target in TARGETS:
        for scheme in CURVES:
            values = np.asarray(samples[(scheme, target)])
            rows.append({"metric": "crossing", "scheme": scheme, "target_bler": target,
                         "estimate_db": estimates[(scheme,target)], "ci95_low_db": np.quantile(values, .025),
                         "ci95_high_db": np.quantile(values, .975), "bootstrap_valid": len(values)})
        for receiver in RECEIVERS:
          for scheme in SCHEMES[:2]:
            values = np.asarray(gains[(scheme, receiver, target)])
            estimate = estimates[(f"{SCHEMES[2]}_{receiver}",target)] - estimates[(f"{scheme}_{receiver}",target)]
            rows.append({"metric": "gain_vs_cycling", "scheme": scheme, "receiver": receiver, "target_bler": target,
                         "estimate_db": estimate, "ci95_low_db": np.quantile(values,.025),
                         "ci95_high_db": np.quantile(values,.975), "bootstrap_valid": len(values)})
    return rows


def plot(points: list[dict[str, Any]], output: Path) -> None:
    import matplotlib.pyplot as plt
    colors = {SCHEMES[0]: "#0072B2", SCHEMES[1]: "#D55E00", SCHEMES[2]: "#009E73"}
    fig, ax = plt.subplots(figsize=(9, 6))
    for scheme in CURVES:
        rows = [row for row in points if row["scheme"] == scheme]
        base, receiver = scheme.rsplit("_", 1); color = colors[base]
        marker, line = (("o", "-") if receiver == "ideal" else ("s", "--"))
        ax.semilogy([r["snr_db"] for r in rows], [max(r["bler"], .5/r["trials"]) for r in rows],
                    label=scheme, color=color, marker=marker, linestyle=line, linewidth=2.5, markersize=8)
    ax.set_xlabel("Reference SNR (dB)", fontsize=16); ax.set_ylabel("TB BLER", fontsize=16)
    ax.grid(True, which="both", alpha=.3); ax.tick_params(labelsize=14); ax.legend(fontsize=12)
    fig.tight_layout(); output.parent.mkdir(parents=True, exist_ok=True); fig.savefig(output, dpi=180); plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("calibration", "pdp", "prescan", "formal"), required=True)
    parser.add_argument("--run-dir", type=Path, action="append")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gate-a-dir", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--pdp-report", type=Path)
    parser.add_argument("--base-config", type=Path)
    parser.add_argument("--formal-config", type=Path)
    parser.add_argument("--bootstrap-repeats", type=int, default=2000)
    args = parser.parse_args()
    if args.stage == "calibration":
        if args.gate_a_dir is None or args.manifest is None:
            parser.error("calibration requires --gate-a-dir and --manifest")
        calibration_plots(args.gate_a_dir, args.manifest, args.output)
        return
    if args.stage == "pdp":
        if args.pdp_report is None:
            parser.error("pdp requires --pdp-report")
        pdp_plots(args.pdp_report, args.output)
        return
    if not args.run_dir:
        parser.error("prescan/formal requires at least one --run-dir")
    curves = aggregate_trials(read_rows(args.run_dir, "trial_metrics.csv"))
    points = point_rows(curves); args.output.mkdir(parents=True, exist_ok=True)
    write_csv(args.output / "bler_points.csv", points)
    plot(points, args.output / "beam8_bler.png")
    if args.stage == "prescan":
        try:
            grid = make_formal_grid(curves)
        except ValueError as error:
            report = {"stage": "prescan", "status": "STOPPED", "reason": str(error),
                      "formal_config_generated": False}
            (args.output / "analysis.json").write_text(
                json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
            )
            print(f"[prescan] status=STOPPED reason={error}", flush=True)
            return
        data = yaml.safe_load(args.base_config.read_text(encoding="utf-8"))
        data["simulation"].update({"snr_points_db": grid, "n_trials_per_snr": 1000,
                                   "max_trials_per_snr": 1000, "absolute_trial_start": 0,
                                   "run_id": "formal_00000_01000",
                                   "output_dir": "outputs/experiment037_cdl_e_32tx_beam8/formal"})
        args.formal_config.write_text(yaml.safe_dump(data, sort_keys=False, allow_unicode=True), encoding="utf-8")
        report = {"stage": "prescan", "status": "READY_FOR_FORMAL_CONFIRMATION",
                  "formal_grid_db": grid, "formal_config": str(args.formal_config)}
    else:
        metrics = bootstrap(curves, args.bootstrap_repeats, 20263704)
        write_csv(args.output / "crossings_and_gains.csv", metrics)
        report = {"stage": "formal", "bootstrap_repeats": args.bootstrap_repeats, "metrics": metrics}
    (args.output / "analysis.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
