"""Execute Plan-039 validation, selection, and frozen formal supplements."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import config_from_dict, dataclass_to_dict, load_config
from cdd_lls.phy.cdl_beam_platform import (
    beam_domain_cdd_precoder, beam_path_statistics, covariance_audit,
    direction_cosines, rows_to_csv, txru_response_matrix,
)
from cdd_lls.phy.channel_cdl_fixed import FixedCDLStatisticsChannel
from cdd_lls.phy.estimators import cdl_spatial_unaware_covariance
from cdd_lls.phy.plan039 import (
    build_angular_full_coverage_codebook, frequency_covariance,
    load_manifest as load_plan039_manifest, search_strict_sidon_top8, write_manifest,
)
from cdd_lls.phy.resource_grid import build_resource_grid
from cdd_lls.sim.orchestrator import CDDLinkLevelOrchestrator

PLAN = "research/plan-039-CDL宽波束CDD与cycling-BLER.md"
SCENARIO = "C_UE2R_DUAL_ASD25"
CODEBOOK_TYPE = "angular_full_coverage_ultrawide"
SCHEMES = ("BEAM8_B0_QC", "BEAM8_SIDON_SELECTED", "BEAM8_PRECODER_CYCLING")
CDD_METHODS = ("PLAN039_COMMON_REFERENCE_PDP", "PLAN039_BEAM_SPECIFIC_PDP_INDEPENDENT")
STAGE1B_CDD_METHODS = (
    "PLAN039_COMMON_REFERENCE_PDP",
)
STAGE1B_CYCLING_METHOD = "PLAN039_PRG_COMMON_REFERENCE_PDP"
TRANSPARENT_CDD_METHOD = "PLAN039_TRANSPARENT_COMMON_REFERENCE_PDP"
AGED_MRT_SCHEME = "PLAN039_AGED_MRT_PRG6"
AGED_MRT_METHOD = "PLAN039_PRG_COMMON_REFERENCE_PDP"
SELECTION_GRID = tuple(float(value) for value in range(0, 22, 2))
SELECTION_TRIALS = 400
FORMAL_BATCH = 1000
FORMAL_CAP = 50000
SNR_MIN = -6.0
SNR_MAX = 26.0


def _json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _stable_seed(base: int, namespace: str) -> int:
    digest = hashlib.sha256(f"{int(base)}|plan039|{namespace}".encode("utf-8")).digest()
    return int.from_bytes(digest[:4], "big")


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _crossing(rows: list[dict[str, Any]], target: float) -> tuple[float, tuple[float, float] | None]:
    points = sorted((float(row["snr_db"]), float(row["bler"])) for row in rows
                    if math.isfinite(float(row["bler"])))
    brackets: list[tuple[float, float, float, float]] = []
    for (s0, b0), (s1, b1) in zip(points[:-1], points[1:]):
        if b0 > 0.0 and b1 > 0.0 and b0 >= target >= b1 and b0 != b1:
            brackets.append((s0, b0, s1, b1))
    if not brackets:
        return float("nan"), None
    s0, b0, s1, b1 = brackets[0]
    fraction = (math.log10(target) - math.log10(b0)) / (math.log10(b1) - math.log10(b0))
    return float(s0 + fraction * (s1 - s0)), (s0, s1)


def _wilson(errors: int, trials: int) -> tuple[float, float]:
    if trials <= 0:
        return float("nan"), float("nan")
    z = 1.959963984540054
    p = float(errors) / float(trials)
    denominator = 1.0 + z*z/trials
    center = (p + z*z/(2.0*trials)) / denominator
    half = z * math.sqrt(p*(1.0-p)/trials + z*z/(4.0*trials*trials)) / denominator
    return max(0.0, center-half), min(1.0, center+half)


def _merge_summary(batch_dirs: list[Path]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, float], list[dict[str, str]]] = {}
    for directory in batch_dirs:
        for row in _read_csv(directory / "summary.csv"):
            groups.setdefault((row["variant_id"], float(row["snr_db"])), []).append(row)
    merged: list[dict[str, Any]] = []
    for (variant, snr), rows in sorted(groups.items()):
        first: dict[str, Any] = dict(rows[0])
        trials = sum(int(row["n_trials"]) for row in rows)
        errors = sum(int(row["tb_errors"]) for row in rows)
        cb_errors = sum(int(row["cb_errors"]) for row in rows)
        ce_sum = sum(float(row["ce_nmse_sum"]) for row in rows)
        ce_sq = sum(float(row["ce_nmse_sum_squares"]) for row in rows)
        ce_only = first["receiver_mode"] == "ce_only"
        wilson_low, wilson_high = (float("nan"), float("nan")) if ce_only else _wilson(errors, trials)
        first.update({"variant_id": variant, "snr_db": snr, "n_trials": trials,
                      "absolute_trial_start": 0, "absolute_trial_stop": trials,
                      "tb_errors": errors, "cb_errors": cb_errors,
                      "bler": float("nan") if ce_only else errors / trials,
                      "wilson95_low": wilson_low, "wilson95_high": wilson_high,
                      "cb_bler": float("nan") if ce_only else cb_errors / max(trials * int(first["n_cbs"]), 1),
                      "ce_nmse_eff": ce_sum / trials, "ce_nmse_sum": ce_sum,
                      "ce_nmse_sum_squares": ce_sq})
        merged.append(first)
    return merged


def _collect_trials(batch_dirs: list[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for directory in batch_dirs:
        rows.extend(_read_csv(directory / "trial_metrics.csv"))
    return rows


def _merge_trial_files(batch_dirs: list[Path], output: Path, variants: set[str] | None = None) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] | None = None
    with output.open("w", encoding="utf-8", newline="") as target:
        for directory in batch_dirs:
            with (directory / "trial_metrics.csv").open("r", encoding="utf-8", newline="") as source:
                reader = csv.DictReader(source)
                current = list(reader.fieldnames or [])
                if fieldnames is None:
                    fieldnames = current
                    writer = csv.DictWriter(target, fieldnames=fieldnames, lineterminator="\n")
                    writer.writeheader()
                elif current != fieldnames:
                    raise RuntimeError(f"Trial CSV schema changed in {directory}.")
                writer = csv.DictWriter(target, fieldnames=fieldnames, lineterminator="\n")
                for row in reader:
                    if variants is None or str(row["variant_id"]) in variants:
                        writer.writerow(row)


def _run_batch(config: dict[str, Any], config_path: Path, expected_dir: Path) -> Path:
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(yaml.safe_dump(config, sort_keys=False, allow_unicode=True), encoding="utf-8")
    config_sha256 = _sha256(config_path)
    complete = expected_dir / "summary.csv"
    if complete.exists() and (expected_dir / "trial_metrics.csv").exists():
        receipt_path = expected_dir / "batch_receipt.json"
        if not receipt_path.exists():
            raise RuntimeError(f"Existing batch has no resume receipt: {expected_dir}.")
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if receipt.get("config_sha256") != config_sha256:
            raise RuntimeError(f"Existing batch configuration changed: {expected_dir}.")
        summary_count = len(_read_csv(complete))
        trial_count = len(_read_csv(expected_dir / "trial_metrics.csv"))
        expected_summary = len(config["variants"]) * len(config["simulation"]["snr_points_db"])
        expected_trials = expected_summary * int(config["simulation"]["n_trials_per_snr"])
        if summary_count == expected_summary and trial_count == expected_trials:
            return expected_dir
        raise RuntimeError(
            f"Incomplete existing batch {expected_dir}: summary {summary_count}/{expected_summary}, "
            f"trials {trial_count}/{expected_trials}. Move it aside before resuming."
        )
    directory = CDDLinkLevelOrchestrator(load_config(config_path)).run()
    _json(directory / "batch_receipt.json", {"schema": "plan039-batch-receipt-v1",
          "config": config_path.resolve().as_posix(), "config_sha256": config_sha256,
          "summary_rows": len(_read_csv(directory / "summary.csv")),
          "trial_rows": len(_read_csv(directory / "trial_metrics.csv"))})
    return directory


def _scenario_payload(raw: dict[str, Any], target_asd_deg: float, selection_only: bool = True) -> dict[str, Any]:
    payload = copy.deepcopy(raw)
    for key in ("plan039", "beam_design", "visualization"):
        payload.pop(key, None)
    payload["channel"]["cdl_profile"] = "C"
    payload["antenna"].update({"n_tx": 32, "n_rx": 2})
    payload["fixed_cdl_statistics"].update({
        "codebook_type": CODEBOOK_TYPE, "profile_native_angles": False,
        "mean_aod_deg": 0.0, "target_aod_asd_deg": float(target_asd_deg),
        "aod_scale": 1.0, "aoa_scale": 1.0, "zod_scale": 1.0, "zoa_scale": 1.0,
        "ue_vertical_elements": 1, "ue_horizontal_elements": 1,
        "ue_polarizations": 2, "ue_polarization_type": "cross",
        "selection_only": bool(selection_only),
    })
    payload["scenarios"] = [{"scenario_id": SCENARIO}]
    payload["variants"] = []
    return payload


def _validate_raw(raw: dict[str, Any]) -> None:
    resource = raw["resource"]
    actual = (int(resource["n_prbs"]), int(resource["n_fft"]), int(resource["pdsch_n_symbols"]),
              list(resource["dmrs_symbol_indices"]), int(resource["dmrs_spacing_sc"]),
              int(resource["prg_size_rb"]))
    if actual != (48, 4096, 10, [2, 7], 6, 6):
        raise ValueError(f"Plan-039 resource mismatch: {actual}.")
    if str(raw["channel"]["cdl_profile"]).upper() != "C":
        raise ValueError("Plan-039 contains only CDL-C.")
    if float(raw["channel"]["delay_spread_ns"]) != 100.0 or float(raw["channel"]["ue_speed_kmh"]) != 3.0:
        raise ValueError("Plan-039 freezes 100 ns and 3 km/h.")
    plan = raw["plan039"]
    if float(plan["target_aod_asd_deg"]) != 25.0 or float(plan["diagnostic_aod_asd_deg"]) != 10.0:
        raise ValueError("Plan-039 freezes ASD25 with ASD10 as diagnostic only.")
    if list(plan["smoke"]["snr_points_db"]) != [4.0, 12.0] or int(plan["smoke"]["trials_per_snr"]) != 20:
        raise ValueError("Plan-039 stage-0 smoke is two SNRs and 20 paired trials per point.")


def _ray_table(channel: FixedCDLStatisticsChannel) -> dict[str, np.ndarray]:
    arrays = channel.profile_arrays()
    aod = np.asarray(arrays["aod_rad"], dtype=np.float64).reshape(-1)
    zod = np.asarray(arrays["zod_rad"], dtype=np.float64).reshape(-1)
    clusters = np.asarray(arrays["powers"], dtype=np.float64).reshape(-1)
    delays = np.asarray(arrays["delays_s"], dtype=np.float64).reshape(-1)
    rays_per_cluster = aod.size // clusters.size
    los = bool(np.asarray(arrays["los_indicator"]).reshape(-1)[0])
    k_factor = float(np.asarray(arrays["k_factor"]).reshape(-1)[0])
    diffuse_scale = 1.0 / (k_factor + 1.0) if los else 1.0
    power = np.repeat(clusters * diffuse_scale / rays_per_cluster, rays_per_cluster)
    delay = np.repeat(delays, rays_per_cluster)
    kind = np.asarray(["diffuse"] * power.size, dtype="U16")
    cluster = np.repeat(np.arange(clusters.size), rays_per_cluster)
    ray = np.tile(np.arange(rays_per_cluster), clusters.size)
    if los:
        aod = np.r_[aod, float(np.asarray(arrays["los_aod_rad"]).reshape(-1)[0])]
        zod = np.r_[zod, float(np.asarray(arrays["los_zod_rad"]).reshape(-1)[0])]
        delay = np.r_[delay, delays[0]]
        power = np.r_[power, k_factor / (k_factor + 1.0)]
        kind = np.r_[kind, "los_specular"]
        cluster, ray = np.r_[cluster, 0], np.r_[ray, -1]
    power /= np.sum(power)
    return {"aod_deg": np.rad2deg(aod), "zod_deg": np.rad2deg(zod), "delay_s": delay,
            "power": power, "kind": kind, "cluster": cluster, "ray": ray}


def _circular_mean_asd(angle_deg: np.ndarray, power: np.ndarray) -> tuple[float, float]:
    mean = float(np.rad2deg(np.angle(np.sum(power * np.exp(1j * np.deg2rad(angle_deg))))))
    offset = (angle_deg - mean + 180.0) % 360.0 - 180.0
    return mean, float(np.sqrt(np.sum(power * offset**2) / np.sum(power)))


def _angle_audit(rays: dict[str, np.ndarray], target_asd: float) -> dict[str, Any]:
    mean, asd = _circular_mean_asd(rays["aod_deg"], rays["power"])
    aod_mask = (rays["aod_deg"] >= -60.0) & (rays["aod_deg"] <= 60.0)
    rectangle = aod_mask & (rays["zod_deg"] >= 90.0) & (rays["zod_deg"] <= 110.0)
    result = {"power_weighted_circular_mean_aod_deg": mean, "power_weighted_rms_asd_deg": asd,
              "target_asd_deg": float(target_asd), "aod_interval_power": float(np.sum(rays["power"][aod_mask])),
              "target_rectangle_power": float(np.sum(rays["power"][rectangle]))}
    result["gate"] = {"mean": abs(mean) <= 0.1, "asd": abs(asd-target_asd) <= 0.1,
                      "aod_coverage": result["aod_interval_power"] >= 0.93,
                      "rectangle_coverage": result["target_rectangle_power"] >= 0.92}
    result["status"] = "PASS" if all(result["gate"].values()) else "FAIL"
    return result


def _normalized_correlation(covariance: np.ndarray) -> tuple[np.ndarray, dict[str, float]]:
    dimension = covariance.shape[0]
    diagonal = np.maximum(np.real(np.diag(covariance)), np.finfo(float).tiny)
    correlation = np.abs(covariance) / np.sqrt(diagonal[:, None] * diagonal[None, :])
    off = correlation[~np.eye(dimension, dtype=bool)]
    eig = np.maximum(np.linalg.eigvalsh(0.5 * (covariance + covariance.conj().T)), 0.0)
    probability = eig / max(float(np.sum(eig)), np.finfo(float).tiny)
    effective_rank = float(np.exp(-np.sum(probability[probability > 0] * np.log(probability[probability > 0]))))
    return correlation, {"off_diagonal_mean": float(np.mean(off)), "off_diagonal_median": float(np.median(off)),
                         "off_diagonal_p95": float(np.percentile(off, 95)), "off_diagonal_max": float(np.max(off)),
                         "effective_rank": effective_rank}


def _correlation(channel: FixedCDLStatisticsChannel, weights: np.ndarray, count: int):
    beam_accum = np.zeros((8, 8), dtype=np.complex128)
    txru_accum = np.zeros((32, 32), dtype=np.complex128)
    sample_count = 0
    for index in range(int(count)):
        h = np.asarray(channel.generate(index).H[0], dtype=np.complex128)
        beam = np.einsum("rtnk,tb->rnkb", h, weights, optimize=True)
        beam_accum += np.einsum("rnkb,rnkc->bc", beam, beam.conj(), optimize=True)
        txru_accum += np.einsum("rtnk,runk->tu", h, h.conj(), optimize=True)
        sample_count += int(beam.shape[0] * beam.shape[1] * beam.shape[2])
    beam_accum /= float(sample_count)
    txru_accum /= float(sample_count)
    beam_correlation, beam_summary = _normalized_correlation(beam_accum)
    txru_correlation, txru_summary = _normalized_correlation(txru_accum)
    return beam_correlation, beam_summary, txru_correlation, txru_summary


def _variants() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for scheme in SCHEMES:
        rows.append({"variant_id": f"{scheme.lower()}__ideal", "transmission": {"tx_scheme": scheme},
                     "channel_estimation": {"ce_method": "IDEAL"}})
        methods = ("PLAN039_PRG_COMMON_REFERENCE_PDP",) if scheme == "BEAM8_PRECODER_CYCLING" else CDD_METHODS
        for method in methods:
            rows.append({"variant_id": f"{scheme.lower()}__{method.lower()}",
                         "transmission": {"tx_scheme": scheme}, "channel_estimation": {"ce_method": method}})
    return rows


def _write_validate_plots(
    validate_dir: Path,
    rays: dict[str, np.ndarray],
    aod_grid: np.ndarray,
    zod_grid: np.ndarray,
    reference_gain: np.ndarray,
    narrow_gain: np.ndarray,
    correlation10: np.ndarray,
    correlation25: np.ndarray,
    cut_angle_deg: np.ndarray,
    reference_horizontal_cut: np.ndarray,
    narrow_horizontal_cuts: np.ndarray,
    reference_vertical_cut: np.ndarray,
    narrow_vertical_cuts: np.ndarray,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    size = 8.0 + 240.0 * rays["power"] / np.max(rays["power"])
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    scatter = ax.scatter(rays["aod_deg"], rays["zod_deg"], s=size, c=rays["power"], cmap="viridis")
    ax.add_patch(plt.Rectangle((-60, 90), 120, 20, fill=False, edgecolor="red", linewidth=1.5))
    ax.set(xlabel="AoD (deg)", ylabel="ZoD (deg)", title="CDL-C ASD25 ray power")
    fig.colorbar(scatter, ax=ax, label="Linear ray power")
    fig.savefig(validate_dir / "angle_power_asd25.png", dpi=160)
    plt.close(fig)

    extent = [float(aod_grid[0]), float(aod_grid[-1]), float(zod_grid[0]), float(zod_grid[-1])]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    for ax, values, title in ((axes[0], reference_gain, "Reference ultrawide"),
                              (axes[1], np.max(narrow_gain, axis=2), "Maximum of 8 narrow beams")):
        db = 10.0 * np.log10(np.maximum(values / np.max(values), 1e-8))
        image = ax.imshow(db, origin="lower", extent=extent, aspect="auto", vmin=-30, vmax=0)
        ax.set(xlabel="AoD (deg)", ylabel="ZoD (deg)", title=title)
        fig.colorbar(image, ax=ax, label="Peak-normalized gain (dB)")
    fig.savefig(validate_dir / "direction_patterns.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, values, title in ((axes[0], correlation10, "ASD10"), (axes[1], correlation25, "ASD25")):
        image = ax.imshow(values, origin="lower", vmin=0, vmax=1, cmap="magma")
        ax.set(xlabel="Beam", ylabel="Beam", title=title)
        fig.colorbar(image, ax=ax, label="|correlation|")
    fig.savefig(validate_dir / "beam_correlation_heatmaps.png", dpi=160)
    plt.close(fig)

    def polar_codebook(path: Path, normalize_each: bool) -> None:
        fig = plt.figure(figsize=(14, 6), constrained_layout=True)
        axes = [fig.add_subplot(1, 2, index + 1, projection="polar") for index in range(2)]
        panels = (
            (axes[0], reference_horizontal_cut, narrow_horizontal_cuts,
             "Horizontal cut: fixed ZoD=90 deg", "AoD (deg)"),
            (axes[1], reference_vertical_cut, narrow_vertical_cuts,
             "Vertical cut: fixed AoD=0 deg", "Elevation from horizon (deg)"),
        )
        common_peak = max(
            float(np.max(reference_horizontal_cut)), float(np.max(narrow_horizontal_cuts)),
            float(np.max(reference_vertical_cut)), float(np.max(narrow_vertical_cuts)),
        )
        for ax, reference_cut, narrow_cuts, title, label in panels:
            reference_scale = float(np.max(reference_cut)) if normalize_each else common_peak
            reference_db = 10.0 * np.log10(np.maximum(reference_cut / reference_scale, 1e-6))
            ax.plot(np.deg2rad(cut_angle_deg), reference_db, color="black", linewidth=3,
                    label="reference ultrawide")
            for beam_index in range(8):
                beam_scale = float(np.max(narrow_cuts[:, beam_index])) if normalize_each else common_peak
                beam_db = 10.0 * np.log10(np.maximum(narrow_cuts[:, beam_index] / beam_scale, 1e-6))
                ax.plot(np.deg2rad(cut_angle_deg), beam_db, linewidth=1.2, label=f"narrow {beam_index}")
            ax.set_thetamin(-90)
            ax.set_thetamax(90)
            ax.set_rlim(-50, 0)
            ax.set_title(title, pad=20)
            ax.set_xlabel(label)
            ax.grid(True, alpha=0.35)
        axes[0].legend(fontsize=8, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.18))
        scale_label = "per-beam peak normalized" if normalize_each else "common absolute peak normalized"
        fig.suptitle(f"Plan-039 beam codebook polar cuts ({scale_label})")
        fig.savefig(path, dpi=180)
        plt.close(fig)

    polar_codebook(validate_dir / "codebook_polar_cuts_peak_normalized.png", normalize_each=True)
    polar_codebook(validate_dir / "codebook_polar_cuts_common_absolute.png", normalize_each=False)


def _merge_equal_delay_pdp(delays_s: np.ndarray, powers: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(np.asarray(delays_s), kind="stable")
    merged_delays: list[float] = []
    merged_powers: list[float] = []
    for delay, power in zip(np.asarray(delays_s)[order], np.asarray(powers)[order]):
        if merged_delays and abs(float(delay) - merged_delays[-1]) <= 1e-15:
            merged_powers[-1] += float(power)
        else:
            merged_delays.append(float(delay))
            merged_powers.append(float(power))
    return np.asarray(merged_delays), np.asarray(merged_powers)


def _write_codebook_pdp(
    validate_dir: Path,
    path_delays_s: np.ndarray,
    reference_path_powers: np.ndarray,
    narrow_path_powers: np.ndarray,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    series: list[tuple[str, np.ndarray, np.ndarray]] = []
    reference_delays, reference_powers = _merge_equal_delay_pdp(path_delays_s, reference_path_powers)
    series.append(("reference_ultrawide", reference_delays, reference_powers))
    for beam_index in range(8):
        delays, powers = _merge_equal_delay_pdp(path_delays_s, narrow_path_powers[beam_index])
        series.append((f"narrow_{beam_index}", delays, powers))
    strongest_tap_power = max(float(np.max(powers)) for _, _, powers in series)
    common_delay_zero = min(float(np.min(delays)) for _, delays, _ in series)
    common_delay_max_ns = max(float(np.max(delays) - common_delay_zero) * 1e9 for _, delays, _ in series)
    rows: list[dict[str, Any]] = []
    fig, axes = plt.subplots(3, 3, figsize=(15, 10), sharex=True, sharey=True, constrained_layout=True)
    for ax, (name, delays, powers) in zip(axes.ravel(), series):
        delay_ns = (delays - common_delay_zero) * 1e9
        markerline, stemlines, baseline = ax.stem(delay_ns, powers, basefmt=" ")
        markerline.set_markersize(3.5)
        stemlines.set_linewidth(1.0)
        ax.set_title(name)
        ax.set_xlim(-0.02 * max(common_delay_max_ns, 1.0), 1.02 * max(common_delay_max_ns, 1.0))
        ax.set_ylim(0.0, 1.05 * strongest_tap_power)
        ax.grid(True, alpha=0.3)
        for tap_index, (delay, delay_from_zero_ns, power) in enumerate(zip(delays, delay_ns, powers)):
            rows.append({"beam": name, "tap_index": tap_index, "delay_s": float(delay),
                         "delay_from_common_zero_ns": float(delay_from_zero_ns),
                         "absolute_power": float(power),
                         "common_y_scale_max_absolute_power": strongest_tap_power})
    fig.supxlabel("Delay from common zero (ns)")
    fig.supylabel("Absolute beam-domain PDP tap power (linear)")
    fig.suptitle(
        "Plan-039 reference and narrow-beam PDPs\n"
        f"Shared y-axis; strongest tap={strongest_tap_power:.6g}")
    fig.savefig(validate_dir / "codebook_pdp_absolute_subplots.png", dpi=180)
    plt.close(fig)
    rows_to_csv(validate_dir / "codebook_pdp_absolute.csv", rows)
    _json(validate_dir / "codebook_pdp_audit.json", {
        "schema": "plan039-codebook-pdp-audit-v1",
        "power_definition": "absolute beam-domain path power before per-beam normalization",
        "equal_delay_components_merged": True,
        "common_delay_zero_s": common_delay_zero,
        "common_x_max_delay_ns": common_delay_max_ns,
        "shared_y_min_absolute_power": 0.0,
        "shared_y_max_absolute_power": 1.05 * strongest_tap_power,
        "strongest_tap_absolute_power": strongest_tap_power,
        "beam_count": len(series),
    })


def audit_smoke(root: Path) -> dict[str, Any]:
    smoke_dir = root / "smoke" / SCENARIO / "stage0_smoke"
    with (smoke_dir / "summary.csv").open("r", encoding="utf-8", newline="") as handle:
        summary = list(csv.DictReader(handle))
    with (smoke_dir / "trial_metrics.csv").open("r", encoding="utf-8", newline="") as handle:
        trials = list(csv.DictReader(handle))
    expected_variants = {row["variant_id"] for row in _variants()}
    actual_variants = {row["variant_id"] for row in trials}
    estimated = [row for row in trials if row["receiver_mode"] == "estimated"]
    group_counts = {(variant, snr): 0 for variant in expected_variants for snr in (4.0, 12.0)}
    paired: dict[tuple[float, int], set[tuple[str, str, str, str]]] = {}
    for row in trials:
        key = (row["variant_id"], float(row["snr_db"]))
        group_counts[key] += 1
        pair_key = (float(row["snr_db"]), int(row["trial"]))
        paired.setdefault(pair_key, set()).add((row["channel_realization_index"], row["payload_seed"],
                                                row["data_noise_seed"], row["pilot_noise_seed"]))
    checks = {
        "summary_rows_16": len(summary) == 16,
        "trial_rows_320": len(trials) == 320,
        "eight_variants": actual_variants == expected_variants,
        "twenty_trials_per_variant_snr": all(value == 20 for value in group_counts.values()),
        "paired_seed_replay": all(len(values) == 1 for values in paired.values()),
        "estimated_nmse_finite": all(math.isfinite(float(row["ce_nmse_eff"])) for row in estimated),
        "filter_condition_finite": all(math.isfinite(float(row["filter_condition_number"])) for row in estimated),
        "unit_precoder_power": all(abs(float(row["precoder_power_min"])-1.0) <= 1e-12
                                     and abs(float(row["precoder_power_max"])-1.0) <= 1e-12 for row in trials),
    }
    report = {"schema": "plan039-stage0-smoke-audit-v1", "status": "PASS" if all(checks.values()) else "FAIL",
              "no_performance_conclusion": True, "checks": checks, "summary_rows": len(summary),
              "trial_rows": len(trials), "estimated_trial_rows": len(estimated)}
    _json(smoke_dir / "stage0_smoke_audit.json", report)
    if report["status"] != "PASS":
        raise RuntimeError(f"Plan-039 smoke audit failed: {checks}.")
    return report


def prepare_stage0(raw: dict[str, Any], config_path: Path, root: Path) -> Path:
    _validate_raw(raw)
    validate_dir = root / "validate" / SCENARIO
    validate_dir.mkdir(parents=True, exist_ok=True)
    configs, channels, audits = {}, {}, {}
    grid = None
    for target in (10.0, 25.0):
        cfg = config_from_dict(_scenario_payload(raw, target, selection_only=True))
        grid = build_resource_grid(cfg.resource)
        channel = FixedCDLStatisticsChannel(cfg, grid, prepare_reference=False)
        rays = _ray_table(channel)
        configs[target], channels[target] = cfg, channel
        audits[f"ASD{int(target)}"] = _angle_audit(rays, target)
        rows_to_csv(validate_dir / f"ray_table_asd{int(target)}.csv", [
            {"cluster_index": int(rays["cluster"][i]), "ray_index": int(rays["ray"][i]),
             "kind": str(rays["kind"][i]), "delay_s": float(rays["delay_s"][i]),
             "aod_deg": float(rays["aod_deg"][i]), "zod_deg": float(rays["zod_deg"][i]),
             "linear_power": float(rays["power"][i])} for i in range(rays["power"].size)])
    _json(validate_dir / "angle_audit.json", audits)
    if audits["ASD25"]["status"] != "PASS":
        raise RuntimeError(f"Plan-039 angle gate failed: {audits['ASD25']}.")
    assert grid is not None

    reference, narrow, codebook = build_angular_full_coverage_codebook(
        tuple(float(value) for value in raw["beam_design"]["regularization_grid"]),
        int(raw["beam_design"]["horizontal_samples"]), int(raw["beam_design"]["vertical_samples"]))
    cfg25, channel25 = configs[25.0], channels[25.0]
    fixed = cfg25.fixed_cdl_statistics
    aod_grid = np.linspace(-60.0, 60.0, int(raw["beam_design"]["horizontal_samples"]))
    zod_grid = np.linspace(90.0, 110.0, int(raw["beam_design"]["vertical_samples"]))
    mesh_aod, mesh_zod = np.meshgrid(aod_grid, zod_grid, indexing="xy")
    uh, uv = direction_cosines(np.deg2rad(mesh_aod.ravel()), np.deg2rad(mesh_zod.ravel()))
    response = txru_response_matrix(uh, uv, fixed.bs_vertical_aes, fixed.bs_horizontal_aes,
                                    fixed.bs_vertical_spacing_lambda, fixed.bs_horizontal_spacing_lambda,
                                    channel25.txru_mapping, fixed.bs_antenna_pattern)
    reference_gain, narrow_gain = np.abs(response @ reference) ** 2, np.abs(response @ narrow) ** 2
    cut_angle_deg = np.linspace(-90.0, 90.0, 361)
    horizontal_uh, horizontal_uv = direction_cosines(
        np.deg2rad(cut_angle_deg), np.full(cut_angle_deg.size, np.deg2rad(90.0)))
    vertical_zod_deg = 90.0 - cut_angle_deg
    vertical_uh, vertical_uv = direction_cosines(
        np.zeros(cut_angle_deg.size), np.deg2rad(vertical_zod_deg))
    horizontal_response = txru_response_matrix(
        horizontal_uh, horizontal_uv, fixed.bs_vertical_aes, fixed.bs_horizontal_aes,
        fixed.bs_vertical_spacing_lambda, fixed.bs_horizontal_spacing_lambda,
        channel25.txru_mapping, fixed.bs_antenna_pattern)
    vertical_response = txru_response_matrix(
        vertical_uh, vertical_uv, fixed.bs_vertical_aes, fixed.bs_horizontal_aes,
        fixed.bs_vertical_spacing_lambda, fixed.bs_horizontal_spacing_lambda,
        channel25.txru_mapping, fixed.bs_antenna_pattern)
    reference_horizontal_cut = np.abs(horizontal_response @ reference) ** 2
    narrow_horizontal_cuts = np.abs(horizontal_response @ narrow) ** 2
    reference_vertical_cut = np.abs(vertical_response @ reference) ** 2
    narrow_vertical_cuts = np.abs(vertical_response @ narrow) ** 2
    codebook.update({"target_grid_shape": [int(zod_grid.size), int(aod_grid.size)],
                     "reference_grid_all_finite": bool(np.all(np.isfinite(reference_gain))),
                     "reference_grid_zero_count": int(np.sum(reference_gain <= 0.0)),
                     "reference_in_band_min_linear": float(np.min(reference_gain)),
                     "reference_in_band_mean_linear": float(np.mean(reference_gain)),
                     "reference_in_band_p5_linear": float(np.percentile(reference_gain, 5)),
                     "reference_in_band_p95_linear": float(np.percentile(reference_gain, 95)),
                     "reference_in_band_p95_p5_ripple_linear": float(np.percentile(reference_gain, 95)-np.percentile(reference_gain, 5))})
    codebook["status"] = "PASS" if (codebook["reference_norm_error"] <= 1e-12
        and codebook["narrow_gram_max_error"] <= 1e-12
        and max(abs(np.asarray(codebook["horizontal_center_aod_deg"]))) <= 61.5
        and codebook["reference_grid_all_finite"] and codebook["reference_grid_zero_count"] == 0) else "FAIL"
    _json(validate_dir / "codebook_audit.json", codebook)
    if codebook["status"] != "PASS":
        raise RuntimeError("Plan-039 codebook gate failed.")
    np.savez_compressed(validate_dir / "direction_pattern_inputs.npz", aod_deg=aod_grid, zod_deg=zod_grid,
                        reference_gain=reference_gain.reshape(zod_grid.size, aod_grid.size),
                        narrow_gain=narrow_gain.reshape(zod_grid.size, aod_grid.size, 8),
                        reference_weight=reference, narrow_weights=narrow,
                        cut_angle_deg=cut_angle_deg,
                        horizontal_cut_fixed_zod_deg=np.asarray([90.0]),
                        reference_horizontal_cut=reference_horizontal_cut,
                        narrow_horizontal_cuts=narrow_horizontal_cuts,
                        vertical_cut_fixed_aod_deg=np.asarray([0.0]),
                        vertical_cut_zod_deg=vertical_zod_deg,
                        reference_vertical_cut=reference_vertical_cut,
                        narrow_vertical_cuts=narrow_vertical_cuts)

    correlation: dict[str, Any] = {}
    beam_correlation_matrices: dict[float, np.ndarray] = {}
    for target in (10.0, 25.0):
        matrix, summary, txru_matrix, txru_summary = _correlation(
            channels[target], narrow, int(raw["plan039"]["correlation_realizations"]))
        np.save(validate_dir / f"beam_correlation_asd{int(target)}.npy", matrix)
        np.save(validate_dir / f"txru_correlation_asd{int(target)}.npy", txru_matrix)
        correlation[f"ASD{int(target)}"] = {"beam_domain": summary, "txru_domain": txru_summary}
        beam_correlation_matrices[target] = matrix
    correlation["gate"] = {"criterion": "ASD25 off-diagonal mean < ASD10 off-diagonal mean",
                           "passed": correlation["ASD25"]["beam_domain"]["off_diagonal_mean"]
                           < correlation["ASD10"]["beam_domain"]["off_diagonal_mean"]}
    _json(validate_dir / "correlation_audit.json", correlation)
    if not correlation["gate"]["passed"]:
        raise RuntimeError("ANGULAR_SCALING_NO_CORRELATION_REDUCTION")
    _write_validate_plots(validate_dir, _ray_table(channel25), aod_grid, zod_grid,
                          reference_gain.reshape(zod_grid.size, aod_grid.size),
                          narrow_gain.reshape(zod_grid.size, aod_grid.size, 8),
                          beam_correlation_matrices[10.0], beam_correlation_matrices[25.0],
                          cut_angle_deg, reference_horizontal_cut, narrow_horizontal_cuts,
                          reference_vertical_cut, narrow_vertical_cuts)

    sidon = search_strict_sidon_top8()
    sidon["selection_status"] = "GEOMETRY_TOP8_ONLY_NO_PERFORMANCE_SELECTION"
    selected_for_smoke = sidon["top8"][0]["delay_indices"]
    _json(validate_dir / "sidon_rule_inference.json", sidon)
    rows_to_csv(validate_dir / "sidon_top8.csv", [
        {"geometric_rank": index+1, "delay_indices": ",".join(map(str, row["delay_indices"])),
         "fold_gap": row["fold_gap"], "pair_gap": row["pair_gap"], "pilot_rank": row["pilot_rank"],
         "pilot_condition_number": row["pilot_condition_number"]} for index, row in enumerate(sidon["top8"])])

    ray_delays, _, ray_covariances = channel25.analytic_ray_covariances()
    beam_power, _ = beam_path_statistics(narrow, ray_covariances)
    reference_power, _ = beam_path_statistics(reference[:, None], ray_covariances)
    _write_codebook_pdp(validate_dir, ray_delays, reference_power[0], beam_power)
    time_covariance = cdl_spatial_unaware_covariance(grid, cfg25.channel).time
    frequencies = np.arange(grid.n_sc, dtype=np.float64) * float(grid.scs_khz) * 1e3
    covariance_reports: dict[str, Any] = {}
    covariance_arrays = {"path_delays_s": ray_delays, "reference_path_powers": reference_power[0],
                         "beam_path_powers": beam_power}
    for name, delays in (("B0_QC", [0,9,18,27,36,45,54,63]),
                         ("SIDON_GEOMETRY_TOP1_SMOKE", selected_for_smoke)):
        alpha = np.full(grid.n_sc, 1.0/np.sqrt(8.0))
        covariance_reports[name] = {method: covariance_audit(
            frequency_covariance(method, frequencies, covariance_arrays, delays, alpha)) for method in CDD_METHODS}
    _json(validate_dir / "covariance_audit.json", covariance_reports)

    cdd, alpha, orthogonal = beam_domain_cdd_precoder(narrow, selected_for_smoke, np.arange(576), 576)
    rng = np.random.default_rng(2026093907)
    h = rng.normal(size=(2,576,32)) + 1j*rng.normal(size=(2,576,32))
    direct = np.einsum("rkt,kt->rk", h, cdd)
    branches = np.einsum("rkt,tb->rkb", h, narrow)
    phase = np.exp(-1j*2*np.pi*np.arange(576)[:,None]*np.asarray(selected_for_smoke)[None,:]/576.0)
    projected = np.einsum("rkb,kb->rk", branches, alpha[:,None]*phase)
    synthesis = {"orthogonal_branches": orthogonal,
                 "unit_power_max_error": float(np.max(np.abs(np.sum(np.abs(cdd)**2, axis=1)-1.0))),
                 "direct_txru_vs_beam_projection_max_error": float(np.max(np.abs(direct-projected)))}
    synthesis["status"] = "PASS" if synthesis["unit_power_max_error"] <= 1e-12 and synthesis[
        "direct_txru_vs_beam_projection_max_error"] <= 1e-10 else "FAIL"
    _json(validate_dir / "synthesis_audit.json", synthesis)
    if synthesis["status"] != "PASS":
        raise RuntimeError(f"Plan-039 synthesis audit failed: {synthesis}.")

    arrays_path = validate_dir / "plan039_numeric.npz"
    np.savez_compressed(arrays_path, narrow_weights=narrow, wide_weight=reference,
                        path_delays_s=ray_delays, reference_path_powers=reference_power[0],
                        beam_path_powers=beam_power,
                        sidon_selected_delay_indices=np.asarray(selected_for_smoke, dtype=np.int64),
                        time_covariance=time_covariance)
    manifest_path = validate_dir / "plan039_manifest.json"
    manifest = {"schema": "plan039-angular-full-coverage-manifest-v1", "status": "FROZEN", "plan": PLAN,
        "stage": "stage0_smoke_only", "scenario_id": SCENARIO,
        "frozen_context": {"cdl_profile": "C", "delay_spread_ns": 100.0,
            "carrier_frequency_hz": float(cfg25.channel.carrier_frequency_hz), "ue_speed_kmh": 3.0,
            "n_sc": int(grid.n_sc), "n_symbols": int(grid.n_symbols),
            "statistics_seed": int(cfg25.fixed_cdl_statistics.statistics_seed)},
        "arrays_file": arrays_path.name, "arrays_sha256": _sha256(arrays_path),
        "selected_ssb_index": 0, "ssb_long_term_powers": [float(np.sum(reference_power[0]))],
        "reference_receive_power": float(np.sum(reference_power[0])),
        "reference_pdp": {"delays_s": ray_delays.tolist(), "powers": reference_power[0].tolist(),
                          "time_covariance_real": time_covariance.real.tolist(),
                          "time_covariance_imag": time_covariance.imag.tolist()},
        "sidon_status": "GEOMETRY_TOP1_FOR_STAGE0_SMOKE_NOT_STAGE1_SELECTED",
        "sidon_delay_indices": selected_for_smoke, "source_config": config_path.resolve().as_posix()}
    digest = write_manifest(manifest_path, manifest)
    rendered = dataclass_to_dict(cfg25)
    rendered["fixed_cdl_statistics"].update({"selection_only": False,
        "frozen_beam_manifest": manifest_path.resolve().as_posix(), "frozen_beam_manifest_sha256": digest})
    rendered["variants"] = _variants()
    rendered["simulation"].update({"snr_points_db": list(raw["plan039"]["smoke"]["snr_points_db"]),
        "n_trials_per_snr": int(raw["plan039"]["smoke"]["trials_per_snr"]),
        "max_trials_per_snr": int(raw["plan039"]["smoke"]["trials_per_snr"]), "min_block_errors": 0,
        "output_dir": (root/"smoke"/SCENARIO).as_posix(), "run_id": "stage0_smoke",
        "save_trial_metrics": True, "common_random_numbers": True, "absolute_trial_start": 0})
    rendered["scenarios"] = [{"scenario_id": SCENARIO}]
    smoke_config = root / "configs" / f"{SCENARIO}_stage0_smoke.yaml"
    smoke_config.parent.mkdir(parents=True, exist_ok=True)
    smoke_config.write_text(yaml.safe_dump(rendered, sort_keys=False, allow_unicode=True), encoding="utf-8")
    _json(validate_dir / "validate_report.json", {"schema": "plan039-stage0-validate-v1", "status": "PASS",
        "scenario": SCENARIO, "no_performance_conclusion": True, "curves": len(_variants()),
        "sidon_search_mode": "rule_inference_not_numerical_simulation", "manifest_sha256": digest,
        "smoke_config": smoke_config.as_posix(),
        "gates": {"angles": audits["ASD25"]["status"], "codebook": codebook["status"],
                  "correlation": correlation["gate"]["passed"], "synthesis": synthesis["status"]}})
    return smoke_config


def _stage0_sources(root: Path) -> tuple[Path, Path, dict[str, Any], dict[str, np.ndarray]]:
    audit_path = root / "smoke" / SCENARIO / "stage0_smoke" / "stage0_smoke_audit.json"
    if not audit_path.exists() or json.loads(audit_path.read_text(encoding="utf-8"))["status"] != "PASS":
        raise RuntimeError("Stage 0 PASS audit is required before stage 1A.")
    source_dir = root / "validate" / SCENARIO
    manifest_path = source_dir / "plan039_manifest.json"
    arrays_path = source_dir / "plan039_numeric.npz"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    with np.load(arrays_path) as loaded:
        arrays = {name: np.asarray(loaded[name]) for name in loaded.files}
    return manifest_path, arrays_path, manifest, arrays


def _candidate_manifests(root: Path, raw: dict[str, Any]) -> list[dict[str, Any]]:
    _, _, source_manifest, source_arrays = _stage0_sources(root)
    top8 = search_strict_sidon_top8()["top8"]
    ranks = [int(value) for value in raw["plan039"].get("selection_candidate_ranks", [1])]
    if not ranks or len(set(ranks)) != len(ranks) or any(rank < 1 or rank > len(top8) for rank in ranks):
        raise ValueError("plan039.selection_candidate_ranks must contain unique ranks in the geometric top-8.")
    candidates = [(rank, top8[rank-1]) for rank in ranks]
    expected = raw["plan039"].get("selection_candidate_delay_indices")
    actual = [candidate["delay_indices"] for _, candidate in candidates]
    if expected is not None and [[int(item) for item in row] for row in expected] != actual:
        raise RuntimeError(f"Configured Sidon prescan candidates do not match geometric search: {actual}.")
    records: list[dict[str, Any]] = []
    for rank, candidate in candidates:
        directory = root / "selection" / "manifests" / f"sidon_{rank:02d}"
        arrays_path = directory / "plan039_numeric.npz"
        manifest_path = directory / "plan039_manifest.json"
        if arrays_path.exists() and manifest_path.exists():
            with np.load(arrays_path) as loaded:
                frozen_delays = np.asarray(loaded["sidon_selected_delay_indices"], dtype=np.int64).tolist()
            frozen_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            if frozen_delays != candidate["delay_indices"] or frozen_manifest.get("arrays_sha256") != _sha256(arrays_path):
                raise RuntimeError(f"Existing candidate artifact is inconsistent: {directory}.")
            digest = str(frozen_manifest.get("manifest_sha256", ""))
            if len(digest) != 64:
                raise RuntimeError(f"Existing candidate manifest has no valid hash: {manifest_path}.")
            load_plan039_manifest(manifest_path, digest)
            records.append({"rank": rank, "candidate": candidate, "manifest": manifest_path.resolve(),
                            "manifest_sha256": digest})
            continue
        arrays = dict(source_arrays)
        arrays["sidon_selected_delay_indices"] = np.asarray(candidate["delay_indices"], dtype=np.int64)
        directory.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(arrays_path, **arrays)
        manifest = dict(source_manifest)
        manifest.update({"stage": "stage1a_selection_candidate", "sidon_status": "STAGE1A_CANDIDATE",
                         "sidon_geometric_rank": rank, "sidon_delay_indices": candidate["delay_indices"],
                         "arrays_file": arrays_path.name, "arrays_sha256": _sha256(arrays_path)})
        digest = write_manifest(manifest_path, manifest)
        records.append({"rank": rank, "candidate": candidate, "manifest": manifest_path.resolve(),
                        "manifest_sha256": digest})
    return records


def _stage_base(raw: dict[str, Any], manifest: Path, digest: str, seed: int,
                realization_seed: int) -> dict[str, Any]:
    cfg = config_from_dict(_scenario_payload(raw, 25.0, selection_only=True))
    rendered = dataclass_to_dict(cfg)
    rendered["fixed_cdl_statistics"].update({
        "selection_only": False, "frozen_beam_manifest": manifest.as_posix(),
        "frozen_beam_manifest_sha256": digest, "realization_seed": int(realization_seed)})
    rendered["simulation"].update({"seed": int(seed), "save_trial_metrics": True,
                                    "common_random_numbers": True, "min_block_errors": 0})
    rendered["scenarios"] = [{"scenario_id": SCENARIO}]
    rendered["plots"] = {"enabled": False}
    return rendered


def _selection_variants(base_manifest: Path, base_digest: str,
                        candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    variants: list[dict[str, Any]] = []
    fixed_base = {"frozen_beam_manifest": base_manifest.as_posix(),
                  "frozen_beam_manifest_sha256": base_digest}
    for scheme, name in (("BEAM8_B0_QC", "b0"), ("BEAM8_PRECODER_CYCLING", "cycling")):
        variants.append({"variant_id": f"{name}__ideal", "transmission": {"tx_scheme": scheme},
                         "channel_estimation": {"ce_method": "IDEAL"},
                         "fixed_cdl_statistics": fixed_base})
    for record in candidates:
        fixed = {"frozen_beam_manifest": record["manifest"].as_posix(),
                 "frozen_beam_manifest_sha256": record["manifest_sha256"]}
        name = f"sidon_{record['rank']:02d}"
        variants.append({"variant_id": f"{name}__ideal", "transmission": {"tx_scheme": "BEAM8_SIDON_SELECTED"},
                         "channel_estimation": {"ce_method": "IDEAL"}, "fixed_cdl_statistics": fixed})
    for method in CDD_METHODS:
        variants.append({"variant_id": f"b0__{method.lower()}", "transmission": {"tx_scheme": "BEAM8_B0_QC"},
                         "channel_estimation": {"ce_method": method}, "fixed_cdl_statistics": fixed_base,
                         "simulation": {"ce_only": True}})
    for record in candidates:
        fixed = {"frozen_beam_manifest": record["manifest"].as_posix(),
                 "frozen_beam_manifest_sha256": record["manifest_sha256"]}
        for method in CDD_METHODS:
            variants.append({"variant_id": f"sidon_{record['rank']:02d}__{method.lower()}",
                             "transmission": {"tx_scheme": "BEAM8_SIDON_SELECTED"},
                             "channel_estimation": {"ce_method": method}, "fixed_cdl_statistics": fixed,
                             "simulation": {"ce_only": True}})
    variants.append({"variant_id": "cycling__plan039_prg_common_reference_pdp",
                     "transmission": {"tx_scheme": "BEAM8_PRECODER_CYCLING"},
                     "channel_estimation": {"ce_method": "PLAN039_PRG_COMMON_REFERENCE_PDP"},
                     "fixed_cdl_statistics": fixed_base, "simulation": {"ce_only": True}})
    return variants


def _rows_for(rows: list[dict[str, Any]], variant: str) -> list[dict[str, Any]]:
    return [row for row in rows if row["variant_id"] == variant]


def _selection_extension(rows: list[dict[str, Any]]) -> list[float]:
    ideal = sorted({str(row["variant_id"]) for row in rows if str(row["variant_id"]).endswith("__ideal")})
    current = sorted({float(row["snr_db"]) for row in rows})
    requested: set[float] = set()
    for variant in ideal:
        curve = _rows_for(rows, variant)
        crossing, _ = _crossing(curve, 0.10)
        if math.isfinite(crossing):
            continue
        values = [float(row["bler"]) for row in curve]
        if values and all(value > 0.10 for value in values) and max(current) < SNR_MAX:
            requested.add(min(SNR_MAX, max(current) + 2.0))
        elif values and all(value < 0.10 for value in values) and min(current) > SNR_MIN:
            requested.add(max(SNR_MIN, min(current) - 2.0))
    return sorted(requested)


def _selection_grid(rows: list[dict[str, Any]], selected_variant: str) -> list[float]:
    crossings: list[float] = []
    selected_curves = ("b0__ideal", selected_variant, "cycling__ideal")
    missing_one_percent = False
    for variant in selected_curves:
        curve = _rows_for(rows, variant)
        for target in (0.10, 0.01):
            value, _ = _crossing(curve, target)
            if math.isfinite(value):
                crossings.append(value)
            elif target == 0.01:
                missing_one_percent = True
    if not crossings:
        raise RuntimeError("Cannot construct evaluation grid without a real selection crossing.")
    lower = max(SNR_MIN, math.floor((min(crossings) - 1.0) * 2.0) / 2.0)
    upper = min(SNR_MAX, math.ceil((max(crossings) + 1.0) * 2.0) / 2.0)
    if missing_one_percent:
        measured = [float(row["snr_db"]) for variant in selected_curves for row in _rows_for(rows, variant)]
        upper = min(SNR_MAX, max(measured))
    grid = set(np.arange(lower, upper + 0.25, 0.5).round(8).tolist())
    for crossing in crossings:
        start = max(lower, math.floor((crossing - 0.5) * 4.0) / 4.0)
        stop = min(upper, math.ceil((crossing + 0.5) * 4.0) / 4.0)
        grid.update(np.arange(start, stop + 0.125, 0.25).round(8).tolist())
    return sorted(float(value) for value in grid)


def run_selection(raw: dict[str, Any], root: Path) -> dict[str, Any]:
    base_manifest, _, base_payload, _ = _stage0_sources(root)
    base_digest = str(base_payload["manifest_sha256"])
    candidates = _candidate_manifests(root, raw)
    base_seed = int(raw["simulation"]["seed"])
    selection_seed = _stable_seed(base_seed, "selection")
    selection_channel_seed = _stable_seed(base_seed, "selection_channel")
    template = _stage_base(raw, base_manifest.resolve(), base_digest, selection_seed, selection_channel_seed)
    template["variants"] = _selection_variants(base_manifest.resolve(), base_digest, candidates)
    batch_dirs: list[Path] = []
    points = list(SELECTION_GRID)
    batch_index = 0
    while points:
        batch_index += 1
        config = copy.deepcopy(template)
        run_id = f"batch_{batch_index:02d}"
        config["simulation"].update({"snr_points_db": points, "n_trials_per_snr": SELECTION_TRIALS,
                                     "max_trials_per_snr": SELECTION_TRIALS, "absolute_trial_start": 0,
                                     "output_dir": (root / "selection" / "batches").as_posix(), "run_id": run_id})
        directory = root / "selection" / "batches" / run_id
        batch_dirs.append(_run_batch(config, root / "configs" / f"stage1a_{run_id}.yaml", directory))
        merged = _merge_summary(batch_dirs)
        points = _selection_extension(merged)
        already = {float(row["snr_db"]) for row in merged}
        points = [point for point in points if point not in already]
    summary = _merge_summary(batch_dirs)
    trials = _collect_trials(batch_dirs)
    selection_dir = root / "selection" / SCENARIO
    _write_csv(selection_dir / "summary.csv", summary)
    _write_csv(selection_dir / "trial_metrics.csv", trials)
    rankings: list[dict[str, Any]] = []
    b0_method2 = _rows_for(summary, "b0__plan039_beam_specific_pdp_independent")
    b0_nmse = {float(row["snr_db"]): float(row["ce_nmse_eff"]) for row in b0_method2}
    for record in candidates:
        name = f"sidon_{record['rank']:02d}"
        c10, bracket10 = _crossing(_rows_for(summary, f"{name}__ideal"), 0.10)
        c01, bracket01 = _crossing(_rows_for(summary, f"{name}__ideal"), 0.01)
        method2 = _rows_for(summary, f"{name}__plan039_beam_specific_pdp_independent")
        penalties = [10.0 * math.log10(float(row["ce_nmse_eff"]) / b0_nmse[float(row["snr_db"])])
                     for row in method2 if float(row["snr_db"]) in b0_nmse]
        penalty = max(penalties) if penalties else float("inf")
        rankings.append({"geometric_rank": record["rank"], "variant_id": name,
                         "delay_indices": record["candidate"]["delay_indices"],
                         "crossing_10pct_db": c10, "crossing_1pct_db": c01,
                         "bracket_10pct_db": bracket10, "bracket_1pct_db": bracket01,
                         "method2_worst_nmse_penalty_db": penalty,
                         "manifest": record["manifest"].as_posix(),
                         "manifest_sha256": record["manifest_sha256"]})
    def finite_or_inf(value: Any) -> float:
        number = float(value)
        return number if math.isfinite(number) else float("inf")
    rankings.sort(key=lambda row: (finite_or_inf(row["crossing_10pct_db"]),
                                   finite_or_inf(row["crossing_1pct_db"]),
                                   finite_or_inf(row["method2_worst_nmse_penalty_db"]),
                                   int(row["geometric_rank"])))
    if not math.isfinite(float(rankings[0]["crossing_10pct_db"])):
        raise RuntimeError("No Sidon candidate obtained a real 10% crossing within [-6,26] dB.")
    selected = rankings[0]
    selected["selection_rank"] = 1
    grid = _selection_grid(summary, f"{selected['variant_id']}__ideal")
    grid_bytes = json.dumps(grid, separators=(",", ":")).encode("utf-8")
    seed_manifest = {"base": base_seed, "selection": selection_seed,
                     "selection_channel": selection_channel_seed,
                     "evaluation": _stable_seed(base_seed, "evaluation"),
                     "evaluation_channel": _stable_seed(base_seed, "evaluation_channel"),
                     "bootstrap": _stable_seed(base_seed, "bootstrap")}
    derived = [value for key, value in seed_manifest.items() if key != "base"]
    if len(set(derived)) != len(derived):
        raise RuntimeError("Plan-039 derived seed namespaces are not disjoint.")
    executed_ranks = {int(record["rank"]) for record in candidates}
    deferred = [{"geometric_rank": rank, "delay_indices": candidate["delay_indices"],
                 "status": "DEFERRED_NOT_RUN"}
                for rank, candidate in enumerate(search_strict_sidon_top8()["top8"], 1)
                if rank not in executed_ranks]
    freeze = {"schema": "plan039-stage1a-freeze-v1", "status": "FROZEN", "selected": selected,
              "executed_candidates": rankings, "deferred_candidates": deferred,
              "evaluation_snr_grid_db": grid,
              "prescan_candidate_ranks": [record["rank"] for record in candidates],
              "scope_note": "Researcher-requested limited prescan; additional Sidon candidates are deferred.",
              "evaluation_snr_grid_sha256": hashlib.sha256(grid_bytes).hexdigest(),
              "seeds": seed_manifest,
              "selection_summary_sha256": _sha256(selection_dir / "summary.csv")}
    _json(selection_dir / "selection_freeze.json", freeze)
    _write_csv(selection_dir / "candidate_ranking.csv", rankings)
    return freeze


def _evaluation_variants(freeze: dict[str, Any]) -> list[dict[str, Any]]:
    selected = freeze["selected"]
    selected_fixed = {"frozen_beam_manifest": selected["manifest"],
                      "frozen_beam_manifest_sha256": selected["manifest_sha256"]}
    variants: list[dict[str, Any]] = []
    for scheme, name, methods, fixed in (
        ("BEAM8_B0_QC", "b0", STAGE1B_CDD_METHODS, None),
        ("BEAM8_SIDON_SELECTED", "sidon_selected", STAGE1B_CDD_METHODS, selected_fixed),
        ("BEAM8_PRECODER_CYCLING", "cycling", (STAGE1B_CYCLING_METHOD,), None),
    ):
        for method in methods:
            row = {"variant_id": f"{name}__{method.lower()}", "transmission": {"tx_scheme": scheme},
                   "channel_estimation": {"ce_method": method}}
            if fixed is not None:
                row["fixed_cdl_statistics"] = fixed
            variants.append(row)
    return variants


def _adaptive_points(summary: list[dict[str, Any]], variants: tuple[str, ...]) -> dict[int, list[float]]:
    requested: dict[int, set[float]] = {}
    for variant in variants:
        curve = _rows_for(summary, variant)
        by_snr = {float(row["snr_db"]): row for row in curve}
        for target, lower, upper in ((0.10, 0.05, 0.20), (0.01, 0.005, 0.02)):
            _, bracket = _crossing(curve, target)
            if bracket is None:
                continue
            for snr in bracket:
                row = by_snr[snr]
                bler = float(row["bler"])
                trials = int(row["n_trials"])
                errors = int(row["tb_errors"])
                if lower <= bler <= upper and errors < 200 and trials < FORMAL_CAP:
                    requested.setdefault(trials, set()).add(snr)
    return {start: sorted(points) for start, points in sorted(requested.items())}


def _stage1b_grid(raw: dict[str, Any], family: str) -> list[float]:
    stage = raw["plan039"]["stage1b"]
    lower = float(stage["snr_min_db"])
    upper = float(stage["snr_max_db"])
    coarse = float(stage["coarse_step_db"])
    fine = float(stage["fine_step_db"])
    window = stage[f"{family}_refinement_db"]
    if not (lower < upper and coarse > 0.0 and fine > 0.0 and len(window) == 2):
        raise ValueError("Invalid plan039.stage1b SNR grid definition.")
    refine_lower, refine_upper = map(float, window)
    if not (lower <= refine_lower < refine_upper <= upper):
        raise ValueError(f"Invalid {family} refinement window {window!r}.")
    coarse_grid = np.arange(lower, upper + coarse / 2.0, coarse)
    fine_grid = np.arange(refine_lower, refine_upper + fine / 2.0, fine)
    return sorted({round(float(value), 8) for value in np.r_[coarse_grid, fine_grid]})


def _aged_mrt_grid(raw: dict[str, Any]) -> list[float]:
    stage = raw["plan039"]["aged_mrt"]
    lower = float(stage["snr_min_db"])
    upper = float(stage["snr_max_db"])
    step = float(stage["snr_step_db"])
    if not (lower == 7.5 and upper == 20.0 and step == 0.5):
        raise ValueError("Plan-039 aged MRT formal grid is frozen to 7.5:0.5:20 dB.")
    return [round(float(value), 8) for value in np.arange(lower, upper + step / 2.0, step)]


def _prepare_aged_mrt_manifest(
    raw: dict[str, Any],
    root: Path,
    speed_kmh: float,
) -> tuple[Path, str, dict[str, Any]]:
    source_path, _, source_manifest, source_arrays = _stage0_sources(root)
    directory = root / "aged_mrt_60kmh" / "manifest"
    manifest_path = directory / "plan039_manifest_60kmh.json"
    if manifest_path.exists():
        candidate = json.loads(manifest_path.read_text(encoding="utf-8"))
        candidate_digest = str(candidate.get("manifest_sha256", ""))
        candidate, _ = load_plan039_manifest(manifest_path, candidate_digest)
        frozen = candidate.get("frozen_context", {})
        if (
            candidate.get("stage") == "supplement_aged_mrt_60kmh_40ms"
            and candidate.get("source_manifest_sha256") == source_manifest["manifest_sha256"]
            and float(frozen.get("ue_speed_kmh", -1.0)) == float(speed_kmh)
        ):
            return manifest_path, candidate_digest, candidate
    config_payload = _scenario_payload(raw, 25.0, selection_only=True)
    config_payload["channel"]["ue_speed_kmh"] = float(speed_kmh)
    cfg = config_from_dict(config_payload)
    grid = build_resource_grid(cfg.resource)
    time_covariance = cdl_spatial_unaware_covariance(grid, cfg.channel).time
    directory.mkdir(parents=True, exist_ok=True)
    arrays_path = directory / "plan039_numeric_60kmh.npz"
    arrays = dict(source_arrays)
    arrays["time_covariance"] = time_covariance
    np.savez_compressed(arrays_path, **arrays)
    manifest = copy.deepcopy(source_manifest)
    manifest.update({
        "stage": "supplement_aged_mrt_60kmh_40ms",
        "source_manifest": source_path.resolve().as_posix(),
        "source_manifest_sha256": str(source_manifest["manifest_sha256"]),
        "arrays_file": arrays_path.name,
        "arrays_sha256": _sha256(arrays_path),
    })
    manifest["frozen_context"] = dict(manifest["frozen_context"])
    manifest["frozen_context"]["ue_speed_kmh"] = float(speed_kmh)
    manifest["reference_pdp"] = dict(manifest["reference_pdp"])
    manifest["reference_pdp"].update({
        "time_covariance_real": time_covariance.real.tolist(),
        "time_covariance_imag": time_covariance.imag.tolist(),
    })
    digest = write_manifest(manifest_path, manifest)
    return manifest_path, digest, manifest


def run_aged_mrt_formal(raw: dict[str, Any], root: Path) -> dict[str, Any]:
    stage = raw["plan039"]["aged_mrt"]
    speed_kmh = float(stage["speed_kmh"])
    age_ms = float(stage["csi_age_ms"])
    trials = int(stage["trials_per_snr"])
    max_workers = int(stage.get("max_parallel_snr", 1))
    if speed_kmh != 60.0 or age_ms != 40.0:
        raise ValueError("Plan-039 aged MRT supplement is frozen to 60 km/h and 40 ms CSI age.")
    if trials != FORMAL_BATCH or bool(stage.get("prescan", False)):
        raise ValueError("Plan-039 aged MRT runs no prescan and exactly 1,000 formal trials per SNR.")
    if max_workers < 1 or max_workers > 3:
        raise ValueError("Plan-039 aged MRT max_parallel_snr must be in [1,3].")
    grid = _aged_mrt_grid(raw)
    manifest_path, manifest_digest, manifest = _prepare_aged_mrt_manifest(raw, root, speed_kmh)
    base_seed = int(raw["simulation"]["seed"])
    evaluation_seed = _stable_seed(base_seed, "aged_mrt_60kmh_40ms_evaluation")
    channel_seed = _stable_seed(base_seed, "aged_mrt_60kmh_40ms_channel")
    template = _stage_base(raw, manifest_path.resolve(), manifest_digest, evaluation_seed, channel_seed)
    template["channel"]["ue_speed_kmh"] = speed_kmh
    variant_id = "aged_mrt_prg6__60kmh__csi_age40ms"
    template["variants"] = [{
        "variant_id": variant_id,
        "transmission": {"tx_scheme": AGED_MRT_SCHEME, "aged_csi_ms": age_ms},
        "channel_estimation": {"ce_method": AGED_MRT_METHOD},
    }]
    grid_digest = hashlib.sha256(
        json.dumps(grid, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    output = root / "aged_mrt_60kmh" / SCENARIO
    freeze = {
        "schema": "plan039-aged-mrt-formal-freeze-v1",
        "status": "FROZEN",
        "speed_kmh": speed_kmh,
        "csi_age_ms": age_ms,
        "csi_time_reference": "old=-40ms,current_slot_start=0ms",
        "mrt_definition": "dominant eigenvector of two-Rx PRG-summed Gram matrix",
        "prg_size_rb": 6,
        "receiver": AGED_MRT_METHOD,
        "prescan_run": False,
        "adaptive_additions": False,
        "snr_grid_db": grid,
        "snr_grid_sha256": grid_digest,
        "trials_per_snr": trials,
        "max_parallel_snr": max_workers,
        "variant_id": variant_id,
        "evaluation_seed": evaluation_seed,
        "channel_seed": channel_seed,
        "manifest": manifest_path.resolve().as_posix(),
        "manifest_sha256": manifest_digest,
        "source_manifest_sha256": manifest["source_manifest_sha256"],
    }
    _json(output / "aged_mrt_formal_freeze.json", freeze)
    tasks: list[tuple[dict[str, Any], Path, Path]] = []
    for snr in grid:
        point_config = copy.deepcopy(template)
        snr_tag = f"{snr:.1f}".replace(".", "p")
        run_id = f"aged_mrt_60kmh_csi40ms_snr{snr_tag}_t1_1000"
        point_config["simulation"].update({
            "snr_points_db": [snr],
            "n_trials_per_snr": trials,
            "max_trials_per_snr": trials,
            "absolute_trial_start": 0,
            "output_dir": (root / "aged_mrt_60kmh" / "batches").as_posix(),
            "run_id": run_id,
        })
        tasks.append((
            point_config,
            root / "configs" / f"{run_id}.yaml",
            root / "aged_mrt_60kmh" / "batches" / run_id,
        ))
    batch_dirs: list[Path] = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(_run_batch, config, config_path, directory): directory
            for config, config_path, directory in tasks
        }
        for future in as_completed(futures):
            batch_dirs.append(future.result())
    summary = _merge_summary(batch_dirs)
    _write_csv(output / "aged_mrt_formal_summary.csv", summary)
    _merge_trial_files(batch_dirs, output / "aged_mrt_formal_trial_metrics.csv", {variant_id})
    if len(summary) != len(grid) or any(int(row["n_trials"]) != trials for row in summary):
        raise RuntimeError("Plan-039 aged MRT formal output coverage is incomplete.")
    report = {
        "schema": "plan039-aged-mrt-formal-report-v1",
        "status": "COMPLETE",
        "prescan_run": False,
        "adaptive_additions": False,
        "variant_id": variant_id,
        "snr_grid_db": grid,
        "trials_per_snr": trials,
        "freeze_sha256": _sha256(output / "aged_mrt_formal_freeze.json"),
        "summary_sha256": _sha256(output / "aged_mrt_formal_summary.csv"),
        "trial_metrics_sha256": _sha256(output / "aged_mrt_formal_trial_metrics.csv"),
    }
    _json(output / "aged_mrt_formal_report.json", report)
    return report


def _run_estimated_family(
    raw: dict[str, Any],
    root: Path,
    template: dict[str, Any],
    family: str,
    variants: list[dict[str, Any]],
    adaptive: bool = True,
) -> list[Path]:
    grid = _stage1b_grid(raw, family)
    variant_ids = tuple(str(row["variant_id"]) for row in variants)
    family_template = copy.deepcopy(template)
    batch_dirs: list[Path] = []

    def run_variant(variant: dict[str, Any], start: int, points: list[float], count: int) -> None:
        config = copy.deepcopy(family_template)
        config["variants"] = [variant]
        variant_id = str(variant["variant_id"])
        short = variant_id.replace("__plan039_", "__").replace("_common_reference_pdp", "")
        points_digest = hashlib.sha256(
            json.dumps(points, separators=(",", ":")).encode("utf-8")
        ).hexdigest()[:8]
        run_id = f"{family}_{short}_t{start+1}_{start+count}_p{points_digest}"
        config["simulation"].update({"snr_points_db": points, "n_trials_per_snr": count,
                                     "max_trials_per_snr": count, "absolute_trial_start": start,
                                     "output_dir": (root / "estimated_confirm" / "batches").as_posix(),
                                     "run_id": run_id})
        directory = root / "estimated_confirm" / "batches" / run_id
        batch_dirs.append(_run_batch(config, root / "configs" / f"stage1b_{run_id}.yaml", directory))

    legacy_initial = root / "estimated_confirm" / "batches" / f"{family}_batch_001_t1_{FORMAL_BATCH}"
    if _completed_batch_covers(legacy_initial, variant_ids, grid, 0, FORMAL_BATCH):
        batch_dirs.append(legacy_initial)
    else:
        for variant in variants:
            run_variant(variant, 0, grid, FORMAL_BATCH)
    if not adaptive:
        return batch_dirs
    while True:
        summary = _merge_summary(batch_dirs)
        requests = _adaptive_points(summary, variant_ids)
        if not requests:
            break
        for current, points in requests.items():
            count = min(FORMAL_BATCH, FORMAL_CAP - current)
            if count > 0:
                for variant in variants:
                    run_variant(variant, current, points, count)
    return batch_dirs


def _completed_batch_covers(
    directory: Path,
    variants: tuple[str, ...],
    points: list[float],
    absolute_start: int,
    count: int,
) -> bool:
    required = {(variant, float(point)) for variant in variants for point in points}
    summary_path = directory / "summary.csv"
    trial_path = directory / "trial_metrics.csv"
    receipt_path = directory / "batch_receipt.json"
    if not (summary_path.exists() and trial_path.exists() and receipt_path.exists()):
        return False
    rows = _read_csv(summary_path)
    covered = {
        (str(row["variant_id"]), float(row["snr_db"]))
        for row in rows
        if int(row["absolute_trial_start"]) == int(absolute_start)
        and int(row["n_trials"]) == int(count)
    }
    return required <= covered


def _completed_method1_sidon_intervals(root: Path) -> dict[tuple[int, int], list[float]]:
    """Return completed formal Sidon method-1 intervals grouped by absolute trial range."""
    variant_id = "sidon_selected__plan039_common_reference_pdp"
    batches = root / "estimated_confirm" / "batches"
    grouped: dict[tuple[int, int], set[float]] = {}
    if not batches.exists():
        raise FileNotFoundError("Run the Sidon method-1 formal scan before the transparent supplement.")
    for directory in sorted(path for path in batches.iterdir() if path.is_dir()):
        summary_path = directory / "summary.csv"
        receipt_path = directory / "batch_receipt.json"
        trial_path = directory / "trial_metrics.csv"
        if not (summary_path.exists() and receipt_path.exists() and trial_path.exists()):
            continue
        for row in _read_csv(summary_path):
            if str(row["variant_id"]) != variant_id:
                continue
            start = int(row["absolute_trial_start"])
            stop = int(row.get("absolute_trial_stop", start + int(row["n_trials"])))
            if stop - start != int(row["n_trials"]) or start < 0 or stop <= start:
                raise RuntimeError(f"Invalid method-1 interval in {summary_path}: {row}.")
            grouped.setdefault((start, stop), set()).add(float(row["snr_db"]))
    if not grouped:
        raise FileNotFoundError("No completed Sidon method-1 formal batches with receipts were found.")
    return {interval: sorted(points) for interval, points in sorted(grouped.items())}


def run_transparent_sidon_formal(raw: dict[str, Any], root: Path) -> dict[str, Any]:
    """Run transparent estimated-CSI Sidon on exactly the completed method-1 coverage."""
    freeze_path = root / "selection" / SCENARIO / "selection_freeze.json"
    if not freeze_path.exists():
        raise FileNotFoundError("Run --stage select before --stage transparent_sidon_formal.")
    selection = json.loads(freeze_path.read_text(encoding="utf-8"))
    if selection.get("status") != "FROZEN":
        raise RuntimeError("Stage 1A selection freeze is not valid.")
    intervals = _completed_method1_sidon_intervals(root)
    base_manifest, _, base_payload, _ = _stage0_sources(root)
    seeds = selection["seeds"]
    template = _stage_base(raw, base_manifest.resolve(), str(base_payload["manifest_sha256"]),
                           int(seeds["evaluation"]), int(seeds["evaluation_channel"]))
    selected = selection["selected"]
    variant_id = "sidon_selected__plan039_transparent_common_reference_pdp"
    variant = {
        "variant_id": variant_id,
        "transmission": {"tx_scheme": "BEAM8_SIDON_SELECTED"},
        "channel_estimation": {"ce_method": TRANSPARENT_CDD_METHOD},
        "fixed_cdl_statistics": {
            "frozen_beam_manifest": selected["manifest"],
            "frozen_beam_manifest_sha256": selected["manifest_sha256"],
        },
    }
    output = root / "transparent_sidon_formal" / SCENARIO
    coverage = [
        {"absolute_trial_start": start, "absolute_trial_stop": stop,
         "trials": stop - start, "snr_points_db": points}
        for (start, stop), points in intervals.items()
    ]
    scope = {
        "schema": "plan039-transparent-sidon-formal-freeze-v1",
        "status": "FROZEN",
        "selection_freeze_sha256": _sha256(freeze_path),
        "source": "completed Sidon method-1 batches with batch receipts",
        "scheme": "BEAM8_SIDON_SELECTED",
        "method": TRANSPARENT_CDD_METHOD,
        "variant_id": variant_id,
        "coverage": coverage,
        "seeds": seeds,
    }
    _json(output / "transparent_sidon_formal_freeze.json", scope)
    batch_dirs: list[Path] = []
    for (start, stop), points in intervals.items():
        count = stop - start
        config = copy.deepcopy(template)
        config["variants"] = [variant]
        points_digest = hashlib.sha256(
            json.dumps(points, separators=(",", ":")).encode("utf-8")
        ).hexdigest()[:8]
        run_id = f"transparent_sidon_t{start+1}_{stop}_p{points_digest}"
        config["simulation"].update({
            "snr_points_db": points,
            "n_trials_per_snr": count,
            "max_trials_per_snr": count,
            "absolute_trial_start": start,
            "output_dir": (root / "transparent_sidon_formal" / "batches").as_posix(),
            "run_id": run_id,
        })
        directory = root / "transparent_sidon_formal" / "batches" / run_id
        batch_dirs.append(_run_batch(
            config, root / "configs" / f"stage1b_{run_id}.yaml", directory
        ))
    summary = _merge_summary(batch_dirs)
    _write_csv(output / "transparent_sidon_summary.csv", summary)
    _merge_trial_files(batch_dirs, output / "transparent_sidon_trial_metrics.csv", {variant_id})
    expected_trials = {snr: 0 for points in intervals.values() for snr in points}
    for (start, stop), points in intervals.items():
        for snr in points:
            expected_trials[snr] += stop - start
    actual_trials = {float(row["snr_db"]): int(row["n_trials"]) for row in summary}
    if actual_trials != expected_trials:
        raise RuntimeError(
            f"Transparent Sidon coverage differs from method 1: {actual_trials} != {expected_trials}."
        )
    report = {
        "schema": "plan039-transparent-sidon-formal-report-v1",
        "status": "COMPLETE",
        "variant_id": variant_id,
        "coverage": coverage,
        "freeze_sha256": _sha256(output / "transparent_sidon_formal_freeze.json"),
        "summary_sha256": _sha256(output / "transparent_sidon_summary.csv"),
        "trial_metrics_sha256": _sha256(output / "transparent_sidon_trial_metrics.csv"),
    }
    _json(output / "transparent_sidon_formal_report.json", report)
    return report


def _bootstrap_gains(batch_dirs: list[Path], seed: int, repeats: int = 1000) -> dict[str, Any]:
    variants = ("b0__ideal", "sidon_selected__ideal", "cycling__ideal")
    collected: dict[tuple[str, float], list[tuple[int, int]]] = {}
    for directory in batch_dirs:
        for row in _read_csv(directory / "trial_metrics.csv"):
            variant = row["variant_id"]
            if variant in variants:
                key = (variant, float(row["snr_db"]))
                collected.setdefault(key, []).append((int(row["trial"]), int(row["tb_error"])))
    snrs = sorted({snr for _, snr in collected})
    arrays: dict[tuple[str, float], np.ndarray] = {}
    for variant in variants:
        for snr in snrs:
            ordered = sorted(collected.get((variant, snr), []))
            arrays[(variant, snr)] = np.asarray([error for _, error in ordered], dtype=np.int8)
    rng = np.random.default_rng(int(seed))
    gains_b0: list[float] = []
    gains_cycling: list[float] = []
    for _ in range(int(repeats)):
        curves = {variant: [] for variant in variants}
        for snr in snrs:
            count = min(arrays[(variant, snr)].size for variant in variants)
            if count == 0:
                continue
            indices = rng.integers(0, count, size=count)
            for variant in variants:
                curves[variant].append({"snr_db": snr, "bler": float(np.mean(arrays[(variant, snr)][:count][indices]))})
        crossings = {variant: _crossing(curves[variant], 0.10)[0] for variant in variants}
        if all(math.isfinite(value) for value in crossings.values()):
            gains_b0.append(crossings["b0__ideal"] - crossings["sidon_selected__ideal"])
            gains_cycling.append(crossings["cycling__ideal"] - crossings["sidon_selected__ideal"])
    def interval(values: list[float]) -> dict[str, Any]:
        if not values:
            return {"effective_repeats": 0, "estimate_db": float("nan"), "ci95_db": [float("nan"), float("nan")]}
        return {"effective_repeats": len(values), "estimate_db": float(np.mean(values)),
                "ci95_db": [float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))]}
    return {"repeats_requested": int(repeats), "b0_minus_sidon_10pct_gain": interval(gains_b0),
            "cycling_minus_sidon_10pct_gain": interval(gains_cycling)}


def run_estimated_confirm(raw: dict[str, Any], root: Path) -> dict[str, Any]:
    freeze_path = root / "selection" / SCENARIO / "selection_freeze.json"
    if not freeze_path.exists():
        raise FileNotFoundError("Run --stage select before --stage estimated_confirm.")
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    if freeze.get("status") != "FROZEN":
        raise RuntimeError("Stage 1A selection freeze is not valid.")
    base_manifest, _, base_payload, _ = _stage0_sources(root)
    seeds = freeze["seeds"]
    template = _stage_base(raw, base_manifest.resolve(), str(base_payload["manifest_sha256"]),
                           int(seeds["evaluation"]), int(seeds["evaluation_channel"]))
    all_variants = _evaluation_variants(freeze)
    cdd_variants = [row for row in all_variants if row["transmission"]["tx_scheme"] != "BEAM8_PRECODER_CYCLING"]
    cycling_variants = [row for row in all_variants if row["transmission"]["tx_scheme"] == "BEAM8_PRECODER_CYCLING"]
    output = root / "estimated_confirm" / SCENARIO
    grids = {family: _stage1b_grid(raw, family) for family in ("cdd", "cycling")}
    grid_digest = hashlib.sha256(
        json.dumps(grids, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    stage_freeze = {
        "schema": "plan039-stage1b-estimated-freeze-v1",
        "status": "FROZEN",
        "selection_freeze_sha256": _sha256(freeze_path),
        "ideal_csi_run": False,
        "methods": {"cdd": list(STAGE1B_CDD_METHODS), "cycling": STAGE1B_CYCLING_METHOD},
        "snr_grids_db": grids,
        "snr_grid_sha256": grid_digest,
        "budget": {"initial_trials_per_snr": FORMAL_BATCH,
                   "batch_trials": FORMAL_BATCH, "max_trials_per_snr": FORMAL_CAP},
        "variant_ids": [row["variant_id"] for row in all_variants],
        "seeds": seeds,
    }
    _json(output / "stage1b_method1_only_freeze.json", stage_freeze)
    batch_dirs = _run_estimated_family(raw, root, template, "cdd", cdd_variants)
    batch_dirs.extend(_run_estimated_family(raw, root, template, "cycling", cycling_variants))
    summary = _merge_summary(batch_dirs)
    active_variant_ids = {str(row["variant_id"]) for row in all_variants}
    summary = [row for row in summary if str(row["variant_id"]) in active_variant_ids]
    _write_csv(output / "method1_summary.csv", summary)
    _merge_trial_files(batch_dirs, output / "method1_trial_metrics.csv", active_variant_ids)
    report = {"schema": "plan039-stage1b-estimated-v1", "status": "COMPLETE",
              "ideal_csi_run": False,
              "methods": {"cdd": list(STAGE1B_CDD_METHODS), "cycling": STAGE1B_CYCLING_METHOD},
              "snr_grids_db": grids,
              "snr_grid_sha256": grid_digest,
              "variant_ids": [row["variant_id"] for row in all_variants],
              "stage1b_freeze_sha256": _sha256(output / "stage1b_method1_only_freeze.json"),
              "summary_sha256": _sha256(output / "method1_summary.csv"),
              "trial_metrics_sha256": _sha256(output / "method1_trial_metrics.csv")}
    _json(output / "stage1b_method1_only_report.json", report)
    return report


def run_cycling_initial(raw: dict[str, Any], root: Path) -> dict[str, Any]:
    freeze_path = root / "selection" / SCENARIO / "selection_freeze.json"
    if not freeze_path.exists():
        raise FileNotFoundError("Run --stage select before --stage cycling_initial.")
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    if freeze.get("status") != "FROZEN":
        raise RuntimeError("Stage 1A selection freeze is not valid.")
    base_manifest, _, base_payload, _ = _stage0_sources(root)
    seeds = freeze["seeds"]
    template = _stage_base(raw, base_manifest.resolve(), str(base_payload["manifest_sha256"]),
                           int(seeds["evaluation"]), int(seeds["evaluation_channel"]))
    cycling_variants = [
        row for row in _evaluation_variants(freeze)
        if row["transmission"]["tx_scheme"] == "BEAM8_PRECODER_CYCLING"
    ]
    output = root / "estimated_confirm" / SCENARIO
    grid = _stage1b_grid(raw, "cycling")
    scope = {
        "schema": "plan039-stage1b-cycling-initial-freeze-v1",
        "status": "FROZEN",
        "selection_freeze_sha256": _sha256(freeze_path),
        "adaptive_additions": False,
        "method": STAGE1B_CYCLING_METHOD,
        "snr_grid_db": grid,
        "snr_grid_sha256": hashlib.sha256(
            json.dumps(grid, separators=(",", ":")).encode("utf-8")
        ).hexdigest(),
        "initial_trials_per_snr": FORMAL_BATCH,
        "variant_ids": [row["variant_id"] for row in cycling_variants],
        "seeds": seeds,
    }
    _json(output / "stage1b_cycling_initial_freeze.json", scope)
    batch_dirs = _run_estimated_family(
        raw, root, template, "cycling", cycling_variants, adaptive=False
    )
    active = {str(row["variant_id"]) for row in cycling_variants}
    summary = [
        row for row in _merge_summary(batch_dirs)
        if str(row["variant_id"]) in active
    ]
    _write_csv(output / "cycling_initial_summary.csv", summary)
    _merge_trial_files(batch_dirs, output / "cycling_initial_trial_metrics.csv", active)
    report = {
        "schema": "plan039-stage1b-cycling-initial-v1",
        "status": "COMPLETE",
        "adaptive_additions": False,
        "snr_grid_db": grid,
        "scope_freeze_sha256": _sha256(output / "stage1b_cycling_initial_freeze.json"),
        "summary_sha256": _sha256(output / "cycling_initial_summary.csv"),
        "trial_metrics_sha256": _sha256(output / "cycling_initial_trial_metrics.csv"),
    }
    _json(output / "stage1b_cycling_initial_report.json", report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--stage", choices=(
        "validate", "smoke", "audit", "select", "estimated_confirm", "cycling_initial",
        "aged_mrt_formal", "transparent_sidon_formal",
    ), required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    raw = yaml.safe_load(args.config.read_text(encoding="utf-8")) or {}
    root = (args.output or ROOT/"outputs"/"experiment039_cdl_beam_bler"/str(raw["plan039"]["run_id"])).resolve()
    smoke_config = root/"configs"/f"{SCENARIO}_stage0_smoke.yaml"
    if args.stage == "validate":
        smoke_config = prepare_stage0(raw, args.config, root)
        print(json.dumps({"status": "PASS", "stage": "validate", "smoke_config": smoke_config.as_posix()}, ensure_ascii=False, indent=2))
    elif args.stage == "smoke":
        if not smoke_config.exists():
            raise FileNotFoundError("Run --stage validate before --stage smoke.")
        CDDLinkLevelOrchestrator(load_config(smoke_config)).run()
        audit_smoke(root)
        print(json.dumps({"status": "PASS", "stage": "smoke", "no_performance_conclusion": True}, indent=2))
    elif args.stage == "audit":
        report = audit_smoke(root)
        print(json.dumps(report, indent=2))
    elif args.stage == "select":
        report = run_selection(raw, root)
        print(json.dumps({"status": report["status"], "stage": "select",
                          "selected": report["selected"],
                          "evaluation_snr_grid_db": report["evaluation_snr_grid_db"]},
                         ensure_ascii=False, indent=2))
    elif args.stage == "estimated_confirm":
        report = run_estimated_confirm(raw, root)
        print(json.dumps(report, ensure_ascii=False, indent=2))
    elif args.stage == "cycling_initial":
        report = run_cycling_initial(raw, root)
        print(json.dumps(report, ensure_ascii=False, indent=2))
    elif args.stage == "aged_mrt_formal":
        report = run_aged_mrt_formal(raw, root)
        print(json.dumps(report, ensure_ascii=False, indent=2))
    else:
        report = run_transparent_sidon_formal(raw, root)
        print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
