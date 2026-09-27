"""Build and visualize the Plan-037 section-13 CDL beam platform artifacts."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import dataclass_to_dict, config_from_dict
from cdd_lls.phy.cdl_beam_platform import (
    ALGORITHM_VERSION,
    AngularRegion,
    CoverageRectangle,
    array_sha256,
    beam_domain_cdd_precoder,
    beam_path_statistics,
    covariance_audit,
    covariance_beam_joint,
    covariance_beam_specific_independent,
    covariance_common_reference_pdp,
    direction_cosines,
    dft_steering_beams,
    maximum_normalized_pattern_capture,
    partition_angular_region,
    regular_dft_codebook,
    rows_to_csv,
    select_beams,
    sha256_json,
    synthesize_region_beam,
    txru_response_matrix,
    validate_regular_txru_mapping,
    validate_cache,
    write_cache,
)
from cdd_lls.phy.channel_cdl_fixed import FixedCDLStatisticsChannel
from cdd_lls.phy.resource_grid import build_resource_grid


def _json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _load(path: Path) -> tuple[dict[str, Any], Any]:
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    platform = {key: value for key, value in raw.items() if key not in {"beam_design", "visualization", "platform_run"}}
    return raw, config_from_dict(platform)


def _profile_config(base: Any, profile: str) -> Any:
    payload = dataclass_to_dict(base)
    payload["channel"]["cdl_profile"] = str(profile).upper()
    payload["fixed_cdl_statistics"]["selection_only"] = True
    payload["transmission"]["tx_scheme"] = "BEAM8_B0_QC"
    payload["channel_estimation"]["ce_method"] = "IDEAL"
    return config_from_dict(payload)


def _ray_table(channel: FixedCDLStatisticsChannel) -> dict[str, np.ndarray]:
    arrays = channel.profile_arrays()
    aod = np.asarray(arrays["aod_rad"], dtype=np.float64).reshape(-1)
    zod = np.asarray(arrays["zod_rad"], dtype=np.float64).reshape(-1)
    clusters = np.asarray(arrays["powers"], dtype=np.float64).reshape(-1)
    cluster_delays = np.asarray(arrays["delays_s"], dtype=np.float64).reshape(-1)
    rays_per_cluster = aod.size // clusters.size
    los = bool(np.asarray(arrays["los_indicator"]).reshape(-1)[0])
    k_factor = float(np.asarray(arrays["k_factor"]).reshape(-1)[0])
    diffuse_scale = 1.0 / (k_factor + 1.0) if los else 1.0
    power = np.repeat(clusters * diffuse_scale / rays_per_cluster, rays_per_cluster)
    delay = np.repeat(cluster_delays, rays_per_cluster)
    kind = np.asarray(["diffuse"] * power.size, dtype="U16")
    if los:
        aod = np.r_[aod, float(np.asarray(arrays["los_aod_rad"]).reshape(-1)[0])]
        zod = np.r_[zod, float(np.asarray(arrays["los_zod_rad"]).reshape(-1)[0])]
        delay = np.r_[delay, cluster_delays[0]]
        power = np.r_[power, k_factor / (k_factor + 1.0)]
        kind = np.concatenate((kind, np.asarray(["los_specular"], dtype="U16")))
    uh, uv = direction_cosines(aod, zod)
    power /= np.sum(power)
    return {"aod_rad": aod, "zod_rad": zod, "delay_s": delay, "power": power, "kind": kind,
            "u_horizontal": uh, "u_vertical": uv}


def _rectangle_dict(rectangle: CoverageRectangle) -> dict[str, float]:
    return {
        "horizontal_min": rectangle.horizontal_min, "horizontal_max": rectangle.horizontal_max,
        "vertical_min": rectangle.vertical_min, "vertical_max": rectangle.vertical_max,
        "area": rectangle.area, "covered_power": rectangle.covered_power, "total_power": rectangle.total_power,
        "coverage": rectangle.coverage,
    }


def _angular_region_dict(region: AngularRegion) -> dict[str, float]:
    return {
        "aod_min_deg": region.aod_min_deg,
        "aod_max_deg": region.aod_max_deg,
        "zod_min_deg": region.zod_min_deg,
        "zod_max_deg": region.zod_max_deg,
        "center_aod_deg": region.center_aod_deg,
        "center_zod_deg": region.center_zod_deg,
    }


def _merge_pdp(delays: np.ndarray, powers: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(delays, kind="stable")
    merged_d: list[float] = []
    merged_p: list[float] = []
    for delay, power in zip(delays[order], powers[order]):
        if merged_d and abs(float(delay) - merged_d[-1]) <= 1e-15:
            merged_p[-1] += float(power)
        else:
            merged_d.append(float(delay)); merged_p.append(float(power))
    return np.asarray(merged_d), np.asarray(merged_p)


def _plot_patterns(output: Path, pattern_file: Path, selected: np.ndarray,
                   ssb_regions: list[AngularRegion], narrow_regions: list[AngularRegion]) -> None:
    with np.load(pattern_file) as data:
        regular = np.asarray(data["regular_selected_power"])
        ssb = np.asarray(data["ssb_power"]); split = np.asarray(data["split_power"])
        selected_ssb_index = int(np.asarray(data["selected_ssb_index"]).reshape(-1)[0])
        pattern_aod_deg = np.asarray(data["pattern_aod_deg"])
        pattern_zod_deg = np.asarray(data["pattern_zod_deg"])
        horizontal_angle_deg = np.asarray(data["horizontal_cut_aod_deg"])
        horizontal_fixed_zod_deg = float(np.asarray(data["horizontal_cut_fixed_zod_deg"]).reshape(-1)[0])
        regular_horizontal = np.asarray(data["regular_horizontal_cut_power"])
        ssb_horizontal = np.asarray(data["ssb_horizontal_cut_power"])
        split_horizontal = np.asarray(data["split_horizontal_cut_power"])
        vertical_angle_deg = np.asarray(data["vertical_cut_zod_deg"])
        vertical_fixed_aod_deg = float(np.asarray(data["vertical_cut_fixed_aod_deg"]).reshape(-1)[0])
        regular_vertical = np.asarray(data["regular_vertical_cut_power"])
        ssb_vertical = np.asarray(data["ssb_vertical_cut_power"])
        split_vertical = np.asarray(data["split_vertical_cut_power"])
    extent = [float(pattern_aod_deg[0]), float(pattern_aod_deg[-1]),
              float(pattern_zod_deg[0]), float(pattern_zod_deg[-1])]
    fig, axes = plt.subplots(1, 3, figsize=(20, 6), constrained_layout=True)
    for ax, values, title in (
        (axes[0], regular, "Regular DFT selected beams"),
        (axes[1], ssb, "8-beam SSB codebook"),
        (axes[2], split, f"8 narrow DFT beams inside SSB {selected_ssb_index}"),
    ):
        maximum = np.max(values, axis=0)
        db = 10 * np.log10(np.maximum(maximum / np.max(maximum), 1e-6))
        image = ax.imshow(db, origin="lower", extent=extent, aspect="auto", vmin=-30, vmax=0, cmap="viridis")
        ax.set(xlabel="AoD (deg)", ylabel="ZoD (deg)", title=title)
        fig.colorbar(image, ax=ax, label="Max normalized response (dB)")
    for region in ssb_regions:
        axes[1].add_patch(plt.Rectangle((region.aod_min_deg, region.zod_min_deg),
                                        region.aod_max_deg-region.aod_min_deg,
                                        region.zod_max_deg-region.zod_min_deg,
                                        fill=False, edgecolor="white", linewidth=0.8))
    for region in narrow_regions:
        axes[2].add_patch(plt.Rectangle((region.aod_min_deg, region.zod_min_deg),
                                        region.aod_max_deg-region.aod_min_deg,
                                        region.zod_max_deg-region.zod_min_deg,
                                        fill=False, edgecolor="white", linewidth=0.8))
    fig.savefig(output / "codebook_direction_pattern_2d.png", dpi=180); plt.close(fig)

    fig = plt.figure(figsize=(16, 14))
    axes = [fig.add_subplot(2, 2, i + 1, projection="polar") for i in range(4)]

    def plot_regular(ax: Any, angle_deg: np.ndarray, cuts: np.ndarray) -> None:
        for beam, cut in zip(selected, cuts):
            ax.plot(np.deg2rad(angle_deg), 10*np.log10(np.maximum(cut/np.max(cut), 1e-5)), label=f"b{beam}")

    def plot_wide_split(ax: Any, angle_deg: np.ndarray, ssb_cuts: np.ndarray, split_cuts: np.ndarray) -> None:
        wide_cut = ssb_cuts[selected_ssb_index]
        ax.plot(np.deg2rad(angle_deg), 10*np.log10(np.maximum(wide_cut/np.max(wide_cut), 1e-5)),
                color="black", linewidth=3, label=f"selected SSB {selected_ssb_index}")
        for index, cut in enumerate(split_cuts):
            ax.plot(np.deg2rad(angle_deg), 10*np.log10(np.maximum(cut/np.max(cut), 1e-5)), label=f"s{index}")

    plot_regular(axes[0], horizontal_angle_deg, regular_horizontal)
    plot_wide_split(axes[1], horizontal_angle_deg, ssb_horizontal, split_horizontal)
    plot_regular(axes[2], vertical_angle_deg, regular_vertical)
    plot_wide_split(axes[3], vertical_angle_deg, ssb_vertical, split_vertical)
    titles = (
        f"Regular DFT horizontal cut (fixed ZoD={horizontal_fixed_zod_deg:.2f} deg)",
        f"Selected SSB/narrow horizontal cut (fixed ZoD={horizontal_fixed_zod_deg:.2f} deg)",
        f"Regular DFT vertical cut (fixed AoD={vertical_fixed_aod_deg:.2f} deg)",
        f"Selected SSB/narrow vertical cut (fixed AoD={vertical_fixed_aod_deg:.2f} deg)",
    )
    for index, (ax, title) in enumerate(zip(axes, titles)):
        ax.set_thetamin(-90); ax.set_thetamax(90); ax.set_rlim(-40, 0); ax.set_title(title, pad=20)
        ax.set_xlabel("AoD" if index < 2 else "ZoD")
        ax.legend(fontsize=7, ncol=2, loc="lower center")
    fig.tight_layout(); fig.savefig(output / "codebook_direction_pattern_polar.png", dpi=180); plt.close(fig)


def _plot_power_series(output: Path, profile: str, rays: dict[str, np.ndarray], regular_rows: list[dict[str, Any]],
                       ssb_rows: list[dict[str, Any]], wide_rows: list[dict[str, Any]],
                       pdp_rows: list[dict[str, Any]]) -> None:
    selected = [row for row in regular_rows if row["selected"]]
    fig, axes = plt.subplots(1, 3, figsize=(21, 5), constrained_layout=True)
    axes[0].bar([row["beam_index"] for row in regular_rows], [row["normalized_db"] for row in regular_rows])
    axes[0].scatter([row["beam_index"] for row in selected], [row["normalized_db"] for row in selected], color="red", marker="*")
    axes[0].set(xlabel="Adaptive DFT beam index", ylabel="Normalized analytic power (dB)", title=f"CDL-{profile}: regular DFT power")
    axes[0].grid(True, alpha=.3)
    ssb_max = max(row["analytic_receive_power"] for row in ssb_rows)
    ssb_db = [10*np.log10(max(row["analytic_receive_power"]/ssb_max, 1e-15)) for row in ssb_rows]
    colors = ["red" if row["selected"] else "steelblue" for row in ssb_rows]
    axes[1].bar([row["beam"] for row in ssb_rows], ssb_db, color=colors)
    axes[1].set(xlabel="SSB beam", ylabel="Normalized analytic receive power (dB)",
                title=f"CDL-{profile}: strongest SSB selection")
    axes[1].grid(True, alpha=.3)
    axes[2].bar([row["beam"] for row in wide_rows], [row["power_relative_to_trace"] for row in wide_rows])
    axes[2].set(xlabel="Selected SSB / narrow beam", ylabel="Power / covariance trace",
                title=f"CDL-{profile}: second-stage narrow beams")
    axes[2].tick_params(axis="x", rotation=45); axes[2].grid(True, alpha=.3)
    fig.savefig(output / "beam_power_distributions.png", dpi=180); plt.close(fig)

    size = 15 + 180 * np.sqrt(rays["power"] / np.max(rays["power"]))
    colors = np.where(rays["kind"] == "los_specular", "red", "royalblue")
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), constrained_layout=True,
                             subplot_kw={"projection": None})
    axes[0].remove(); axes[0] = fig.add_subplot(1, 2, 1, projection="polar")
    axes[0].scatter(rays["aod_rad"], rays["power"], s=size, c=colors, alpha=.65)
    axes[0].set_title(f"CDL-{profile}: raw ray AoD-power")
    scatter = axes[1].scatter(np.rad2deg(rays["aod_rad"]), np.rad2deg(rays["zod_rad"]), s=size,
                              c=rays["power"], cmap="viridis", alpha=.75)
    axes[1].set(xlabel="AoD (deg)", ylabel="ZoD (deg)", title=f"CDL-{profile}: raw ray AoD-ZoD power")
    fig.colorbar(scatter, ax=axes[1], label="Linear raw ray power")
    fig.savefig(output / "raw_ray_angle_power.png", dpi=180); plt.close(fig)

    groups: dict[str, list[dict[str, Any]]] = {}
    for row in pdp_rows:
        groups.setdefault(str(row["beam"]), []).append(row)
    global_zero_s = min(float(row["delay_s"]) for row in pdp_rows)
    common_max_delay_ns = max((float(row["delay_s"]) - global_zero_s) * 1e9 for row in pdp_rows)

    def stem_group(names: list[str], path: Path, title: str, columns: int, power_field: str,
                   power_scale: float, common_y_max: float, ylabel: str) -> None:
        rows_count = int(np.ceil(len(names) / columns))
        fig, axes = plt.subplots(rows_count, columns, figsize=(4.2 * columns, 3.0 * rows_count),
                                 sharex=True, sharey=True, squeeze=False, constrained_layout=True)
        for ax, name in zip(axes.ravel(), names):
            entries = groups[name]
            x = (np.asarray([row["delay_s"] for row in entries]) - global_zero_s) * 1e9
            y = np.asarray([row[power_field] for row in entries]) / power_scale
            markerline, stemlines, baseline = ax.stem(x, y, basefmt=" ")
            markerline.set_markersize(4); stemlines.set_linewidth(1.1)
            ax.set(title=name, xlim=(-0.02 * max(common_max_delay_ns, 1.0), 1.02 * max(common_max_delay_ns, 1.0)),
                   ylim=(0.0, 1.05 * max(common_y_max, np.finfo(float).eps)))
            ax.grid(True, alpha=.25)
        for ax in axes.ravel()[len(names):]:
            ax.set_visible(False)
        fig.suptitle(title)
        fig.supxlabel("Delay from common global zero (ns)")
        fig.supylabel(ylabel)
        fig.savefig(path, dpi=180); plt.close(fig)

    raw_y_max = max(float(row["normalized_power"]) for row in groups["raw"])
    stem_group(["raw"], output / "raw_pdp.png", f"CDL-{profile}: raw PDP", 1,
               "normalized_power", 1.0, raw_y_max, "Normalized tap power")
    regular_names = [name for name in groups if name.startswith("regular_")]
    wide_names = ["wide"] + [name for name in groups if name.startswith("split_")]
    beam_names = regular_names + wide_names
    beam_total_power = {
        name: sum(float(row["absolute_power"]) for row in groups[name])
        for name in beam_names
    }
    strongest_beam = max(beam_names, key=lambda name: beam_total_power[name])
    strongest_beam_power = beam_total_power[strongest_beam]
    if not np.isfinite(strongest_beam_power) or strongest_beam_power <= 0.0:
        raise RuntimeError("Per-beam PDPs have no positive finite total power.")
    beam_y_max = max(
        float(row["absolute_power"]) / strongest_beam_power
        for name in beam_names for row in groups[name]
    )
    beam_ylabel = f"Tap power / total power of strongest beam ({strongest_beam})"
    stem_group(regular_names, output / "regular_codebook_pdp.png",
               f"CDL-{profile}: selected regular-DFT beam PDPs", min(4, len(regular_names)),
               "absolute_power", strongest_beam_power, beam_y_max, beam_ylabel)
    stem_group(wide_names, output / "wide_split_codebook_pdp.png",
               f"CDL-{profile}: wide/split-codebook beam PDPs", min(3, len(wide_names)),
               "absolute_power", strongest_beam_power, beam_y_max, beam_ylabel)
    (output / "raw_and_per_beam_pdp.png").unlink(missing_ok=True)


def run_profile(raw: dict[str, Any], base: Any, profile: str, root: Path) -> dict[str, Any]:
    config = _profile_config(base, profile)
    grid = build_resource_grid(config.resource)
    channel = FixedCDLStatisticsChannel(config, grid, prepare_reference=False)
    design = raw.get("beam_design", {})
    branches = int(design.get("num_branches", 8))
    ssb_grid = design.get("ssb_grid", {})
    ssb_aod_range = ssb_grid.get("aod_range_deg", [-60.0, 60.0])
    ssb_zod_range = ssb_grid.get("zod_range_deg", [90.0, 110.0])
    ssb_vertical = int(ssb_grid.get("vertical_beams", 2))
    ssb_horizontal = int(ssb_grid.get("horizontal_beams", 4))
    narrow_grid = design.get("narrow_grid", {})
    narrow_vertical = int(narrow_grid.get("vertical_beams", 2))
    narrow_horizontal = int(narrow_grid.get("horizontal_beams", 4))
    if ssb_vertical * ssb_horizontal != 8:
        raise ValueError("The first-stage SSB grid must contain exactly 8 beams.")
    if branches != 8 or narrow_vertical * narrow_horizontal != branches:
        raise ValueError("The selected SSB region must be split into exactly 8 narrow beams.")
    oversampling = design.get("dft_oversampling", {"vertical": 1, "horizontal": 1})
    if isinstance(oversampling, int):
        ov = oh = int(oversampling)
    else:
        ov, oh = int(oversampling.get("vertical", 1)), int(oversampling.get("horizontal", 1))
    rays = _ray_table(channel)
    ray_delays, _, ray_covariances = channel.analytic_ray_covariances()
    if ray_covariances.shape[0] != rays["power"].size or not np.allclose(ray_delays, rays["delay_s"]):
        raise RuntimeError("Expanded ray table and analytic ray covariance components are not aligned.")
    covariance = np.sum(ray_covariances, axis=0)
    fixed = config.fixed_cdl_statistics
    if any(abs(float(value)) > 1e-12 for value in fixed.bs_orientation_deg):
        raise ValueError("REGULAR_DFT_UNSUPPORTED_GEOMETRY: nonzero BS orientation is not yet supported by the grid synthesizer.")
    validate_regular_txru_mapping(channel.txru_mapping, fixed.bs_vertical_aes, fixed.bs_horizontal_aes,
                                  fixed.bs_vertical_txrus_per_pol, fixed.bs_horizontal_txrus_per_pol,
                                  fixed.bs_polarizations)
    regular, regular_meta = regular_dft_codebook(fixed.bs_vertical_txrus_per_pol, fixed.bs_horizontal_txrus_per_pol,
                                                  fixed.bs_polarizations, ov, oh)
    regular_power, selected = select_beams(covariance, regular, branches, design.get("explicit_beam_indices"))
    ssb_regions = partition_angular_region(ssb_aod_range, ssb_zod_range, ssb_vertical, ssb_horizontal)
    grid_size = int(design.get("synthesis_grid_size", 81))
    synthesis_aod = np.linspace(-90.0, 90.0, grid_size)
    synthesis_zod = np.linspace(70.0, 130.0, grid_size)
    mesh_aod, mesh_zod = np.meshgrid(synthesis_aod, synthesis_zod, indexing="xy")
    synthesis_uh, synthesis_uv = direction_cosines(np.deg2rad(mesh_aod.ravel()), np.deg2rad(mesh_zod.ravel()))
    response = txru_response_matrix(synthesis_uh, synthesis_uv, fixed.bs_vertical_aes, fixed.bs_horizontal_aes,
                                    fixed.bs_vertical_spacing_lambda, fixed.bs_horizontal_spacing_lambda, channel.txru_mapping,
                                    fixed.bs_antenna_pattern)
    def synthesis(region: AngularRegion) -> np.ndarray:
        rectangle = CoverageRectangle(region.aod_min_deg, region.aod_max_deg,
                                      region.zod_min_deg, region.zod_max_deg, 0.0, 1.0)
        return synthesize_region_beam(rectangle, response, mesh_aod.ravel(), mesh_zod.ravel(),
                                      float(design.get("regularization", 1e-3)),
                                      design.get("weight_constraint", "unit_norm"))

    ssb = np.column_stack([synthesis(region) for region in ssb_regions])
    ssb_power, strongest_ssb_array = select_beams(covariance, ssb, 1)
    strongest_ssb = int(strongest_ssb_array[0])
    selected_ssb_region = ssb_regions[strongest_ssb]
    regions = partition_angular_region(
        [selected_ssb_region.aod_min_deg, selected_ssb_region.aod_max_deg],
        [selected_ssb_region.zod_min_deg, selected_ssb_region.zod_max_deg],
        narrow_vertical,
        narrow_horizontal,
    )
    center_aod = np.deg2rad([region.center_aod_deg for region in regions])
    center_zod = np.deg2rad([region.center_zod_deg for region in regions])
    center_uh, center_uv = direction_cosines(center_aod, center_zod)
    center_response = txru_response_matrix(
        center_uh, center_uv, fixed.bs_vertical_aes, fixed.bs_horizontal_aes,
        fixed.bs_vertical_spacing_lambda, fixed.bs_horizontal_spacing_lambda,
        channel.txru_mapping, fixed.bs_antenna_pattern,
    )
    split = dft_steering_beams(regions, center_response)
    wide = ssb[:, strongest_ssb]
    selected_regular = regular[:, selected]
    ray_response = txru_response_matrix(
        rays["u_horizontal"], rays["u_vertical"], fixed.bs_vertical_aes, fixed.bs_horizontal_aes,
        fixed.bs_vertical_spacing_lambda, fixed.bs_horizontal_spacing_lambda, channel.txru_mapping,
        fixed.bs_antenna_pattern,
    )
    wide_pattern_capture, wide_pattern_peak = maximum_normalized_pattern_capture(
        wide, ray_response, rays["power"], response
    )
    split_pattern_metrics = [
        maximum_normalized_pattern_capture(split[:, index], ray_response, rays["power"], response)
        for index in range(branches)
    ]

    profile_arrays = channel.profile_arrays()
    profile_hashes = {name: array_sha256(value) for name, value in profile_arrays.items()}
    key_payload = {
        "algorithm_version": ALGORITHM_VERSION, "profile": profile, "delay_spread_ns": config.channel.delay_spread_ns,
        "carrier_frequency_hz": config.channel.carrier_frequency_hz, "direction": config.channel.cdl_direction,
        "fixed_cdl_statistics": dataclass_to_dict(config)["fixed_cdl_statistics"], "beam_design": design,
        "profile_array_sha256": profile_hashes, "mapping_sha256": array_sha256(channel.txru_mapping),
    }
    cache_key = sha256_json(key_payload)
    cache_root = ROOT / str(design.get("cache_dir", "outputs/experiment037_cdl_platform/cache")) / profile.lower() / cache_key
    arrays = {"regular_weights": regular, "selected_regular_weights": selected_regular, "ssb_weights": ssb,
              "wide_weight": wide, "split_weights": split}
    cache_hit = validate_cache(cache_root, cache_key)
    if cache_hit is None:
        cache_manifest = {"schema": "plan037-section13-beam-cache-v2", "cache_key": cache_key, "key_inputs": key_payload,
                          "selected_regular_indices": selected.tolist(), "selected_ssb_index": strongest_ssb,
                          "ssb_regions": [_angular_region_dict(region) for region in ssb_regions],
                          "narrow_regions": [_angular_region_dict(region) for region in regions]}
        cache_manifest_sha256 = write_cache(cache_root, cache_manifest, arrays)
        cache_status = "MISS_WRITTEN"
    else:
        cache_manifest, cached = cache_hit
        arrays = cached
        regular, selected_regular, ssb, wide, split = (
            arrays["regular_weights"], arrays["selected_regular_weights"], arrays["ssb_weights"],
            arrays["wide_weight"], arrays["split_weights"])
        cache_manifest_sha256 = str(cache_manifest["manifest_sha256"])
        cache_status = "HIT"

    all_weights = np.column_stack((selected_regular, wide, split))
    all_powers, all_joint = beam_path_statistics(all_weights, ray_covariances)
    regular_path_power = all_powers[:branches]
    wide_path_power = all_powers[branches]
    split_path_power = all_powers[branches + 1:]
    trace_power = float(np.trace(covariance).real)
    wide_ratio = float(np.sum(wide_path_power) / trace_power)
    output = root / f"cdl_{profile.lower()}"
    output.mkdir(parents=True, exist_ok=True)
    (output / "resolved_config.yaml").write_text(yaml.safe_dump({**dataclass_to_dict(config), "beam_design": design,
                                                                  "visualization": raw.get("visualization", {})}, sort_keys=False),
                                                         encoding="utf-8")
    np.savez_compressed(output / "profile_arrays.npz", **profile_arrays)
    raw_delays = np.asarray(profile_arrays["profile_raw_delays"], dtype=float).reshape(-1)
    raw_powers = np.asarray(profile_arrays["profile_raw_powers"], dtype=float).reshape(-1)
    rows_to_csv(output / "standard_profile_table.csv", [
        {"cluster_index": index, "normalized_delay": raw_delays[index], "linear_power": raw_powers[index]}
        for index in range(raw_delays.size)])
    rows_to_csv(output / "ray_table.csv", [
        {"ray_index": i, "kind": rays["kind"][i], "delay_s": rays["delay_s"][i], "aod_deg": np.rad2deg(rays["aod_rad"][i]),
         "zod_deg": np.rad2deg(rays["zod_rad"][i]), "u_horizontal": rays["u_horizontal"][i],
         "u_vertical": rays["u_vertical"][i], "raw_linear_power": rays["power"][i]}
        for i in range(rays["power"].size)])
    maximum_regular = float(np.max(regular_power))
    regular_rows = [{**regular_meta[i], "analytic_power": float(regular_power[i]),
                     "normalized_db": float(10*np.log10(max(regular_power[i]/maximum_regular, 1e-15))),
                     "selected": int(i in set(selected.tolist()))} for i in range(regular.shape[1])]
    rows_to_csv(output / "regular_dft_power.csv", regular_rows)
    ssb_rows = [{"beam": f"ssb_{index}", **_angular_region_dict(region),
                 "analytic_receive_power": float(ssb_power[index]),
                 "selected": int(index == strongest_ssb)}
                for index, region in enumerate(ssb_regions)]
    rows_to_csv(output / "ssb_codebook_power.csv", ssb_rows)
    wide_rows = [{"beam": f"ssb_{strongest_ssb}", "absolute_power": float(np.sum(wide_path_power)),
                  "power_relative_to_trace": wide_ratio,
                  "maximum_normalized_pattern_capture": wide_pattern_capture,
                  "pattern_peak_linear": wide_pattern_peak,
                  **_angular_region_dict(selected_ssb_region), "selection_status": "STRONGEST_SSB"}]
    for index, region in enumerate(regions):
        wide_rows.append({"beam": f"split_{index}", "absolute_power": float(np.sum(split_path_power[index])),
                          "power_relative_to_trace": float(np.sum(split_path_power[index])/trace_power),
                          "maximum_normalized_pattern_capture": split_pattern_metrics[index][0],
                          "pattern_peak_linear": split_pattern_metrics[index][1],
                          **_angular_region_dict(region), "selection_status": "NARROW_DFT"})
    rows_to_csv(output / "wide_split_power.csv", wide_rows)

    raw_d, raw_p = _merge_pdp(rays["delay_s"], rays["power"])
    pdp_rows: list[dict[str, Any]] = [{"beam": "raw", "delay_s": d, "absolute_power": p,
                                      "normalized_power": p/np.sum(raw_p)} for d, p in zip(raw_d, raw_p)]
    names = [f"regular_{int(value)}" for value in selected] + ["wide"] + [f"split_{i}" for i in range(branches)]
    for name, powers in zip(names, all_powers):
        delays, mass = _merge_pdp(ray_delays, powers)
        for delay, power in zip(delays, mass):
            pdp_rows.append({"beam": name, "delay_s": delay, "absolute_power": power,
                             "normalized_power": power/np.sum(mass)})
    rows_to_csv(output / "pdp_tables.csv", pdp_rows)

    plot_size = int(raw.get("visualization", {}).get("pattern_grid_size", 181))
    pattern_aod_deg = np.linspace(-90.0, 90.0, plot_size)
    pattern_zod_deg = np.linspace(70.0, 130.0, plot_size)
    pattern_aod, pattern_zod = np.meshgrid(pattern_aod_deg, pattern_zod_deg, indexing="xy")
    pattern_uh, pattern_uv = direction_cosines(np.deg2rad(pattern_aod.ravel()), np.deg2rad(pattern_zod.ravel()))
    pattern_response = txru_response_matrix(pattern_uh, pattern_uv, fixed.bs_vertical_aes, fixed.bs_horizontal_aes,
                                            fixed.bs_vertical_spacing_lambda, fixed.bs_horizontal_spacing_lambda,
                                            channel.txru_mapping, fixed.bs_antenna_pattern)
    pattern = lambda weight: np.abs(pattern_response @ weight) ** 2
    fixed_zod_rad = np.deg2rad(selected_ssb_region.center_zod_deg)
    fixed_aod_rad = np.deg2rad(selected_ssb_region.center_aod_deg)
    cut_angle_deg = np.linspace(-90.0, 90.0, plot_size)
    horizontal_aod_rad = np.deg2rad(cut_angle_deg)
    horizontal_h = np.sin(fixed_zod_rad) * np.sin(horizontal_aod_rad)
    horizontal_v = np.full(plot_size, np.cos(fixed_zod_rad))
    vertical_zod_deg = np.linspace(70.0, 130.0, plot_size)
    vertical_zod_rad = np.deg2rad(vertical_zod_deg)
    vertical_h = np.sin(vertical_zod_rad) * np.sin(fixed_aod_rad)
    vertical_v = np.cos(vertical_zod_rad)
    horizontal_response = txru_response_matrix(
        horizontal_h, horizontal_v, fixed.bs_vertical_aes, fixed.bs_horizontal_aes,
        fixed.bs_vertical_spacing_lambda, fixed.bs_horizontal_spacing_lambda, channel.txru_mapping,
        fixed.bs_antenna_pattern,
    )
    vertical_response = txru_response_matrix(
        vertical_h, vertical_v, fixed.bs_vertical_aes, fixed.bs_horizontal_aes,
        fixed.bs_vertical_spacing_lambda, fixed.bs_horizontal_spacing_lambda, channel.txru_mapping,
        fixed.bs_antenna_pattern,
    )
    cut_pattern = lambda response_matrix, weight: np.abs(response_matrix @ weight) ** 2
    pattern_file = output / "direction_pattern_data.npz"
    np.savez_compressed(pattern_file, pattern_aod_deg=pattern_aod_deg, pattern_zod_deg=pattern_zod_deg,
                        selected_ssb_index=np.asarray([strongest_ssb]),
                        regular_selected_power=np.asarray([pattern(selected_regular[:, i]).reshape(plot_size, plot_size)
                                                           for i in range(branches)]),
                        wide_power=pattern(wide).reshape(plot_size, plot_size),
                        ssb_power=np.asarray([pattern(ssb[:, i]).reshape(plot_size, plot_size) for i in range(8)]),
                        split_power=np.asarray([pattern(split[:, i]).reshape(plot_size, plot_size) for i in range(branches)]),
                        horizontal_cut_aod_deg=cut_angle_deg,
                        horizontal_cut_fixed_zod_deg=np.asarray([np.rad2deg(fixed_zod_rad)]),
                        regular_horizontal_cut_power=np.asarray([
                            cut_pattern(horizontal_response, selected_regular[:, i]) for i in range(branches)]),
                        ssb_horizontal_cut_power=np.asarray([
                            cut_pattern(horizontal_response, ssb[:, i]) for i in range(8)]),
                        split_horizontal_cut_power=np.asarray([
                            cut_pattern(horizontal_response, split[:, i]) for i in range(branches)]),
                        vertical_cut_zod_deg=vertical_zod_deg,
                        vertical_cut_fixed_aod_deg=np.asarray([np.rad2deg(fixed_aod_rad)]),
                        regular_vertical_cut_power=np.asarray([
                            cut_pattern(vertical_response, selected_regular[:, i]) for i in range(branches)]),
                        ssb_vertical_cut_power=np.asarray([
                            cut_pattern(vertical_response, ssb[:, i]) for i in range(8)]),
                        split_vertical_cut_power=np.asarray([
                            cut_pattern(vertical_response, split[:, i]) for i in range(branches)]))
    rows_to_csv(cache_root / "ray_table.csv", [
        {"ray_index": i, "kind": rays["kind"][i], "delay_s": rays["delay_s"][i], "aod_deg": np.rad2deg(rays["aod_rad"][i]),
         "zod_deg": np.rad2deg(rays["zod_rad"][i]), "u_horizontal": rays["u_horizontal"][i],
         "u_vertical": rays["u_vertical"][i], "raw_linear_power": rays["power"][i]}
        for i in range(rays["power"].size)])
    rows_to_csv(cache_root / "beam_power.csv", regular_rows + [
        {"beam_index": row["beam"], "q_vertical": "", "q_horizontal": "", "spatial_frequency_vertical_cycles": "",
         "spatial_frequency_horizontal_cycles": "", "analytic_power": row["absolute_power"], "normalized_db": "",
         "selected": "wide_or_split"} for row in wide_rows])
    (cache_root / "direction_pattern_data.npz").write_bytes(pattern_file.read_bytes())
    cache_payload = json.loads((cache_root / "manifest.json").read_text(encoding="utf-8"))
    cache_payload.pop("manifest_sha256", None)
    cache_payload["artifact_sha256"] = {
        name: hashlib.sha256((cache_root / name).read_bytes()).hexdigest()
        for name in ("ray_table.csv", "beam_power.csv", "direction_pattern_data.npz")
    }
    cache_manifest_sha256 = sha256_json(cache_payload)
    cache_payload["manifest_sha256"] = cache_manifest_sha256
    _json(cache_root / "manifest.json", cache_payload)
    _plot_patterns(output, pattern_file, selected, ssb_regions, regions)
    _plot_power_series(output, profile, rays, regular_rows, ssb_rows, wide_rows, pdp_rows)

    delay_indices = np.asarray(design.get("cdd_delay_indices", [0, 1, 3, 7, 12, 20, 30, 65]), dtype=float)
    if delay_indices.size != branches:
        raise ValueError("cdd_delay_indices must match num_branches.")
    full_precoder, full_alpha, orthogonal = beam_domain_cdd_precoder(split, delay_indices, np.arange(grid.n_sc), grid.n_sc)
    rng = np.random.default_rng(37013)
    h = rng.normal(size=(2, grid.n_sc, split.shape[0])) + 1j*rng.normal(size=(2, grid.n_sc, split.shape[0]))
    direct = np.einsum("rkt,kt->rk", h, full_precoder)
    branch_h = np.einsum("rkt,tm->rkm", h, split)
    phase = np.exp(-1j*2*np.pi*np.arange(grid.n_sc)[:, None]*delay_indices[None, :]/grid.n_sc)
    projected = np.einsum("rkm,km->rk", branch_h, full_alpha[:, None]*phase)
    sample_index = np.linspace(0, grid.n_sc-1, min(48, grid.n_sc), dtype=int)
    frequencies = np.asarray(grid.subcarrier_indices[sample_index], dtype=float) * grid.scs_khz * 1e3
    artificial_s = delay_indices / (grid.n_sc * grid.scs_khz * 1e3)
    reference_index = int(np.argmax(np.sum(regular_path_power, axis=1)))
    explicit_reference = design.get("reference_beam_index")
    if explicit_reference is not None:
        matches = np.flatnonzero(selected == int(explicit_reference))
        if matches.size != 1:
            raise ValueError("reference_beam_index must identify one selected regular DFT beam.")
        reference_index = int(matches[0])
    common = covariance_common_reference_pdp(frequencies, ray_delays, regular_path_power[reference_index], artificial_s)
    independent = covariance_beam_specific_independent(frequencies, ray_delays, split_path_power, artificial_s)
    split_joint = all_joint[:, branches+1:, branches+1:]
    joint = covariance_beam_joint(frequencies, ray_delays, split_joint, artificial_s, full_alpha[sample_index])
    pilot_positions = np.arange(0, frequencies.size, 6, dtype=int)
    loading = max(float(np.mean(np.real(np.diag(joint)))) * 1e-3, np.finfo(float).eps)
    lmmse_filter = np.linalg.solve(
        joint[np.ix_(pilot_positions, pilot_positions)] + loading * np.eye(pilot_positions.size),
        joint[pilot_positions, :],
    ).conj().T
    synthetic_ls = np.ones(pilot_positions.size, dtype=np.complex128)
    synthetic_lmmse = lmmse_filter @ synthetic_ls
    if synthetic_lmmse.shape != (frequencies.size,) or not np.all(np.isfinite(synthetic_lmmse)):
        raise RuntimeError("Synthetic frequency-LMMSE shape smoke failed.")
    covariance_report = {"common_reference_pdp": covariance_audit(common),
                         "beam_specific_pdp_independent": covariance_audit(independent),
                         "beam_joint_covariance": covariance_audit(joint),
                         "synthetic_ls_lmmse_shape_smoke": {
                             "status": "PASS", "pilot_count": int(pilot_positions.size),
                             "ls_shape": list(synthetic_ls.shape), "lmmse_filter_shape": list(lmmse_filter.shape),
                             "estimate_shape": list(synthetic_lmmse.shape), "diagonal_loading": loading,
                         },
                         "input_hash": sha256_json({"ray_covariances": array_sha256(ray_covariances),
                                                    "split_weights": array_sha256(split), "delay_indices": delay_indices.tolist()})}
    np.savez_compressed(output / "covariance_audit_data.npz", frequencies_hz=frequencies, common=common,
                        beam_specific_independent=independent, beam_joint=joint)
    _json(output / "covariance_audit.json", covariance_report)
    smoke_count = int(design.get("smoke_realizations", 2))
    smoke_powers = np.asarray([channel.realization_beam_statistics(index, regular)[1] for index in range(smoke_count)])
    smoke_mean = np.mean(smoke_powers, axis=0)
    analytic_normalized = regular_power / np.max(regular_power)
    smoke_normalized = smoke_mean / np.max(smoke_mean)
    smoke_db_error = np.max(np.abs(10*np.log10(np.maximum(analytic_normalized, 1e-15))
                                       - 10*np.log10(np.maximum(smoke_normalized, 1e-15))))
    covariance_direct = channel.analytic_transmit_covariance()
    geometry_audit = {
        "mapping_shape": list(channel.txru_mapping.shape),
        "mapping_column_norm_max_error": float(np.max(np.abs(np.linalg.norm(channel.txru_mapping, axis=0)-1))),
        "covariance_hermitian_max_error": float(np.max(np.abs(covariance-covariance.conj().T))),
        "covariance_minimum_eigenvalue": float(np.min(np.linalg.eigvalsh(covariance))),
        "covariance_trace": trace_power,
        "ray_sum_vs_direct_covariance_max_error": float(np.max(np.abs(covariance-covariance_direct))),
        "smoke_realizations": smoke_count,
        "smoke_normalized_power_max_error_db": float(smoke_db_error),
        "smoke_analytic_top8": np.argsort(-regular_power, kind="stable")[:8].tolist(),
        "smoke_sample_top8": np.argsort(-smoke_mean, kind="stable")[:8].tolist(),
    }
    _json(output / "geometry_profile_audit.json", geometry_audit)
    report = {
        "schema": "plan037-section13-platform-run-v3", "profile": profile,
        "status": "PASS", "algorithm_version": ALGORITHM_VERSION,
        "cache_status": cache_status, "cache_key": cache_key, "cache_manifest_sha256": cache_manifest_sha256,
        "selected_regular_indices": selected.tolist(), "selected_ssb_index": strongest_ssb,
        "ssb_regions": [_angular_region_dict(region) for region in ssb_regions],
        "selected_ssb_region": _angular_region_dict(selected_ssb_region),
        "narrow_regions": [_angular_region_dict(region) for region in regions],
        "selected_ssb_power_relative_to_covariance_trace": wide_ratio,
        "selected_ssb_maximum_normalized_pattern_capture": wide_pattern_capture,
        "selected_ssb_pattern_peak_linear": wide_pattern_peak,
        "regular_codebook_shape": list(regular.shape), "ssb_codebook_shape": list(ssb.shape),
        "split_codebook_shape": list(split.shape),
        "cdd_split_branches_orthogonal": orthogonal,
        "cdd_unit_power_max_error": float(np.max(np.abs(np.sum(np.abs(full_precoder)**2, axis=1)-1))),
        "cdd_direct_projection_max_error": float(np.max(np.abs(direct-projected))), "covariance": covariance_report,
        "reference_regular_beam_index": int(selected[reference_index]), "geometry_profile_audit": geometry_audit,
        "coordinate_semantics": {"aod": "Sionna GCS azimuth", "zod": "Sionna GCS zenith angle",
                                 "u_horizontal": "sin(ZoD)*sin(AoD)", "u_vertical": "cos(ZoD)"},
        "outputs": ["profile_arrays.npz", "standard_profile_table.csv", "ray_table.csv", "regular_dft_power.csv",
                    "ssb_codebook_power.csv", "wide_split_power.csv", "pdp_tables.csv",
                    "direction_pattern_data.npz", "codebook_direction_pattern_2d.png", "codebook_direction_pattern_polar.png",
                    "beam_power_distributions.png", "raw_ray_angle_power.png", "raw_pdp.png",
                    "regular_codebook_pdp.png", "wide_split_codebook_pdp.png",
                    "geometry_profile_audit.json", "covariance_audit.json", "covariance_audit_data.npz"],
    }
    _json(output / "run_report.json", report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--stage", choices=("validate", "run"), required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    raw, base = _load(args.config)
    profiles = [str(value).upper() for value in raw.get("platform_run", {}).get("profiles", [base.channel.cdl_profile])]
    design = raw.get("beam_design", {})
    if str(design.get("method", "regular_dft")) not in {"regular_dft", "wide_beam_split", "both"}:
        raise ValueError("beam_design.method must be regular_dft, wide_beam_split, or both.")
    ssb_grid = design.get("ssb_grid", {})
    narrow_grid = design.get("narrow_grid", {})
    partition_angular_region(ssb_grid.get("aod_range_deg", [-60.0, 60.0]),
                             ssb_grid.get("zod_range_deg", [90.0, 110.0]),
                             int(ssb_grid.get("vertical_beams", 2)), int(ssb_grid.get("horizontal_beams", 4)))
    if int(ssb_grid.get("vertical_beams", 2)) * int(ssb_grid.get("horizontal_beams", 4)) != 8:
        raise ValueError("beam_design.ssb_grid must define 2 x 4 = 8 beams (or another factorization totaling 8).")
    if int(narrow_grid.get("vertical_beams", 2)) * int(narrow_grid.get("horizontal_beams", 4)) != 8:
        raise ValueError("beam_design.narrow_grid must define 2 x 4 = 8 beams (or another factorization totaling 8).")
    if int(design.get("num_branches", 8)) != 8:
        raise ValueError("beam_design.num_branches must be 8 for the two-stage SSB/DFT design.")
    if args.stage == "validate":
        print(json.dumps({"status": "PASS", "profiles": profiles, "num_branches": int(design.get("num_branches", 8)),
                          "methods": ["regular_dft", "two_stage_ssb_dft"], "no_bler": True}, indent=2))
        return
    output = (args.output or (ROOT / raw.get("platform_run", {}).get("output_dir", "outputs/experiment037_cdl_platform/beam_patterns_e_c"))).resolve()
    reports = [run_profile(raw, base, profile, output) for profile in profiles]
    _json(output / "run_summary.json", {"schema": "plan037-section13-platform-summary-v3", "reports": reports,
                                         "reproduce": f"python tools/run_plan037_cdl_platform.py --config {args.config.as_posix()} --stage run --output {output.as_posix()}"})
    print(json.dumps({"output": output.as_posix(), "profiles": profiles,
                      "status": [report["status"] for report in reports]}, indent=2))


if __name__ == "__main__":
    main()
