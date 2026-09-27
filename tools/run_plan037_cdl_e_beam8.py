"""Execute Plan-037 section 13 gates for 32T2R CDL-E Beam8.

Long-term statistics are a one-time calibration.  ``freeze`` writes a
hash-locked manifest which all link stages reuse without reselecting beams.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import dataclass_to_dict, load_config
from cdd_lls.phy.beam8 import (
    CODEBOOK_TYPE,
    DELAYS,
    SCHEMES,
    build_beam8_precoder,
    build_dft_2x8_same_pol_codebook,
    circular_pdp_moments,
    frequency_covariance,
    equivalent_channel_two_paths,
    load_frozen_manifest,
    merge_discrete_pdp,
    pdp_moments,
    shifted_reference_pdp,
    write_manifest,
)
from cdd_lls.phy.channel_cdl_fixed import FixedCDLStatisticsChannel
from cdd_lls.phy.resource_grid import build_resource_grid
from cdd_lls.sim.orchestrator import CDDLinkLevelOrchestrator


SCHEMA = "plan037-cdl-e-32tx-beam8-v1"
SCENARIO = "E100_NT32_NR2_V60_BEAM8"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else [])
        if rows:
            writer.writeheader()
            writer.writerows(rows)


def _context(cfg) -> dict[str, Any]:
    grid = build_resource_grid(cfg.resource)
    return {
        "cdl_profile": str(cfg.channel.cdl_profile).upper(),
        "delay_spread_ns": float(cfg.channel.delay_spread_ns),
        "carrier_frequency_hz": float(cfg.channel.carrier_frequency_hz),
        "ue_speed_kmh": float(cfg.channel.ue_speed_kmh),
        "n_sc": int(grid.n_sc),
        "n_symbols": int(grid.n_symbols),
        "statistics_seed": int(cfg.fixed_cdl_statistics.statistics_seed),
    }


def validate_frozen_scene(cfg, require_manifest: bool) -> None:
    fixed = cfg.fixed_cdl_statistics
    checks = {
        "scenario dimensions": (cfg.antenna.n_tx, cfg.antenna.n_rx) == (32, 2),
        "CDL-E": str(cfg.channel.cdl_profile).upper() == "E",
        "100 ns": float(cfg.channel.delay_spread_ns) == 100.0,
        "4 GHz": float(cfg.channel.carrier_frequency_hz) == 4.0e9,
        "60 km/h": float(cfg.channel.ue_speed_kmh) == 60.0,
        "profile-native angles": bool(fixed.profile_native_angles),
        "DFT2x8 codebook": str(fixed.codebook_type).lower() == CODEBOOK_TYPE,
        "48 PRB": int(cfg.resource.n_prbs) == 48,
        "576 active SC": int(cfg.resource.n_prbs) * 12 == 576,
        "FFT/CP": (int(cfg.resource.n_fft), int(cfg.resource.cyclic_prefix_length)) == (4096, 288),
        "PDSCH/DMRS": (int(cfg.resource.pdsch_n_symbols), list(cfg.resource.dmrs_symbol_indices), int(cfg.resource.dmrs_spacing_sc)) == (10, [2, 7], 6),
        "MCS": (str(cfg.mcs.table), int(cfg.mcs.index)) == ("nr_256qam", 8),
        "ideal CSI": str(cfg.channel_estimation.ce_method).upper() == "IDEAL",
        "statistics/formal seeds separated": int(fixed.statistics_seed) != int(fixed.realization_seed),
    }
    failed = [label for label, passed in checks.items() if not passed]
    if failed:
        raise ValueError("Plan-037 section 13 scene mismatch: " + ", ".join(failed))
    if require_manifest:
        manifest = load_frozen_manifest(fixed.frozen_beam_manifest, fixed.frozen_beam_manifest_sha256)
        if manifest["frozen_context"] != _context(cfg):
            raise ValueError("Frozen manifest context does not match the link configuration.")


def gate_a(cfg, output: Path) -> None:
    validate_frozen_scene(cfg, require_manifest=False)
    grid = build_resource_grid(cfg.resource)
    channel = FixedCDLStatisticsChannel(cfg, grid, prepare_reference=False)
    codebook = build_dft_2x8_same_pol_codebook()
    mapping = channel.txru_mapping
    gram_b = codebook.conj().T @ codebook
    gram_t = mapping.conj().T @ mapping
    analytic = channel.analytic_transmit_covariance()
    analytic_powers = np.einsum("tb,tu,ub->b", codebook.conj(), analytic, codebook, optimize=True).real
    eig = np.linalg.eigvalsh(analytic)
    rng = np.random.default_rng(20260922)
    h = rng.normal(size=(2, 576, 32)) + 1j * rng.normal(size=(2, 576, 32))
    direct_errors = {}
    for scheme in SCHEMES:
        precoder = build_beam8_precoder(grid, scheme, list(range(8)), codebook, 6)
        direct, reconstructed = equivalent_channel_two_paths(h, precoder.C, codebook[:, :8], scheme)
        direct_errors[scheme] = float(np.max(np.abs(direct - reconstructed)))
    single_ray_cases = {
        "boresight": 0,
        "horizontal_positive_spatial_frequency": 7,
        "horizontal_negative_spatial_frequency": 1,
        "vertical_alternating_phase": 8,
    }
    single_ray_results = {}
    for label, expected in single_ray_cases.items():
        ray = codebook[:, expected].conj()[None, :]
        response = np.abs(ray @ codebook).reshape(-1)
        strongest = np.flatnonzero(np.isclose(response, np.max(response), atol=1e-12)).tolist()
        single_ray_results[label] = {"expected_beam_id": expected, "strongest_beam_ids": strongest,
                                     "peak_response": float(np.max(response))}
    profile = channel.profile_arrays()
    profile_hashes = {name: hashlib.sha256(np.ascontiguousarray(value).view(np.uint8)).hexdigest() for name, value in profile.items()}
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output / "gate_a_numeric.npz",
        txru_mapping=mapping,
        codebook=codebook,
        ae_codebook=channel.txru_mapping @ codebook,
        analytic_covariance=analytic,
        analytic_beam_powers=analytic_powers,
        **profile,
    )
    report = {
        "schema": SCHEMA,
        "stage": "gate_a",
        "status": "PASS" if (
            np.max(np.abs(gram_b - np.eye(16))) <= 1e-12
            and np.max(np.abs(gram_t - np.eye(32))) <= 1e-12
            and float(np.min(eig)) >= -1e-10 * max(float(np.max(eig)), 1.0)
            and max(direct_errors.values()) <= 1e-10
            and all(item["strongest_beam_ids"] == [item["expected_beam_id"]] for item in single_ray_results.values())
        ) else "FAIL",
        "coordinate_semantics": {
            "aod_aoa": "Sionna TR 38.901 GCS azimuth in radians; downlink uses BS departure and UT arrival",
            "zod_zoa": "zenith angle in radians, not elevation",
            "profile_native_angles": True,
            "angle_shift_or_scale_applied": False,
            "txru_flatten_order": "polarization-major, then vertical-major/horizontal-minor",
            "beam_id": "q_v*8+q_h",
            "dft_sign": "negative; H@b transmit semantics",
        },
        "max_codebook_gram_error": float(np.max(np.abs(gram_b - np.eye(16)))),
        "max_mapping_gram_error": float(np.max(np.abs(gram_t - np.eye(32)))),
        "mapping_nonzero_amplitudes": sorted(set(np.abs(mapping[np.nonzero(mapping)]).tolist())),
        "analytic_covariance_hermitian_error": float(np.max(np.abs(analytic - analytic.conj().T))),
        "analytic_covariance_min_eigenvalue": float(np.min(eig)),
        "analytic_covariance_trace": float(np.trace(analytic).real),
        "beam_direct_equivalence_max_errors": direct_errors,
        "single_ray_calibration": single_ray_results,
        "profile_array_sha256": profile_hashes,
        "profile_array_shapes": {name: list(value.shape) for name, value in profile.items()},
        "angle_statistics": channel.angle_statistics,
        "frozen_context": _context(cfg),
        "numeric_npz_sha256": "pending",
    }
    report["numeric_npz_sha256"] = _sha256(output / "gate_a_numeric.npz")
    _json(output / "gate_a_report.json", report)
    if report["status"] != "PASS":
        raise RuntimeError("Gate A failed; inspect gate_a_report.json.")


def statistics(cfg, gate_a_dir: Path, output: Path, start: int, count: int) -> None:
    validate_frozen_scene(cfg, require_manifest=False)
    gate = json.loads((gate_a_dir / "gate_a_report.json").read_text(encoding="utf-8"))
    if gate.get("status") != "PASS" or gate.get("frozen_context") != _context(cfg):
        raise ValueError("Statistics requires a matching PASS Gate A report.")
    if start < 0 or count <= 0 or start + count > 10000:
        raise ValueError("Statistics interval must satisfy 0 <= start < start+count <= 10000.")
    grid = build_resource_grid(cfg.resource)
    channel = FixedCDLStatisticsChannel(cfg, grid, prepare_reference=False)
    powers = np.empty((count, 16), dtype=np.float64)
    covariance = np.empty((count, 32, 32), dtype=np.complex128)
    for local, realization in enumerate(range(start, start + count)):
        covariance[local], powers[local] = channel.realization_beam_statistics(realization)
        if (local + 1) % 50 == 0 or local + 1 == count:
            print(f"[statistics] {local + 1}/{count} interval=[{start},{start + count})", flush=True)
    output.mkdir(parents=True, exist_ok=True)
    name = f"statistics_{start:05d}_{start + count:05d}.npz"
    np.savez_compressed(output / name, realization_index=np.arange(start, start + count), powers=powers, covariance=covariance)
    _json(output / name.replace(".npz", ".json"), {
        "schema": SCHEMA, "stage": "statistics", "absolute_interval": [start, start + count],
        "statistics_seed": int(cfg.fixed_cdl_statistics.statistics_seed), "frozen_context": _context(cfg),
        "npz_sha256": _sha256(output / name),
    })


def _load_statistics(directory: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    pieces = []
    sources = []
    for path in sorted(directory.glob("statistics_*.npz")):
        data = np.load(path)
        indices = np.asarray(data["realization_index"], dtype=np.int64)
        powers = np.asarray(data["powers"], dtype=np.float64)
        pieces.extend((int(index), row.copy()) for index, row in zip(indices, powers))
        sources.append(path.name)
    pieces.sort(key=lambda item: item[0])
    if not pieces or [item[0] for item in pieces] != list(range(len(pieces))):
        raise ValueError("Statistics batches must form one exact absolute interval [0,D) without gaps or duplicates.")
    return np.asarray([item[0] for item in pieces]), np.vstack([item[1] for item in pieces]), sources


def _bootstrap(powers: np.ndarray, repeats: int, seed: int) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    means = powers.mean(axis=0)
    order = np.lexsort((np.arange(16), -means))
    selected = set(order[:8].tolist())
    membership = np.zeros((repeats, 16), dtype=np.float64)
    gaps = np.empty(repeats, dtype=np.float64)
    normalized = np.empty((repeats, 16), dtype=np.float64)
    for repeat in range(repeats):
        sample = powers[rng.integers(0, powers.shape[0], powers.shape[0])].mean(axis=0)
        ranking = np.lexsort((np.arange(16), -sample))
        membership[repeat, ranking[:8]] = 1.0
        gaps[repeat] = sample[ranking[7]] - sample[ranking[8]]
        normalized[repeat] = 10.0 * np.log10(sample / np.max(sample))
    probability = membership.mean(axis=0)
    return {
        "order": order.tolist(), "selected": sorted(selected), "selection_probability": probability.tolist(),
        "gap_99_ci": np.quantile(gaps, [0.005, 0.995]).tolist(),
        "normalized_power_99_ci_db": np.quantile(normalized, [0.005, 0.995], axis=0).T.tolist(),
    }


def freeze(cfg, gate_a_dir: Path, statistics_dir: Path, output: Path, repeats: int) -> None:
    validate_frozen_scene(cfg, require_manifest=False)
    indices, powers, sources = _load_statistics(statistics_dir)
    if powers.shape[0] < 2000 or powers.shape[0] > 10000 or powers.shape[0] % 500:
        raise ValueError("Freeze requires D in {2000,2500,...,10000}.")
    checkpoints = []
    for count in range(2000, powers.shape[0] + 1, 500):
        mean = powers[:count].mean(axis=0)
        order = np.lexsort((np.arange(16), -mean))
        checkpoints.append({"D": count, "top8": sorted(order[:8].tolist()), "order": order.tolist()})
    boot = _bootstrap(powers, repeats, int(cfg.fixed_cdl_statistics.statistics_seed) ^ 0xB00757A9)
    selected = set(boot["selected"])
    stable = len(checkpoints) >= 3 and len({tuple(item["top8"]) for item in checkpoints[-3:]}) == 1
    probabilities = np.asarray(boot["selection_probability"])
    identifiable = (
        stable and all(probabilities[index] >= 0.99 for index in selected)
        and all(probabilities[index] <= 0.01 for index in set(range(16)) - selected)
        and float(boot["gap_99_ci"][0]) > 0.0
    )
    numeric = np.load(gate_a_dir / "gate_a_numeric.npz")
    analytic = np.asarray(numeric["analytic_beam_powers"], dtype=np.float64)
    mean = powers.mean(axis=0)
    analytic_db = 10.0 * np.log10(analytic / np.max(analytic))
    mc_db = 10.0 * np.log10(mean / np.max(mean))
    analytic_order = np.lexsort((np.arange(16), -analytic))
    rank = float(np.corrcoef(np.argsort(np.argsort(-analytic)), np.argsort(np.argsort(-mean)))[0, 1])
    equivalent_top8 = set(analytic_order[:8]) == selected
    interval = np.asarray(boot["normalized_power_99_ci_db"])
    point_ok = np.logical_or(np.abs(analytic_db - mc_db) <= 0.15, (analytic_db >= interval[:, 0]) & (analytic_db <= interval[:, 1]))
    analytic_match = equivalent_top8 and rank >= 0.99 and bool(np.all(point_ok))
    status = "FROZEN" if identifiable and analytic_match else (
        "TOP8_NOT_IDENTIFIABLE" if not identifiable and powers.shape[0] >= 10000 else
        "CDL_ANGLE_POWER_MISMATCH" if identifiable and not analytic_match else "NEEDS_MORE_STATISTICS"
    )
    manifest = {
        "schema": SCHEMA, "status": status, "codebook_type": CODEBOOK_TYPE,
        "selected_beam_indices": boot["selected"], "ranking": boot["order"], "mean_rsrp": mean.tolist(),
        "reference_receive_power": float(np.max(mean)), "D": int(powers.shape[0]),
        "statistics_seed": int(cfg.fixed_cdl_statistics.statistics_seed), "bootstrap_repeats": int(repeats),
        "selection_probability": boot["selection_probability"], "gap_99_ci": boot["gap_99_ci"],
        "normalized_mc_power_db": mc_db.tolist(), "normalized_analytic_power_db": analytic_db.tolist(),
        "normalized_power_99_ci_db": boot["normalized_power_99_ci_db"], "spearman_rank_correlation": rank,
        "analytic_top8_matches": equivalent_top8, "analytic_point_acceptance": point_ok.tolist(),
        "stable_last_three_checkpoints": stable, "checkpoints": checkpoints, "statistics_files": sources,
        "frozen_context": _context(cfg), "gate_a_report_sha256": _sha256(gate_a_dir / "gate_a_report.json"),
    }
    output.mkdir(parents=True, exist_ok=True)
    digest = write_manifest(output / "beam8_frozen_manifest.json", manifest)
    _csv(output / "beam8_rsrp.csv", [
        {"beam_id": index, "q_v": index // 8, "q_h": index % 8, "mean_rsrp": mean[index],
         "normalized_mc_db": mc_db[index], "normalized_analytic_db": analytic_db[index],
         "selection_probability": probabilities[index], "selected": int(index in selected)} for index in range(16)
    ])
    print(f"[freeze] status={status} D={powers.shape[0]} manifest_sha256={digest}", flush=True)
    if status in {"TOP8_NOT_IDENTIFIABLE", "CDL_ANGLE_POWER_MISMATCH"}:
        raise RuntimeError(f"Gate B terminal stop: {status}")


def _pdp_rows(delays: np.ndarray, powers: np.ndarray, **labels: Any) -> list[dict[str, Any]]:
    mean, rms = pdp_moments(delays, powers)
    return [{**labels, "tap": index, "delay_s": float(delay), "delay_ns": float(delay * 1e9),
             "power": float(power), "mean_delay_ns": mean * 1e9, "rms_delay_ns": rms * 1e9}
            for index, (delay, power) in enumerate(zip(delays, powers))]


def pdp(cfg, manifest_path: Path, output: Path) -> None:
    """Complete Gate B with standard, per-beam, reference, and CDD PDP artifacts."""
    validate_frozen_scene(cfg, require_manifest=False)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    embedded = str(manifest.get("manifest_sha256", ""))
    manifest = load_frozen_manifest(manifest_path, embedded)
    if manifest.get("frozen_context") != _context(cfg):
        raise ValueError("PDP extraction requires a matching frozen manifest context.")
    grid = build_resource_grid(cfg.resource)
    channel = FixedCDLStatisticsChannel(cfg, grid, prepare_reference=False)
    codebook = build_dft_2x8_same_pol_codebook()
    delays, path_covariances, component = channel.analytic_delay_covariances()
    beam_component_power = np.einsum(
        "tb,ctu,ub->bc", codebook.conj(), path_covariances, codebook, optimize=True
    ).real
    if np.any(beam_component_power < -1e-10):
        raise RuntimeError("Analytic beam PDP contains negative component power.")
    beam_component_power = np.maximum(beam_component_power, 0.0)
    analytic = np.einsum("tb,tu,ub->b", codebook.conj(), np.sum(path_covariances, axis=0), codebook).real
    totals = beam_component_power.sum(axis=1)
    if not np.allclose(totals, analytic, atol=1e-9, rtol=1e-10):
        raise RuntimeError("Beam PDP totals do not equal analytic beam powers.")

    raw_delays = np.asarray(channel._profile_raw_arrays["delays"], dtype=np.float64).reshape(-1) * float(cfg.channel.delay_spread_ns) * 1e-9
    raw_powers = np.asarray(channel._profile_raw_arrays["powers"], dtype=np.float64).reshape(-1)
    raw_powers = np.power(10.0, raw_powers / 10.0) if np.any(raw_powers < 0.0) else raw_powers
    raw_delays, raw_powers = merge_discrete_pdp(raw_delays, raw_powers)
    expanded_cluster_power = np.asarray(channel.powers, dtype=np.float64).reshape(-1)
    expanded_delays = np.asarray(channel.delays_s, dtype=np.float64).reshape(-1)
    if bool(component["los"]):
        k = float(component["k_factor_linear"])
        final_delays = np.r_[expanded_delays, expanded_delays[0]]
        final_power = np.r_[expanded_cluster_power / (k + 1.0), k / (k + 1.0)]
    else:
        final_delays, final_power = expanded_delays, expanded_cluster_power
    final_delays, final_power = merge_discrete_pdp(final_delays, final_power)

    selected = [int(value) for value in manifest["selected_beam_indices"]]
    strongest = int(next(value for value in manifest["ranking"] if value in selected))
    rows = _pdp_rows(raw_delays, raw_powers, pdp_type="standard_raw", beam_id="")
    beam_tables: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for beam in selected:
        merged_d, merged_p = merge_discrete_pdp(delays, beam_component_power[beam])
        beam_tables[beam] = (merged_d, merged_p)
        rows.extend(_pdp_rows(merged_d, merged_p, pdp_type="beam_normalized", beam_id=beam))
    ref_delays, ref_power = beam_tables[strongest]
    rows.extend(_pdp_rows(ref_delays, ref_power, pdp_type="reference", beam_id=strongest))

    ray_delays, ray_dopplers, ray_covariances = channel.analytic_ray_covariances()
    reference_vector = codebook[:, strongest]
    reference_ray_power = np.einsum(
        "t,ctu,u->c", reference_vector.conj(), ray_covariances, reference_vector, optimize=True
    ).real
    reference_ray_power = np.maximum(reference_ray_power, 0.0)
    reference_ray_power /= np.sum(reference_ray_power)
    symbol_times = np.arange(grid.n_symbols, dtype=np.float64) * float(grid.ofdm_symbol_duration_s)
    time_delta = symbol_times[:, None] - symbol_times[None, :]
    reference_time_covariance = np.sum(
        reference_ray_power[None, None, :]
        * np.exp(1j * 2.0 * np.pi * time_delta[:, :, None] * ray_dopplers[None, None, :]), axis=2
    )

    scs_hz = float(grid.scs_khz) * 1e3
    covariance_checks = {}
    cdd_tables = {}
    ref_cov = frequency_covariance(ref_delays, ref_power, grid.n_sc, scs_hz)
    for scheme, indices in DELAYS.items():
        eff_d, eff_p = shifted_reference_pdp(ref_delays, ref_power, indices, grid.n_sc, scs_hz)
        cdd_tables[scheme] = (eff_d, eff_p)
        direct = frequency_covariance(eff_d, eff_p, grid.n_sc, scs_hz)
        delta = np.arange(grid.n_sc)[:, None] - np.arange(grid.n_sc)[None, :]
        a = np.mean(np.exp(-1j * 2.0 * np.pi * delta[:, :, None] * np.asarray(indices)[None, None, :] / grid.n_sc), axis=2)
        covariance_checks[scheme] = float(np.max(np.abs(direct - ref_cov * a)))
        rows.extend(_pdp_rows(eff_d, eff_p, pdp_type="cdd_effective", beam_id=scheme))
    if max(covariance_checks.values()) > 1e-10:
        raise RuntimeError("Shifted-PDP and Hadamard frequency covariance disagree.")

    output.mkdir(parents=True, exist_ok=True)
    _csv(output / "pdp_tables.csv", rows)
    np.savez_compressed(output / "pdp_numeric.npz", component_delays_s=delays,
                       component_covariances=path_covariances, beam_component_power=beam_component_power,
                       standard_delays_s=final_delays, standard_power=final_power,
                       reference_delays_s=ref_delays, reference_power=ref_power,
                       reference_ray_delays_s=ray_delays, reference_ray_dopplers_hz=ray_dopplers,
                       reference_ray_power=reference_ray_power, reference_time_covariance=reference_time_covariance)
    ref_payload = {"beam_id": strongest, "delays_s": ref_delays.tolist(), "powers": ref_power.tolist(),
                   "ray_delays_s": ray_delays.tolist(), "ray_dopplers_hz": ray_dopplers.tolist(),
                   "ray_powers": reference_ray_power.tolist(),
                   "time_covariance_real": reference_time_covariance.real.tolist(),
                   "time_covariance_imag": reference_time_covariance.imag.tolist()}
    ref_hash = hashlib.sha256(json.dumps(ref_payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    report = {
        "schema": SCHEMA, "stage": "pdp", "status": "PASS", "selected_beam_indices": selected,
        "reference_beam_id": strongest, "reference_pdp_sha256": ref_hash,
        "standard_profile_raw": {"delays_s": raw_delays.tolist(), "powers": raw_powers.tolist()},
        "sionna_expanded": {"delays_s": expanded_delays.tolist(), "powers": expanded_cluster_power.tolist()},
        "standard_final": {"delays_s": final_delays.tolist(), "powers": final_power.tolist()},
        "beam_total_power": {str(beam): float(totals[beam]) for beam in selected},
        "beam_pdp": {str(beam): {"delays_s": beam_tables[beam][0].tolist(), "powers": beam_tables[beam][1].tolist()}
                     for beam in selected},
        "reference_pdp": ref_payload, "frequency_covariance_max_errors": covariance_checks,
        "moments": {}, "numeric_npz_sha256": _sha256(output / "pdp_numeric.npz"),
    }
    period = 1.0 / scs_hz
    for label, table in [("standard", (final_delays, final_power)), ("reference", (ref_delays, ref_power)), *cdd_tables.items()]:
        mean, rms = pdp_moments(*table); cmean, crms = circular_pdp_moments(*table, period)
        report["moments"][label] = {"mean_delay_ns": mean * 1e9, "rms_delay_ns": rms * 1e9,
                                      "circular_mean_ns": cmean * 1e9, "circular_rms_ns": crms * 1e9}
    _json(output / "pdp_report.json", report)
    manifest.update({"gate_b_pdp_status": "PASS", "reference_beam_id": strongest,
                     "reference_pdp": ref_payload, "reference_pdp_sha256": ref_hash,
                     "pdp_report_sha256": _sha256(output / "pdp_report.json")})
    digest = write_manifest(output / "beam8_frozen_manifest.json", manifest)
    print(f"[pdp] status=PASS reference_beam={strongest} manifest_sha256={digest}", flush=True)


def validate_link(cfg) -> None:
    validate_frozen_scene(cfg, require_manifest=True)
    scenario_ids = [str(item.get("scenario_id", item.get("id", ""))) for item in cfg.scenarios]
    if scenario_ids != [SCENARIO]:
        raise ValueError(f"Beam8 link config must contain exactly scenarios: [{{scenario_id: {SCENARIO}}}].")
    grid = build_resource_grid(cfg.resource)
    manifest = load_frozen_manifest(
        cfg.fixed_cdl_statistics.frozen_beam_manifest,
        cfg.fixed_cdl_statistics.frozen_beam_manifest_sha256,
    )
    if manifest.get("gate_b_pdp_status") != "PASS" or not manifest.get("reference_pdp_sha256"):
        raise ValueError("Beam8 link requires the section 13.5 Gate-B PDP-enriched manifest.")
    for scheme in SCHEMES:
        result = build_beam8_precoder(grid, scheme, manifest["selected_beam_indices"], prg_size_rb=6)
        power = np.sum(np.abs(result.C) ** 2, axis=1)
        if not np.allclose(power, 1.0, atol=1e-12, rtol=0.0):
            raise RuntimeError(f"Unit-power audit failed for {scheme}.")
    print(f"[validate] {SCENARIO}: frozen_manifest={manifest['manifest_sha256']} schemes={','.join(SCHEMES)}", flush=True)


def render_link_configs(cfg, manifest_path: Path, output: Path) -> None:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    digest = str(manifest.get("manifest_sha256", ""))
    load_frozen_manifest(manifest_path, digest)
    if manifest.get("gate_b_pdp_status") != "PASS" or not manifest.get("reference_pdp_sha256"):
        raise ValueError("Link configs require the section 13.5 Gate-B PDP-enriched manifest.")
    if manifest.get("frozen_context") != _context(cfg):
        raise ValueError("Cannot render link configs from a manifest with a different frozen context.")
    base = dataclass_to_dict(cfg)
    base["fixed_cdl_statistics"].update({
        "selection_only": False,
        "frozen_beam_manifest": manifest_path.resolve().as_posix(),
        "frozen_beam_manifest_sha256": digest,
    })
    base["simulation"].update({"common_random_numbers": True, "save_trial_metrics": True})
    base["variants"] = []
    for scheme in SCHEMES:
        base["variants"].append({
            "variant_id": f"{scheme.lower()}_ideal",
            "transmission": {"tx_scheme": scheme},
            "channel_estimation": {"ce_method": "IDEAL"},
        })
        estimated_method = "BEAM8_PRG_LMMSE" if scheme == "BEAM8_PRECODER_CYCLING" else "BEAM8_CDD_AWARE_LMMSE"
        base["variants"].append({
            "variant_id": f"{scheme.lower()}_estimated",
            "transmission": {"tx_scheme": scheme},
            "channel_estimation": {"ce_method": estimated_method},
        })
    base["scenarios"] = [{"scenario_id": SCENARIO}]
    definitions = {
        "smoke": {"snr_points_db": [4.0], "n_trials_per_snr": 2, "max_trials_per_snr": 2,
                  "min_block_errors": 0, "run_id": "smoke", "output_dir": "outputs/experiment037_cdl_e_32tx_beam8/smoke"},
        "prescan": {"snr_points_db": [-4,-2,0,2,4,6,8,10,12,14], "n_trials_per_snr": 400,
                    "max_trials_per_snr": 400, "min_block_errors": 0, "run_id": "prescan",
                    "output_dir": "outputs/experiment037_cdl_e_32tx_beam8/prescan"},
        "prescan_extend22": {"snr_points_db": [16,18,20,22], "n_trials_per_snr": 400,
                             "max_trials_per_snr": 400, "min_block_errors": 0,
                             "run_id": "prescan_extend22",
                             "output_dir": "outputs/experiment037_cdl_e_32tx_beam8/prescan_extend22"},
    }
    output.mkdir(parents=True, exist_ok=True)
    receipt = {}
    for name, simulation in definitions.items():
        data = json.loads(json.dumps(base))
        data["simulation"].update(simulation)
        path = output / f"plan037_cdl_e_32tx_2rx_beam8_{name}.yaml"
        path.write_text(yaml.safe_dump(data, sort_keys=False, allow_unicode=True), encoding="utf-8")
        receipt[path.name] = _sha256(path)
    _json(output / "plan037_cdl_e_32tx_2rx_beam8_config_receipt.json", {
        "schema": SCHEMA, "manifest": manifest_path.resolve().as_posix(), "manifest_sha256": digest,
        "config_sha256": receipt,
    })


def render_batch_config(cfg, output: Path, start: int, count: int) -> None:
    validate_frozen_scene(cfg, require_manifest=True)
    if start < 0 or count <= 0 or count % 1000 or start % 1000 or start + count > 50000:
        raise ValueError("Formal batch requires 1000-aligned start/count and stop <= 50000.")
    data = dataclass_to_dict(cfg)
    data["simulation"].update({
        "absolute_trial_start": int(start),
        "n_trials_per_snr": int(count),
        "max_trials_per_snr": int(count),
        "min_block_errors": 0,
        "run_id": f"formal_{start:05d}_{start + count:05d}",
    })
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(yaml.safe_dump(data, sort_keys=False, allow_unicode=True), encoding="utf-8")
    print(f"[render-batch] {output} sha256={_sha256(output)} interval=[{start},{start + count})", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--stage", choices=("gate-a", "statistics", "freeze", "pdp", "render-configs", "render-batch", "validate", "run"), required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--gate-a-dir", type=Path)
    parser.add_argument("--statistics-dir", type=Path)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--count", type=int, default=500)
    parser.add_argument("--bootstrap-repeats", type=int, default=2000)
    parser.add_argument("--manifest", type=Path)
    args = parser.parse_args()
    cfg = load_config(args.config.resolve())
    if args.stage == "gate-a":
        gate_a(cfg, args.output.resolve())
    elif args.stage == "statistics":
        statistics(cfg, args.gate_a_dir.resolve(), args.output.resolve(), args.start, args.count)
    elif args.stage == "freeze":
        freeze(cfg, args.gate_a_dir.resolve(), args.statistics_dir.resolve(), args.output.resolve(), args.bootstrap_repeats)
    elif args.stage == "pdp":
        pdp(cfg, args.manifest.resolve(), args.output.resolve())
    elif args.stage == "render-configs":
        render_link_configs(cfg, args.manifest.resolve(), args.output.resolve())
    elif args.stage == "render-batch":
        render_batch_config(cfg, args.output.resolve(), args.start, args.count)
    elif args.stage == "validate":
        validate_link(cfg)
    else:
        validate_link(cfg)
        output = CDDLinkLevelOrchestrator(cfg).run()
        print(f"[run] output={output}", flush=True)


if __name__ == "__main__":
    main()
