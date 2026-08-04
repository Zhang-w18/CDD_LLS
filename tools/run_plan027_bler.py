"""Run the manifest-gated Plan-027 multi-candidate estimated-CSI BLER stages."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import time
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from scipy.optimize import minimize
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.config import ChannelConfig, ResourceConfig
from cdd_lls.core.mcs import build_tb_layout, get_mcs
from cdd_lls.phy.channel_tdl import generate_sionna_tdl_channel
from cdd_lls.phy.estimators import (
    build_frequency_rmmse_filter,
    tdl_active_frequency_covariance,
)
from cdd_lls.phy.ldpc import SionnaLDPCAdapter
from cdd_lls.phy.precoding import build_active_dft_grid_precoder, equivalent_channel
from cdd_lls.phy.qam import qam_demapper_maxlog, qam_modulate
from cdd_lls.phy.resource_grid import (
    build_resource_grid,
    local_indices_for_subcarriers,
)
from tools.run_plan025_delay_matched_tdl import (
    read_csv_rows,
    stable_seed,
    wilson,
    write_csv_rows,
)
from tools.run_v_design_piecewise_tradeoff import decode_same_tb_batch


DEFAULT_OUTPUT_ROOT = ROOT / "outputs" / "experiment027_meff_sidon"
SCENE_ORDER = ("A5", "A1", "A10", "A30", "A100", "A300")
ORIGINAL_BASELINE_FAMILIES = ("B0_QC", "AP_RMS_T1", "AP_TEPS_T1")
SYSTEM_BASELINE_FAMILIES = ("AP_TU_NT", "AP_TU_NTM1", "AP_TALIAS_NT")
BASELINE_FAMILIES = (*ORIGINAL_BASELINE_FAMILIES, *SYSTEM_BASELINE_FAMILIES)
BASELINE_SHORT_NAMES = {
    "B0_QC": "b0",
    "AP_RMS_T1": "ap_rms",
    "AP_TEPS_T1": "ap_teps",
    "AP_TU_NT": "ap_tu_nt",
    "AP_TU_NTM1": "ap_tu_ntm1",
    "AP_TALIAS_NT": "ap_talias_nt",
}
SEED = 20260727
TARGETS = (("10pct", 0.10, "outage10_snr_db"), ("1pct", 0.01, "outage1_snr_db"))


def validate_approved_manifest(
    manifest_path: Path,
    approved_sha256: str,
    approval_path: Path,
) -> dict:
    payload = manifest_path.read_bytes()
    actual = hashlib.sha256(payload).hexdigest()
    if actual.lower() != str(approved_sha256).strip().lower():
        raise RuntimeError("Approved manifest SHA-256 does not match the file.")
    if not approval_path.exists():
        raise RuntimeError("E4 researcher approval receipt is missing.")
    approval = json.loads(approval_path.read_text(encoding="utf-8"))
    if approval.get("approved") is not True:
        raise RuntimeError("E4 approval receipt is not affirmative.")
    if str(approval.get("manifest_sha256", "")).lower() != actual.lower():
        raise RuntimeError("E4 approval receipt does not match the manifest SHA-256.")
    manifest = json.loads(payload.decode("utf-8"))
    schema = manifest.get("schema")
    curve_count = int(manifest.get("curve_count", -1))
    accepted = {
        "plan027-e4-envelope-link-manifest-v2": 47,
        "plan027-e4-envelope-link-manifest-v3": 49,
        "plan027-e4-envelope-link-manifest-v4": 53,
    }
    dense_schemas = {
        "plan027-e5-a30-dense-dmrs-manifest-v1": "A30",
        "plan027-e6-a100-dense-dmrs-manifest-v1": "A100",
        "plan028-a300-comb6-manifest-v1": "A300",
    }
    if schema in accepted and curve_count != accepted[schema]:
        raise RuntimeError(
            "The approved manifest is not a supported 47/49/53-curve E4 link schema."
        )
    if schema not in accepted and schema not in dense_schemas:
        raise RuntimeError("The approved manifest schema is not supported.")
    spacing = int(manifest.get("dmrs_spacing_subcarriers", -1))
    if schema in accepted and spacing != 24:
        raise RuntimeError("Plan-027 E4 link execution is fixed to comb 24.")
    if schema in dense_schemas:
        dense_scenario = dense_schemas[schema]
        if spacing not in (6, 12):
            raise RuntimeError("Plan-027 E5 dense-DMRS execution requires comb 12 or 6.")
        if set(manifest.get("scene_counts", {})) != {dense_scenario}:
            raise RuntimeError(
                "Plan-027 dense-DMRS manifest scenario does not match its schema."
            )
        if dense_scenario in ("A100", "A300") and spacing != 6:
            raise RuntimeError(f"{dense_scenario} dense-DMRS execution is fixed to comb 6.")
        if not 8 <= curve_count <= 13:
            raise RuntimeError(
                "Plan-027 dense-DMRS execution requires 8 to 13 physical curves."
            )
    physical = manifest.get("physical_definition", {})
    if int(physical.get("phase_denominator", -1)) != 576:
        raise RuntimeError("Plan-027 E4 phase denominator must be 576.")
    if physical.get("delay_input_field") != "delay_grid_coordinates":
        raise RuntimeError("Plan-027 E4 must consume delay_grid_coordinates.")
    return manifest


def validate_prior_scene_approval(
    manifest: dict,
    scenario_id: str,
    approval_path: Path | None,
) -> None:
    execution_order = tuple(manifest.get("execution_order", SCENE_ORDER))
    if scenario_id not in execution_order:
        raise RuntimeError(
            f"{scenario_id} is not authorized by this manifest execution order."
        )
    index = execution_order.index(scenario_id)
    if index == 0:
        return
    previous = execution_order[index - 1]
    if approval_path is None or not approval_path.exists():
        raise RuntimeError(
            f"{scenario_id} requires researcher confirmation of completed scene {previous}."
        )
    approval = json.loads(approval_path.read_text(encoding="utf-8"))
    if approval.get("approved") is not True or approval.get("scenario_id") != previous:
        raise RuntimeError(
            f"Prior-scene approval must affirm completed scene {previous}."
        )


def resource_config(dmrs_spacing_subcarriers: int = 24) -> ResourceConfig:
    return ResourceConfig(
        carrier_bandwidth_mhz=100.0,
        scs_khz=30,
        n_fft=4096,
        n_prbs=48,
        pdsch_n_symbols=10,
        dmrs_symbol_indices=[2, 7],
        dmrs_spacing_sc=int(dmrs_spacing_subcarriers),
        dmrs_offset_sc=0,
        cyclic_prefix_length=288,
    )


def scenario_candidates(manifest: dict, scenario_id: str) -> list[dict]:
    if scenario_id not in manifest.get("scene_counts", {}):
        raise ValueError(f"Unsupported scenario {scenario_id}.")
    rows = [
        dict(row)
        for row in manifest["candidates"]
        if row["scenario_id"] == scenario_id
    ]
    expected = int(manifest["scene_counts"][scenario_id])
    if len(rows) != expected:
        raise RuntimeError(
            f"Manifest contains {len(rows)} curves for {scenario_id}, expected {expected}."
        )
    identifiers = [str(row["candidate_id"]) for row in rows]
    if len(identifiers) != len(set(identifiers)):
        raise RuntimeError(f"Duplicate candidate IDs in scene {scenario_id}.")
    for row in rows:
        coordinates = row.get("delay_grid_coordinates")
        if not isinstance(coordinates, list) or len(coordinates) != 8:
            raise RuntimeError(
                f"{scenario_id}/{row['candidate_id']} has invalid delay coordinates."
            )
        if not all(math.isfinite(float(value)) for value in coordinates):
            raise RuntimeError(
                f"{scenario_id}/{row['candidate_id']} has non-finite delays."
            )
        expected_spacing = int(manifest["dmrs_spacing_subcarriers"])
        if int(row.get("dmrs_spacing_subcarriers", -1)) != expected_spacing:
            raise RuntimeError(
                f"{scenario_id}/{row['candidate_id']} is not comb "
                f"{expected_spacing}."
            )
    return rows


def build_scene(
    manifest: dict,
    candidates: Sequence[dict],
    scenario_id: str,
):
    spread_ns = float(scenario_id.removeprefix("A"))
    resource = resource_config(int(manifest["dmrs_spacing_subcarriers"]))
    channel = ChannelConfig(
        backend="sionna_tdl",
        model="3gpp_tr38901_tdl",
        tdl_profile="A",
        delay_spread_ns=spread_ns,
        carrier_frequency_hz=3.5e9,
        ue_speed_kmh=0.0,
        num_sinusoids=20,
    )
    grid = build_resource_grid(resource)
    precoders = {}
    for row in candidates:
        candidate_id = str(row["candidate_id"])
        coordinates = [float(value) for value in row["delay_grid_coordinates"]]
        precoder = build_active_dft_grid_precoder(
            grid,
            coordinates,
            n_tx=8,
            normalize=False,
        )
        actual = np.asarray(
            precoder.metadata["cdd_delay_indices"], dtype=np.float64
        )
        if not np.array_equal(actual, np.asarray(coordinates, dtype=np.float64)):
            raise RuntimeError(f"Delay coordinates changed for {candidate_id}.")
        if int(precoder.metadata["phase_denominator"]) != 576:
            raise RuntimeError(f"Phase denominator changed for {candidate_id}.")
        precoders[candidate_id] = precoder
    base_covariance = tdl_active_frequency_covariance(grid, channel)
    covariances = {
        candidate_id: base_covariance * (precoder.C @ precoder.C.conj().T)
        for candidate_id, precoder in precoders.items()
    }
    return channel, grid, precoders, covariances


def _snr_key(snr_db: float) -> str:
    return f"{float(snr_db):+.2f}".replace("+", "p").replace("-", "m").replace(".", "p")


def _flag_path(stage_dir: Path, snr_db: float, candidate_id: str) -> Path:
    safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", candidate_id)
    return stage_dir / "error_flags" / _snr_key(snr_db) / f"{safe_id}.npy"


def _write_rows_atomic(rows: Sequence[dict], path: Path) -> None:
    if not rows:
        return
    temporary = path.with_suffix(path.suffix + ".tmp")
    write_csv_rows(rows, temporary)
    temporary.replace(path)


def _pilot_geometry(precoder, pilot_local: np.ndarray) -> dict[str, float | int]:
    matrix = precoder.C[pilot_local]
    singular = np.linalg.svd(matrix, compute_uv=False)
    return {
        "pilot_matrix_rank": int(np.linalg.matrix_rank(matrix, tol=1e-10)),
        "pilot_matrix_condition_number": float(singular[0] / singular[-1]),
        "pilot_matrix_min_singular_value": float(singular[-1]),
    }


def _pairwise_counts(stage_dir: Path, rows: Sequence[dict]) -> None:
    output = []
    grouped: dict[tuple[float, int], list[dict]] = {}
    for row in rows:
        grouped.setdefault((float(row["snr_db"]), int(row["trials"])), []).append(row)
    for (snr_db, trials), values in sorted(grouped.items()):
        flags = {}
        for row in values:
            candidate_id = str(row["candidate_id"])
            path = _flag_path(stage_dir, snr_db, candidate_id)
            if path.exists():
                flags[candidate_id] = np.load(path).astype(bool)
        identifiers = sorted(flags)
        for left_index, left in enumerate(identifiers):
            for right in identifiers[left_index + 1 :]:
                left_flags = flags[left]
                right_flags = flags[right]
                if len(left_flags) != trials or len(right_flags) != trials:
                    raise RuntimeError("Stored error flags do not match the trial budget.")
                output.append(
                    {
                        "snr_db": snr_db,
                        "trials": trials,
                        "candidate_a": left,
                        "candidate_b": right,
                        "both_error": int(np.sum(left_flags & right_flags)),
                        "a_only_error": int(np.sum(left_flags & ~right_flags)),
                        "b_only_error": int(np.sum(~left_flags & right_flags)),
                        "both_success": int(np.sum(~left_flags & ~right_flags)),
                        "mcnemar_exact_pvalue_two_sided": float(
                            binomtest(
                                min(
                                    int(np.sum(left_flags & ~right_flags)),
                                    int(np.sum(~left_flags & right_flags)),
                                ),
                                int(np.sum(left_flags & ~right_flags))
                                + int(np.sum(~left_flags & right_flags)),
                                p=0.5,
                                alternative="two-sided",
                            ).pvalue
                        )
                        if int(np.sum(left_flags ^ right_flags)) > 0
                        else 1.0,
                    }
                )
    if output:
        _write_rows_atomic(output, stage_dir / "paired_error_counts.csv")


def simulate_schedule(
    manifest: dict,
    manifest_sha256: str,
    scenario_id: str,
    schedule: dict[float, set[str]],
    trials: int,
    batch_size: int,
    stage_dir: Path,
    stage_name: str,
) -> list[dict]:
    phase = (
        "E5"
        if str(manifest.get("schema", "")).startswith("plan027-e5-")
        else "E4"
    )
    candidates = scenario_candidates(manifest, scenario_id)
    by_id = {str(row["candidate_id"]): row for row in candidates}
    requested = set().union(*schedule.values()) if schedule else set()
    unknown = requested - set(by_id)
    if unknown:
        raise RuntimeError(f"Unknown candidates requested: {sorted(unknown)}")
    channel, grid, precoders, covariances = build_scene(
        manifest,
        candidates,
        scenario_id,
    )
    mcs = get_mcs("nr_256qam", 8, None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(
        tb.cb_k_values,
        tb.cb_e_values,
        num_iter=8,
        llr_clip=50.0,
    )
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    geometry = {
        candidate_id: _pilot_geometry(precoder, pilot_local)
        for candidate_id, precoder in precoders.items()
    }
    result_path = stage_dir / "bler_points.csv"
    rows = [dict(row) for row in read_csv_rows(result_path)]
    stage_dir.mkdir(parents=True, exist_ok=True)
    receipt = {
        "experiment": (
            "plan-027 E5"
            if str(manifest.get("schema", "")).startswith("plan027-e5-")
            else "plan-027 E4"
        ),
        "stage": stage_name,
        "scenario_id": scenario_id,
        "profile": "A",
        "delay_spread_ns": float(scenario_id.removeprefix("A")),
        "dmrs_spacing_subcarriers": int(
            manifest["dmrs_spacing_subcarriers"]
        ),
        "pilot_re": int(grid.n_dmrs_re),
        "data_re": int(grid.n_data_re),
        "trials_per_snr": int(trials),
        "batch_size": int(batch_size),
        "seed": SEED,
        "manifest_sha256": manifest_sha256,
        "phase_denominator": 576,
        "phase_reference": "first_active_subcarrier",
        "candidate_delay_grid_coordinates": {
            candidate_id: by_id[candidate_id]["delay_grid_coordinates"]
            for candidate_id in sorted(by_id)
        },
        "receiver": "two-DMRS averaged, V-aware matched frequency LMMSE",
    }
    (stage_dir / "resolved_experiment.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    for snr_db in sorted(schedule):
        candidate_ids = sorted(schedule[snr_db])
        completed = {
            str(row["candidate_id"])
            for row in rows
            if math.isclose(float(row["snr_db"]), float(snr_db), abs_tol=1e-12)
            and int(row["trials"]) == int(trials)
            and _flag_path(
                stage_dir, snr_db, str(row["candidate_id"])
            ).exists()
        }
        pending = [candidate_id for candidate_id in candidate_ids if candidate_id not in completed]
        if not pending:
            print(
                f"[{phase} {scenario_id} {stage_name}] resume skip snr={snr_db:g} "
                f"candidates={len(candidate_ids)}",
                flush=True,
            )
            continue
        rows = [
            row
            for row in rows
            if not (
                math.isclose(float(row["snr_db"]), float(snr_db), abs_tol=1e-12)
                and str(row["candidate_id"]) in pending
            )
        ]
        snr_linear = 10.0 ** (float(snr_db) / 10.0)
        noise_variance = 8.0 / snr_linear
        ls_noise_variance = noise_variance / 2.0
        filters = {
            candidate_id: build_frequency_rmmse_filter(
                covariances[candidate_id],
                pilot_local,
                ls_noise_variance,
                diagonal_loading=1e-10,
            )
            for candidate_id in pending
        }
        errors = {candidate_id: 0 for candidate_id in pending}
        error_flags = {
            candidate_id: np.zeros(int(trials), dtype=bool)
            for candidate_id in pending
        }
        nmse_sum = {candidate_id: 0.0 for candidate_id in pending}
        started = time.time()
        completed_trials = 0
        for batch_start in range(1, int(trials) + 1, int(batch_size)):
            current = min(int(batch_size), int(trials) - batch_start + 1)
            payload_rng = np.random.default_rng(
                stable_seed(SEED, scenario_id, snr_db, batch_start, "payload")
            )
            payload = [
                payload_rng.integers(0, 2, size=int(k), dtype=np.int8)
                for k in tb.cb_k_values
            ]
            coded = np.concatenate(adapter.encode(payload))
            symbols = qam_modulate(coded, int(mcs.qm))
            realization = generate_sionna_tdl_channel(
                grid,
                channel,
                n_tx=8,
                n_rx=1,
                batch_size=current,
                seed=stable_seed(
                    SEED, scenario_id, snr_db, batch_start, "channel"
                ),
            )
            noise_rng = np.random.default_rng(
                stable_seed(SEED, scenario_id, snr_db, batch_start, "noise")
            )
            averaged_ls_noise = math.sqrt(ls_noise_variance / 2.0) * (
                noise_rng.normal(size=(current, 1, grid.pilot_count))
                + 1j * noise_rng.normal(size=(current, 1, grid.pilot_count))
            )
            data_noise = math.sqrt(noise_variance / 2.0) * (
                noise_rng.normal(size=(current, 1, grid.n_data_re))
                + 1j * noise_rng.normal(size=(current, 1, grid.n_data_re))
            )
            for candidate_id in pending:
                effective = equivalent_channel(
                    realization.H, precoders[candidate_id].C
                )[:, 0:1]
                pilot0 = effective[
                    :, :, int(grid.pilot_symbol_indices[0]), pilot_local
                ]
                pilot1 = effective[
                    :, :, int(grid.pilot_symbol_indices[-1]), pilot_local
                ]
                true_pilot_average = 0.5 * (pilot0 + pilot1)
                estimate_full = filters[candidate_id].estimate_full_band(
                    true_pilot_average + averaged_ls_noise
                )
                estimate_data = estimate_full[:, :, data_local]
                true_data = effective[
                    :, :, grid.data_symbol_indices, data_local
                ]
                per_trial_nmse = np.sum(
                    np.abs(estimate_data - true_data) ** 2, axis=(1, 2)
                ) / np.maximum(
                    np.sum(np.abs(true_data) ** 2, axis=(1, 2)), 1e-30
                )
                nmse_sum[candidate_id] += float(np.sum(per_trial_nmse))
                received = true_data * symbols[None, None, :] + data_noise
                denominator = np.maximum(
                    np.sum(np.abs(estimate_data) ** 2, axis=1), 1e-10
                )
                equalized = np.sum(
                    np.conj(estimate_data) * received, axis=1
                ) / denominator
                effective_noise = noise_variance / denominator
                llrs = [
                    qam_demapper_maxlog(
                        equalized[index],
                        effective_noise[index],
                        int(mcs.qm),
                    )
                    for index in range(current)
                ]
                decoded = decode_same_tb_batch(adapter, llrs, payload)
                flags = np.asarray(
                    [not item.tb_success for item in decoded], dtype=bool
                )
                start = batch_start - 1
                error_flags[candidate_id][start : start + current] = flags
                errors[candidate_id] += int(np.sum(flags))
            completed_trials += current
            if completed_trials % 100 < current or completed_trials == int(trials):
                print(
                    f"[{phase} {scenario_id} {stage_name}] snr={snr_db:g} "
                    f"trial={completed_trials}/{trials} candidates={len(pending)} "
                    f"errors={errors} elapsed={time.time()-started:.1f}s",
                    flush=True,
                )

        for candidate_id in pending:
            flag_path = _flag_path(stage_dir, snr_db, candidate_id)
            flag_path.parent.mkdir(parents=True, exist_ok=True)
            np.save(flag_path, error_flags[candidate_id])
            low, high = wilson(errors[candidate_id], int(trials))
            mean_nmse = nmse_sum[candidate_id] / float(trials)
            filt = filters[candidate_id]
            row = {
                "scenario_id": scenario_id,
                "stage": stage_name,
                "candidate_id": candidate_id,
                "family": by_id[candidate_id]["family"],
                "snr_db": float(snr_db),
                "trials": int(trials),
                "tb_errors": int(errors[candidate_id]),
                "bler": float(errors[candidate_id]) / float(trials),
                "bler_wilson95_lo": low,
                "bler_wilson95_hi": high,
                "ce_nmse_mean": mean_nmse,
                "ce_nmse_mean_db": 10.0 * math.log10(mean_nmse),
                "estimator_condition_number": filt.condition_number,
                "estimator_min_singular_value": filt.minimum_singular_value,
                "numerical_jitter": filt.numerical_jitter,
                "noise_variance": noise_variance,
                "ls_noise_variance": ls_noise_variance,
                **geometry[candidate_id],
            }
            rows.append(row)
        rows.sort(key=lambda row: (float(row["snr_db"]), str(row["candidate_id"])))
        _write_rows_atomic(rows, result_path)
        _pairwise_counts(stage_dir, rows)
    receipt["executed_snrs_by_candidate"] = {
        candidate_id: sorted(
            {
                float(row["snr_db"])
                for row in rows
                if str(row["candidate_id"]) == candidate_id
            }
        )
        for candidate_id in sorted(by_id)
    }
    (stage_dir / "resolved_experiment.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return rows


def _round_half(value: float) -> float:
    return round(float(value) * 2.0) / 2.0


def find_crossing_interval(
    rows: Iterable[dict],
    candidate_id: str,
    target: float,
    lower_bound: float | None = None,
    upper_bound: float | None = None,
) -> tuple[float, float] | None:
    selected = sorted(
        (
            row
            for row in rows
            if str(row["candidate_id"]) == candidate_id
            and (lower_bound is None or float(row["snr_db"]) >= lower_bound - 1e-12)
            and (upper_bound is None or float(row["snr_db"]) <= upper_bound + 1e-12)
        ),
        key=lambda row: float(row["snr_db"]),
    )
    for left, right in zip(selected, selected[1:]):
        left_probability = float(left["bler"])
        right_probability = float(right["bler"])
        if left_probability >= target and right_probability <= target:
            return float(left["snr_db"]), float(right["snr_db"])
    return None


def _read_stage_rows(stage_dir: Path) -> list[dict]:
    return [dict(row) for row in read_csv_rows(stage_dir / "bler_points.csv")]


def run_smoke(
    manifest: dict,
    digest: str,
    scenario_id: str,
    batch_size: int,
    scene_root: Path,
) -> None:
    execution_order = tuple(manifest.get("execution_order", SCENE_ORDER))
    if not execution_order or scenario_id != execution_order[0]:
        raise RuntimeError(
            "The integration preflight must use the first authorized scene."
        )
    candidates = scenario_candidates(manifest, scenario_id)
    schedule = {}
    for row in candidates:
        snr_db = _round_half(float(row["outage10_snr_db"]) + 5.0)
        schedule.setdefault(snr_db, set()).add(str(row["candidate_id"]))
    rows = simulate_schedule(
        manifest,
        digest,
        scenario_id,
        schedule,
        20,
        batch_size,
        scene_root / "smoke",
        "smoke",
    )
    if len({str(row["candidate_id"]) for row in rows}) != len(candidates):
        raise RuntimeError(
            "Smoke did not produce one point for every scene candidate."
        )
    if any(not math.isfinite(float(row["ce_nmse_mean_db"])) for row in rows):
        raise RuntimeError("Smoke produced a non-finite CE result.")


def run_prescan(
    manifest: dict,
    digest: str,
    scenario_id: str,
    batch_size: int,
    scene_root: Path,
) -> None:
    candidates = scenario_candidates(manifest, scenario_id)
    budget = manifest["trial_budget"]
    gap = float(budget["initial_bler_gap_from_outage_db"])
    step = float(budget["prescan_step_db"])
    maximum_extension = float(budget["maximum_extension_each_target_db"])
    trials = int(budget["prescan_trials_per_snr"])
    bounds = {}
    schedule: dict[float, set[str]] = {}
    for row in candidates:
        candidate_id = str(row["candidate_id"])
        for label, target, field in TARGETS:
            center = _round_half(float(row[field]) + gap)
            lower = center - step
            upper = center + step
            bounds[(candidate_id, label)] = {
                "target_probability": target,
                "initial_center_db": center,
                "lower_db": lower,
                "upper_db": upper,
                "extensions": 0,
            }
            for snr_db in np.arange(lower, upper + 0.01, step):
                schedule.setdefault(float(round(snr_db, 6)), set()).add(candidate_id)
    stage_dir = scene_root / "prescan"
    simulate_schedule(
        manifest,
        digest,
        scenario_id,
        schedule,
        trials,
        batch_size,
        stage_dir,
        "prescan",
    )
    maximum_extensions = int(round(maximum_extension / step))
    while True:
        rows = _read_stage_rows(stage_dir)
        additions: dict[float, set[str]] = {}
        for (candidate_id, label), spec in bounds.items():
            interval = find_crossing_interval(
                rows,
                candidate_id,
                float(spec["target_probability"]),
                float(spec["lower_db"]),
                float(spec["upper_db"]),
            )
            if interval is not None or int(spec["extensions"]) >= maximum_extensions:
                continue
            candidate_rows = [
                row
                for row in rows
                if str(row["candidate_id"]) == candidate_id
                and float(spec["lower_db"]) - 1e-12
                <= float(row["snr_db"])
                <= float(spec["upper_db"]) + 1e-12
            ]
            at_lower = min(candidate_rows, key=lambda row: float(row["snr_db"]))
            at_upper = max(candidate_rows, key=lambda row: float(row["snr_db"]))
            target = float(spec["target_probability"])
            if float(at_lower["bler"]) < target:
                new_snr = float(spec["lower_db"]) - step
                spec["lower_db"] = new_snr
            elif float(at_upper["bler"]) > target:
                new_snr = float(spec["upper_db"]) + step
                spec["upper_db"] = new_snr
            else:
                break
            spec["extensions"] = int(spec["extensions"]) + 1
            additions.setdefault(float(round(new_snr, 6)), set()).add(candidate_id)
        if not additions:
            break
        simulate_schedule(
            manifest,
            digest,
            scenario_id,
            additions,
            trials,
            batch_size,
            stage_dir,
            "prescan",
        )
    rows = _read_stage_rows(stage_dir)
    bracket_summary = []
    for (candidate_id, label), spec in sorted(bounds.items()):
        interval = find_crossing_interval(
            rows,
            candidate_id,
            float(spec["target_probability"]),
            float(spec["lower_db"]),
            float(spec["upper_db"]),
        )
        bracket_summary.append(
            {
                "candidate_id": candidate_id,
                "target": label,
                **spec,
                "available": interval is not None,
                "crossing_interval_db": list(interval) if interval else None,
            }
        )
    (stage_dir / "prescan_brackets.json").write_text(
        json.dumps(
            {
                "scenario_id": scenario_id,
                "method": "candidate-specific 0.5 dB adaptive bracket",
                "brackets": bracket_summary,
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )


def run_refine(
    manifest: dict,
    digest: str,
    scenario_id: str,
    target_label: str,
    batch_size: int,
    scene_root: Path,
) -> None:
    target_lookup = {label: probability for label, probability, _ in TARGETS}
    if target_label not in target_lookup:
        raise ValueError(f"Unknown target label {target_label}.")
    bracket_path = scene_root / "prescan" / "prescan_brackets.json"
    if not bracket_path.exists():
        raise RuntimeError("Prescan brackets are missing.")
    brackets = json.loads(bracket_path.read_text(encoding="utf-8"))["brackets"]
    schedule: dict[float, set[str]] = {}
    candidate_ids = set()
    for spec in brackets:
        if spec["target"] != target_label or not spec["available"]:
            continue
        candidate_id = str(spec["candidate_id"])
        candidate_ids.add(candidate_id)
        lower, upper = [float(value) for value in spec["crossing_interval_db"]]
        for snr_db in np.arange(lower, upper + 0.001, 0.25):
            schedule.setdefault(float(round(snr_db, 6)), set()).add(candidate_id)
    if not schedule:
        raise RuntimeError(f"No {target_label} candidate brackets are available.")
    trials = int(manifest["trial_budget"]["refine_trials_per_snr"])
    simulate_schedule(
        manifest,
        digest,
        scenario_id,
        schedule,
        trials,
        batch_size,
        scene_root / f"refine_{target_label}",
        f"refine_{target_label}",
    )
    stage_dir = scene_root / f"refine_{target_label}"
    step = 0.25
    maximum_extension = float(
        manifest["trial_budget"]["maximum_extension_each_target_db"]
    )
    maximum_extensions = int(round(maximum_extension / step))
    extension_counts = {candidate_id: 0 for candidate_id in candidate_ids}
    while True:
        rows = _read_stage_rows(stage_dir)
        additions: dict[float, set[str]] = {}
        for candidate_id in sorted(candidate_ids):
            selected = sorted(
                (
                    row
                    for row in rows
                    if str(row["candidate_id"]) == candidate_id
                ),
                key=lambda row: float(row["snr_db"]),
            )
            interval = find_crossing_interval(
                selected,
                candidate_id,
                float(target_lookup[target_label]),
            )
            if interval is not None:
                fitted = logistic_target(
                    selected,
                    candidate_id,
                    float(target_lookup[target_label]),
                )
                fitted_snr = float(fitted["target_snr_db"])
                sampled_lower = float(selected[0]["snr_db"])
                sampled_upper = float(selected[-1]["snr_db"])
                if sampled_lower - 1e-12 <= fitted_snr <= sampled_upper + 1e-12:
                    continue
            if extension_counts[candidate_id] >= maximum_extensions:
                continue
            target = float(target_lookup[target_label])
            at_lower = selected[0]
            at_upper = selected[-1]
            if interval is not None and fitted_snr < sampled_lower:
                new_snr = sampled_lower - step
            elif interval is not None and fitted_snr > sampled_upper:
                new_snr = sampled_upper + step
            elif float(at_lower["bler"]) < target:
                new_snr = float(at_lower["snr_db"]) - step
            elif float(at_upper["bler"]) > target:
                new_snr = float(at_upper["snr_db"]) + step
            else:
                continue
            extension_counts[candidate_id] += 1
            additions.setdefault(float(round(new_snr, 6)), set()).add(candidate_id)
        if not additions:
            break
        simulate_schedule(
            manifest,
            digest,
            scenario_id,
            additions,
            trials,
            batch_size,
            stage_dir,
            f"refine_{target_label}",
        )
    rows = _read_stage_rows(stage_dir)
    bracket_summary = []
    for candidate_id in sorted(candidate_ids):
        selected = sorted(
            (
                row
                for row in rows
                if str(row["candidate_id"]) == candidate_id
            ),
            key=lambda row: float(row["snr_db"]),
        )
        interval = find_crossing_interval(
            selected,
            candidate_id,
            float(target_lookup[target_label]),
        )
        bracket_summary.append(
            {
                "candidate_id": candidate_id,
                "target": target_label,
                "target_probability": float(target_lookup[target_label]),
                "extensions": extension_counts[candidate_id],
                "available": interval is not None,
                "crossing_interval_db": list(interval) if interval else None,
                "sampled_snr_db": [
                    float(row["snr_db"])
                    for row in selected
                ],
            }
        )
    (stage_dir / "refine_brackets.json").write_text(
        json.dumps(
            {
                "scenario_id": scenario_id,
                "method": "formal 3000-trial 0.25 dB adaptive bracket audit",
                "brackets": bracket_summary,
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )


def logistic_target(
    rows: Sequence[dict],
    candidate_id: str,
    target: float,
) -> dict[str, object]:
    selected = sorted(
        (row for row in rows if str(row["candidate_id"]) == candidate_id),
        key=lambda row: float(row["snr_db"]),
    )
    x = np.asarray([float(row["snr_db"]) for row in selected], dtype=np.float64)
    n = np.asarray([int(row["trials"]) for row in selected], dtype=np.float64)
    y = np.asarray([int(row["tb_errors"]) for row in selected], dtype=np.float64)
    if len(x) < 2 or not np.any(y > 0) or not np.any(y < n):
        raise ValueError(f"Insufficient finite-error points for {candidate_id}.")
    matrix = np.column_stack((np.ones_like(x), x))

    def negative_log_likelihood(theta):
        linear = matrix @ theta
        return np.sum(n * np.logaddexp(0.0, linear) - y * linear)

    fit = minimize(negative_log_likelihood, [15.0, -1.0], method="BFGS")
    if not fit.success and not np.all(np.isfinite(fit.x)):
        raise RuntimeError(f"Logistic fit failed for {candidate_id}: {fit.message}")
    intercept, slope = fit.x
    fitted_probability = 1.0 / (1.0 + np.exp(-(matrix @ fit.x)))
    weight = n * fitted_probability * (1.0 - fitted_probability)
    covariance = np.linalg.inv(matrix.T @ (weight[:, None] * matrix))
    logit_target = math.log(target / (1.0 - target))
    snr_db = float((logit_target - intercept) / slope)
    gradient = np.asarray([-1.0 / slope, -snr_db / slope])
    standard_error = float(np.sqrt(gradient @ covariance @ gradient))
    near_errors = sum(
        int(row["tb_errors"])
        for row in selected
        if 0.5 * target <= float(row["bler"]) <= 2.0 * target
    )
    return {
        "target_snr_db": snr_db,
        "standard_error_db": standard_error,
        "ci95_low_db": snr_db - 1.96 * standard_error,
        "ci95_high_db": snr_db + 1.96 * standard_error,
        "logit_intercept": float(intercept),
        "logit_slope_per_db": float(slope),
        "point_count": len(selected),
        "total_trials": int(np.sum(n)),
        "total_errors": int(np.sum(y)),
        "errors_in_half_to_twice_target_band": int(near_errors),
        "sufficient_error_count": bool(near_errors >= 30),
    }


def interpolate_ce_nmse_db(
    rows: Sequence[dict],
    candidate_id: str,
    target_snr_db: float,
) -> float:
    selected = sorted(
        (row for row in rows if str(row["candidate_id"]) == candidate_id),
        key=lambda row: float(row["snr_db"]),
    )
    x = np.asarray([float(row["snr_db"]) for row in selected], dtype=np.float64)
    y = np.asarray(
        [float(row["ce_nmse_mean_db"]) for row in selected],
        dtype=np.float64,
    )
    if not len(x) or target_snr_db < x[0] - 1e-12 or target_snr_db > x[-1] + 1e-12:
        raise ValueError(
            f"CE interpolation point is outside sampled support for {candidate_id}."
        )
    return float(np.interp(float(target_snr_db), x, y))


def curve_style(candidate_id: str, family: str) -> dict[str, object]:
    colors = {
        "B0_QC": "#111111",
        "AP_RMS_T1": "#7f7f7f",
        "AP_TEPS_T1": "#bcbd22",
        "AP_TU_NT": "#d62728",
        "AP_TU_NTM1": "#8c564b",
        "AP_TALIAS_NT": "#17becf",
        "S0_SIDON": "#9467bd",
        "AP_T2_CTRL": "#ff7f0e",
        "GEO_T1_CTRL": "#2ca02c",
        "MEFF_T2_CAND": "#1f77b4",
    }
    variants = (
        ("-", "o"),
        ("--", "s"),
        ("-.", "^"),
        (":", "D"),
        ((0, (5, 1)), "v"),
        ((0, (3, 1, 1, 1)), "P"),
        ((0, (1, 1)), "X"),
    )
    match = re.search(r"_(\d+)$", str(candidate_id))
    variant = int(match.group(1)) if match else 0
    linestyle, marker = variants[variant % len(variants)]
    return {
        "color": colors.get(str(family), "#8c564b"),
        "linestyle": linestyle,
        "marker": marker,
    }


def curve_label(candidate_id: str) -> str:
    return re.sub(r"^A\d+_", "", str(candidate_id))


def _format_signed(value: object) -> str:
    number = float(value)
    return f"{number:+.3f}" if math.isfinite(number) else "未闭合"


def _plot_scene(
    summaries: Sequence[dict],
    rows_by_target: dict[str, list[dict]],
    path: Path,
    dmrs_spacing_subcarriers: int = 24,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    identifiers = sorted({str(row["candidate_id"]) for row in summaries})
    family_by_id = {
        str(row["candidate_id"]): str(row["family"]) for row in summaries
    }
    figure, axes = plt.subplots(1, 2, figsize=(15.0, 6.0))
    for axis, (label, target, _) in zip(axes, TARGETS):
        rows = rows_by_target.get(label, [])
        for candidate_id in identifiers:
            selected = sorted(
                (
                    row
                    for row in rows
                    if str(row["candidate_id"]) == candidate_id
                ),
                key=lambda row: float(row["snr_db"]),
            )
            if not selected:
                continue
            x = np.asarray([float(row["snr_db"]) for row in selected])
            y = np.asarray(
                [
                    max(float(row["bler"]), 0.5 / int(row["trials"]))
                    for row in selected
                ]
            )
            axis.plot(
                x,
                y,
                linewidth=2.5,
                markersize=8,
                label=curve_label(candidate_id),
                **curve_style(candidate_id, family_by_id[candidate_id]),
            )
        axis.axhline(target, color="#555555", linestyle="--", linewidth=2.5)
        axis.set_yscale("log")
        axis.set_xlabel("SNR (dB)", fontsize=16)
        axis.set_ylabel("estimated-CSI BLER", fontsize=16)
        axis.set_title(
            f"comb-{int(dmrs_spacing_subcarriers)} {label} refinement "
            "(closed targets)"
            , fontsize=18
        )
        axis.tick_params(axis="both", which="both", labelsize=14)
        axis.grid(True, which="both", alpha=0.25)
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(
        handles, labels, loc="center left", bbox_to_anchor=(0.99, 0.5), fontsize=16
    )
    figure.tight_layout(rect=(0.0, 0.0, 0.82, 1.0))
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def _plot_prescan(
    rows: Sequence[dict],
    scenario_id: str,
    path: Path,
    dmrs_spacing_subcarriers: int = 24,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    identifiers = sorted({str(row["candidate_id"]) for row in rows})
    family_by_id = {
        str(row["candidate_id"]): str(row["family"]) for row in rows
    }
    figure, axis = plt.subplots(figsize=(12.0, 7.0))
    for candidate_id in identifiers:
        selected = sorted(
            (
                row
                for row in rows
                if str(row["candidate_id"]) == candidate_id
            ),
            key=lambda row: float(row["snr_db"]),
        )
        x = np.asarray([float(row["snr_db"]) for row in selected])
        y = np.asarray(
            [
                max(float(row["bler"]), 0.5 / int(row["trials"]))
                for row in selected
            ]
        )
        axis.plot(
            x,
            y,
            linewidth=2.5,
            markersize=8,
            label=curve_label(candidate_id),
            **curve_style(candidate_id, family_by_id[candidate_id]),
        )
    axis.axhline(0.10, color="#555555", linestyle="--", linewidth=2.5)
    axis.axhline(0.01, color="#555555", linestyle=":", linewidth=2.5)
    axis.set_yscale("log")
    axis.set_xlabel("SNR (dB)", fontsize=16)
    axis.set_ylabel("estimated-CSI BLER", fontsize=16)
    axis.set_title(
        f"{scenario_id} comb-{int(dmrs_spacing_subcarriers)} "
        "candidate-specific prescan (all physical curves)"
        , fontsize=18
    )
    axis.tick_params(axis="both", which="both", labelsize=14)
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=16)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def analyze_scene(
    manifest: dict,
    digest: str,
    scenario_id: str,
    scene_root: Path,
) -> None:
    candidates = scenario_candidates(manifest, scenario_id)
    by_id = {str(row["candidate_id"]): row for row in candidates}
    summaries = []
    rows_by_target = {}
    failures = []
    for label, target, _ in TARGETS:
        stage_dir = scene_root / f"refine_{label}"
        rows = _read_stage_rows(stage_dir)
        _pairwise_counts(stage_dir, rows)
        rows_by_target[label] = rows
        bracket_path = stage_dir / "refine_brackets.json"
        if not bracket_path.exists():
            raise RuntimeError(f"Formal bracket audit is missing for {label}.")
        bracket_rows = json.loads(
            bracket_path.read_text(encoding="utf-8")
        )["brackets"]
        available = {
            str(row["candidate_id"])
            for row in bracket_rows
            if bool(row["available"])
        }
        for candidate_id in sorted(by_id):
            if candidate_id not in available:
                failures.append(
                    {
                        "candidate_id": candidate_id,
                        "target": label,
                        "reason": (
                            "No formal 3000-trial two-sided target bracket; "
                            "target SNR was not extrapolated."
                        ),
                    }
                )
                continue
            try:
                fit = logistic_target(rows, candidate_id, target)
                target_ce_nmse_db = interpolate_ce_nmse_db(
                    rows,
                    candidate_id,
                    float(fit["target_snr_db"]),
                )
            except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
                failures.append(
                    {
                        "candidate_id": candidate_id,
                        "target": label,
                        "reason": str(error),
                    }
                )
                continue
            summaries.append(
                {
                    "scenario_id": scenario_id,
                    "target": label,
                    "target_probability": target,
                    "candidate_id": candidate_id,
                    "family": by_id[candidate_id]["family"],
                    "outage_target_snr_db": float(
                        by_id[candidate_id][
                            "outage10_snr_db" if label == "10pct" else "outage1_snr_db"
                        ]
                    ),
                    "target_ce_nmse_db": target_ce_nmse_db,
                    **fit,
                }
            )
    final_dir = scene_root / "final"
    final_dir.mkdir(parents=True, exist_ok=True)
    if summaries:
        _write_rows_atomic(summaries, final_dir / "target_summary.csv")
    comparisons = []
    for label, _, _ in TARGETS:
        target_rows = [row for row in summaries if row["target"] == label]
        baselines = {
            family: next(
                (
                    row
                    for row in target_rows
                    if row["family"] == family
                ),
                None,
            )
            for family in BASELINE_FAMILIES
        }
        available_baselines = {
            family: row for family, row in baselines.items() if row is not None
        }
        if not available_baselines:
            continue
        original_available = {
            family: row
            for family, row in available_baselines.items()
            if family in ORIGINAL_BASELINE_FAMILIES
        }
        best_three_family, best_three_row = min(
            original_available.items(),
            key=lambda item: float(item[1]["target_snr_db"]),
        )
        best_five_family, best_five_row = min(
            available_baselines.items(),
            key=lambda item: float(item[1]["target_snr_db"]),
        )
        for candidate_row in target_rows:
            comparison = {
                "scenario_id": scenario_id,
                "target": label,
                "candidate_id": candidate_row["candidate_id"],
                "family": candidate_row["family"],
                "candidate_target_snr_db": candidate_row["target_snr_db"],
                "candidate_target_ce_nmse_db": candidate_row[
                    "target_ce_nmse_db"
                ],
                "estimated_csi_best_three_family": best_three_family,
                "estimated_csi_best_three_candidate_id": best_three_row[
                    "candidate_id"
                ],
                "estimated_csi_best_five_family": best_five_family,
                "estimated_csi_best_five_candidate_id": best_five_row[
                    "candidate_id"
                ],
            }
            for family, short in BASELINE_SHORT_NAMES.items():
                comparison[f"gain_vs_{short}_db"] = math.nan
                comparison[f"gain_vs_{short}_ci95_low_db"] = math.nan
                comparison[f"gain_vs_{short}_ci95_high_db"] = math.nan
                comparison[f"ce_delta_vs_{short}_db"] = math.nan
            for family, baseline_row in available_baselines.items():
                short = BASELINE_SHORT_NAMES[family]
                gain = float(baseline_row["target_snr_db"]) - float(
                    candidate_row["target_snr_db"]
                )
                gain_se = math.sqrt(
                    float(baseline_row["standard_error_db"]) ** 2
                    + float(candidate_row["standard_error_db"]) ** 2
                )
                comparison[f"gain_vs_{short}_db"] = gain
                comparison[f"gain_vs_{short}_ci95_low_db"] = gain - 1.96 * gain_se
                comparison[f"gain_vs_{short}_ci95_high_db"] = gain + 1.96 * gain_se
                comparison[f"ce_delta_vs_{short}_db"] = float(
                    candidate_row["target_ce_nmse_db"]
                ) - float(baseline_row["target_ce_nmse_db"])
            best_gain = float(best_three_row["target_snr_db"]) - float(
                candidate_row["target_snr_db"]
            )
            best_gain_se = math.sqrt(
                float(best_three_row["standard_error_db"]) ** 2
                + float(candidate_row["standard_error_db"]) ** 2
            )
            comparison["gain_vs_best_three_db"] = best_gain
            comparison["gain_vs_best_three_ci95_low_db"] = (
                best_gain - 1.96 * best_gain_se
            )
            comparison["gain_vs_best_three_ci95_high_db"] = (
                best_gain + 1.96 * best_gain_se
            )
            comparison["ce_delta_vs_best_three_db"] = float(
                candidate_row["target_ce_nmse_db"]
            ) - float(best_three_row["target_ce_nmse_db"])
            best_five_gain = float(best_five_row["target_snr_db"]) - float(
                candidate_row["target_snr_db"]
            )
            best_five_gain_se = math.sqrt(
                float(best_five_row["standard_error_db"]) ** 2
                + float(candidate_row["standard_error_db"]) ** 2
            )
            comparison["gain_vs_best_five_db"] = best_five_gain
            comparison["gain_vs_best_five_ci95_low_db"] = (
                best_five_gain - 1.96 * best_five_gain_se
            )
            comparison["gain_vs_best_five_ci95_high_db"] = (
                best_five_gain + 1.96 * best_five_gain_se
            )
            comparison["ce_delta_vs_best_five_db"] = float(
                candidate_row["target_ce_nmse_db"]
            ) - float(best_five_row["target_ce_nmse_db"])
            comparison["estimated_csi_best_all_family"] = best_five_family
            comparison["estimated_csi_best_all_candidate_id"] = best_five_row[
                "candidate_id"
            ]
            comparison["gain_vs_best_all_db"] = best_five_gain
            comparison["gain_vs_best_all_ci95_low_db"] = (
                best_five_gain - 1.96 * best_five_gain_se
            )
            comparison["gain_vs_best_all_ci95_high_db"] = (
                best_five_gain + 1.96 * best_five_gain_se
            )
            comparison["ce_delta_vs_best_all_db"] = comparison[
                "ce_delta_vs_best_five_db"
            ]
            comparisons.append(comparison)
    if comparisons:
        _write_rows_atomic(comparisons, final_dir / "baseline_comparison.csv")
    dmrs_spacing = int(manifest["dmrs_spacing_subcarriers"])
    _plot_scene(
        summaries,
        rows_by_target,
        final_dir / "bler_curves.png",
        dmrs_spacing,
    )
    _plot_prescan(
        _read_stage_rows(scene_root / "prescan"),
        scenario_id,
        final_dir / "prescan_all_curves.png",
        dmrs_spacing,
    )
    schema = str(manifest.get("schema", ""))
    if schema.startswith("plan027-e6-"):
        experiment_section = "E6"
    elif schema.startswith("plan027-e5-"):
        experiment_section = "E5"
    else:
        experiment_section = "E4"
    report_lines = [
        f"# Plan-027 {experiment_section} 场景 {scenario_id} "
        f"comb-{dmrs_spacing}",
        "",
        f"- manifest SHA-256：`{digest}`",
        f"- DMRS：comb {dmrs_spacing}。",
        f"- 物理曲线数：{len(candidates)}。",
        "- 状态：场景计算完成，等待研究者确认；确认前不得运行下一场景。",
        "",
        "## 目标 SNR",
        "",
        "| 目标 | candidate | family | ideal-outage SNR | BLER SNR | BLER-outage | 95%区间 | 目标带错误数 | 充分 |",
        "|---|---|---|---:|---:|---:|---|---:|---|",
    ]
    for row in summaries:
        report_lines.append(
            f"| {row['target']} | {row['candidate_id']} | {row['family']} | "
            f"{float(row['outage_target_snr_db']):.3f} | "
            f"{float(row['target_snr_db']):.3f} | "
            f"{float(row['target_snr_db']) - float(row['outage_target_snr_db']):+.3f} | "
            f"[{float(row['ci95_low_db']):.3f},"
            f"{float(row['ci95_high_db']):.3f}] | "
            f"{row['errors_in_half_to_twice_target_band']} | "
            f"{'是' if row['sufficient_error_count'] else '否'} |"
        )
    summary_by_candidate_target = {
        (str(row["candidate_id"]), str(row["target"])): row for row in summaries
    }
    report_lines.extend(
        [
            "",
            "## 各方案 estimated-CSI BLER SNR 总览",
            "",
            "数值为正式拟合的目标SNR点估计；未形成正式双侧 bracket 的目标标为“未闭合”，"
            "不作外推。95%区间见上一表。",
            "",
            "| candidate | family | 10% BLER SNR (dB) | 1% BLER SNR (dB) |",
            "|---|---|---:|---:|",
        ]
    )
    for candidate in candidates:
        candidate_id = str(candidate["candidate_id"])
        family = str(candidate["family"])
        target_10 = summary_by_candidate_target.get((candidate_id, "10pct"))
        target_1 = summary_by_candidate_target.get((candidate_id, "1pct"))
        snr_10 = (
            f"{float(target_10['target_snr_db']):.3f}"
            if target_10 is not None
            else "未闭合"
        )
        snr_1 = (
            f"{float(target_1['target_snr_db']):.3f}"
            if target_1 is not None
            else "未闭合"
        )
        report_lines.append(
            f"| {candidate_id} | {family} | {snr_10} | {snr_1} |"
        )
    report_lines.extend(
        [
            "",
            "## 相对各条基线及全部基线最优",
            "",
            "BLER增益 = 基线目标SNR - 候选目标SNR，正值表示候选所需SNR更低；"
            "全部基线最优按同一场景、同一BLER目标下已闭合基线中目标SNR最低者确定；"
            "三者最优字段仍保留原 B0/AP_RMS/AP_TEPS 口径。",
            "",
            "| 目标 | candidate | 目标SNR | vs B0 | vs AP_RMS | vs AP_TEPS | vs TU/Nt | vs TU/(Nt-1) | vs Talias/Nt | vs 全部基线最优 |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in comparisons:
        report_lines.append(
            f"| {row['target']} | {row['candidate_id']} | "
            f"{float(row['candidate_target_snr_db']):.3f} | "
            f"{_format_signed(row['gain_vs_b0_db'])} | "
            f"{_format_signed(row['gain_vs_ap_rms_db'])} | "
            f"{_format_signed(row['gain_vs_ap_teps_db'])} | "
            f"{_format_signed(row['gain_vs_ap_tu_nt_db'])} | "
            f"{_format_signed(row['gain_vs_ap_tu_ntm1_db'])} | "
            f"{_format_signed(row['gain_vs_ap_talias_nt_db'])} | "
            f"{_format_signed(row['gain_vs_best_all_db'])} |"
        )
    report_lines.extend(
        [
            "",
            "## 各自工作点的 CE NMSE",
            "",
            "CE差值 = 候选在自身BLER目标SNR处的CE NMSE - 基线在自身BLER目标SNR处"
            "的CE NMSE；负值表示候选CE更低。CE NMSE由相邻正式SNR点在线性dB尺度内插。",
            "",
            "| 目标 | candidate | own SNR | own CE(dB) | ΔB0 | ΔAP_RMS | ΔAP_TEPS | ΔTU/Nt | ΔTU/(Nt-1) | ΔTalias/Nt | Δ全部基线最优 |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in comparisons:
        report_lines.append(
            f"| {row['target']} | {row['candidate_id']} | "
            f"{float(row['candidate_target_snr_db']):.3f} | "
            f"{float(row['candidate_target_ce_nmse_db']):.3f} | "
            f"{_format_signed(row['ce_delta_vs_b0_db'])} | "
            f"{_format_signed(row['ce_delta_vs_ap_rms_db'])} | "
            f"{_format_signed(row['ce_delta_vs_ap_teps_db'])} | "
            f"{_format_signed(row['ce_delta_vs_ap_tu_nt_db'])} | "
            f"{_format_signed(row['ce_delta_vs_ap_tu_ntm1_db'])} | "
            f"{_format_signed(row['ce_delta_vs_ap_talias_nt_db'])} | "
            f"{_format_signed(row['ce_delta_vs_best_all_db'])} |"
        )
    report_lines.extend(
        [
            "",
            "## 未闭合或拟合失败",
            "",
        ]
    )
    if failures:
        for failure in failures:
            report_lines.append(
                f"- `{failure['target']}` `{failure['candidate_id']}`："
                f"{failure['reason']}"
            )
    else:
        report_lines.append("- 无。")
    report_lines.extend(
        [
            "",
            "## BLER 曲线",
            "",
            "粗扫与细扫对同一 candidate 使用同一颜色、线型和 marker。"
            "细扫图只绘制正式闭合并成功拟合的目标。",
            "",
            "### 自适应粗扫",
            "",
            f"![{scenario_id} estimated-CSI BLER 自适应粗扫](prescan_all_curves.png)",
            "",
            "### 10%/1% 正式细扫",
            "",
            f"![{scenario_id} estimated-CSI BLER 正式细扫](bler_curves.png)",
            "",
            "机器可读比较：`final/baseline_comparison.csv`；逐点数据分别位于 "
            "`refine_10pct/bler_points.csv` 和 `refine_1pct/bler_points.csv`；"
            f"全{len(candidates)}条物理曲线的预扫图为 `final/prescan_all_curves.png`，"
            "正式闭合目标细扫图为 `final/bler_curves.png`。",
            "",
        ]
    )
    report_path = final_dir / "scenario_report.md"
    report_path.write_text(
        "\n".join(report_lines),
        encoding="utf-8",
    )
    result_hash = hashlib.sha256(report_path.read_bytes()).hexdigest()
    (final_dir / "scenario_result.sha256").write_text(
        f"{result_hash}  scenario_report.md\n",
        encoding="utf-8",
    )
    (final_dir / "scene_status.json").write_text(
        json.dumps(
            {
                "scenario_id": scenario_id,
                "status": "awaiting_researcher_confirmation",
                "manifest_sha256": digest,
                "scenario_report_sha256": result_hash,
                "fit_count": len(summaries),
                "failure_count": len(failures),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage",
        choices=("smoke", "prescan", "refine-10pct", "refine-1pct", "analyze"),
        required=True,
    )
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--scenario", choices=SCENE_ORDER, required=True)
    parser.add_argument("--approved-manifest", type=Path, required=True)
    parser.add_argument("--approved-sha256", required=True)
    parser.add_argument("--approval-receipt", type=Path, required=True)
    parser.add_argument("--prior-scene-approval", type=Path)
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    args = parser.parse_args()
    manifest = validate_approved_manifest(
        args.approved_manifest,
        args.approved_sha256,
        args.approval_receipt,
    )
    validate_prior_scene_approval(
        manifest,
        args.scenario,
        args.prior_scene_approval,
    )
    output_subdir = str(manifest.get("output_subdir", "e4_link"))
    scene_root = (
        args.output_root / str(args.run_id) / output_subdir / str(args.scenario)
    )
    if args.stage == "smoke":
        run_smoke(
            manifest,
            args.approved_sha256,
            args.scenario,
            int(args.batch_size),
            scene_root,
        )
    elif args.stage == "prescan":
        run_prescan(
            manifest,
            args.approved_sha256,
            args.scenario,
            int(args.batch_size),
            scene_root,
        )
    elif args.stage == "refine-10pct":
        run_refine(
            manifest,
            args.approved_sha256,
            args.scenario,
            "10pct",
            int(args.batch_size),
            scene_root,
        )
    elif args.stage == "refine-1pct":
        run_refine(
            manifest,
            args.approved_sha256,
            args.scenario,
            "1pct",
            int(args.batch_size),
            scene_root,
        )
    else:
        analyze_scene(
            manifest,
            args.approved_sha256,
            args.scenario,
            scene_root,
        )


if __name__ == "__main__":
    main()
