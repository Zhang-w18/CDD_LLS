"""Configuration-driven, resumable BLER-curve runner for static TDL-A studies.

The runner can start a new curve set or append non-overlapping trial intervals to
existing point-level CSV data.  Estimated-CSI and ideal-CSI receivers use the
same absolute-trial seed derivation, so samples remain paired wherever their
trial intervals overlap.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

try:
    import yaml
except ModuleNotFoundError as exc:  # pragma: no cover
    raise RuntimeError("PyYAML is required by tools/run_bler_curves.py") from exc

from cdd_lls.core.mcs import build_tb_layout, get_mcs
from cdd_lls.phy.channel_tdl import generate_sionna_tdl_channel
from cdd_lls.phy.estimators import build_frequency_rmmse_filter
from cdd_lls.phy.ldpc import SionnaLDPCAdapter
from cdd_lls.phy.precoding import equivalent_channel
from cdd_lls.phy.qam import qam_demapper_maxlog, qam_modulate
from cdd_lls.phy.resource_grid import local_indices_for_subcarriers
from tools import run_plan027_bler as bler027
from tools import run_plan028_csi_curves as curves028
from tools.run_plan025_delay_matched_tdl import (
    read_csv_rows,
    stable_seed,
    wilson,
    write_csv_rows,
)
from tools.run_v_design_piecewise_tradeoff import decode_same_tb_batch


SCHEMA = "bler-curve-runner-v1"
RECEIVERS = ("estimated", "ideal")


def _resolve(path: str | Path, base: Path = ROOT) -> Path:
    value = Path(path)
    return value if value.is_absolute() else base / value


def _repo_relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return str(path.resolve())


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_rows_atomic(rows: Sequence[dict], path: Path) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    write_csv_rows(rows, temporary)
    temporary.replace(path)


def load_runner_config(path: Path) -> dict:
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if payload.get("schema") != SCHEMA:
        raise ValueError(f"config.schema must be {SCHEMA!r}")
    if not isinstance(payload.get("scenes"), list) or not payload["scenes"]:
        raise ValueError("config.scenes must be a non-empty list")
    batch_size = int(payload.get("batch_size", 100))
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    payload["batch_size"] = batch_size
    payload["seed"] = int(payload.get("seed", bler027.SEED))
    payload["max_passes"] = int(payload.get("max_passes", 3))
    payload["output_dir"] = str(_resolve(payload["output_dir"]))
    payload["config_path"] = str(path.resolve())
    return payload


def _load_manifest(scene: dict) -> tuple[dict, str]:
    bundle = scene.get("manifest") or {}
    path = _resolve(bundle["path"])
    digest = _sha256(path)
    if bundle.get("sha256_path"):
        recorded = _resolve(bundle["sha256_path"]).read_text(encoding="utf-8").split()[0]
        if recorded.lower() != digest.lower():
            raise RuntimeError(f"Manifest SHA-256 mismatch: {path}")
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if bundle.get("approval_path"):
        approval_path = _resolve(bundle["approval_path"])
        approval = json.loads(approval_path.read_text(encoding="utf-8"))
        if approval.get("approved") is not True:
            raise RuntimeError(f"Manifest approval is not affirmative: {approval_path}")
        if str(approval.get("manifest_sha256", "")).lower() != digest.lower():
            raise RuntimeError(f"Manifest approval hash mismatch: {approval_path}")
    physical = manifest.get("physical_definition", {})
    if int(physical.get("phase_denominator", -1)) != 576:
        raise RuntimeError("BLER runner requires phase_denominator=576")
    if physical.get("delay_input_field") != "delay_grid_coordinates":
        raise RuntimeError("BLER runner requires delay_grid_coordinates")
    if int(manifest.get("dmrs_spacing_subcarriers", 0)) <= 0:
        raise RuntimeError("Manifest DMRS spacing must be positive")
    scenario_id = str(scene["scenario_id"])
    try:
        float(scenario_id.removeprefix("A"))
    except ValueError as exc:
        raise ValueError("scenario_id must encode TDL-A delay spread as A<number>") from exc
    candidates = bler027.scenario_candidates(manifest, scenario_id)
    requested = scene.get("candidate_ids", "all")
    if requested != "all":
        wanted = {str(value) for value in requested}
        actual = {str(row["candidate_id"]) for row in candidates}
        if not wanted <= actual:
            raise ValueError(f"Unknown candidate_ids: {sorted(wanted-actual)}")
    return manifest, digest


def _selected_candidates(scene: dict, manifest: dict) -> list[dict]:
    candidates = bler027.scenario_candidates(manifest, str(scene["scenario_id"]))
    requested = scene.get("candidate_ids", "all")
    if requested == "all":
        return candidates
    wanted = {str(value) for value in requested}
    return [row for row in candidates if str(row["candidate_id"]) in wanted]


def _base_rows(scene: dict, receiver: str, candidates: Sequence[dict]) -> list[dict]:
    receiver_cfg = scene["receivers"][receiver]
    base_csv = receiver_cfg.get("base_csv")
    if base_csv:
        rows = [dict(row) for row in read_csv_rows(_resolve(base_csv))]
        wanted = {str(row["candidate_id"]) for row in candidates}
        rows = [row for row in rows if str(row["candidate_id"]) in wanted]
        if not rows:
            raise RuntimeError(f"No selected rows in base_csv={base_csv}")
        return rows
    snrs = [float(value) for value in receiver_cfg.get("snr_db", [])]
    if not snrs:
        raise ValueError(f"{receiver} requires base_csv or snr_db")
    rows = []
    for candidate in candidates:
        for snr_db in snrs:
            rows.append(
                {
                    "scenario_id": scene["scenario_id"],
                    "candidate_id": candidate["candidate_id"],
                    "family": candidate["family"],
                    "snr_db": snr_db,
                    "trials": 0,
                    "tb_errors": 0,
                    "bler": 0.0,
                }
            )
    return rows


def _scene_output(config: dict, scene: dict) -> Path:
    return Path(config["output_dir"]) / str(scene["scenario_id"]).lower()


def _supplement_path(config: dict, scene: dict, receiver: str) -> Path:
    return _scene_output(config, scene) / receiver / "supplemental_points.csv"


def _flag_path(
    config: dict,
    scene: dict,
    receiver: str,
    snr_db: float,
    candidate_id: str,
    trial_start: int,
    trial_end: int,
) -> Path:
    snr_key = f"{float(snr_db):+.2f}".replace("+", "p").replace("-", "m").replace(".", "p")
    safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", candidate_id)
    name = f"{safe_id}_t{trial_start:06d}_{trial_end:06d}.npy"
    return _scene_output(config, scene) / receiver / "error_flags" / snr_key / name


def _read_supplements(config: dict, scene: dict, receiver: str) -> list[dict]:
    path = _supplement_path(config, scene, receiver)
    return [dict(row) for row in read_csv_rows(path)] if path.exists() else []


def _validate_intervals(base_rows: Sequence[dict], supplement_rows: Sequence[dict]) -> None:
    base = {
        (str(row["candidate_id"]), float(row["snr_db"])): int(row["trials"])
        for row in base_rows
    }
    grouped: dict[tuple[str, float], list[tuple[int, int]]] = defaultdict(list)
    for row in supplement_rows:
        key = (str(row["candidate_id"]), float(row["snr_db"]))
        grouped[key].append((int(row["trial_start"]), int(row["trial_end"])))
    for key, intervals in grouped.items():
        expected = base.get(key, 0) + 1
        for start, end in sorted(intervals):
            if start != expected or end < start:
                raise RuntimeError(
                    f"Non-contiguous or overlapping supplemental interval for {key}: "
                    f"expected start {expected}, got [{start}, {end}]"
                )
            expected = end + 1


def merge_rows(base_rows: Sequence[dict], supplement_rows: Sequence[dict]) -> list[dict]:
    _validate_intervals(base_rows, supplement_rows)
    additions: dict[tuple[str, float], list[dict]] = defaultdict(list)
    for row in supplement_rows:
        additions[(str(row["candidate_id"]), float(row["snr_db"]))].append(row)
    merged = []
    for source in base_rows:
        row = dict(source)
        key = (str(row["candidate_id"]), float(row["snr_db"]))
        extra = additions.get(key, [])
        base_trials = int(row["trials"])
        base_errors = int(row["tb_errors"])
        trials = base_trials + sum(int(item["trials"]) for item in extra)
        errors = base_errors + sum(int(item["tb_errors"]) for item in extra)
        low, high = wilson(errors, trials)
        row.update(
            {
                "base_trials": base_trials,
                "supplemental_trials": trials - base_trials,
                "trials": trials,
                "tb_errors": errors,
                "bler": errors / trials,
                "bler_wilson95_lo": low,
                "bler_wilson95_hi": high,
            }
        )
        if "ce_nmse_mean" in row:
            total_nmse = float(row["ce_nmse_mean"]) * base_trials
            total_nmse += sum(float(item.get("ce_nmse_sum", 0.0)) for item in extra)
            mean_nmse = total_nmse / trials
            row["ce_nmse_mean"] = mean_nmse
            row["ce_nmse_mean_db"] = 10.0 * math.log10(mean_nmse)
        merged.append(row)
    merged.sort(key=lambda row: (str(row["candidate_id"]), float(row["snr_db"])))
    return merged


def _plot_keys(rows: Sequence[dict], stop_below_bler: float | None) -> set[tuple[str, float]]:
    if stop_below_bler is None:
        return {(str(row["candidate_id"]), float(row["snr_db"])) for row in rows}
    output = set()
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[str(row["candidate_id"])].append(row)
    for candidate_id, values in grouped.items():
        for row in sorted(values, key=lambda item: float(item["snr_db"])):
            output.add((candidate_id, float(row["snr_db"])))
            if float(row["bler"]) < float(stop_below_bler):
                break
    return output


def _candidate_policy(policy: dict, candidate_id: str) -> dict:
    """Return the receiver policy after applying an optional candidate override."""
    effective = {
        key: value for key, value in policy.items() if key != "candidate_overrides"
    }
    overrides = policy.get("candidate_overrides", {})
    if overrides:
        override = overrides.get(candidate_id, {})
        if not isinstance(override, dict):
            raise ValueError(f"candidate_overrides[{candidate_id!r}] must be a mapping")
        effective.update(override)
    return effective


def _policy_plot_keys(rows: Sequence[dict], policy: dict) -> set[tuple[str, float]]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[str(row["candidate_id"])].append(row)
    output = set()
    for candidate_id, values in grouped.items():
        effective = _candidate_policy(policy, candidate_id)
        stop = effective.get("stop_below_bler")
        output.update(_plot_keys(values, None if stop is None else float(stop)))
    return output


def build_tasks(
    rows: Sequence[dict], policy: dict, batch_size: int
) -> list[dict]:
    selected = _policy_plot_keys(rows, policy)
    tasks = []
    for row in rows:
        key = (str(row["candidate_id"]), float(row["snr_db"]))
        if key not in selected:
            continue
        effective = _candidate_policy(policy, key[0])
        target_errors = int(effective.get("target_errors", 0))
        maximum = int(effective["max_total_trials"])
        minimum = int(effective.get("min_total_trials", 0))
        fixed = effective.get("fixed_total_trials")
        current = int(row["trials"])
        errors = int(row["tb_errors"])
        if current % batch_size:
            raise ValueError(f"Existing trials must align to batch_size for {key}")
        if fixed is not None:
            target = int(fixed)
        elif target_errors > 0 and errors >= target_errors:
            target = current
        elif errors > 0:
            target = math.ceil((target_errors * current / errors) / batch_size) * batch_size
        else:
            target = maximum
        target = max(current, minimum, target)
        target = min(maximum, math.ceil(target / batch_size) * batch_size)
        if target > current:
            tasks.append(
                {
                    "candidate_id": key[0],
                    "snr_db": key[1],
                    "trial_start": current + 1,
                    "trial_end": target,
                    "trials": target - current,
                }
            )
    return tasks


def _decode_flags(
    adapter: SionnaLDPCAdapter,
    equalized: np.ndarray,
    effective_noise: np.ndarray,
    qm: int,
    payload: Sequence[np.ndarray],
) -> np.ndarray:
    llrs = [
        qam_demapper_maxlog(equalized[index], effective_noise[index], qm)
        for index in range(len(equalized))
    ]
    decoded = decode_same_tb_batch(adapter, llrs, payload)
    return np.asarray([not item.tb_success for item in decoded], dtype=bool)


def _run_task_groups(
    config: dict,
    scene: dict,
    manifest: dict,
    candidates: Sequence[dict],
    receiver: str,
    tasks: Sequence[dict],
) -> None:
    if not tasks:
        return
    scenario_id = str(scene["scenario_id"])
    by_id = {str(row["candidate_id"]): row for row in candidates}
    channel, grid, precoders, covariances = bler027.build_scene(
        manifest, candidates, scenario_id
    )
    mcs = get_mcs("nr_256qam", 8, None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(tb.cb_k_values, tb.cb_e_values, num_iter=8, llr_clip=50.0)
    pilot_local = local_indices_for_subcarriers(grid, grid.pilot_subcarriers)
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    stored = _read_supplements(config, scene, receiver)
    groups: dict[float, list[dict]] = defaultdict(list)
    for task in tasks:
        groups[float(task["snr_db"])].append(dict(task))
    for snr_db, snr_tasks in sorted(groups.items()):
        task_by_id = {str(task["candidate_id"]): task for task in snr_tasks}
        candidate_ids = sorted(task_by_id)
        flags = {
            candidate_id: np.zeros(int(task_by_id[candidate_id]["trials"]), dtype=bool)
            for candidate_id in candidate_ids
        }
        nmse_sum = {candidate_id: 0.0 for candidate_id in candidate_ids}
        first_trial = min(int(task["trial_start"]) for task in snr_tasks)
        last_trial = max(int(task["trial_end"]) for task in snr_tasks)
        snr_linear = 10.0 ** (snr_db / 10.0)
        noise_variance = 8.0 / snr_linear
        ls_noise_variance = noise_variance / 2.0
        filters = (
            {
                candidate_id: build_frequency_rmmse_filter(
                    covariances[candidate_id], pilot_local, ls_noise_variance, diagonal_loading=1e-10
                )
                for candidate_id in candidate_ids
            }
            if receiver == "estimated"
            else {}
        )
        started = time.time()
        for absolute_start in range(first_trial, last_trial + 1, int(config["batch_size"])):
            active = [
                candidate_id
                for candidate_id in candidate_ids
                if int(task_by_id[candidate_id]["trial_start"]) <= absolute_start
                <= int(task_by_id[candidate_id]["trial_end"])
            ]
            if not active:
                continue
            current = min(int(config["batch_size"]), last_trial - absolute_start + 1)
            payload_rng = np.random.default_rng(
                stable_seed(config["seed"], scenario_id, snr_db, absolute_start, "payload")
            )
            payload = [payload_rng.integers(0, 2, size=int(k), dtype=np.int8) for k in tb.cb_k_values]
            symbols = qam_modulate(np.concatenate(adapter.encode(payload)), int(mcs.qm))
            realization = generate_sionna_tdl_channel(
                grid,
                channel,
                n_tx=8,
                n_rx=1,
                batch_size=current,
                seed=stable_seed(config["seed"], scenario_id, snr_db, absolute_start, "channel"),
            )
            noise_rng = np.random.default_rng(
                stable_seed(config["seed"], scenario_id, snr_db, absolute_start, "noise")
            )
            averaged_ls_noise = math.sqrt(ls_noise_variance / 2.0) * (
                noise_rng.normal(size=(current, 1, grid.pilot_count))
                + 1j * noise_rng.normal(size=(current, 1, grid.pilot_count))
            )
            data_noise = math.sqrt(noise_variance / 2.0) * (
                noise_rng.normal(size=(current, 1, grid.n_data_re))
                + 1j * noise_rng.normal(size=(current, 1, grid.n_data_re))
            )
            decode_batches = []
            for candidate_id in active:
                task = task_by_id[candidate_id]
                used = min(current, int(task["trial_end"]) - absolute_start + 1)
                offset = absolute_start - int(task["trial_start"])
                effective = equivalent_channel(realization.H, precoders[candidate_id].C)[:used, 0:1]
                true_data = effective[:, :, grid.data_symbol_indices, data_local]
                received = true_data * symbols[None, None, :] + data_noise[:used]
                if receiver == "estimated":
                    pilot0 = effective[:, :, int(grid.pilot_symbol_indices[0]), pilot_local]
                    pilot1 = effective[:, :, int(grid.pilot_symbol_indices[-1]), pilot_local]
                    estimate_full = filters[candidate_id].estimate_full_band(
                        0.5 * (pilot0 + pilot1) + averaged_ls_noise[:used]
                    )
                    estimate_data = estimate_full[:, :, data_local]
                    ratios = np.sum(np.abs(estimate_data - true_data) ** 2, axis=(1, 2)) / np.maximum(
                        np.sum(np.abs(true_data) ** 2, axis=(1, 2)), 1e-30
                    )
                    nmse_sum[candidate_id] += float(np.sum(ratios))
                    denominator = np.maximum(np.sum(np.abs(estimate_data) ** 2, axis=1), 1e-10)
                    equalized = np.sum(np.conj(estimate_data) * received, axis=1) / denominator
                    effective_noise = noise_variance / denominator
                else:
                    equalized, effective_noise = curves028.ideal_csi_equalize(
                        true_data, received, noise_variance
                    )
                decode_batches.append(
                    (candidate_id, offset, used, equalized, effective_noise)
                )
            combined_equalized = np.concatenate(
                [item[3] for item in decode_batches], axis=0
            )
            combined_noise = np.concatenate(
                [item[4] for item in decode_batches], axis=0
            )
            combined_flags = _decode_flags(
                adapter,
                combined_equalized,
                combined_noise,
                int(mcs.qm),
                payload,
            )
            cursor = 0
            for candidate_id, offset, used, _, _ in decode_batches:
                flags[candidate_id][offset : offset + used] = combined_flags[
                    cursor : cursor + used
                ]
                cursor += used
            completed = absolute_start + current - 1
            if absolute_start == first_trial or completed % 500 < current or completed == last_trial:
                print(
                    f"[BLER {scenario_id} {receiver}] snr={snr_db:g} "
                    f"absolute_trials={completed}/{last_trial} active={len(active)} "
                    f"elapsed={time.time()-started:.1f}s",
                    flush=True,
                )
        for candidate_id in candidate_ids:
            task = task_by_id[candidate_id]
            trial_start = int(task["trial_start"])
            trial_end = int(task["trial_end"])
            trial_count = int(task["trials"])
            flag_path = _flag_path(
                config, scene, receiver, snr_db, candidate_id, trial_start, trial_end
            )
            flag_path.parent.mkdir(parents=True, exist_ok=True)
            np.save(flag_path, flags[candidate_id])
            errors = int(np.sum(flags[candidate_id]))
            low, high = wilson(errors, trial_count)
            stored.append(
                {
                    "scenario_id": scenario_id,
                    "candidate_id": candidate_id,
                    "family": by_id[candidate_id]["family"],
                    "receiver": receiver,
                    "snr_db": snr_db,
                    "trial_start": trial_start,
                    "trial_end": trial_end,
                    "trials": trial_count,
                    "tb_errors": errors,
                    "bler": errors / trial_count,
                    "bler_wilson95_lo": low,
                    "bler_wilson95_hi": high,
                    "ce_nmse_sum": nmse_sum[candidate_id] if receiver == "estimated" else "",
                    "seed": int(config["seed"]),
                    "error_flags": _repo_relative(flag_path),
                }
            )
        stored.sort(
            key=lambda row: (
                str(row["candidate_id"]),
                float(row["snr_db"]),
                int(row["trial_start"]),
            )
        )
        _write_rows_atomic(stored, _supplement_path(config, scene, receiver))


def _write_styles_and_mapping(output: Path, candidates: Sequence[dict]) -> tuple[dict, dict, dict]:
    ordered = sorted(candidates, key=curves028._candidate_sort_key)
    styles = {}
    candidate_labels = {}
    abbreviations = {}
    for index, candidate in enumerate(ordered, start=1):
        candidate_id = str(candidate["candidate_id"])
        style = bler027.curve_style(candidate_id, str(candidate["family"]))
        styles[candidate_id] = {
            "color": style["color"],
            "linestyle": style["linestyle"],
            "marker": style["marker"],
        }
        candidate_labels[candidate_id] = f"candidate {index:02d}"
        abbreviations[candidate_id] = bler027.curve_label(candidate_id)
    output.mkdir(parents=True, exist_ok=True)
    (output / "curve_styles.json").write_text(
        json.dumps(styles, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return styles, candidate_labels, abbreviations


def merge_and_plot_scene(config: dict, scene: dict, manifest: dict, candidates: Sequence[dict]) -> None:
    output = _scene_output(config, scene) / "final"
    styles, candidate_labels, abbreviations = _write_styles_and_mapping(output, candidates)
    scenario_id = str(scene["scenario_id"])
    merged_by_receiver = {}
    for receiver in RECEIVERS:
        if receiver not in scene.get("receivers", {}):
            continue
        base = _base_rows(scene, receiver, candidates)
        merged = merge_rows(base, _read_supplements(config, scene, receiver))
        policy = scene["receivers"][receiver]
        included = _policy_plot_keys(merged, policy)
        for row in merged:
            row["plot_included"] = (
                str(row["candidate_id"]), float(row["snr_db"])
            ) in included
        merged_by_receiver[receiver] = merged
        _write_rows_atomic(merged, output / f"{receiver}_csi_bler_points.csv")
        plotted = [row for row in merged if bool(row["plot_included"])]
        for mode, mapping in (("candidate", candidate_labels), ("abbrev", abbreviations)):
            curves028._plot_metric(
                plotted,
                mapping,
                styles,
                "bler",
                f"{receiver.capitalize()}-CSI BLER",
                f"{scenario_id} comb-{manifest['dmrs_spacing_subcarriers']} {receiver}-CSI BLER",
                output / f"{scenario_id.lower()}_comb{manifest['dmrs_spacing_subcarriers']}_{mode}_{receiver}_csi_bler.png",
                True,
            )
    if "estimated" in merged_by_receiver:
        plotted = [row for row in merged_by_receiver["estimated"] if bool(row["plot_included"])]
        for mode, mapping in (("candidate", candidate_labels), ("abbrev", abbreviations)):
            curves028._plot_metric(
                plotted,
                mapping,
                styles,
                "ce_nmse_mean_db",
                "CE NMSE (dB)",
                f"{scenario_id} comb-{manifest['dmrs_spacing_subcarriers']} matched CE NMSE",
                output / f"{scenario_id.lower()}_comb{manifest['dmrs_spacing_subcarriers']}_{mode}_ce_nmse.png",
                False,
            )


def validate_config(config: dict) -> None:
    seen = set()
    for scene in config["scenes"]:
        scenario_id = str(scene["scenario_id"])
        if scenario_id in seen:
            raise ValueError(f"Duplicate scenario_id: {scenario_id}")
        seen.add(scenario_id)
        manifest, _ = _load_manifest(scene)
        candidates = _selected_candidates(scene, manifest)
        for receiver in scene.get("receivers", {}):
            if receiver not in RECEIVERS:
                raise ValueError(f"Unsupported receiver: {receiver}")
            base = _base_rows(scene, receiver, candidates)
            _validate_intervals(base, _read_supplements(config, scene, receiver))
            policy = scene["receivers"][receiver]
            if int(policy["max_total_trials"]) <= 0:
                raise ValueError("max_total_trials must be positive")
            overrides = policy.get("candidate_overrides", {})
            if not isinstance(overrides, dict):
                raise ValueError("candidate_overrides must be a mapping")
            selected_ids = {str(row["candidate_id"]) for row in candidates}
            unknown = set(overrides) - selected_ids
            if unknown:
                raise ValueError(f"Unknown candidate_overrides: {sorted(unknown)}")
            for candidate_id in overrides:
                effective = _candidate_policy(policy, candidate_id)
                if int(effective["max_total_trials"]) <= 0:
                    raise ValueError("max_total_trials must be positive")
        print(
            f"[validate] {scenario_id}: candidates={len(candidates)} "
            f"receivers={sorted(scene.get('receivers', {}))}",
            flush=True,
        )


def run_config(config: dict) -> None:
    for scene in config["scenes"]:
        manifest, digest = _load_manifest(scene)
        candidates = _selected_candidates(scene, manifest)
        scene_output = _scene_output(config, scene)
        scene_output.mkdir(parents=True, exist_ok=True)
        resolved = {
            "schema": SCHEMA,
            "scenario_id": scene["scenario_id"],
            "manifest_sha256": digest,
            "seed": config["seed"],
            "batch_size": config["batch_size"],
            "config_path": config["config_path"],
            "config_sha256": _sha256(Path(config["config_path"])),
            "receivers": scene.get("receivers", {}),
            "absolute_trial_seeded": True,
        }
        (scene_output / "resolved_run.json").write_text(
            json.dumps(resolved, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        for receiver in RECEIVERS:
            if receiver not in scene.get("receivers", {}):
                continue
            for _ in range(int(config["max_passes"])):
                base = _base_rows(scene, receiver, candidates)
                supplements = _read_supplements(config, scene, receiver)
                merged = merge_rows(base, supplements)
                tasks = build_tasks(merged, scene["receivers"][receiver], int(config["batch_size"]))
                if not tasks:
                    break
                _run_task_groups(config, scene, manifest, candidates, receiver, tasks)
        merge_and_plot_scene(config, scene, manifest, candidates)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--stage", choices=("validate", "run", "merge", "all"), default="all")
    args = parser.parse_args()
    config = load_runner_config(args.config.resolve())
    validate_config(config)
    if args.stage in ("run", "all"):
        run_config(config)
    elif args.stage == "merge":
        for scene in config["scenes"]:
            manifest, _ = _load_manifest(scene)
            merge_and_plot_scene(config, scene, manifest, _selected_candidates(scene, manifest))


if __name__ == "__main__":
    main()
