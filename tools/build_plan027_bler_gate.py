"""Build the immutable, deduplicated Plan-027 E4 candidate gate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.design import (
    canonical_delay_set,
    delay_indices_to_ns,
    fold_min_gap_indices,
    pair_sum_min_gap_indices,
    residues_and_lifts,
)
from tools.run_plan025_delay_matched_tdl import write_csv_rows
from tools.run_plan027_meff_design import (
    DEFAULT_OUTPUT_ROOT,
    DEFAULT_RUN_ID,
    GRID_NS,
    K,
    PILOT_PERIOD,
    SCS_HZ,
    decode,
    encode,
    read_csv,
    save_json,
)


SOURCE026 = ROOT / "outputs" / "experiment026_cdd_design" / "20260724_main"
LINK_SCENARIO_ORDER = ("A5", "A1", "A10", "A30", "A100")
LINK_BASELINE_FAMILIES = ("B0_QC", "AP_RMS_T1", "AP_TEPS_T1")
LINK_SYSTEM_BASELINE_FAMILIES = ("AP_TU_NT", "AP_TALIAS_NT")
LINK_SYSTEM_BASELINE_SCENARIOS = ("A5", "A10", "A30")
LINK_ALL_BASELINE_FAMILIES = (
    *LINK_BASELINE_FAMILIES,
    *LINK_SYSTEM_BASELINE_FAMILIES,
)
LINK_NEW_FAMILIES = ("S0_SIDON", "AP_T2_CTRL", "GEO_T1_CTRL", "MEFF_T2_CAND")
LINK_EXPECTED_COUNTS = {"A1": 10, "A5": 11, "A10": 12, "A30": 11, "A100": 9}


def positive_027(output: Path) -> list[dict[str, object]]:
    targets = read_csv(output / "e1_outage" / "e1_outage_targets.csv")
    triggers = {
        (row["scenario_id"], row["candidate_id"]): row
        for row in read_csv(output / "e3_ce_density" / "e3_density_trigger.csv")
    }
    ce_rows = read_csv(output / "e3_ce_density" / "e3_ce_at_targets.csv")
    ce_lookup = {}
    for row in ce_rows:
        if (
            row["dmrs_spacing_subcarriers"] == "24"
            and row["working_point"] == "common_baseline_target"
        ):
            ce_lookup[
                (
                    row["scenario_id"],
                    row["candidate_id"],
                    float(row["target_outage_probability"]),
                )
            ] = row
    by_id = {(row["scenario_id"], row["candidate_id"]): row for row in targets}
    output_rows = []
    for row in targets:
        if row["family"] not in ("AP_T2_CTRL", "GEO_T1_CTRL", "MEFF_T2_CAND"):
            continue
        gain10 = float(row["gain_vs_ap_t1_base_10_db"])
        gain1 = float(row["gain_vs_ap_t1_base_1_db"])
        if gain10 <= 0.0 and gain1 <= 0.0:
            continue
        baseline10 = by_id[(row["scenario_id"], row["ap_t1_base_10_candidate_id"])]
        baseline1 = by_id[(row["scenario_id"], row["ap_t1_base_1_candidate_id"])]
        indices = tuple(int(value) for value in decode(row["delay_indices"]))
        key = canonical_delay_set(indices, K)
        trigger = triggers[(row["scenario_id"], row["candidate_id"])]
        output_rows.append(
            {
                "source_experiment": "027",
                "scenario_id": row["scenario_id"],
                "profile": "A",
                "delay_spread_ns": float(
                    row["scenario_id"].removeprefix("A").replace("p", ".")
                ),
                "original_candidate_id": row["candidate_id"],
                "family": row["family"],
                "selection_rule": row["selection_rule"],
                "delay_indices": encode(indices),
                "delay_grid_coordinates": row["delay_grid_coordinates"],
                "delay_ns": row["delay_ns"],
                "residues": row["residues"],
                "lifts": row["lifts"],
                "canonical_key": encode(key),
                "baseline_10_candidate_id": baseline10["candidate_id"],
                "baseline_10_delay_grid_coordinates": baseline10[
                    "delay_grid_coordinates"
                ],
                "baseline_1_candidate_id": baseline1["candidate_id"],
                "baseline_1_delay_grid_coordinates": baseline1[
                    "delay_grid_coordinates"
                ],
                "outage10_snr_db": float(row["outage10_snr_db"]),
                "outage1_snr_db": float(row["outage1_snr_db"]),
                "gain10_db": gain10,
                "gain10_ci95_low_db": float(
                    row["gain_vs_ap_t1_base_10_ci95_low_db"]
                ),
                "gain10_ci95_high_db": float(
                    row["gain_vs_ap_t1_base_10_ci95_high_db"]
                ),
                "gain1_db": gain1,
                "gain1_ci95_low_db": float(
                    row["gain_vs_ap_t1_base_1_ci95_low_db"]
                ),
                "gain1_ci95_high_db": float(
                    row["gain_vs_ap_t1_base_1_ci95_high_db"]
                ),
                "ce_comb24_10_nmse_db": float(
                    ce_lookup[(row["scenario_id"], row["candidate_id"], 0.10)][
                        "ce_nmse_db"
                    ]
                ),
                "ce_comb24_1_nmse_db": float(
                    ce_lookup[(row["scenario_id"], row["candidate_id"], 0.01)][
                        "ce_nmse_db"
                    ]
                ),
                "comb12_triggered": trigger["comb12_triggered"] == "True",
                "comb6_triggered": trigger["comb6_triggered"] == "True",
                "strict_sidon": pair_sum_min_gap_indices(indices, K) >= 1,
                "hard_thick_sidon": row["hard_thick_sidon"] == "True",
                "pair_sum_gap_indices": int(row["pair_sum_gap_indices"]),
                "fold_gap_indices": int(row["fold_gap_indices"]),
                "m2_eff": float(row["m2_eff"]),
                "m4_eff": float(row["m4_eff"]),
                "provenance": [f"027:{row['candidate_id']}:{row['family']}"],
                "recommended_dmrs_spacings": (
                    [24, 12, 6]
                    if trigger["comb6_triggered"] == "True"
                    else [24, 12]
                    if trigger["comb12_triggered"] == "True"
                    else [24]
                ),
            }
        )
    return output_rows


def baseline_026(rows: list[dict[str, str]], scenario: str, stage: str) -> dict[str, str]:
    selected = [row for row in rows if row["scenario_id"] == scenario]
    if stage == "e1":
        baselines = [row for row in selected if "B1" in row["family"]]
    else:
        baselines = [row for row in selected if row["family"] == "B1_T1"]
    if not baselines:
        raise RuntimeError(f"Could not locate the 026 baseline for {scenario}.")
    return min(baselines, key=lambda row: float(row["outage10_snr_db"]))


def positive_026() -> list[dict[str, object]]:
    sources = (
        ("e1", SOURCE026 / "e1_outage" / "outage_targets.csv"),
        ("e2", SOURCE026 / "e2_outage" / "outage_targets.csv"),
    )
    metric_paths = {
        "e1": SOURCE026 / "e1_search" / "e1_metrics.csv",
        "e2": SOURCE026 / "e2_thick_sidon" / "e2_metrics.csv",
    }
    output_rows = []
    for stage, path in sources:
        rows = read_csv(path)
        metrics = {
            (row.get("scenario_id", "E1-A5"), row["candidate_id"]): row
            for row in read_csv(metric_paths[stage])
        }
        for scenario in sorted({row["scenario_id"] for row in rows}):
            baseline = baseline_026(rows, scenario, stage)
            for row in [value for value in rows if value["scenario_id"] == scenario]:
                if (
                    "B0" in row["family"]
                    or "B1" in row["family"]
                    or row["family"] == "S0"
                ):
                    continue
                gain10 = float(baseline["outage10_snr_db"]) - float(
                    row["outage10_snr_db"]
                )
                gain1 = float(baseline["outage1_snr_db"]) - float(
                    row["outage1_snr_db"]
                )
                if gain10 <= 0.0 and gain1 <= 0.0:
                    continue
                indices = canonical_delay_set(
                    [int(value) for value in json.loads(row["delay_indices"])],
                    K,
                )
                residues, lifts = residues_and_lifts(indices, PILOT_PERIOD)
                metric = metrics.get((scenario, row["candidate_id"]), {})
                profile_spread = scenario.split("-")[-1]
                profile = profile_spread[0]
                spread = float(profile_spread[1:].replace("p", "."))
                output_rows.append(
                    {
                        "source_experiment": "026",
                        "scenario_id": scenario,
                        "profile": profile,
                        "delay_spread_ns": spread,
                        "original_candidate_id": row["candidate_id"],
                        "family": row["family"],
                        "selection_rule": "result026_frozen_positive_outage",
                        "delay_indices": encode(indices),
                        "delay_grid_coordinates": encode(indices),
                        "delay_ns": encode(delay_indices_to_ns(indices, K, SCS_HZ)),
                        "residues": encode(residues),
                        "lifts": encode(lifts),
                        "canonical_key": encode(indices),
                        "baseline_10_candidate_id": baseline["candidate_id"],
                        "baseline_10_delay_grid_coordinates": baseline["delay_indices"],
                        "baseline_1_candidate_id": baseline["candidate_id"],
                        "baseline_1_delay_grid_coordinates": baseline["delay_indices"],
                        "outage10_snr_db": float(row["outage10_snr_db"]),
                        "outage1_snr_db": float(row["outage1_snr_db"]),
                        "gain10_db": gain10,
                        "gain10_ci95_low_db": "",
                        "gain10_ci95_high_db": "",
                        "gain1_db": gain1,
                        "gain1_ci95_low_db": "",
                        "gain1_ci95_high_db": "",
                        "ce_comb24_10_nmse_db": metric.get("ce_nmse_db_16db", ""),
                        "ce_comb24_1_nmse_db": metric.get("ce_nmse_db_16db", ""),
                        "comb12_triggered": False,
                        "comb6_triggered": False,
                        "strict_sidon": pair_sum_min_gap_indices(indices, K) >= 1,
                        "hard_thick_sidon": metric.get("hard_thick_sidon", "False")
                        == "True",
                        "pair_sum_gap_indices": pair_sum_min_gap_indices(indices, K),
                        "fold_gap_indices": fold_min_gap_indices(indices, PILOT_PERIOD),
                        "m2_eff": metric.get("m2_eff_16db", ""),
                        "m4_eff": metric.get("m4_eff_16db", ""),
                        "provenance": [f"026:{stage}:{row['candidate_id']}:{row['family']}"],
                        "recommended_dmrs_spacings": [24],
                    }
                )
    return output_rows


def deduplicate(rows: Iterable[dict[str, object]]) -> list[dict[str, object]]:
    merged: dict[tuple[str, float, str], dict[str, object]] = {}
    for row in rows:
        identity = (
            str(row["profile"]),
            float(row["delay_spread_ns"]),
            str(row["canonical_key"]),
        )
        if identity not in merged:
            merged[identity] = dict(row)
            merged[identity]["duplicate_candidate_ids"] = []
            continue
        current = merged[identity]
        current["provenance"] = sorted(
            set(current["provenance"]) | set(row["provenance"])
        )
        current["duplicate_candidate_ids"].append(row["original_candidate_id"])
        current["family"] = "+".join(
            sorted(set(str(current["family"]).split("+")) | set(str(row["family"]).split("+")))
        )
        if row["source_experiment"] == "027" and current["source_experiment"] != "027":
            preserved_provenance = current["provenance"]
            duplicates = current["duplicate_candidate_ids"]
            merged[identity] = dict(row)
            merged[identity]["provenance"] = preserved_provenance
            merged[identity]["duplicate_candidate_ids"] = duplicates
    result = list(merged.values())
    result.sort(
        key=lambda row: (
            str(row["profile"]),
            float(row["delay_spread_ns"]),
            str(row["canonical_key"]),
        )
    )
    for index, row in enumerate(result):
        row["gate_id"] = f"E4_{index:04d}"
        row["duplicate_candidate_ids"] = sorted(set(row["duplicate_candidate_ids"]))
        row["provenance"] = sorted(set(row["provenance"]))
    return result


def _decode_optional_array(value: object) -> list[object] | None:
    text = str(value).strip()
    if not text:
        return None
    decoded = json.loads(text)
    if not isinstance(decoded, list):
        raise ValueError(f"Expected JSON array, got {type(decoded).__name__}.")
    return decoded


def envelope_link_candidates(
    targets: list[dict[str, str]],
) -> list[dict[str, object]]:
    """Select the union of 10%/1% family winners plus all fixed references."""
    output: list[dict[str, object]] = []
    target_specs = (("10pct", "outage10_snr_db"), ("1pct", "outage1_snr_db"))
    for scenario_id in LINK_SCENARIO_ORDER:
        scenario_rows = [row for row in targets if row["scenario_id"] == scenario_id]
        if not scenario_rows:
            raise RuntimeError(f"No E1 outage rows found for {scenario_id}.")
        selected: dict[str, dict[str, object]] = {}

        for family in LINK_BASELINE_FAMILIES:
            matches = [row for row in scenario_rows if row["family"] == family]
            if len(matches) != 1:
                raise RuntimeError(
                    f"Expected one {family} row for {scenario_id}, got {len(matches)}."
                )
            selected[matches[0]["candidate_id"]] = {
                "row": matches[0],
                "roles": ["baseline"],
                "outage_winner_for": [],
            }

        if scenario_id in LINK_SYSTEM_BASELINE_SCENARIOS:
            for family in LINK_SYSTEM_BASELINE_FAMILIES:
                matches = [row for row in scenario_rows if row["family"] == family]
                if len(matches) != 1:
                    raise RuntimeError(
                        f"Expected one {family} row for {scenario_id}, "
                        f"got {len(matches)}."
                    )
                selected[matches[0]["candidate_id"]] = {
                    "row": matches[0],
                    "roles": ["system_geometry_baseline"],
                    "outage_winner_for": [],
                }

        s0_matches = [row for row in scenario_rows if row["family"] == "S0_SIDON"]
        if len(s0_matches) != 1:
            raise RuntimeError(
                f"Expected one S0_SIDON row for {scenario_id}, got {len(s0_matches)}."
            )
        selected[s0_matches[0]["candidate_id"]] = {
            "row": s0_matches[0],
            "roles": ["historical_reference"],
            "outage_winner_for": ["10pct", "1pct"],
        }

        for family in LINK_NEW_FAMILIES[1:]:
            matches = [row for row in scenario_rows if row["family"] == family]
            if not matches:
                raise RuntimeError(f"No {family} rows found for {scenario_id}.")
            for target_label, field in target_specs:
                winner = min(matches, key=lambda row: float(row[field]))
                candidate_id = winner["candidate_id"]
                if candidate_id not in selected:
                    selected[candidate_id] = {
                        "row": winner,
                        "roles": ["outage_family_winner"],
                        "outage_winner_for": [],
                    }
                selected[candidate_id]["outage_winner_for"].append(target_label)

        for candidate_id, item in selected.items():
            row = item["row"]
            coordinates = _decode_optional_array(row["delay_grid_coordinates"])
            if coordinates is None or len(coordinates) != 8:
                raise RuntimeError(
                    f"{scenario_id}/{candidate_id} does not have eight delay coordinates."
                )
            delay_ns = _decode_optional_array(row["delay_ns"])
            indices = _decode_optional_array(row.get("delay_indices", ""))
            entry = {
                "curve_id": f"{scenario_id}:{candidate_id}",
                "scenario_id": scenario_id,
                "profile": "A",
                "delay_spread_ns": float(scenario_id.removeprefix("A")),
                "candidate_id": candidate_id,
                "family": row["family"],
                "roles": item["roles"],
                "outage_winner_for": sorted(set(item["outage_winner_for"])),
                "selection_rule": row["selection_rule"],
                "delay_grid_coordinates": [float(value) for value in coordinates],
                "delay_indices": (
                    [int(value) for value in indices] if indices is not None else None
                ),
                "delay_ns": [float(value) for value in delay_ns or []],
                "grid_aligned": row["grid_aligned"] == "True",
                "outage10_snr_db": float(row["outage10_snr_db"]),
                "outage1_snr_db": float(row["outage1_snr_db"]),
                "phase_denominator_active_subcarriers": K,
                "phase_reference": "first_active_subcarrier",
                "precoder_normalized": False,
                "dmrs_spacing_subcarriers": 24,
            }
            output.append(entry)

    output.sort(
        key=lambda row: (
            LINK_SCENARIO_ORDER.index(str(row["scenario_id"])),
            LINK_ALL_BASELINE_FAMILIES.index(str(row["family"]))
            if row["family"] in LINK_ALL_BASELINE_FAMILIES
            else len(LINK_ALL_BASELINE_FAMILIES)
            + LINK_NEW_FAMILIES.index(str(row["family"])),
            str(row["candidate_id"]),
        )
    )
    actual_counts = {
        scenario_id: sum(row["scenario_id"] == scenario_id for row in output)
        for scenario_id in LINK_SCENARIO_ORDER
    }
    if actual_counts != LINK_EXPECTED_COUNTS:
        raise RuntimeError(
            f"Envelope-link curve counts changed: {actual_counts}; "
            f"expected {LINK_EXPECTED_COUNTS}."
        )
    if len(output) != 53:
        raise RuntimeError(f"Expected 53 link curves, got {len(output)}.")
    for index, row in enumerate(output):
        row["gate_id"] = f"E4L_{index:04d}"
    return output


def link_markdown_table(
    rows: list[dict[str, object]],
    digest: str,
    manifest: dict[str, object],
) -> str:
    lines = [
        "# Plan-027 E4 53曲线链路执行清单（A10/A30系统基线增补）",
        "",
        f"- manifest SHA-256：`{digest}`",
        f"- 场景×物理候选数：{len(rows)}",
        f"- 场景顺序：`{' → '.join(manifest['scene_order'])}`",
        "- DMRS：固定 comb 24。",
        "- 状态：等待研究者确认确切 SHA-256；确认前不得启动正式链路。",
        "",
        "| gate | 场景 | candidate | family | 角色 | 胜出目标 | delay-grid coordinates | out10 | out1 |",
        "|---|---|---|---|---|---|---|---:|---:|",
    ]
    for row in rows:
        roles = ",".join(row["roles"])
        targets = ",".join(row["outage_winner_for"])
        delays = json.dumps(
            row["delay_grid_coordinates"], separators=(",", ":")
        )
        lines.append(
            f"| {row['gate_id']} | {row['scenario_id']} | "
            f"{row['candidate_id']} | {row['family']} | `{roles}` | "
            f"`{targets}` | `{delays}` | "
            f"{float(row['outage10_snr_db']):.3f} | "
            f"{float(row['outage1_snr_db']):.3f} |"
        )
    lines.extend(
        [
            "",
            "预算：A5一次20-trial集成preflight；prescan每点400 trials、0.5 dB步长、目标中心为ideal-outage+5 dB、最大扩展2 dB；10%/1% refine每点3,000 trials、0.25 dB步长。",
            "",
        ]
    )
    return "\n".join(lines)


def build_envelope_link_manifest(output: Path) -> tuple[Path, str, list[dict[str, object]]]:
    stage = output / "e4_link_gate"
    stage.mkdir(parents=True, exist_ok=True)
    targets = read_csv(output / "e1_outage" / "e1_outage_targets.csv")
    rows = envelope_link_candidates(targets)
    scene_counts = {
        scenario_id: sum(row["scenario_id"] == scenario_id for row in rows)
        for scenario_id in LINK_SCENARIO_ORDER
    }
    manifest: dict[str, object] = {
        "schema": "plan027-e4-envelope-link-manifest-v4",
        "status": "awaiting_researcher_confirmation",
        "run_id": output.name,
        "selection_scope": (
            "Previous 49-curve addendum plus A10/A30 AP_TU_NT and "
            "AP_TALIAS_NT system-geometry baselines"
        ),
        "curve_count": len(rows),
        "scene_counts": scene_counts,
        "scene_order": list(LINK_SCENARIO_ORDER),
        "execution_order": ["A5", "A10", "A30"],
        "dmrs_spacing_subcarriers": 24,
        "physical_definition": {
            "K_active_subcarriers": K,
            "subcarrier_spacing_hz": SCS_HZ,
            "phase_reference": "first_active_subcarrier",
            "phase_denominator": K,
            "delay_input_field": "delay_grid_coordinates",
            "continuous_coordinates_must_not_be_rounded": True,
            "precoder_normalized": False,
        },
        "trial_budget": {
            "smoke_first_scene_only": True,
            "smoke_trials_per_snr": 20,
            "prescan_trials_per_snr": 400,
            "prescan_step_db": 0.5,
            "initial_bler_gap_from_outage_db": 5.0,
            "maximum_extension_each_target_db": 2.0,
            "refine_trials_per_snr": 3000,
            "refine_step_db": 0.25,
            "minimum_errors_in_half_to_twice_target_band": 30,
        },
        "required_approvals": {
            "manifest_sha256": True,
            "first_scene_start": True,
            "next_scene_after_prior_report": True,
        },
        "amendment": {
            "authorized_request": (
                "Run A10 then A30 with useful-symbol/Nt and "
                "pilot-alias-period/Nt included in the baseline set"
            ),
            "previous_curve_count": 49,
            "added_curve_count": 4,
            "added_scenarios": ["A10", "A30"],
            "smoke_repeated": False,
        },
        "candidates": rows,
    }
    manifest_path = stage / "e4_link_manifest.json"
    historical_path = stage / "e4_link_manifest_pre_system_baseline_addendum.json"
    if manifest_path.exists() and not historical_path.exists():
        historical_path.write_bytes(manifest_path.read_bytes())
    previous_addendum_path = (
        stage / "e4_link_manifest_pre_a10_a30_system_baseline_addendum.json"
    )
    if manifest_path.exists() and not previous_addendum_path.exists():
        previous_addendum_path.write_bytes(manifest_path.read_bytes())
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    digest = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    (stage / "e4_link_manifest.sha256").write_text(
        f"{digest}  e4_link_manifest.json\n",
        encoding="utf-8",
    )
    csv_rows = []
    for row in rows:
        encoded = dict(row)
        for key in (
            "roles",
            "outage_winner_for",
            "delay_grid_coordinates",
            "delay_indices",
            "delay_ns",
        ):
            encoded[key] = json.dumps(row[key], separators=(",", ":"))
        csv_rows.append(encoded)
    write_csv_rows(csv_rows, stage / "e4_link_gate.csv")
    (stage / "e4_link_gate.md").write_text(
        link_markdown_table(rows, digest, manifest),
        encoding="utf-8",
    )
    save_json(
        {
            "manifest_sha256": digest,
            "status": "awaiting_researcher_confirmation",
            "curve_count": len(rows),
            "scene_counts": scene_counts,
            "scene_order": list(LINK_SCENARIO_ORDER),
        },
        stage / "link_gate_summary.json",
    )
    save_json(
        {
            "approved": False,
            "manifest_sha256": digest,
            "approved_curve_count": 53,
            "approved_dmrs_spacing_subcarriers": 24,
            "approved_trial_budget": manifest["trial_budget"],
            "approved_next_scenes": ["A10", "A30"],
            "researcher_confirmation": "",
        },
        stage / "e4_link_approval_template.json",
    )
    return manifest_path, digest, rows


def markdown_table(rows: list[dict[str, object]], digest: str) -> str:
    lines = [
        "# Plan-027 E4 候选门控清单",
        "",
        f"- manifest SHA-256：`{digest}`",
        f"- 去重后场景×候选数：{len(rows)}",
        "- 状态：未确认；不得启动 E4 estimated-CSI BLER。",
        "",
        "| gate | 来源 | 场景 | 原候选 | family | delay index | 10%增益(dB) | 1%增益(dB) | DMRS建议 |",
        "|---|---|---|---|---|---|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            "| {gate_id} | {source_experiment} | {scenario_id} | "
            "{original_candidate_id} | {family} | `{delay_indices}` | "
            "{gain10_db:.3f} | {gain1_db:.3f} | `{dmrs}` |".format(
                **row,
                dmrs=json.dumps(row["recommended_dmrs_spacings"]),
            )
        )
    lines.extend(
        [
            "",
            "确认 E4 时必须同时确认上述 manifest SHA-256、场景、候选、DMRS 变体和粗扫/加密 trial 预算。",
            "",
        ]
    )
    return "\n".join(lines)


def run(output: Path) -> None:
    stage = output / "e4_gate"
    stage.mkdir(parents=True, exist_ok=True)
    rows027 = positive_027(output)
    rows026 = positive_026()
    rows = deduplicate([*rows027, *rows026])
    csv_rows = []
    for row in rows:
        encoded = dict(row)
        for key in ("provenance", "duplicate_candidate_ids", "recommended_dmrs_spacings"):
            encoded[key] = json.dumps(row[key], separators=(",", ":"))
        csv_rows.append(encoded)
    write_csv_rows(csv_rows, stage / "e4_candidate_gate.csv")
    manifest = {
        "schema": "plan027-e4-candidate-manifest-v1",
        "status": "awaiting_researcher_confirmation",
        "run_id": output.name,
        "source_counts_before_deduplication": {
            "027": len(rows027),
            "026": len(rows026),
        },
        "candidate_count_after_deduplication": len(rows),
        "candidates": rows,
        "required_second_confirmation": {
            "manifest_sha256": True,
            "trial_budget": True,
            "dense_dmrs_variants": True,
        },
    }
    manifest_path = stage / "e4_candidate_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    digest = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    (stage / "e4_candidate_manifest.sha256").write_text(
        f"{digest}  e4_candidate_manifest.json\n",
        encoding="utf-8",
    )
    (stage / "e4_candidate_gate.md").write_text(
        markdown_table(rows, digest),
        encoding="utf-8",
    )
    save_json(
        {
            "manifest_sha256": digest,
            "status": "awaiting_researcher_confirmation",
            "candidate_count": len(rows),
            "source_027_before_deduplication": len(rows027),
            "source_026_before_deduplication": len(rows026),
        },
        stage / "gate_summary.json",
    )
    print(json.dumps({"sha256": digest, "candidate_count": len(rows)}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", default=DEFAULT_RUN_ID)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--scope",
        choices=("positive-all", "envelope-link"),
        default="positive-all",
    )
    args = parser.parse_args()
    output = args.output_root / str(args.run_id)
    if args.scope == "envelope-link":
        manifest_path, digest, rows = build_envelope_link_manifest(output)
        print(
            json.dumps(
                {
                    "manifest": str(manifest_path),
                    "sha256": digest,
                    "curve_count": len(rows),
                },
                indent=2,
            )
        )
    else:
        run(output)


if __name__ == "__main__":
    main()
