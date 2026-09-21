"""Analyze Plan-027 E1--E3, gate E4, and write paired result documents."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.run_plan025_delay_matched_tdl import stable_seed, write_csv_rows
from tools.run_plan027_ce_density import build_metrics
from tools.run_plan027_meff_design import (
    BOOTSTRAP_REPEATS,
    DEFAULT_OUTPUT_ROOT,
    DEFAULT_RUN_ID,
    SCENARIOS,
    paired_bootstrap_gain,
    read_csv,
    save_json,
    decode,
)


FAMILIES = ("AP_T2_CTRL", "GEO_T1_CTRL", "MEFF_T2_CAND")
NEW_SCHEMES = ("S0_SIDON", *FAMILIES)
COMPARISON_BASELINES = ("B0_QC", "AP_RMS_T1", "AP_TEPS_T1")


def target_field(probability: float) -> str:
    return "outage10_snr_db" if probability == 0.10 else "outage1_snr_db"


def gain_field(probability: float) -> str:
    return (
        "gain_vs_ap_t1_base_10_db"
        if probability == 0.10
        else "gain_vs_ap_t1_base_1_db"
    )


def best_family(
    rows: list[dict[str, str]],
    family: str,
    probability: float,
) -> dict[str, str]:
    return min(
        [row for row in rows if row["family"] == family],
        key=lambda row: float(row[target_field(probability)]),
    )


def curve_lookup(path: Path) -> dict[tuple[str, str], tuple[np.ndarray, np.ndarray]]:
    rows = read_csv(path)
    grouped: dict[tuple[str, str], list[tuple[float, float]]] = {}
    for row in rows:
        grouped.setdefault((row["scenario_id"], row["candidate_id"]), []).append(
            (float(row["snr_db"]), float(row["outage_probability"]))
        )
    output = {}
    for key, values in grouped.items():
        ordered = sorted(values)
        output[key] = (
            np.asarray([value[0] for value in ordered]),
            np.asarray([value[1] for value in ordered]),
        )
    return output


def thresholds_for(output: Path, scenario_id: str) -> tuple[dict[str, np.ndarray], np.ndarray]:
    payload = np.load(
        output / "e1_outage" / f"{scenario_id}_paired_threshold_snr_db.npz"
    )
    ids = [str(value) for value in payload["candidate_ids"]]
    thresholds = np.asarray(payload["threshold_snr_db"], dtype=np.float32)
    return {candidate_id: thresholds[index] for index, candidate_id in enumerate(ids)}, np.asarray(
        payload["snr_grid_db"], dtype=np.float64
    )


def family_comparisons(output: Path, targets: list[dict[str, str]]) -> list[dict[str, object]]:
    curves = curve_lookup(output / "e1_outage" / "e1_outage_curves.csv")
    result = []
    for scenario_id, _ in SCENARIOS:
        rows = [row for row in targets if row["scenario_id"] == scenario_id]
        thresholds, snr_grid = thresholds_for(output, scenario_id)
        for probability in (0.10, 0.01):
            base_id = rows[0][
                "ap_t1_base_10_candidate_id"
                if probability == 0.10
                else "ap_t1_base_1_candidate_id"
            ]
            by_id = {row["candidate_id"]: row for row in rows}
            comparisons = (
                ("AP_T2_CTRL_vs_AP_T1_BASE", best_family(rows, "AP_T2_CTRL", probability), by_id[base_id]),
                ("MEFF_T2_vs_GEO_T1", best_family(rows, "MEFF_T2_CAND", probability), best_family(rows, "GEO_T1_CTRL", probability)),
                ("MEFF_T2_vs_AP_T2", best_family(rows, "MEFF_T2_CAND", probability), best_family(rows, "AP_T2_CTRL", probability)),
            )
            for label, candidate, baseline in comparisons:
                _, candidate_curve = curves[(scenario_id, candidate["candidate_id"])]
                _, baseline_curve = curves[(scenario_id, baseline["candidate_id"])]
                gains = paired_bootstrap_gain(
                    thresholds[candidate["candidate_id"]],
                    thresholds[baseline["candidate_id"]],
                    candidate_curve,
                    baseline_curve,
                    snr_grid,
                    probability,
                    BOOTSTRAP_REPEATS,
                    stable_seed(
                        20260727,
                        scenario_id,
                        probability,
                        label,
                        "analysis_bootstrap",
                    ),
                )
                low, high = np.percentile(gains, [2.5, 97.5])
                result.append(
                    {
                        "scenario_id": scenario_id,
                        "target_outage_probability": probability,
                        "comparison": label,
                        "candidate_id": candidate["candidate_id"],
                        "candidate_family": candidate["family"],
                        "baseline_candidate_id": baseline["candidate_id"],
                        "baseline_family": baseline["family"],
                        "gain_snr_db": float(baseline[target_field(probability)])
                        - float(candidate[target_field(probability)]),
                        "gain_ci95_low_db": float(low),
                        "gain_ci95_high_db": float(high),
                        "selection_note": "post-outage family envelope among pre-frozen representatives",
                    }
                )
    return result


def rank(values: list[float]) -> np.ndarray:
    order = np.argsort(np.asarray(values, dtype=np.float64), kind="mergesort")
    ranks = np.empty(len(order), dtype=np.float64)
    ranks[order] = np.arange(len(order), dtype=np.float64)
    return ranks


def ranking_rows(
    targets: list[dict[str, str]],
    ce_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    ce_lookup = {}
    for row in ce_rows:
        if (
            row["dmrs_spacing_subcarriers"] == "24"
            and row["working_point"] == "common_baseline_target"
            and float(row["target_outage_probability"]) == 0.10
        ):
            ce_lookup.setdefault((row["scenario_id"], row["candidate_id"]), float(row["ce_nmse_db"]))
    output = []
    for scenario_id, _ in SCENARIOS:
        rows = [
            row
            for row in targets
            if row["scenario_id"] == scenario_id
            and row["family"] in ("B0_QC", "S0_SIDON", *FAMILIES)
        ]
        m2 = [float(row["m2_eff"]) for row in rows]
        m4 = [float(row["m4_eff"]) for row in rows]
        outage = [float(row["outage10_snr_db"]) for row in rows]
        ce = [ce_lookup[(scenario_id, row["candidate_id"])] for row in rows]
        output.append(
            {
                "scenario_id": scenario_id,
                "candidate_count": len(rows),
                "spearman_m2_vs_outage10": float(np.corrcoef(rank(m2), rank(outage))[0, 1]),
                "spearman_m4_vs_outage10": float(np.corrcoef(rank(m4), rank(outage))[0, 1]),
                "spearman_m2_vs_ce_nmse_db": float(np.corrcoef(rank(m2), rank(ce))[0, 1]),
                "spearman_m4_vs_ce_nmse_db": float(np.corrcoef(rank(m4), rank(ce))[0, 1]),
            }
        )
    return output


def density_rows(
    targets: list[dict[str, str]],
    ce_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    by_id = {(row["scenario_id"], row["candidate_id"]): row for row in targets}
    grouped: dict[tuple[str, str, int], list[float]] = {}
    lookup = {}
    for row in ce_rows:
        if row["working_point"] != "common_baseline_target":
            continue
        lookup[
            (
                row["scenario_id"],
                row["candidate_id"],
                int(row["dmrs_spacing_subcarriers"]),
                float(row["target_outage_probability"]),
            )
        ] = float(row["ce_nmse_db"])
    for row in targets:
        if row["family"] not in FAMILIES:
            continue
        if float(row["gain_vs_ap_t1_base_10_db"]) <= 0.0 and float(
            row["gain_vs_ap_t1_base_1_db"]
        ) <= 0.0:
            continue
        for spacing in (24, 12, 6):
            deltas = []
            for probability in (0.10, 0.01):
                key = (row["scenario_id"], row["candidate_id"], spacing, probability)
                if key not in lookup:
                    continue
                baseline_id = row[
                    "ap_t1_base_10_candidate_id"
                    if probability == 0.10
                    else "ap_t1_base_1_candidate_id"
                ]
                deltas.append(
                    lookup[key]
                    - lookup[(row["scenario_id"], baseline_id, spacing, probability)]
                )
            if deltas:
                grouped[(row["scenario_id"], row["candidate_id"], spacing)] = deltas
    output = []
    for (scenario_id, candidate_id, spacing), deltas in sorted(grouped.items()):
        comb24 = grouped[(scenario_id, candidate_id, 24)]
        output.append(
            {
                "scenario_id": scenario_id,
                "candidate_id": candidate_id,
                "family": by_id[(scenario_id, candidate_id)]["family"],
                "dmrs_spacing_subcarriers": spacing,
                "pilot_count_per_symbol": 576 // spacing,
                "maximum_ce_degradation_vs_baseline_db": max(deltas),
                "ce_degradation_reduction_vs_comb24_db": max(comb24) - max(deltas),
                "still_above_0p5db_trigger": max(deltas) > 0.5,
            }
        )
    return output


def best_summary(targets: list[dict[str, str]]) -> list[dict[str, object]]:
    output = []
    for scenario_id, _ in SCENARIOS:
        rows = [row for row in targets if row["scenario_id"] == scenario_id]
        by_id = {row["candidate_id"]: row for row in rows}
        for probability in (0.10, 0.01):
            baseline_id = rows[0][
                "ap_t1_base_10_candidate_id"
                if probability == 0.10
                else "ap_t1_base_1_candidate_id"
            ]
            baseline = by_id[baseline_id]
            item: dict[str, object] = {
                "scenario_id": scenario_id,
                "target_outage_probability": probability,
                "ap_t1_base_candidate_id": baseline_id,
                "ap_t1_base_snr_db": float(baseline[target_field(probability)]),
            }
            for family in ("AP_T2_CTRL", "GEO_T1_CTRL", "MEFF_T2_CAND"):
                best = best_family(rows, family, probability)
                prefix = {
                    "AP_T2_CTRL": "ap_t2",
                    "GEO_T1_CTRL": "geo_t1",
                    "MEFF_T2_CAND": "meff_t2",
                }[family]
                item[f"{prefix}_candidate_id"] = best["candidate_id"]
                item[f"{prefix}_snr_db"] = float(best[target_field(probability)])
                item[f"{prefix}_gain_vs_ap_t1_db"] = float(best[gain_field(probability)])
                item[f"{prefix}_delay_indices"] = best["delay_indices"]
                item[f"{prefix}_selection_rule"] = best["selection_rule"]
            s0 = next(row for row in rows if row["family"] == "S0_SIDON")
            item["s0_snr_db"] = float(s0[target_field(probability)])
            item["s0_gain_vs_ap_t1_db"] = float(s0[gain_field(probability)])
            output.append(item)
    return output


def comprehensive_outage_rows(
    targets: list[dict[str, str]],
) -> list[dict[str, object]]:
    output = []
    for scenario_id, _ in SCENARIOS:
        rows = [row for row in targets if row["scenario_id"] == scenario_id]
        by_family: dict[str, list[dict[str, str]]] = {}
        for row in rows:
            by_family.setdefault(row["family"], []).append(row)
        for probability in (0.10, 0.01):
            field = target_field(probability)
            baselines = {
                family: min(by_family[family], key=lambda row: float(row[field]))
                for family in COMPARISON_BASELINES
            }
            best_three_family, best_three = min(
                baselines.items(),
                key=lambda item: float(item[1][field]),
            )
            for scheme in NEW_SCHEMES:
                if scheme == "S0_SIDON":
                    candidate = by_family[scheme][0]
                    scope = "fixed_single_candidate"
                else:
                    candidate = min(
                        by_family[scheme],
                        key=lambda row: float(row[field]),
                    )
                    scope = "post_outage_envelope_of_prefrozen_representatives"
                candidate_snr = float(candidate[field])
                output.append(
                    {
                        "scenario_id": scenario_id,
                        "target_outage_probability": probability,
                        "scheme": scheme,
                        "representative_candidate_id": candidate["candidate_id"],
                        "representative_delay_indices": candidate["delay_indices"],
                        "representative_selection_rule": candidate["selection_rule"],
                        "family_frozen_representative_count": len(by_family[scheme]),
                        "representative_scope": scope,
                        "candidate_target_snr_db": candidate_snr,
                        "b0_candidate_id": baselines["B0_QC"]["candidate_id"],
                        "b0_target_snr_db": float(baselines["B0_QC"][field]),
                        "gain_vs_b0_db": float(baselines["B0_QC"][field])
                        - candidate_snr,
                        "ap_rms_candidate_id": baselines["AP_RMS_T1"]["candidate_id"],
                        "ap_rms_target_snr_db": float(
                            baselines["AP_RMS_T1"][field]
                        ),
                        "gain_vs_ap_rms_db": float(
                            baselines["AP_RMS_T1"][field]
                        )
                        - candidate_snr,
                        "ap_teps_candidate_id": baselines["AP_TEPS_T1"][
                            "candidate_id"
                        ],
                        "ap_teps_target_snr_db": float(
                            baselines["AP_TEPS_T1"][field]
                        ),
                        "gain_vs_ap_teps_db": float(
                            baselines["AP_TEPS_T1"][field]
                        )
                        - candidate_snr,
                        "best_three_baseline_family": best_three_family,
                        "best_three_baseline_candidate_id": best_three["candidate_id"],
                        "best_three_target_snr_db": float(best_three[field]),
                        "gain_vs_best_three_db": float(best_three[field])
                        - candidate_snr,
                    }
                )
    return output


def comprehensive_ce_rows(
    targets: list[dict[str, str]],
    outage_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    by_id = {(row["scenario_id"], row["candidate_id"]): row for row in targets}
    output = []
    for scenario_id, spread_ns in SCENARIOS:
        metric = build_metrics(spread_ns)[24]
        scenario_rows = [
            row for row in outage_rows if row["scenario_id"] == scenario_id
        ]
        for summary in scenario_rows:
            probability = float(summary["target_outage_probability"])
            candidate = by_id[
                (scenario_id, str(summary["representative_candidate_id"]))
            ]
            candidate_coordinates = decode(candidate["delay_grid_coordinates"])
            own_snr = float(summary["candidate_target_snr_db"])
            own_ce = metric.evaluate_grid_coordinates(
                candidate_coordinates,
                own_snr,
            )
            row: dict[str, object] = {
                "scenario_id": scenario_id,
                "target_outage_probability": probability,
                "scheme": summary["scheme"],
                "representative_candidate_id": candidate["candidate_id"],
                "representative_selection_rule": summary[
                    "representative_selection_rule"
                ],
                "candidate_own_target_snr_db": own_snr,
                "candidate_own_target_ce_nmse_db": float(own_ce["nmse_db"]),
                "best_three_baseline_family": summary[
                    "best_three_baseline_family"
                ],
            }
            baseline_specs = (
                ("b0", summary["b0_candidate_id"]),
                ("ap_rms", summary["ap_rms_candidate_id"]),
                ("ap_teps", summary["ap_teps_candidate_id"]),
                ("best_three", summary["best_three_baseline_candidate_id"]),
            )
            for prefix, baseline_id in baseline_specs:
                baseline = by_id[(scenario_id, str(baseline_id))]
                baseline_snr = float(baseline[target_field(probability)])
                candidate_at_baseline = metric.evaluate_grid_coordinates(
                    candidate_coordinates,
                    baseline_snr,
                )
                baseline_at_baseline = metric.evaluate_grid_coordinates(
                    decode(baseline["delay_grid_coordinates"]),
                    baseline_snr,
                )
                row[f"{prefix}_candidate_id"] = baseline_id
                row[f"{prefix}_target_snr_db"] = baseline_snr
                row[f"candidate_ce_at_{prefix}_target_nmse_db"] = float(
                    candidate_at_baseline["nmse_db"]
                )
                row[f"{prefix}_ce_at_own_target_nmse_db"] = float(
                    baseline_at_baseline["nmse_db"]
                )
                row[f"ce_delta_vs_{prefix}_db"] = float(
                    candidate_at_baseline["nmse_db"]
                ) - float(baseline_at_baseline["nmse_db"])
            output.append(row)
    return output


def own_workpoint_ce_rows(
    comprehensive_ce: list[dict[str, object]],
) -> list[dict[str, object]]:
    output = []
    for row in comprehensive_ce:
        candidate_ce = float(row["candidate_own_target_ce_nmse_db"])
        compact: dict[str, object] = {
            "scenario_id": row["scenario_id"],
            "target_outage_probability": row["target_outage_probability"],
            "scheme": row["scheme"],
            "representative_candidate_id": row["representative_candidate_id"],
            "representative_selection_rule": row["representative_selection_rule"],
            "candidate_own_target_snr_db": row["candidate_own_target_snr_db"],
            "candidate_own_target_ce_nmse_db": candidate_ce,
            "best_three_baseline_family": row["best_three_baseline_family"],
        }
        for prefix in ("b0", "ap_rms", "ap_teps", "best_three"):
            baseline_ce = float(row[f"{prefix}_ce_at_own_target_nmse_db"])
            compact[f"{prefix}_candidate_id"] = row[f"{prefix}_candidate_id"]
            compact[f"{prefix}_own_target_snr_db"] = row[
                f"{prefix}_target_snr_db"
            ]
            compact[f"{prefix}_own_target_ce_nmse_db"] = baseline_ce
            compact[f"own_ce_delta_vs_{prefix}_db"] = candidate_ce - baseline_ce
        output.append(compact)
    return output


def md_table(headers: list[str], rows: list[list[object]]) -> list[str]:
    output = [
        "| " + " | ".join(headers) + " |",
        "|" + "|".join("---" for _ in headers) + "|",
    ]
    output.extend("| " + " | ".join(str(value) for value in row) + " |" for row in rows)
    return output


def render_result(
    output: Path,
    support: list[dict[str, str]],
    best: list[dict[str, object]],
    comparisons: list[dict[str, object]],
    rankings: list[dict[str, object]],
    densities: list[dict[str, object]],
    comprehensive_outage: list[dict[str, object]],
    comprehensive_ce: list[dict[str, object]],
    own_workpoint_ce: list[dict[str, object]],
    gate: dict[str, object],
    text_only: bool,
) -> str:
    example_ce = next(
        row
        for row in comprehensive_ce
        if row["scenario_id"] == "A1"
        and row["scheme"] == "S0_SIDON"
        and math.isclose(float(row["target_outage_probability"]), 0.1)
    )
    title = (
        "# result-027-text：协方差加权有效矩、Sidon 适用条件与 E4 门控（无图版）"
        if text_only
        else "# result-027：协方差加权有效矩、Sidon 适用条件与 E4 门控"
    )
    lines = [
        title,
        "",
        "> 对应 `research/plan-027-有效矩-Sidon-DMRS.md`。E1–E3 和 E4 候选门控已完成；E4 estimated-CSI BLER 因缺少第二次 manifest 确认而按计划停止。",
        "",
        "## 1. 结论与状态",
        "",
        "事实：五个场景、两个 outage 目标均形成双侧 bracket。`AP_TEPS_T1` 在 10% 和 1% 两个目标上始终是 `AP_T1_BASE` 的胜出组成项。完整协方差搜索与解除等差约束的联合设计相对该直接物理公差基线有明确 ideal-CSI outage 价值，但优势随物理展宽增加而减小。",
        "",
        "事实：匹配非等差候选空间后，`MEFF_T2_CAND` 与 `GEO_T1_CTRL` 的最佳冻结代表基本持平；不能据本轮声称完整协方差知识具有独立且有实际幅度的 outage 增益。T2 的主要可见价值来自相对直接物理公差等差基线的候选空间放宽和设计优化总体作用。",
        "",
        "事实：E3 显示多数正 outage 候选相对主基线有显著 CE 劣化。提高导频密度能降低部分劣化，但改变了导频 RE 数、总导频能量和开销，且多数场景 comb 6 后仍超过 0.5 dB 门槛。因此 ideal-CSI 潜力不能替代 estimated-CSI BLER。",
        "",
        f"门控事实：027 的 56 个正点估计条目与 026 的 12 个条目按场景和规范 delay key 去重后得到 {gate['candidate_count']} 个场景×候选；manifest SHA-256 为 `{gate['manifest_sha256']}`。当前状态是 `awaiting_researcher_confirmation`，未运行 E4。",
        "",
        "## 2. 仿真条件与复现",
        "",
        "- 系统：48 PRB、$K=576$、30 kHz SCS、4096 FFT、288 CP、10 OFDM symbols、8 Tx / 1 Rx / 1 layer。",
        "- 默认 DMRS：symbols `[2,7]`、comb 24、48 pilot RE、5712 data RE；E3 comb 12/6 分别为 96/192 pilot RE，共同数据交集为 5568 RE。",
        "- 信道：Sionna 1.0.2 TDL-A，RMS delay spread `[1,5,10,30,100]` ns，0 km/h，3.5 GHz。",
        "- CDD：$V_{k,n}=\\exp(-j2\\pi k\\Delta f\\tau_n)$；网格 index 分辨率 57.870370 ns。`AP_RMS_T1`/`AP_TEPS_T1` 直接使用连续物理时延，不做 index 舍入。",
        "- 搜索：每场景 200,000 个唯一 T2 状态、200,000 个 T1 几何状态；有效矩容差 `[0.10,0.25,0.50] dB`。",
        "- outage：每场景 200,000 个共同信道样本，SNR 0–24 dB、0.5 dB 步长及预授权 4 dB 扩展，1,000 次成对 bootstrap。",
        "- 目标谱效率：16QAM Gray BICM，$4\\times553/1024=2.1602$ bit/RE。",
        f"- 展开配置：`outputs/experiment027_meff_sidon/{output.name}/e1_outage/resolved_experiment.json`；环境：`outputs/experiment027_meff_sidon/{output.name}/validation/environment.json`。",
        "",
        "核心命令：",
        "",
        "```powershell",
        f"& D:\\venvs\\cdd-s102\\Scripts\\python.exe tools\\run_plan027_meff_design.py --stage validate --run-id {output.name}",
        f"& D:\\venvs\\cdd-s102\\Scripts\\python.exe tools\\run_plan027_meff_design.py --stage search --run-id {output.name}",
        f"& D:\\venvs\\cdd-s102\\Scripts\\python.exe tools\\run_plan027_meff_design.py --stage outage --run-id {output.name}",
        f"& D:\\venvs\\cdd-s102\\Scripts\\python.exe tools\\run_plan027_ce_density.py --stage ce --run-id {output.name}",
        f"& D:\\venvs\\cdd-s102\\Scripts\\python.exe tools\\build_plan027_bler_gate.py --run-id {output.name}",
        f"& D:\\venvs\\cdd-s102\\Scripts\\python.exe tools\\enrich_plan027_candidate_outputs.py --run-id {output.name}",
        f"& D:\\venvs\\cdd-s102\\Scripts\\python.exe tools\\analyze_plan027.py --run-id {output.name}",
        "```",
        "",
        "## 3. Phase 0 与 Sidon 几何",
        "",
        "Phase 0 的三组直接/lag 有效矩误差均小于 $2\\times10^{-11}$。A100 协方差 Hermitian，最小特征值 `-3.62e-13`（数值误差量级），对角为 1。搜索同 seed 可重放；comb 24/12/6 的 pilot 数、fold 周期与共同数据集合通过验证。",
        "",
    ]
    lines.extend(
        md_table(
            ["场景", "$T_\\epsilon$ (ns)", "要求 pair/fold", "硬样本数/200k", "硬装填判定", "S0硬厚"],
            [
                [
                    row["scenario_id"],
                    f"{float(row['support_width_ns']):.3f}",
                    f"{row['required_pair_gap_indices']}/{row['required_fold_gap_indices']}",
                    row.get("hard_candidates_found", ""),
                    row["hard_feasible_by_packing_bound"],
                    row.get("s0_hard_thick_sidon", ""),
                ]
                for row in support
            ],
        )
    )
    lines.extend(["", "几何层级事实：A1/A5/A10/A30 的装填必要条件可行；A100 的要求 `17/9` 超过上界 `16/3`，硬不可行。严格 Sidon 的加性构造始终成立，但 S0 只在 A1/A5 同时满足本轮厚保护距离。几何成立本身不等价于 outage 或 BLER 成立。", "", "两条 T1 连续物理公差基线及包络胜出者：", ""])
    target_source = {
        row["scenario_id"]: row
        for row in read_csv(output / "e1_outage" / "e1_outage_targets.csv")
        if row["family"] == "AP_RMS_T1"
    }
    teps_source = {
        row["scenario_id"]: row
        for row in read_csv(output / "e1_outage" / "e1_outage_targets.csv")
        if row["family"] == "AP_TEPS_T1"
    }
    lines.extend(
        md_table(
            ["场景", "AP_RMS_T1 delay(ns)", "AP_TEPS_T1 delay(ns)", "10%胜出", "1%胜出"],
            [
                [
                    scenario_id,
                    f"`{target_source[scenario_id]['delay_ns']}`",
                    f"`{teps_source[scenario_id]['delay_ns']}`",
                    target_source[scenario_id]["ap_t1_base_10_candidate_id"],
                    target_source[scenario_id]["ap_t1_base_1_candidate_id"],
                ]
                for scenario_id, _ in SCENARIOS
            ],
        )
    )
    lines.extend(["", "## 4. ideal-CSI outage", ""])
    lines.extend(
        md_table(
            ["场景", "目标", "AP_T1基线(dB)", "AP_T2增益", "GEO增益", "MEFF增益", "S0增益"],
            [
                [
                    row["scenario_id"],
                    f"{100*float(row['target_outage_probability']):.0f}%",
                    f"{float(row['ap_t1_base_snr_db']):.3f}",
                    f"{float(row['ap_t2_gain_vs_ap_t1_db']):+.3f}",
                    f"{float(row['geo_t1_gain_vs_ap_t1_db']):+.3f}",
                    f"{float(row['meff_t2_gain_vs_ap_t1_db']):+.3f}",
                    f"{float(row['s0_gain_vs_ap_t1_db']):+.3f}",
                ]
                for row in best
            ],
        )
    )
    lines.extend(["", "表中各 family 取预冻结代表的后验包络，仅用于汇总，不能忽略 family 内多重比较。逐候选点估计和相对 AP_T1 的预定成对区间见 `e1_outage/e1_outage_targets.csv`；块级配对阈值见同目录五个 `*_paired_threshold_snr_db.npz`。", "", "匹配比较：", ""])
    lines.extend(
        md_table(
            ["场景", "目标", "比较", "增益(dB)", "95%区间(dB)"],
            [
                [
                    row["scenario_id"],
                    f"{100*float(row['target_outage_probability']):.0f}%",
                    row["comparison"],
                    f"{float(row['gain_snr_db']):+.4f}",
                    f"[{float(row['gain_ci95_low_db']):+.4f}, {float(row['gain_ci95_high_db']):+.4f}]",
                ]
                for row in comparisons
            ],
        )
    )
    lines.extend(["", "事实：`MEFF_T2_vs_GEO_T1` 的差值均很小；即使部分区间不跨 0，其幅度也远低于相对 AP_T1 的总体增益，且比较使用 family 后验包络。因此本轮不把它解释为协方差知识的独立实用价值。", "", "各场景最佳构造与 delay：", ""])
    lines.extend(
        md_table(
            ["场景", "目标", "AP_T2", "GEO_T1", "MEFF_T2"],
            [
                [
                    row["scenario_id"],
                    f"{100*float(row['target_outage_probability']):.0f}%",
                    f"{row['ap_t2_candidate_id']} `{row['ap_t2_delay_indices']}`",
                    f"{row['geo_t1_candidate_id']} `{row['geo_t1_delay_indices']}`",
                    f"{row['meff_t2_candidate_id']} `{row['meff_t2_delay_indices']}`",
                ]
                for row in best
            ],
        )
    )
    lines.extend(["", "上述最佳代表的冻结规则：", ""])
    lines.extend(
        md_table(
            ["场景", "目标", "AP_T2规则", "GEO规则", "MEFF规则"],
            [
                [
                    row["scenario_id"],
                    f"{100*float(row['target_outage_probability']):.0f}%",
                    row["ap_t2_selection_rule"],
                    row["geo_t1_selection_rule"],
                    row["meff_t2_selection_rule"],
                ]
                for row in best
            ],
        )
    )
    lines.extend(["", "候选的最低 $M_2^{\\rm eff}$、受限最低 $M_4^{\\rm eff}$、最低 $M_4^{\\rm eff}$ 消融和 Pareto provenance 均保存在 `e1_search/e1_frozen_candidates.csv` 的 `selection_rule` 字段。完整 1,000,000 个 T2 状态与 1,000,000 个几何状态分别在 `e1_all_candidates.csv` 和 `e2_sidon_candidates.csv`。"])
    lines.extend(
        [
            "",
            "冻结与后验包络的精确定义：",
            "",
            "- `AP_T2_CTRL`：每场景先只按完整协方差下较低的 $M_2^{\\rm eff}/M_4^{\\rm eff}$ 冻结 2 个代表，再把这 2 个都运行 outage。",
            "- `GEO_T1_CTRL`：只按 pair/fold gap 与硬/软厚 Sidon 几何冻结；A1/A5/A10 各 3 个，A30/A100 各 1 个，再全部运行 outage。",
            "- `MEFF_T2_CAND`：每场景先按最低 $M_2^{\\rm eff}$、受限最低 $M_4^{\\rm eff}$、最低 $M_4^{\\rm eff}$ 消融与 $M_2/M_4$ Pareto 规则冻结 7 个，再全部运行 outage。",
            "- `S0_SIDON`：固定单一候选，不存在 family 内后验挑选。",
            "- `pareto_representative` 只表示在“两个指标都越低越好”的 $(M_2^{\\rm eff},M_4^{\\rm eff})$ 平面上不被另一搜索状态同时支配，并从该前沿取代表；该标签在 outage 运行前确定，未读取 outage。",
            "",
            "因此，下表的 `AP_T2/GEO/MEFF` 数字是各自预冻结代表全部完成 outage 后，在同一 family 内取目标 SNR 最低者形成的后验包络；它适合浏览 family 潜力，但不是单个代理指标预先指定候选的无选择偏差估计。",
            "",
            "四类新方案相对三条基线及三者最优者的全面对比（正增益表示新方案所需 SNR 更低）：",
            "",
        ]
    )
    lines.extend(
        md_table(
            [
                "场景",
                "目标",
                "新方案",
                "代表/冻结数",
                "方案SNR",
                "vs B0",
                "vs AP_RMS",
                "vs AP_TEPS",
                "vs 三者最优",
                "三者最优",
            ],
            [
                [
                    row["scenario_id"],
                    f"{100*float(row['target_outage_probability']):.0f}%",
                    row["scheme"],
                    f"{row['representative_candidate_id']}/{row['family_frozen_representative_count']}",
                    f"{float(row['candidate_target_snr_db']):.3f}",
                    f"{float(row['gain_vs_b0_db']):+.3f}",
                    f"{float(row['gain_vs_ap_rms_db']):+.3f}",
                    f"{float(row['gain_vs_ap_teps_db']):+.3f}",
                    f"{float(row['gain_vs_best_three_db']):+.3f}",
                    row["best_three_baseline_family"],
                ]
                for row in comprehensive_outage
            ],
        )
    )
    lines.extend(
        [
            "",
            "**各自工作点的 comb-24 CE 全面对比。** 下表中，新方案和每条基线都在各自达到同一 outage 目标 $p$ 的 SNR 处计算 CE，即",
            "",
            "$$\\Delta_b^{\\mathrm{own}}(s;p)=L_{\\mathrm{CE,dB}}(s;\\gamma_p(s))-L_{\\mathrm{CE,dB}}(b;\\gamma_p(b)).$$",
            "",
            "这里 $\\Delta^{\\mathrm{own}}>0$ 表示新方案在各自工作点上的 CE NMSE 更大。由于两端 SNR 通常不同，该差值同时包含“工作点 SNR 不同”和“方案 CE 结构不同”的影响；它回答端到端各自工作点的 CE 水平，不等价于第 5 节后面的同 SNR CE 结构对比。`三者最优` 仍按三条基线中 outage 目标 SNR 最低者选择。",
            "",
        ]
    )
    lines.extend(
        md_table(
            [
                "场景",
                "目标",
                "新方案",
                "代表",
                "新SNR",
                "新CE",
                "B0 SNR",
                "B0 CE",
                "Δown-B0",
                "RMS SNR",
                "RMS CE",
                "Δown-RMS",
                "TEPS SNR",
                "TEPS CE",
                "Δown-TEPS",
                "三者最优",
                "最优SNR",
                "最优CE",
                "Δown-最优",
            ],
            [
                [
                    row["scenario_id"],
                    f"{100*float(row['target_outage_probability']):.0f}%",
                    row["scheme"],
                    row["representative_candidate_id"],
                    f"{float(row['candidate_own_target_snr_db']):.3f}",
                    f"{float(row['candidate_own_target_ce_nmse_db']):.3f}",
                    f"{float(row['b0_own_target_snr_db']):.3f}",
                    f"{float(row['b0_own_target_ce_nmse_db']):.3f}",
                    f"{float(row['own_ce_delta_vs_b0_db']):+.3f}",
                    f"{float(row['ap_rms_own_target_snr_db']):.3f}",
                    f"{float(row['ap_rms_own_target_ce_nmse_db']):.3f}",
                    f"{float(row['own_ce_delta_vs_ap_rms_db']):+.3f}",
                    f"{float(row['ap_teps_own_target_snr_db']):.3f}",
                    f"{float(row['ap_teps_own_target_ce_nmse_db']):.3f}",
                    f"{float(row['own_ce_delta_vs_ap_teps_db']):+.3f}",
                    row["best_three_baseline_family"],
                    f"{float(row['best_three_own_target_snr_db']):.3f}",
                    f"{float(row['best_three_own_target_ce_nmse_db']):.3f}",
                    f"{float(row['own_ce_delta_vs_best_three_db']):+.3f}",
                ]
                for row in own_workpoint_ce
            ],
        )
    )
    lines.extend(
        [
            "",
            "注意：plan 预定主基线 `AP_T1_BASE` 只在 `AP_RMS_T1` 与 `AP_TEPS_T1` 中取优，不包含历史参考 `B0_QC`。上表新增的“三者最优”是本次对齐理解所增加的补充口径，不改写原 plan 主判据。",
            "",
            "## 5. 指标排序与 CE",
            "",
        ]
    )
    lines.extend(
        md_table(
            ["场景", "候选数", "rho(M2,out10)", "rho(M4,out10)", "rho(M2,CE)", "rho(M4,CE)"],
            [
                [
                    row["scenario_id"],
                    row["candidate_count"],
                    f"{float(row['spearman_m2_vs_outage10']):+.3f}",
                    f"{float(row['spearman_m4_vs_outage10']):+.3f}",
                    f"{float(row['spearman_m2_vs_ce_nmse_db']):+.3f}",
                    f"{float(row['spearman_m4_vs_ce_nmse_db']):+.3f}",
                ]
                for row in rankings
            ],
        )
    )
    lines.extend(
        [
            "",
            "这里的 $\\rho(x,y)$ 是 Spearman 等级相关系数，范围为 $[-1,1]$。它比较候选按两个量排序是否一致，不是物理频率相关系数。因为 $M_2/M_4$、outage 目标 SNR 和 CE NMSE 都是越低越好，正相关表示排序方向较一致，负相关表示一个量降低时另一个量倾向升高。",
            "",
            "上表“候选数”包含 `B0_QC`、`S0_SIDON`、全部预冻结 `AP_T2_CTRL`、`GEO_T1_CTRL` 和 `MEFF_T2_CAND`；不包含连续基线 `AP_RMS_T1/AP_TEPS_T1`。所以 A1/A5/A10 为 $2+2+3+7=14$，A30/A100 为 $2+2+1+7=12$。",
            "",
        ]
    )
    density_summary = {}
    for row in densities:
        key = (row["scenario_id"], int(row["dmrs_spacing_subcarriers"]))
        density_summary.setdefault(key, []).append(
            float(row["maximum_ce_degradation_vs_baseline_db"])
        )
    lines.extend(["", "CE 密度汇总（正 outage 候选中相对基线的最大劣化范围）：", ""])
    lines.extend(
        md_table(
            ["场景", "comb", "候选数", "劣化最小(dB)", "劣化最大(dB)", "仍>0.5数"],
            [
                [
                    scenario_id,
                    spacing,
                    len(values),
                    f"{min(values):+.3f}",
                    f"{max(values):+.3f}",
                    sum(value > 0.5 for value in values),
                ]
                for (scenario_id, spacing), values in sorted(density_summary.items())
            ],
        )
    )
    lines.extend(
        [
            "",
            "“仍 >0.5 数”是该场景、该 comb 下，正 outage 的 `AP_T2/GEO/MEFF` 候选中，候选 CE NMSE 相对对应 `AP_T1_BASE` 在共同基线目标 SNR 处仍劣化超过 `0.5 dB` 的候选个数；不是 outage、相关系数或 BLER 大于 0.5。S0 在原 E3 中是固定诊断参考，不参与密度触发计数。",
            "",
            "第 4 节四类新方案最佳代表的 comb-24 CE：这里的“最优”首先按 **outage 目标所需 SNR 最低**选择，不按 CE NMSE 最小选择。令 $\\gamma_p(x)$ 为方案 $x$ 达到 outage 目标 $p$ 所需的 SNR，三条基线集合为 $\\mathcal B=\\{\\mathrm{B0},\\mathrm{AP\\_RMS},\\mathrm{AP\\_TEPS}\\}$，则",
            "",
            "$$b_p^\\star=\\arg\\min_{b\\in\\mathcal B}\\gamma_p(b).$$",
            "",
            "令 $L_{\\mathrm{CE,dB}}(x;\\gamma)$ 为方案 $x$ 在 SNR $\\gamma$ 处的 comb-24 matched CE NMSE（dB）。`own CE` 为 $L_{\\mathrm{CE,dB}}(s;\\gamma_p(s))$；`候选CE@三者最优` 为 $L_{\\mathrm{CE,dB}}(s;\\gamma_p(b_p^\\star))$；`三者最优CE` 为 $L_{\\mathrm{CE,dB}}(b_p^\\star;\\gamma_p(b_p^\\star))$。注意 `own CE` 与后两列通常不在同一 SNR，不能直接相减。",
            "",
            "对任一基线 $b$，CE 差值统一在该基线自己的 outage 目标 SNR 处计算：",
            "",
            "$$\\Delta_b(s;p)=L_{\\mathrm{CE,dB}}(s;\\gamma_p(b))-L_{\\mathrm{CE,dB}}(b;\\gamma_p(b))=10\\log_{10}\\frac{\\mathrm{NMSE}_s(\\gamma_p(b))}{\\mathrm{NMSE}_b(\\gamma_p(b))}.$$",
            "",
            "因此 `Δ三者最优` 取 $b=b_p^\\star$，`ΔB0`、`ΔAP_RMS`、`ΔAP_TEPS` 分别取对应基线。$\\Delta>0$ 表示候选 NMSE 更大、CE 更差；$\\Delta<0$ 表示候选 CE 更好。",
            "",
        ]
    )
    lines.extend(
        md_table(
            [
                "场景",
                "目标",
                "方案",
                "代表",
                "own CE",
                "候选CE@三者最优",
                "三者最优CE",
                "Δ三者最优",
                "ΔB0",
                "ΔAP_RMS",
                "ΔAP_TEPS",
            ],
            [
                [
                    row["scenario_id"],
                    f"{100*float(row['target_outage_probability']):.0f}%",
                    row["scheme"],
                    row["representative_candidate_id"],
                    f"{float(row['candidate_own_target_ce_nmse_db']):.3f}",
                    f"{float(row['candidate_ce_at_best_three_target_nmse_db']):.3f}",
                    f"{float(row['best_three_ce_at_own_target_nmse_db']):.3f}",
                    f"{float(row['ce_delta_vs_best_three_db']):+.3f}",
                    f"{float(row['ce_delta_vs_b0_db']):+.3f}",
                    f"{float(row['ce_delta_vs_ap_rms_db']):+.3f}",
                    f"{float(row['ce_delta_vs_ap_teps_db']):+.3f}",
                ]
                for row in comprehensive_ce
            ],
        )
    )
    lines.extend(
        [
            "",
            "**A1、10%、S0_SIDON 示例。** 三条基线达到 10% outage 所需的 SNR 分别为 "
            f"B0 `{float(example_ce['b0_target_snr_db']):.3f} dB`、"
            f"AP_TEPS `{float(example_ce['ap_teps_target_snr_db']):.3f} dB`、"
            f"AP_RMS `{float(example_ce['ap_rms_target_snr_db']):.3f} dB`，"
            "所以“三者最优”是 B0，而不是 CE NMSE 最小的 AP_RMS。",
            "",
            f"- `own CE = {float(example_ce['candidate_own_target_ce_nmse_db']):.3f} dB`：S0 在自己的 10% outage SNR `{float(example_ce['candidate_own_target_snr_db']):.3f} dB` 处计算。",
            f"- `Δ三者最优 = ΔB0 = {float(example_ce['candidate_ce_at_b0_target_nmse_db']):.3f} - ({float(example_ce['b0_ce_at_own_target_nmse_db']):.3f}) = {float(example_ce['ce_delta_vs_b0_db']):+.3f} dB`。",
            f"- `ΔAP_RMS = {float(example_ce['candidate_ce_at_ap_rms_target_nmse_db']):.3f} - ({float(example_ce['ap_rms_ce_at_own_target_nmse_db']):.3f}) = {float(example_ce['ce_delta_vs_ap_rms_db']):+.3f} dB`。",
            f"- `ΔAP_TEPS = {float(example_ce['candidate_ce_at_ap_teps_target_nmse_db']):.3f} - ({float(example_ce['ap_teps_ce_at_own_target_nmse_db']):.3f}) = {float(example_ce['ce_delta_vs_ap_teps_db']):+.3f} dB`。",
            "",
            "A1 的 `AP_RMS_T1` 在数值容差 $10^{-10}$ 下 pilot rank 为 7 且条件数极大；它不是任何目标的 AP_T1 包络胜出者。其他冻结候选的 comb-24 pilot rank 为 8。E3 固定每个 pilot RE 功率，所以 comb 12/6 分别把 pilot RE 增至 96/192；结果只说明增加观测后的 CE 可恢复性。",
            "",
            "## 6. E4 门控与未完成项",
            "",
            f"- gate CSV：`outputs/experiment027_meff_sidon/{output.name}/e4_gate/e4_candidate_gate.csv`。",
            f"- manifest：`outputs/experiment027_meff_sidon/{output.name}/e4_gate/e4_candidate_manifest.json`。",
            f"- SHA-256：`{gate['manifest_sha256']}`。",
            f"- 去重前：027 `{gate['source_027_before_deduplication']}` 条、026 `{gate['source_026_before_deduplication']}` 条；去重后 `{gate['candidate_count']}` 条。",
            "- E4 未执行：没有研究者对确切 manifest、hash、DMRS 变体和粗扫/加密 trial 预算的第二次确认。`tools/run_plan027_bler.py` 会拒绝缺少确认回执或 hash 不匹配的请求。",
            "",
            "因此，Q4“estimated-CSI BLER 是否保持增益”仍未回答。当前结果只能评价 ideal-CSI outage 与 matched CE，不能作为最终链路胜出判据。",
            "",
            "## 7. 异常、解释边界与可进入知识的候选结论",
            "",
            "已确认事实：直接/lag 有效矩一致；五个场景 outage 目标闭合；相对直接物理公差 AP_T1 基线存在正 ideal-CSI outage 增益；A100 硬厚 Sidon 几何不可行；多数候选存在 CE 代价；E4 manifest 已冻结但未确认。",
            "",
            "推断：相对 AP_T1 的大增益主要反映候选空间放宽和协方差感知优化的总体作用，而不是完整协方差知识的独立作用。该推断由 MEFF/GEO 近似持平支持，但 family 包络是后验汇总，仍需 E4 或预先指定单一代表进一步验证。",
            "",
            "不能进入已验证知识的内容：任何 estimated-CSI BLER 增益、高密度导频的净链路收益、对 TDL-C/D/E、CDL、空间相关、移动性或协方差失配的推广。",
            "",
            "建议在研究者确认 result 后，可把“有效矩实现验证”“AP_TEPS 主基线胜出”“协方差知识未显示独立实用 outage 幅度”“A100 硬厚 Sidon 不可行”和“CE 是主要门控风险”作为带限定条件的结论更新到 `KNOWLEDGE.md`；在 E4 完成前不得更新最终 BLER 阶段验收。",
            "",
            "## 8. 证据路径",
            "",
            f"- `outputs/experiment027_meff_sidon/{output.name}/validation/`：Phase 0、复用回执、环境。",
            f"- `outputs/experiment027_meff_sidon/{output.name}/e1_search/`：全候选、冻结、Pareto 与搜索摘要。",
            f"- `outputs/experiment027_meff_sidon/{output.name}/e1_outage/`：曲线、目标、成对 bootstrap 与块级阈值。",
            f"- `outputs/experiment027_meff_sidon/{output.name}/e2_sidon_map/`：支撑、几何候选与适用性。",
            f"- `outputs/experiment027_meff_sidon/{output.name}/e3_ce_density/`：目标 SNR CE、密度触发与矩阵诊断。",
            f"- `outputs/experiment027_meff_sidon/{output.name}/e4_gate/`：去重清单、manifest 与 SHA-256。",
            f"- `outputs/experiment027_meff_sidon/{output.name}/final/`：本文表格的机器可读汇总。",
            "",
        ]
    )
    return "\n".join(lines)


def run(output: Path) -> None:
    final = output / "final"
    final.mkdir(parents=True, exist_ok=True)
    targets = read_csv(output / "e1_outage" / "e1_outage_targets.csv")
    support = read_csv(output / "e2_sidon_map" / "e2_support_geometry.csv")
    applicability = {
        row["scenario_id"]: row
        for row in read_csv(output / "e2_sidon_map" / "e2_applicability.csv")
    }
    for row in support:
        row.update(applicability[row["scenario_id"]])
    ce_rows = read_csv(output / "e3_ce_density" / "e3_ce_at_targets.csv")
    best = best_summary(targets)
    comparisons = family_comparisons(output, targets)
    rankings = ranking_rows(targets, ce_rows)
    densities = density_rows(targets, ce_rows)
    comprehensive_outage = comprehensive_outage_rows(targets)
    comprehensive_ce = comprehensive_ce_rows(targets, comprehensive_outage)
    own_workpoint_ce = own_workpoint_ce_rows(comprehensive_ce)
    gate = json.loads(
        (output / "e4_gate" / "gate_summary.json").read_text(encoding="utf-8")
    )
    write_csv_rows(best, final / "best_family_outage.csv")
    write_csv_rows(comparisons, final / "family_comparisons.csv")
    write_csv_rows(rankings, final / "metric_rank_correlations.csv")
    write_csv_rows(densities, final / "ce_density_summary.csv")
    write_csv_rows(
        comprehensive_outage,
        final / "new_scheme_baseline_outage_comparison.csv",
    )
    write_csv_rows(
        comprehensive_ce,
        final / "best_representative_ce_comparison.csv",
    )
    write_csv_rows(
        own_workpoint_ce,
        final / "own_workpoint_ce_comparison.csv",
    )
    summary = {
        "run_id": output.name,
        "status": "E1-E3_complete_E4_awaiting_second_confirmation",
        "support_geometry": support,
        "best_family_outage": best,
        "family_comparisons": comparisons,
        "metric_rank_correlations": rankings,
        "ce_density_summary": densities,
        "new_scheme_baseline_outage_comparison": comprehensive_outage,
        "best_representative_ce_comparison": comprehensive_ce,
        "own_workpoint_ce_comparison": own_workpoint_ce,
        "e4_gate": gate,
        "estimated_csi_bler_executed": False,
    }
    save_json(summary, final / "final_summary.json")
    commands = [
        [
            sys.executable,
            "tools/run_plan027_meff_design.py",
            "--stage",
            stage,
            "--run-id",
            output.name,
        ]
        for stage in ("validate", "search", "outage")
    ]
    commands.extend(
        [
            [
                sys.executable,
                "tools/run_plan027_ce_density.py",
                "--stage",
                "ce",
                "--run-id",
                output.name,
            ],
            [
                sys.executable,
                "tools/build_plan027_bler_gate.py",
                "--run-id",
                output.name,
            ],
            [
                sys.executable,
                "tools/enrich_plan027_candidate_outputs.py",
                "--run-id",
                output.name,
            ],
            [
                sys.executable,
                "tools/analyze_plan027.py",
                "--run-id",
                output.name,
            ],
        ]
    )
    save_json({"commands": commands}, output / "commands.json")
    save_json(
        {
            "run_id": output.name,
            "plan": "research/plan-027-有效矩-Sidon-DMRS.md",
            "K_active_subcarriers": 576,
            "subcarrier_spacing_hz": 30000.0,
            "n_tx": 8,
            "n_rx": 1,
            "tdl_profile": "A",
            "rms_delay_spreads_ns": [1.0, 5.0, 10.0, 30.0, 100.0],
            "search_states_per_scenario": 200000,
            "geometry_states_per_scenario": 200000,
            "outage_samples_per_scenario": 200000,
            "paired_bootstrap_repeats": 1000,
            "e4_status": "awaiting_researcher_confirmation",
        },
        output / "resolved_experiment.json",
    )
    environment = json.loads(
        (output / "validation" / "environment.json").read_text(encoding="utf-8")
    )
    save_json(environment, output / "environment.json")
    (output / "execution.log").write_text(
        "\n".join(
            [
                "Plan-027 execution summary",
                "Phase 0: passed",
                "Search/freeze: completed before outage",
                "E1 outage and paired bootstrap: completed",
                "E2 Sidon applicability: completed",
                "E3 CE density diagnostics: completed",
                "E4 manifest: generated; awaiting second researcher confirmation",
                "E4 estimated-CSI BLER: not executed",
                "",
            ]
        ),
        encoding="utf-8",
    )
    (ROOT / "research" / "result-027-有效矩-Sidon-DMRS.md").write_text(
        render_result(
            output,
            support,
            best,
            comparisons,
            rankings,
            densities,
            comprehensive_outage,
            comprehensive_ce,
            own_workpoint_ce,
            gate,
            False,
        ),
        encoding="utf-8",
    )
    (ROOT / "research" / "result-027-有效矩-Sidon-DMRS-text.md").write_text(
        render_result(
            output,
            support,
            best,
            comparisons,
            rankings,
            densities,
            comprehensive_outage,
            comprehensive_ce,
            own_workpoint_ce,
            gate,
            True,
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", default=DEFAULT_RUN_ID)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    args = parser.parse_args()
    run(args.output_root / str(args.run_id))


if __name__ == "__main__":
    main()
