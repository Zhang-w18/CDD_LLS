from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable

import yaml


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs" / "reports" / "result027_onward_target_bler_snr_gain.csv"


FIELDS = [
    "experiment",
    "comparison_group",
    "result_source",
    "result_status",
    "waveform",
    "channel_model",
    "delay_spread_ns",
    "carrier_frequency_ghz",
    "ue_speed_kmh",
    "allocated_rb",
    "occupied_rb",
    "symbols",
    "n_tx",
    "n_rx",
    "layers",
    "scs_khz",
    "fft_size",
    "cp_samples",
    "payload_bits",
    "modulation",
    "mcs_or_code",
    "dmrs",
    "csi_mode",
    "channel_estimator",
    "demodulator",
    "decoder",
    "tx_power_normalization",
    "snr_definition",
    "candidate_id",
    "candidate_family",
    "receiver_knowledge",
    "delay_coordinates",
    "delay_ns",
    "precoder_description",
    "target_bler",
    "target_snr_db",
    "target_ci95_low_db",
    "target_ci95_high_db",
    "target_status",
    "target_method",
    "baseline_candidate_id",
    "baseline_rule",
    "gain_vs_baseline_db",
    "gain_definition",
    "data_source",
    "notes",
]


def read_csv(path: str | Path) -> list[dict[str, str]]:
    with (ROOT / path).open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path: str | Path) -> Any:
    with (ROOT / path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def read_yaml(path: str | Path) -> Any:
    with (ROOT / path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def number(value: Any, digits: int = 6) -> str:
    if value in (None, ""):
        return ""
    value = float(value)
    if not math.isfinite(value):
        return ""
    return f"{value:.{digits}f}".rstrip("0").rstrip(".")


def seq(values: Iterable[Any] | None, digits: int = 6) -> str:
    if values is None:
        return ""
    return "[" + ",".join(number(v, digits) for v in values) + "]"


def local_log_crossing(points: list[dict[str, Any]], target: float) -> tuple[str, str]:
    clean = sorted(
        (float(row["snr_db"]), float(row["bler"]))
        for row in points
        if row.get("snr_db") not in (None, "") and row.get("bler") not in (None, "")
    )
    for (x0, y0), (x1, y1) in zip(clean, clean[1:]):
        if y0 >= target >= y1 and y0 > 0 and y1 > 0 and x1 > x0:
            if y0 == y1:
                return number((x0 + x1) / 2), f"raw bracket [{number(x0)},{number(x1)}]"
            fraction = (math.log10(target) - math.log10(y0)) / (math.log10(y1) - math.log10(y0))
            return number(x0 + fraction * (x1 - x0)), f"raw bracket [{number(x0)},{number(x1)}]"
    return "", "unbracketed"


def base_row(**kwargs: Any) -> dict[str, str]:
    row = {field: "" for field in FIELDS}
    for key, value in kwargs.items():
        if key not in row:
            raise KeyError(key)
        row[key] = str(value)
    return row


def pds_common(**kwargs: Any) -> dict[str, str]:
    return base_row(
        waveform="PDSCH",
        channel_model="Sionna 1.0.2 TDL-A",
        allocated_rb=48,
        occupied_rb=48,
        symbols=10,
        layers=1,
        scs_khz=30,
        fft_size=4096,
        cp_samples=288,
        payload_bits=12024,
        modulation="16QAM",
        mcs_or_code="NR 256QAM table MCS 8; target rate 553/1024",
        dmrs="symbols [2,7], comb-6, 192 pilot RE",
        demodulator="soft 16QAM demapper; coherent single-layer equalization/MRC",
        decoder="Sionna NR LDPC, max 8 iterations",
        **kwargs,
    )


def pdcch_common(**kwargs: Any) -> dict[str, str]:
    return base_row(
        waveform="PDCCH",
        allocated_rb=48,
        layers=1,
        scs_khz=30,
        fft_size=4096,
        cp_samples=288,
        modulation="QPSK",
        demodulator="coherent QPSK LLR; MRC when n_rx>1",
        decoder="CRC24C/RNTI mask; Sionna Polar hybrid CA-SCL list 8",
        tx_power_normalization="unit total transmit power per active RE",
        snr_definition="Es/N0: total transmit power per active RE divided by noise variance of each Rx branch; no n_rx scaling",
        **kwargs,
    )


def manifest_index(path: str | Path) -> dict[str, dict[str, Any]]:
    data = read_json(path)
    candidates = data["candidates"] if isinstance(data, dict) else data
    return {str(item["candidate_id"]): item for item in candidates}


def add_027_028(rows: list[dict[str, str]]) -> None:
    datasets = [
        ("027", "A5", 5, 24, "e4_link/A5/final/target_summary.csv", "e4_link_gate/e4_link_manifest.json"),
        ("027", "A10", 10, 24, "e4_link/A10/final/target_summary.csv", "e4_link_gate/e4_link_manifest.json"),
        ("027", "A30", 30, 24, "e4_link/A30/final/target_summary.csv", "e4_link_gate/e4_link_manifest.json"),
        ("027", "A30", 30, 12, "e5_dense_dmrs/comb12/link/A30/final/target_summary.csv", "e5_dense_dmrs/comb12/manifest/e5_link_manifest.json"),
        ("027", "A30", 30, 6, "e5_dense_dmrs/comb6/link/A30/final/target_summary.csv", "e5_dense_dmrs/comb6/manifest/e5_link_manifest.json"),
        ("027", "A100", 100, 6, "e6_a100_dense_dmrs/comb6/link/A100/final/target_summary.csv", "e6_a100_dense_dmrs/comb6/manifest/e6_link_manifest.json"),
        ("028", "A300", 300, 6, "a300_comb6/comb6/link/A300/final/target_summary.csv", "a300_comb6/comb6/manifest/p28_link_manifest.json"),
    ]
    root027 = Path("outputs/experiment027_meff_sidon/20260726_main")
    root028 = Path("outputs/experiment028_csi_curves/20260803_main")
    for experiment, scenario, ds, comb, csv_rel, manifest_rel in datasets:
        root = root027 if experiment == "027" else root028
        index = manifest_index(root / manifest_rel)
        source_path = root / csv_rel
        for item in read_csv(source_path):
            candidate = index.get(item["candidate_id"], {})
            status = "bracketed" if item.get("target_snr_db") else "unbracketed"
            rows.append(
                pds_common(
                    experiment=experiment,
                    comparison_group=f"PDSCH_{scenario}_static_8T1R_comb{comb}_estimated",
                    result_source=f"result-{experiment}",
                    result_status="executed; awaiting researcher confirmation",
                    delay_spread_ns=ds,
                    carrier_frequency_ghz=3.5,
                    ue_speed_kmh=0,
                    n_tx=8,
                    n_rx=1,
                    csi_mode="estimated",
                    channel_estimator="two static DMRS observations averaged, then matched frequency LMMSE",
                    tx_power_normalization="precoder squared norm 8; noise variance 8/SNR",
                    snr_definition="total Tx power per RE divided by single-Rx noise variance (equivalent to unit-power convention)",
                    candidate_id=item["candidate_id"],
                    candidate_family=item.get("family", candidate.get("family", "")),
                    receiver_knowledge="non-transparent; matched effective covariance",
                    delay_coordinates=seq(candidate.get("delay_grid_coordinates")),
                    delay_ns=seq(candidate.get("delay_ns")),
                    precoder_description=f"CDD exp(-j2pi*k*j/576); DMRS comb-{comb}",
                    target_bler=item.get("target_probability", ""),
                    target_snr_db=number(item.get("target_snr_db")),
                    target_ci95_low_db=number(item.get("ci95_low_db")),
                    target_ci95_high_db=number(item.get("ci95_high_db")),
                    target_status=status,
                    target_method="local logistic fit used by result-027/028 target_summary",
                    data_source=str(source_path).replace("\\", "/"),
                    notes="A30/A100 data are reused by result-028 rather than independently rerun" if experiment == "027" and scenario in {"A30", "A100"} and comb == 6 else "",
                )
            )

    # Result-028 A100 supplemental transparent PRG and small-delay receiver comparisons.
    supplemental = [
        ("estimated", "outputs/experiment028_csi_curves/20260803_main/transparent_prg_baselines/a100/final/estimated_csi_bler_points.csv"),
        ("estimated", "outputs/experiment028_csi_curves/20260803_main/transparent_prg_6rb_1pct_extension/a100/final/estimated_csi_bler_points.csv"),
        ("ideal", "outputs/experiment028_csi_curves/20260803_main/transparent_prg_baselines/a100/final/ideal_csi_bler_points.csv"),
        ("estimated", "outputs/experiment028_csi_curves/20260803_main/transparent_cdd_small_delay/a100/final/estimated_csi_bler_points.csv"),
        ("estimated", "outputs/experiment028_csi_curves/20260803_main/small_delay_cdd_matched/a100/final/estimated_csi_bler_points.csv"),
        ("ideal", "outputs/experiment028_csi_curves/20260803_main/small_delay_cdd_ideal/a100/final/ideal_csi_bler_points.csv"),
    ]
    for csi_mode, path in supplemental:
        full = ROOT / path
        if not full.exists():
            continue
        raw = read_csv(path)
        by_candidate: dict[str, list[dict[str, str]]] = {}
        for item in raw:
            by_candidate.setdefault(item["candidate_id"], []).append(item)
        for candidate_id, points in by_candidate.items():
            if "PRG" in candidate_id:
                family = "PRECODER_CYCLING"
                knowledge = "transparent; physical covariance within each PRG"
                delays = ""
                precoder = "DFT8 precoder cycling by 4-RB or 6-RB PRG"
            else:
                family = "SMALL_DELAY_CDD"
                knowledge = "ideal CSI" if csi_mode == "ideal" else ("matched effective covariance" if "MATCHED" in candidate_id else "transparent; physical fullband covariance")
                delays = seq([0, 14.467593, 28.935185, 43.402778, 57.87037, 72.337963, 86.805556, 101.273148])
                precoder = "CDD coordinates [0,0.25,...,1.75], phase denominator 576"
            for target in (0.1, 0.01):
                crossing, bracket = local_log_crossing(points, target)
                rows.append(
                    pds_common(
                        experiment="028",
                        comparison_group=f"PDSCH_A100_static_8T1R_comb6_{csi_mode}",
                        result_source="result-028",
                        result_status="executed; awaiting researcher confirmation",
                        delay_spread_ns=100,
                        carrier_frequency_ghz=3.5,
                        ue_speed_kmh=0,
                        n_tx=8,
                        n_rx=1,
                        csi_mode=csi_mode,
                        channel_estimator="true data-RE effective channel" if csi_mode == "ideal" else "static frequency LMMSE; receiver covariance depends on candidate mode",
                        tx_power_normalization="precoder squared norm 8; noise variance 8/SNR",
                        snr_definition="total Tx power per RE divided by single-Rx noise variance (equivalent to unit-power convention)",
                        candidate_id=candidate_id,
                        candidate_family=family,
                        receiver_knowledge=knowledge,
                        delay_coordinates="[0,0.25,0.5,0.75,1,1.25,1.5,1.75]" if delays else "",
                        delay_ns=delays,
                        precoder_description=precoder,
                        target_bler=target,
                        target_snr_db=crossing,
                        target_status="bracketed" if crossing else "unbracketed",
                        target_method="audit recomputation: adjacent raw points, linear interpolation in log10(BLER)",
                        data_source=path,
                        notes=bracket,
                    )
                )


def yaml_candidate(config_path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    config = read_yaml(config_path)
    candidate = config.get("candidates", [{}])[0]
    return config, candidate


def delay_ns_from_config(config: dict[str, Any], candidate: dict[str, Any]) -> list[float] | None:
    coords = candidate.get("delay_grid_coordinates")
    if not coords:
        return None
    resource = config["resource"]
    al = int(resource["aggregation_level"])
    duration = int(resource["duration_symbols"])
    occupied_rb = 6 * al // duration
    k_active = occupied_rb * 12
    q_ns = 1e9 / (k_active * float(resource["scs_khz"]) * 1e3)
    return [float(v) * q_ns for v in coords]


def add_simple_pdcch_curve(
    rows: list[dict[str, str]],
    experiment: str,
    group: str,
    result_source: str,
    status: str,
    points_path: str,
    meta: dict[str, Any],
    candidate_id: str,
    family: str,
    precoder: str,
) -> None:
    points = [row for row in read_csv(points_path) if row.get("candidate_id", candidate_id) == candidate_id]
    for target in (0.1, 0.01):
        crossing, bracket = local_log_crossing(points, target)
        rows.append(
            pdcch_common(
                experiment=experiment,
                comparison_group=group,
                result_source=result_source,
                result_status=status,
                candidate_id=candidate_id,
                candidate_family=family,
                precoder_description=precoder,
                target_bler=target,
                target_snr_db=crossing,
                target_status="bracketed" if crossing else "unbracketed",
                target_method="audit recomputation: adjacent raw points, linear interpolation in log10(BLER)",
                data_source=points_path,
                notes=bracket,
                **meta,
            )
        )


def add_029_030_034(rows: list[dict[str, str]]) -> None:
    common = dict(
        channel_model="Sionna 1.0.2 TDL-A",
        delay_spread_ns=100,
        carrier_frequency_ghz=3.5,
        ue_speed_kmh=0,
        symbols=2,
        n_tx=8,
        n_rx=1,
        payload_bits=41,
        dmrs="PDCCH DMRS, 3 RE/REG; bundle-local across two symbols",
        csi_mode="estimated",
        channel_estimator="bundle LS across two symbols, then frequency-linear interpolation",
        receiver_knowledge="transparent to DFT index; bundle-local estimate",
    )
    add_simple_pdcch_curve(
        rows,
        "029",
        "PDCCH_A100_8T1R_2sym_AL4_est_A41",
        "result-029",
        "executed; awaiting researcher confirmation",
        "outputs/experiment029_pdcch/20260902_main/bler_points.csv",
        {**common, "occupied_rb": 12, "mcs_or_code": "A=41, CRC24C, Polar K=65 N=512 E=432"},
        "PDCCH_AL4_DFT8",
        "PRECODER_CYCLING",
        "unit-norm DFT8 bundle cycling order [0,1,2,3]",
    )
    add_simple_pdcch_curve(
        rows,
        "030",
        "PDCCH_A100_8T1R_2sym_AL8_est_A41",
        "result-030",
        "executed; awaiting researcher confirmation",
        "outputs/experiment030_pdcch_al8/20260903_main/bler_points.csv",
        {**common, "occupied_rb": 24, "mcs_or_code": "A=41, CRC24C, Polar K=65 N=512 E=864 with repetition"},
        "PDCCH_AL8_DFT8",
        "PRECODER_CYCLING",
        "unit-norm DFT8 bundle cycling order [0..7]",
    )
    # Result-034 ran a single fixed-DFT BLER curve and explicitly did not publish a target crossing.
    paths = [
        "outputs/experiment034_rsrp_pdcch_4t4r/20260916_plan034/bler/low/bler_points.csv",
        "outputs/experiment034_rsrp_pdcch_4t4r/20260916_plan034/bler/mid/bler_points.csv",
        "outputs/experiment034_rsrp_pdcch_4t4r/20260916_plan034/bler/high/bler_points.csv",
    ]
    points: list[dict[str, str]] = []
    for path in paths:
        points.extend(read_csv(path))
    unique = {(float(p["snr_db"]), p["candidate_id"]): p for p in points}
    merged_path = "outputs/experiment034_rsrp_pdcch_4t4r/20260916_plan034/analysis/bler_combined.csv"
    for target in (0.1, 0.01):
        crossing, bracket = local_log_crossing(list(unique.values()), target)
        rows.append(
            pdcch_common(
                experiment="034",
                comparison_group="PDCCH_C300_4T4R_2sym_AL1_est_A40",
                result_source="result-034",
                result_status="executed; awaiting researcher confirmation",
                channel_model="Sionna 1.0.2 TDL-C",
                delay_spread_ns=300,
                carrier_frequency_ghz=4,
                ue_speed_kmh=3,
                occupied_rb=3,
                symbols=2,
                n_tx=4,
                n_rx=4,
                payload_bits=40,
                mcs_or_code="A=40, CRC24C, Polar E=108",
                dmrs="comb-4 on both CORESET symbols; 18 DMRS RE",
                csi_mode="estimated",
                channel_estimator="per-Rx 2D time-frequency LMMSE, then coherent 4Rx MRC",
                candidate_id="C300_AL1_FIXED_DFT0_NR4_EST",
                candidate_family="FIXED_DFT",
                receiver_knowledge="transparent physical-PRG covariance; AL1 has one fixed DFT0 bundle",
                precoder_description="fixed DFT4 index 0, [1,1,1,1]/2",
                target_bler=target,
                target_snr_db=crossing,
                target_status="derived_in_audit" if crossing else "unbracketed",
                target_method="audit recomputation only; result-034 explicitly did not publish a threshold interpolation",
                data_source=merged_path,
                notes=bracket,
            )
        )


def add_031_table(
    rows: list[dict[str, str]],
    csv_path: str,
    experiment_tag: str,
    config_resolver,
    metadata_resolver,
) -> None:
    for item in read_csv(csv_path):
        meta = metadata_resolver(item)
        config_path = config_resolver(item)
        config: dict[str, Any] = {}
        candidate: dict[str, Any] = {}
        if config_path and (ROOT / config_path).exists():
            config, candidate = yaml_candidate(Path(config_path))
        target_value = item.get("target_snr_db", "")
        rows.append(
            pdcch_common(
                experiment="031",
                result_source="result-031",
                result_status="executed; awaiting researcher confirmation" if "partial" not in experiment_tag else "partial; AL4 incomplete",
                candidate_id=item["candidate_id"],
                candidate_family=candidate.get("family", meta.get("candidate_family", "")),
                receiver_knowledge=meta.get("receiver_knowledge", candidate.get("receiver_covariance_mode", "")),
                delay_coordinates=seq(candidate.get("delay_grid_coordinates")),
                delay_ns=seq(delay_ns_from_config(config, candidate)) if config else meta.get("delay_ns", ""),
                precoder_description=meta.get("precoder_description", candidate.get("scheme", "")),
                target_bler=item.get("target_bler", ""),
                target_snr_db=number(target_value),
                target_ci95_low_db=number(item.get("target_snr_ci95_lo_db") or item.get("ci95_lo_db")),
                target_ci95_high_db=number(item.get("target_snr_ci95_hi_db") or item.get("ci95_hi_db")),
                target_status="bracketed" if target_value else item.get("status", "unbracketed"),
                target_method=item.get("method", "Jeffreys-smoothed weighted decreasing isotonic; local log10-BLER interpolation"),
                data_source=csv_path,
                **{key: value for key, value in meta.items() if key in FIELDS and key not in {"candidate_family", "receiver_knowledge", "delay_ns", "precoder_description"}},
            )
        )


def add_031(rows: list[dict[str, str]]) -> None:
    base_root = "outputs/experiment031_pdcch_cdd/20260903_main"

    def base_meta(item: dict[str, str]) -> dict[str, Any]:
        al = int(item["aggregation_level"])
        return dict(
            comparison_group=f"PDCCH_A100_8T1R_1sym_AL{al}_estimated_A41",
            channel_model="Sionna 1.0.2 TDL-A",
            delay_spread_ns=100,
            carrier_frequency_ghz=3.5,
            ue_speed_kmh=0,
            occupied_rb=6 * al,
            symbols=1,
            n_tx=8,
            n_rx=1,
            payload_bits=41,
            mcs_or_code=f"A=41, CRC24C, Polar E={108*al}; AL8 repetition when E=864",
            dmrs=f"PDCCH DMRS, 3 RE/REG; {36*al} DMRS RE",
            csi_mode="estimated",
            channel_estimator="matched fullband frequency LMMSE for CDD; physical covariance per PRG for cycling",
        )

    add_031_table(
        rows,
        f"{base_root}/analysis/target_snr.csv",
        "base",
        lambda i: f"{base_root}/formal/al{i['aggregation_level']}/{i['candidate_id']}/resolved_config.yaml",
        base_meta,
    )

    strict_root = "outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search"

    def strict_meta(item: dict[str, str]) -> dict[str, Any]:
        profile, al_text = item["scene"].split("_al")
        al = int(al_text)
        a100 = profile == "a100"
        return dict(
            comparison_group=f"PDCCH_{'A100_8T1R_1sym' if a100 else 'C300_4T1R_2sym'}_AL{al}_estimated_A41",
            channel_model=f"Sionna 1.0.2 TDL-{'A' if a100 else 'C'}",
            delay_spread_ns=100 if a100 else 300,
            carrier_frequency_ghz=3.5 if a100 else 4,
            ue_speed_kmh=0 if a100 else 3,
            occupied_rb=6 * al if a100 else 3 * al,
            symbols=1 if a100 else 2,
            n_tx=8 if a100 else 4,
            n_rx=1,
            payload_bits=41,
            mcs_or_code=f"A=41, CRC24C, Polar E={108*al}",
            dmrs="PDCCH DMRS; comb-3/REG for 1-symbol A100 or comb-4 on both symbols for C300",
            csi_mode="estimated",
            channel_estimator="candidate-matched frequency or 2D time-frequency LMMSE",
            receiver_knowledge="non-transparent; matched effective covariance",
        )

    def strict_config(item: dict[str, str]) -> str:
        primary = ROOT / strict_root / "fine" / item["scene"] / item["candidate_id"] / "resolved_config.yaml"
        if primary.exists():
            return str(primary.relative_to(ROOT))
        matches = list((ROOT / strict_root).glob(f"**/{item['candidate_id']}/resolved_config.yaml"))
        return str(matches[0].relative_to(ROOT)) if matches else ""

    add_031_table(rows, f"{strict_root}/fine_analysis/target_snr.csv", "strict", strict_config, strict_meta)

    c300_root = "outputs/experiment031_pdcch_cdd/20260907_c300_4tx_2sym"

    def c300_meta(item: dict[str, str]) -> dict[str, Any]:
        al = int(item["aggregation_level"])
        return dict(
            comparison_group=f"PDCCH_C300_4T1R_2sym_AL{al}_estimated_A41",
            channel_model="Sionna 1.0.2 TDL-C",
            delay_spread_ns=300,
            carrier_frequency_ghz=4,
            ue_speed_kmh=3,
            occupied_rb=3 * al,
            symbols=2,
            n_tx=4,
            n_rx=1,
            payload_bits=41,
            mcs_or_code=f"A=41, CRC24C, Polar E={108*al}",
            dmrs="comb-4 on both CORESET symbols",
            csi_mode="estimated",
            channel_estimator="matched/physical 2D time-frequency LMMSE; per-PRG for cycling",
        )

    add_031_table(
        rows,
        f"{c300_root}/analysis/target_snr.csv",
        "c300",
        lambda i: f"{c300_root}/formal/al{i['aggregation_level']}/{i['candidate_id']}/resolved_config.yaml",
        c300_meta,
    )

    ideal_root = "outputs/experiment031_pdcch_cdd/20260910_c300_4tx_2sym_ideal_csi"

    def ideal_meta(item: dict[str, str]) -> dict[str, Any]:
        meta = c300_meta(item)
        al = int(item["aggregation_level"])
        meta.update(
            comparison_group=f"PDCCH_C300_4T1R_2sym_AL{al}_ideal_A41",
            csi_mode="ideal",
            channel_estimator="true data-RE effective channel",
            receiver_knowledge="ideal CSI",
        )
        return meta

    add_031_table(
        rows,
        f"{ideal_root}/analysis/target_snr.csv",
        "ideal",
        lambda i: f"{ideal_root}/formal/al{i['aggregation_level']}/{i['candidate_id']}/resolved_config.yaml",
        ideal_meta,
    )

    rx2_root = "outputs/experiment031_pdcch_cdd/20260911_c300_4tx_2rx_2sym"

    def rx2_meta(item: dict[str, str]) -> dict[str, Any]:
        al = int(item["aggregation_level"])
        csi = item["csi_mode"]
        return dict(
            comparison_group=f"PDCCH_C300_4T2R_2sym_AL{al}_{csi}_A41",
            channel_model="Sionna 1.0.2 TDL-C; independent Rx branches",
            delay_spread_ns=300,
            carrier_frequency_ghz=4,
            ue_speed_kmh=3,
            occupied_rb=3 * al,
            symbols=2,
            n_tx=4,
            n_rx=2,
            payload_bits=41,
            mcs_or_code=f"A=41, CRC24C, Polar E={108*al}",
            dmrs="comb-4 on both CORESET symbols",
            csi_mode=csi,
            channel_estimator="true data-RE channel" if csi == "ideal" else "independent per-Rx 2D time-frequency LMMSE, then coherent MRC",
        )

    add_031_table(
        rows,
        f"{rx2_root}/analysis/partial_al1_al2/target_snr.csv",
        "partial_2rx",
        lambda i: f"{rx2_root}/formal/{i['csi_mode']}/al{i['aggregation_level']}/{i['candidate_id']}/resolved_config.yaml",
        rx2_meta,
    )

    supplement = "outputs/experiment031_pdcch_cdd/20260914_c300_al1_cdd911/analysis/target_snr.csv"
    for item in read_csv(supplement):
        csi = item["receiver"]
        n_rx = int(item["n_rx"])
        cid = item["candidate_id"]
        delay_map = {
            "CDD911": [0, 0, 911, 911],
            "CDD911_transparent": [0, 0, 911, 911],
            "CDD130_transparent": [0, 0, 130, 130],
            "CDD130": [0, 0, 130, 130],
            "C300_S0_SIDON": [0, 925.925926, 2777.777778, 6481.481481],
            "C300_S0_SIDON_TRANSPARENT": [0, 925.925926, 2777.777778, 6481.481481],
        }
        value = item.get("target_snr_db", "")
        rows.append(
            pdcch_common(
                experiment="031",
                comparison_group=f"PDCCH_C300_4T{n_rx}R_2sym_AL1_{csi}_A41",
                result_source="result-031",
                result_status="executed; awaiting researcher confirmation",
                channel_model="Sionna 1.0.2 TDL-C; independent Rx branches",
                delay_spread_ns=300,
                carrier_frequency_ghz=4,
                ue_speed_kmh=3,
                occupied_rb=3,
                symbols=2,
                n_tx=4,
                n_rx=n_rx,
                payload_bits=41,
                mcs_or_code="A=41, CRC24C, Polar E=108",
                dmrs="comb-4 on both CORESET symbols; 18 DMRS RE",
                csi_mode=csi,
                channel_estimator="true data-RE channel" if csi == "ideal" else "per-Rx 2D time-frequency LMMSE, then coherent MRC",
                candidate_id=cid,
                candidate_family="FIXED_CODEBOOK" if "codebook" in cid.lower() else ("PRECODER_CYCLING" if "cycling" in item["curve_id"].lower() else "CDD"),
                receiver_knowledge="ideal CSI" if csi == "ideal" else ("transparent physical covariance" if "transparent" in cid.lower() or "codebook" in cid.lower() else "non-transparent matched covariance"),
                delay_ns=seq(delay_map.get(cid)),
                precoder_description="fixed codebook/cycling vector" if cid not in delay_map else "CDD with explicit physical delays",
                target_bler=item["target_bler"],
                target_snr_db=number(value),
                target_ci95_low_db=number(item.get("target_snr_ci95_lo_db")),
                target_ci95_high_db=number(item.get("target_snr_ci95_hi_db")),
                target_status="bracketed" if value else item.get("status", "unbracketed"),
                target_method=item.get("method", "Jeffreys/isotonic/local log interpolation"),
                data_source=supplement,
                notes="Rows with blank target failed the <=0.25 dB bracket rule and are not extrapolated.",
            )
        )


def add_032_033_035(rows: list[dict[str, str]]) -> None:
    path032 = "outputs/experiment032_tdl_mobility/20260906_main/final/target_summary.csv"
    manifest032 = read_json("outputs/experiment032_tdl_mobility/20260906_main/formal/resolved_run.json")
    index032 = {item["candidate_id"]: item for item in manifest032["candidates"]}
    for item in read_csv(path032):
        candidate = index032.get(item["candidate_id"], {})
        value = item.get("target_snr_db", "") if item.get("bracketed", "").lower() == "true" else ""
        rows.append(
            pds_common(
                experiment="032",
                comparison_group="PDSCH_A100_60kmh_8T1R_comb6_estimated",
                result_source="result-032",
                result_status="executed; awaiting researcher confirmation",
                delay_spread_ns=100,
                carrier_frequency_ghz=3.5,
                ue_speed_kmh=60,
                n_tx=8,
                n_rx=1,
                csi_mode="estimated",
                channel_estimator="two-DMRS joint 2D time-frequency RMMSE; matched CDD or physical per-PRG covariance",
                tx_power_normalization="precoder squared norm 8; noise variance 8/SNR",
                snr_definition="total Tx power per RE divided by single-Rx noise variance",
                candidate_id=item["candidate_id"],
                candidate_family=candidate.get("family", ""),
                receiver_knowledge="matched CDD covariance or transparent PRG covariance as identified by candidate",
                delay_coordinates=seq(candidate.get("delay_grid_coordinates")),
                delay_ns=seq(candidate.get("delay_ns")),
                precoder_description="CDD exp(-j2pi*k*j/576), DFT8 PRG6 cycling, or 5ms-aged PRG MRT",
                target_bler=item["target_bler"],
                target_snr_db=number(value),
                target_status="bracketed" if value else "unbracketed",
                target_method="adjacent raw bracket; linear interpolation in log10(BLER); no smoothing/extrapolation",
                data_source=path032,
            )
        )

    path033 = "outputs/experiment033_tdl_mobility_mimo/20260915_main/analysis_partial/target_summary_initial.csv"
    manifest_cache: dict[str, dict[str, Any]] = {}
    for item in read_csv(path033):
        scenario = item["scenario_id"]
        if scenario not in manifest_cache:
            path = f"outputs/experiment033_tdl_mobility_mimo/20260915_main/formal/{scenario}/candidate_receiver_manifest.json"
            manifest_cache[scenario] = manifest_index(path)
        candidate = manifest_cache[scenario].get(item["candidate_id"], {})
        n_tx = int(candidate.get("n_tx", scenario.split("_NT")[1].split("_")[0]))
        n_rx = int(candidate.get("n_rx", 4))
        speed = float(candidate.get("speed_kmh", scenario.rsplit("V", 1)[1]))
        value = item.get("crossing_snr_db", "") if item.get("status") == "bracketed" else ""
        rows.append(
            pds_common(
                experiment="033",
                comparison_group=f"PDSCH_A100_{number(speed)}kmh_{n_tx}T{n_rx}R_comb6_estimated",
                result_source="result-033",
                result_status="partial formal experiment; endpoint append and paired bootstrap incomplete",
                delay_spread_ns=100,
                carrier_frequency_ghz=3.5,
                ue_speed_kmh=number(speed),
                n_tx=n_tx,
                n_rx=n_rx,
                csi_mode="estimated",
                channel_estimator="independent per-Rx two-DMRS 2D time-frequency RMMSE, then coherent MRC",
                tx_power_normalization="unit squared-norm precoder; noise variance 1/SNR per Rx branch",
                snr_definition="total Tx power per RE divided by each Rx branch noise variance; no n_rx scaling",
                candidate_id=item["candidate_id"],
                candidate_family=candidate.get("family", ""),
                receiver_knowledge=candidate.get("receiver_knowledge", ""),
                delay_coordinates=seq(candidate.get("delay_grid_coordinates")),
                delay_ns=seq([float(v) * 1e9 / (576 * 30000) for v in candidate.get("delay_grid_coordinates", [])]),
                precoder_description="unit-power CDD, DFT PRG6 cycling, or 5ms-aged PRG MRT",
                target_bler=item["target_bler"],
                target_snr_db=number(value),
                target_status="preliminary_bracketed" if value else item.get("status", "unbracketed"),
                target_method="initial adjacent raw bracket; linear interpolation in log10(BLER)",
                data_source=path033,
                notes="Do not treat as final: result-033 reports insufficient endpoint errors and incomplete append/bootstrap.",
            )
        )

    root035 = "outputs/experiment035_pdsch_2rx/20260917_fix1/formal"
    for scenario in ("A100_NT4_NR2_V60", "A100_NT8_NR2_V60"):
        path = f"{root035}/{scenario}/analysis/target_crossings_bootstrap.csv"
        index = manifest_index(f"{root035}/{scenario}/candidate_receiver_manifest.json")
        for item in read_csv(path):
            candidate = index[item["candidate_id"]]
            n_tx = int(candidate["n_tx"])
            csi = item["receiver"]
            rows.append(
                pds_common(
                    experiment="035",
                    comparison_group=f"PDSCH_A100_60kmh_{n_tx}T2R_comb6_{csi}",
                    result_source="result-035",
                    result_status="executed; awaiting researcher confirmation",
                    delay_spread_ns=100,
                    carrier_frequency_ghz=3.5,
                    ue_speed_kmh=60,
                    n_tx=n_tx,
                    n_rx=2,
                    csi_mode=csi,
                    channel_estimator="true data-RE channel" if csi == "ideal" else "independent per-Rx two-DMRS 2D time-frequency RMMSE, then coherent MRC",
                    tx_power_normalization="unit squared-norm precoder; noise variance 1/SNR per Rx branch",
                    snr_definition="total Tx power per RE divided by each Rx branch noise variance; no n_rx scaling",
                    candidate_id=item["candidate_id"],
                    candidate_family=candidate.get("family", ""),
                    receiver_knowledge=candidate.get("receiver_knowledge", ""),
                    delay_coordinates=seq(candidate.get("delay_grid_coordinates")),
                    delay_ns=seq([float(v) * 1e9 / (576 * 30000) for v in candidate.get("delay_grid_coordinates", [])]),
                    precoder_description="unit-power CDD, DFT PRG6 cycling, or 5ms-aged PRG MRT",
                    target_bler=item["target_bler"],
                    target_snr_db=number(item["crossing_snr_db"]),
                    target_ci95_low_db=number(item["crossing_ci95_lo_db"]),
                    target_ci95_high_db=number(item["crossing_ci95_hi_db"]),
                    target_status="bracketed",
                    target_method="adjacent raw bracket, log10(BLER) interpolation; 500 paired bootstrap replicates",
                    data_source=path,
                )
            )


def add_036(rows: list[dict[str, str]]) -> None:
    datasets = [
        (4, 1, "outputs/experiment036_pdcch_4t4r_1symbol/20260918_formal/al1_fd5/bler_points.csv"),
        (4, 2, "outputs/experiment036_pdcch_4t4r_1symbol/20260918_formal/al2_fd5/bler_points.csv"),
        (2, 1, "outputs/experiment036_pdcch_4t2r_1symbol/20260919_formal/al1_fd5/bler_points.csv"),
        (2, 2, "outputs/experiment036_pdcch_4t2r_1symbol/20260919_formal/al2_fd5/bler_points.csv"),
    ]
    delay_coords = {
        1: {"SIDON": [0, 1, 3, 7], "QC_DELAY_NT": [0, .24192, .48384, .72576], "QC_DELAY_TRANSPARENT": [0, .24192, .48384, .72576]},
        2: {"SIDON": [0, 11, 19, 64], "QC_DELAY_NT": [0, .41472, .82944, 1.24416], "QC_DELAY_TRANSPARENT": [0, .41472, .82944, 1.24416]},
    }
    for n_rx, al, path in datasets:
        raw = read_csv(path)
        by_curve: dict[tuple[str, str], list[dict[str, str]]] = {}
        for item in raw:
            by_curve.setdefault((item["csi_mode"], item["candidate_id"]), []).append(item)
        occupied_rb = 6 * al
        q_ns = 1e9 / (occupied_rb * 12 * 30000)
        for (csi, candidate_id), points in by_curve.items():
            coords = delay_coords[al].get(candidate_id)
            for target in (0.1, 0.01):
                crossing, bracket = local_log_crossing(points, target)
                rows.append(
                    pdcch_common(
                        experiment="036",
                        comparison_group=f"PDCCH_C300_4T{n_rx}R_1sym_AL{al}_{csi}_A40_fd5",
                        result_source="result-036",
                        result_status="5 Hz subset executed; 1100 Hz and paired bootstrap not executed",
                        channel_model="Sionna 1.0.2 TDL-C",
                        delay_spread_ns=300,
                        carrier_frequency_ghz=4,
                        ue_speed_kmh=1.349066061,
                        occupied_rb=occupied_rb,
                        symbols=1,
                        n_tx=4,
                        n_rx=n_rx,
                        payload_bits=40,
                        mcs_or_code=f"A=40, CRC24C, Polar E={108*al}",
                        dmrs="single-symbol comb-4 PDCCH DMRS",
                        csi_mode=csi,
                        channel_estimator="true data-RE channel" if csi == "ideal" else "per-Rx frequency LMMSE, then coherent MRC",
                        candidate_id=candidate_id,
                        candidate_family="PRECODER_CYCLING" if "DFT" in candidate_id else "CDD",
                        receiver_knowledge="ideal CSI" if csi == "ideal" else ("transparent physical covariance" if "TRANSPARENT" in candidate_id or "DFT" in candidate_id else "non-transparent matched covariance"),
                        delay_coordinates=seq(coords),
                        delay_ns=seq([v * q_ns for v in coords]) if coords else "",
                        precoder_description="DFT4 cycling" if "DFT" in candidate_id else "CDD with phase denominator equal to occupied subcarriers",
                        target_bler=target,
                        target_snr_db=crossing,
                        target_status="bracketed" if crossing else "unbracketed",
                        target_method="adjacent raw 1 dB points, log10(BLER) interpolation; point estimate only",
                        data_source=path,
                        notes=bracket,
                    )
                )


def assign_baselines(rows: list[dict[str, str]]) -> None:
    groups: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        groups.setdefault((row["comparison_group"], row["target_bler"]), []).append(row)
    for group_rows in groups.values():
        candidates = [row for row in group_rows if row["target_snr_db"]]
        def is_cycling(row: dict[str, str]) -> bool:
            token = (row["candidate_id"] + " " + row["candidate_family"]).upper()
            return row["candidate_family"] == "PRECODER_CYCLING" or "PRG_DFT" in token or "TRANSPARENT_PRG" in token or "CYCLING" in token

        all_cycling = [row for row in group_rows if is_cycling(row)]
        cycling = [row for row in all_cycling if row["target_snr_db"]]
        if all_cycling:
            baseline = sorted(
                cycling or all_cycling,
                key=lambda row: ("PRG6" not in row["candidate_id"].upper() and "6RB" not in row["candidate_id"].upper(), row["candidate_id"]),
            )[0]
            rule = "precoder cycling present" if baseline["target_snr_db"] else "precoder cycling present, but its target is unavailable"
        else:
            all_b0 = [row for row in group_rows if "b0" in (row["candidate_id"] + row["candidate_family"]).lower()]
            b0 = [row for row in all_b0 if row["target_snr_db"]]
            baseline = (b0 or all_b0)[0] if (b0 or all_b0) else None
            rule = ("B0QC because no precoder cycling is present" if baseline and baseline["target_snr_db"] else "B0QC present, but its target is unavailable") if baseline else "NA: neither precoder cycling nor B0QC is present"
        for row in group_rows:
            row["baseline_rule"] = rule
            row["gain_definition"] = "baseline target SNR - candidate target SNR; positive means candidate needs less SNR"
            if baseline:
                row["baseline_candidate_id"] = baseline["candidate_id"]
                if row["target_snr_db"] and baseline["target_snr_db"]:
                    row["gain_vs_baseline_db"] = number(float(baseline["target_snr_db"]) - float(row["target_snr_db"]))
                    if row["target_method"] != baseline["target_method"]:
                        warning = "gain combines source target estimates produced by different crossing methods"
                        row["notes"] = f"{row['notes']}; {warning}" if row["notes"] else warning


def deduplicate(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    # Later supplemental result-031 tables deliberately reuse earlier formal curves.
    # Preserve one row per exact group/candidate/target and concatenate provenance.
    merged: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in rows:
        key = (row["comparison_group"], row["candidate_id"], row["target_bler"])
        if key not in merged:
            merged[key] = row
            continue
        old = merged[key]
        if not old["target_snr_db"] and row["target_snr_db"]:
            merged[key] = row
            old = row
        sources = {value for value in (old["data_source"], row["data_source"]) if value}
        old["data_source"] = ";".join(sorted(sources))
        notes = {value for value in (old["notes"], row["notes"], "duplicate curve reused; not an independent rerun") if value}
        old["notes"] = "; ".join(sorted(notes))
    return list(merged.values())


def main() -> None:
    rows: list[dict[str, str]] = []
    add_027_028(rows)
    add_029_030_034(rows)
    add_031(rows)
    add_032_033_035(rows)
    add_036(rows)
    rows = deduplicate(rows)
    assign_baselines(rows)
    rows.sort(key=lambda row: (int(row["experiment"]), row["comparison_group"], float(row["target_bler"] or 9), row["candidate_id"]))
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows to {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
