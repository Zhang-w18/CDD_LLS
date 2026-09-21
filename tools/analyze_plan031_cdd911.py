"""Audit and analyze the plan-031 AL1 CDD911 supplement."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import platform
import subprocess

import numpy as np
import yaml

import analyze_plan031 as common
import analyze_plan031_c300 as c300


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs" / "experiment031_pdcch_cdd" / "20260914_c300_al1_cdd911"
ANALYSIS = OUTPUT / "analysis"
FIGURES = ROOT / "docs" / "figures" / "result-031"
TARGETS = (0.10, 0.01)
NEW_CONFIGS = {
    ("estimated", 1): ROOT / "configs" / "pdcch_result031_cdd911_1rx_estimated_formal.yaml",
    ("estimated", 2): ROOT / "configs" / "pdcch_result031_cdd911_2rx_estimated_formal.yaml",
    ("ideal", 1): ROOT / "configs" / "pdcch_result031_cdd911_1rx_ideal_formal.yaml",
    ("ideal", 2): ROOT / "configs" / "pdcch_result031_cdd911_2rx_ideal_formal.yaml",
}
TRANSPARENT_CONFIGS = {
    1: ROOT / "configs" / "pdcch_result031_cdd_transparent_1rx_formal.yaml",
    2: ROOT / "configs" / "pdcch_result031_cdd_transparent_2rx_formal.yaml",
}
TRANSPARENT_DELAYS = {
    "CDD911_transparent": np.asarray([0.0, 0.0, 0.98388, 0.98388]),
    "CDD130_transparent": np.asarray([0.0, 0.0, 0.1404, 0.1404]),
}
SUPPLEMENT_20260915_CONFIGS = (
    (
        "estimated",
        1,
        ROOT / "configs" / "pdcch_result031_sidon_transparent_1rx_formal.yaml",
        "C300_S0_SIDON_TRANSPARENT",
    ),
    (
        "ideal",
        1,
        ROOT / "configs" / "pdcch_result031_cdd130_1rx_ideal_formal.yaml",
        "CDD130",
    ),
    (
        "ideal",
        2,
        ROOT / "configs" / "pdcch_result031_cdd130_2rx_ideal_formal.yaml",
        "CDD130",
    ),
)
SUPPLEMENT_20260915_DELAYS = {
    "C300_S0_SIDON_TRANSPARENT": np.asarray([0.0, 1.0, 3.0, 7.0]),
    "CDD130": np.asarray([0.0, 0.0, 0.1404, 0.1404]),
}
CODEBOOK_CONFIGS = {
    "estimated": ROOT / "configs" / "pdcch_result031_codebooks_2rx_estimated_formal.yaml",
    "ideal": ROOT / "configs" / "pdcch_result031_codebooks_2rx_ideal_formal.yaml",
}
CODEBOOK_TARGET_SNRS = (5.0, 6.0, 7.0, 7.5, 8.0, 8.5, 9.0, 10.0)
HISTORICAL_DFT_2RX_POINTS = {
    receiver: ROOT
    / "outputs"
    / "experiment031_pdcch_cdd"
    / "20260911_c300_4tx_2rx_2sym"
    / "formal"
    / receiver
    / "al1"
    / "C300_PRG_DFT4_6REG"
    / "bler_points.csv"
    for receiver in ("estimated", "ideal")
}
HISTORICAL_POINTS = {
    "estimated": ROOT
    / "outputs"
    / "experiment031_pdcch_cdd"
    / "20260907_c300_4tx_2sym"
    / "analysis"
    / "formal_points.csv",
    "ideal": ROOT
    / "outputs"
    / "experiment031_pdcch_cdd"
    / "20260910_c300_4tx_2sym_ideal_csi"
    / "analysis"
    / "formal_points.csv",
}
HISTORICAL_2RX_POINTS = (
    ROOT
    / "outputs"
    / "experiment031_pdcch_cdd"
    / "20260911_c300_4tx_2rx_2sym"
    / "analysis"
    / "partial_al1_al2"
    / "formal_points.csv"
)
BASELINES = ("C300_S0_SIDON", "C300_PRG_DFT4_6REG")
CURVE_STYLES = {
    "CDD911_1Rx": {
        "label": "CDD911, 1Rx",
        "estimated_label": "CDD911 non-transparent, 1Rx",
        "ideal_label": "CDD911 ideal-CSI, 1Rx",
        "color": "#d62728",
        "linestyle": "-",
        "marker": "o",
    },
    "CDD911_2Rx": {
        "label": "CDD911, 2Rx",
        "estimated_label": "CDD911 non-transparent, 2Rx",
        "ideal_label": "CDD911 ideal-CSI, 2Rx",
        "color": "#d62728",
        "linestyle": "--",
        "marker": "s",
    },
    "CDD911_transparent_1Rx": {
        "label": "CDD911 transparent, 1Rx",
        "color": "#d62728",
        "linestyle": ":",
        "marker": "P",
    },
    "CDD911_transparent_2Rx": {
        "label": "CDD911 transparent, 2Rx",
        "color": "#d62728",
        "linestyle": "-.",
        "marker": "X",
    },
    "CDD130_transparent_1Rx": {
        "label": "CDD130 transparent, 1Rx",
        "color": "#ff7f0e",
        "linestyle": ":",
        "marker": "P",
    },
    "CDD130_transparent_2Rx": {
        "label": "CDD130 transparent, 2Rx",
        "color": "#ff7f0e",
        "linestyle": "-.",
        "marker": "X",
    },
    "C300_S0_SIDON_1Rx": {
        "label": "Sidon, 1Rx",
        "estimated_label": "Sidon non-transparent, 1Rx",
        "ideal_label": "Sidon ideal-CSI, 1Rx",
        "color": "#1f77b4",
        "linestyle": "-.",
        "marker": "^",
    },
    "C300_S0_SIDON_2Rx": {
        "label": "Sidon, 2Rx",
        "estimated_label": "Sidon non-transparent, 2Rx",
        "ideal_label": "Sidon ideal-CSI, 2Rx",
        "color": "#1f77b4",
        "linestyle": "--",
        "marker": "v",
    },
    "C300_PRG_DFT4_6REG_1Rx": {
        "label": "Precoder cycling, 1Rx",
        "estimated_label": "Precoder cycling transparent, 1Rx",
        "ideal_label": "Precoder cycling ideal-CSI, 1Rx",
        "color": "#2ca02c",
        "linestyle": ":",
        "marker": "D",
    },
    "C300_S0_SIDON_TRANSPARENT_1Rx": {
        "label": "Sidon transparent, 1Rx",
        "color": "#1f77b4",
        "linestyle": ":",
        "marker": "P",
    },
    "CDD130_1Rx": {
        "label": "CDD130 ideal-CSI, 1Rx",
        "color": "#ff7f0e",
        "linestyle": "-",
        "marker": "X",
    },
    "CDD130_2Rx": {
        "label": "CDD130 ideal-CSI, 2Rx",
        "color": "#ff7f0e",
        "linestyle": "--",
        "marker": "X",
    },
    "DFTcodebook_2Rx": {
        "label": "DFTcodebook, 2Rx",
        "estimated_label": "DFTcodebook transparent, 2Rx",
        "ideal_label": "DFTcodebook ideal-CSI, 2Rx",
        "color": "#2ca02c",
        "linestyle": "--",
        "marker": ">",
    },
    "LTEcodebook_2Rx": {
        "label": "LTEcodebook, 2Rx",
        "estimated_label": "LTEcodebook transparent, 2Rx",
        "ideal_label": "LTEcodebook ideal-CSI, 2Rx",
        "color": "#9467bd",
        "linestyle": "-.",
        "marker": "D",
    },
}


def _read_csv(path: Path) -> list[dict[str, str]]:
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _read_yaml(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def _relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def _normalized_point(
    row: dict[str, str], *, receiver: str, n_rx: int, curve_id: str, source_csv: Path
) -> dict[str, object]:
    out: dict[str, object] = dict(row)
    out.pop("csi_mode_label", None)
    for key in (
        "snr_db",
        "bler",
        "bler_wilson95_lo",
        "bler_wilson95_hi",
        "ce_nmse_db",
        "pilot_condition_number",
        "ce_floor_nmse_db",
    ):
        out[key] = float(row[key])
    for key in (
        "trials",
        "errors",
        "aggregation_level",
        "occupied_rb",
        "k_active",
        "data_re",
        "dmrs_re",
        "coded_bits",
        "pilot_rank",
    ):
        out[key] = int(float(row[key]))
    out["receiver"] = receiver
    out["n_rx"] = n_rx
    out["curve_id"] = curve_id
    out["source_csv"] = _relative(source_csv)
    candidate_id = str(row["candidate_id"])
    if candidate_id == "CDD911":
        out["source_kind"] = "new_cdd911"
    elif candidate_id in TRANSPARENT_DELAYS:
        out["source_kind"] = "new_transparent_supplement"
    elif candidate_id in SUPPLEMENT_20260915_DELAYS:
        out["source_kind"] = "new_20260915_supplement"
    elif candidate_id in ("DFTcodebook", "LTEcodebook"):
        out["source_kind"] = "new_20260915_codebook_supplement"
    else:
        out["source_kind"] = "historical_reuse"
    return out


def _wilson(errors: int, trials: int) -> tuple[float, float]:
    z = 1.959963984540054
    p = errors / trials
    scale = 1.0 + z * z / trials
    center = (p + z * z / (2.0 * trials)) / scale
    half = z * math.sqrt(p * (1.0 - p) / trials + z * z / (4.0 * trials * trials)) / scale
    return center - half, center + half


def _audit_source_config(
    config: dict, *, receiver: str, n_rx: int, candidate_id: str, is_new: bool
) -> None:
    expected_estimator = "frequency_lmmse" if receiver == "estimated" else "ideal"
    assert config["schema"] == "pdcch-cdd-bler-v2"
    assert int(config["antenna"]["n_tx"]) == 4
    assert int(config["antenna"]["n_rx"]) == n_rx
    assert int(config["resource"]["duration_symbols"]) == 2
    assert int(config["resource"]["aggregation_level"]) == 1
    assert int(config["pdcch"]["coded_bits"]) == 108
    assert str(config["channel"]["tdl_profile"]).upper() == "C"
    assert math.isclose(float(config["channel"]["delay_spread_ns"]), 300.0)
    assert math.isclose(float(config["channel"]["carrier_frequency_hz"]), 4.0e9)
    assert math.isclose(float(config["channel"]["ue_speed_kmh"]), 3.0)
    assert str(config["receiver"]["channel_estimation"]).lower() == expected_estimator
    candidates = [item for item in config["candidates"] if str(item["candidate_id"]) == candidate_id]
    assert len(candidates) == 1
    if is_new:
        candidate = candidates[0]
        expected_duplicates = candidate_id != "C300_S0_SIDON_TRANSPARENT"
        assert candidate.get("allow_duplicate_delays", False) is expected_duplicates
        expected_mode = (
            "physical_fullband"
            if candidate_id in TRANSPARENT_DELAYS
            or candidate_id == "C300_S0_SIDON_TRANSPARENT"
            else "matched_effective"
        )
        expected_delays = TRANSPARENT_DELAYS.get(candidate_id)
        if expected_delays is None:
            expected_delays = SUPPLEMENT_20260915_DELAYS.get(
                candidate_id, np.asarray([0.0, 0.0, 0.98388, 0.98388])
            )
        assert str(candidate["receiver_covariance_mode"]) == expected_mode
        assert np.allclose(
            np.asarray(candidate["delay_grid_coordinates"], dtype=np.float64),
            expected_delays,
            rtol=0.0,
            atol=1e-12,
        )


def _load_new_points() -> tuple[list[dict[str, object]], list[dict[str, object]], list[Path]]:
    points: list[dict[str, object]] = []
    sources: list[dict[str, object]] = []
    configs_used: list[Path] = []
    for (receiver, n_rx), config_path in NEW_CONFIGS.items():
        config = _read_yaml(config_path)
        _audit_source_config(
            config, receiver=receiver, n_rx=n_rx, candidate_id="CDD911", is_new=True
        )
        candidate_dir = ROOT / str(config["output_dir"]) / "CDD911"
        csv_path = candidate_dir / "bler_points.csv"
        resolved_path = candidate_dir / "resolved_config.yaml"
        metadata_path = candidate_dir / "run_metadata.json"
        if not csv_path.exists() or not resolved_path.exists() or not metadata_path.exists():
            raise FileNotFoundError(f"Incomplete CDD911 shard: {candidate_dir}")
        resolved = _read_yaml(resolved_path)
        _audit_source_config(
            resolved, receiver=receiver, n_rx=n_rx, candidate_id="CDD911", is_new=True
        )
        rows = _read_csv(csv_path)
        expected_snr = sorted(float(value) for value in config["simulation"]["snr_points_db"])
        actual_snr = sorted(float(row["snr_db"]) for row in rows)
        if actual_snr != expected_snr or len(actual_snr) != len(set(actual_snr)):
            raise ValueError(f"Incomplete or duplicate SNR grid: {csv_path}")
        curve_id = f"CDD911_{n_rx}Rx"
        points.extend(
            _normalized_point(
                row, receiver=receiver, n_rx=n_rx, curve_id=curve_id, source_csv=csv_path
            )
            for row in rows
        )
        sources.append(
            {
                "receiver": receiver,
                "n_rx": n_rx,
                "candidate_id": "CDD911",
                "curve_id": curve_id,
                "source_csv": _relative(csv_path),
                "resolved_config": _relative(resolved_path),
                "run_metadata": _relative(metadata_path),
                "config_sha256": common.sha256(config_path),
                "source_kind": "new_cdd911",
            }
        )
        configs_used.append(config_path)
    return points, sources, configs_used


def _load_transparent_points() -> tuple[
    list[dict[str, object]], list[dict[str, object]], list[Path]
]:
    points: list[dict[str, object]] = []
    sources: list[dict[str, object]] = []
    configs_used: list[Path] = []
    for n_rx, config_path in TRANSPARENT_CONFIGS.items():
        config = _read_yaml(config_path)
        expected_snr = sorted(float(value) for value in config["simulation"]["snr_points_db"])
        for candidate_id in TRANSPARENT_DELAYS:
            _audit_source_config(
                config,
                receiver="estimated",
                n_rx=n_rx,
                candidate_id=candidate_id,
                is_new=True,
            )
            candidate_dir = ROOT / str(config["output_dir"]) / candidate_id
            csv_path = candidate_dir / "bler_points.csv"
            resolved_path = candidate_dir / "resolved_config.yaml"
            metadata_path = candidate_dir / "run_metadata.json"
            if not csv_path.exists() or not resolved_path.exists() or not metadata_path.exists():
                raise FileNotFoundError(f"Incomplete transparent shard: {candidate_dir}")
            resolved = _read_yaml(resolved_path)
            _audit_source_config(
                resolved,
                receiver="estimated",
                n_rx=n_rx,
                candidate_id=candidate_id,
                is_new=True,
            )
            rows = _read_csv(csv_path)
            actual_snr = sorted(float(row["snr_db"]) for row in rows)
            if actual_snr != expected_snr or len(actual_snr) != len(set(actual_snr)):
                raise ValueError(f"Incomplete or duplicate SNR grid: {csv_path}")
            curve_id = f"{candidate_id}_{n_rx}Rx"
            points.extend(
                _normalized_point(
                    row,
                    receiver="estimated",
                    n_rx=n_rx,
                    curve_id=curve_id,
                    source_csv=csv_path,
                )
                for row in rows
            )
            sources.append(
                {
                    "receiver": "estimated",
                    "n_rx": n_rx,
                    "candidate_id": candidate_id,
                    "curve_id": curve_id,
                    "source_csv": _relative(csv_path),
                    "resolved_config": _relative(resolved_path),
                    "run_metadata": _relative(metadata_path),
                    "config_sha256": common.sha256(config_path),
                    "source_kind": "new_transparent_supplement",
                }
            )
        configs_used.append(config_path)
    return points, sources, configs_used


def _load_20260915_points() -> tuple[
    list[dict[str, object]], list[dict[str, object]], list[Path]
]:
    points: list[dict[str, object]] = []
    sources: list[dict[str, object]] = []
    configs_used: list[Path] = []
    for receiver, n_rx, config_path, candidate_id in SUPPLEMENT_20260915_CONFIGS:
        config = _read_yaml(config_path)
        _audit_source_config(
            config,
            receiver=receiver,
            n_rx=n_rx,
            candidate_id=candidate_id,
            is_new=True,
        )
        candidate_dir = ROOT / str(config["output_dir"]) / candidate_id
        csv_path = candidate_dir / "bler_points.csv"
        resolved_path = candidate_dir / "resolved_config.yaml"
        metadata_path = candidate_dir / "run_metadata.json"
        if not csv_path.exists() or not resolved_path.exists() or not metadata_path.exists():
            raise FileNotFoundError(f"Incomplete 20260915 supplement shard: {candidate_dir}")
        resolved = _read_yaml(resolved_path)
        _audit_source_config(
            resolved,
            receiver=receiver,
            n_rx=n_rx,
            candidate_id=candidate_id,
            is_new=True,
        )
        rows = _read_csv(csv_path)
        expected_snr = sorted(float(value) for value in config["simulation"]["snr_points_db"])
        actual_snr = sorted(float(row["snr_db"]) for row in rows)
        if actual_snr != expected_snr or len(actual_snr) != len(set(actual_snr)):
            raise ValueError(f"Incomplete or duplicate SNR grid: {csv_path}")
        curve_id = f"{candidate_id}_{n_rx}Rx"
        points.extend(
            _normalized_point(
                row,
                receiver=receiver,
                n_rx=n_rx,
                curve_id=curve_id,
                source_csv=csv_path,
            )
            for row in rows
        )
        sources.append(
            {
                "receiver": receiver,
                "n_rx": n_rx,
                "candidate_id": candidate_id,
                "curve_id": curve_id,
                "source_csv": _relative(csv_path),
                "resolved_config": _relative(resolved_path),
                "run_metadata": _relative(metadata_path),
                "config_sha256": common.sha256(config_path),
                "source_kind": "new_20260915_supplement",
            }
        )
        configs_used.append(config_path)
    return points, sources, configs_used


def _audit_codebook_config(
    config: dict, *, receiver: str, candidate_id: str, expected_snrs: list[float]
) -> None:
    _audit_source_config(
        config,
        receiver=receiver,
        n_rx=2,
        candidate_id=candidate_id,
        is_new=False,
    )
    candidate = next(
        item for item in config["candidates"] if str(item["candidate_id"]) == candidate_id
    )
    assert str(candidate["scheme"]).lower() == "reg_bundle_dft_cycling"
    assert str(candidate["receiver_covariance_mode"]).lower() == "physical_prg"
    assert [int(value) for value in candidate["cycling_order"]] == [
        0 if candidate_id == "DFTcodebook" else 2
    ]
    assert sorted(float(value) for value in candidate["snr_points_db"]) == sorted(
        expected_snrs
    )


def _load_codebook_points() -> tuple[
    list[dict[str, object]], list[dict[str, object]], list[Path]
]:
    points: list[dict[str, object]] = []
    sources: list[dict[str, object]] = []
    configs_used: list[Path] = []
    target_snrs = set(CODEBOOK_TARGET_SNRS)
    for receiver, config_path in CODEBOOK_CONFIGS.items():
        config = _read_yaml(config_path)
        new_by_candidate: dict[str, list[dict[str, str]]] = {}
        for candidate_id in ("DFTcodebook", "LTEcodebook"):
            candidate = next(
                item
                for item in config["candidates"]
                if str(item["candidate_id"]) == candidate_id
            )
            expected_new_snrs = [float(value) for value in candidate["snr_points_db"]]
            _audit_codebook_config(
                config,
                receiver=receiver,
                candidate_id=candidate_id,
                expected_snrs=expected_new_snrs,
            )
            candidate_dir = ROOT / str(config["output_dir"]) / candidate_id
            csv_path = candidate_dir / "bler_points.csv"
            resolved_path = candidate_dir / "resolved_config.yaml"
            metadata_path = candidate_dir / "run_metadata.json"
            if not csv_path.exists() or not resolved_path.exists() or not metadata_path.exists():
                raise FileNotFoundError(f"Incomplete codebook shard: {candidate_dir}")
            resolved = _read_yaml(resolved_path)
            _audit_codebook_config(
                resolved,
                receiver=receiver,
                candidate_id=candidate_id,
                expected_snrs=expected_new_snrs,
            )
            rows = _read_csv(csv_path)
            actual_snrs = sorted(float(row["snr_db"]) for row in rows)
            if actual_snrs != sorted(expected_new_snrs) or len(actual_snrs) != len(
                set(actual_snrs)
            ):
                raise ValueError(f"Incomplete or duplicate codebook SNR grid: {csv_path}")
            new_by_candidate[candidate_id] = rows
            curve_id = f"{candidate_id}_2Rx"
            points.extend(
                _normalized_point(
                    row,
                    receiver=receiver,
                    n_rx=2,
                    curve_id=curve_id,
                    source_csv=csv_path,
                )
                for row in rows
            )
            sources.append(
                {
                    "receiver": receiver,
                    "n_rx": 2,
                    "candidate_id": candidate_id,
                    "curve_id": curve_id,
                    "source_csv": _relative(csv_path),
                    "resolved_config": _relative(resolved_path),
                    "run_metadata": _relative(metadata_path),
                    "config_sha256": common.sha256(config_path),
                    "source_kind": "new_20260915_codebook_supplement",
                }
            )

        historical_csv = HISTORICAL_DFT_2RX_POINTS[receiver]
        historical_rows = [
            row
            for row in _read_csv(historical_csv)
            if float(row["snr_db"]) in target_snrs
            and float(row["snr_db"])
            not in {float(item["snr_db"]) for item in new_by_candidate["DFTcodebook"]}
        ]
        historical_dir = historical_csv.parent
        historical_resolved = historical_dir / "resolved_config.yaml"
        historical_metadata = historical_dir / "run_metadata.json"
        _audit_source_config(
            _read_yaml(historical_resolved),
            receiver=receiver,
            n_rx=2,
            candidate_id="C300_PRG_DFT4_6REG",
            is_new=False,
        )
        points.extend(
            _normalized_point(
                row,
                receiver=receiver,
                n_rx=2,
                curve_id="DFTcodebook_2Rx",
                source_csv=historical_csv,
            )
            for row in historical_rows
        )
        sources.append(
            {
                "receiver": receiver,
                "n_rx": 2,
                "candidate_id": "C300_PRG_DFT4_6REG",
                "curve_id": "DFTcodebook_2Rx",
                "source_csv": _relative(historical_csv),
                "resolved_config": _relative(historical_resolved),
                "run_metadata": _relative(historical_metadata),
                "source_kind": "historical_dft_2rx_reuse",
            }
        )
        for curve_id in ("DFTcodebook_2Rx", "LTEcodebook_2Rx"):
            actual = sorted(
                float(row["snr_db"])
                for row in points
                if row["receiver"] == receiver and row["curve_id"] == curve_id
            )
            if actual != list(CODEBOOK_TARGET_SNRS):
                raise ValueError(f"Incomplete combined grid for {receiver}/{curve_id}: {actual}")
        configs_used.append(config_path)
    return points, sources, configs_used


def _load_historical_points() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    points: list[dict[str, object]] = []
    sources: list[dict[str, object]] = []
    for receiver, summary_path in HISTORICAL_POINTS.items():
        selected = [
            row
            for row in _read_csv(summary_path)
            if int(row["aggregation_level"]) == 1 and row["candidate_id"] in BASELINES
        ]
        for candidate_id in BASELINES:
            candidate_rows = [row for row in selected if row["candidate_id"] == candidate_id]
            if not candidate_rows:
                raise ValueError(f"Missing historical {receiver} {candidate_id}")
            curve_id = f"{candidate_id}_1Rx"
            for row in candidate_rows:
                csv_path = ROOT / row["source_csv"]
                points.append(
                    _normalized_point(
                        row,
                        receiver=receiver,
                        n_rx=1,
                        curve_id=curve_id,
                        source_csv=csv_path,
                    )
                )
            source_paths = sorted({str(row["source_csv"]) for row in candidate_rows})
            for source in source_paths:
                csv_path = ROOT / source
                candidate_dir = csv_path.parent
                resolved_path = candidate_dir / "resolved_config.yaml"
                metadata_path = candidate_dir / "run_metadata.json"
                resolved = _read_yaml(resolved_path)
                _audit_source_config(
                    resolved,
                    receiver=receiver,
                    n_rx=1,
                    candidate_id=candidate_id,
                    is_new=False,
                )
                sources.append(
                    {
                        "receiver": receiver,
                        "n_rx": 1,
                        "candidate_id": candidate_id,
                        "curve_id": curve_id,
                        "source_csv": _relative(csv_path),
                        "resolved_config": _relative(resolved_path),
                        "run_metadata": _relative(metadata_path),
                        "source_kind": "historical_reuse",
                    }
                )
    two_rx_rows = _read_csv(HISTORICAL_2RX_POINTS)
    for receiver, csi_mode in (("estimated", "frequency_lmmse"), ("ideal", "ideal")):
        selected = [
            row
            for row in two_rx_rows
            if int(row["aggregation_level"]) == 1
            and row["candidate_id"] == "C300_S0_SIDON"
            and str(row["csi_mode"]).lower() == csi_mode
            and int(row["n_rx"]) == 2
        ]
        if not selected:
            raise ValueError(f"Missing historical 2Rx Sidon {receiver} points")
        curve_id = "C300_S0_SIDON_2Rx"
        for row in selected:
            csv_path = ROOT / row["source_csv"]
            points.append(
                _normalized_point(
                    row,
                    receiver=receiver,
                    n_rx=2,
                    curve_id=curve_id,
                    source_csv=csv_path,
                )
            )
        source_paths = sorted({str(row["source_csv"]) for row in selected})
        for source in source_paths:
            csv_path = ROOT / source
            candidate_dir = csv_path.parent
            resolved_path = candidate_dir / "resolved_config.yaml"
            metadata_path = candidate_dir / "run_metadata.json"
            resolved = _read_yaml(resolved_path)
            _audit_source_config(
                resolved,
                receiver=receiver,
                n_rx=2,
                candidate_id="C300_S0_SIDON",
                is_new=False,
            )
            sources.append(
                {
                    "receiver": receiver,
                    "n_rx": 2,
                    "candidate_id": "C300_S0_SIDON",
                    "curve_id": curve_id,
                    "source_csv": _relative(csv_path),
                    "resolved_config": _relative(resolved_path),
                    "run_metadata": _relative(metadata_path),
                    "source_summary": _relative(HISTORICAL_2RX_POINTS),
                    "source_kind": "historical_2rx_reuse",
                }
            )
    return points, sources


def _audit_points(points: list[dict[str, object]]) -> list[dict[str, object]]:
    diagnostics: list[dict[str, object]] = []
    for receiver in ("estimated", "ideal"):
        curve_ids = [
            "CDD911_1Rx",
            "CDD911_2Rx",
            "C300_S0_SIDON_1Rx",
            "C300_S0_SIDON_2Rx",
            "C300_PRG_DFT4_6REG_1Rx",
        ]
        if receiver == "estimated":
            curve_ids.extend(
                [
                    "CDD911_transparent_1Rx",
                    "CDD911_transparent_2Rx",
                    "CDD130_transparent_1Rx",
                    "CDD130_transparent_2Rx",
                    "C300_S0_SIDON_TRANSPARENT_1Rx",
                ]
            )
        else:
            curve_ids.extend(["CDD130_1Rx", "CDD130_2Rx"])
        curve_ids.extend(["DFTcodebook_2Rx", "LTEcodebook_2Rx"])
        for curve_id in curve_ids:
            rows = sorted(
                [
                    row
                    for row in points
                    if row["receiver"] == receiver and row["curve_id"] == curve_id
                ],
                key=lambda row: float(row["snr_db"]),
            )
            if not rows:
                raise ValueError(f"Missing curve {receiver}/{curve_id}")
            max_flag_error = 0
            max_wilson_error = 0.0
            max_nmse_error_db = 0.0
            increases = 0
            previous = None
            for row in rows:
                if int(row["aggregation_level"]) != 1:
                    raise ValueError("Non-AL1 point entered CDD911 analysis.")
                if previous is not None and float(row["bler"]) > previous + 1e-15:
                    increases += 1
                previous = float(row["bler"])
                csv_path = ROOT / str(row["source_csv"])
                candidate_dir = csv_path.parent
                stem = common._flag_stem(str(row["candidate_id"]), float(row["snr_db"]))
                flags = np.load(candidate_dir / "trial_error_flags" / f"{stem}_error_flags.npy")
                ce = np.load(candidate_dir / "trial_error_flags" / f"{stem}_ce_nmse.npy")
                max_flag_error = max(
                    max_flag_error,
                    abs(len(flags) - int(row["trials"])),
                    abs(int(np.sum(flags)) - int(row["errors"])),
                )
                lo, hi = _wilson(int(row["errors"]), int(row["trials"]))
                max_wilson_error = max(
                    max_wilson_error,
                    abs(lo - float(row["bler_wilson95_lo"])),
                    abs(hi - float(row["bler_wilson95_hi"])),
                )
                ce_db = 10.0 * math.log10(max(float(np.mean(ce)), 1e-30))
                max_nmse_error_db = max(
                    max_nmse_error_db, abs(ce_db - float(row["ce_nmse_db"]))
                )
            diagnostics.append(
                {
                    "receiver": receiver,
                    "curve_id": curve_id,
                    "candidate_id": rows[0]["candidate_id"],
                    "n_rx": rows[0]["n_rx"],
                    "formal_points": len(rows),
                    "total_trials": sum(int(row["trials"]) for row in rows),
                    "total_errors": sum(int(row["errors"]) for row in rows),
                    "adjacent_bler_increases": increases,
                    "max_flag_count_or_sum_error": max_flag_error,
                    "max_wilson_recompute_error": max_wilson_error,
                    "max_ce_recompute_error_db": max_nmse_error_db,
                    "reached_max_without_200_errors": sum(
                        int(row["trials"]) == 50000 and int(row["errors"]) < 200
                        for row in rows
                    ),
                    "pilot_rank": rows[0]["pilot_rank"],
                    "pilot_condition_number": rows[0]["pilot_condition_number"],
                    "ce_floor_nmse_db": rows[0]["ce_floor_nmse_db"],
                }
            )
    return diagnostics


def _targets(
    points: list[dict[str, object]], repeats: int
) -> tuple[list[dict[str, object]], dict[tuple[str, str, float], np.ndarray]]:
    rows_out: list[dict[str, object]] = []
    replicate_map: dict[tuple[str, str, float], np.ndarray] = {}
    for receiver in ("estimated", "ideal"):
        for curve_id in CURVE_STYLES:
            rows = sorted(
                [
                    row
                    for row in points
                    if row["receiver"] == receiver and row["curve_id"] == curve_id
                ],
                key=lambda row: float(row["snr_db"]),
            )
            if not rows:
                continue
            for target in TARGETS:
                key = (receiver, curve_id, target)
                raw_bracket_lo = raw_bracket_hi = None
                for left, right in zip(rows, rows[1:]):
                    left_delta = float(left["bler"]) - target
                    right_delta = float(right["bler"]) - target
                    if left_delta * right_delta <= 0.0:
                        raw_bracket_lo = float(left["snr_db"])
                        raw_bracket_hi = float(right["snr_db"])
                        break
                try:
                    bracket_lo, bracket_hi = c300.qualified_raw_bracket(rows, target)
                    estimate, ci_lo, ci_hi, replicates = common.bootstrap_targets(
                        rows,
                        target,
                        repeats,
                        common.stable_seed(
                            20260914, receiver, curve_id, target, "cdd911-bootstrap"
                        ),
                    )
                    status = "estimated"
                    method = "Jeffreys-smoothed weighted decreasing isotonic curve; local log10-BLER interpolation; raw adjacent bracket <=0.25 dB"
                except ValueError as error:
                    bracket_lo = raw_bracket_lo
                    bracket_hi = raw_bracket_hi
                    estimate = ci_lo = ci_hi = None
                    replicates = np.asarray([], dtype=np.float64)
                    status = "unqualified"
                    method = f"not estimated: {error}"
                replicate_map[key] = replicates
                rows_out.append(
                    {
                        "receiver": receiver,
                        "curve_id": curve_id,
                        "candidate_id": rows[0]["candidate_id"],
                        "n_rx": rows[0]["n_rx"],
                        "target_bler": target,
                        "target_snr_db": estimate,
                        "target_snr_ci95_lo_db": ci_lo,
                        "target_snr_ci95_hi_db": ci_hi,
                        "raw_bracket_lo_db": bracket_lo,
                        "raw_bracket_hi_db": bracket_hi,
                        "raw_bracket_width_db": (
                            None if bracket_lo is None else float(bracket_hi) - float(bracket_lo)
                        ),
                        "status": status,
                        "method": method,
                        "ci_method": (
                            "independent per-point Bernoulli bootstrap"
                            if status == "estimated"
                            else "not applicable"
                        ),
                        "bootstrap_repeats_requested": repeats,
                        "bootstrap_repeats_valid": len(replicates),
                    }
                )
    return rows_out, replicate_map


def _differences(
    targets: list[dict[str, object]],
    replicates: dict[tuple[str, str, float], np.ndarray],
) -> list[dict[str, object]]:
    lookup = {
        (str(row["receiver"]), str(row["curve_id"]), float(row["target_bler"])): row
        for row in targets
    }
    comparisons: list[tuple[str, str, str, str, str]] = []
    for receiver in ("estimated", "ideal"):
        for baseline in ("C300_S0_SIDON_1Rx", "C300_PRG_DFT4_6REG_1Rx"):
            comparisons.append(
                (
                    f"gain_cdd911_vs_{baseline.removesuffix('_1Rx')}",
                    receiver,
                    baseline,
                    receiver,
                    "CDD911_1Rx",
                )
            )
        comparisons.append(
            ("delta_rx_cdd911", receiver, "CDD911_1Rx", receiver, "CDD911_2Rx")
        )
    for n_rx in (1, 2):
        comparisons.append(
            (
                f"delta_ce_cdd911_{n_rx}rx",
                "estimated",
                f"CDD911_{n_rx}Rx",
                "ideal",
                f"CDD911_{n_rx}Rx",
            )
        )
        comparisons.extend(
            [
                (
                    f"delta_covariance_cdd911_{n_rx}rx",
                    "estimated",
                    f"CDD911_transparent_{n_rx}Rx",
                    "estimated",
                    f"CDD911_{n_rx}Rx",
                ),
                (
                    f"delta_delay_transparent_{n_rx}rx",
                    "estimated",
                    f"CDD911_transparent_{n_rx}Rx",
                    "estimated",
                    f"CDD130_transparent_{n_rx}Rx",
                ),
            ]
        )
    for candidate_id in TRANSPARENT_DELAYS:
        comparisons.append(
            (
                f"delta_rx_{candidate_id.lower()}",
                "estimated",
                f"{candidate_id}_1Rx",
                "estimated",
                f"{candidate_id}_2Rx",
            )
        )
    for n_rx in (1, 2):
        comparisons.append(
            (
                f"delta_ce_cdd130_{n_rx}rx",
                "estimated",
                f"CDD130_transparent_{n_rx}Rx",
                "ideal",
                f"CDD130_{n_rx}Rx",
            )
        )
    comparisons.append(
        (
            "delta_rx_cdd130_ideal",
            "ideal",
            "CDD130_1Rx",
            "ideal",
            "CDD130_2Rx",
        )
    )
    for receiver in ("estimated", "ideal"):
        comparisons.append(
            (
                f"delta_codebook_2rx_{receiver}",
                receiver,
                "DFTcodebook_2Rx",
                receiver,
                "LTEcodebook_2Rx",
            )
        )
    out: list[dict[str, object]] = []
    for comparison_id, left_receiver, left_curve, right_receiver, right_curve in comparisons:
        for target in TARGETS:
            left_key = (left_receiver, left_curve, target)
            right_key = (right_receiver, right_curve, target)
            left_rep = replicates[left_key]
            right_rep = replicates[right_key]
            count = min(len(left_rep), len(right_rep))
            if count:
                samples = left_rep[:count] - right_rep[:count]
                lo, hi = (float(value) for value in np.quantile(samples, [0.025, 0.975]))
                value = float(lookup[left_key]["target_snr_db"]) - float(
                    lookup[right_key]["target_snr_db"]
                )
                status = "estimated"
            else:
                value = lo = hi = None
                status = "unavailable_unqualified_target"
            out.append(
                {
                    "comparison_id": comparison_id,
                    "target_bler": target,
                    "left_receiver": left_receiver,
                    "left_curve_id": left_curve,
                    "right_receiver": right_receiver,
                    "right_curve_id": right_curve,
                    "left_minus_right_db": value,
                    "ci95_lo_db": lo,
                    "ci95_hi_db": hi,
                    "status": status,
                    "bootstrap_repeats_valid": count,
                }
            )
    return out


def _same_snr(points: list[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for receiver in ("estimated", "ideal"):
        cdd_rows = {
            float(row["snr_db"]): row
            for row in points
            if row["receiver"] == receiver and row["curve_id"] == "CDD911_1Rx"
        }
        for baseline_curve in ("C300_S0_SIDON_1Rx", "C300_PRG_DFT4_6REG_1Rx"):
            baseline_rows = {
                float(row["snr_db"]): row
                for row in points
                if row["receiver"] == receiver and row["curve_id"] == baseline_curve
            }
            for snr in sorted(set(cdd_rows) & set(baseline_rows)):
                cdd_bler = float(cdd_rows[snr]["bler"])
                baseline_bler = float(baseline_rows[snr]["bler"])
                out.append(
                    {
                        "receiver": receiver,
                        "baseline_curve_id": baseline_curve,
                        "snr_db": snr,
                        "cdd911_bler": cdd_bler,
                        "baseline_bler": baseline_bler,
                        "cdd911_minus_baseline_bler": cdd_bler - baseline_bler,
                        "cdd911_over_baseline_bler": (
                            cdd_bler / baseline_bler if baseline_bler > 0.0 else None
                        ),
                    }
                )
    return out


def _same_snr_cdd911(points: list[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    comparisons = [
        ("rx_effect_estimated", "estimated", "CDD911_1Rx", "estimated", "CDD911_2Rx"),
        ("rx_effect_ideal", "ideal", "CDD911_1Rx", "ideal", "CDD911_2Rx"),
        ("ce_effect_1rx", "estimated", "CDD911_1Rx", "ideal", "CDD911_1Rx"),
        ("ce_effect_2rx", "estimated", "CDD911_2Rx", "ideal", "CDD911_2Rx"),
        (
            "covariance_effect_1rx",
            "estimated",
            "CDD911_transparent_1Rx",
            "estimated",
            "CDD911_1Rx",
        ),
        (
            "covariance_effect_2rx",
            "estimated",
            "CDD911_transparent_2Rx",
            "estimated",
            "CDD911_2Rx",
        ),
        (
            "delay_effect_transparent_1rx",
            "estimated",
            "CDD911_transparent_1Rx",
            "estimated",
            "CDD130_transparent_1Rx",
        ),
        (
            "delay_effect_transparent_2rx",
            "estimated",
            "CDD911_transparent_2Rx",
            "estimated",
            "CDD130_transparent_2Rx",
        ),
        (
            "ce_effect_cdd130_1rx",
            "estimated",
            "CDD130_transparent_1Rx",
            "ideal",
            "CDD130_1Rx",
        ),
        (
            "ce_effect_cdd130_2rx",
            "estimated",
            "CDD130_transparent_2Rx",
            "ideal",
            "CDD130_2Rx",
        ),
        (
            "covariance_effect_sidon_1rx",
            "estimated",
            "C300_S0_SIDON_TRANSPARENT_1Rx",
            "estimated",
            "C300_S0_SIDON_1Rx",
        ),
        (
            "codebook_effect_estimated_2rx",
            "estimated",
            "DFTcodebook_2Rx",
            "estimated",
            "LTEcodebook_2Rx",
        ),
        (
            "codebook_effect_ideal_2rx",
            "ideal",
            "DFTcodebook_2Rx",
            "ideal",
            "LTEcodebook_2Rx",
        ),
    ]
    for comparison_id, left_receiver, left_curve, right_receiver, right_curve in comparisons:
        left_rows = {
            float(row["snr_db"]): row
            for row in points
            if row["receiver"] == left_receiver and row["curve_id"] == left_curve
        }
        right_rows = {
            float(row["snr_db"]): row
            for row in points
            if row["receiver"] == right_receiver and row["curve_id"] == right_curve
        }
        for snr in sorted(set(left_rows) & set(right_rows)):
            left_bler = float(left_rows[snr]["bler"])
            right_bler = float(right_rows[snr]["bler"])
            out.append(
                {
                    "comparison_id": comparison_id,
                    "snr_db": snr,
                    "left_receiver": left_receiver,
                    "left_curve_id": left_curve,
                    "left_bler": left_bler,
                    "right_receiver": right_receiver,
                    "right_curve_id": right_curve,
                    "right_bler": right_bler,
                    "left_minus_right_bler": left_bler - right_bler,
                    "left_over_right_bler": (
                        left_bler / right_bler if right_bler > 0.0 else None
                    ),
                }
            )
    return out


def _plot(points: list[dict[str, object]]) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator

    FIGURES.mkdir(parents=True, exist_ok=True)
    for receiver in ("estimated", "ideal"):
        fig, axis = plt.subplots(figsize=(17.0, 9.0) if receiver == "estimated" else (15.0, 8.5))
        for curve_id, style in CURVE_STYLES.items():
            rows = sorted(
                [
                    row
                    for row in points
                    if row["receiver"] == receiver and row["curve_id"] == curve_id
                ],
                key=lambda row: float(row["snr_db"]),
            )
            if not rows:
                continue
            x = np.asarray([float(row["snr_db"]) for row in rows])
            y = np.asarray([max(float(row["bler"]), 0.5 / int(row["trials"])) for row in rows])
            lo = np.asarray([float(row["bler_wilson95_lo"]) for row in rows])
            hi = np.asarray([float(row["bler_wilson95_hi"]) for row in rows])
            axis.errorbar(
                x,
                y,
                yerr=np.vstack((y - lo, hi - y)),
                label=style.get(f"{receiver}_label", style["label"]),
                color=style["color"],
                linestyle=style["linestyle"],
                marker=style["marker"],
                linewidth=2.5,
                markersize=8,
                capsize=2.5,
                elinewidth=1.3,
            )
        axis.axhline(0.10, color="black", linestyle=":", linewidth=1.8, label="10% BLER")
        axis.axhline(0.01, color="black", linestyle="--", linewidth=1.8, label="1% BLER")
        axis.set_yscale("log")
        axis.set_ylim(1e-4, 1.0)
        axis.set_xlabel("SNR (dB)", fontsize=16)
        axis.set_ylabel(f"{receiver.capitalize()}-CSI DCI BLER", fontsize=16)
        axis.set_title(f"C300 4Tx two-symbol PDCCH AL1: {receiver}-CSI", fontsize=17)
        axis.xaxis.set_major_locator(MultipleLocator(1.0))
        axis.xaxis.set_minor_locator(MultipleLocator(0.5))
        axis.tick_params(labelsize=14)
        axis.grid(True, which="major", axis="x", alpha=0.42, linewidth=0.9)
        axis.grid(True, which="minor", axis="x", alpha=0.22, linewidth=0.55)
        axis.grid(True, which="both", axis="y", alpha=0.35)
        axis.legend(fontsize=16, loc="center left", bbox_to_anchor=(1.01, 0.5))
        fig.tight_layout()
        figure_path = FIGURES / f"c300_4tx_2sym_al1_cdd911_{receiver}_bler.png"
        fig.savefig(
            figure_path,
            dpi=180,
            bbox_inches="tight",
        )
        plt.close(fig)
        from PIL import Image

        preview_dir = ANALYSIS / "previews_13cm"
        preview_dir.mkdir(parents=True, exist_ok=True)
        with Image.open(figure_path) as source:
            preview_width_px = 1535
            preview_height_px = round(source.height * preview_width_px / source.width)
            preview = source.resize((preview_width_px, preview_height_px), Image.Resampling.LANCZOS)
            preview.save(preview_dir / figure_path.name)


def _environment() -> dict[str, object]:
    import matplotlib
    import scipy
    import sionna
    import tensorflow

    return {
        "platform": platform.platform(),
        "python": platform.python_version(),
        "packages": {
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
            "tensorflow": tensorflow.__version__,
            "sionna": sionna.__version__,
            "pyyaml": yaml.__version__,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bootstrap-repeats", type=int, default=4000)
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Regenerate figures without rewriting analysis tables or metadata.",
    )
    args = parser.parse_args()
    new_points, new_sources, configs_used = _load_new_points()
    transparent_points, transparent_sources, transparent_configs = _load_transparent_points()
    supplement_points, supplement_sources, supplement_configs = _load_20260915_points()
    codebook_points, codebook_sources, codebook_configs = _load_codebook_points()
    historical_points, historical_sources = _load_historical_points()
    points = (
        new_points
        + transparent_points
        + supplement_points
        + codebook_points
        + historical_points
    )
    if args.plot_only:
        _plot(points)
        return
    diagnostics = _audit_points(points)
    targets, replicate_map = _targets(points, int(args.bootstrap_repeats))
    differences = _differences(targets, replicate_map)
    same_snr = _same_snr(points)
    same_snr_cdd911 = _same_snr_cdd911(points)
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    common.write_csv(ANALYSIS / "formal_points.csv", points)
    common.write_csv(ANALYSIS / "diagnostics.csv", diagnostics)
    common.write_csv(ANALYSIS / "target_snr.csv", targets)
    common.write_csv(ANALYSIS / "target_differences.csv", differences)
    common.write_csv(ANALYSIS / "same_snr_comparisons.csv", same_snr)
    common.write_csv(ANALYSIS / "same_snr_cdd911_comparisons.csv", same_snr_cdd911)
    (ANALYSIS / "curve_styles.json").write_text(
        json.dumps(CURVE_STYLES, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    sources = (
        new_sources
        + transparent_sources
        + supplement_sources
        + codebook_sources
        + historical_sources
    )
    (ANALYSIS / "source_receipt.json").write_text(
        json.dumps(sources, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    _plot(points)
    git_head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True
    ).stdout.strip()
    metadata = {
        "schema": "plan031-cdd911-analysis-v4",
        "plan": "research/plan-031-PDCCH-CDD时延-BLER.md (2026-09-14/15 AL1 CDD supplements)",
        "bootstrap_repeats": int(args.bootstrap_repeats),
        "git_head": git_head,
        "working_tree_note": "CDD911, transparent CDD, Sidon transparent, CDD130 ideal, and codebook supplements are uncommitted",
        "analysis_script_sha256": common.sha256(Path(__file__)),
        "config_sha256": {
            _relative(path): common.sha256(path)
            for path in configs_used
            + transparent_configs
            + supplement_configs
            + codebook_configs
        },
        "formal_point_count": len(points),
        "new_cdd911_point_count": len(new_points),
        "historical_reused_point_count": len(historical_points),
        "total_trials_all_plotted_curves": sum(int(row["trials"]) for row in points),
        "total_errors_all_plotted_curves": sum(int(row["errors"]) for row in points),
        "new_cdd911_trials": sum(int(row["trials"]) for row in new_points),
        "new_cdd911_errors": sum(int(row["errors"]) for row in new_points),
        "new_transparent_point_count": len(transparent_points),
        "new_transparent_trials": sum(int(row["trials"]) for row in transparent_points),
        "new_transparent_errors": sum(int(row["errors"]) for row in transparent_points),
        "new_20260915_point_count": len(supplement_points),
        "new_20260915_trials": sum(int(row["trials"]) for row in supplement_points),
        "new_20260915_errors": sum(int(row["errors"]) for row in supplement_points),
        "new_codebook_point_count": sum(
            row["source_kind"] == "new_20260915_codebook_supplement"
            for row in codebook_points
        ),
        "historical_dft_codebook_reused_point_count": sum(
            row["source_kind"] == "historical_reuse"
            for row in codebook_points
        ),
        "new_codebook_trials": sum(
            int(row["trials"])
            for row in codebook_points
            if row["source_kind"] == "new_20260915_codebook_supplement"
        ),
        "new_codebook_errors": sum(
            int(row["errors"])
            for row in codebook_points
            if row["source_kind"] == "new_20260915_codebook_supplement"
        ),
        "new_20260915_curves_reaching_raw_bler_le_1pct": sum(
            min(
                float(row["bler"])
                for row in supplement_points
                if row["curve_id"] == curve_id
            )
            <= 0.01
            for curve_id in (
                "C300_S0_SIDON_TRANSPARENT_1Rx",
                "CDD130_1Rx",
                "CDD130_2Rx",
            )
        ),
        "transparent_curves_reaching_raw_bler_le_1pct": sum(
            min(
                float(row["bler"])
                for row in transparent_points
                if row["curve_id"] == curve_id
            )
            <= 0.01
            for curve_id in (
                "CDD911_transparent_1Rx",
                "CDD911_transparent_2Rx",
                "CDD130_transparent_1Rx",
                "CDD130_transparent_2Rx",
            )
        ),
        "qualified_target_rows": sum(row["status"] == "estimated" for row in targets),
        "target_rows": len(targets),
        "max_flag_count_or_sum_error": max(
            int(row["max_flag_count_or_sum_error"]) for row in diagnostics
        ),
        "max_wilson_recompute_error": max(
            float(row["max_wilson_recompute_error"]) for row in diagnostics
        ),
        "max_ce_recompute_error_db": max(
            float(row["max_ce_recompute_error_db"]) for row in diagnostics
        ),
    }
    (ANALYSIS / "analysis_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (OUTPUT / "environment_receipt.json").write_text(
        json.dumps(_environment(), indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
