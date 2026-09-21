from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import yaml
from scipy.stats import ks_2samp

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.run_plan034_rsrp_cdf import (  # noqa: E402
    PER_RX_CANDIDATES,
    PER_RX_FIELDNAMES,
    PER_RX_METADATA,
)


DEFAULT_ROOT = ROOT / "outputs" / "experiment034_rsrp_pdcch_4t4r" / "20260916_plan034"
BLER_CONFIG_NAMES = {
    "low": "pdcch_plan034_fixed_dft0_4rx_bler_low.yaml",
    "mid": "pdcch_plan034_fixed_dft0_4rx_bler_mid.yaml",
    "high": "pdcch_plan034_fixed_dft0_4rx_bler_high.yaml",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    fields = list(rows[0])
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _quantile(values: np.ndarray, probability: float) -> float:
    return float(np.quantile(values, probability, method="inverted_cdf"))


def _summary(values: np.ndarray) -> dict[str, float]:
    q01 = _quantile(values, 0.01)
    q10 = _quantile(values, 0.10)
    q50 = _quantile(values, 0.50)
    q90 = _quantile(values, 0.90)
    return {
        "q01_db": q01,
        "q10_db": q10,
        "q50_db": q50,
        "q90_db": q90,
        "central80_width_db": q90 - q10,
        "mean_db": float(np.mean(values)),
        "std_db": float(np.std(values, ddof=1)),
    }


def _bootstrap_quantiles(values: np.ndarray) -> tuple[float, float, float]:
    n = len(values)
    indices = tuple(max(int(np.ceil(probability * n)) - 1, 0) for probability in (0.01, 0.10, 0.90))
    partitioned = np.partition(values, indices)
    return tuple(float(partitioned[index]) for index in indices)


def paired_bootstrap(
    fixed: np.ndarray,
    freq: np.ndarray,
    *,
    repeats: int,
    seed: int,
) -> dict[str, dict[str, float]]:
    if fixed.shape != freq.shape or fixed.ndim != 1:
        raise ValueError("Paired bootstrap inputs must be aligned one-dimensional arrays.")
    rng = np.random.default_rng(int(seed))
    n = len(fixed)
    draws = {"delta_q01_db": [], "delta_q10_db": [], "delta_b80_db": []}
    for _ in range(int(repeats)):
        indices = rng.integers(0, n, size=n)
        fixed_sample = fixed[indices]
        freq_sample = freq[indices]
        fixed_q01, fixed_q10, fixed_q90 = _bootstrap_quantiles(fixed_sample)
        freq_q01, freq_q10, freq_q90 = _bootstrap_quantiles(freq_sample)
        draws["delta_q01_db"].append(freq_q01 - fixed_q01)
        draws["delta_q10_db"].append(freq_q10 - fixed_q10)
        draws["delta_b80_db"].append((freq_q90 - freq_q10) - (fixed_q90 - fixed_q10))
    return {
        name: {
            "estimate": float(
                (_summary(freq)["central80_width_db"] - _summary(fixed)["central80_width_db"])
                if name == "delta_b80_db"
                else _quantile(freq, 0.01 if name == "delta_q01_db" else 0.10)
                - _quantile(fixed, 0.01 if name == "delta_q01_db" else 0.10)
            ),
            "percentile95_lo": float(np.quantile(values, 0.025)),
            "percentile95_hi": float(np.quantile(values, 0.975)),
        }
        for name, values in draws.items()
    }


def _recursive_differences(left: Any, right: Any, prefix: str = "") -> list[str]:
    if isinstance(left, dict) and isinstance(right, dict):
        out = []
        for key in sorted(set(left) | set(right)):
            path = f"{prefix}.{key}" if prefix else str(key)
            if key not in left or key not in right:
                out.append(path)
            else:
                out.extend(_recursive_differences(left[key], right[key], path))
        return out
    if left != right:
        return [prefix]
    return []


def _config_diff_receipt(input_root: Path | None = None) -> dict[str, object]:
    config_root = (input_root / "configs") if input_root is not None else (ROOT / "configs")
    loaded = {}
    for tier, name in BLER_CONFIG_NAMES.items():
        path = config_root / name
        with open(path, "r", encoding="utf-8") as handle:
            loaded[tier] = yaml.safe_load(handle)
    allowed = {
        "output_dir",
        "random_stream_namespace",
        "simulation.snr_points_db",
        "simulation.min_trials_per_snr",
        "simulation.max_trials_per_snr",
    }
    receipt: dict[str, object] = {"allowed_differences": sorted(allowed), "comparisons": {}}
    for tier in ("mid", "high"):
        differences = _recursive_differences(loaded["low"], loaded[tier])
        unexpected = sorted(set(differences) - allowed)
        receipt["comparisons"][f"low_vs_{tier}"] = {
            "differences": differences,
            "unexpected": unexpected,
        }
        if unexpected:
            raise RuntimeError(f"Unexpected BLER config differences low vs {tier}: {unexpected}")
    return receipt


def analyze(root: Path, bootstrap_repeats: int) -> None:
    analysis = root / "analysis"
    analysis.mkdir(parents=True, exist_ok=True)
    figure_dir = ROOT / "docs" / "figures" / "result-034"
    figure_dir.mkdir(parents=True, exist_ok=True)
    rows = _read_csv(root / "rsrp" / "paired_rsrp_trials.csv")
    trials = np.asarray([int(row["absolute_trial"]) for row in rows], dtype=np.int64)
    if not np.array_equal(trials, np.arange(1, len(rows) + 1, dtype=np.int64)):
        raise RuntimeError("RSRP absolute trials are not contiguous from 1.")
    keys = [row["trial_key"] for row in rows]
    if len(set(keys)) != len(keys):
        raise RuntimeError("RSRP trial keys are not unique.")
    fixed = np.asarray([float(row["fixed_block_power_db"]) for row in rows])
    freq = np.asarray([float(row["freq_block_power_db"]) for row in rows])
    fixed_single = np.asarray([float(row["fixed_single_re_power_db"]) for row in rows])
    freq_single = np.asarray([float(row["freq_single_re_power_db"]) for row in rows])
    if any(np.any(~np.isfinite(values)) for values in (fixed, freq, fixed_single, freq_single)):
        raise RuntimeError("RSRP input contains non-finite values.")

    summaries = {"FIXED_DFT0": _summary(fixed), "FREQ_SIDON_0137": _summary(freq)}
    fixed_linear = np.asarray([float(row["fixed_block_power_linear"]) for row in rows])
    freq_linear = np.asarray([float(row["freq_block_power_linear"]) for row in rows])
    summaries["FIXED_DFT0"].update(
        {
            "mean_linear": float(np.mean(fixed_linear)),
            "std_linear": float(np.std(fixed_linear, ddof=1)),
        }
    )
    summaries["FREQ_SIDON_0137"].update(
        {
            "mean_linear": float(np.mean(freq_linear)),
            "std_linear": float(np.std(freq_linear, ddof=1)),
        }
    )
    bootstrap = paired_bootstrap(
        fixed, freq, repeats=int(bootstrap_repeats), seed=20260916
    )
    single_quantile_differences = {
        f"q{int(probability * 100):02d}_difference_db": _quantile(freq_single, probability)
        - _quantile(fixed_single, probability)
        for probability in (0.01, 0.10, 0.50)
    }
    single_ks = float(ks_2samp(fixed_single, freq_single, method="auto").statistic)
    single_control_pass = bool(
        all(abs(value) <= 0.10 for value in single_quantile_differences.values())
        and single_ks <= 0.01
    )
    intervals = bootstrap
    concentration = bool(
        intervals["delta_q01_db"]["estimate"] > 0.0
        and intervals["delta_q10_db"]["estimate"] > 0.0
        and intervals["delta_b80_db"]["estimate"] < 0.0
        and intervals["delta_q01_db"]["percentile95_lo"] > 0.0
        and intervals["delta_q10_db"]["percentile95_lo"] > 0.0
        and intervals["delta_b80_db"]["percentile95_hi"] < 0.0
    )
    dkw_epsilon = float(np.sqrt(np.log(2.0 / 0.05) / (2.0 * len(rows))))

    rsrp_summary = {
        "schema": "plan034-rsrp-analysis-v1",
        "trials": len(rows),
        "quantile_definition": "left-continuous empirical inverse CDF (numpy inverted_cdf)",
        "candidates": summaries,
        "paired_bootstrap": {
            "seed": 20260916,
            "repeats": int(bootstrap_repeats),
            "differences": bootstrap,
        },
        "dkw95_epsilon": dkw_epsilon,
        "single_re_control": {
            **single_quantile_differences,
            "maximum_absolute_empirical_cdf_difference": single_ks,
            "pass": single_control_pass,
        },
        "concentration_criterion_pass": concentration,
    }
    (analysis / "rsrp_summary.json").write_text(
        json.dumps(rsrp_summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    cdf_rows: list[dict[str, object]] = []
    for candidate_id, values in (("FIXED_DFT0", fixed), ("FREQ_SIDON_0137", freq)):
        ordered = np.sort(values)
        probabilities = np.arange(1, len(ordered) + 1, dtype=np.float64) / len(ordered)
        for value, probability in zip(ordered, probabilities):
            cdf_rows.append(
                {
                    "candidate_id": candidate_id,
                    "power_db": float(value),
                    "empirical_cdf": float(probability),
                    "dkw95_lo": float(max(0.0, probability - dkw_epsilon)),
                    "dkw95_hi": float(min(1.0, probability + dkw_epsilon)),
                }
            )
    _write_csv(analysis / "rsrp_cdf_input.csv", cdf_rows)

    bler_rows: list[dict[str, object]] = []
    bler_array_checks: list[dict[str, object]] = []
    for tier in ("low", "mid", "high"):
        for row in _read_csv(root / "bler" / tier / "bler_points.csv"):
            converted: dict[str, object] = {"tier": tier}
            for key, value in row.items():
                converted[key] = value
            bler_rows.append(converted)
            candidate_id = "".join(
                character if character.isalnum() or character in "-_" else "_"
                for character in str(row["candidate_id"])
            )
            snr_tag = str(float(row["snr_db"])).replace("-", "m").replace(".", "p")
            stem = f"{candidate_id}_snr_{snr_tag}"
            arrays = root / "bler" / tier / "trial_error_flags"
            flags = np.load(arrays / f"{stem}_error_flags.npy")
            ce_nmse = np.load(arrays / f"{stem}_ce_nmse.npy")
            mrc = np.load(arrays / f"{stem}_mrc_denominator_trial_mean.npy")
            trials_expected = int(row["trials"])
            errors_expected = int(row["errors"])
            lengths_ok = len(flags) == len(ce_nmse) == len(mrc) == trials_expected
            errors_ok = int(np.sum(flags)) == errors_expected
            ce_ok = np.isclose(
                10.0 * np.log10(max(float(np.mean(ce_nmse)), 1e-30)),
                float(row["ce_nmse_db"]),
                rtol=0.0,
                atol=1e-12,
            )
            mrc_ok = np.isclose(
                float(np.mean(mrc)),
                float(row["mrc_denominator_mean"]),
                rtol=0.0,
                atol=1e-12,
            )
            finite_positive_ok = bool(
                np.all(np.isfinite(ce_nmse))
                and np.all(ce_nmse >= 0.0)
                and np.all(np.isfinite(mrc))
                and np.all(mrc > 0.0)
            )
            check = {
                "snr_db": float(row["snr_db"]),
                "lengths_ok": bool(lengths_ok),
                "error_recount_ok": bool(errors_ok),
                "ce_linear_mean_to_db_ok": bool(ce_ok),
                "mrc_mean_ok": bool(mrc_ok),
                "finite_positive_ok": finite_positive_ok,
            }
            bler_array_checks.append(check)
            if not all(value for key, value in check.items() if key != "snr_db"):
                raise RuntimeError(f"BLER array audit failed: {check}")
    bler_rows.sort(key=lambda row: float(row["snr_db"]))
    if [float(row["snr_db"]) for row in bler_rows] != [float(value) for value in range(-3, 6)]:
        raise RuntimeError("BLER inputs do not contain exactly the frozen -3..5 dB grid.")
    _write_csv(analysis / "bler_combined.csv", bler_rows)

    config_receipt = _config_diff_receipt(root)
    (analysis / "config_diff_receipt.json").write_text(
        json.dumps(config_receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    style = {
        "figure_width_cm": 13.0,
        "fixed_color": "#1f77b4",
        "freq_color": "#d62728",
        "fixed_linestyle": "-",
        "freq_linestyle": "--",
        "fixed_marker": "o",
        "freq_marker": "s",
        "line_width": 2.5,
        "marker_size": 8.0,
        "dkw_band_alpha": 0.12,
        "axis_label_fontsize": 16,
        "legend_fontsize": 16,
        "tick_fontsize": 14,
        "dpi": 200,
    }
    (analysis / "plot_style.json").write_text(
        json.dumps(style, indent=2) + "\n", encoding="utf-8"
    )

    preview_width = style["figure_width_cm"] / 2.54
    fig, axis = plt.subplots(figsize=(8.6, 6.2))
    for candidate_id, label, values, color, linestyle, marker in (
        (
            "FIXED_DFT0",
            "Fixed DFT0",
            fixed,
            style["fixed_color"],
            style["fixed_linestyle"],
            style["fixed_marker"],
        ),
        (
            "FREQ_SIDON_0137",
            "Freq Sidon [0,1,3,7]",
            freq,
            style["freq_color"],
            style["freq_linestyle"],
            style["freq_marker"],
        ),
    ):
        ordered = np.sort(values)
        probability = np.arange(1, len(ordered) + 1) / len(ordered)
        axis.fill_between(
            ordered,
            np.maximum(0.0, probability - dkw_epsilon),
            np.minimum(1.0, probability + dkw_epsilon),
            color=color,
            alpha=style["dkw_band_alpha"],
            linewidth=0.0,
        )
        axis.plot(
            ordered,
            probability,
            label=label,
            color=color,
            linestyle=linestyle,
            linewidth=style["line_width"],
            marker=marker,
            markersize=style["marker_size"],
            markevery=5000,
        )
    axis.set_xlabel(
        "Normalized block-average receive power (dB)",
        fontsize=style["axis_label_fontsize"],
    )
    axis.set_ylabel("Empirical CDF", fontsize=style["axis_label_fontsize"])
    axis.grid(True, alpha=0.3)
    axis.tick_params(labelsize=style["tick_fontsize"])
    axis.legend(fontsize=style["legend_fontsize"])
    fig.tight_layout()
    fig.savefig(analysis / "rsrp_cdf.png", dpi=style["dpi"])
    fig.savefig(figure_dir / "rsrp_cdf.png", dpi=style["dpi"])
    fig.set_size_inches(preview_width, preview_width * 0.72)
    fig.tight_layout()
    fig.savefig(analysis / "rsrp_cdf_preview.png", dpi=style["dpi"])
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(8.6, 6.2))
    snr = np.asarray([float(row["snr_db"]) for row in bler_rows])
    bler = np.asarray([float(row["bler"]) for row in bler_rows])
    trials_per_point = np.asarray([int(row["trials"]) for row in bler_rows])
    display = np.where(bler > 0.0, bler, 0.5 / trials_per_point)
    axis.semilogy(
        snr,
        display,
        color=style["fixed_color"],
        linestyle=style["fixed_linestyle"],
        marker=style["fixed_marker"],
        linewidth=style["line_width"],
        markersize=style["marker_size"],
        label="Fixed DFT0",
    )
    axis.set_xlabel("SNR (dB)", fontsize=style["axis_label_fontsize"])
    axis.set_ylabel("PDCCH BLER", fontsize=style["axis_label_fontsize"])
    axis.grid(True, which="both", alpha=0.3)
    axis.tick_params(labelsize=style["tick_fontsize"])
    axis.legend(fontsize=style["legend_fontsize"])
    fig.tight_layout()
    fig.savefig(analysis / "fixed_dft0_bler.png", dpi=style["dpi"])
    fig.savefig(figure_dir / "fixed_dft0_bler.png", dpi=style["dpi"])
    fig.set_size_inches(preview_width, preview_width * 0.72)
    fig.tight_layout()
    fig.savefig(analysis / "fixed_dft0_bler_preview.png", dpi=style["dpi"])
    plt.close(fig)

    audit = {
        "schema": "plan034-analysis-audit-v1",
        "rsrp_trial_continuity": True,
        "rsrp_trial_key_uniqueness": True,
        "rsrp_candidate_pairing": True,
        "single_re_control_pass": single_control_pass,
        "bler_snr_grid": [float(row["snr_db"]) for row in bler_rows],
        "bler_array_checks": bler_array_checks,
        "config_diff_receipt": "analysis/config_diff_receipt.json",
        "rsrp_summary": "analysis/rsrp_summary.json",
        "rsrp_cdf_input": "analysis/rsrp_cdf_input.csv",
        "bler_input": "analysis/bler_combined.csv",
    }
    (analysis / "audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def _read_per_rx_trials(root: Path) -> tuple[np.ndarray, list[str], np.ndarray, np.ndarray]:
    csv_path = root / "rsrp" / "per_rx_rsrp_trials.csv"
    receipt_path = root / "rsrp" / "run_receipt.json"
    if not receipt_path.exists():
        raise FileNotFoundError(receipt_path)
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    trial_count = int(receipt["trials"])
    linear = np.empty((trial_count, 4, len(PER_RX_CANDIDATES)), dtype=np.float64)
    db = np.empty_like(linear)
    trial_keys = [""] * trial_count
    rows_per_trial = 4 * len(PER_RX_CANDIDATES)
    row_count = 0
    with open(csv_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != PER_RX_FIELDNAMES:
            raise RuntimeError("Per-Rx CSV field list does not match the frozen schema.")
        for row_index, row in enumerate(reader):
            trial_index = row_index // rows_per_trial
            position = row_index % rows_per_trial
            candidate_index = position // 4
            rx_index = position % 4
            if trial_index >= trial_count:
                raise RuntimeError("Per-Rx CSV contains more rows than its run receipt.")
            expected_trial = trial_index + 1
            candidate_id = PER_RX_CANDIDATES[candidate_index]
            metadata = PER_RX_METADATA[candidate_id]
            if int(row["absolute_trial"]) != expected_trial:
                raise RuntimeError(f"Expected absolute trial {expected_trial} at data row {row_index + 1}.")
            if row["candidate_id"] != candidate_id or int(row["rx_index"]) != rx_index:
                raise RuntimeError(f"Candidate/Rx ordering mismatch at data row {row_index + 1}.")
            if row["trial_key"] != row["channel_key"]:
                raise RuntimeError(f"Trial/channel key mismatch at data row {row_index + 1}.")
            if json.loads(row["delay_grid_coordinates"]) != metadata["delay_grid_coordinates"]:
                raise RuntimeError(f"Delay-grid metadata mismatch for {candidate_id}.")
            if json.loads(row["delay_ns"]) != metadata["delay_ns"]:
                raise RuntimeError(f"Delay metadata mismatch for {candidate_id}.")
            expected_denominator = metadata["phase_denominator"]
            actual_denominator = None if row["phase_denominator"] == "" else int(row["phase_denominator"])
            if actual_denominator != expected_denominator:
                raise RuntimeError(f"Phase denominator mismatch for {candidate_id}.")
            if (row["allow_duplicate_delays"].lower() == "true") != bool(
                metadata["allow_duplicate_delays"]
            ):
                raise RuntimeError(f"Duplicate-delay flag mismatch for {candidate_id}.")
            if float(row["total_transmit_power"]) != 1.0:
                raise RuntimeError("Per-Rx total transmit power must be exactly 1.0.")
            value = float(row["block_power_linear"])
            value_db = float(row["block_power_db"])
            if not np.isfinite(value) or value <= 0.0 or not np.isfinite(value_db):
                raise RuntimeError(f"Invalid per-Rx power at data row {row_index + 1}.")
            if not np.isclose(value_db, 10.0 * np.log10(value), rtol=0.0, atol=1e-12):
                raise RuntimeError(f"Linear/dB mismatch at data row {row_index + 1}.")
            linear[trial_index, rx_index, candidate_index] = value
            db[trial_index, rx_index, candidate_index] = value_db
            if not trial_keys[trial_index]:
                trial_keys[trial_index] = row["trial_key"]
            elif trial_keys[trial_index] != row["trial_key"]:
                raise RuntimeError(f"Candidate channel pairing failed for trial {expected_trial}.")
            row_count += 1
    if row_count != trial_count * rows_per_trial:
        raise RuntimeError(
            f"Per-Rx row count mismatch: found {row_count}, expected {trial_count * rows_per_trial}."
        )
    trials = np.arange(1, trial_count + 1, dtype=np.int64)
    return trials, trial_keys, linear, db


def _paired_bootstrap_per_rx(
    values: np.ndarray, *, repeats: int, seed: int
) -> dict[str, dict[str, dict[str, float]]]:
    if values.ndim != 2 or values.shape[1] != len(PER_RX_CANDIDATES):
        raise ValueError("Per-Rx bootstrap values must have shape [trial, candidate].")
    rng = np.random.default_rng(int(seed))
    draws = {
        candidate_id: {
            "delta_q01_db": [],
            "delta_q10_db": [],
            "delta_b80_db": [],
        }
        for candidate_id in PER_RX_CANDIDATES[1:]
    }
    for _ in range(int(repeats)):
        indices = rng.integers(0, values.shape[0], size=values.shape[0])
        fixed_q01, fixed_q10, fixed_q90 = _bootstrap_quantiles(values[indices, 0])
        for candidate_index, candidate_id in enumerate(PER_RX_CANDIDATES[1:], start=1):
            q01, q10, q90 = _bootstrap_quantiles(values[indices, candidate_index])
            draws[candidate_id]["delta_q01_db"].append(q01 - fixed_q01)
            draws[candidate_id]["delta_q10_db"].append(q10 - fixed_q10)
            draws[candidate_id]["delta_b80_db"].append(
                (q90 - q10) - (fixed_q90 - fixed_q10)
            )
    fixed_summary = _summary(values[:, 0])
    result: dict[str, dict[str, dict[str, float]]] = {}
    for candidate_index, candidate_id in enumerate(PER_RX_CANDIDATES[1:], start=1):
        candidate_summary = _summary(values[:, candidate_index])
        estimates = {
            "delta_q01_db": candidate_summary["q01_db"] - fixed_summary["q01_db"],
            "delta_q10_db": candidate_summary["q10_db"] - fixed_summary["q10_db"],
            "delta_b80_db": candidate_summary["central80_width_db"]
            - fixed_summary["central80_width_db"],
        }
        result[candidate_id] = {
            metric: {
                "estimate": float(estimates[metric]),
                "percentile95_lo": float(np.quantile(samples, 0.025)),
                "percentile95_hi": float(np.quantile(samples, 0.975)),
            }
            for metric, samples in draws[candidate_id].items()
        }
    return result


def _original_reconstruction_audit(
    root: Path, trial_keys: list[str], linear: np.ndarray
) -> dict[str, object]:
    original_path = root.parent / "rsrp" / "paired_rsrp_trials.csv"
    if not original_path.exists():
        return {"status": "missing_original_input", "path": str(original_path)}
    original = _read_csv(original_path)
    if len(original) != linear.shape[0]:
        raise RuntimeError(
            f"Original RSRP trial count {len(original)} differs from per-Rx {linear.shape[0]}."
        )
    checks: dict[str, object] = {
        "status": "checked",
        "path": str(original_path),
        "trial_keys_equal": [row["trial_key"] for row in original] == trial_keys,
        "candidates": {},
    }
    if not checks["trial_keys_equal"]:
        raise RuntimeError("Original and per-Rx trial keys differ.")
    for candidate_id, candidate_index, field in (
        ("FIXED_DFT0", 0, "fixed_block_power_linear"),
        ("FREQ_SIDON_0137", 1, "freq_block_power_linear"),
    ):
        expected = np.asarray([float(row[field]) for row in original], dtype=np.float64)
        reconstructed = np.mean(linear[:, :, candidate_index], axis=1)
        maximum_absolute_difference = float(np.max(np.abs(reconstructed - expected)))
        exact = bool(np.array_equal(reconstructed, expected))
        within_tolerance = bool(
            np.allclose(reconstructed, expected, rtol=0.0, atol=1e-12)
        )
        checks["candidates"][candidate_id] = {
            "array_equal": exact,
            "linear_atol_1e-12": within_tolerance,
            "maximum_absolute_difference": maximum_absolute_difference,
        }
        if not within_tolerance:
            raise RuntimeError(
                f"{candidate_id} per-Rx reconstruction differs from original by "
                f"{maximum_absolute_difference:.3e}."
            )
    return checks


def analyze_per_rx(root: Path, bootstrap_repeats: int) -> None:
    analysis = root / "analysis"
    analysis.mkdir(parents=True, exist_ok=True)
    trials, trial_keys, linear, db = _read_per_rx_trials(root)
    if len(set(trial_keys)) != len(trial_keys):
        raise RuntimeError("Per-Rx trial keys are not unique.")
    summaries: list[dict[str, object]] = []
    summary_lookup: dict[tuple[int, str], dict[str, float]] = {}
    for rx_index in range(4):
        for candidate_index, candidate_id in enumerate(PER_RX_CANDIDATES):
            summary = _summary(db[:, rx_index, candidate_index])
            summary.update(
                {
                    "mean_linear": float(np.mean(linear[:, rx_index, candidate_index])),
                    "std_linear": float(
                        np.std(linear[:, rx_index, candidate_index], ddof=1)
                    ),
                }
            )
            summary_lookup[(rx_index, candidate_id)] = summary
            summaries.append({"rx_index": rx_index, "candidate_id": candidate_id, **summary})
    _write_csv(analysis / "per_rx_summary.csv", summaries)

    comparisons: list[dict[str, object]] = []
    bootstrap_json: dict[str, object] = {}
    for rx_index in range(4):
        result = _paired_bootstrap_per_rx(
            db[:, rx_index, :], repeats=bootstrap_repeats, seed=20260916
        )
        bootstrap_json[f"rx{rx_index}"] = result
        for candidate_id, metrics in result.items():
            for metric, interval in metrics.items():
                comparisons.append(
                    {
                        "rx_index": rx_index,
                        "candidate_id": candidate_id,
                        "baseline_id": "FIXED_DFT0",
                        "metric": metric,
                        **interval,
                    }
                )
    _write_csv(analysis / "per_rx_paired_bootstrap.csv", comparisons)

    extrema_rows: list[dict[str, object]] = []
    metric_names = (
        "q01_db",
        "q10_db",
        "q50_db",
        "q90_db",
        "central80_width_db",
        "mean_db",
        "std_db",
        "mean_linear",
        "std_linear",
    )
    for candidate_id in PER_RX_CANDIDATES:
        for metric in metric_names:
            values = np.asarray(
                [summary_lookup[(rx, candidate_id)][metric] for rx in range(4)]
            )
            extrema_rows.append(
                {
                    "candidate_id": candidate_id,
                    "metric": metric,
                    "minimum": float(np.min(values)),
                    "maximum": float(np.max(values)),
                    "range": float(np.max(values) - np.min(values)),
                }
            )
    _write_csv(analysis / "per_rx_extrema.csv", extrema_rows)

    dkw_epsilon = float(
        np.sqrt(np.log(2.0 * 16.0 / 0.05) / (2.0 * len(trials)))
    )
    consistency_rows: list[dict[str, object]] = []
    all_consistent = True
    for candidate_index, candidate_id in enumerate(PER_RX_CANDIDATES):
        for left_rx in range(4):
            for right_rx in range(left_rx + 1, 4):
                statistic = float(
                    ks_2samp(
                        db[:, left_rx, candidate_index],
                        db[:, right_rx, candidate_index],
                        method="auto",
                    ).statistic
                )
                consistent = statistic <= 2.0 * dkw_epsilon
                all_consistent = all_consistent and consistent
                consistency_rows.append(
                    {
                        "candidate_id": candidate_id,
                        "left_rx": left_rx,
                        "right_rx": right_rx,
                        "maximum_absolute_empirical_cdf_difference": statistic,
                        "bonferroni_dkw_pair_threshold": 2.0 * dkw_epsilon,
                        "consistent": consistent,
                    }
                )
    _write_csv(analysis / "per_rx_consistency.csv", consistency_rows)
    reconstruction = _original_reconstruction_audit(root, trial_keys, linear)

    cdf_path = analysis / "per_rx_cdf_input.csv"
    with open(cdf_path, "w", encoding="utf-8", newline="") as handle:
        fields = (
            "rx_index",
            "candidate_id",
            "power_db",
            "empirical_cdf",
            "bonferroni_dkw_lo",
            "bonferroni_dkw_hi",
        )
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for rx_index in range(4):
            for candidate_index, candidate_id in enumerate(PER_RX_CANDIDATES):
                ordered = np.sort(db[:, rx_index, candidate_index])
                probabilities = np.arange(1, len(ordered) + 1) / len(ordered)
                for value, probability in zip(ordered, probabilities):
                    writer.writerow(
                        {
                            "rx_index": rx_index,
                            "candidate_id": candidate_id,
                            "power_db": float(value),
                            "empirical_cdf": float(probability),
                            "bonferroni_dkw_lo": float(
                                max(0.0, probability - dkw_epsilon)
                            ),
                            "bonferroni_dkw_hi": float(
                                min(1.0, probability + dkw_epsilon)
                            ),
                        }
                    )

    style = {
        "FIXED_DFT0": {"color": "#1f77b4", "linestyle": "-", "label": "Fixed DFT0; no delay"},
        "FREQ_SIDON_0137": {
            "color": "#d62728",
            "linestyle": "--",
            "label": "Sidon [0,925.926,2777.778,6481.481] ns",
        },
        "CDD911": {
            "color": "#2ca02c",
            "linestyle": "-.",
            "label": "CDD911 [0,0,911,911] ns",
        },
        "CDD130": {
            "color": "#9467bd",
            "linestyle": ":",
            "label": "CDD130 [0,0,130,130] ns",
        },
    }
    (analysis / "per_rx_plot_style.json").write_text(
        json.dumps(style, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    fig, axes = plt.subplots(2, 2, figsize=(15.0, 11.0), sharex=True, sharey=True)
    for rx_index, axis in enumerate(axes.flat):
        for candidate_index, candidate_id in enumerate(PER_RX_CANDIDATES):
            ordered = np.sort(db[:, rx_index, candidate_index])
            probability = np.arange(1, len(ordered) + 1) / len(ordered)
            candidate_style = style[candidate_id]
            axis.fill_between(
                ordered,
                np.maximum(0.0, probability - dkw_epsilon),
                np.minimum(1.0, probability + dkw_epsilon),
                color=candidate_style["color"],
                alpha=0.08,
                linewidth=0.0,
            )
            axis.plot(
                ordered,
                probability,
                color=candidate_style["color"],
                linestyle=candidate_style["linestyle"],
                linewidth=2.5,
                label=candidate_style["label"],
            )
        axis.set_title(f"Rx{rx_index}", fontsize=16)
        axis.grid(True, alpha=0.3)
        axis.tick_params(labelsize=14)
    axes[1, 0].set_xlabel("Normalized block-average receive power (dB)", fontsize=16)
    axes[1, 1].set_xlabel("Normalized block-average receive power (dB)", fontsize=16)
    axes[0, 0].set_ylabel("Empirical CDF", fontsize=16)
    axes[1, 0].set_ylabel("Empirical CDF", fontsize=16)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, fontsize=14)
    fig.tight_layout(rect=(0.0, 0.07, 1.0, 1.0))
    figure_path = analysis / "per_rx_4waveform_cdf.png"
    fig.savefig(figure_path, dpi=200)
    fig.set_size_inches(13.0 / 2.54, 10.0 / 2.54)
    fig.tight_layout(rect=(0.0, 0.16, 1.0, 1.0))
    preview_path = analysis / "per_rx_4waveform_cdf_preview.png"
    fig.savefig(preview_path, dpi=200)
    plt.close(fig)

    summary_payload = {
        "schema": "plan034-per-rx-rsrp-analysis-v1",
        "trials": int(len(trials)),
        "candidate_ids": list(PER_RX_CANDIDATES),
        "rx_indices": [0, 1, 2, 3],
        "bootstrap_seed": 20260916,
        "bootstrap_repeats": int(bootstrap_repeats),
        "bonferroni_dkw_fwer": 0.95,
        "bonferroni_dkw_half_width": dkw_epsilon,
        "all_rx_pairs_consistent": bool(all_consistent),
        "original_average_reconstruction": reconstruction,
        "paired_bootstrap": bootstrap_json,
    }
    summary_path = analysis / "per_rx_analysis_summary.json"
    summary_path.write_text(
        json.dumps(summary_payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    metadata = {
        "schema": "plan034-per-rx-analysis-metadata-v1",
        "input": {
            "path": "rsrp/per_rx_rsrp_trials.csv",
            "sha256": _sha256(root / "rsrp" / "per_rx_rsrp_trials.csv"),
        },
        "outputs": {
            path.name: _sha256(path)
            for path in (
                analysis / "per_rx_summary.csv",
                analysis / "per_rx_paired_bootstrap.csv",
                analysis / "per_rx_extrema.csv",
                analysis / "per_rx_consistency.csv",
                cdf_path,
                summary_path,
                figure_path,
                preview_path,
            )
        },
    }
    (analysis / "per_rx_analysis_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze Plan-034 RSRP and PDCCH BLER outputs.")
    parser.add_argument("--input-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--bootstrap-repeats", type=int, default=1000)
    parser.add_argument("--mode", choices=("original", "per-rx"), default="original")
    args = parser.parse_args()
    if args.bootstrap_repeats <= 0:
        raise ValueError("--bootstrap-repeats must be positive.")
    if args.mode == "per-rx":
        analyze_per_rx(args.input_root.resolve(), args.bootstrap_repeats)
    else:
        analyze(args.input_root.resolve(), args.bootstrap_repeats)


if __name__ == "__main__":
    main()
