"""Run and plot the plan-033 8Tx/4Rx ideal-CSI Sidon/PRG6 diagnostic."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.mcs import build_tb_layout, get_mcs
from cdd_lls.phy.channel_tdl import generate_sionna_tdl_channel_active
from cdd_lls.phy.ldpc import SionnaLDPCAdapter
from cdd_lls.phy.precoding import build_prg_dft_precoder_batch, equivalent_channel
from cdd_lls.phy.qam import qam_modulate
from cdd_lls.phy.resource_grid import local_indices_for_subcarriers
from tools import run_plan028_csi_curves as curves028
from tools.run_bler_curves import _decode_flags
from tools.run_plan025_delay_matched_tdl import stable_seed, wilson
from tools.run_plan032_tdl_mobility import _repo_relative, _safe, _snr_key
from tools.run_plan033_tdl_mobility_mimo import (
    N_RX,
    PRG_SIZE_RB,
    _build_scene,
    awgn_variance,
    candidate_ids,
)


SCHEMA = "plan033-ideal-sidon-prg6-diagnostic-v1"
SCENARIO_ID = "A100_NT8_NR4_V3"
N_TX = 8
SPEED_KMH = 3.0
SEED = 20260727
BATCH_SIZE = 25
TRIALS = 1000
SNR_DB = (3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0)
CURRENT_SYMBOL_INDICES = tuple(range(140, 150))
DEFAULT_OUTPUT = (
    ROOT
    / "outputs"
    / "experiment033_tdl_mobility_mimo"
    / "20260916_ideal_sidon_prg6"
    / SCENARIO_ID
)

LABELS = {
    "A100_NT8_S0_SIDON": "Sidon CDD",
    "A100_NT8_TRANSPARENT_PRG6": "Transparent PRG6",
}
STYLES = {
    "A100_NT8_S0_SIDON": {"color": "#009E73", "linestyle": "-", "marker": "^"},
    "A100_NT8_TRANSPARENT_PRG6": {"color": "#000000", "linestyle": "-", "marker": "o"},
}


def _write_rows(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _read_rows(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _resolved_config(output: Path) -> dict:
    return {
        "schema": SCHEMA,
        "scenario_id": SCENARIO_ID,
        "n_tx": N_TX,
        "n_rx": N_RX,
        "speed_kmh": SPEED_KMH,
        "tdl_profile": "A",
        "delay_spread_ns": 100.0,
        "receiver": "ideal CSI per data RE with 4Rx MRC",
        "precoder_normalize": True,
        "snr_definition": "per-Rx average received SNR with unit total transmit power",
        "noise_variance": "1 / SNR_linear",
        "seed": SEED,
        "batch_size": BATCH_SIZE,
        "trials_per_snr": TRIALS,
        "snr_db": list(SNR_DB),
        "current_symbol_indices": list(CURRENT_SYMBOL_INDICES),
        "paired_fields": ["physical channel", "payload", "data noise", "absolute trial"],
        "candidates": [
            {
                "candidate_id": "A100_NT8_S0_SIDON",
                "delay_grid_coordinates": [0, 1, 3, 7, 12, 20, 30, 65],
            },
            {
                "candidate_id": "A100_NT8_TRANSPARENT_PRG6",
                "prg_size_rb": 6,
                "prg_vector_indices": list(range(8)),
            },
        ],
        "output_dir": _repo_relative(output),
    }


def _scene_config(output: Path) -> dict:
    return {
        "n_tx": N_TX,
        "n_rx": N_RX,
        "speed_kmh": SPEED_KMH,
        "precoder_normalize": True,
        "output_dir": str(output),
    }


def _candidate_channels(realization_h: np.ndarray, grid, precoders: dict) -> dict[str, np.ndarray]:
    ids = candidate_ids(N_TX)
    order = np.tile(np.arange(8, dtype=np.int64), (len(realization_h), 1))
    prg = build_prg_dft_precoder_batch(
        grid,
        n_tx=N_TX,
        prg_size_rb=PRG_SIZE_RB,
        prg_vector_indices=order,
        normalize=True,
    )
    return {
        ids["sidon"]: equivalent_channel(realization_h, precoders["S0_SIDON"].C),
        ids["prg"]: equivalent_channel(realization_h, prg),
    }


def _merge(output: Path, candidate_order: list[str]) -> list[dict]:
    intervals = _read_rows(output / "intervals.csv")
    merged = []
    for candidate_id in candidate_order:
        for snr_db in SNR_DB:
            selected = [
                row
                for row in intervals
                if row["candidate_id"] == candidate_id
                and math.isclose(float(row["snr_db"]), snr_db, abs_tol=1e-12)
            ]
            selected.sort(key=lambda row: int(row["trial_start"]))
            if not selected:
                continue
            expected_start = 1
            for row in selected:
                if int(row["trial_start"]) != expected_start:
                    raise RuntimeError(f"Non-contiguous intervals for {candidate_id} at {snr_db:g} dB")
                expected_start = int(row["trial_end"]) + 1
            trials = sum(int(row["trials"]) for row in selected)
            errors = sum(int(row["tb_errors"]) for row in selected)
            low, high = wilson(errors, trials)
            merged.append(
                {
                    "scenario_id": SCENARIO_ID,
                    "n_tx": N_TX,
                    "n_rx": N_RX,
                    "speed_kmh": SPEED_KMH,
                    "candidate_id": candidate_id,
                    "receiver": "ideal",
                    "snr_db": snr_db,
                    "trials": trials,
                    "tb_errors": errors,
                    "bler": errors / trials,
                    "bler_wilson95_lo": low,
                    "bler_wilson95_hi": high,
                }
            )
    if merged:
        _write_rows(output / "ideal_csi_bler_points.csv", merged)
    return merged


def run(output: Path) -> list[dict]:
    output.mkdir(parents=True, exist_ok=True)
    resolved = _resolved_config(output)
    (output / "resolved_config.json").write_text(
        json.dumps(resolved, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (output / "curve_styles.json").write_text(
        json.dumps(STYLES, indent=2) + "\n", encoding="utf-8"
    )

    channel, grid, precoders, _ = _build_scene(_scene_config(output))
    ids = candidate_ids(N_TX)
    candidate_order = [ids["sidon"], ids["prg"]]
    if set(candidate_order) != set(LABELS):
        raise RuntimeError("Plan-033 candidate IDs changed.")
    if grid.n_data_re != 5568 or list(np.unique(grid.pilot_symbol_indices)) != [2, 7]:
        raise RuntimeError("Plan-033 resource grid changed.")
    if not np.allclose(
        np.sum(np.abs(precoders["S0_SIDON"].C) ** 2, axis=1), 1.0, atol=1e-12
    ):
        raise RuntimeError("Sidon precoder normalization changed.")

    mcs = get_mcs("nr_256qam", 8, None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(tb.cb_k_values, tb.cb_e_values, num_iter=8, llr_clip=50.0)
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    intervals_path = output / "intervals.csv"
    stored = _read_rows(intervals_path)

    for snr_db in SNR_DB:
        counts = {
            candidate_id: sum(
                int(row["trials"])
                for row in stored
                if row["candidate_id"] == candidate_id
                and math.isclose(float(row["snr_db"]), snr_db, abs_tol=1e-12)
            )
            for candidate_id in candidate_order
        }
        if len(set(counts.values())) != 1:
            raise RuntimeError(f"Paired trial counts differ at {snr_db:g} dB: {counts}")
        current_total = next(iter(counts.values()))
        if current_total >= TRIALS:
            print(f"[resume] {SCENARIO_ID} snr={snr_db:g} trials={current_total}", flush=True)
            continue

        trial_start = current_total + 1
        trial_end = TRIALS
        trial_count = trial_end - trial_start + 1
        flags = {candidate_id: np.zeros(trial_count, dtype=bool) for candidate_id in candidate_order}
        noise_variance = awgn_variance(snr_db)
        started = time.time()

        for absolute_start in range(trial_start, trial_end + 1, BATCH_SIZE):
            batch = min(BATCH_SIZE, trial_end - absolute_start + 1)
            offset = absolute_start - trial_start
            payload_rng = np.random.default_rng(
                stable_seed(SEED, SCENARIO_ID, snr_db, absolute_start, "payload")
            )
            payload = [
                payload_rng.integers(0, 2, size=int(k), dtype=np.int8)
                for k in tb.cb_k_values
            ]
            symbols = qam_modulate(np.concatenate(adapter.encode(payload)), int(mcs.qm))
            realization = generate_sionna_tdl_channel_active(
                grid,
                channel,
                n_tx=N_TX,
                n_rx=N_RX,
                batch_size=batch,
                seed=stable_seed(SEED, SCENARIO_ID, snr_db, absolute_start, "channel"),
                time_sample_indices=CURRENT_SYMBOL_INDICES,
            )
            effective = _candidate_channels(realization.H, grid, precoders)
            noise_rng = np.random.default_rng(
                stable_seed(SEED, SCENARIO_ID, snr_db, absolute_start, "noise")
            )
            data_noise = math.sqrt(noise_variance / 2.0) * (
                noise_rng.normal(size=(batch, N_RX, grid.n_data_re))
                + 1j * noise_rng.normal(size=(batch, N_RX, grid.n_data_re))
            )

            decoded = []
            for candidate_id in candidate_order:
                true_data = effective[candidate_id][
                    :, :, grid.data_symbol_indices, data_local
                ]
                received = true_data * symbols[None, None, :] + data_noise
                equalized, effective_noise = curves028.ideal_csi_equalize(
                    true_data, received, noise_variance
                )
                decoded.append((candidate_id, equalized, effective_noise))
            combined = _decode_flags(
                adapter,
                np.concatenate([item[1] for item in decoded], axis=0),
                np.concatenate([item[2] for item in decoded], axis=0),
                int(mcs.qm),
                payload,
            )
            cursor = 0
            for candidate_id, _, _ in decoded:
                flags[candidate_id][offset : offset + batch] = combined[cursor : cursor + batch]
                cursor += batch
            completed = absolute_start + batch - 1
            print(
                f"[ideal {SCENARIO_ID}] snr={snr_db:g} "
                f"absolute_trials={completed}/{trial_end} elapsed={time.time()-started:.1f}s",
                flush=True,
            )

        interval_rows = []
        for candidate_id in candidate_order:
            flag_path = (
                output
                / "error_flags"
                / _snr_key(snr_db)
                / f"{_safe(candidate_id)}_t{trial_start:06d}_{trial_end:06d}.npy"
            )
            flag_path.parent.mkdir(parents=True, exist_ok=True)
            np.save(flag_path, flags[candidate_id])
            errors = int(np.sum(flags[candidate_id]))
            low, high = wilson(errors, trial_count)
            interval_rows.append(
                {
                    "scenario_id": SCENARIO_ID,
                    "candidate_id": candidate_id,
                    "receiver": "ideal",
                    "snr_db": snr_db,
                    "trial_start": trial_start,
                    "trial_end": trial_end,
                    "trials": trial_count,
                    "tb_errors": errors,
                    "bler": errors / trial_count,
                    "bler_wilson95_lo": low,
                    "bler_wilson95_hi": high,
                    "error_flags": _repo_relative(flag_path),
                    "seed": SEED,
                }
            )
        stored.extend(interval_rows)
        stored.sort(
            key=lambda row: (
                str(row["candidate_id"]),
                float(row["snr_db"]),
                int(row["trial_start"]),
            )
        )
        _write_rows(intervals_path, stored)

    return _merge(output, candidate_order)


def plot(output: Path, rows: list[dict]) -> None:
    if not rows:
        rows = _read_rows(output / "ideal_csi_bler_points.csv")
    expected = len(LABELS) * len(SNR_DB)
    if len(rows) != expected:
        raise RuntimeError(f"Expected {expected} merged points, found {len(rows)}")
    if any(int(row["trials"]) != TRIALS for row in rows):
        raise RuntimeError("Plot requires the frozen 1000 trials at every point.")

    plot_rows = []
    for row in rows:
        copied = dict(row)
        copied["plot_bler"] = max(float(row["bler"]), 0.5 / int(row["trials"]))
        plot_rows.append(copied)
    _write_rows(output / "plot_data.csv", plot_rows)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def draw(path: Path, preview: bool) -> None:
        size = (5.12, 4.1) if preview else (9.5, 6.6)
        figure, axis = plt.subplots(figsize=size)
        for candidate_id in LABELS:
            selected = sorted(
                (row for row in plot_rows if row["candidate_id"] == candidate_id),
                key=lambda row: float(row["snr_db"]),
            )
            x = np.asarray([float(row["snr_db"]) for row in selected])
            y = np.asarray([float(row["plot_bler"]) for row in selected])
            low = np.asarray([float(row["bler_wilson95_lo"]) for row in selected])
            high = np.asarray([float(row["bler_wilson95_hi"]) for row in selected])
            low = np.minimum(low, y)
            high = np.maximum(high, y)
            style = STYLES[candidate_id]
            axis.errorbar(
                x,
                y,
                yerr=np.vstack((y - low, high - y)),
                label=LABELS[candidate_id],
                color=style["color"],
                linestyle=style["linestyle"],
                marker=style["marker"],
                linewidth=2.5,
                markersize=8,
                markeredgewidth=1.0,
                capsize=3,
            )
        axis.set_yscale("log")
        axis.set_ylim(3e-4, 1.0)
        axis.axhline(0.1, color="#666666", linewidth=1.2, linestyle="--")
        axis.axhline(0.01, color="#666666", linewidth=1.2, linestyle=":")
        axis.set_xlabel("SNR (dB), unit total transmit power", fontsize=16)
        axis.set_ylabel("Transport-block error rate", fontsize=16)
        axis.set_title("8Tx/4Rx, 3 km/h: ideal-CSI BLER", fontsize=16)
        axis.grid(True, which="both", alpha=0.28)
        axis.tick_params(axis="both", labelsize=14)
        axis.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.20),
            ncol=2,
            fontsize=10 if preview else 13,
            frameon=False,
        )
        figure.tight_layout()
        figure.savefig(path, dpi=200, bbox_inches="tight")
        plt.close(figure)

    draw(output / "ideal_sidon_vs_prg6.png", preview=False)
    draw(output / "ideal_sidon_vs_prg6_preview_13cm.png", preview=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    rows = run(output)
    plot(output, rows)
    print(f"[done] data and figures: {output}", flush=True)


if __name__ == "__main__":
    main()
