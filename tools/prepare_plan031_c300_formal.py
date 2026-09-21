"""Freeze plan-031 C300 formal grids and write reproducibility receipts."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys

import numpy as np
import yaml


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

OUTPUT = ROOT / "outputs" / "experiment031_pdcch_cdd" / "20260907_c300_4tx_2sym"
TARGETS = (0.10, 0.01)


def isotonic_decreasing(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    blocks: list[list[float | int]] = []
    for index, (value, weight) in enumerate(zip(values, weights)):
        blocks.append([index, index + 1, float(value), float(weight)])
        while len(blocks) >= 2 and float(blocks[-2][2]) < float(blocks[-1][2]):
            right = blocks.pop()
            left = blocks.pop()
            total = float(left[3]) + float(right[3])
            mean = (float(left[2]) * float(left[3]) + float(right[2]) * float(right[3])) / total
            blocks.append([int(left[0]), int(right[1]), mean, total])
    out = np.empty(len(values), dtype=np.float64)
    for start, stop, value, _ in blocks:
        out[int(start) : int(stop)] = float(value)
    return out


def interpolate_target(snr: np.ndarray, probability: np.ndarray, target: float) -> float:
    for index in range(len(snr) - 1):
        if probability[index] >= target >= probability[index + 1]:
            if probability[index] == probability[index + 1]:
                return float((snr[index] + snr[index + 1]) / 2.0)
            fraction = (math.log10(target) - math.log10(probability[index])) / (
                math.log10(probability[index + 1]) - math.log10(probability[index])
            )
            return float(snr[index] + fraction * (snr[index + 1] - snr[index]))
    raise ValueError(f"Target {target} is not bracketed by the prescan.")


def freeze_grid(estimate: float) -> list[float]:
    center = round(float(estimate) * 4.0) / 4.0
    return [center + offset for offset in (-0.5, -0.25, 0.0, 0.25, 0.5)]


def read_rows(path: Path) -> list[dict[str, str]]:
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def candidate_receipt(al: int, config: dict, candidate: dict) -> dict[str, object]:
    k_active = 36 * int(al)
    n_p = k_active // 4
    base = {
        "aggregation_level": al,
        "candidate_id": str(candidate["candidate_id"]),
        "scheme": str(candidate["scheme"]),
        "receiver_covariance_mode": str(
            candidate.get(
                "receiver_covariance_mode",
                "physical_prg" if candidate["scheme"] == "reg_bundle_dft_cycling" else "matched_effective",
            )
        ),
        "k_active": k_active,
        "n_p": n_p,
        "q_ns": 1e9 / (k_active * float(config["resource"]["scs_khz"]) * 1e3),
    }
    if candidate["scheme"] == "reg_bundle_dft_cycling":
        return {**base, "cycling_order": list(candidate["cycling_order"])}
    delays = np.asarray(candidate["delay_grid_coordinates"], dtype=np.float64)
    residues = np.mod(delays, float(n_p))
    circular_distances = []
    pair_sums = []
    for left in range(len(delays)):
        for right in range(left + 1, len(delays)):
            delta = abs(float(residues[left] - residues[right]))
            circular_distances.append(min(delta, float(n_p) - delta))
        for right in range(left, len(delays)):
            pair_sums.append(float(np.mod(delays[left] + delays[right], k_active)))
    ordered_pair_sums = np.sort(np.asarray(pair_sums, dtype=np.float64))
    pair_gaps = np.diff(ordered_pair_sums)
    pair_wrap = ordered_pair_sums[0] + k_active - ordered_pair_sums[-1]
    minimum_pair_gap = float(min(np.min(pair_gaps), pair_wrap))
    return {
        **base,
        "delay_grid_coordinates": delays.tolist(),
        "delay_ns": (delays * float(base["q_ns"])).tolist(),
        "folded_residues": residues.tolist(),
        "minimum_folded_circular_distance": min(circular_distances),
        "minimum_pair_sum_circular_distance": minimum_pair_gap,
        "pair_sum_count": len(pair_sums),
        "unique_pair_sum_count": len(set(round(value, 12) for value in pair_sums)),
        "strict_sidon_pair_sums": len(set(round(value, 12) for value in pair_sums)) == len(pair_sums),
    }


def write_pdp_receipt() -> None:
    from sionna.phy.channel.tr38901 import TDL

    tdl = TDL(
        model="C",
        delay_spread=300e-9,
        carrier_frequency=4.0e9,
        num_sinusoids=20,
        min_speed=3.0 / 3.6,
        max_speed=3.0 / 3.6,
        num_rx_ant=1,
        num_tx_ant=1,
        precision="double",
    )
    delays_ns = np.asarray(tdl.delays.numpy(), dtype=np.float64).reshape(-1) * 1e9
    powers = np.asarray(tdl.mean_powers.numpy(), dtype=np.float64).reshape(-1)
    powers /= np.sum(powers)
    order = np.argsort(delays_ns)
    cumulative = np.cumsum(powers[order])
    index = int(np.flatnonzero(cumulative >= 0.99)[0])
    receipt = {
        "profile": "TDL-C",
        "delay_spread_ns": 300.0,
        "tap_count": int(len(delays_ns)),
        "support_definition": "smallest sorted tap-delay interval from zero containing at least 99% power",
        "t_0p01_ns": float(delays_ns[order][index]),
        "contained_power": float(cumulative[index]),
        "tap_delays_ns": delays_ns.tolist(),
        "normalized_tap_powers": powers.tolist(),
        "speed_mps": 3.0 / 3.6,
        "maximum_doppler_hz": (3.0 / 3.6) * 4.0e9 / 299792458.0,
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "pdp_t0p01_receipt.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    freeze_rows: list[dict[str, object]] = []
    candidate_rows: list[dict[str, object]] = []
    for al in (1, 2, 4):
        prescan_path = ROOT / "configs" / f"pdcch_result031_c300_al{al}_prescan.yaml"
        with open(prescan_path, "r", encoding="utf-8") as handle:
            config = yaml.safe_load(handle) or {}
        formal = json.loads(json.dumps(config))
        formal["output_dir"] = (
            f"outputs/experiment031_pdcch_cdd/20260907_c300_4tx_2sym/formal/al{al}"
        )
        formal["simulation"].update(
            {
                "snr_points_db": [0.0],
                "min_trials_per_snr": 10000,
                "target_errors": 200,
                "max_trials_per_snr": 50000,
                "resume": True,
            }
        )
        for candidate, formal_candidate in zip(config["candidates"], formal["candidates"]):
            candidate_id = str(candidate["candidate_id"])
            rows = read_rows(ROOT / config["output_dir"] / candidate_id / "bler_points.csv")
            snr = np.asarray([float(row["snr_db"]) for row in rows], dtype=np.float64)
            trials = np.asarray([int(row["trials"]) for row in rows], dtype=np.float64)
            errors = np.asarray([int(row["errors"]) for row in rows], dtype=np.float64)
            fitted = isotonic_decreasing((errors + 0.5) / (trials + 1.0), trials)
            grids: list[float] = []
            for target in TARGETS:
                try:
                    estimate = interpolate_target(snr, fitted, target)
                    local_grid = freeze_grid(estimate)
                    status = "prescan_bracketed"
                except ValueError:
                    if target != 0.01 or fitted[-1] < target:
                        raise
                    estimate = None
                    local_grid = [float(snr[-1]) - 2.0 + 0.25 * index for index in range(9)]
                    status = "prescan_unbracketed_high_snr_floor_probe"
                grids.extend(local_grid)
                freeze_rows.append(
                    {
                        "aggregation_level": al,
                        "candidate_id": candidate_id,
                        "target_bler": target,
                        "prescan_estimate_snr_db": estimate,
                        "frozen_local_grid_db": local_grid,
                        "status": status,
                    }
                )
            formal_candidate["snr_points_db"] = sorted(set(grids))
            candidate_rows.append(candidate_receipt(al, config, candidate))
        formal_path = ROOT / "configs" / f"pdcch_result031_c300_al{al}_formal.yaml"
        formal_path.write_text(
            yaml.safe_dump(formal, allow_unicode=True, sort_keys=False), encoding="utf-8"
        )
    (OUTPUT / "prescan_freeze_receipt.json").write_text(
        json.dumps(freeze_rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (OUTPUT / "candidate_receipt.json").write_text(
        json.dumps(candidate_rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    write_pdp_receipt()
    print(json.dumps({"frozen_target_rows": len(freeze_rows), "candidate_rows": len(candidate_rows)}))


if __name__ == "__main__":
    main()
