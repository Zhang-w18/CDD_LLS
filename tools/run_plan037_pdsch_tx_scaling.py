"""Run plan-037 4/16Tx, 1Rx PDSCH scaling curves and audit frozen geometry.

Strict 32Tx Sidon is analytically infeasible for K=576.  This entry point
therefore rejects every 32Tx run config; use ``--stage audit`` to emit the
proof and all frozen 4/16Tx geometry without starting Sionna.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.phy.precoding import build_prg_dft_precoder
from tools import run_plan033_tdl_mobility_mimo as core


SCHEMA = "plan037-pdsch-tx-scaling-v1"
PLAN_PATH = "research/plan-037-PDSCH-Tx数影响.md"
K = 576
FOLD_PERIOD = 96

FROZEN_DELAYS = {
    4: {
        "B0_QC": [0, 18, 36, 54],
        "S0_SIDON": [0, 24, 72, 240],
    },
    16: {
        "B0_QC_GROUPED": [0, 9, 18, 27, 36, 45, 54, 63] * 2,
        "S0_SIDON": [0, 9, 56, 61, 71, 174, 243, 280, 303, 307, 388, 407, 419, 427, 451, 509],
    },
    32: {
        "B0_QC_GROUPED": [0, 9, 18, 27, 36, 45, 54, 63] * 4,
    },
}


def prg_vector_indices(n_tx: int) -> list[int]:
    if int(n_tx) == 4:
        return [0, 1, 2, 3, 0, 1, 2, 3]
    if int(n_tx) in {16, 32}:
        return list(range(8))
    raise ValueError("Plan-037 n_tx must be 4, 16, or 32.")


def candidate_definitions(n_tx: int) -> list[dict]:
    n_tx = int(n_tx)
    if n_tx == 32:
        raise ValueError(
            "Plan-037 rejects 32Tx runs: strict Sidon is infeasible because "
            "32*31=992 distinct nonzero directed differences exceed 575."
        )
    if n_tx not in {4, 16}:
        raise ValueError("Plan-037 runnable n_tx must be 4 or 16.")
    prefix = f"A100_NT{n_tx}"
    delays = FROZEN_DELAYS[n_tx]
    rows = []
    for label, values in delays.items():
        rows.append({
            "candidate_id": f"{prefix}_{label}", "family": "CDD_MATCHED",
            "label": label, "delay_grid_coordinates": list(values),
            "receiver_knowledge": "matched CDD effective-channel covariance",
        })
    rows.append({
        "candidate_id": f"{prefix}_TRANSPARENT_PRG6",
        "family": "TRANSPARENT_PRG_DFT", "label": "transparent PRG6",
        "prg_vector_indices": prg_vector_indices(n_tx),
        "receiver_knowledge": "physical covariance within each PRG",
    })
    return rows


def load_config(path: Path) -> dict:
    config = core.load_config(path, expected_schema=SCHEMA)
    config["plan_path"] = PLAN_PATH
    return config


def _cyclic_minimum_gap(values: Sequence[int], modulus: int) -> int:
    unique = sorted({int(value) % modulus for value in values})
    if len(unique) < 2:
        return 0
    return min((unique[(index + 1) % len(unique)] - unique[index]) % modulus for index in range(len(unique)))


def audit_delay_set(n_tx: int, label: str, values: Sequence[int]) -> dict:
    delays = [int(value) for value in values]
    if len(delays) != int(n_tx):
        raise ValueError(f"{label} must contain exactly {n_tx} delays.")
    pair_sums = sorted((delays[i] + delays[j]) % K for i in range(n_tx) for j in range(i, n_tx))
    residues = [value % FOLD_PERIOD for value in delays]
    pilot = np.exp(-1j * 2.0 * np.pi * np.arange(0, K, 6)[:, None] * np.asarray(delays)[None, :] / K)
    singular_values = np.linalg.svd(pilot, compute_uv=False)
    rank = int(np.linalg.matrix_rank(pilot, tol=1e-10))
    return {
        "n_tx": n_tx, "label": label, "delay_grid_coordinates": delays,
        "delay_ns": [value / (K * 30e3) * 1e9 for value in delays],
        "delay_fft_samples": [value * 4096.0 / K for value in delays],
        "comb6_residues_mod_96": residues,
        "unordered_pair_sums_mod_576": pair_sums,
        "unique_pair_sum_count": len(set(pair_sums)),
        "strict_sidon": len(set(pair_sums)) == len(pair_sums),
        "minimum_pair_sum_gap": _cyclic_minimum_gap(pair_sums, K),
        "minimum_fold_gap": _cyclic_minimum_gap(residues, FOLD_PERIOD),
        "pilot_rank": rank,
        "pilot_singular_values": singular_values.tolist(),
        "pilot_condition_number": float(np.linalg.cond(pilot)) if rank == n_tx else "inf",
    }


def geometry_audit() -> dict:
    delay_rows = [
        audit_delay_set(n_tx, label, values)
        for n_tx, sets in FROZEN_DELAYS.items()
        for label, values in sets.items()
    ]
    return {
        "schema": "plan037-geometry-audit-v1",
        "delay_sets": delay_rows,
        "prg_vector_indices": {str(n_tx): prg_vector_indices(n_tx) for n_tx in (4, 16, 32)},
        "sidon_32tx": {
            "status": "NOT_FOUND_PROVEN_INFEASIBLE",
            "required_distinct_nonzero_directed_differences": 32 * 31,
            "available_nonzero_residues_mod_576": K - 1,
            "inequality_holds": 32 * 31 <= K - 1,
        },
    }


def validate_config(config: dict):
    if config["n_tx"] == 32:
        candidate_definitions(32)
    if config["n_tx"] not in {4, 16} or config["n_rx"] != 1 or config["speed_kmh"] != 60.0:
        raise ValueError("Plan-037 runnable scope is 4Tx/16Tx, 1Rx, 60 km/h.")
    if config["receiver_modes"] != ["estimated", "ideal"]:
        raise ValueError("Plan-037 requires receiver_modes: [estimated, ideal].")
    if config["channel_model"] != "tdl_a":
        raise ValueError("Plan-037 requires TDL-A.")
    if config["scenario_id"] != core.scenario_id(config["n_tx"], 60.0, 1):
        raise ValueError("Plan-037 scenario_id does not match n_tx/n_rx/speed.")
    candidates = candidate_definitions(config["n_tx"])
    config["_candidate_manifest"] = candidates
    # Reuse the common numeric/budget checks while keeping the 037 antenna scope local.
    if config["run_kind"] not in {"smoke", "prescan", "formal"}:
        raise ValueError("run_kind must be smoke, prescan, or formal")
    targets = [core._trial_target(config, value) for value in config["snr_db"]]
    expected = 20 if config["run_kind"] == "smoke" else 400 if config["run_kind"] == "prescan" else None
    if expected is not None and targets != [expected] * len(targets):
        raise ValueError(f"Plan-037 {config['run_kind']} requires {expected} trials per point.")
    if not config["snr_db"] or any(b <= a for a, b in zip(config["snr_db"], config["snr_db"][1:])):
        raise ValueError("snr_db must be non-empty and strictly increasing.")
    if min(config["snr_db"]) < -2.0 or max(config["snr_db"]) > 22.0:
        raise ValueError("Plan-037 SNR must remain within [-2,22] dB.")
    if config["batch_size"] <= 0 or any(value <= 0 for value in targets):
        raise ValueError("Batch size and trial targets must be positive.")
    if config["run_kind"] != "smoke" and any(value % config["batch_size"] for value in targets):
        raise ValueError("Non-smoke trial targets must align to batch_size.")
    channel, grid, precoders, base = core._build_scene(config)
    for row in candidates:
        if row["family"] == "CDD_MATCHED":
            power = np.sum(np.abs(precoders[row["candidate_id"]].C) ** 2, axis=1)
            if not np.allclose(power, 1.0, atol=1e-12):
                raise RuntimeError(f"Precoder power changed for {row['candidate_id']}.")
        else:
            prg = build_prg_dft_precoder(grid, config["n_tx"], 6, row["prg_vector_indices"], normalize=True)
            if not np.allclose(np.sum(np.abs(prg.C) ** 2, axis=1), 1.0, atol=1e-12):
                raise RuntimeError("PRG6 unit-power audit failed.")
    print(f"[validate] {config['scenario_id']}: curves=3 snr_points={len(config['snr_db'])} trials={targets[0]} batch={config['batch_size']}", flush=True)
    return candidates, channel, grid, precoders, base


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--stage", choices=("audit", "validate", "run", "merge", "adaptive"), default="run")
    parser.add_argument("--audit-output", type=Path)
    args = parser.parse_args()
    if args.stage == "audit":
        payload = geometry_audit()
        text = json.dumps(payload, indent=2, allow_nan=False) + "\n"
        if args.audit_output:
            output = args.audit_output.resolve()
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(text, encoding="utf-8")
        else:
            print(text, end="")
        return
    if args.config is None:
        parser.error("--config is required unless --stage audit is selected")
    config = load_config(args.config.resolve())
    candidates, channel, grid, precoders, base = validate_config(config)
    if args.stage == "run":
        core.run(config, candidates, channel, grid, precoders, base)
    elif args.stage == "adaptive":
        if "adaptive_policy" not in config:
            raise ValueError("--stage adaptive requires adaptive_policy in the config.")
        core.run_adaptive(config, candidates, channel, grid, precoders, base)
    elif args.stage == "merge":
        core._merge_intervals(Path(config["output_dir"]), candidates, config["scenario_id"])


if __name__ == "__main__":
    main()
