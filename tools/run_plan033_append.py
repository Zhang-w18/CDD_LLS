"""Advance one recoverable plan-033 endpoint-append interval.

The frozen initial endpoint list is read from the partial analysis.  One call
advances every still-pending SNR in one scenario by at most 1,000 paired
absolute trials.  Repeat the same command until the reported status is
``complete``.  Existing trials are never restarted or copied.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.run_plan033_tdl_mobility_mimo import (  # noqa: E402
    load_config,
    run,
    validate_config,
)


TARGET_ERRORS = 200
MAX_TOTAL_TRIALS = 50_000
INTERVAL_TRIALS = 1_000
EXPECTED_REASON = "1pct_bracket_endpoint_in_0p5_to_2pct_with_lt_200_errors"
BASELINE_FILES = (
    "original_config.yaml",
    "resolved_run.json",
    "resolved_run.sha256",
    "candidate_receiver_manifest.json",
    "delay_audit.json",
    "seed_pairing_audit.json",
    "code_state.json",
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def endpoint_requirements(path: Path, scenario_id: str) -> list[dict[str, object]]:
    selected: list[dict[str, object]] = []
    seen: set[tuple[str, float]] = set()
    for row in _read_csv(path):
        if row["scenario_id"] != scenario_id:
            continue
        if row.get("reason") != EXPECTED_REASON:
            raise ValueError(f"Unexpected append reason for {scenario_id}: {row.get('reason')!r}")
        key = (row["candidate_id"], float(row["snr_db"]))
        if key in seen:
            raise ValueError(f"Duplicate append endpoint: {key}")
        seen.add(key)
        selected.append(
            {
                "scenario_id": scenario_id,
                "candidate_id": key[0],
                "snr_db": key[1],
                "initial_bler": float(row["bler"]),
                "initial_errors": int(row["errors"]),
                "initial_trials": int(row["trials"]),
                "reason": row["reason"],
            }
        )
    if not selected:
        raise ValueError(f"No frozen append endpoints found for {scenario_id} in {path}")
    return sorted(selected, key=lambda row: (float(row["snr_db"]), str(row["candidate_id"])))


def plan_next_step(
    requirements: list[dict[str, object]],
    final_rows: list[dict[str, str]],
) -> dict[str, object]:
    by_key = {
        (row["candidate_id"], float(row["snr_db"])): row
        for row in final_rows
    }
    trials_by_snr: dict[float, set[int]] = {}
    for row in final_rows:
        trials_by_snr.setdefault(float(row["snr_db"]), set()).add(int(row["trials"]))

    endpoints: list[dict[str, object]] = []
    pending_snrs: set[float] = set()
    for requirement in requirements:
        key = (str(requirement["candidate_id"]), float(requirement["snr_db"]))
        if key not in by_key:
            raise ValueError(f"Missing final point for frozen append endpoint: {key}")
        row = by_key[key]
        errors = int(row["tb_errors"])
        trials = int(row["trials"])
        if trials > MAX_TOTAL_TRIALS:
            raise ValueError(f"Endpoint exceeds {MAX_TOTAL_TRIALS} trials: {key} -> {trials}")
        complete = errors >= TARGET_ERRORS or trials >= MAX_TOTAL_TRIALS
        if not complete:
            pending_snrs.add(key[1])
        endpoints.append(
            {
                **requirement,
                "current_errors": errors,
                "current_trials": trials,
                "complete": complete,
                "completion_reason": (
                    "target_errors" if errors >= TARGET_ERRORS else
                    "max_total_trials" if trials >= MAX_TOTAL_TRIALS else
                    "pending"
                ),
            }
        )

    targets: dict[float, int] = {}
    for snr in sorted(pending_snrs):
        counts = trials_by_snr.get(snr, set())
        if len(counts) != 1:
            raise ValueError(f"Paired trial counts differ at {snr:g} dB: {sorted(counts)}")
        current = next(iter(counts))
        targets[snr] = min(current + INTERVAL_TRIALS, MAX_TOTAL_TRIALS)

    return {
        "schema": "plan033-endpoint-append-state-v1",
        "status": "complete" if not targets else "pending",
        "target_errors": TARGET_ERRORS,
        "max_total_trials": MAX_TOTAL_TRIALS,
        "next_target_total_trials_by_snr": targets,
        "endpoints": endpoints,
    }


def _request_config(base_path: Path, state: dict[str, object]) -> dict[str, object]:
    config = yaml.safe_load(base_path.read_text(encoding="utf-8")) or {}
    targets = state["next_target_total_trials_by_snr"]
    if not isinstance(targets, dict) or not targets:
        raise ValueError("No pending append interval to schedule.")
    config.pop("target_total_trials", None)
    config["snr_db"] = sorted(float(value) for value in targets)
    config["target_total_trials_by_snr"] = {
        float(snr): int(targets[snr]) for snr in config["snr_db"]
    }
    return config


def _preserve_baseline(output_dir: Path) -> Path:
    baseline_dir = output_dir / "append_scheduler" / "baseline"
    baseline_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = baseline_dir / "manifest.json"
    if manifest_path.exists():
        return manifest_path
    files = []
    for name in BASELINE_FILES:
        source = output_dir / name
        if not source.exists():
            raise FileNotFoundError(f"Cannot preserve missing pre-append evidence: {source}")
        target = baseline_dir / name
        shutil.copy2(source, target)
        files.append({"name": name, "sha256": _sha256(target), "bytes": target.stat().st_size})
    manifest_path.write_text(
        json.dumps({"schema": "plan033-pre-append-baseline-v1", "files": files}, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest_path


def _write_request(output_dir: Path, request: dict[str, object]) -> Path:
    payload = yaml.safe_dump(request, sort_keys=False, allow_unicode=True)
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    path = output_dir / "append_scheduler" / "requests" / f"request_{digest[:16]}.yaml"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.read_text(encoding="utf-8") != payload:
        raise RuntimeError(f"Append request hash collision: {path}")
    path.write_text(payload, encoding="utf-8")
    return path


def _state(base_path: Path, requirements_path: Path) -> tuple[dict[str, object], dict[str, object]]:
    base = load_config(base_path.resolve())
    output_dir = Path(base["output_dir"])
    final_path = output_dir / "final" / "estimated_csi_bler_points.csv"
    requirements = endpoint_requirements(requirements_path.resolve(), str(base["scenario_id"]))
    state = plan_next_step(requirements, _read_csv(final_path))
    state.update(
        {
            "scenario_id": base["scenario_id"],
            "base_config": str(base_path.resolve()),
            "base_config_sha256": _sha256(base_path.resolve()),
            "requirements_path": str(requirements_path.resolve()),
            "requirements_sha256": _sha256(requirements_path.resolve()),
        }
    )
    return base, state


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", required=True, type=Path)
    parser.add_argument("--requirements", required=True, type=Path)
    parser.add_argument("--stage", choices=("validate", "run"), default="validate")
    args = parser.parse_args()

    base, state = _state(args.base_config, args.requirements)
    output_dir = Path(base["output_dir"])
    if args.stage == "run" and state["status"] == "pending":
        baseline_manifest = _preserve_baseline(output_dir)
        request_path = _write_request(output_dir, _request_config(args.base_config, state))
        request = load_config(request_path)
        candidates, channel, grid, precoders, covariance = validate_config(request)
        run(request, candidates, channel, grid, precoders, covariance)
        _, state = _state(args.base_config, args.requirements)
        state["executed_request"] = str(request_path)
        state["baseline_manifest"] = str(baseline_manifest)

    state_path = output_dir / "append_scheduler" / "state.json"
    if args.stage == "run":
        state_path.parent.mkdir(parents=True, exist_ok=True)
        state_path.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(state, indent=2))
    if args.stage == "run":
        print(f"[written] {state_path}")


if __name__ == "__main__":
    main()
