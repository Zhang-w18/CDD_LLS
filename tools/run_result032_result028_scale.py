"""One-command, resumable orchestration for the result-032 trial extension."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import socket
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "outputs/experiment032_tdl_mobility/20260906_main"
STATE_DIR = RUN_ROOT / "result028_scale_orchestration"
CONFIG_PATHS = (
    ROOT / "configs/bler_curves_result032_a100_60kmh_result028_scale_formal.yaml",
    ROOT / "configs/bler_curves_result032_a100_60kmh_result028_scale_shard1.yaml",
    ROOT / "configs/bler_curves_result032_a100_60kmh_result028_scale_shard2.yaml",
    ROOT / "configs/bler_curves_result032_a100_60kmh_result028_scale_shard3.yaml",
)
EXPECTED_CANDIDATES = 12
SIMULATION_SUMMARY = {
    "scenario": "A100_V60",
    "channel": "Sionna TDL-A, RMS delay spread 100 ns",
    "mobility": "60 km/h at 3.5 GHz; no ICI/CFO",
    "link": "8 Tx / 1 Rx, single layer, 48 PRB, 30 kHz SCS",
    "receiver": "two-DMRS [2,7] joint 2D time-frequency RMMSE",
    "curves": "10 matched CDD + transparent PRG6 + aged-CSI MRT PRG6",
    "budget": "14-20 dB >=10000 trials; 0.5%-2% BLER targets 200 errors; cap 50000",
}
_PRINT_LOCK = threading.Lock()
_CHILD_LOCK = threading.Lock()
_ACTIVE_CHILDREN: dict[str, subprocess.Popen[str]] = {}


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


class _SingleInstanceLock:
    """An OS-released lock; a stale file after a crash does not block resume."""

    def __init__(self, path: Path):
        self.path = path
        self.handle = None

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = self.path.open("a+b")
        self.handle.seek(0, os.SEEK_END)
        if self.handle.tell() == 0:
            self.handle.write(b"\0")
            self.handle.flush()
        self.handle.seek(0)
        try:
            if os.name == "nt":
                import msvcrt

                msvcrt.locking(self.handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl

                fcntl.flock(self.handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            self.handle.close()
            raise RuntimeError(
                "Another result-032 orchestration process is already running."
            ) from error
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.handle is None:
            return
        self.handle.seek(0)
        if os.name == "nt":
            import msvcrt

            msvcrt.locking(self.handle.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl

            fcntl.flock(self.handle.fileno(), fcntl.LOCK_UN)
        self.handle.close()


def read_completed_trials(output_dir: Path) -> dict[float, int]:
    """Return paired completed trials per SNR and reject gaps or partial curves."""
    interval_path = output_dir / "intervals.csv"
    if not interval_path.exists():
        return {}
    with interval_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    grouped: dict[tuple[str, float], list[tuple[int, int]]] = {}
    for row in rows:
        key = (str(row["candidate_id"]), float(row["snr_db"]))
        grouped.setdefault(key, []).append(
            (int(row["trial_start"]), int(row["trial_end"]))
        )
    totals: dict[float, list[int]] = {}
    for (_, snr_db), intervals in grouped.items():
        expected_start = 1
        for trial_start, trial_end in sorted(intervals):
            if trial_start != expected_start or trial_end < trial_start:
                raise RuntimeError(
                    f"Trial gap/overlap in {interval_path} at {snr_db:g} dB"
                )
            expected_start = trial_end + 1
        totals.setdefault(snr_db, []).append(expected_start - 1)
    completed = {}
    for snr_db, candidate_totals in totals.items():
        if len(candidate_totals) != EXPECTED_CANDIDATES:
            raise RuntimeError(
                f"Partial candidate checkpoint in {interval_path} at {snr_db:g} dB"
            )
        if len(set(candidate_totals)) != 1:
            raise RuntimeError(
                f"Paired trial counts differ in {interval_path} at {snr_db:g} dB"
            )
        completed[snr_db] = candidate_totals[0]
    return completed


def load_jobs() -> list[dict]:
    jobs = []
    for config_path in CONFIG_PATHS:
        config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
        if config.get("schema") != "plan032-tdl-mobility-v1":
            raise RuntimeError(f"Unexpected schema in {config_path}")
        snr_values = [float(value) for value in config["snr_db"]]
        raw_targets = config.get("target_total_trials_by_snr")
        if not isinstance(raw_targets, dict):
            raise RuntimeError(f"Per-SNR trial targets are required in {config_path}")
        targets = {float(key): int(value) for key, value in raw_targets.items()}
        if set(targets) != set(snr_values):
            raise RuntimeError(f"SNR and target keys differ in {config_path}")
        output_dir = Path(config["output_dir"])
        if not output_dir.is_absolute():
            output_dir = ROOT / output_dir
        completed = read_completed_trials(output_dir)
        points = []
        for snr_db in snr_values:
            current = int(completed.get(snr_db, 0))
            target = targets[snr_db]
            if current > target:
                raise RuntimeError(
                    f"Existing trials exceed target at {snr_db:g} dB in {config_path}"
                )
            points.append(
                {
                    "snr_db": snr_db,
                    "current_trials": current,
                    "target_trials": target,
                    "remaining_trials": target - current,
                }
            )
        jobs.append(
            {
                "name": config_path.stem.replace(
                    "bler_curves_result032_a100_60kmh_result028_scale_", ""
                ),
                "config_path": config_path,
                "config_sha256": _sha256(config_path),
                "output_dir": output_dir,
                "points": points,
            }
        )
    return jobs


def build_status(jobs: list[dict], status: str, started_at: str, **extra) -> dict:
    total_current = sum(
        point["current_trials"] for job in jobs for point in job["points"]
    )
    total_target = sum(
        point["target_trials"] for job in jobs for point in job["points"]
    )
    payload = {
        "schema": "result032-result028-scale-orchestration-v1",
        "status": status,
        "started_at_utc": started_at,
        "updated_at_utc": _utc_now(),
        "host": socket.gethostname(),
        "python": sys.executable,
        "checkpoint_granularity": "completed SNR point",
        "candidate_count": EXPECTED_CANDIDATES,
        "simulation": SIMULATION_SUMMARY,
        "common_trials": {
            "current": total_current,
            "target": total_target,
            "remaining": total_target - total_current,
        },
        "candidate_trials_remaining": (total_target - total_current)
        * EXPECTED_CANDIDATES,
        "jobs": [
            {
                "name": job["name"],
                "config_path": str(job["config_path"].relative_to(ROOT)),
                "config_sha256": job["config_sha256"],
                "output_dir": str(job["output_dir"].relative_to(ROOT)),
                "points": job["points"],
            }
            for job in jobs
        ],
    }
    payload.update(extra)
    return payload


def print_summary(jobs: list[dict], max_workers: int) -> None:
    total_current = sum(
        point["current_trials"] for job in jobs for point in job["points"]
    )
    total_target = sum(
        point["target_trials"] for job in jobs for point in job["points"]
    )
    print("result-032 60 km/h trial extension", flush=True)
    print(
        f"scenario={SIMULATION_SUMMARY['scenario']}  "
        f"channel={SIMULATION_SUMMARY['channel']}",
        flush=True,
    )
    print(
        f"link={SIMULATION_SUMMARY['link']}  receiver={SIMULATION_SUMMARY['receiver']}",
        flush=True,
    )
    print(f"curves={SIMULATION_SUMMARY['curves']}", flush=True)
    print(f"budget={SIMULATION_SUMMARY['budget']}", flush=True)
    print(
        f"checkpoint=SNR point  workers={max_workers}  candidates={EXPECTED_CANDIDATES}",
        flush=True,
    )
    print(
        f"common trials: current={total_current:,} target={total_target:,} "
        f"remaining={total_target-total_current:,}",
        flush=True,
    )
    print(
        f"equivalent remaining candidate-trials="
        f"{(total_target-total_current)*EXPECTED_CANDIDATES:,}",
        flush=True,
    )
    print("job       SNR points  complete  current/target common trials", flush=True)
    for job in jobs:
        complete = sum(
            point["remaining_trials"] == 0 for point in job["points"]
        )
        current = sum(point["current_trials"] for point in job["points"])
        target = sum(point["target_trials"] for point in job["points"])
        print(
            f"{job['name']:<10} {len(job['points']):>10} {complete:>9}  "
            f"{current:>8,}/{target:<8,}",
            flush=True,
        )
    print(
        "Progress below is emitted by each shard; completed SNR points are checkpointed.",
        flush=True,
    )


def _run_job(job: dict, stamp: str) -> tuple[str, int, str]:
    log_dir = STATE_DIR / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{stamp}_{job['name']}.log"
    command = [
        sys.executable,
        str(ROOT / "tools/run_plan032_tdl_mobility.py"),
        "--config",
        str(job["config_path"]),
        "--stage",
        "run",
    ]
    environment = dict(os.environ)
    environment.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
    matplotlib_dir = STATE_DIR / "matplotlib"
    matplotlib_dir.mkdir(parents=True, exist_ok=True)
    environment.setdefault("MPLCONFIGDIR", str(matplotlib_dir))
    with log_path.open("a", encoding="utf-8", buffering=1) as log:
        log.write(f"[{_utc_now()}] command={json.dumps(command)}\n")
        process = subprocess.Popen(
            command,
            cwd=ROOT,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )
        with _CHILD_LOCK:
            _ACTIVE_CHILDREN[job["name"]] = process
        try:
            assert process.stdout is not None
            for line in process.stdout:
                log.write(line)
                with _PRINT_LOCK:
                    print(f"[{job['name']}] {line}", end="", flush=True)
            return_code = process.wait()
        finally:
            with _CHILD_LOCK:
                _ACTIVE_CHILDREN.pop(job["name"], None)
        log.write(f"[{_utc_now()}] exit_code={return_code}\n")
    return job["name"], return_code, str(log_path.relative_to(ROOT))


def _terminate_children() -> None:
    with _CHILD_LOCK:
        children = list(_ACTIVE_CHILDREN.values())
    for process in children:
        if process.poll() is None:
            process.terminate()
    for process in children:
        try:
            process.wait(timeout=15)
        except subprocess.TimeoutExpired:
            process.kill()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--max-workers",
        type=int,
        default=3,
        help="Maximum concurrent shard processes (default: 3).",
    )
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="Print and save the resumable progress summary without running trials.",
    )
    args = parser.parse_args()
    if args.max_workers < 1 or args.max_workers > 3:
        parser.error("--max-workers must be between 1 and 3")

    started_at = _utc_now()
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    status_path = STATE_DIR / "status.json"
    with _SingleInstanceLock(STATE_DIR / "run.lock"):
        jobs = load_jobs()
        print_summary(jobs, args.max_workers)
        _write_json_atomic(
            status_path,
            build_status(
                jobs,
                "ready" if args.summary_only else "running",
                started_at,
                max_workers=args.max_workers,
            ),
        )
        if args.summary_only:
            print(f"Summary saved to {status_path.relative_to(ROOT)}", flush=True)
            return 0
        pending = [
            job
            for job in jobs
            if any(point["remaining_trials"] > 0 for point in job["points"])
        ]
        if not pending:
            _write_json_atomic(
                status_path,
                build_status(jobs, "complete", started_at, max_workers=args.max_workers),
            )
            print("All requested SNR points are already complete.", flush=True)
            return 0

        results = []
        executor = ThreadPoolExecutor(max_workers=args.max_workers)
        futures = {executor.submit(_run_job, job, stamp): job for job in pending}
        try:
            for future in as_completed(futures):
                results.append(future.result())
        except KeyboardInterrupt:
            print("\nInterrupt received; terminating active shards...", flush=True)
            _terminate_children()
            executor.shutdown(wait=True, cancel_futures=True)
            refreshed = load_jobs()
            _write_json_atomic(
                status_path,
                build_status(
                    refreshed,
                    "interrupted",
                    started_at,
                    max_workers=args.max_workers,
                    completed_processes=results,
                ),
            )
            print(
                "Completed SNR checkpoints are preserved. Run the same command to resume.",
                flush=True,
            )
            return 130
        finally:
            executor.shutdown(wait=True)

        refreshed = load_jobs()
        failures = [result for result in results if result[1] != 0]
        incomplete = any(
            point["remaining_trials"] > 0
            for job in refreshed
            for point in job["points"]
        )
        final_status = "failed" if failures or incomplete else "complete"
        _write_json_atomic(
            status_path,
            build_status(
                refreshed,
                final_status,
                started_at,
                max_workers=args.max_workers,
                completed_processes=results,
            ),
        )
        print_summary(refreshed, args.max_workers)
        if final_status == "complete":
            print(
                "Simulation complete. Raw intervals, error flags, CE arrays, and shard CSVs "
                "are ready for analysis.",
                flush=True,
            )
            return 0
        print(
            f"Run incomplete ({len(failures)} failed shard process(es)). "
            "Run the same command to resume.",
            flush=True,
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
