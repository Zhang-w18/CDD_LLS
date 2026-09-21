"""One-command validation and resumable prescan for plan-031 strict Sidon sets."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
import platform
import socket
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search"
STATE_DIR = RUN_ROOT / "prescan_orchestration"
SHORTLIST_PATH = RUN_ROOT / "sidon_shortlist.json"
RUNNER = ROOT / "tools/run_pdcch_bler_curves.py"
PRESCAN_CONFIG_PATHS = (
    ROOT / "configs/pdcch_result031_strict_sidon_a100_al2_prescan.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_a100_al4_prescan.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_a100_al8_prescan.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_c300_al1_prescan.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_c300_al2_prescan.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_c300_al4_prescan.yaml",
)
EXTENSION_CONFIG_PATHS = (
    ROOT / "configs/pdcch_result031_strict_sidon_c300_al1_prescan_extension.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_c300_al2_prescan_extension.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_c300_al4_prescan_extension.yaml",
)
COARSE_CONFIG_PATHS = (
    ROOT / "configs/pdcch_result031_strict_sidon_a100_al2_coarse.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_a100_al4_coarse.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_a100_al8_coarse.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_c300_al1_coarse.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_c300_al2_coarse.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_c300_al4_coarse.yaml",
)
FINE_CONFIG_PATHS = (
    ROOT / "configs/pdcch_result031_strict_sidon_a100_al2_fine.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_a100_al4_fine.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_a100_al8_fine.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_c300_al1_fine.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_c300_al2_fine.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_c300_al4_fine.yaml",
)
FINE_SUPPLEMENT_CONFIG_PATHS = (
    ROOT / "configs/pdcch_result031_strict_sidon_a100_al4_fine_supplement.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_a100_al8_fine_supplement.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_c300_al2_fine_supplement.yaml",
    ROOT / "configs/pdcch_result031_strict_sidon_c300_al4_fine_supplement.yaml",
)
CONFIG_PATHS = PRESCAN_CONFIG_PATHS
EXPECTED_CANDIDATES_PER_SCENE = 8
EXPECTED_TOTAL_CANDIDATES = 48
EXPECTED_TRIALS_PER_POINT = 300
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
                "Another plan-031 strict-Sidon prescan is already running."
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


def _scene_name(config: dict) -> str:
    profile = str(config["channel"]["tdl_profile"]).upper()
    family = "a100" if profile == "A" else "c300"
    return f"{family}_al{int(config['resource']['aggregation_level'])}"


def _expected_shortlist() -> dict[tuple[str, int], list[list[int]]]:
    receipt = json.loads(SHORTLIST_PATH.read_text(encoding="utf-8"))
    if int(receipt["shortlist_size_per_scenario"]) != EXPECTED_CANDIDATES_PER_SCENE:
        raise RuntimeError("Unexpected shortlist size in sidon_shortlist.json")
    grouped: dict[tuple[str, int], list[tuple[int, list[int]]]] = {}
    for row in receipt["candidates"]:
        key = (str(row["family"]).upper(), int(row["aggregation_level"]))
        grouped.setdefault(key, []).append(
            (int(row["candidate_index"]), [int(value) for value in row["delay_indices"]])
        )
    return {
        key: [delays for _, delays in sorted(rows)] for key, rows in grouped.items()
    }


def _completed_points(
    output_dir: Path,
    candidate_id: str,
    snrs: list[float],
    minimum_trials: int,
    maximum_trials: int,
    target_errors: int,
) -> tuple[int, int]:
    csv_path = output_dir / candidate_id / "bler_points.csv"
    if not csv_path.exists():
        return 0, 0
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    seen = set()
    target_seen = set()
    for row in rows:
        if str(row["candidate_id"]) != candidate_id:
            raise RuntimeError(f"Unexpected candidate in {csv_path}")
        snr = float(row["snr_db"])
        if snr in seen:
            raise RuntimeError(f"Duplicate SNR in {csv_path}: {snr:g}")
        trials = int(row["trials"])
        errors = int(row["errors"])
        if not minimum_trials <= trials <= maximum_trials:
            raise RuntimeError(f"Unexpected trial count in {csv_path}: {trials}")
        if trials < maximum_trials and errors < target_errors:
            raise RuntimeError(
                f"Premature adaptive stop in {csv_path}: {errors}/{trials}"
            )
        seen.add(snr)
        if snr in snrs:
            target_seen.add(snr)
    completed_trials = sum(
        int(row["trials"]) for row in rows if float(row["snr_db"]) in snrs
    )
    return len(target_seen), completed_trials


def load_jobs(
    config_paths: tuple[Path, ...] = CONFIG_PATHS,
    *,
    expected_total_candidates: int = EXPECTED_TOTAL_CANDIDATES,
    require_full_shortlist: bool = True,
) -> tuple[list[dict], list[dict]]:
    expected = _expected_shortlist()
    scenes = []
    jobs = []
    for config_path in config_paths:
        config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
        scene = _scene_name(config)
        family = scene.split("_")[0].upper()
        aggregation = int(config["resource"]["aggregation_level"])
        candidates = list(config["candidates"])
        expected_delays = expected[(family, aggregation)]
        actual_delays = [row["delay_grid_coordinates"] for row in candidates]
        if require_full_shortlist:
            matches_shortlist = (
                len(candidates) == EXPECTED_CANDIDATES_PER_SCENE
                and actual_delays == expected_delays
            )
        else:
            matches_shortlist = bool(candidates)
            for candidate in candidates:
                try:
                    candidate_index = int(str(candidate["candidate_id"]).rsplit("_", 1)[1])
                except (IndexError, ValueError) as error:
                    raise RuntimeError(
                        f"Invalid strict-Sidon candidate ID in {config_path.name}"
                    ) from error
                matches_shortlist = matches_shortlist and (
                    1 <= candidate_index <= len(expected_delays)
                    and candidate["delay_grid_coordinates"]
                    == expected_delays[candidate_index - 1]
                )
        if not matches_shortlist:
            raise RuntimeError(f"{config_path.name} does not match sidon_shortlist.json")
        if any(not bool(row.get("strict_sidon")) for row in candidates):
            raise RuntimeError(f"{config_path.name} contains a non-strict candidate")
        default_snrs = [float(value) for value in config["simulation"]["snr_points_db"]]
        minimum_trials = int(config["simulation"]["min_trials_per_snr"])
        maximum_trials = int(config["simulation"]["max_trials_per_snr"])
        target_errors = int(config["simulation"]["target_errors"])
        output_dir = Path(config["output_dir"])
        if not output_dir.is_absolute():
            output_dir = ROOT / output_dir
        scene_jobs = []
        for candidate in candidates:
            candidate_id = str(candidate["candidate_id"])
            snrs = [
                float(value)
                for value in candidate.get("snr_points_db", default_snrs)
            ]
            complete, completed_trials = _completed_points(
                output_dir,
                candidate_id,
                snrs,
                minimum_trials,
                maximum_trials,
                target_errors,
            )
            job = {
                "name": f"{scene}:{candidate_id}",
                "scene": scene,
                "candidate_id": candidate_id,
                "config_path": config_path,
                "output_dir": output_dir / candidate_id,
                "snr_count": len(snrs),
                "snr_min_db": min(snrs),
                "snr_max_db": max(snrs),
                "trials_per_point": minimum_trials,
                "minimum_trials_per_point": minimum_trials,
                "maximum_trials_per_point": maximum_trials,
                "target_errors": target_errors,
                "completed_trials": completed_trials,
                "completed_points": complete,
                "remaining_points": len(snrs) - complete,
            }
            jobs.append(job)
            scene_jobs.append(job)
        scenes.append(
            {
                "name": scene,
                "config_path": config_path,
                "config_sha256": _sha256(config_path),
                "n_tx": int(config["antenna"]["n_tx"]),
                "duration_symbols": int(config["resource"]["duration_symbols"]),
                "channel": (
                    f"TDL-{str(config['channel']['tdl_profile']).upper()} "
                    f"{float(config['channel']['delay_spread_ns']):g} ns, "
                    f"{float(config['channel']['ue_speed_kmh']):g} km/h"
                ),
                "candidate_count": len(candidates),
                "snr_count_min": min(job["snr_count"] for job in scene_jobs),
                "snr_count_max": max(job["snr_count"] for job in scene_jobs),
                "snr_min_db": min(job["snr_min_db"] for job in scene_jobs),
                "snr_max_db": max(job["snr_max_db"] for job in scene_jobs),
                "trials_per_point": minimum_trials,
                "minimum_trials_per_point": minimum_trials,
                "maximum_trials_per_point": maximum_trials,
                "target_errors": target_errors,
                "completed_trials": sum(job["completed_trials"] for job in scene_jobs),
                "completed_points": sum(job["completed_points"] for job in scene_jobs),
                "total_points": sum(job["snr_count"] for job in scene_jobs),
            }
        )
    if len(jobs) != expected_total_candidates:
        raise RuntimeError(
            f"Expected exactly {expected_total_candidates} strict Sidon candidate jobs"
        )
    return scenes, jobs


def _serializable_scene(scene: dict) -> dict:
    return {
        **scene,
        "config_path": str(scene["config_path"].relative_to(ROOT)),
    }


def build_status(
    scenes: list[dict], jobs: list[dict], status: str, started_at: str, **extra
) -> dict:
    current = sum(job["completed_points"] for job in jobs)
    target = sum(job["snr_count"] for job in jobs)
    current_trials = sum(job["completed_trials"] for job in jobs)
    minimum_trials = sum(
        job["snr_count"] * job["minimum_trials_per_point"] for job in jobs
    )
    target_trials = sum(
        job["snr_count"] * job["maximum_trials_per_point"] for job in jobs
    )
    payload = {
        "schema": "plan031-strict-sidon-prescan-orchestration-v1",
        "status": status,
        "started_at_utc": started_at,
        "updated_at_utc": _utc_now(),
        "host": socket.gethostname(),
        "python": sys.executable,
        "checkpoint_granularity": "completed candidate/SNR point",
        "candidate_count": len(jobs),
        "point_progress": {"current": current, "target": target, "remaining": target - current},
        "candidate_trials": {
            "current": current_trials,
            "minimum_budget": minimum_trials,
            "maximum_budget": target_trials,
            "remaining_to_maximum": target_trials - current_trials,
        },
        "scenes": [_serializable_scene(scene) for scene in scenes],
    }
    payload.update(extra)
    return payload


def print_summary(
    scenes: list[dict], jobs: list[dict], max_workers: int, phase: str = "prescan"
) -> None:
    current = sum(job["completed_points"] for job in jobs)
    target = sum(job["snr_count"] for job in jobs)
    current_trials = sum(job["completed_trials"] for job in jobs)
    minimum_trials = sum(
        job["snr_count"] * job["minimum_trials_per_point"] for job in jobs
    )
    target_trials = sum(
        job["snr_count"] * job["maximum_trials_per_point"] for job in jobs
    )
    policies = sorted(
        {
            (
                job["minimum_trials_per_point"],
                job["target_errors"],
                job["maximum_trials_per_point"],
            )
            for job in jobs
        }
    )
    print(f"plan-031 strict Sidon candidate {phase}", flush=True)
    print(
        f"scope={len(scenes)} scene(s); {len(jobs)} non-transparent matched-LMMSE candidates",
        flush=True,
    )
    print(
        "budget="
        + ",".join(
            f"{minimum} min/{errors} errors/{maximum} max"
            for minimum, errors, maximum in policies
        )
        + " per SNR; "
        "BLER and per-trial data-RE CE NMSE saved",
        flush=True,
    )
    print(
        f"workers={max_workers}  checkpoint=candidate/SNR point  "
        f"points={current}/{target}  remaining={target-current}",
        flush=True,
    )
    print(
        f"candidate-trials={current_trials:,}  minimum-budget={minimum_trials:,}  "
        f"maximum-budget={target_trials:,}",
        flush=True,
    )
    print("scene       Tx sym candidates  SNR range/points  completed/total", flush=True)
    for scene in scenes:
        print(
            f"{scene['name']:<11} {scene['n_tx']:>2} {scene['duration_symbols']:>3} "
            f"{scene['candidate_count']:>10}  "
            f"{scene['snr_min_db']:g}..{scene['snr_max_db']:g}/"
            f"{scene['snr_count_min']}-{scene['snr_count_max']}  "
            f"{scene['completed_points']:>9}/{scene['total_points']:<5}",
            flush=True,
        )


def _environment_payload() -> dict:
    packages = {}
    for name in ("numpy", "scipy", "tensorflow", "sionna", "pyyaml"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    return {
        "created_at_utc": _utc_now(),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "python": sys.version,
        "executable": sys.executable,
        "packages": packages,
        "source_sha256": {
            str(path.relative_to(ROOT)): _sha256(path)
            for path in (
                Path(__file__).resolve(),
                RUNNER,
                ROOT / "cdd_lls/sim/pdcch_cdd.py",
                SHORTLIST_PATH,
            )
        },
    }


def validate_scenes(scenes: list[dict], stamp: str) -> None:
    validation_dir = STATE_DIR / "validation"
    log_dir = STATE_DIR / "logs"
    validation_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    print("Running frozen-input, pilot, CE-floor, and noiseless-codec validation...", flush=True)
    for scene in scenes:
        receipt = validation_dir / f"{scene['name']}.json"
        command = [
            sys.executable,
            str(RUNNER),
            "--config",
            str(scene["config_path"]),
            "--stage",
            "validate",
            "--validation-output",
            str(receipt),
        ]
        environment = dict(os.environ)
        environment.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
        completed = subprocess.run(
            command,
            cwd=ROOT,
            env=environment,
            text=True,
            encoding="utf-8",
            errors="replace",
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        (log_dir / f"{stamp}_validate_{scene['name']}.log").write_text(
            completed.stdout, encoding="utf-8"
        )
        if completed.returncode != 0 or not receipt.exists():
            print(completed.stdout.rstrip(), flush=True)
            raise RuntimeError(f"Validation failed for {scene['name']}")
        print(f"validated {scene['name']} -> {receipt.relative_to(ROOT)}", flush=True)


def _run_job(job: dict, stamp: str) -> tuple[str, int, str]:
    log_dir = STATE_DIR / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    safe_name = job["name"].replace(":", "_")
    log_path = log_dir / f"{stamp}_{safe_name}.log"
    command = [
        sys.executable,
        str(RUNNER),
        "--config",
        str(job["config_path"]),
        "--candidate",
        job["candidate_id"],
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
    parser.add_argument("--max-workers", type=int, default=3)
    parser.add_argument("--summary-only", action="store_true")
    parser.add_argument(
        "--phase",
        choices=(
            "prescan",
            "extension",
            "coarse-confirmation",
            "fine",
            "fine-supplement",
        ),
        default="prescan",
        help=(
            "Run the complete original prescan, frozen high-SNR extension, "
            "0.5 dB/3000-trial coarse confirmation, or final 0.25 dB fine scan."
        ),
    )
    args = parser.parse_args()
    if args.max_workers < 1 or args.max_workers > 3:
        parser.error("--max-workers must be between 1 and 3")

    started_at = _utc_now()
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    phase_settings = {
        "prescan": (RUN_ROOT / "prescan_orchestration", PRESCAN_CONFIG_PATHS, 48, True),
        "extension": (
            RUN_ROOT / "prescan_extension_orchestration",
            EXTENSION_CONFIG_PATHS,
            5,
            False,
        ),
        "coarse-confirmation": (
            RUN_ROOT / "coarse_confirmation_orchestration",
            COARSE_CONFIG_PATHS,
            43,
            False,
        ),
        "fine": (
            RUN_ROOT / "fine_orchestration",
            FINE_CONFIG_PATHS,
            18,
            False,
        ),
        "fine-supplement": (
            RUN_ROOT / "fine_supplement_orchestration",
            FINE_SUPPLEMENT_CONFIG_PATHS,
            10,
            False,
        ),
    }
    state_dir, config_paths, expected_candidates, require_full_shortlist = phase_settings[
        args.phase
    ]
    global STATE_DIR
    STATE_DIR = state_dir
    status_path = state_dir / "status.json"
    with _SingleInstanceLock(state_dir / "run.lock"):
        scenes, jobs = load_jobs(
            config_paths,
            expected_total_candidates=expected_candidates,
            require_full_shortlist=require_full_shortlist,
        )
        print_summary(scenes, jobs, args.max_workers, args.phase)
        _write_json_atomic(state_dir / "environment.json", _environment_payload())
        _write_json_atomic(
            status_path,
            build_status(
                scenes,
                jobs,
                "ready" if args.summary_only else "validating",
                started_at,
                max_workers=args.max_workers,
            ),
        )
        if args.summary_only:
            print(f"Summary saved to {status_path.relative_to(ROOT)}", flush=True)
            return 0

        validate_scenes(scenes, stamp)
        pending = [job for job in jobs if job["remaining_points"] > 0]
        if not pending:
            _write_json_atomic(
                status_path,
                build_status(scenes, jobs, "complete", started_at, max_workers=args.max_workers),
            )
            print("All prescan points are already complete.", flush=True)
            return 0

        _write_json_atomic(
            status_path,
            build_status(scenes, jobs, "running", started_at, max_workers=args.max_workers),
        )
        results = []
        executor = ThreadPoolExecutor(max_workers=args.max_workers)
        futures = {executor.submit(_run_job, job, stamp): job for job in pending}
        try:
            for future in as_completed(futures):
                results.append(future.result())
        except KeyboardInterrupt:
            print("\nInterrupt received; terminating active candidate processes...", flush=True)
            _terminate_children()
            executor.shutdown(wait=True, cancel_futures=True)
            refreshed_scenes, refreshed_jobs = load_jobs(
                config_paths,
                expected_total_candidates=expected_candidates,
                require_full_shortlist=require_full_shortlist,
            )
            _write_json_atomic(
                status_path,
                build_status(
                    refreshed_scenes,
                    refreshed_jobs,
                    "interrupted",
                    started_at,
                    max_workers=args.max_workers,
                    completed_processes=results,
                ),
            )
            print("Completed SNR points are preserved. Run the same command to resume.", flush=True)
            return 130
        finally:
            executor.shutdown(wait=True)

        refreshed_scenes, refreshed_jobs = load_jobs(
            config_paths,
            expected_total_candidates=expected_candidates,
            require_full_shortlist=require_full_shortlist,
        )
        failures = [result for result in results if result[1] != 0]
        incomplete = any(job["remaining_points"] > 0 for job in refreshed_jobs)
        final_status = "failed" if failures or incomplete else "complete"
        _write_json_atomic(
            status_path,
            build_status(
                refreshed_scenes,
                refreshed_jobs,
                final_status,
                started_at,
                max_workers=args.max_workers,
                completed_processes=results,
            ),
        )
        print_summary(refreshed_scenes, refreshed_jobs, args.max_workers, args.phase)
        if final_status == "complete":
            print(
                f"{args.phase} complete. BLER CSVs, trial error flags, CE NMSE arrays, "
                "resolved configs, validation receipts, logs, and environment metadata are saved.",
                flush=True,
            )
            return 0
        print(
            f"{args.phase} incomplete ({len(failures)} failed process(es)). "
            "Run the same command to resume.",
            flush=True,
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
