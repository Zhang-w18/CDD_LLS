"""Export compact, copyable BLER evidence from one or more completed runner configs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import yaml


ROOT = Path(__file__).resolve().parents[1]
PACKAGE_NAMES = ("numpy", "scipy", "matplotlib", "PyYAML", "tensorflow", "sionna")
CORE_CODE_PATHS = (
    "tools/run_bler_curves.py",
    "tools/run_v_design_piecewise_tradeoff.py",
    "cdd_lls/phy/channel_tdl.py",
    "cdd_lls/phy/ldpc.py",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _typed(value: str) -> Any:
    text = value.strip()
    if text == "":
        return None
    if text.lower() in {"true", "false"}:
        return text.lower() == "true"
    if re.fullmatch(r"[-+]?\d+", text):
        return int(text)
    try:
        return float(text)
    except ValueError:
        return text


def _read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [
            {key: _typed(value) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _artifact_receipt(paths: Iterable[Path], root: Path) -> dict[str, Any]:
    rows = []
    for path in sorted({item.resolve() for item in paths if item.is_file()}):
        try:
            display = path.relative_to(root.resolve()).as_posix()
        except ValueError:
            display = str(path)
        rows.append({"path": display, "bytes": path.stat().st_size, "sha256": _sha256(path)})
    aggregate = hashlib.sha256()
    for row in rows:
        aggregate.update(
            f"{row['path']}\0{row['bytes']}\0{row['sha256']}\n".encode("utf-8")
        )
    return {
        "file_count": len(rows),
        "total_bytes": sum(int(row["bytes"]) for row in rows),
        "aggregate_sha256": aggregate.hexdigest(),
        "files": rows,
    }


def _compact_raw_receipt(scene_output: Path, root: Path) -> dict[str, Any]:
    paths = sorted(scene_output.rglob("*.npy"))
    receipt = _artifact_receipt(paths, root)
    return {
        "file_count": receipt["file_count"],
        "total_bytes": receipt["total_bytes"],
        "aggregate_sha256": receipt["aggregate_sha256"],
    }


def _interval_audit(
    final_rows: list[dict[str, Any]], supplemental_rows: list[dict[str, Any]]
) -> dict[str, Any]:
    grouped: dict[tuple[str, float], list[tuple[int, int]]] = defaultdict(list)
    for row in supplemental_rows:
        grouped[(str(row["candidate_id"]), float(row["snr_db"]))].append(
            (int(row["trial_start"]), int(row["trial_end"]))
        )
    failures = []
    audited = []
    for row in final_rows:
        key = (str(row["candidate_id"]), float(row["snr_db"]))
        base_trials = int(row.get("base_trials") or 0)
        expected = base_trials + 1
        intervals = sorted(grouped.get(key, []))
        for start, end in intervals:
            if start != expected or end < start:
                failures.append(
                    {
                        "candidate_id": key[0],
                        "snr_db": key[1],
                        "expected_start": expected,
                        "actual_interval": [start, end],
                    }
                )
                break
            expected = end + 1
        observed_end = expected - 1
        if observed_end != int(row["trials"]):
            failures.append(
                {
                    "candidate_id": key[0],
                    "snr_db": key[1],
                    "expected_total_trials": int(row["trials"]),
                    "observed_interval_end": observed_end,
                }
            )
        audited.append(
            {
                "candidate_id": key[0],
                "snr_db": key[1],
                "base_trials": base_trials,
                "interval_count": len(intervals),
                "trial_end": observed_end,
            }
        )
    return {"ok": not failures, "failures": failures, "points": audited}


def environment_receipt() -> dict[str, Any]:
    packages = {}
    for name in PACKAGE_NAMES:
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    receipt: dict[str, Any] = {
        "python": sys.version,
        "executable": sys.executable,
        "platform": platform.platform(),
        "packages": packages,
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
        "tensorflow_device_policy": (
            "automatic placement; run_bler_curves.py has no explicit tf.device scope "
            "or multi-GPU distribution strategy"
        ),
    }
    try:
        import tensorflow as tf

        physical = tf.config.list_physical_devices()
        logical = tf.config.list_logical_devices()
        gpu_details = []
        for device in tf.config.list_physical_devices("GPU"):
            details = tf.config.experimental.get_device_details(device)
            gpu_details.append(
                {
                    "name": device.name,
                    "device_name": details.get("device_name"),
                    "compute_capability": details.get("compute_capability"),
                }
            )
        probe = tf.linalg.matmul(tf.ones((64, 64)), tf.ones((64, 64)))
        receipt["tensorflow"] = {
            "built_with_cuda": bool(tf.test.is_built_with_cuda()),
            "physical_devices": [device.name for device in physical],
            "logical_devices": [device.name for device in logical],
            "gpu_details": gpu_details,
            "gpu_visible": bool(gpu_details),
            "matmul_placement_probe": probe.device,
        }
    except Exception as exc:  # pragma: no cover - depends on host TensorFlow installation
        receipt["tensorflow"] = {"probe_error": f"{type(exc).__name__}: {exc}"}
    return receipt


def _git_head(root: Path) -> str | None:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=root,
        capture_output=True,
        text=True,
        check=False,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _progress_rows(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            rows.append(json.loads(line))
    return rows


def collect_config_results(config_path: Path, root: Path = ROOT) -> dict[str, Any]:
    config_path = config_path.resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    result: dict[str, Any] = {
        "config_path": str(config_path),
        "config_sha256": _sha256(config_path),
        "config": config,
        "scenes": [],
    }
    output_root = Path(config["output_dir"])
    if not output_root.is_absolute():
        output_root = root / output_root
    for scene in config["scenes"]:
        scene_output = output_root / str(scene["scenario_id"]).lower()
        resolved_path = scene_output / "resolved_run.json"
        if not resolved_path.exists():
            raise FileNotFoundError(f"Missing completed-run receipt: {resolved_path}")
        scene_result: dict[str, Any] = {
            "scenario_id": scene["scenario_id"],
            "n_rx": int(scene.get("n_rx", 1)),
            "mode": scene.get("mode", "cdd_manifest"),
            "resolved_run": _read_json(resolved_path),
            "progress": _progress_rows(scene_output / "run_progress.jsonl"),
            "receivers": {},
        }
        artifact_paths = [resolved_path, config_path]
        for manifest_name in ("transparent_prg_manifest.json", "transparent_cdd_manifest.json"):
            manifest_path = scene_output / "manifest" / manifest_name
            if manifest_path.exists():
                scene_result["transparent_manifest"] = _read_json(manifest_path)
                artifact_paths.append(manifest_path)
        for receiver in scene.get("receivers", {}):
            final_path = scene_output / "final" / f"{receiver}_csi_bler_points.csv"
            if not final_path.exists():
                raise FileNotFoundError(
                    f"Missing {final_path}; finish the run or run --stage merge before export"
                )
            supplemental_path = scene_output / receiver / "supplemental_points.csv"
            final_rows = _read_csv(final_path)
            supplemental_rows = _read_csv(supplemental_path) if supplemental_path.exists() else []
            scene_result["receivers"][receiver] = {
                "points": final_rows,
                "interval_audit": _interval_audit(final_rows, supplemental_rows),
            }
            artifact_paths.append(final_path)
            if supplemental_path.exists():
                artifact_paths.append(supplemental_path)
        scene_result["run_elapsed_seconds"] = sum(
            float(row.get("elapsed_seconds", 0.0)) for row in scene_result["progress"]
        )
        scene_result["raw_array_receipt"] = _compact_raw_receipt(scene_output, root)
        scene_result["key_artifact_receipt"] = _artifact_receipt(artifact_paths, root)
        result["scenes"].append(scene_result)
    return result


def build_bundle(config_paths: list[Path], root: Path = ROOT) -> dict[str, Any]:
    code_paths = [root / relative for relative in CORE_CODE_PATHS]
    return {
        "schema": "bler-copy-bundle-v1",
        "purpose": "compact evidence for result-document analysis; raw NPY arrays remain at source",
        "environment": environment_receipt(),
        "repository": {
            "root": str(root.resolve()),
            "git_head": _git_head(root),
            "core_code_receipt": _artifact_receipt(code_paths, root),
        },
        "runs": [collect_config_results(path, root) for path in config_paths],
    }


def _format_number(value: Any, digits: int = 6) -> str:
    if value is None:
        return "-"
    return f"{float(value):.{digits}g}"


def render_text(bundle: dict[str, Any]) -> str:
    tf_info = bundle["environment"].get("tensorflow", {})
    lines = [
        "BLER copy bundle",
        f"schema: {bundle['schema']}",
        f"git_head: {bundle['repository'].get('git_head')}",
        f"platform: {bundle['environment']['platform']}",
        f"python: {bundle['environment']['python'].splitlines()[0]}",
        f"tensorflow_gpu_visible: {tf_info.get('gpu_visible')}",
        f"tensorflow_built_with_cuda: {tf_info.get('built_with_cuda')}",
        f"tensorflow_matmul_probe: {tf_info.get('matmul_placement_probe')}",
        f"tensorflow_devices: {', '.join(tf_info.get('physical_devices', []))}",
        "",
        "Columns: run | scenario | nRx | receiver | candidate | SNR_dB | trials | errors | BLER | Wilson95_lo | Wilson95_hi | CE_NMSE_dB",
    ]
    for run_index, run in enumerate(bundle["runs"], start=1):
        for scene in run["scenes"]:
            lines.extend(
                [
                    "",
                    f"run {run_index}: {run['config_path']}",
                    f"mode: {scene['mode']}; elapsed_seconds: {scene['run_elapsed_seconds']:.3f}; interval_audits_ok: "
                    + str(
                        all(
                            receiver["interval_audit"]["ok"]
                            for receiver in scene["receivers"].values()
                        )
                    ),
                    f"raw_arrays: count={scene['raw_array_receipt']['file_count']}, bytes={scene['raw_array_receipt']['total_bytes']}, aggregate_sha256={scene['raw_array_receipt']['aggregate_sha256']}",
                ]
            )
            for receiver, receiver_data in scene["receivers"].items():
                for row in receiver_data["points"]:
                    lines.append(
                        " | ".join(
                            [
                                str(run_index),
                                str(scene["scenario_id"]),
                                str(scene["n_rx"]),
                                receiver,
                                str(row["candidate_id"]),
                                _format_number(row["snr_db"]),
                                str(row["trials"]),
                                str(row["tb_errors"]),
                                _format_number(row["bler"]),
                                _format_number(row.get("bler_wilson95_lo")),
                                _format_number(row.get("bler_wilson95_hi")),
                                _format_number(row.get("ce_nmse_mean_db")),
                            ]
                        )
                    )
    lines.extend(
        [
            "",
            "Interpretation note: 1000 fixed trials are suitable for a preliminary waterfall comparison. At BLER=1%, the expected error count is 10; a zero-error point has a Wilson 95% upper bound of about 0.00383.",
            "The JSON file is authoritative and includes configs, manifests, point-level confidence intervals, interval audits, environment/device data, timing, and hashes.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", action="append", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--output-txt", required=True, type=Path)
    args = parser.parse_args()
    bundle = build_bundle(args.config)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_txt.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(bundle, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    args.output_txt.write_text(render_text(bundle), encoding="utf-8")
    print(f"Wrote {args.output_json}")
    print(f"Wrote {args.output_txt}")


if __name__ == "__main__":
    main()
