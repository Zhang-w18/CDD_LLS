"""Migrate pre-freeze plan-033 curve IDs without changing numeric arrays."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ALLOWED_ROOT = (ROOT / "outputs" / "experiment033_tdl_mobility_mimo").resolve()
RECEIPT_NAME = "curve_id_migration_receipt.json"
TEXT_SUFFIXES = {".csv", ".json", ".jsonl", ".md", ".sha256", ".txt", ".yaml", ".yml"}
REPLACEMENTS = {
    "A100_NT4_AGED_CSI_MRT_PRG6_SLOTS10": "A100_NT4_AGED_MRT_PRG6",
    "A100_NT8_AGED_CSI_MRT_PRG6_SLOTS10": "A100_NT8_AGED_MRT_PRG6",
    "A100_NT4_SMALL_CDD_QSTEP0P25_NONTRANSPARENT": "A100_NT4_SMALL_CDD_QSTEP0P25_MATCHED",
    "A100_NT8_SMALL_CDD_QSTEP0P25_NONTRANSPARENT": "A100_NT8_SMALL_CDD_QSTEP0P25_MATCHED",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _replace(value: str) -> str:
    result = value
    for old, new in REPLACEMENTS.items():
        result = result.replace(old, new)
    return result


def _within(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def _atomic_write(path: Path, content: str) -> None:
    temporary = path.with_suffix(path.suffix + ".migration-tmp")
    temporary.write_text(content, encoding="utf-8")
    temporary.replace(path)


def _inventory(target: Path) -> tuple[list[dict], list[dict]]:
    text_changes = []
    renames = []
    for path in sorted(item for item in target.rglob("*") if item.is_file()):
        if path.name == RECEIPT_NAME:
            continue
        if path.suffix.lower() in TEXT_SUFFIXES:
            content = path.read_text(encoding="utf-8")
            replaced = _replace(content)
            if replaced != content:
                text_changes.append(
                    {
                        "path": str(path.relative_to(ROOT)),
                        "replacement_count": sum(content.count(old) for old in REPLACEMENTS),
                    }
                )
        new_name = _replace(path.name)
        if new_name != path.name:
            destination = path.with_name(new_name)
            if destination.exists():
                raise FileExistsError(f"Migration destination already exists: {destination}")
            if not _within(path, target) or not _within(destination, target):
                raise RuntimeError("Refusing to rename outside the requested migration root.")
            renames.append(
                {
                    "source": str(path.relative_to(ROOT)),
                    "destination": str(destination.relative_to(ROOT)),
                    "sha256_before": _sha256(path),
                    "size_bytes": path.stat().st_size,
                }
            )
    return text_changes, renames


def _refresh_hashes(target: Path) -> list[dict]:
    refreshed = []
    for manifest_path in sorted(target.rglob("candidate_receiver_manifest.json")):
        resolved_path = manifest_path.parent / "resolved_run.json"
        if not resolved_path.exists():
            continue
        payload = json.loads(resolved_path.read_text(encoding="utf-8"))
        payload["candidate_manifest_sha256"] = _sha256(manifest_path)
        _atomic_write(resolved_path, json.dumps(payload, indent=2) + "\n")
        digest = _sha256(resolved_path)
        sha_path = resolved_path.parent / "resolved_run.sha256"
        _atomic_write(sha_path, f"{digest}  {resolved_path.name}\n")
        refreshed.append(
            {
                "resolved_run": str(resolved_path.relative_to(ROOT)),
                "resolved_run_sha256": digest,
                "candidate_manifest_sha256": payload["candidate_manifest_sha256"],
            }
        )
    return refreshed


def _validate_outputs(target: Path, renamed: list[dict]) -> dict:
    stale_text = []
    stale_names = []
    for path in sorted(item for item in target.rglob("*") if item.is_file()):
        if path.name == RECEIPT_NAME:
            continue
        if _replace(path.name) != path.name:
            stale_names.append(str(path.relative_to(ROOT)))
        if path.suffix.lower() in TEXT_SUFFIXES:
            content = path.read_text(encoding="utf-8")
            if _replace(content) != content:
                stale_text.append(str(path.relative_to(ROOT)))
    if stale_names or stale_text:
        raise RuntimeError(f"Stale curve IDs remain: names={stale_names}, text={stale_text}")
    verified_binary_files = []
    for item in renamed:
        destination = ROOT / item["destination"]
        digest = _sha256(destination)
        if digest != item["sha256_before"] or destination.stat().st_size != item["size_bytes"]:
            raise RuntimeError(f"Renamed file content changed: {destination}")
        verified_binary_files.append(
            {
                "path": item["destination"],
                "sha256_after": digest,
                "size_bytes": item["size_bytes"],
            }
        )
    interval_files = list(target.rglob("intervals.csv"))
    referenced_arrays = 0
    for interval_path in interval_files:
        with interval_path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                for key in ("error_flags", "ce_nmse_trial_path"):
                    referenced = Path(row[key])
                    referenced = referenced if referenced.is_absolute() else ROOT / referenced
                    if not referenced.exists():
                        raise FileNotFoundError(f"Missing migrated array reference: {referenced}")
                    referenced_arrays += 1
    return {
        "stale_text_count": 0,
        "stale_filename_count": 0,
        "verified_renamed_file_count": len(verified_binary_files),
        "verified_renamed_files": verified_binary_files,
        "interval_file_count": len(interval_files),
        "verified_array_reference_count": referenced_arrays,
    }


def migrate(target: Path, apply: bool) -> dict:
    target = target.resolve()
    if not target.is_dir() or not _within(target, ALLOWED_ROOT):
        raise ValueError(f"Migration root must be inside {ALLOWED_ROOT}")
    text_changes, renames = _inventory(target)
    plan = {
        "schema": "plan033-curve-id-migration-v1",
        "mode": "apply" if apply else "dry-run",
        "root": str(target.relative_to(ROOT)),
        "replacements": REPLACEMENTS,
        "text_changes": text_changes,
        "renames": renames,
    }
    if not apply:
        return plan
    for item in text_changes:
        path = ROOT / item["path"]
        _atomic_write(path, _replace(path.read_text(encoding="utf-8")))
    for item in renames:
        source = ROOT / item["source"]
        destination = ROOT / item["destination"]
        source.replace(destination)
    plan["refreshed_hashes"] = _refresh_hashes(target)
    plan["validation"] = _validate_outputs(target, renames)
    plan["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
    receipt = target / RECEIPT_NAME
    _atomic_write(receipt, json.dumps(plan, indent=2) + "\n")
    plan["receipt"] = str(receipt.relative_to(ROOT))
    return plan


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    print(json.dumps(migrate(args.root, args.apply), indent=2))


if __name__ == "__main__":
    main()
