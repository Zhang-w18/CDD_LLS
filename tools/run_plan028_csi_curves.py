"""Run Plan-028 CSI-curve supplements and the A300 comb-6 comparison."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
import sys
import time
from pathlib import Path
from typing import Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cdd_lls.core.mcs import build_tb_layout, get_mcs
from cdd_lls.phy.channel_tdl import generate_sionna_tdl_channel
from cdd_lls.phy.ldpc import SionnaLDPCAdapter
from cdd_lls.phy.precoding import equivalent_channel
from cdd_lls.phy.qam import qam_demapper_maxlog, qam_modulate
from cdd_lls.phy.resource_grid import local_indices_for_subcarriers
from tools import run_plan027_bler as bler027
from tools import run_plan027_dense_dmrs as dense027
from tools.run_plan025_delay_matched_tdl import (
    read_csv_rows,
    stable_seed,
    wilson,
    write_csv_rows,
)
from tools.run_v_design_piecewise_tradeoff import decode_same_tb_batch


DEFAULT_OUTPUT_ROOT = ROOT / "outputs" / "experiment028_csi_curves"
DEFAULT_RUN_ID = "20260803_main"
SOURCE_RUN = ROOT / "outputs" / "experiment027_meff_sidon" / "20260726_main"
SOURCE_SPECS = {
    "A30": {
        "root": SOURCE_RUN / "e5_dense_dmrs" / "comb6",
        "manifest": Path("manifest/e5_link_manifest.json"),
        "sha": Path("manifest/e5_link_manifest.sha256"),
    },
    "A100": {
        "root": SOURCE_RUN / "e6_a100_dense_dmrs" / "comb6",
        "manifest": Path("manifest/e6_link_manifest.json"),
        "sha": Path("manifest/e6_link_manifest.sha256"),
    },
}
STAGE_PRIORITY = {"prescan": 0, "refine_10pct": 1, "refine_1pct": 2}
FAMILY_ORDER = (
    "B0_QC",
    "AP_RMS_T1",
    "AP_TEPS_T1",
    "AP_TU_NT",
    "AP_TU_NTM1",
    "AP_TALIAS_NT",
    "S0_SIDON",
    "AP_T2_CTRL",
    "GEO_T1_CTRL",
    "MEFF_T2_CAND",
)
PLOT_FONT = 16
TICK_FONT = 14
TITLE_FONT = 18
LINE_WIDTH = 2.5
MARKER_SIZE = 8.0


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repo_relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT.resolve()))


def scenario_root(output: Path, scenario_id: str) -> Path:
    return output / f"{scenario_id.lower()}_comb6"


def a300_dense_root(output: Path) -> Path:
    dense027.configure_scenario("A300", 6)
    return dense027.dense_root(output, 6)


def manifest_paths(output: Path, scenario_id: str) -> tuple[Path, Path, Path]:
    if scenario_id in SOURCE_SPECS:
        root = scenario_root(output, scenario_id) / "manifest"
        return (
            root / "source_manifest.json",
            root / "source_manifest.sha256",
            root / "source_approval.json",
        )
    root = a300_dense_root(output) / "manifest"
    return (
        root / "p28_link_manifest.json",
        root / "p28_link_manifest.sha256",
        root / "p28_link_approval.json",
    )


def load_manifest(output: Path, scenario_id: str) -> tuple[dict, str, Path, Path]:
    manifest_path, sha_path, approval_path = manifest_paths(output, scenario_id)
    if not manifest_path.exists() or not sha_path.exists() or not approval_path.exists():
        raise RuntimeError(f"Manifest bundle is incomplete for {scenario_id}.")
    digest = file_sha256(manifest_path)
    recorded = sha_path.read_text(encoding="utf-8").split()[0]
    if digest != recorded:
        raise RuntimeError(f"Manifest SHA-256 mismatch for {scenario_id}.")
    manifest = bler027.validate_approved_manifest(
        manifest_path, digest, approval_path
    )
    return manifest, digest, manifest_path, approval_path


def _flag_path(stage_dir: Path, snr_db: float, candidate_id: str) -> Path:
    key = f"{float(snr_db):+.2f}".replace("+", "p").replace("-", "m").replace(".", "p")
    safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", candidate_id)
    return stage_dir / "error_flags" / key / f"{safe_id}.npy"


def _canonical_curve_rows(link_root: Path, scenario_id: str) -> list[dict]:
    selected: dict[tuple[str, float], dict] = {}
    for stage, priority in STAGE_PRIORITY.items():
        stage_dir = link_root / stage
        points = stage_dir / "bler_points.csv"
        if not points.exists():
            continue
        for source in read_csv_rows(points):
            if str(source["scenario_id"]) != scenario_id:
                raise RuntimeError(f"Unexpected scenario in {points}.")
            row = dict(source)
            candidate_id = str(row["candidate_id"])
            snr_db = float(row["snr_db"])
            flag = _flag_path(stage_dir, snr_db, candidate_id)
            if not flag.exists() or len(np.load(flag, mmap_mode="r")) != int(row["trials"]):
                raise RuntimeError(f"Missing or invalid error flags: {flag}")
            row["source_stage"] = stage
            row["source_priority"] = priority
            row["source_points_csv"] = repo_relative(points)
            row["source_error_flags"] = repo_relative(flag)
            key = (candidate_id, snr_db)
            previous = selected.get(key)
            rank = (int(row["trials"]), priority)
            old_rank = (
                (int(previous["trials"]), int(previous["source_priority"]))
                if previous is not None
                else (-1, -1)
            )
            if rank > old_rank:
                selected[key] = row
    rows = list(selected.values())
    rows.sort(key=lambda row: (str(row["candidate_id"]), float(row["snr_db"])))
    if not rows:
        raise RuntimeError(f"No estimated-CSI points found for {scenario_id}.")
    return rows


def _candidate_sort_key(row: dict) -> tuple[int, int, str]:
    family = str(row["family"])
    try:
        family_index = FAMILY_ORDER.index(family)
    except ValueError:
        family_index = len(FAMILY_ORDER)
    match = re.search(r"_(\d+)$", str(row["candidate_id"]))
    suffix = int(match.group(1)) if match else -1
    return family_index, suffix, str(row["candidate_id"])


def write_candidate_assets(output: Path, scenario_id: str, manifest: dict) -> None:
    final = scenario_root(output, scenario_id) / "final"
    final.mkdir(parents=True, exist_ok=True)
    candidates = sorted(
        bler027.scenario_candidates(manifest, scenario_id), key=_candidate_sort_key
    )
    rows = []
    styles = {}
    for index, candidate in enumerate(candidates, start=1):
        candidate_id = str(candidate["candidate_id"])
        family = str(candidate["family"])
        style = bler027.curve_style(candidate_id, family)
        styles[candidate_id] = {
            "color": style["color"],
            "linestyle": style["linestyle"],
            "marker": style["marker"],
        }
        rows.append(
            {
                "scenario_id": scenario_id,
                "candidate_number": f"candidate {index:02d}",
                "candidate_id": candidate_id,
                "family": family,
                "represented_families": json.dumps(
                    candidate.get("represented_families", [family]), separators=(",", ":")
                ),
                "roles": json.dumps(candidate.get("roles", []), separators=(",", ":")),
                "delay_grid_coordinates": json.dumps(
                    candidate["delay_grid_coordinates"], separators=(",", ":")
                ),
                "delay_ns": json.dumps(candidate["delay_ns"], separators=(",", ":")),
                "pilot_rank": int(candidate["pilot_rank"]),
                "pilot_condition_number": float(candidate["pilot_condition_number"]),
            }
        )
    write_csv_rows(rows, final / "candidate_number_delay_table.csv")
    (final / "curve_styles.json").write_text(
        json.dumps(styles, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def freeze_curve_grid(output: Path, scenario_id: str) -> list[dict]:
    manifest, digest, manifest_path, _ = load_manifest(output, scenario_id)
    if scenario_id in SOURCE_SPECS:
        link_root = SOURCE_SPECS[scenario_id]["root"] / "link" / scenario_id
    else:
        link_root = a300_dense_root(output) / "link" / scenario_id
    rows = _canonical_curve_rows(link_root, scenario_id)
    expected = {row["candidate_id"] for row in manifest["candidates"]}
    actual = {row["candidate_id"] for row in rows}
    if actual != expected:
        raise RuntimeError(
            f"Canonical curve candidates differ from manifest for {scenario_id}: "
            f"missing={sorted(expected-actual)}, extra={sorted(actual-expected)}"
        )
    destination = scenario_root(output, scenario_id) / "final"
    destination.mkdir(parents=True, exist_ok=True)
    write_csv_rows(rows, destination / "estimated_csi_curve_points.csv")
    schedule = {
        "scenario_id": scenario_id,
        "dmrs_spacing_subcarriers": 6,
        "manifest_sha256": digest,
        "manifest_path": repo_relative(manifest_path),
        "selection_rule": "max trials; ties refine_1pct > refine_10pct > prescan",
        "points": [
            {
                "candidate_id": row["candidate_id"],
                "snr_db": float(row["snr_db"]),
                "trials": int(row["trials"]),
                "source_stage": row["source_stage"],
            }
            for row in rows
        ],
    }
    (destination / "curve_schedule.json").write_text(
        json.dumps(schedule, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    write_candidate_assets(output, scenario_id, manifest)
    return rows


def import_existing(output: Path, scenario_id: str) -> None:
    if scenario_id not in SOURCE_SPECS:
        raise ValueError("Only A30 and A100 are imported from result-027.")
    spec = SOURCE_SPECS[scenario_id]
    source_manifest = spec["root"] / spec["manifest"]
    source_sha = spec["root"] / spec["sha"]
    source_approval = source_manifest.with_name(source_manifest.stem.replace("manifest", "approval") + ".json")
    recorded = source_sha.read_text(encoding="utf-8").split()[0]
    if file_sha256(source_manifest) != recorded:
        raise RuntimeError(f"Source manifest hash mismatch for {scenario_id}.")
    destination = scenario_root(output, scenario_id) / "manifest"
    destination.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_manifest, destination / "source_manifest.json")
    (destination / "source_manifest.sha256").write_text(
        f"{recorded}  source_manifest.json\n", encoding="utf-8"
    )
    shutil.copy2(source_approval, destination / "source_approval.json")
    (destination / "source_provenance.json").write_text(
        json.dumps(
            {
                "source_manifest": repo_relative(source_manifest),
                "source_sha256": recorded,
                "estimated_csi_reused_without_rerun": True,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    freeze_curve_grid(output, scenario_id)


def run_ideal(output: Path, scenario_id: str, batch_size: int) -> None:
    manifest, digest, _, _ = load_manifest(output, scenario_id)
    canonical_path = scenario_root(output, scenario_id) / "final" / "estimated_csi_curve_points.csv"
    if not canonical_path.exists():
        freeze_curve_grid(output, scenario_id)
    schedule_rows = [dict(row) for row in read_csv_rows(canonical_path)]
    candidates = bler027.scenario_candidates(manifest, scenario_id)
    by_id = {str(row["candidate_id"]): row for row in candidates}
    channel, grid, precoders, _ = bler027.build_scene(manifest, candidates, scenario_id)
    mcs = get_mcs("nr_256qam", 8, None, None)
    tb = build_tb_layout(grid.n_data_re, mcs)
    adapter = SionnaLDPCAdapter(
        tb.cb_k_values, tb.cb_e_values, num_iter=8, llr_clip=50.0
    )
    data_local = local_indices_for_subcarriers(grid, grid.data_subcarrier_indices)
    stage = scenario_root(output, scenario_id) / "ideal_link"
    stage.mkdir(parents=True, exist_ok=True)
    result_path = stage / "ideal_csi_bler_points.csv"
    existing = [dict(row) for row in read_csv_rows(result_path)] if result_path.exists() else []
    complete = {
        (str(row["candidate_id"]), float(row["snr_db"]), int(row["trials"]))
        for row in existing
        if _flag_path(stage, float(row["snr_db"]), str(row["candidate_id"])).exists()
    }
    grouped: dict[float, list[dict]] = {}
    for row in schedule_rows:
        grouped.setdefault(float(row["snr_db"]), []).append(row)
    for snr_db in sorted(grouped):
        pending = [
            row
            for row in grouped[snr_db]
            if (str(row["candidate_id"]), snr_db, int(row["trials"])) not in complete
        ]
        if not pending:
            continue
        maximum_trials = max(int(row["trials"]) for row in pending)
        snr_linear = 10.0 ** (snr_db / 10.0)
        noise_variance = 8.0 / snr_linear
        flags = {
            str(row["candidate_id"]): np.zeros(int(row["trials"]), dtype=bool)
            for row in pending
        }
        started = time.time()
        for batch_start in range(1, maximum_trials + 1, int(batch_size)):
            current = min(int(batch_size), maximum_trials - batch_start + 1)
            payload_rng = np.random.default_rng(
                stable_seed(bler027.SEED, scenario_id, snr_db, batch_start, "payload")
            )
            payload = [
                payload_rng.integers(0, 2, size=int(k), dtype=np.int8)
                for k in tb.cb_k_values
            ]
            coded = np.concatenate(adapter.encode(payload))
            symbols = qam_modulate(coded, int(mcs.qm))
            realization = generate_sionna_tdl_channel(
                grid,
                channel,
                n_tx=8,
                n_rx=1,
                batch_size=current,
                seed=stable_seed(
                    bler027.SEED, scenario_id, snr_db, batch_start, "channel"
                ),
            )
            noise_rng = np.random.default_rng(
                stable_seed(bler027.SEED, scenario_id, snr_db, batch_start, "noise")
            )
            # Preserve the Plan-027 RNG stream exactly: estimated-CSI consumes
            # the averaged pilot-noise draw before drawing data noise.
            _ = noise_rng.normal(size=(current, 1, grid.pilot_count))
            _ = noise_rng.normal(size=(current, 1, grid.pilot_count))
            data_noise = math.sqrt(noise_variance / 2.0) * (
                noise_rng.normal(size=(current, 1, grid.n_data_re))
                + 1j * noise_rng.normal(size=(current, 1, grid.n_data_re))
            )
            for row in pending:
                candidate_id = str(row["candidate_id"])
                trials = int(row["trials"])
                if batch_start > trials:
                    continue
                used = min(current, trials - batch_start + 1)
                effective = equivalent_channel(
                    realization.H, precoders[candidate_id].C
                )[:used, 0:1]
                true_data = effective[:, :, grid.data_symbol_indices, data_local]
                received = true_data * symbols[None, None, :] + data_noise[:used]
                equalized, effective_noise = ideal_csi_equalize(
                    true_data, received, noise_variance
                )
                llrs = [
                    qam_demapper_maxlog(
                        equalized[index], effective_noise[index], int(mcs.qm)
                    )
                    for index in range(used)
                ]
                decoded = decode_same_tb_batch(adapter, llrs, payload)
                first = batch_start - 1
                flags[candidate_id][first : first + used] = np.asarray(
                    [not item.tb_success for item in decoded], dtype=bool
                )
            completed_trials = batch_start + current - 1
            if (
                batch_start == 1
                or completed_trials % 500 < current
                or completed_trials == maximum_trials
            ):
                print(
                    f"[P28 ideal {scenario_id}] snr={snr_db:g} "
                    f"trial={completed_trials}/{maximum_trials} "
                    f"candidates={len(pending)} elapsed={time.time()-started:.1f}s",
                    flush=True,
                )
        keys = {(str(row["candidate_id"]), snr_db) for row in pending}
        existing = [
            row
            for row in existing
            if (str(row["candidate_id"]), float(row["snr_db"])) not in keys
        ]
        for row in pending:
            candidate_id = str(row["candidate_id"])
            trials = int(row["trials"])
            errors = int(np.sum(flags[candidate_id]))
            low, high = wilson(errors, trials)
            flag_path = _flag_path(stage, snr_db, candidate_id)
            flag_path.parent.mkdir(parents=True, exist_ok=True)
            np.save(flag_path, flags[candidate_id])
            existing.append(
                {
                    "scenario_id": scenario_id,
                    "candidate_id": candidate_id,
                    "family": by_id[candidate_id]["family"],
                    "snr_db": snr_db,
                    "trials": trials,
                    "tb_errors": errors,
                    "bler": errors / trials,
                    "bler_wilson95_lo": low,
                    "bler_wilson95_hi": high,
                    "noise_variance": noise_variance,
                    "receiver": "true data-RE effective CSI equalization and LDPC decoding",
                }
            )
        existing.sort(key=lambda row: (str(row["candidate_id"]), float(row["snr_db"])))
        write_csv_rows(existing, result_path)
    (stage / "resolved_experiment.json").write_text(
        json.dumps(
            {
                "scenario_id": scenario_id,
                "manifest_sha256": digest,
                "seed": bler027.SEED,
                "receiver": "true data-RE effective CSI",
                "real_link_components": [
                    "TDL-A channel",
                    "payload",
                    "LDPC encode",
                    "256QAM",
                    "AWGN",
                    "true-CSI equalization",
                    "soft demap",
                    "LDPC decode",
                ],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def _style_for_plot(style: dict) -> dict:
    linestyle = style["linestyle"]
    if isinstance(linestyle, list):
        linestyle = tuple(
            tuple(value) if isinstance(value, list) else value for value in linestyle
        )
    return {"color": style["color"], "linestyle": linestyle, "marker": style["marker"]}


def ideal_csi_equalize(
    true_data: np.ndarray, received: np.ndarray, noise_variance: float
) -> tuple[np.ndarray, np.ndarray]:
    """Equalize data REs with the exact effective channel used to form ``received``."""
    denominator = np.maximum(np.sum(np.abs(true_data) ** 2, axis=1), 1e-10)
    equalized = np.sum(np.conj(true_data) * received, axis=1) / denominator
    return equalized, float(noise_variance) / denominator


def _plot_metric(
    rows: Sequence[dict],
    mapping: dict[str, str],
    styles: dict,
    metric: str,
    ylabel: str,
    title: str,
    path: Path,
    log_y: bool,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(16.0, 9.0))
    for candidate_id, label in mapping.items():
        selected = sorted(
            (row for row in rows if str(row["candidate_id"]) == candidate_id),
            key=lambda row: float(row["snr_db"]),
        )
        if not selected:
            continue
        x = [float(row["snr_db"]) for row in selected]
        if metric == "bler":
            y = [
                max(float(row[metric]), 0.5 / int(row["trials"]))
                for row in selected
            ]
        else:
            y = [float(row[metric]) for row in selected]
        axis.plot(
            x,
            y,
            label=label,
            linewidth=LINE_WIDTH,
            markersize=MARKER_SIZE,
            **_style_for_plot(styles[candidate_id]),
        )
    if log_y:
        axis.set_yscale("log")
    axis.set_xlabel("SNR (dB)", fontsize=PLOT_FONT)
    axis.set_ylabel(ylabel, fontsize=PLOT_FONT)
    axis.set_title(title, fontsize=TITLE_FONT)
    axis.tick_params(axis="both", which="both", labelsize=TICK_FONT)
    axis.grid(True, which="both", alpha=0.3)
    axis.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=2,
        fontsize=PLOT_FONT,
        frameon=True,
    )
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def plot_curves(output: Path, scenario_id: str) -> None:
    final = scenario_root(output, scenario_id) / "final"
    estimated = [dict(row) for row in read_csv_rows(final / "estimated_csi_curve_points.csv")]
    ideal = [
        dict(row)
        for row in read_csv_rows(
            scenario_root(output, scenario_id) / "ideal_link" / "ideal_csi_bler_points.csv"
        )
    ]
    number_rows = [dict(row) for row in read_csv_rows(final / "candidate_number_delay_table.csv")]
    styles = json.loads((final / "curve_styles.json").read_text(encoding="utf-8"))
    estimated_grid = {
        (str(row["candidate_id"]), float(row["snr_db"]), int(row["trials"]))
        for row in estimated
    }
    ideal_grid = {
        (str(row["candidate_id"]), float(row["snr_db"]), int(row["trials"]))
        for row in ideal
    }
    if estimated_grid != ideal_grid:
        raise RuntimeError(f"Ideal/estimated curve grids differ for {scenario_id}.")
    numbered = {
        str(row["candidate_id"]): str(row["candidate_number"]) for row in number_rows
    }
    abbreviated = {
        str(row["candidate_id"]): str(row["candidate_id"]).removeprefix(f"{scenario_id}_")
        for row in number_rows
    }
    label_modes = {"candidate": numbered}
    if scenario_id == "A300":
        label_modes = {"abbrev": abbreviated, "candidate": numbered}
    for mode, mapping in label_modes.items():
        prefix = f"{scenario_id.lower()}_comb6_{mode}"
        _plot_metric(
            estimated,
            mapping,
            styles,
            "bler",
            "Estimated-CSI BLER",
            f"{scenario_id} comb-6 estimated-CSI BLER",
            final / f"{prefix}_estimated_csi_bler.png",
            True,
        )
        _plot_metric(
            ideal,
            mapping,
            styles,
            "bler",
            "Ideal-CSI BLER",
            f"{scenario_id} comb-6 ideal-CSI BLER (decoded link)",
            final / f"{prefix}_ideal_csi_bler.png",
            True,
        )
        _plot_metric(
            estimated,
            mapping,
            styles,
            "ce_nmse_mean_db",
            "Matched CE NMSE (dB)",
            f"{scenario_id} comb-6 matched CE NMSE",
            final / f"{prefix}_ce_nmse.png",
            False,
        )


def run_a300_estimated(output: Path, stage: str, batch_size: int) -> None:
    manifest, digest, _, _ = load_manifest(output, "A300")
    scene_root = a300_dense_root(output) / "link" / "A300"
    if stage == "estimated-smoke":
        bler027.run_smoke(manifest, digest, "A300", batch_size, scene_root)
    elif stage == "estimated-prescan":
        bler027.run_prescan(manifest, digest, "A300", batch_size, scene_root)
    elif stage == "estimated-refine-10pct":
        bler027.run_refine(manifest, digest, "A300", "10pct", batch_size, scene_root)
    elif stage == "estimated-refine-1pct":
        bler027.run_refine(manifest, digest, "A300", "1pct", batch_size, scene_root)
    else:
        raise ValueError(stage)


def run_analyze(output: Path, scenario_id: str) -> None:
    manifest, digest, _, _ = load_manifest(output, scenario_id)
    if scenario_id == "A300":
        scene_root = a300_dense_root(output) / "link" / "A300"
        bler027.analyze_scene(manifest, digest, "A300", scene_root)
        dense027.run_analyze(output, 6)
    freeze_curve_grid(output, scenario_id)
    plot_curves(output, scenario_id)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage",
        choices=(
            "import-existing",
            "search",
            "outage",
            "manifest",
            "estimated-smoke",
            "estimated-prescan",
            "estimated-refine-10pct",
            "estimated-refine-1pct",
            "ideal",
            "analyze",
        ),
        required=True,
    )
    parser.add_argument("--scenario", choices=("A30", "A100", "A300"), required=True)
    parser.add_argument("--run-id", default=DEFAULT_RUN_ID)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--search-states", type=int, default=200_000)
    parser.add_argument("--geometry-states", type=int, default=200_000)
    parser.add_argument("--outage-samples", type=int, default=200_000)
    parser.add_argument("--bootstrap-repeats", type=int, default=1_000)
    parser.add_argument("--batch-size", type=int, default=100)
    args = parser.parse_args()
    output = args.output_root / str(args.run_id)
    if args.stage == "import-existing":
        import_existing(output, str(args.scenario))
        return
    if args.stage in ("search", "outage", "manifest"):
        if args.scenario != "A300":
            raise ValueError("Search/outage/manifest are only run for A300.")
        dense027.configure_scenario("A300", 6)
        if args.stage == "search":
            dense027.run_search(
                output, 6, int(args.search_states), int(args.geometry_states)
            )
        elif args.stage == "outage":
            dense027.run_outage(
                output, 6, int(args.outage_samples), int(args.bootstrap_repeats)
            )
        else:
            dense027.run_manifest(output, 6)
        return
    if args.stage.startswith("estimated-"):
        if args.scenario != "A300":
            raise ValueError("Existing A30/A100 estimated-CSI links are not rerun.")
        run_a300_estimated(output, str(args.stage), int(args.batch_size))
    elif args.stage == "ideal":
        run_ideal(output, str(args.scenario), int(args.batch_size))
    else:
        run_analyze(output, str(args.scenario))


if __name__ == "__main__":
    main()
