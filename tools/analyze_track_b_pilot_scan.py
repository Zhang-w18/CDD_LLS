"""plan-023 analysis: H1 dominance / H2 collapse / H3 pairwise / QC gap + figures.

Reads the scan outputs produced by tools/run_track_b_pilot_scan.py and writes
analysis CSVs plus all figures required by plan-023 §7.2.
"""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import sys
from typing import Dict, List, Optional

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/cdd_lls_matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tools.run_track_b_pilot_scan import SNR_GRID_DB, write_csv

FAM_COLOR = {
    "B1_cdd": "#4878a8",
    "B2_prg_cycling": "#e69f00",
    "B3_cc": "#d62728",
    "B4_cn": "#9467bd",
    "B5_qc": "#111111",
    "B6_appendix": "#2a9d8f",
}
FAM_LABEL = {
    "B1_cdd": "B1 transparent CDD",
    "B2_prg_cycling": "B2 PRG cycling",
    "B3_cc": "B3 CC (continuous cycling)",
    "B4_cn": "B4 CN (capped slope)",
    "B5_qc": "B5 QC benchmark",
}


def fnum(x: object) -> float:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float("nan")
    return v


def read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


# ----------------------------------------------------------------------------
# H1: dominance vs baseline (B1 ∪ B2) Pareto front on (nmse_R1_16dB, outage10)
# ----------------------------------------------------------------------------
def h1_analysis(rows: List[Dict[str, str]]) -> Dict[str, object]:
    base = [r for r in rows if r["family"] in ("B1_cdd", "B2_prg_cycling")
            and np.isfinite(fnum(r["nmse_R1_16dB"])) and np.isfinite(fnum(r["outage10_dB"]))]
    challengers = [r for r in rows if r["family"] in ("B3_cc", "B4_cn")]

    def base_front() -> List[Dict[str, str]]:
        front = []
        for r in base:
            x, y = fnum(r["nmse_R1_16dB"]), fnum(r["outage10_dB"])
            dominated = any(
                fnum(o["nmse_R1_16dB"]) <= x and fnum(o["outage10_dB"]) <= y
                and (fnum(o["nmse_R1_16dB"]) < x or fnum(o["outage10_dB"]) < y)
                for o in base if o is not r)
            if not dominated:
                front.append(r)
        return sorted(front, key=lambda r: fnum(r["nmse_R1_16dB"]))

    front = base_front()
    out_rows = []
    for c in challengers:
        cn, co = fnum(c["nmse_R1_16dB"]), fnum(c["outage10_dB"])
        if not (np.isfinite(cn) and np.isfinite(co)):
            continue
        # (i) baselines at least as good in outage (within +0.05 dB): NMSE advantage
        s1 = [fnum(b["nmse_R1_16dB"]) for b in base if fnum(b["outage10_dB"]) <= co + 0.05]
        margin_i = (min(s1) - cn) if s1 else float("inf")
        # (ii) baselines at least as good in NMSE (within +0.2 dB): outage advantage
        s2 = [fnum(b["outage10_dB"]) for b in base if fnum(b["nmse_R1_16dB"]) <= cn + 0.2]
        margin_ii = (min(s2) - co) if s2 else float("inf")
        out_rows.append({
            "id": c["id"], "family": c["family"],
            "nmse_R1_16dB": cn, "outage10_dB": co, "outage1_dB": fnum(c["outage1_dB"]),
            "margin_i_nmse_dB": round(margin_i, 4) if np.isfinite(margin_i) else "inf",
            "margin_ii_outage_dB": round(margin_ii, 4) if np.isfinite(margin_ii) else "inf",
            "pass_i": int(margin_i >= 1.0),
            "pass_ii": int(margin_ii >= 0.3),
            "pass_h1": int(margin_i >= 1.0 or margin_ii >= 0.3),
        })
    out_rows.sort(key=lambda r: -(1e9 if r["margin_i_nmse_dB"] == "inf" else fnum(r["margin_i_nmse_dB"])))
    return {"front": front, "challengers": out_rows}


# ----------------------------------------------------------------------------
# H2: collapse — bucket by tau_span (±10%), range of R2 NMSE (tau_w = 1389 ns)
# ----------------------------------------------------------------------------
def h2_analysis(rows: List[Dict[str, str]]) -> List[Dict[str, object]]:
    pts = [r for r in rows if np.isfinite(fnum(r.get("nmse_R2_tw1389_16dB")))
           and np.isfinite(fnum(r["tau_span_ns"]))]
    pts.sort(key=lambda r: fnum(r["tau_span_ns"]))
    buckets: List[List[Dict[str, str]]] = []
    cur: List[Dict[str, str]] = []
    lo = None
    for r in pts:
        t = fnum(r["tau_span_ns"])
        if not cur:
            cur = [r]
            lo = t
            continue
        # ±10% band around center c = lo/0.9 -> upper edge = lo * 11/9
        if t <= max(lo * 11.0 / 9.0, lo + 1e-9):
            cur.append(r)
        else:
            buckets.append(cur)
            cur = [r]
            lo = t
    if cur:
        buckets.append(cur)
    out = []
    for b in buckets:
        vals = [fnum(r["nmse_R2_tw1389_16dB"]) for r in b]
        taus = [fnum(r["tau_span_ns"]) for r in b]
        fams = sorted({r["family"] for r in b})
        vmin, vmax = min(vals), max(vals)
        out.append({
            "tau_lo_ns": round(min(taus), 1), "tau_hi_ns": round(max(taus), 1),
            "n": len(b), "families": "+".join(f.split("_")[0] for f in fams),
            "nmse_min_dB": round(vmin, 3), "nmse_max_dB": round(vmax, 3),
            "range_dB": round(vmax - vmin, 3),
            "collapse_ok": int(vmax - vmin < 0.5) if len(b) > 1 else "",
            "id_min": min(b, key=lambda r: fnum(r["nmse_R2_tw1389_16dB"]))["id"],
            "id_max": max(b, key=lambda r: fnum(r["nmse_R2_tw1389_16dB"]))["id"],
        })
    return out


# ----------------------------------------------------------------------------
# H3: pairwise outage among B6 (48 PRB only)
# ----------------------------------------------------------------------------
def h3_analysis(rows: List[Dict[str, str]]) -> List[Dict[str, object]]:
    b6 = {r["id"]: r for r in rows if r["family"] == "B6_appendix"}
    out = []
    arith = {int(r["sigma"]): r for r in b6.values() if r.get("sigma") not in ("", None)}
    sig = sorted(arith)
    for i, s1 in enumerate(sig):
        for s2 in sig[i + 1:]:
            d10 = fnum(arith[s2]["outage10_dB"]) - fnum(arith[s1]["outage10_dB"])
            d1 = fnum(arith[s2]["outage1_dB"]) - fnum(arith[s1]["outage1_dB"])
            out.append({
                "pair": f"sigma{s2} - sigma{s1}", "d_outage10_dB": round(d10, 4),
                "d_outage1_dB": round(d1, 4),
                "same_orbit_ok": int(abs(d10) <= 0.05 and abs(d1) <= 0.05),
            })
    if "B6_sidon" in b6 and arith:
        best1 = min(fnum(r["outage1_dB"]) for r in arith.values())
        best10 = min(fnum(r["outage10_dB"]) for r in arith.values())
        out.append({
            "pair": "sidon - best_arith",
            "d_outage10_dB": round(fnum(b6["B6_sidon"]["outage10_dB"]) - best10, 4),
            "d_outage1_dB": round(fnum(b6["B6_sidon"]["outage1_dB"]) - best1, 4),
            "same_orbit_ok": "",
        })
    return out


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
def fig_pareto(rows: List[Dict[str, str]], front: List[Dict[str, str]],
               outage_key: str, n_prbs: int, path: Path, top_ids: List[str]) -> None:
    fig, ax = plt.subplots(figsize=(9.2, 6.6))
    for fam in ("B4_cn", "B1_cdd", "B2_prg_cycling", "B3_cc"):
        xs = [fnum(r["nmse_R1_16dB"]) for r in rows if r["family"] == fam]
        ys = [fnum(r[outage_key]) for r in rows if r["family"] == fam]
        size = 14 if fam == "B4_cn" else 42
        marker = {"B1_cdd": "o", "B2_prg_cycling": "s", "B3_cc": "^", "B4_cn": "."}[fam]
        ax.scatter(xs, ys, s=size, c=FAM_COLOR[fam], marker=marker, label=FAM_LABEL[fam],
                   alpha=0.75, zorder=3)
    # staircase front (B1 ∪ B2), sorted by NMSE; steps go down-right
    fx = [fnum(r["nmse_R1_16dB"]) for r in front]
    fy = [fnum(r[outage_key]) for r in front]
    order = np.argsort(fx)
    fx = [fx[i] for i in order]
    fy = [fy[i] for i in order]
    sx, sy = [], []
    for i in range(len(fx)):
        if i > 0:
            sx.append(fx[i])
            sy.append(fy[i - 1])
        sx.append(fx[i])
        sy.append(fy[i])
    ax.plot(sx, sy, "-", c="#666666", lw=1.4, label="baseline front (B1∪B2)", zorder=2)
    qc = [r for r in rows if r["family"] == "B5_qc"]
    if qc and np.isfinite(fnum(qc[0].get("nmse_R3_16dB"))):
        ax.scatter([fnum(qc[0]["nmse_R3_16dB"])], [fnum(qc[0][outage_key])],
                   marker="*", s=380, c="#111111", label="B5 QC benchmark (R3 matched)", zorder=5)
    for r in rows:
        if r["id"] in top_ids:
            ax.annotate(r["id"].replace("B3_cc_", "").replace("B4_cn_", "cn "),
                        (fnum(r["nmse_R1_16dB"]), fnum(r[outage_key])),
                        textcoords="offset points", xytext=(6, 5), fontsize=7.5,
                        color=FAM_COLOR[r["family"]])
    level = "10%" if outage_key == "outage10_dB" else "1%"
    ax.set_xlabel("R1 sliding-window mismatched NMSE @ snr̄=16 dB (dB)")
    ax.set_ylabel(f"{level} outage SNR (dB)")
    ax.set_title(f"{n_prbs} PRB — CE cost vs diversity ({level} outage), transparent R1 receiver")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_collapse(rows_by_prb: Dict[int, List[Dict[str, str]]], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4), sharey=True)
    for ax, (n_prbs, rows) in zip(axes, sorted(rows_by_prb.items(), reverse=True)):
        for fam in ("B1_cdd", "B2_prg_cycling", "B3_cc", "B4_cn", "B5_qc"):
            pts = [r for r in rows if r["family"] == fam
                   and np.isfinite(fnum(r.get("nmse_R2_tw1389_16dB")))]
            ax.scatter([fnum(r["tau_span_ns"]) for r in pts],
                       [fnum(r["nmse_R2_tw1389_16dB"]) for r in pts],
                       s=22, c=FAM_COLOR[fam], label=FAM_LABEL[fam], alpha=0.8)
        ax.set_xlabel("tau_span (95% energy, ns)")
        ax.set_title(f"{n_prbs} PRB")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("R2 wideband NMSE @16 dB, tau_w=1389 ns (dB)")
    axes[0].legend(fontsize=8)
    fig.suptitle("H2 collapse check: does CE cost depend only on composite delay span?")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_b6(curve_rows: List[Dict[str, str]], path: Path) -> None:
    b6 = [r for r in curve_rows if r["family"] == "B6_appendix"]
    if not b6:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2))
    ref = None
    styles = {"B6_arith_s1": "-", "B6_arith_s2": "--", "B6_arith_s4": "-.",
              "B6_arith_s9": ":", "B6_sidon": "-"}
    colors = {"B6_arith_s1": "#26547c", "B6_arith_s2": "#ef476f", "B6_arith_s4": "#06a77d",
              "B6_arith_s9": "#f77f00", "B6_sidon": "#111111"}
    curves = {}
    for r in b6:
        p = np.array([fnum(r[f"p{s:g}"]) for s in SNR_GRID_DB])
        curves[r["id"]] = p
        axes[0].semilogy(SNR_GRID_DB, np.maximum(p, 1e-6), styles[r["id"]],
                         c=colors[r["id"]], label=r["id"], lw=1.6)
    axes[0].axhline(0.10, c="k", ls=":", lw=0.8)
    axes[0].axhline(0.01, c="k", ls=":", lw=0.8)
    axes[0].set_xlabel("snr̄ (dB)")
    axes[0].set_ylabel("P_out")
    axes[0].set_ylim(1e-4, 1.0)
    axes[0].set_xlim(4, 18)
    axes[0].set_title("B6 outage curves (48 PRB, common h samples)")
    axes[0].grid(alpha=0.3, which="both")
    axes[0].legend(fontsize=8)

    # Delta SNR at each outage level, relative to sigma=1
    def snr_at(p: np.ndarray, target: float) -> float:
        p = np.maximum(p, 1e-9)
        for a in range(len(p) - 1):
            if p[a] >= target > p[a + 1]:
                l0, l1 = math.log10(p[a]), math.log10(p[a + 1])
                return float(SNR_GRID_DB[a] + (l0 - math.log10(target)) / (l0 - l1) * 0.5)
        return float("nan")

    levels = np.logspace(math.log10(0.005), math.log10(0.5), 40)
    ref_curve = curves["B6_arith_s1"]
    for cid, p in curves.items():
        if cid == "B6_arith_s1":
            continue
        d = [snr_at(p, lv) - snr_at(ref_curve, lv) for lv in levels]
        axes[1].semilogx(levels, d, styles[cid], c=colors[cid], label=f"{cid} − arith_s1", lw=1.6)
    axes[1].axhline(0.0, c="k", lw=0.8)
    axes[1].axhline(0.05, c="k", ls=":", lw=0.8)
    axes[1].axhline(-0.05, c="k", ls=":", lw=0.8)
    axes[1].set_xlabel("outage level P_out")
    axes[1].set_ylabel("Δ outage SNR vs arith σ=1 (dB)")
    axes[1].set_title("pairwise outage-SNR differences (±0.05 dB gate dotted)")
    axes[1].grid(alpha=0.3, which="both")
    axes[1].legend(fontsize=8)
    axes[1].invert_xaxis()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_calib0a(calib_dir: Path, path: Path) -> None:
    rows = read_rows(calib_dir / "calib0a_nmse.csv")
    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    idx = np.arange(len(rows))
    cf = [fnum(r["nmse_closed_form_dB"]) for r in rows]
    mc = [fnum(r["nmse_mc_dB"]) for r in rows]
    ax.bar(idx - 0.18, cf, width=0.36, label="closed-form R3", color="#26547c")
    ax.bar(idx + 0.18, mc, width=0.36, label="MC (Alg1 estimator, 300 trials)", color="#ef476f")
    ax.set_xticks(idx)
    ax.set_xticklabels([f'{r["id"].replace("B3_cc_nseg12_T12_seq", "CC(12,12,seq)").replace("B5_qc", "QC")}\n{r["snr_db"]} dB'
                        for r in rows], fontsize=8)
    for i, r in enumerate(rows):
        ax.text(i, max(cf[i], mc[i]) + 0.5, f'Δ={fnum(r["diff_dB"]):+.2f}', ha="center", fontsize=8)
    ax.set_ylabel("matched NMSE (dB)")
    ax.set_title("Phase 0a: closed-form mismatched-MSE vs MC estimator (gate ±0.3 dB)")
    ax.legend()
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_calib0b(calib_dir: Path, path: Path) -> None:
    curves = read_rows(calib_dir / "calib0b_outage_curves.csv")
    bler = read_rows(calib_dir / "calib0b_bler_curves.csv")
    gaps = {r["id"]: r for r in read_rows(calib_dir / "calib0b_gaps.csv")}
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), sharey=True)
    for ax, crow in zip(axes, curves):
        cid = crow["id"]
        p = np.array([fnum(crow[f"p{s:g}"]) for s in SNR_GRID_DB])
        ax.semilogy(SNR_GRID_DB, np.maximum(p, 1e-6), "-", c="#26547c", lw=1.8,
                    label="MI outage P_out")
        pts = sorted([(fnum(r["snr_db"]), fnum(r["bler"]), fnum(r["trials"])) for r in bler
                      if r["id"] == cid])
        ax.semilogy([p0 for p0, _, _ in pts],
                    [max(b, 0.5 / t) for _, b, t in pts], "o-", c="#ef476f", lw=1.8,
                    ms=4, label="ideal-CSI BLER")
        ax.axhline(0.10, c="k", ls=":", lw=0.8)
        ax.axhline(0.01, c="k", ls=":", lw=0.8)
        g = gaps.get(cid, {})
        ax.set_title(f'{cid}\ngap@10%={fnum(g.get("gap10_dB")):.2f} dB, gap@1%={fnum(g.get("gap1_dB")):.2f} dB, '
                     f'|Δ|={fnum(g.get("gap_diff_dB")):.2f} dB')
        ax.set_xlabel("snr̄ (dB)")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=9)
        ax.set_xlim(4, 22)
        ax.set_ylim(1e-3, 1.0)
    axes[0].set_ylabel("probability")
    fig.suptitle("Phase 0b: outage vs ideal-CSI BLER parallelism (gate: |gap10 − gap1| < 0.4 dB)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------------
def main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/track_b_pilot_scan/20260709_main")
    args = parser.parse_args()
    out_dir = Path(args.out)
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    summary: Dict[str, object] = {}
    rows_by_prb: Dict[int, List[Dict[str, str]]] = {}
    for n_prbs in (48, 24):
        path = out_dir / f"candidates_{n_prbs}prb.csv"
        if not path.exists():
            print(f"[analyze] missing {path}, skipping {n_prbs} PRB")
            continue
        rows = read_rows(path)
        rows_by_prb[n_prbs] = rows

        h1 = h1_analysis(rows)
        write_csv(h1["challengers"], out_dir / f"analysis_h1_{n_prbs}prb.csv")
        front_rows = [{k: r.get(k, "") for k in ("id", "family", "nmse_R1_16dB", "outage10_dB", "outage1_dB")}
                      for r in h1["front"]]
        write_csv(front_rows, out_dir / f"baseline_front_{n_prbs}prb.csv")
        n_pass = sum(r["pass_h1"] for r in h1["challengers"])
        top = h1["challengers"][:5]
        summary[f"h1_{n_prbs}prb"] = {
            "n_challengers": len(h1["challengers"]),
            "n_pass": n_pass,
            "top5": [r["id"] for r in top],
        }

        # QC gap (criterion-2 preview)
        qc = [r for r in rows if r["family"] == "B5_qc"]
        gap_rows = []
        if qc:
            qn = fnum(qc[0].get("nmse_R3_16dB"))
            qo = fnum(qc[0]["outage10_dB"])
            for r in top:
                gap_rows.append({
                    "id": r["id"],
                    "d_nmse_vs_qc_dB": round(fnum(r["nmse_R1_16dB"]) - qn, 3),
                    "d_outage10_vs_qc_dB": round(fnum(r["outage10_dB"]) - qo, 3),
                    "near_qc_ref": int(fnum(r["nmse_R1_16dB"]) - qn <= 2.0
                                       and fnum(r["outage10_dB"]) - qo <= 0.5),
                })
            write_csv(gap_rows, out_dir / f"analysis_qc_gap_{n_prbs}prb.csv")
            summary[f"qc_{n_prbs}prb"] = {"nmse_R3_16dB": qn, "outage10_dB": qo}

        h2 = h2_analysis(rows)
        write_csv(h2, out_dir / f"analysis_h2_buckets_{n_prbs}prb.csv")
        multi = [b for b in h2 if b["n"] > 1]
        summary[f"h2_{n_prbs}prb"] = {
            "n_buckets_multi": len(multi),
            "max_range_dB": max((b["range_dB"] for b in multi), default=float("nan")),
            "n_buckets_violating": sum(1 for b in multi if b["range_dB"] >= 0.5),
        }

        top_ids = [r["id"] for r in top]
        fig_pareto(rows, h1["front"], "outage10_dB", n_prbs,
                   fig_dir / f"pareto_{n_prbs}prb_out10.png", top_ids)
        fig_pareto(rows, h1["front"], "outage1_dB", n_prbs,
                   fig_dir / f"pareto_{n_prbs}prb_out1.png", top_ids)

        if n_prbs == 48:
            h3 = h3_analysis(rows)
            write_csv(h3, out_dir / "analysis_h3_pairs.csv")
            summary["h3"] = h3
            curve_rows = read_rows(out_dir / "outage_curves_48prb.csv")
            fig_b6(curve_rows, fig_dir / "b6_outage_curves.png")

    if rows_by_prb:
        fig_collapse(rows_by_prb, fig_dir / "collapse_tw1389.png")
    calib_dir = out_dir / "calibration"
    if (calib_dir / "calib0a_nmse.csv").exists():
        fig_calib0a(calib_dir, fig_dir / "calib0a_compare.png")
    if (calib_dir / "calib0b_gaps.csv").exists():
        fig_calib0b(calib_dir, fig_dir / "calib0b_parallelism.png")

    with (out_dir / "analysis_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
