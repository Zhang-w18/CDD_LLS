from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List
import csv
import json
import math


def save_json(data: object, path: str | Path) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def save_csv(rows: List[Dict[str, object]], path: str | Path) -> None:
    if not rows:
        return
    fields = sorted(set().union(*(row.keys() for row in rows)))
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def snr_values(range_db: Iterable[float]) -> List[float]:
    start, stop, step = [float(x) for x in range_db]
    vals = []
    x = start
    if step == 0:
        raise ValueError("SNR step must be nonzero.")
    if step > 0:
        while x <= stop + 1e-9:
            vals.append(float(round(x, 10)))
            x += step
    else:
        while x >= stop - 1e-9:
            vals.append(float(round(x, 10)))
            x += step
    return vals


def interpolate_target_snr(rows: List[Dict[str, object]], target_bler: float = 0.10) -> float:
    pts = sorted((float(r["snr_db"]), float(r["bler"])) for r in rows if "snr_db" in r and "bler" in r)
    if not pts:
        return float("nan")
    for (s0, b0), (s1, b1) in zip(pts[:-1], pts[1:]):
        if (b0 - target_bler) == 0:
            return float(s0)
        if (b0 - target_bler) * (b1 - target_bler) <= 0 and b0 != b1:
            y0 = math.log10(max(b0, 1e-4))
            y1 = math.log10(max(b1, 1e-4))
            yt = math.log10(max(float(target_bler), 1e-4))
            alpha = (yt - y0) / (y1 - y0)
            return float(s0 + alpha * (s1 - s0))
    return float("nan")
