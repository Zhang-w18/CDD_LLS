"""Digitize the CDF curves in two specific screenshots and draw an overlay.

The axis calibrations are based on the visible tick/grid locations in the
screenshots.  Consequently, the exported values are approximate rather than
the original numerical samples.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from matplotlib.ticker import MultipleLocator


@dataclass(frozen=True)
class Curve:
    label: str
    color: tuple[int, int, int]
    plot_color: str
    linestyle: str = "-"


@dataclass(frozen=True)
class Calibration:
    left: int
    right: int
    top: int
    bottom: int
    x_at_left: float
    x_at_right: float

    def x(self, pixel_x: np.ndarray) -> np.ndarray:
        return self.x_at_left + (pixel_x - self.left) * (
            (self.x_at_right - self.x_at_left) / (self.right - self.left)
        )

    def y(self, pixel_y: np.ndarray) -> np.ndarray:
        return (self.bottom - pixel_y) / (self.bottom - self.top)


IMAGE1_CURVES = (
    Curve("Scheme1: fix precoder", (31, 119, 180), "tab:blue"),
    Curve("Scheme2: trans smallCDD", (255, 127, 14), "tab:orange"),
    Curve("Scheme3: non-trans CDD911", (44, 160, 44), "tab:green"),
    Curve("Scheme4: non-trans CDD", (255, 0, 51), "crimson"),
)

IMAGE2_CURVES = (
    Curve("Fixed DFT0; no delay", (76, 114, 176), "#4c72b0"),
    Curve("Sidon [0,925.926,2777.778,6481.481] ns", (196, 78, 68), "#c44e44", "--"),
    Curve("CDD911 [0,0,911,911] ns", (89, 161, 79), "#59a14f", "-."),
    Curve("CDD130 [0,0,130,130] ns", (129, 114, 179), "#8172b3", ":"),
)

# Pixel locations were fitted from the visible grid/tick marks.
IMAGE1_CAL = Calibration(50, 702, 33, 382, -11.94, 7.08)
IMAGE2_CAL = Calibration(102, 1075, 38, 601, -25.0, 12.56)


def digitize(
    image_path: Path,
    calibration: Calibration,
    curves: tuple[Curve, ...],
    color_tolerance: float,
    ignore_boxes: tuple[tuple[int, int, int, int], ...] = (),
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    rgb = np.asarray(Image.open(image_path).convert("RGB"), dtype=np.int16)
    result: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    for curve in curves:
        target = np.asarray(curve.color, dtype=np.int16)
        distance = np.linalg.norm(rgb - target, axis=2)
        for left, top, right, bottom in ignore_boxes:
            distance[top : bottom + 1, left : right + 1] = np.inf
        pixel_xs: list[int] = []
        pixel_ys: list[float] = []
        for pixel_x in range(calibration.left, calibration.right + 1):
            rows = np.flatnonzero(
                distance[calibration.top : calibration.bottom + 1, pixel_x]
                <= color_tolerance
            )
            if rows.size:
                pixel_xs.append(pixel_x)
                pixel_ys.append(float(np.median(rows + calibration.top)))

        px = np.asarray(pixel_xs, dtype=float)
        py = np.asarray(pixel_ys, dtype=float)
        if px.size < 10:
            raise RuntimeError(f"Too few pixels found for {curve.label!r}: {px.size}")

        # Use a one-pixel grid and interpolate short dashed/occluded sections.
        dense_px = np.arange(int(px.min()), int(px.max()) + 1, dtype=float)
        dense_py = np.interp(dense_px, px, py)
        x = calibration.x(dense_px)
        y = np.clip(calibration.y(dense_py), 0.0, 1.0)
        y = np.maximum.accumulate(y)
        result[curve.label] = (x, y)

    return result


def write_csv(
    path: Path,
    source: str,
    values: dict[str, tuple[np.ndarray, np.ndarray]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("source", "curve", "x_db", "cdf"))
        for label, (x, y) in values.items():
            writer.writerows(
                (source, label, f"{x_value:.6f}", f"{y_value:.6f}")
                for x_value, y_value in zip(x, y, strict=True)
            )


def draw_overlay(
    output_path: Path,
    image1_values: dict[str, tuple[np.ndarray, np.ndarray]],
    image2_values: dict[str, tuple[np.ndarray, np.ndarray]],
) -> None:
    pair_indices = ((0, 0), (1, 3), (2, 2), (3, 1))
    panel_titles = (
        "Scheme1 vs Fixed DFT0",
        "Scheme2 vs CDD130",
        "Scheme3 vs CDD911",
        "Scheme4 vs Sidon",
    )
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True, sharey=True, constrained_layout=True)
    for ax, (image1_index, image2_index), title in zip(
        axes.flat, pair_indices, panel_titles, strict=True
    ):
        curve1 = IMAGE1_CURVES[image1_index]
        curve2 = IMAGE2_CURVES[image2_index]
        x1, y1 = image1_values[curve1.label]
        x2, y2 = image2_values[curve2.label]
        ax.plot(x1, y1, color=curve1.plot_color, linewidth=2.2, label=curve1.label)
        ax.plot(
            x2,
            y2,
            color=curve2.plot_color,
            linestyle="--",
            linewidth=2.5,
            label=curve2.label,
        )
        ax.set_title(title)
        ax.set_xlim(-25, 12.5)
        ax.set_ylim(-0.01, 1.01)
        ax.grid(True, alpha=0.3)
        ax.xaxis.set_minor_locator(MultipleLocator(2.5))
        ax.grid(True, which="minor", axis="x", alpha=0.16, linewidth=0.8)
        ax.legend(loc="upper left", fontsize=9, frameon=True)

    for ax in axes[-1, :]:
        ax.set_xlabel("Power (dB)")
    for ax in axes[:, 0]:
        ax.set_ylabel("CDF")
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("image1", type=Path)
    parser.add_argument("image2", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    image1_values = digitize(
        args.image1,
        IMAGE1_CAL,
        IMAGE1_CURVES,
        color_tolerance=80.0,
        ignore_boxes=((515, 313, 700, 381),),
    )
    image2_values = digitize(args.image2, IMAGE2_CAL, IMAGE2_CURVES, color_tolerance=65.0)
    write_csv(args.output_dir / "image1_digitized.csv", "image1", image1_values)
    write_csv(args.output_dir / "image2_digitized.csv", "image2", image2_values)
    draw_overlay(args.output_dir / "cdf_overlay.png", image1_values, image2_values)


if __name__ == "__main__":
    main()
