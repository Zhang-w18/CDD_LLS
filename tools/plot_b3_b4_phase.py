"""Plot representative B3 and B4 frequency-domain phase laws.

The construction matches tools/run_track_b_pilot_scan.py exactly. The plotted
phase is unwrapped along the subcarrier axis so a physically irrelevant 2*pi
representation jump is not mistaken for a discontinuity.
"""

from __future__ import annotations

from pathlib import Path
import html
import math

import numpy as np
from PIL import Image, ImageDraw, ImageFont, JpegImagePlugin  # noqa: F401


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "figures" / "b3_b4_phase"
K = 576
NT = 8
BITREV = [0, 4, 2, 6, 1, 5, 3, 7]

WIDTH, HEIGHT = 1800, 2300
MARGIN_X, TOP, BOTTOM = 110, 170, 90
COL_GAP, ROW_GAP = 105, 82
COLS, ROWS = 2, 4
PLOT_W = (WIDTH - 2 * MARGIN_X - COL_GAP) // COLS
PLOT_H = (HEIGHT - TOP - BOTTOM - (ROWS - 1) * ROW_GAP) // ROWS


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    roots = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold
        else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/DejaVuSans-Bold.ttf" if bold
        else "/System/Library/Fonts/Supplemental/DejaVuSans.ttf",
    ]
    for path in roots:
        if Path(path).exists():
            return ImageFont.truetype(path, size=size)
    return ImageFont.load_default()


def wrap_pi(x: np.ndarray) -> np.ndarray:
    return -(np.mod(-x + math.pi, 2.0 * math.pi) - math.pi)


def build_b3_phase(
    k_count: int = K,
    n_seg: int = 8,
    transition_width: int = 24,
    schedule: str = "seq",
) -> np.ndarray:
    seg_len = k_count // n_seg
    n = np.arange(NT, dtype=np.float64)
    if schedule == "seq":
        r = np.array([s % 8 for s in range(n_seg)], dtype=np.float64)
    else:
        r = np.array([BITREV[s % 8] for s in range(n_seg)], dtype=np.float64)
    plateau = 2.0 * math.pi * n[None, :] * r[:, None] / 8.0
    k = np.arange(k_count)
    phase = plateau[k // seg_len, :].copy()
    half = transition_width // 2
    for s in range(1, n_seg):
        boundary = s * seg_len
        transition_k = np.arange(boundary - half, boundary + half)
        delta = wrap_pi(plateau[s, :] - plateau[s - 1, :])
        fraction = (transition_k - (boundary - half)) / float(transition_width)
        phase[transition_k, :] = (
            plateau[s - 1, :][None, :]
            + delta[None, :] * fraction[:, None]
        )
    return np.unwrap(phase, axis=0)


def build_b4_phase(
    k_count: int = K,
    slope_cap: int = 2,
    seed_offset: int = 21,
) -> tuple[np.ndarray, np.ndarray]:
    n_seg = 8
    seg_len = k_count // n_seg
    rng = np.random.default_rng([23, slope_cap, seed_offset])
    slope_index = rng.integers(0, slope_cap + 1, size=(n_seg, NT))
    slope = slope_index[np.arange(k_count) // seg_len, :].astype(np.float64)
    phase = np.zeros((k_count, NT), dtype=np.float64)
    phase[1:, :] = -2.0 * math.pi / k_count * np.cumsum(slope[:-1, :], axis=0)
    return phase, slope_index


def nice_limits(values: np.ndarray) -> tuple[float, float]:
    lo, hi = float(values.min()), float(values.max())
    if hi - lo < 0.25:
        mid = 0.5 * (lo + hi)
        return mid - 0.5, mid + 0.5
    pad = 0.08 * (hi - lo)
    return lo - pad, hi + pad


def panel_box(index: int) -> tuple[int, int, int, int]:
    row, col = divmod(index, COLS)
    left = MARGIN_X + col * (PLOT_W + COL_GAP)
    top = TOP + row * (PLOT_H + ROW_GAP)
    return left, top, left + PLOT_W, top + PLOT_H


def dashed_vertical(draw: ImageDraw.ImageDraw, x: int, y0: int, y1: int) -> None:
    for y in range(y0, y1, 15):
        draw.line((x, y, x, min(y + 8, y1)), fill=(95, 95, 95), width=2)


def make_png(
    phase: np.ndarray,
    title: str,
    subtitle: str,
    color: tuple[int, int, int],
    kind: str,
    slopes: np.ndarray | None = None,
) -> Image.Image:
    data = phase / math.pi
    ymin, ymax = nice_limits(data)
    image = Image.new("RGB", (WIDTH, HEIGHT), "white")
    draw = ImageDraw.Draw(image, "RGBA")
    draw.text((WIDTH // 2, 38), title, font=font(36, True), fill=(25, 25, 25), anchor="ma")
    draw.text((WIDTH // 2, 92), subtitle, font=font(24), fill=(70, 70, 70), anchor="ma")
    seg_len = K // 8
    half = 12

    for ant in range(NT):
        left, top, right, bottom = panel_box(ant)
        # Background and neutral grid.
        draw.rectangle((left, top, right, bottom), fill=(252, 252, 252), outline=(80, 80, 80), width=2)
        for t in (0.0, 0.5, 1.0):
            y = int(bottom - t * (bottom - top))
            draw.line((left, y, right, y), fill=(215, 215, 215), width=1)
            value = ymin + t * (ymax - ymin)
            draw.text((left - 12, y), f"{value:.2f}", font=font(18), fill=(70, 70, 70), anchor="rm")
        for xtick in (0, 144, 288, 432, 575):
            x = int(left + xtick / (K - 1) * (right - left))
            draw.line((x, top, x, bottom), fill=(225, 225, 225), width=1)
            draw.text((x, bottom + 8), str(xtick), font=font(17), fill=(70, 70, 70), anchor="ma")

        if kind == "b3":
            for s in range(1, 8):
                boundary = s * seg_len
                x0 = int(left + (boundary - half) / (K - 1) * (right - left))
                x1 = int(left + (boundary + half) / (K - 1) * (right - left))
                draw.rectangle((x0, top, x1, bottom), fill=(249, 168, 37, 30))
                xb = int(left + boundary / (K - 1) * (right - left))
                dashed_vertical(draw, xb, top, bottom)
        else:
            for s in range(1, 8):
                xb = int(left + s * seg_len / (K - 1) * (right - left))
                dashed_vertical(draw, xb, top, bottom)

        pts = []
        for k, value in enumerate(data[:, ant]):
            x = left + k / (K - 1) * (right - left)
            y = bottom - (value - ymin) / (ymax - ymin) * (bottom - top)
            pts.append((x, y))
        draw.line(pts, fill=(*color, 255), width=4, joint="curve")
        draw.text((left + 10, top + 9), f"Antenna {ant}", font=font(22, True), fill=(30, 30, 30))
        if slopes is not None:
            labels = " ".join(str(int(v)) for v in slopes[:, ant])
            draw.text((left + 10, bottom - 34), f"j_s = [{labels}]", font=font(17), fill=(50, 50, 50))

    draw.text((WIDTH // 2, HEIGHT - 42), "Frequency-domain subcarrier index k", font=font(25), fill=(35, 35, 35), anchor="ms")
    return image


def svg_text(x: float, y: float, value: str, size: int = 18, anchor: str = "start", weight: int = 400) -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" '
        f'fill="#333">{html.escape(value)}</text>'
    )


def make_svg(
    phase: np.ndarray,
    title: str,
    subtitle: str,
    color: str,
    kind: str,
    slopes: np.ndarray | None = None,
) -> str:
    data = phase / math.pi
    ymin, ymax = nice_limits(data)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">',
        '<rect width="100%" height="100%" fill="white"/>',
        svg_text(WIDTH / 2, 55, title, 36, "middle", 700),
        svg_text(WIDTH / 2, 105, subtitle, 24, "middle"),
    ]
    seg_len, half = K // 8, 12
    for ant in range(NT):
        left, top, right, bottom = panel_box(ant)
        parts.append(f'<rect x="{left}" y="{top}" width="{right-left}" height="{bottom-top}" fill="#fcfcfc" stroke="#555" stroke-width="2"/>')
        for t in (0.0, 0.5, 1.0):
            y = bottom - t * (bottom - top)
            parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}" stroke="#d7d7d7"/>')
            parts.append(svg_text(left - 12, y + 6, f"{ymin + t*(ymax-ymin):.2f}", 18, "end"))
        for xtick in (0, 144, 288, 432, 575):
            x = left + xtick / (K - 1) * (right - left)
            parts.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{bottom}" stroke="#e1e1e1"/>')
            parts.append(svg_text(x, bottom + 27, str(xtick), 17, "middle"))
        if kind == "b3":
            for s in range(1, 8):
                boundary = s * seg_len
                x0 = left + (boundary - half) / (K - 1) * (right - left)
                x1 = left + (boundary + half) / (K - 1) * (right - left)
                xb = left + boundary / (K - 1) * (right - left)
                parts.append(f'<rect x="{x0:.1f}" y="{top}" width="{x1-x0:.1f}" height="{bottom-top}" fill="#f9a825" opacity="0.13"/>')
                parts.append(f'<line x1="{xb:.1f}" y1="{top}" x2="{xb:.1f}" y2="{bottom}" stroke="#666" stroke-width="2" stroke-dasharray="8 7"/>')
        else:
            for s in range(1, 8):
                xb = left + s * seg_len / (K - 1) * (right - left)
                parts.append(f'<line x1="{xb:.1f}" y1="{top}" x2="{xb:.1f}" y2="{bottom}" stroke="#666" stroke-width="2" stroke-dasharray="8 7"/>')
        points = []
        for k, value in enumerate(data[:, ant]):
            x = left + k / (K - 1) * (right - left)
            y = bottom - (value - ymin) / (ymax - ymin) * (bottom - top)
            points.append(f"{x:.1f},{y:.1f}")
        parts.append(f'<polyline points="{" ".join(points)}" fill="none" stroke="{color}" stroke-width="4" stroke-linejoin="round"/>')
        parts.append(svg_text(left + 10, top + 29, f"Antenna {ant}", 22, "start", 700))
        if slopes is not None:
            labels = " ".join(str(int(v)) for v in slopes[:, ant])
            parts.append(svg_text(left + 10, bottom - 13, f"j_s = [{labels}]", 17))
    parts.append(svg_text(WIDTH / 2, HEIGHT - 28, "Frequency-domain subcarrier index k", 25, "middle"))
    parts.append("</svg>")
    return "\n".join(parts)


def save_all() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    b3_phase = build_b3_phase()
    b4_phase, b4_slopes = build_b4_phase()

    specs = [
        (
            "b3-cc-phase",
            b3_phase,
            "B3 / CC: phase-continuous cycling",
            "Vertical axis: unwrapped phase / pi; K=576, Nseg=8, T=24; dashed = segment boundary; shaded = transition band",
            (21, 101, 192),
            "#1565c0",
            "b3",
            None,
        ),
        (
            "b4-cn-phase",
            b4_phase,
            "B4 / CN: continuous phase with capped piecewise-constant slope",
            "Vertical axis: unwrapped phase / pi; K=576, a=2, seed=21; dashed = slope-change boundary",
            (123, 31, 162),
            "#7b1fa2",
            "b4",
            b4_slopes,
        ),
    ]
    png_images = []
    for stem, phase, title, subtitle, rgb, hex_color, kind, slopes in specs:
        image = make_png(phase, title, subtitle, rgb, kind, slopes)
        image.save(OUT / f"{stem}.png", dpi=(220, 220))
        (OUT / f"{stem}.svg").write_text(
            make_svg(phase, title, subtitle, hex_color, kind, slopes),
            encoding="utf-8",
        )
        png_images.append(image.convert("RGB"))

    png_images[0].save(
        OUT / "b3-b4-phase-combined.pdf",
        save_all=True,
        append_images=png_images[1:],
        resolution=220.0,
    )
    np.savez(
        OUT / "b3-b4-phase-data.npz",
        b3_phase=b3_phase,
        b4_phase=b4_phase,
        b4_slope_index=b4_slopes,
    )


if __name__ == "__main__":
    save_all()
