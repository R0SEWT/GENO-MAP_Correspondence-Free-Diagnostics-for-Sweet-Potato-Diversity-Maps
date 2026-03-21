#!/usr/bin/env python3
"""Generate a high-quality pipeline flow diagram for the GENO-MAP poster.

Horizontal pipeline:
  Genotype Matrix → Mode Imputation → PCA-30D → kNN (k=15) → UMAP-2D

Style: Okabe-Ito palette, rounded boxes, clean arrows, annotations.
Output: docs/poster/figures/fig_pipeline.pdf + .png
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
from pathlib import Path

# ── Palette (Okabe-Ito) ──────────────────────────────────
BLUE   = "#0072B2"
TEAL   = "#009E73"
VERM   = "#D55E00"
ROSE   = "#CC79A7"
GOLD   = "#E69F00"
DARK   = "#333333"
LIGHT  = "#F7F7F7"
MIDGRAY = "#BBBBBB"
WHITE  = "#FFFFFF"

# ── Layout parameters ────────────────────────────────────
N_STEPS = 5
BOX_W = 1.8
BOX_H = 1.1
GAP = 0.65   # horizontal gap between boxes
ARROW_LEN = GAP
Y_CENTER = 0.0
TOTAL_W = N_STEPS * BOX_W + (N_STEPS - 1) * GAP
X_START = 0.0

# Step definitions: (label_line1, label_line2, color, annotation, icon_type)
STEPS = [
    ("Genotype\nMatrix",    None,         MIDGRAY, "$n \\times p$\nDArT / DArTSeq", "matrix"),
    ("Mode\nImputation",    None,         TEAL,    "NaN → mode\nper marker",         "clean"),
    ("PCA",                 "30D",        BLUE,    "Analytic\nspace",                "reduce"),
    ("kNN Graph",           "$k\\!=\\!15$", VERM,  "Cosine\ndistance",               "graph"),
    ("UMAP",                "2D",         ROSE,    "Rendering\nonly",                 "scatter"),
]


def draw_mini_icon(ax, cx, cy, icon_type, color, size=0.18):
    """Draw a tiny icon inside each pipeline box."""
    s = size
    if icon_type == "matrix":
        # 3x4 grid of dots
        for r in range(3):
            for c in range(4):
                ax.plot(cx - s*0.6 + c * s*0.4, cy - s*0.3 + r * s*0.3,
                        'o', color=color, markersize=2.5, alpha=0.5)
    elif icon_type == "clean":
        # Sparkle / broom
        ax.plot([cx - s*0.3, cx + s*0.3], [cy, cy],
                color=color, lw=1.5, alpha=0.5)
        ax.plot([cx, cx], [cy - s*0.3, cy + s*0.3],
                color=color, lw=1.5, alpha=0.5)
        ax.plot(cx, cy, '*', color=color, markersize=8, alpha=0.6)
    elif icon_type == "reduce":
        # Converging arrows (many dims → few)
        for dy in [-s*0.35, -s*0.15, 0, s*0.15, s*0.35]:
            ax.annotate("", xy=(cx + s*0.15, cy),
                       xytext=(cx - s*0.35, cy + dy),
                       arrowprops=dict(arrowstyle="-", color=color,
                                      alpha=0.35, lw=0.8))
        ax.plot(cx + s*0.25, cy, 'o', color=color, markersize=4, alpha=0.6)
    elif icon_type == "graph":
        # Small network: 4 nodes + edges
        nodes = [(cx-s*0.25, cy+s*0.2), (cx+s*0.25, cy+s*0.2),
                 (cx-s*0.15, cy-s*0.2), (cx+s*0.15, cy-s*0.2)]
        edges = [(0,1), (0,2), (1,3), (2,3), (0,3)]
        for i, j in edges:
            ax.plot([nodes[i][0], nodes[j][0]],
                    [nodes[i][1], nodes[j][1]],
                    color=color, lw=0.8, alpha=0.4)
        for nx, ny in nodes:
            ax.plot(nx, ny, 'o', color=color, markersize=3.5, alpha=0.7)
    elif icon_type == "scatter":
        # Small scatter
        rng = np.random.RandomState(42)
        xs = cx + rng.randn(12) * s * 0.25
        ys = cy + rng.randn(12) * s * 0.2
        ax.scatter(xs, ys, s=4, color=color, alpha=0.5, zorder=5)


def main():
    fig, ax = plt.subplots(figsize=(14, 3.2))
    ax.set_xlim(-0.5, TOTAL_W + 0.5)
    ax.set_ylim(-1.6, 1.4)
    ax.set_aspect("equal")
    ax.axis("off")

    # ── Title ─────────────────────────────────────────────
    ax.text(TOTAL_W / 2, 1.25, "Analytical Pipeline",
            ha="center", va="center", fontsize=16, fontweight="bold",
            color=DARK, family="sans-serif")

    # ── Draw steps ────────────────────────────────────────
    box_centers_x = []
    for i, (label, sublabel, color, annotation, icon_type) in enumerate(STEPS):
        cx = X_START + i * (BOX_W + GAP) + BOX_W / 2
        cy = Y_CENTER
        box_centers_x.append(cx)

        # ── Box ───────────────────────────────────────────
        box = FancyBboxPatch(
            (cx - BOX_W/2, cy - BOX_H/2), BOX_W, BOX_H,
            boxstyle="round,pad=0.08",
            facecolor=WHITE,
            edgecolor=color,
            linewidth=2.5,
            zorder=3,
        )
        ax.add_patch(box)

        # ── Step number badge ─────────────────────────────
        badge_r = 0.16
        badge_x = cx - BOX_W/2 + 0.08
        badge_y = cy + BOX_H/2 - 0.08
        badge = plt.Circle((badge_x, badge_y), badge_r,
                           color=color, zorder=5)
        ax.add_patch(badge)
        ax.text(badge_x, badge_y, str(i + 1),
                ha="center", va="center", fontsize=9, fontweight="bold",
                color=WHITE, zorder=6)

        # ── Icon (top half of box) ────────────────────────
        draw_mini_icon(ax, cx, cy + 0.15, icon_type, color, size=0.22)

        # ── Label (bottom half of box) ────────────────────
        label_y = cy - 0.25
        ax.text(cx, label_y, label,
                ha="center", va="center", fontsize=10, fontweight="bold",
                color=DARK, family="sans-serif", linespacing=1.1)

        # ── Sublabel ──────────────────────────────────────
        if sublabel:
            ax.text(cx + BOX_W/2 - 0.12, cy + BOX_H/2 - 0.12, sublabel,
                    ha="center", va="center", fontsize=9,
                    color=color, family="sans-serif",
                    fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.08",
                              facecolor=color, alpha=0.12, edgecolor="none"))

        # ── Annotation below box ──────────────────────────
        ax.text(cx, cy - BOX_H/2 - 0.22, annotation,
                ha="center", va="top", fontsize=8,
                color=DARK, alpha=0.6, family="sans-serif",
                linespacing=1.2, style="italic")

    # ── Arrows between steps ──────────────────────────────
    for i in range(N_STEPS - 1):
        x_from = box_centers_x[i] + BOX_W / 2 + 0.04
        x_to = box_centers_x[i + 1] - BOX_W / 2 - 0.04
        arrow = FancyArrowPatch(
            (x_from, Y_CENTER), (x_to, Y_CENTER),
            arrowstyle="-|>",
            mutation_scale=18,
            color=DARK,
            lw=2.0,
            alpha=0.5,
            zorder=2,
        )
        ax.add_patch(arrow)

    # ── Bracket: "Deterministic analytic layer" over steps 2–4 ──
    bx_left = box_centers_x[1] - BOX_W/2
    bx_right = box_centers_x[3] + BOX_W/2
    by = Y_CENTER + BOX_H/2 + 0.25
    bracket_h = 0.12
    # Draw bracket lines
    ax.plot([bx_left, bx_left, bx_right, bx_right],
            [by, by + bracket_h, by + bracket_h, by],
            color=BLUE, lw=1.5, alpha=0.6)
    ax.text((bx_left + bx_right) / 2, by + bracket_h + 0.06,
            "Deterministic analytic layer",
            ha="center", va="bottom", fontsize=9, color=BLUE,
            fontweight="bold", alpha=0.7, family="sans-serif")

    # ── "Rendering only" label for UMAP ───────────────────
    ux = box_centers_x[4]
    uy = Y_CENTER + BOX_H/2 + 0.25
    ax.annotate("Perceptual\nrendering",
                xy=(ux, uy - 0.05),
                ha="center", va="bottom", fontsize=8,
                color=ROSE, fontweight="bold", alpha=0.7,
                family="sans-serif")

    # ── Save ──────────────────────────────────────────────
    out_dir = Path(__file__).resolve().parent.parent / "docs" / "poster" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out_dir / f"fig_pipeline.{ext}",
                    dpi=300, bbox_inches="tight",
                    facecolor="white", edgecolor="none")
    print(f"✓ Pipeline diagram saved to {out_dir}/fig_pipeline.{{pdf,png}}")
    plt.close(fig)


if __name__ == "__main__":
    main()
