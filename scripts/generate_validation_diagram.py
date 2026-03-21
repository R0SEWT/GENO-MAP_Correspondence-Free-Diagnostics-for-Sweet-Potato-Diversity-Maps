#!/usr/bin/env python3
"""Generate a compact 2:1 horizontal validation framework diagram for the GENO-MAP poster.

Layout:
    Baseline graph  →  perturbation family  →  rebuilt variants  →  compare + metrics

Three perturbation branches are stacked vertically and then merged into one
reconstruction stage to reduce visual repetition while preserving methodology.
Designed for poster use: low cognitive load, no overlapping, clear flow.

Output:
    docs/poster/figures/fig_validation_framework_wide.pdf + .png
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import FancyBboxPatch
import numpy as np
from pathlib import Path

# ── Palette ───────────────────────────────────────────────
BLUE    = "#0072B2"
TEAL    = "#009E73"
VERM    = "#D55E00"
ROSE    = "#CC79A7"
DARK    = "#333333"
WHITE   = "#FFFFFF"
MIDGRAY = "#AAAAAA"

ITALIC_TEXT_KW = dict(
    fontsize=6.3,
    color="#262626",
    alpha=0.88,
    family="sans-serif",
    style="italic",
    fontweight="semibold",
    zorder=5,
    path_effects=[pe.withStroke(linewidth=1.25, foreground="white", alpha=0.92)],
)

SUPPORT_TEXT_KW = dict(
    fontsize=6.0,
    color="#2B2B2B",
    alpha=0.84,
    family="sans-serif",
    zorder=5,
    path_effects=[pe.withStroke(linewidth=1.0, foreground="white", alpha=0.88)],
)


def draw_mini_graph(ax, cx, cy, color, seed=42, n_nodes=8, alpha=0.62, r=0.16):
    rng = np.random.RandomState(seed)
    angles = np.linspace(0, 2 * np.pi, n_nodes, endpoint=False) + rng.uniform(0, 0.25)
    radii = r * (0.75 + 0.25 * rng.rand(n_nodes))
    xs = cx + radii * np.cos(angles)
    ys = cy + radii * np.sin(angles)

    for i in range(n_nodes):
        neigh = rng.choice([j for j in range(n_nodes) if j != i],
                           size=min(3, n_nodes - 1), replace=False)
        for j in neigh:
            ax.plot([xs[i], xs[j]], [ys[i], ys[j]],
                    color=color, lw=0.7, alpha=0.28, zorder=2)

    ax.scatter(xs, ys, s=14, color=color, alpha=alpha, zorder=3,
               edgecolors="white", linewidth=0.25)


def draw_box(ax, cx, cy, w, h, color, label, sublabel=None,
             fontsize=8.2, lw=1.6, facealpha=0.05, rounding=0.08):
    box = FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle=f"round,pad=0.05,rounding_size={rounding}",
        facecolor=color, alpha=facealpha,
        edgecolor=color, linewidth=lw, zorder=3
    )
    ax.add_patch(box)

    ax.text(
        cx, cy + (0.05 if sublabel else 0.0),
        label,
        ha="center", va="center",
        fontsize=fontsize, fontweight="bold",
        color=DARK, family="sans-serif", zorder=5,
        linespacing=1.0
    )

    if sublabel:
        ax.text(
            cx, cy - 0.18,
            sublabel,
            ha="center", va="center",
            **ITALIC_TEXT_KW,
        )


def draw_pill(ax, cx, cy, w, h, color, label, sublabel=None, fontsize=6.8):
    pill = FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.04,rounding_size=0.10",
        facecolor=color, alpha=0.10,
        edgecolor=color, linewidth=1.3, zorder=4
    )
    ax.add_patch(pill)

    ax.text(
        cx, cy + (0.03 if sublabel else 0.0),
        label,
        ha="center", va="center",
        fontsize=fontsize, fontweight="bold",
        color=color, family="sans-serif", zorder=5,
        linespacing=1.0
    )

    if sublabel:
        ax.text(
            cx, cy - 0.17,
            sublabel,
            ha="center", va="center",
            **ITALIC_TEXT_KW,
        )


def draw_arrow(ax, x0, y0, x1, y1, color=DARK, lw=1.2, alpha=0.45, rad=0.0):
    ax.annotate(
        "",
        xy=(x1, y1),
        xytext=(x0, y0),
        arrowprops=dict(
            arrowstyle="-|>",
            color=color,
            lw=lw,
            mutation_scale=11,
            alpha=alpha,
            connectionstyle=f"arc3,rad={rad}"
        ),
        zorder=2,
    )


def draw_metric_label(ax, cx, cy, color, label, sublabel):
    ax.text(
        cx, cy + 0.07,
        label,
        ha="center", va="center",
        fontsize=10.2, fontweight="bold",
        color=color, family="sans-serif", zorder=5,
        path_effects=[pe.withStroke(linewidth=1.3, foreground="white", alpha=0.92)],
    )
    ax.text(
        cx, cy - 0.10,
        sublabel,
        ha="center", va="center",
        fontsize=5.8, fontweight="semibold",
        color="#2B2B2B", family="sans-serif", zorder=5,
        linespacing=0.95,
        path_effects=[pe.withStroke(linewidth=1.0, foreground="white", alpha=0.88)],
    )


def draw_umap_snapshot(ax, cx, cy, w, h, title, seed=0, mode="baseline"):
    rng = np.random.RandomState(seed)
    n = 42

    if mode == "local":
        c1, c2, sigma = (-0.42, 0.22), (0.46, -0.20), 0.17
    elif mode == "global":
        c1, c2, sigma = (-0.18, 0.10), (0.20, -0.08), 0.29
    else:
        c1, c2, sigma = (-0.24, 0.12), (0.28, -0.10), 0.23

    pts1 = rng.normal(loc=c1, scale=sigma, size=(n, 2))
    pts2 = rng.normal(loc=c2, scale=sigma, size=(n, 2))

    sx = w * 0.34
    sy = h * 0.34

    ax.scatter(
        cx + sx * pts1[:, 0], cy + sy * pts1[:, 1],
        s=8.5, color=TEAL, alpha=0.60, edgecolors="white", linewidth=0.20, zorder=4
    )
    ax.scatter(
        cx + sx * pts2[:, 0], cy + sy * pts2[:, 1],
        s=8.5, color=ROSE, alpha=0.60, edgecolors="white", linewidth=0.20, zorder=4
    )

    ax.text(
        cx, cy + h / 2 - 0.10,
        title,
        ha="center", va="center",
        **SUPPORT_TEXT_KW,
    )


def main():
    fig, ax = plt.subplots(figsize=(13.6, 5.6))
    ax.set_xlim(0, 14.0)
    ax.set_ylim(-1.75, 1.75)
    ax.axis("off")

    # ── Geometry anchors ──────────────────────────────────
    x_baseline = 1.6
    x_perturb  = 5.0
    x_rebuilt  = 8.2
    x_compare  = 11.4

    ys = [0.95, 0.0, -0.95]

    # ══════════════════════════════════════════════════════
    # Baseline block
    # ══════════════════════════════════════════════════════
    draw_box(
        ax, x_baseline, 0.0, 2.15, 1.05, BLUE,
        "Reference\ngraph",
        sublabel="full-marker PCA-30D · kNN (k=15)",
        fontsize=9.0, lw=1.8, facealpha=0.06
    )
    draw_mini_graph(ax, x_baseline, 0.33, BLUE, seed=42, r=0.15)

    # subtle title above center area
    ax.text(
        7.0, 1.50,
        "Correspondence-free validation flow",
        ha="center", va="center",
        fontsize=10.5, fontweight="bold",
        color=DARK, family="sans-serif"
    )

    # ══════════════════════════════════════════════════════
    # Perturbation branches
    # ══════════════════════════════════════════════════════
    branches = [
        ("Marker\nsubsampling", "5–80 %", TEAL, 11),
        ("MCAR\ninjection", "+0–20 %", VERM, 22),
        ("PCA dim.\nvariation", "5–50 PCs", ROSE, 33),
    ]

    for (label, sublabel, color, seed), y in zip(branches, ys):
        # baseline -> perturbation
        draw_arrow(ax, x_baseline + 1.15, 0.0, x_perturb - 0.9, y,
                   color=color, lw=1.3, alpha=0.55,
                   rad=-0.08 * np.sign(y) if y != 0 else 0.0)

        # perturbation pill
        draw_pill(ax, x_perturb, y, 1.65, 0.62, color, label, sublabel)

        # perturbation -> merged rebuilt variants
        draw_arrow(ax, x_perturb + 0.90, y, x_rebuilt - 1.15, 0.40 * y,
                   color=color, lw=1.3, alpha=0.55)

    # merged rebuilt variants block
    rebuilt_box = FancyBboxPatch(
        (x_rebuilt - 1.075, -0.725), 2.15, 1.45,
        boxstyle="round,pad=0.05,rounding_size=0.08",
        facecolor=BLUE, alpha=0.04,
        edgecolor=BLUE, linewidth=1.6, zorder=3
    )
    ax.add_patch(rebuilt_box)
    ax.text(
        x_rebuilt, 0.58,
        "Reconstructed\nvariants",
        ha="center", va="center",
        fontsize=8.1, fontweight="bold",
        color=DARK, family="sans-serif", zorder=5
    )
    ax.text(
        x_rebuilt, -0.60,
        "same PCA-30D + kNN protocol",
        ha="center", va="center",
        **ITALIC_TEXT_KW,
    )
    for x, color, seed in zip(
        [x_rebuilt - 0.42, x_rebuilt, x_rebuilt + 0.42],
        [TEAL, VERM, ROSE],
        [11, 22, 33],
    ):
        draw_mini_graph(ax, x, 0.02, color, seed=seed, r=0.10)

    # rebuilt variants -> compare
    draw_arrow(
        ax, x_rebuilt + 1.12, 0.0, x_compare - 1.05, 0.0,
        color=DARK, lw=1.2, alpha=0.45
    )

    # annotation under branches
    ax.text(
        6.5, -1.48,
        "Additional controls: mode/median imputation and contiguous block removal",
        ha="center", va="center",
        **ITALIC_TEXT_KW,
    )

    # ══════════════════════════════════════════════════════
    # Compare panel (illustrative ADR-0014 layouts + compact metrics)
    # ══════════════════════════════════════════════════════
    panel = FancyBboxPatch(
        (x_compare - 1.50, -1.00), 3.00, 2.00,
        boxstyle="round,pad=0.06,rounding_size=0.08",
        facecolor=WHITE, edgecolor=DARK,
        linewidth=1.4, linestyle="--", alpha=0.85, zorder=2
    )
    ax.add_patch(panel)

    ax.text(
        x_compare, 0.83,
        "UMAP renderings",
        ha="center", va="center",
        fontsize=8.5, fontweight="bold",
        color=DARK, family="sans-serif", zorder=5
    )

    ax.text(
        x_compare, 0.65,
        "reference / local / global layouts",
        ha="center", va="center",
        **ITALIC_TEXT_KW,
    )

    draw_umap_snapshot(
        ax, x_compare - 0.98, 0.10, 0.82, 0.78,
        "reference", seed=101, mode="baseline"
    )
    draw_umap_snapshot(
        ax, x_compare, 0.10, 0.82, 0.78,
        "local stress", seed=202, mode="local"
    )
    draw_umap_snapshot(
        ax, x_compare + 0.98, 0.10, 0.82, 0.78,
        "global stress", seed=303, mode="global"
    )

    ax.text(
        x_compare, -0.18,
        "analytic structure remains unchanged",
        ha="center", va="center",
        **ITALIC_TEXT_KW,
    )

    draw_metric_label(ax, x_compare - 0.86, -0.56, BLUE, r"$J_{\mathrm{edge}}$", "edge\noverlap")
    draw_metric_label(ax, x_compare,        -0.56, TEAL, r"$J_{\mathrm{nbr}}$", "neighbor\noverlap")
    draw_metric_label(ax, x_compare + 0.86, -0.56, VERM, "SS", "stress\nshift")

    out_dir = Path(__file__).resolve().parent.parent / "docs" / "poster" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    for ext in ("pdf", "png"):
        fig.savefig(
            out_dir / f"fig_validation_framework_wide.{ext}",
            dpi=300,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="none"
        )

    print(f"✓ Saved to {out_dir}/fig_validation_framework_wide.{{pdf,png}}")
    plt.close(fig)


if __name__ == "__main__":
    main()
