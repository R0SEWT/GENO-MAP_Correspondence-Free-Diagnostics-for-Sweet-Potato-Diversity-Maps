"""
Problem statement figure:
Disjoint identifier namespaces prevent cross-panel alignment.

Key message:
The same biological accession may exist in both panels, but direct
sample-to-sample matching is not identifiable a priori because the
ID systems are incompatible and no lookup table is available.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent

# Okabe-Ito-inspired palette
BLUE = "#0072B2"    # Global SNP
TEAL = "#009E73"    # Global SilicoDArT
VERM = "#D55E00"    # LD SNP
ROSE = "#CC79A7"    # LD SilicoDArT
GREY = "#6F6F6F"
LIGHT = "#F4F4F4"
MID = "#D9D9D9"
GOLD = "#E69F00"
DARK = "#222222"

rng = np.random.default_rng(42)


def _scatter_in_ellipse(ax, cx, cy, rx, ry, n, marker, color, size=40, alpha=0.75):
    pts = []
    while len(pts) < n:
        x = rng.uniform(cx - 0.78 * rx, cx + 0.78 * rx)
        y = rng.uniform(cy - 0.72 * ry, cy + 0.72 * ry)
        if ((x - cx) / rx) ** 2 + ((y - cy) / ry) ** 2 < 0.62:
            pts.append((x, y))
    pts = np.array(pts)
    ax.scatter(
        pts[:, 0], pts[:, 1],
        marker=marker, s=size, color=color,
        edgecolors="white", linewidths=0.5,
        alpha=alpha, zorder=3
    )


def _blocked_link(ax, x1, y1, x2, y2, color=GOLD, lw=1.6):
    # dashed connection
    ax.plot([x1, x2], [y1, y2], linestyle=(0, (4, 3)), color=color, lw=lw, zorder=4, alpha=0.9)

    # central "X" meaning not matchable / not directly identifiable
    xm, ym = (x1 + x2) / 2, (y1 + y2) / 2
    dx, dy = 0.09, 0.09
    ax.plot([xm - dx, xm + dx], [ym - dy, ym + dy], color=color, lw=lw + 0.2, zorder=5)
    ax.plot([xm - dx, xm + dx], [ym + dy, ym - dy], color=color, lw=lw + 0.2, zorder=5)


def _labeled_point(ax, x, y, marker, facecolor, label, label_dx, label_dy,
                   edgecolor=GOLD, size=90, textcolor=DARK):
    ax.scatter(
        [x], [y], marker=marker, s=size, color=facecolor,
        edgecolors=edgecolor, linewidths=1.6, zorder=6
    )
    ax.text(
        x + label_dx, y + label_dy, label,
        fontsize=7.4, color=textcolor, ha="center", va="center",
        bbox=dict(boxstyle="round,pad=0.22", fc="white", ec=MID, lw=0.8),
        zorder=7
    )


def fig_problem_disjoint(out):
    fig, ax = plt.subplots(figsize=(8.0, 4.4))
    ax.set_xlim(-6.0, 6.0)
    ax.set_ylim(-2.9, 2.9)
    ax.set_aspect("equal")
    ax.axis("off")

    # ------------------------------------------------------------------
    # Central blocking band
    # ------------------------------------------------------------------
    
    band = mpatches.FancyBboxPatch(
    (-1.1, -1.7), 2.2, 3.4,
        boxstyle="round,pad=0.02,rounding_size=0.12",
        facecolor=LIGHT, edgecolor=MID, linewidth=1.0, zorder=0
    )
    ax.add_patch(band)

    ax.text(
        0,1.30,
        "Disjoint ID namespaces",
        fontsize=7.5,
        fontweight="bold",
        ha="center",
        linespacing=1.4
    )

    ax.text(
        0, -0.25,
        "No cross-reference table\nis available to link\nsamples across panels",
        fontsize=7.5,
        ha="center",
        linespacing=1.4,
        color=GREY,
        style="italic"
    )

    # ------------------------------------------------------------------
    # Left: Global panel
    # ------------------------------------------------------------------
    ell_l = mpatches.Ellipse(
        (-3.6, 0.0), 3.4, 3.7,
        facecolor=BLUE, alpha=0.06,
        edgecolor=BLUE, linewidth=2.0, zorder=1
    )
    ax.add_patch(ell_l)

    ax.text(-3.6, 2.15, "Global panel", fontsize=11,
            color=BLUE, fontweight="bold", ha="center")
    ax.text(-3.6, 1.9, "IDs: CIP xxxxxx.x", fontsize=7.8,
            color=GREY, ha="center")

    _scatter_in_ellipse(ax, -3.6, 0.0, 1.7, 1.85, 8, "o", BLUE, size=42, alpha=0.75)
    _scatter_in_ellipse(ax, -3.6, 0.0, 1.7, 1.85, 6, "s", TEAL, size=40, alpha=0.75)

    # Highlighted examples
    _labeled_point(ax, -2.25, 0.70, "o", BLUE, "CIP 400363", -0.02, 0.36)
    _labeled_point(ax, -2.35, -0.65, "s", TEAL, "CIP 401025", -0.02, -0.38)

    # ------------------------------------------------------------------
    # Right: Low-density panel
    # ------------------------------------------------------------------
    ell_r = mpatches.Ellipse(
        (3.6, 0.0), 3.4, 3.7,
        facecolor=VERM, alpha=0.06,
        edgecolor=VERM, linewidth=2.0, zorder=1
    )
    ax.add_patch(ell_r)

    ax.text(3.6, 2.15, "Low-density panel", fontsize=11,
            color=VERM, fontweight="bold", ha="center")
    ax.text(3.6, 1.9, "IDs: barcode_plate_well", fontsize=7.8,
            color=GREY, ha="center")

    _scatter_in_ellipse(ax, 3.6, 0.0, 1.7, 1.85, 7, "D", VERM, size=40, alpha=0.75)
    _scatter_in_ellipse(ax, 3.6, 0.0, 1.7, 1.85, 6, "^", ROSE, size=44, alpha=0.75)

    _labeled_point(ax, 2.25, 0.70, "D", VERM, "908625126003_A_7", 0.02, 0.36,
                   textcolor=DARK, size=90)
    _labeled_point(ax, 2.35, -0.65, "^", ROSE, "908625126009_G_3", 0.02, -0.38,
                   textcolor=DARK, size=90)

    # ------------------------------------------------------------------
    # Blocked correspondence attempts
    # ------------------------------------------------------------------
    _blocked_link(ax, -2.25, 0.70, 2.25, 0.70)
    _blocked_link(ax, -2.35, -0.65, 2.35, -0.65)

    ax.text(0, -1.1, "same accession?", fontsize=8.5,
            color=GOLD, ha="center", fontweight="bold")

    # ------------------------------------------------------------------
    # Legend — one entry per dataset
    # ------------------------------------------------------------------
    legend_items = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=BLUE,
                markersize=7.5, markeredgecolor="white", markeredgewidth=0.6,
                label="Global SNP"),
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=TEAL,
                markersize=7.0, markeredgecolor="white", markeredgewidth=0.6,
                label="Global SilicoDArT"),
        plt.Line2D([0], [0], marker="D", color="w", markerfacecolor=VERM,
                markersize=7.0, markeredgecolor="white", markeredgewidth=0.6,
                label="LD SNP"),
        plt.Line2D([0], [0], marker="^", color="w", markerfacecolor=ROSE,
                markersize=7.5, markeredgecolor="white", markeredgewidth=0.6,
                label="LD SilicoDArT"),
    ]

    leg = ax.legend(
        handles=legend_items,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.975),
        ncol=2,
        fontsize=6.6,
        frameon=False,
        markerscale=0.85,
        handletextpad=0.4,
        columnspacing=1.3,
        labelspacing=0.45
    )

    # Reduce saliencia del texto
    for txt in leg.get_texts():
        txt.set_alpha(0.9)
    leg.get_frame().set_facecolor("white")
    leg.get_frame().set_linewidth(0.9)
    plt.tight_layout(pad=0.4)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f" -> {out}")


if __name__ == "__main__":
    out_dir = ROOT / "docs" / "poster" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_problem_disjoint(out_dir / "fig_problem_disjoint_improved.pdf")
    print("Done.")
