#!/usr/bin/env python3
"""
generate_beta_figures.py — Additional figures for poster_a1_beta.tex

New figures from ADR experiments:
  Fig 7  — PC sensitivity curve (ADR-0012)
  Fig 8  — Block vs Random removal paired bars (ADR-0013)

Uses same Frontiers style as generate_poster_figures.py.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

ROOT = Path(__file__).resolve().parent.parent

# ── Frontiers house style ───────────────────────────────────────────
DPI = 300
SINGLE_COL_MM = 85
DOUBLE_COL_MM = 170
MM = 1 / 25.4

# Okabe-Ito palette
C = {
    "global_snp":    "#0072B2",
    "global_silico": "#009E73",
    "ld_snp":        "#D55E00",
    "ld_silico":     "#CC79A7",
}
DATASET_ORDER = ["global_snp", "global_silico", "ld_snp", "ld_silico"]
DS_SHORT = ["G-SNP", "G-Silico", "LD-SNP", "LD-Silico"]
DS_LABEL = {
    "global_snp":    "Global SNP",
    "global_silico": "Global SilicoDArT",
    "ld_snp":        "LowDensity SNP",
    "ld_silico":     "LowDensity SilicoDArT",
}
MARKERS = ["o", "s", "D", "^"]


def _style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "lines.linewidth": 1.4,
        "lines.markersize": 5,
        "legend.frameon": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })


def _lbl(ax, s, x=-0.12, y=1.08):
    ax.text(x, y, s, transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top", ha="left")


# ════════════════════════════════════════════════════════════════════
# DATA from experiments
# ════════════════════════════════════════════════════════════════════

# ADR-0012: mean J_edge per PCA dim (seed mean)
PC_DIMS = [5, 10, 15, 20, 30, 50]
PC_JEDGE = {
    "global_snp":    [0.338, 0.472, 0.592, 0.699, 0.930, 0.653],
    "global_silico": [0.340, 0.474, 0.589, 0.694, 0.949, 0.667],
    "ld_snp":        [0.343, 0.471, 0.616, 0.690, 0.955, 0.650],
    "ld_silico":     [0.346, 0.503, 0.610, 0.705, 0.957, 0.667],
}
PC_SS = {
    "global_snp":    [1.000, 1.000, 0.999, 0.997, 0.988, 0.997],
    "global_silico": [1.000, 1.000, 1.000, 0.999, 0.996, 0.997],
    "ld_snp":        [1.000, 1.000, 1.000, 1.000, 0.997, 0.998],
    "ld_silico":     [1.000, 1.000, 1.000, 1.000, 0.998, 0.999],
}

# ADR-0013: block vs random at 10% removal  (mean over seeds/positions)
BLOCK_DATA = {
    "global_snp":    {"random": 0.875, "block": 0.884},
    "global_silico": {"random": 0.911, "block": 0.896},
    "ld_snp":        {"random": 0.918, "block": 0.912},
    "ld_silico":     {"random": 0.913, "block": 0.883},
}

# ADR-0013: block position at 10% (seed 42)
BLOCK_POS = {
    "global_snp":    {"start": 0.861, "center": 0.875, "end": 0.917},
    "global_silico": {"start": 0.872, "center": 0.889, "end": 0.927},
    "ld_snp":        {"start": 0.937, "center": 0.913, "end": 0.887},
    "ld_silico":     {"start": 0.850, "center": 0.896, "end": 0.902},
}

# ADR-0014: UMAP invariance proof
UMAP_PCA_INVARIANCE = 1.0000  # across all configs
UMAP_JEDGE = {
    "global_snp": 0.497,
    "global_silico": 0.492,
    "ld_snp": 0.565,
    "ld_silico": 0.575,
}

# ADR-0015: stability gap
STAB_GAP = {
    "global_snp":    {"pca": 0.904, "ae": 0.524, "gap": 0.380},
    "global_silico": {"pca": 0.919, "ae": 0.498, "gap": 0.421},
    "ld_snp":        {"pca": 0.935, "ae": 0.163, "gap": 0.772},
    "ld_silico":     {"pca": 0.932, "ae": 0.207, "gap": 0.725},
}


# ════════════════════════════════════════════════════════════════════
# FIGURES
# ════════════════════════════════════════════════════════════════════

def fig7_pc_sensitivity(out: Path):
    """ADR-0012: J_edge vs PCA dimensionality — shows 30D is the sweet spot."""
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(DOUBLE_COL_MM * MM, 80 * MM))

    # (A) J_edge vs PCs
    for i, ds in enumerate(DATASET_ORDER):
        ax1.plot(PC_DIMS, PC_JEDGE[ds],
                 marker=MARKERS[i], color=C[ds],
                 ms=6, lw=1.8,
                 markeredgecolor="white", markeredgewidth=0.5)

    # Highlight 30D sweet spot
    ax1.axvline(30, ls="--", lw=1.0, color="#888", alpha=0.7, zorder=0)
    ax1.annotate("sweet spot\n(30D)", xy=(30, 0.96),
                 xytext=(38, 0.50), fontsize=7, color="#555",
                 ha="left", va="center",
                 arrowprops=dict(arrowstyle="->", color="#999", lw=0.7))

    # Semantic zones
    ax1.axhspan(0.80, 1.05, color="#009E73", alpha=0.06, zorder=0)
    ax1.axhspan(0.60, 0.80, color="#E69F00", alpha=0.06, zorder=0)
    ax1.axhspan(0.00, 0.60, color="#D55E00", alpha=0.06, zorder=0)

    ax1.set_xlabel("PCA dimensions")
    ax1.set_ylabel("$J_{\\mathrm{edge}}$ vs baseline (PCA-30D)")
    ax1.set_ylim(0.25, 1.05)
    ax1.set_xlim(0, 55)
    ax1.xaxis.set_major_locator(mticker.FixedLocator(PC_DIMS))
    ax1.grid(axis="y", lw=0.3, alpha=0.4)
    _lbl(ax1, "A")

    # (B) SS vs PCs — all near 1.0
    for i, ds in enumerate(DATASET_ORDER):
        ax2.plot(PC_DIMS, PC_SS[ds],
                 marker=MARKERS[i], color=C[ds],
                 ms=6, lw=1.8,
                 markeredgecolor="white", markeredgewidth=0.5)

    ax2.set_xlabel("PCA dimensions")
    ax2.set_ylabel("Subspace Similarity (SS)")
    ax2.set_ylim(0.97, 1.005)
    ax2.set_xlim(0, 55)
    ax2.xaxis.set_major_locator(mticker.FixedLocator(PC_DIMS))
    ax2.grid(axis="y", lw=0.3, alpha=0.4)
    ax2.axhline(0.99, ls=":", lw=0.5, color="#999")
    _lbl(ax2, "B")

    fig.tight_layout(w_pad=2.0)
    fig.savefig(out / "fig7_pc_sensitivity.png")
    fig.savefig(out / "fig7_pc_sensitivity.pdf")
    plt.close(fig)
    print(f"  ✓ fig7_pc_sensitivity")


def fig8_block_removal(out: Path):
    """ADR-0013: Block vs random paired bars + position breakdown."""
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(DOUBLE_COL_MM * MM, 80 * MM))

    x = np.arange(4)
    w = 0.35
    panel_colors = [C[ds] for ds in DATASET_ORDER]

    # Lighter versions for block bars
    def _light(hex_c, t=0.35):
        r, g, b = int(hex_c[1:3], 16), int(hex_c[3:5], 16), int(hex_c[5:7], 16)
        r2 = int(r + (255 - r) * t)
        g2 = int(g + (255 - g) * t)
        b2 = int(b + (255 - b) * t)
        return f"#{r2:02x}{g2:02x}{b2:02x}"
    block_colors = [_light(c) for c in panel_colors]

    # (A) Random vs Block at 10%
    rand_vals = [BLOCK_DATA[ds]["random"] for ds in DATASET_ORDER]
    block_vals = [BLOCK_DATA[ds]["block"] for ds in DATASET_ORDER]

    ax1.bar(x - w/2, rand_vals, w, color=panel_colors,
            edgecolor="white", linewidth=0.4, label="Random")
    ax1.bar(x + w/2, block_vals, w, color=block_colors,
            edgecolor=[c for c in panel_colors], linewidth=0.8,
            hatch="//", label="Contiguous block")

    for i in range(4):
        ax1.text(x[i] - w/2, rand_vals[i] + 0.004, f"{rand_vals[i]:.3f}",
                 ha="center", va="bottom", fontsize=5.5, color="#333")
        ax1.text(x[i] + w/2, block_vals[i] + 0.004, f"{block_vals[i]:.3f}",
                 ha="center", va="bottom", fontsize=5.5, color="#333")

    ax1.set_ylabel("$J_{\\mathrm{edge}}$ (10% markers removed)")
    ax1.set_xticks(x)
    ax1.set_xticklabels(DS_SHORT, rotation=30, ha="right")
    ax1.set_ylim(0.80, 0.97)
    ax1.yaxis.set_major_locator(mticker.MultipleLocator(0.05))
    ax1.axhline(0.85, ls=":", lw=0.6, color="#999", alpha=0.6)

    from matplotlib.patches import Patch
    ax1.legend(handles=[
        Patch(facecolor="#888", edgecolor="white", label="Random"),
        Patch(facecolor="#ccc", edgecolor="#888", hatch="//", label="Contiguous block"),
    ], loc="lower left", fontsize=7)
    _lbl(ax1, "A")

    # (B) Block position breakdown
    positions = ["start", "center", "end"]
    pos_x = np.arange(3)
    bar_w = 0.18

    for i, ds in enumerate(DATASET_ORDER):
        vals = [BLOCK_POS[ds][p] for p in positions]
        offset = (i - 1.5) * bar_w
        ax2.bar(pos_x + offset, vals, bar_w, color=C[ds],
                edgecolor="white", linewidth=0.4, label=DS_SHORT[i])

    ax2.set_ylabel("$J_{\\mathrm{edge}}$ (block removal, 10%)")
    ax2.set_xticks(pos_x)
    ax2.set_xticklabels(["Start", "Center", "End"])
    ax2.set_ylim(0.80, 0.97)
    ax2.yaxis.set_major_locator(mticker.MultipleLocator(0.05))
    ax2.axhline(0.85, ls=":", lw=0.6, color="#999", alpha=0.6)
    ax2.legend(loc="lower right", fontsize=7, ncol=2)
    _lbl(ax2, "B")

    fig.tight_layout(w_pad=2.0)
    fig.savefig(out / "fig8_block_removal.png")
    fig.savefig(out / "fig8_block_removal.pdf")
    plt.close(fig)
    print(f"  ✓ fig8_block_removal")


def fig9_unified_stability(out: Path):
    """ADR-0015: Unified stability comparison PCA vs AE — horizontal bars."""
    fig, ax = plt.subplots(figsize=(SINGLE_COL_MM * 1.1 * MM,
                                    SINGLE_COL_MM * 0.9 * MM))

    y = np.arange(4)
    h = 0.35
    panel_colors = [C[ds] for ds in DATASET_ORDER]

    pca_vals = [STAB_GAP[ds]["pca"] for ds in DATASET_ORDER]
    ae_vals = [STAB_GAP[ds]["ae"] for ds in DATASET_ORDER]

    ax.barh(y - h/2, pca_vals, h, color=panel_colors,
            edgecolor="white", linewidth=0.4, label="PCA-30D")
    ax.barh(y + h/2, ae_vals, h, color=[C[ds] for ds in DATASET_ORDER],
            edgecolor=[C[ds] for ds in DATASET_ORDER], linewidth=0.8,
            alpha=0.35, hatch="//", label="AE-64D")

    # Gap annotations
    for i, ds in enumerate(DATASET_ORDER):
        gap = STAB_GAP[ds]["gap"]
        mid_x = (pca_vals[i] + ae_vals[i]) / 2
        ax.annotate("", xy=(ae_vals[i] + 0.01, y[i]),
                     xytext=(pca_vals[i] - 0.01, y[i]),
                     arrowprops=dict(arrowstyle="<->", color="#555", lw=0.7))
        ax.text(mid_x, y[i] + 0.38, f"Δ={gap:.2f}",
                fontsize=6, color="#555", ha="center", va="bottom")

    ax.set_yticks(y)
    ax.set_yticklabels(DS_SHORT)
    ax.set_xlabel("Inter-seed $J_{\\mathrm{edge}}$")
    ax.set_xlim(0, 1.05)
    ax.xaxis.set_major_locator(mticker.MultipleLocator(0.2))
    ax.grid(axis="x", lw=0.3, alpha=0.4)
    ax.invert_yaxis()

    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor="#888", edgecolor="white", label="PCA-30D"),
        Patch(facecolor="#888", edgecolor="#888", alpha=0.35, hatch="//", label="AE-64D"),
    ], loc="lower right", fontsize=7)

    fig.tight_layout()
    fig.savefig(out / "fig9_unified_stability.png")
    fig.savefig(out / "fig9_unified_stability.pdf")
    plt.close(fig)
    print(f"  ✓ fig9_unified_stability")


def main():
    ap = argparse.ArgumentParser(
        description="Generate beta poster figures")
    ap.add_argument("--outdir", type=Path,
                    default=ROOT / "docs" / "poster" / "figures",
                    help="Output directory")
    args = ap.parse_args()

    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)

    _style()

    print("Generating beta poster figures …")
    fig7_pc_sensitivity(out)
    fig8_block_removal(out)
    fig9_unified_stability(out)
    print(f"\nDone → {out}/")


if __name__ == "__main__":
    main()
