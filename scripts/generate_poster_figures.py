#!/usr/bin/env python3
"""
generate_poster_figures.py – Poster figure set for GENO-MAP (Frontiers style).

Generates 5 publication-ready data figures:
  Fig 3  — Panel geometry diagnostics (dot-plot, 3 metrics × 4 panels)
  Fig 4  — Robustness curves HERO (marker subsampling + missing injection + SS inset)
  Fig 5  — UMAP diversity scatter (2×2 grid, 4 panels)
  Fig 6  — PCA vs AE quality–stability trade-off (grouped bars)
  Fig 10 — Stability frontier (trust crossover + stability gap)

Style: Frontiers in Plant Science / Genetics
  • Arial / DejaVu Sans, min 8 pt
  • 300 DPI, colour-blind-safe Okabe-Ito palette
  • White bg, left+bottom spines only
  • Panel labels: bold (A), (B), … top-left

Usage:
    python3 scripts/generate_poster_figures.py [--outdir docs/figures/poster]
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

# ── Frontiers house style ───────────────────────────────────────────
DPI = 300
SINGLE_COL_MM = 85
DOUBLE_COL_MM = 170
MM = 1 / 25.4  # mm → inch

# Okabe-Ito colour-blind-safe palette
C = {
    "global_snp":    "#0072B2",   # blue
    "global_silico": "#009E73",   # teal
    "ld_snp":        "#D55E00",   # vermilion
    "ld_silico":     "#CC79A7",   # rose
    "pca":           "#0072B2",
    "ae":            "#D55E00",
}

DATASET_ORDER = ["global_snp", "global_silico", "ld_snp", "ld_silico"]
DS_SHORT = ["G-SNP", "G-Silico", "LD-SNP", "LD-Silico"]
DS_LABEL = {
    "global_snp":    "Global SNP",
    "global_silico": "Global SilicoDArT",
    "ld_snp":        "LowDensity SNP",
    "ld_silico":     "LowDensity SilicoDArT",
}


def _style():
    """Set Frontiers-compatible rcParams."""
    plt.rcParams.update({
        "font.family":       "sans-serif",
        "font.sans-serif":   ["Arial", "DejaVu Sans"],
        "font.size":         9,
        "axes.titlesize":    10,
        "axes.labelsize":    10,
        "xtick.labelsize":   8,
        "ytick.labelsize":   8,
        "legend.fontsize":   8,
        "figure.dpi":        DPI,
        "savefig.dpi":       DPI,
        "savefig.bbox":      "tight",
        "savefig.pad_inches": 0.05,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "axes.linewidth":    0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size":  3,
        "ytick.major.size":  3,
        "lines.linewidth":   1.4,
        "lines.markersize":  5,
        "legend.frameon":    False,
        "figure.facecolor":  "white",
        "axes.facecolor":    "white",
    })


def _lbl(ax, s, x=-0.12, y=1.08):
    """Bold panel label (A), (B), …"""
    ax.text(x, y, s, transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top", ha="left")


# ════════════════════════════════════════════════════════════════════
# DATA (hardcoded from experiments/ to avoid runtime JSON dependency)
# ════════════════════════════════════════════════════════════════════

# -- Fig 3: geometry diagnostics --
DIAG = {  # (PC1%, eff_rank, reciprocity)
    "global_snp":    (9.49,  10.4, 0.648),
    "global_silico": (10.20, 10.0, 0.647),
    "ld_snp":        (8.87,  16.2, 0.699),
    "ld_silico":     (9.46,  16.1, 0.686),
}

# -- Fig 4: robustness --
MK_F = [0.05, 0.10, 0.20, 0.50, 0.80]
MK_J = {
    "global_snp":    [0.429, 0.511, 0.605, 0.745, 0.844],
    "global_silico": [0.563, 0.646, 0.727, 0.837, 0.901],
    "ld_snp":        [0.590, 0.665, 0.728, 0.826, 0.891],
    "ld_silico":     [0.527, 0.606, 0.680, 0.802, 0.884],
}
MK_SS = {
    "global_snp":    [0.917, 0.953, 0.981, 0.996, 0.999],
    "global_silico": [0.942, 0.979, 0.991, 0.998, 0.9995],
    "ld_snp":        [0.969, 0.988, 0.994, 0.998, 0.9996],
    "ld_silico":     [0.932, 0.967, 0.986, 0.997, 0.9994],
}
MS_R = [0.00, 0.05, 0.10, 0.20]
MS_J = {
    "global_snp":    [0.992, 0.876, 0.834, 0.772],
    "global_silico": [0.993, 0.925, 0.898, 0.858],
    "ld_snp":        [0.990, 0.895, 0.868, 0.828],
    "ld_silico":     [0.990, 0.887, 0.854, 0.808],
}

# -- Fig 6: PCA vs AE --
TRUST = {  # (PCA, AE-v1)
    "global_snp":    (0.976, 0.988),
    "global_silico": (0.979, 0.989),
    "ld_snp":        (0.933, 0.836),
    "ld_silico":     (0.911, 0.829),
}
STAB = {  # (PCA, AE-v1) edge Jaccard
    "global_snp":    (0.879, 0.524),
    "global_silico": (0.891, 0.498),
    "ld_snp":        (0.907, 0.163),
    "ld_silico":     (0.908, 0.207),
}

# -- Fig 10: frontier (global_snp) --
FR_N   = [100, 200, 500, 1000, 2000, 3000, 3500, 4000,
          4100, 4200, 4300, 4500, 5000, 5970]
FR_PT  = [.939, .932, .962, .971, .978, .982, .983, .984,
          .984, .984, .984, .985, .985, .986]
FR_AT  = [.886, .869, .875, .961, .927, .981, .983, .978,
          .984, .984, .983, .986, .987, .988]
FR_PS  = [.805, .861, .897, .985, .986, .987, .986, .988,
          .987, .987, .987, .990, .988, .987]
FR_AS  = [.394, .309, .208, .458, .173, .499, .500, .397,
          .503, .451, .468, .535, .532, .524]


# ════════════════════════════════════════════════════════════════════
# FIGURES
# ════════════════════════════════════════════════════════════════════

# QA flag thresholds (from panel_diagnostics.py)
THRESH = {
    "PC1 variance (%)": {"val": 50, "label": "1D-dominant", "side": "right"},
    "Effective rank":   {"val": 5,  "label": "Low-rank",    "side": "left"},
    "kNN reciprocity":  {"val": 0.40, "label": "Low-reciprocity", "side": "left"},
}


def fig3_geometry(out: Path):
    """Dot-plot: 3 geometry metrics × 4 panels + QA threshold lines.

    Poster improvements:
    - Distinct marker shapes per panel (consistent with fig4)
    - Larger markers with numeric value annotations
    - "Safe" zone shading (green) vs "alert" zone (red tint)
    - Horizontal connecting lines for visual grouping
    - Panel labels on every subplot for standalone readability
    """
    MARKERS = ["o", "s", "D", "^"]  # same as fig4
    MS = 10  # larger for poster distance

    metrics = ["PC1 variance (%)", "Effective rank", "kNN reciprocity"]
    data = np.array([DIAG[d] for d in DATASET_ORDER])

    fig, axes = plt.subplots(1, 3,
                             figsize=(DOUBLE_COL_MM * MM, 80 * MM))
    colors = [C[d] for d in DATASET_ORDER]

    for j, (ax, m) in enumerate(zip(axes, metrics)):
        y = np.arange(4)
        t = THRESH[m]

        # Safe / alert zone shading
        if t["side"] == "right":
            # Values BELOW threshold are safe (PC1 < 50%)
            ax.axvspan(ax.get_xlim()[0] if j == 0 else 0, t["val"],
                       color="#009E73", alpha=0.06, zorder=0)
            ax.axvspan(t["val"], 100,
                       color="#D55E00", alpha=0.06, zorder=0)
        else:
            # Values ABOVE threshold are safe (eff_rank > 5, recip > 0.40)
            ax.axvspan(t["val"], 100 if j == 1 else 1.0,
                       color="#009E73", alpha=0.06, zorder=0)
            ax.axvspan(0, t["val"],
                       color="#D55E00", alpha=0.06, zorder=0)

        # Horizontal connecting lines (subtle)
        for i in range(4):
            ax.hlines(y[i], 0, data[i, j], color=colors[i],
                      alpha=0.2, lw=1.5, zorder=1)

        # Markers
        for i in range(4):
            ax.plot(data[i, j], y[i], marker=MARKERS[i], color=colors[i],
                    markersize=MS, zorder=3,
                    markeredgecolor="white", markeredgewidth=0.7)
            # Value annotation
            fmt = f"{data[i,j]:.1f}" if j < 2 else f"{data[i,j]:.3f}"
            ax.text(data[i, j], y[i] - 0.30, fmt,
                    fontsize=6, ha="center", va="bottom", color="#333",
                    fontweight="bold")

        ax.set_yticks(y)
        ax.set_yticklabels(DS_SHORT if j == 0 else [])
        ax.set_xlabel(m)
        ax.invert_yaxis()
        ax.grid(axis="x", lw=0.3, alpha=0.4)
        ax.set_axisbelow(True)
        _lbl(ax, chr(65 + j))

        # QA threshold reference line
        ax.axvline(t["val"], ls="--", lw=1.0, color="#CC3333", alpha=0.7,
                   zorder=2)
        txt_ha = "right" if t["side"] == "right" else "left"
        ax.text(t["val"], -0.6, f"  {t['label']}  ", fontsize=6.5,
                color="#CC3333", ha=txt_ha, va="center", style="italic")

    fig.tight_layout(w_pad=1.2)
    fig.savefig(out / "fig3_panel_geometry.png")
    fig.savefig(out / "fig3_panel_geometry.pdf")
    plt.close(fig)
    print(f"  ✓ fig3_panel_geometry")


def fig4_robustness(out: Path):
    """HERO: marker subsampling + missing injection + SS inset.

    Improvements for poster readability:
    - Distinct marker shapes per panel (visible from 1–2 m)
    - Semantic background zones (safe / caution / degraded)
    - Key threshold annotations
    - Heavier lines and larger markers
    """
    MARKERS = ["o", "s", "D", "^"]  # per panel
    MS = 6  # marker size for main curves

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(DOUBLE_COL_MM * MM, 90 * MM), sharey=True)

    # ── Semantic zones (both panels) ──
    for ax in (ax1, ax2):
        ax.axhspan(0.80, 1.02, color="#009E73", alpha=0.06, zorder=0)   # safe
        ax.axhspan(0.60, 0.80, color="#E69F00", alpha=0.06, zorder=0)   # caution
        ax.axhspan(0.00, 0.60, color="#D55E00", alpha=0.06, zorder=0)   # degraded
        ax.axhline(0.80, ls=":", lw=0.5, color="#999", zorder=0)
        ax.axhline(0.60, ls=":", lw=0.5, color="#999", zorder=0)

    # (A) Marker subsampling
    fracs = [f * 100 for f in MK_F]
    for i, ds in enumerate(DATASET_ORDER):
        ax1.plot(fracs, MK_J[ds],
                 marker=MARKERS[i], color=C[ds], label=DS_LABEL[ds],
                 ms=MS, lw=1.8, markeredgecolor="white", markeredgewidth=0.5)
    ax1.set_xlabel("Markers retained (%)")
    ax1.set_ylabel("Jaccard neighbor overlap ($J_{\\mathrm{nbr}}$)")
    ax1.set_ylim(0.35, 1.02)
    ax1.set_xlim(0, 85)
    ax1.xaxis.set_major_locator(mticker.MultipleLocator(20))
    ax1.grid(axis="y", lw=0.3, alpha=0.4)
    _lbl(ax1, "A")

    # Zone labels (right edge of panel A)
    ax1.text(83, 0.91, "safe", fontsize=6, color="#009E73",
             ha="right", va="center", fontstyle="italic", alpha=0.7)
    ax1.text(83, 0.70, "caution", fontsize=6, color="#E69F00",
             ha="right", va="center", fontstyle="italic", alpha=0.7)
    ax1.text(83, 0.48, "degraded", fontsize=6, color="#D55E00",
             ha="right", va="center", fontstyle="italic", alpha=0.7)

    # SS inset (repositioned to avoid data overlap)
    ins = ax1.inset_axes([0.55, 0.12, 0.42, 0.30])
    for i, ds in enumerate(DATASET_ORDER):
        ins.plot([f * 100 for f in MK_F], MK_SS[ds],
                 marker=MARKERS[i], color=C[ds], markersize=3, lw=1.0,
                 markeredgecolor="white", markeredgewidth=0.3)
    ins.set_ylabel("SS", fontsize=7, labelpad=2)
    ins.set_xlabel("Markers (%)", fontsize=7, labelpad=2)
    ins.set_ylim(0.90, 1.005)
    ins.tick_params(labelsize=6)
    ins.axhline(0.91, ls="--", lw=0.5, color="grey", alpha=0.6)
    ins.set_title("PCA subspace sim.", fontsize=7, pad=2)
    for sp in ["top", "right"]:
        ins.spines[sp].set_visible(False)
    # Annotation: "≥ 0.91 even at 5%"
    ins.annotate("≥ 0.91 at 5%", xy=(5, 0.917), xytext=(30, 0.925),
                 fontsize=5.5, color="#555", ha="center",
                 arrowprops=dict(arrowstyle="->", color="#999", lw=0.5))

    # (B) Missing injection
    rates = [r * 100 for r in MS_R]
    for i, ds in enumerate(DATASET_ORDER):
        ax2.plot(rates, MS_J[ds],
                 marker=MARKERS[i], color=C[ds], label=DS_LABEL[ds],
                 ms=MS, lw=1.8, markeredgecolor="white", markeredgewidth=0.5)
    ax2.set_xlabel("Additional MCAR missing (%)")
    ax2.set_xlim(-1, 22)
    ax2.xaxis.set_major_locator(mticker.MultipleLocator(5))
    ax2.grid(axis="y", lw=0.3, alpha=0.4)
    _lbl(ax2, "B")

    # Shared legend
    h, l = ax1.get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, -0.02), fontsize=8, frameon=False)

    fig.tight_layout(rect=[0, 0.07, 1, 1], w_pad=1.5)
    fig.savefig(out / "fig4_robustness_hero.png")
    fig.savefig(out / "fig4_robustness_hero.pdf")
    plt.close(fig)
    print(f"  ✓ fig4_robustness_hero")


def fig6_pca_vs_ae(out: Path):
    """Grouped bars: trust + stability, PCA vs AE, per-panel colours."""
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(DOUBLE_COL_MM * MM, 80 * MM))

    x = np.arange(4)
    w = 0.35
    panel_colors = [C[d] for d in DATASET_ORDER]
    # Lighter tint for AE bars (alpha blend with white)
    def _light(hex_c, t=0.45):
        r, g, b = int(hex_c[1:3], 16), int(hex_c[3:5], 16), int(hex_c[5:7], 16)
        r2 = int(r + (255 - r) * t)
        g2 = int(g + (255 - g) * t)
        b2 = int(b + (255 - b) * t)
        return f"#{r2:02x}{g2:02x}{b2:02x}"
    ae_colors = [_light(c) for c in panel_colors]

    # (A) Trustworthiness
    pt = [TRUST[d][0] for d in DATASET_ORDER]
    at = [TRUST[d][1] for d in DATASET_ORDER]
    pca_bars_a = ax1.bar(x - w / 2, pt, w, color=panel_colors,
                         edgecolor="white", linewidth=0.4, label="PCA-30D")
    ae_bars_a = ax1.bar(x + w / 2, at, w, color=ae_colors,
                        edgecolor=[c for c in panel_colors], linewidth=0.8,
                        hatch="//", label="AE-64D")
    ax1.set_ylabel("Trustworthiness ($k$=15)")
    ax1.set_xticks(x)
    ax1.set_xticklabels(DS_SHORT, rotation=30, ha="right")
    ax1.set_ylim(0, 1.05)
    ax1.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
    # Custom legend: solid = PCA, hatched = AE
    from matplotlib.patches import Patch
    ax1.legend(handles=[Patch(facecolor="#888", edgecolor="white", label="PCA-30D"),
                        Patch(facecolor="#ccc", edgecolor="#888", hatch="//", label="AE-64D")],
               loc="lower left", fontsize=8)
    _lbl(ax1, "A")

    # Bar-top numeric values
    for i in range(4):
        ax1.text(x[i] - w / 2, pt[i] + 0.004, f"{pt[i]:.3f}",
                 ha="center", va="bottom", fontsize=5.5, color="#333")
        ax1.text(x[i] + w / 2, at[i] + 0.004, f"{at[i]:.3f}",
                 ha="center", va="bottom", fontsize=5.5, color="#333")

    # (B) Stability
    ps = [STAB[d][0] for d in DATASET_ORDER]
    ae_s = [STAB[d][1] for d in DATASET_ORDER]
    ax2.bar(x - w / 2, ps, w, color=panel_colors,
            edgecolor="white", linewidth=0.4, label="PCA-30D")
    ax2.bar(x + w / 2, ae_s, w, color=ae_colors,
            edgecolor=[c for c in panel_colors], linewidth=0.8,
            hatch="//", label="AE-64D")
    ax2.set_ylabel("Edge Jaccard (inter-seed stability)")
    ax2.set_xticks(x)
    ax2.set_xticklabels(DS_SHORT, rotation=30, ha="right")
    ax2.set_ylim(0.0, 1.05)
    ax2.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
    ax2.legend(handles=[Patch(facecolor="#888", edgecolor="white", label="PCA-30D"),
                        Patch(facecolor="#ccc", edgecolor="#888", hatch="//", label="AE-64D")],
               loc="upper right", fontsize=8)
    _lbl(ax2, "B")

    # Bar-top numeric values (panel B)
    for i in range(4):
        ax2.text(x[i] - w / 2, ps[i] + 0.015, f"{ps[i]:.3f}",
                 ha="center", va="bottom", fontsize=5.5, color="#333")
        ax2.text(x[i] + w / 2, ae_s[i] + 0.015, f"{ae_s[i]:.3f}",
                 ha="center", va="bottom", fontsize=5.5, color="#333")

    # Gap annotations (LD panels) — offset further right
    for i in [2, 3]:
        gap = ps[i] - ae_s[i]
        mid = (ps[i] + ae_s[i]) / 2
        ax2.annotate("", xy=(x[i] + w / 2 + 0.05, ae_s[i] + 0.01),
                     xytext=(x[i] + w / 2 + 0.05, ps[i] - 0.01),
                     arrowprops=dict(arrowstyle="<->", color="#555", lw=0.8))
        ax2.text(x[i] + w / 2 + 0.18, mid,
                 f"Δ={gap:.2f}", fontsize=6, color="#555", va="center")

    fig.tight_layout(w_pad=2.0)
    fig.savefig(out / "fig6_pca_vs_ae.png")
    fig.savefig(out / "fig6_pca_vs_ae.pdf")
    plt.close(fig)
    print(f"  ✓ fig6_pca_vs_ae")


def fig10_frontier(out: Path):
    """Trust crossover + stability gap (global_snp subsample frontier)."""
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(DOUBLE_COL_MM * MM, 120 * MM), sharex=True)
    n = np.array(FR_N)

    # (A) Trustworthiness
    ax1.plot(n, FR_PT, "o-", color=C["pca"], label="PCA-30D", ms=4)
    ax1.plot(n, FR_AT, "s-", color=C["ae"],  label="AE-64D",  ms=4)
    ax1.axvspan(4000, 4500, alpha=0.10, color="#999", zorder=0)
    ax1.annotate("$n^*$ ≈ 4 100", xy=(4100, 0.984),
                 xytext=(4100, 1.005), fontsize=7.5,
                 ha="center", va="bottom", color="#555",
                 arrowprops=dict(arrowstyle="->", color="#999", lw=0.7))
    ax1.set_ylabel("Trustworthiness")
    ax1.set_ylim(0.85, 1.02)
    ax1.legend(loc="lower right", fontsize=8)
    _lbl(ax1, "A", x=-0.10)

    # (B) Stability
    ax2.plot(n, FR_PS, "o-", color=C["pca"], label="PCA-30D", ms=4)
    ax2.plot(n, FR_AS, "s-", color=C["ae"],  label="AE-64D",  ms=4)
    ax2.annotate("", xy=(5970, FR_AS[-1] + 0.01),
                 xytext=(5970, FR_PS[-1] - 0.01),
                 arrowprops=dict(arrowstyle="<->", color="#555", lw=0.8))
    ax2.text(5400, 0.75, "gap ≈ 0.46", fontsize=7.5, color="#555", ha="right")
    ax2.set_xlabel("Number of samples ($n$)")
    ax2.set_ylabel("Edge Jaccard stability")
    ax2.set_ylim(0.05, 1.05)
    ax2.legend(loc="center right", fontsize=8)
    _lbl(ax2, "B", x=-0.10)
    ax2.xaxis.set_major_locator(mticker.MultipleLocator(1000))
    ax2.set_xlim(0, 6200)

    fig.tight_layout(h_pad=0.8)
    fig.savefig(out / "fig10_stability_frontier.png")
    fig.savefig(out / "fig10_stability_frontier.pdf")
    plt.close(fig)
    print(f"  ✓ fig10_stability_frontier")


def fig_stability_regime(out: Path):
    """Stability frontier vs n/p — single-panel conceptual figure.

    Poster improvements (consistent with fig3/fig4/fig6):
    - Distinct marker shapes (o=PCA, s=AE) with larger size
    - Semantic zones: stable (green) / unstable (red)
    - Heavier lines, white marker edges
    - Shaded region between curves showing the persistent gap
    - Trust crossover marker with cleaner annotation
    """
    p = 20_069  # Global SNP marker count
    n = np.array(FR_N, dtype=float)
    np_ratio = n / p

    fig, ax = plt.subplots(figsize=(SINGLE_COL_MM * 1.1 * MM,
                                    SINGLE_COL_MM * 0.85 * MM))

    # ── Semantic zones ──
    ax.axhspan(0.80, 1.10, color="#009E73", alpha=0.06, zorder=0)   # stable
    ax.axhspan(0.00, 0.80, color="#D55E00", alpha=0.06, zorder=0)   # unstable
    ax.axhline(0.80, ls=":", lw=0.5, color="#999", zorder=0)
    ax.text(0.005, 0.83, "stable", fontsize=6, color="#009E73",
            fontstyle="italic", alpha=0.7, va="bottom")
    ax.text(0.005, 0.77, "unstable", fontsize=6, color="#D55E00",
            fontstyle="italic", alpha=0.7, va="top")

    # ── Shaded gap between PCA and AE ──
    ax.fill_between(np_ratio, FR_AS, FR_PS,
                    color="#D55E00", alpha=0.08, zorder=1,
                    label="Stability gap")

    # ── PCA: nearly flat high line ──
    ax.plot(np_ratio, FR_PS, "o-", color=C["pca"], label="PCA-30D",
            ms=5, lw=1.8, zorder=3,
            markeredgecolor="white", markeredgewidth=0.5)
    # ── AE: growing but always below ──
    ax.plot(np_ratio, FR_AS, "s-", color=C["ae"], label="AE-64D",
            ms=5, lw=1.8, zorder=3,
            markeredgecolor="white", markeredgewidth=0.5)

    # ── Trust crossover vertical line at n*/p ≈ 0.20 ──
    n_star_ratio = 4100 / p  # ≈ 0.204
    ax.axvline(n_star_ratio, ls="--", lw=1.0, color="#888", alpha=0.7,
               zorder=2)
    ax.annotate(f"$n^*/p \\approx {n_star_ratio:.2f}$\ntrust crossover",
                xy=(n_star_ratio, 0.30), xytext=(n_star_ratio + 0.04, 0.22),
                fontsize=7, color="#555", ha="left", va="center",
                arrowprops=dict(arrowstyle="->", color="#999", lw=0.7))

    # ── Stability gap annotation at max n/p ──
    max_np = np_ratio[-1]
    gap = FR_PS[-1] - FR_AS[-1]
    mid_y = (FR_PS[-1] + FR_AS[-1]) / 2
    ax.annotate("", xy=(max_np + 0.005, FR_AS[-1] + 0.01),
                xytext=(max_np + 0.005, FR_PS[-1] - 0.01),
                arrowprops=dict(arrowstyle="<->", color="#555", lw=0.9))
    ax.text(max_np - 0.01, mid_y,
            f"Δ ≈ {gap:.2f}", fontsize=7.5, color="#555",
            ha="right", va="center", fontweight="bold")

    ax.set_xlabel("$n / p$")
    ax.set_ylabel("Inter-seed kNN Edge Jaccard")
    ax.set_ylim(0.05, 1.08)
    ax.set_xlim(-0.01, 0.32)
    ax.legend(loc="center right", fontsize=7, framealpha=0.9)
    ax.xaxis.set_major_locator(mticker.MultipleLocator(0.05))
    ax.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
    ax.grid(axis="y", lw=0.3, alpha=0.4)

    fig.tight_layout()
    fig.savefig(out / "fig_stability_regime.png")
    fig.savefig(out / "fig_stability_regime.pdf")
    plt.close(fig)
    print(f"  ✓ fig_stability_regime")


# ════════════════════════════════════════════════════════════════════

# Node JSON paths (poster-v2 run, seed 42)
NODE_FILES = {
    "global_snp":    "experiments/global_snp/poster-v2/seed42_global_snp_nodes.json",
    "global_silico": "experiments/global_silico/poster-v2/seed42_global_silico_nodes.json",
    "ld_snp":        "experiments/lowdensity_snp/poster-v2/seed42_lowdensity_snp_nodes.json",
    "ld_silico":     "experiments/lowdensity_silico/poster-v2/seed42_lowdensity_silico_nodes.json",
}


def fig5_umap_scatter(out: Path):
    """2×2 UMAP diversity scatter — one panel per dataset."""
    fig, axes = plt.subplots(
        2, 2, figsize=(DOUBLE_COL_MM * MM, 130 * MM))
    axes = axes.ravel()

    panel_colors = [C[d] for d in DATASET_ORDER]
    n_samples = {"global_snp": 5970, "global_silico": 5970,
                 "ld_snp": 630, "ld_silico": 635}
    panel_titles = {
        "global_snp":    "Global SNP",
        "global_silico": "Global SilicoDArT",
        "ld_snp":        "LowDensity SNP",
        "ld_silico":     "LowDensity SilicoDArT",
    }

    for i, ds in enumerate(DATASET_ORDER):
        ax = axes[i]
        fpath = Path(NODE_FILES[ds])
        if not fpath.exists():
            ax.text(0.5, 0.5, "[data not found]",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=10, color="#999")
            _lbl(ax, chr(65 + i))
            continue

        nodes = json.load(open(fpath))
        coords = np.array([n["embedding"] for n in nodes])
        x, y = coords[:, 0], coords[:, 1]

        ax.scatter(x, y, s=1.5, alpha=0.45, color=panel_colors[i],
                   edgecolors="none", rasterized=True)

        # Remove numeric ticks (UMAP coords are not interpretable)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("UMAP-1", fontsize=8)
        ax.set_ylabel("UMAP-2", fontsize=8)

        # Panel info
        title = f"{panel_titles[ds]}  ($n$ = {n_samples[ds]:,})"
        ax.set_title(title, fontsize=9, pad=4)
        _lbl(ax, chr(65 + i))

    fig.tight_layout(w_pad=1.5, h_pad=2.0)
    fig.savefig(out / "fig5_umap_scatter.png")
    fig.savefig(out / "fig5_umap_scatter.pdf")
    plt.close(fig)
    print(f"  ✓ fig5_umap_scatter")


def fig5_umap_single(out: Path):
    """Single-panel UMAP scatter (Global SNP) — A1 poster inset version."""
    fpath = Path(NODE_FILES["global_snp"])
    if not fpath.exists():
        print("  ⚠ fig5_umap_single skipped (data not found)")
        return

    nodes = json.load(open(fpath))
    coords = np.array([n["embedding"] for n in nodes])
    x, y = coords[:, 0], coords[:, 1]

    fig, ax = plt.subplots(
        figsize=(SINGLE_COL_MM * MM, SINGLE_COL_MM * MM))

    ax.scatter(x, y, s=1.2, alpha=0.40, color=C["global_snp"],
               edgecolors="none", rasterized=True)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP-1", fontsize=8)
    ax.set_ylabel("UMAP-2", fontsize=8)
    ax.set_title("Global SNP  ($n$ = 5,970)", fontsize=9, pad=4)

    fig.tight_layout()
    fig.savefig(out / "fig5_umap_single.png")
    fig.savefig(out / "fig5_umap_single.pdf")
    plt.close(fig)
    print(f"  ✓ fig5_umap_single")


def main():
    ap = argparse.ArgumentParser(
        description="Generate poster figures (Frontiers style)")
    ap.add_argument("--outdir", default="docs/figures/poster",
                    help="Output directory (default: docs/figures/poster)")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    _style()

    print("Generating poster figures (Frontiers style) …")
    fig3_geometry(out)
    fig4_robustness(out)
    fig5_umap_scatter(out)
    fig5_umap_single(out)
    fig6_pca_vs_ae(out)
    fig10_frontier(out)
    fig_stability_regime(out)
    print(f"\nDone → {out}/")


if __name__ == "__main__":
    main()
