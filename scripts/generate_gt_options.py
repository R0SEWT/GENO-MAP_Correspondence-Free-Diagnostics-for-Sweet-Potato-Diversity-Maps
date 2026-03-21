#!/usr/bin/env python3
"""
Generate 3 candidate figures for ADR-0016 ground-truth validation.
The user will compare and pick the best one for the poster.
"""

import json
import pathlib

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Okabe-Ito palette
BLUE   = "#0072B2"
TEAL   = "#009E73"
VERM   = "#D55E00"
GOLD   = "#E69F00"
GREY   = "#6F6F6F"
LIGHT  = "#F4F4F4"
MID    = "#D9D9D9"
DARK   = "#222222"
ROSE   = "#CC79A7"

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "experiments" / "adr" / "adr0016_ground_truth.json"
OUT  = ROOT / "docs" / "poster" / "figures"


def load_results():
    with open(DATA) as f:
        return json.load(f)["results"]


def _avg(seeds, key):
    vals = [s[key] for s in seeds if key in s]
    return float(np.mean(vals)) if vals else float("nan")


# ══════════════════════════════════════════════════════════
#  Option A: Grouped vertical bars (log inset)
# ══════════════════════════════════════════════════════════

def option_a(out):
    """Grouped vertical bars with log-scale inset."""
    results = load_results()
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.5), sharey=False)
    fig.subplots_adjust(wspace=0.35)

    for ax_idx, (pair_key, title) in enumerate([
        ("global", "Global (n = 5 930)"),
        ("lowdensity", "LowDensity (n = 630)")
    ]):
        ax = axes[ax_idx]
        seeds = results[pair_key]["seeds"]

        vals = [
            _avg(seeds, "J_edge_intra"),
            _avg(seeds, "J_edge_gt"),
            _avg(seeds, "J_edge_cf_shuffled"),
        ]
        labels = ["INTRA\n(same tech)", "GT\n(aligned)", "NULL\n(shuffled)"]
        colors = [BLUE, TEAL, VERM]

        bars = ax.bar(range(3), vals, color=colors, alpha=0.8,
                      edgecolor="white", linewidth=1.2, width=0.65)

        for bar, v in zip(bars, vals):
            if v > 0.05:
                ax.text(bar.get_x() + bar.get_width()/2, v + 0.02,
                        f"{v:.3f}", ha="center", fontsize=9, fontweight="bold")
            else:
                ax.text(bar.get_x() + bar.get_width()/2, v + 0.02,
                        f"{v:.4f}", ha="center", fontsize=8, color=VERM,
                        fontweight="bold")

        # ×450 annotation
        factor = vals[1] / vals[2] if vals[2] > 0 else 0
        ax.annotate(
            f"×{factor:.0f}", xy=(2, vals[2] + 0.01),
            xytext=(1.5, vals[1] * 0.5),
            fontsize=12, fontweight="bold", color=GOLD,
            arrowprops=dict(arrowstyle="->", color=GOLD, lw=1.5),
            ha="center",
        )

        ax.set_xticks(range(3))
        ax.set_xticklabels(labels, fontsize=7.5)
        ax.set_ylim(0, 1.08)
        ax.set_ylabel("$J_{edge}$", fontsize=10)
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", alpha=0.15)

    fig.suptitle("Ground-truth validation: $J_{edge}$ across conditions",
                 fontsize=11, fontweight="bold", y=1.01)
    plt.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  A → {out}")


# ══════════════════════════════════════════════════════════
#  Option B: Number line / ruler with regions
# ══════════════════════════════════════════════════════════

def option_b(out):
    """Ruler-style: 0→1 scale with positioned markers and regions."""
    results = load_results()

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 3.0), sharex=True)
    fig.subplots_adjust(hspace=0.55)

    for ax_idx, (pair_key, title) in enumerate([
        ("global", "Global collection (n = 5 930)"),
        ("lowdensity", "LowDensity collection (n = 630)")
    ]):
        ax = axes[ax_idx]
        seeds = results[pair_key]["seeds"]

        intra = _avg(seeds, "J_edge_intra")
        gt    = _avg(seeds, "J_edge_gt")
        shuf  = _avg(seeds, "J_edge_cf_shuffled")

        # Background regions
        ax.axhspan(-0.5, 0.5, xmin=0, xmax=shuf + 0.015,
                   color=VERM, alpha=0.07)
        ax.axhspan(-0.5, 0.5, xmin=gt - 0.02, xmax=gt + 0.02,
                   color=TEAL, alpha=0.07)
        ax.axhspan(-0.5, 0.5, xmin=intra - 0.015, xmax=1.0,
                   color=BLUE, alpha=0.07)

        # Horizontal baseline
        ax.plot([0, 1], [0, 0], color=GREY, lw=1.5, zorder=1)
        ax.set_xlim(-0.02, 1.05)
        ax.set_ylim(-0.35, 0.45)

        # Markers
        for val, color, marker, label, yoff in [
            (shuf,  VERM, "X",  f"Shuffled\n{shuf:.4f}", -0.22),
            (gt,    TEAL, "o",  f"Ground truth\n{gt:.3f}", 0.22),
            (intra, BLUE, "s",  f"Intra-panel\n{intra:.3f}", 0.22),
        ]:
            ax.scatter([val], [0], marker=marker, s=120, color=color,
                       edgecolors="white", linewidths=1, zorder=5)
            ax.plot([val, val], [0, yoff * 0.6], color=color, lw=1,
                    alpha=0.5, zorder=3)
            ax.text(val, yoff, label, ha="center", va="center",
                    fontsize=7, color=color, fontweight="bold",
                    bbox=dict(fc="white", ec=color, alpha=0.8,
                              boxstyle="round,pad=0.25", lw=0.8))

        # ×N annotation
        factor = gt / shuf if shuf > 0 else 0
        mid_x = (shuf + gt) / 2
        ax.annotate(
            f"×{factor:.0f}", xy=(mid_x, -0.02),
            fontsize=13, fontweight="bold", color=GOLD,
            ha="center", va="top",
            bbox=dict(fc=GOLD, ec="none", alpha=0.12,
                      boxstyle="round,pad=0.3"),
        )

        ax.set_title(title, fontsize=9.5, fontweight="bold", loc="left")
        ax.set_xlabel("$J_{edge}$", fontsize=10)
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

    fig.suptitle("Ground-truth validation",
                 fontsize=11, fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  B → {out}")


# ══════════════════════════════════════════════════════════
#  Option C: Heatmap (metrics × conditions)
# ══════════════════════════════════════════════════════════

def option_c(out):
    """Heatmap: rows = metrics, columns = conditions × collections."""
    results = load_results()

    metrics = [
        ("$J_{edge}$", "J_edge"),
        ("$J_{nbr}$",  "J_nbr"),
        ("SS",         "SS"),
    ]
    conditions = [
        ("Shuffled",  "cf_shuffled"),
        ("GT",        "gt"),
        ("Intra",     "intra"),
    ]
    collections = [
        ("Global", "global"),
        ("LD",     "lowdensity"),
    ]

    # Build matrix: rows = metrics, cols = cond × collection
    n_rows = len(metrics)
    n_cols = len(conditions) * len(collections)
    mat = np.full((n_rows, n_cols), np.nan)
    col_labels = []

    for ci, (coll_label, coll_key) in enumerate(collections):
        seeds = results[coll_key]["seeds"]
        for cj, (cond_label, cond_suffix) in enumerate(conditions):
            col_idx = ci * len(conditions) + cj
            col_labels.append(f"{cond_label}\n({coll_label})")
            for ri, (_, m_key) in enumerate(metrics):
                key = f"{m_key}_{cond_suffix}"
                mat[ri, col_idx] = _avg(seeds, key)

    fig, ax = plt.subplots(figsize=(7.5, 2.8))

    # Custom colormap: red(0) → white(0.5) → blue(1)
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(
        "gt_cmap", [(0, VERM), (0.15, "#FDEADC"), (0.5, "#F7F7F7"),
                     (0.7, "#D4E8F0"), (1.0, BLUE)]
    )

    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1)

    # Annotate cells
    for i in range(n_rows):
        for j in range(n_cols):
            v = mat[i, j]
            if np.isnan(v):
                continue
            txt = f"{v:.4f}" if v < 0.05 else f"{v:.3f}"
            color = "white" if v > 0.8 or v < 0.02 else DARK
            ax.text(j, i, txt, ha="center", va="center",
                    fontsize=8.5, fontweight="bold", color=color)

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(col_labels, fontsize=7.5)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([m[0] for m in metrics], fontsize=10)

    # Separator between collections
    ax.axvline(x=2.5, color=DARK, lw=1.5, alpha=0.5)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label("Score", fontsize=9)

    ax.set_title("Ground-truth validation: metrics × conditions",
                 fontsize=11, fontweight="bold", pad=10)

    plt.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  C → {out}")


# ══════════════════════════════════════════════════════════

if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    option_a(OUT / "fig_gt_option_a_bars.pdf")
    option_b(OUT / "fig_gt_option_b_ruler.pdf")
    option_c(OUT / "fig_gt_option_c_heatmap.pdf")
    print("\nDone — 3 options generated. Compare them side by side.")
