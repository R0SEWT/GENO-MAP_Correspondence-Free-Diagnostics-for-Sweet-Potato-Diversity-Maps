#!/usr/bin/env python3
"""
Generate ADR-0016 ground-truth validation figure for the poster.

Visual argument:
  - INTRA (same technology, diff seed) → upper bound ≈ 0.92
  - GT (SNP ↔ Silico, aligned by accession) → real signal ≈ 0.60–0.70
  - SHUFFLED (random alignment) → null ≈ 0.001–0.01

The ×450 gap between GT and SHUFFLED proves the method detects real structure.
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
DARK   = "#222222"

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "experiments" / "adr" / "adr0016_ground_truth.json"


def load_results():
    with open(DATA) as f:
        return json.load(f)["results"]


def fig_ground_truth(out):
    results = load_results()

    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.2), sharey=True)
    fig.subplots_adjust(wspace=0.08)

    conditions = ["INTRA\n(same tech,\ndiff seed)", "GT\n(aligned by\naccession)", "SHUFFLED\n(random\npermutation)"]
    colors = [BLUE, TEAL, VERM]
    markers = ["s", "o", "X"]

    for ax_idx, (pair_key, pair_label) in enumerate([
        ("global", "Global collection (n = 5 930)"),
        ("lowdensity", "LowDensity collection (n = 630)")
    ]):
        ax = axes[ax_idx]
        pair = results[pair_key]
        seeds = pair["seeds"]

        # Extract J_edge values per condition across seeds
        intra_vals  = [s["J_edge_intra"] for s in seeds]
        gt_vals     = [s["J_edge_gt"] for s in seeds]
        shuf_vals   = [s["J_edge_cf_shuffled"] for s in seeds]

        all_vals = [intra_vals, gt_vals, shuf_vals]

        for i, (vals, color, marker) in enumerate(zip(all_vals, colors, markers)):
            mean_v = np.mean(vals)
            # Individual seed dots
            jitter = np.linspace(-0.08, 0.08, len(vals))
            ax.scatter([i + j for j in jitter], vals,
                       color=color, marker=marker, s=60, zorder=5,
                       edgecolors="white", linewidths=0.5, alpha=0.85)
            # Mean bar
            ax.barh(i, mean_v, height=0.45, color=color, alpha=0.15, zorder=1)
            # Mean label
            if mean_v > 0.05:
                ax.text(mean_v + 0.02, i, f"{mean_v:.3f}",
                        va="center", ha="left", fontsize=8.5, color=DARK,
                        fontweight="bold")
            else:
                ax.text(mean_v + 0.02, i, f"{mean_v:.4f}",
                        va="center", ha="left", fontsize=8, color=VERM,
                        fontweight="bold")

        # Annotation: multiplier between GT and SHUFFLED
        gt_mean = np.mean(gt_vals)
        shuf_mean = np.mean(shuf_vals)
        factor = gt_mean / shuf_mean if shuf_mean > 0 else float("inf")

        ax.annotate(
            f"×{factor:.0f}",
            xy=(shuf_mean, 2), xytext=(gt_mean * 0.55, 1.55),
            fontsize=11, fontweight="bold", color=GOLD,
            ha="center",
            arrowprops=dict(arrowstyle="-", color=GOLD, lw=1.2, ls="--"),
        )

        ax.set_yticks(range(3))
        ax.set_yticklabels(conditions, fontsize=7.5)
        ax.set_xlim(-0.02, 1.05)
        ax.set_xlabel("$J_{edge}$", fontsize=10)
        ax.set_title(pair_label, fontsize=9.5, fontweight="bold", pad=8)
        ax.axvline(x=0, color=GREY, lw=0.5, alpha=0.3)

        # Light grid
        ax.grid(axis="x", alpha=0.15, lw=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Shared super-title
    fig.suptitle(
        "Ground-truth validation: kNN graphs from different marker technologies\n"
        "capture the same diversity structure",
        fontsize=10, fontweight="bold", y=1.02
    )

    plt.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  → {out}")


if __name__ == "__main__":
    out_dir = ROOT / "docs" / "poster" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_ground_truth(out_dir / "fig_ground_truth_validation.pdf")
    print("Done.")
