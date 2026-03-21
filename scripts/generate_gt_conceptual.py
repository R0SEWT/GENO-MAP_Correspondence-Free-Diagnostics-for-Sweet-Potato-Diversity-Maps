#!/usr/bin/env python3
"""
Ground-truth validation figure (ADR-0016) — v3.

Three columns comparing the 3 conditions side-by-side:
  1. SHUFFLED  (null) — random identity mapping  → J ≈ 0
  2. GROUND TRUTH     — correct identity         → J ≈ 0.59
  3. INTRA-PANEL      — same data, diff seed     → J ≈ 0.92

Each column: two mini kNN-graphs with correspondence lines.
Bottom: ruler placing the 3 values on a 0→1 scale.
"""

import pathlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Palette ───────────────────────────────────────────────
SLATE = "#4C5C68"
STONE = "#8A817C"
WINE  = "#A44A3F"
OLIVE = "#5E7C4D"
AMBER = "#B98918"
GOLD  = "#D8A62A"
GREY  = "#777777"
DARK  = "#222222"
BG_NO = "#F7ECE9"
BG_GT = "#EEF4E8"
BG_IN = "#F5F0E6"

OUT = pathlib.Path(__file__).resolve().parent.parent / "docs" / "poster" / "figures"

N = 8


def circle_pos(n, cx, cy, r=0.20):
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False) - np.pi / 2
    return [(cx + r * np.cos(a), cy + r * np.sin(a)) for a in angles]


# ── Edge sets ─────────────────────────────────────────────
# Base ring + diagonals = 12 edges each
ring = [(i, (i + 1) % N) for i in range(N)]           # 8 ring edges

diag_shared = [(0, 3), (1, 5), (2, 6), (4, 7)]        # 4 shared diags
diag_snp_only = [(0, 4), (2, 5), (3, 7), (1, 6)]      # 4 SNP-only diags
diag_sil_only = [(0, 5), (1, 4), (3, 6), (2, 7)]      # 4 SilicoDArT-only diags

edges_snp = ring + diag_shared + diag_snp_only       # 16
edges_sil = ring + diag_shared + diag_sil_only       # 16

# For INTRA: same panel, different seed → 92 % overlap
# We reuse edges_snp but drop 1 edge and add 1 different → ~94 %
edges_snp_seed2 = ring + diag_shared + diag_snp_only[:-1] + [(1, 3)]  # 15 same + 1 diff

shared_gt_set = set(tuple(sorted(e)) for e in ring + diag_shared)  # 12 shared across techs
shared_intra_set = set(tuple(sorted(e)) for e in ring + diag_shared + diag_snp_only[:-1])

identity_map = list(range(N))
shuffled_map = [5, 3, 7, 1, 6, 0, 4, 2]


def draw_mini_graph(ax, pos, edges, node_color, shared_set=None,
                     highlight=False, edge_alpha=0.25):
    """Draw a compact mini-graph."""
    for i, j in edges:
        key = tuple(sorted((i, j)))
        is_shared = shared_set and key in shared_set
        if highlight and is_shared:
            ax.plot([pos[i][0], pos[j][0]], [pos[i][1], pos[j][1]],
                    color=GOLD, lw=5, alpha=0.22, zorder=0,
                    solid_capstyle="round")
        lw_ = 1.4 if (highlight and is_shared) else 1.0
        al_ = 0.50 if (highlight and is_shared) else edge_alpha
        ax.plot([pos[i][0], pos[j][0]], [pos[i][1], pos[j][1]],
                color=node_color, lw=lw_, alpha=al_, zorder=1)
    for idx, (x, y) in enumerate(pos):
        ax.scatter(x, y, s=130, color=node_color, edgecolors="white",
                   linewidths=1.2, zorder=5)
        ax.text(x, y, str(idx + 1), ha="center", va="center",
                fontsize=5.5, color="white", fontweight="bold", zorder=6)


def draw_corr(ax, pos_a, pos_b, mapping, color, style="-", lw=0.8, alpha=0.25):
    for i, j in enumerate(mapping):
        ax.plot([pos_a[i][0], pos_b[j][0]],
                [pos_a[i][1], pos_b[j][1]],
                color=color, ls=style, lw=lw, alpha=alpha, zorder=2)


# ══════════════════════════════════════════════════════════
#  FIGURE — 3 columns + bottom ruler
# ══════════════════════════════════════════════════════════

fig = plt.figure(figsize=(9.5, 5.2))

columns = [
    # (title, subtitle, edges_top, edges_bot, color_top, color_bot,
    #  mapping, shared_set, bg, accent, j_val, highlight)
    {
        "title": "Shuffled (null)",
        "sub": "identities scrambled\n→ no edge overlap",
        "e_top": edges_snp, "e_bot": edges_sil,
        "c_top": SLATE, "c_bot": STONE,
        "mapping": shuffled_map,
        "shared": None, "bg": BG_NO, "accent": WINE,
        "j_val": 0.001, "highlight": False,
        "label_top": "Panel A", "label_bot": "Panel B",
    },
    {
        "title": "Ground truth",
        "sub": "correct identity\n→ 59 % edge overlap",
        "e_top": edges_snp, "e_bot": edges_sil,
        "c_top": SLATE, "c_bot": STONE,
        "mapping": identity_map,
        "shared": shared_gt_set, "bg": BG_GT, "accent": OLIVE,
        "j_val": 0.593, "highlight": True,
        "label_top": "Panel A", "label_bot": "Panel B",
    },
    {
        "title": "Intra-panel (ceiling)",
        "sub": "same data, different seed\n→ 92 % edge overlap",
        "e_top": edges_snp, "e_bot": edges_snp_seed2,
        "c_top": SLATE, "c_bot": STONE,
        "mapping": identity_map,
        "shared": shared_intra_set, "bg": BG_IN, "accent": AMBER,
        "j_val": 0.917, "highlight": True,
        "label_top": "Panel A  seed 42", "label_bot": "Panel A  seed 52",
    },
]

for ci, col in enumerate(columns):
    x0 = 0.02 + ci * 0.33
    ax = fig.add_axes([x0, 0.24, 0.30, 0.66])
    ax.set_xlim(-0.08, 1.08)
    ax.set_ylim(-0.22, 1.22)

    # Background card
    rect = mpatches.FancyBboxPatch(
        (-0.06, -0.18), 1.12, 1.36,
        boxstyle="round,pad=0.03", fc=col["bg"], ec="none",
        alpha=0.35, zorder=-1)
    ax.add_patch(rect)

    cx = 0.5
    pos_top = circle_pos(N, cx, 0.80, r=0.20)
    pos_bot = circle_pos(N, cx, 0.25, r=0.20)

    # Graphs
    draw_mini_graph(ax, pos_top, col["e_top"], col["c_top"],
                     col["shared"], col["highlight"])
    draw_mini_graph(ax, pos_bot, col["e_bot"], col["c_bot"],
                     col["shared"], col["highlight"])

    # Correspondence
    corr_style = "--" if ci == 0 else "-"
    draw_corr(ax, pos_bot, pos_top, col["mapping"],
              col["accent"], corr_style,
              lw=0.7 if ci == 0 else 0.9,
              alpha=0.18 if ci == 0 else 0.25)

    # Graph labels
    ax.text(cx, 1.07, col["label_top"], ha="center", fontsize=8,
            fontweight="bold", color=col["c_top"])
    ax.text(cx, -0.02, col["label_bot"], ha="center", fontsize=8,
            fontweight="bold", color=col["c_bot"])

    # Subtitle between graphs
    ax.text(cx, 0.52, col["sub"], ha="center", fontsize=7,
            color=col["accent"], style="italic", alpha=0.85,
            linespacing=1.3)

    # Title
    ax.set_title(col["title"], fontsize=10.5, fontweight="bold",
                 color=col["accent"], pad=6)

    # J_edge badge at bottom
    j_txt = f"$J_{{edge}}$ = {col['j_val']:.3f}"
    ax.text(cx, -0.15, j_txt, ha="center", fontsize=11,
            fontweight="bold", color=col["accent"],
            bbox=dict(fc="white", ec=col["accent"], lw=1.3,
                      boxstyle="round,pad=0.35", alpha=0.92))

    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])

# ── Gold highlight legend ─────────────────────────────────
fig.patches.append(mpatches.FancyBboxPatch(
    (0.35, 0.895), 0.30, 0.03,
    boxstyle="round,pad=0.005", fc="white", ec=GREY,
    alpha=0.8, transform=fig.transFigure, figure=fig, lw=0.5))
# line + text drawn via fig.text
fig.text(0.40, 0.91, "━━", fontsize=9, color=GOLD, fontweight="bold",
         ha="center", va="center")
fig.text(0.43, 0.91, "= shared edge", fontsize=7.5, color=GREY,
         ha="left", va="center")

# ── Bottom ruler ─────────────────────────────────────────
ax_r = fig.add_axes([0.08, 0.04, 0.84, 0.14])
ax_r.set_xlim(-0.03, 1.03)
ax_r.set_ylim(-0.5, 0.9)

# Gradient background zones
ax_r.axvspan(-0.02, 0.05, color=WINE, alpha=0.06)
ax_r.axvspan(0.85, 1.02, color=AMBER, alpha=0.06)

# Baseline
ax_r.plot([0, 1], [0, 0], color=GREY, lw=2.5, zorder=1, solid_capstyle="round")
for v in np.arange(0, 1.01, 0.1):
    ax_r.plot([v, v], [-0.06, 0.06], color=GREY, lw=0.8, zorder=1)

# Label
ax_r.text(0.5, -0.40, "$J_{edge}$  (edge overlap)",
          ha="center", fontsize=9, color=DARK)
ax_r.text(0.0, -0.22, "0", ha="center", fontsize=7, color=GREY)
ax_r.text(1.0, -0.22, "1", ha="center", fontsize=7, color=GREY)

# Markers
pts = [
    (0.001, WINE,  "X",  "Shuffled",     0.55),
    (0.593, OLIVE, "o",  "Ground truth", 0.60),
    (0.917, AMBER, "D",  "Intra-panel",  0.55),
]
for val, color, marker, label, y_top in pts:
    ax_r.scatter([val], [0], marker=marker, s=160, color=color,
                 edgecolors="white", linewidths=1.3, zorder=5)
    ax_r.plot([val, val], [0.08, y_top - 0.15], color=color,
              lw=0.9, alpha=0.5, zorder=3)
    ax_r.text(val, y_top, f"{label}\n{val:.3f}", ha="center", va="center",
              fontsize=7.5, color=color, fontweight="bold",
              bbox=dict(fc="white", ec=color, alpha=0.85,
                        boxstyle="round,pad=0.2", lw=0.8))

# Zone labels
ax_r.text(0.025, -0.35, "noise", ha="center", fontsize=6, color=WINE, alpha=0.7)
ax_r.text(0.935, -0.35, "ceiling", ha="center", fontsize=6, color=AMBER, alpha=0.7)

for sp in ax_r.spines.values():
    sp.set_visible(False)
ax_r.set_xticks([])
ax_r.set_yticks([])


# ── Title ─────────────────────────────────────────────────
fig.suptitle("Ground-truth validation",
             fontsize=14, fontweight="bold", y=0.98, color=DARK)

# ── Save ──────────────────────────────────────────────────
OUT.mkdir(parents=True, exist_ok=True)
for ext in ["pdf", "png"]:
    fig.savefig(OUT / f"fig_gt_conceptual.{ext}", dpi=300, bbox_inches="tight")
plt.close()
print(f"Done → {OUT}/fig_gt_conceptual.pdf|png")
