#!/usr/bin/env python3
"""Generate real ADR-0014 UMAP comparison figures from the analytic pipeline.

Outputs
-------
  fig5_umap_layouts_global.{png,pdf}
      Representative triptych for Global SNP:
      baseline / local-stress / global-stress.

  fig5_umap_layouts_all.{png,pdf}
      Full 4x3 grid across all panels and representative layouts.

  fig5_umap_layouts_summary.json
      Metrics and parameters used for the plotted layouts.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors

import warnings

warnings.filterwarnings("ignore", message=".*n_jobs.*overridden.*")
warnings.filterwarnings("ignore", message=".*Graph is not fully connected.*")

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_adr_experiments import (  # noqa: E402
    ADR0014_LAYOUTS,
    DATASETS,
    K,
    _impute,
    _load_panel,
    _pca_knn,
    jaccard_edges,
    jaccard_neighbors,
)


MM = 1 / 25.4
DPI = 300
DATASET_ORDER = [
    "global_snp",
    "global_silico",
    "lowdensity_snp",
    "lowdensity_silico",
]
DATASET_COLORS = {
    "global_snp": "#0072B2",
    "global_silico": "#009E73",
    "lowdensity_snp": "#D55E00",
    "lowdensity_silico": "#CC79A7",
}
LAYOUT_LABELS = {
    "baseline": "Baseline",
    "local_stress": "Local stress",
    "global_stress": "Global stress",
}
PARTITION_COLORS = ["#48B89F", "#D89AC6"]


def _style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.04,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.6,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })


def _compute_layouts(ds_key: str) -> Dict[str, Any]:
    import umap

    X_num = _load_panel(ds_key)
    X_imp = _impute(X_num)
    baseline = _pca_knn(X_imp, n_pca=30, seed=42)
    X_pca30 = baseline["X_pca"][:, :min(30, baseline["X_pca"].shape[1])]
    partition = KMeans(n_clusters=2, random_state=42, n_init=20).fit_predict(X_pca30)
    layouts = []

    for cfg in ADR0014_LAYOUTS:
        reducer = umap.UMAP(
            n_components=2,
            n_neighbors=cfg["n_neighbors"],
            min_dist=cfg["min_dist"],
            random_state=cfg["seed"],
            metric="cosine",
        )
        X_umap = reducer.fit_transform(X_pca30)

        k_eff = min(K, X_umap.shape[0] - 1)
        nbrs = NearestNeighbors(n_neighbors=k_eff + 1, metric="euclidean")
        nbrs.fit(X_umap)
        _, inds_umap = nbrs.kneighbors(X_umap)
        inds_umap = inds_umap[:, 1:]

        layouts.append({
            "name": cfg["name"],
            "n_neighbors": cfg["n_neighbors"],
            "min_dist": cfg["min_dist"],
            "seed": cfg["seed"],
            "J_edge_pca_vs_umap": round(
                jaccard_edges(baseline["knn_inds"], inds_umap), 4
            ),
            "J_nbr_pca_vs_umap": round(
                jaccard_neighbors(baseline["knn_inds"], inds_umap), 4
            ),
            "J_edge_pca_invariance": 1.0,
            "coords": X_umap,
        })

    return {
        "dataset": ds_key,
        "label": DATASETS[ds_key]["label"],
        "n": int(X_imp.shape[0]),
        "p": int(X_imp.shape[1]),
        "partition": partition,
        "layouts": layouts,
    }


def _panel_limits(layouts: list[Dict[str, Any]]) -> tuple[float, float, float, float]:
    xs = np.concatenate([l["coords"][:, 0] for l in layouts])
    ys = np.concatenate([l["coords"][:, 1] for l in layouts])
    dx = max(xs.max() - xs.min(), 1e-6)
    dy = max(ys.max() - ys.min(), 1e-6)
    pad_x = dx * 0.05
    pad_y = dy * 0.05
    return xs.min() - pad_x, xs.max() + pad_x, ys.min() - pad_y, ys.max() + pad_y


def _plot_layout(ax, coords: np.ndarray, partition: np.ndarray, title: str, metric_text: str,
                 xlim: tuple[float, float], ylim: tuple[float, float], point_size: float) -> None:
    for idx, color in enumerate(PARTITION_COLORS):
        mask = partition == idx
        if np.any(mask):
            ax.scatter(
                coords[mask, 0], coords[mask, 1],
                s=point_size, alpha=0.30, color=color,
                edgecolors="none", rasterized=True,
            )
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, pad=2, color="#444444")
    ax.text(
        0.03, 0.03, metric_text,
        transform=ax.transAxes,
        ha="left", va="bottom",
        fontsize=6.2, color="#666666",
    )


def _save(fig: plt.Figure, out: Path, name: str) -> None:
    fig.savefig(out / f"{name}.png")
    fig.savefig(out / f"{name}.pdf")
    plt.close(fig)
    print(f"  [saved] {out / (name + '.png')}")
    print(f"  [saved] {out / (name + '.pdf')}")


def _write_summary(out: Path, data: Dict[str, Dict[str, Any]]) -> None:
    summary = {
        ds_key: {
            "label": ds["label"],
            "n": ds["n"],
            "p": ds["p"],
            "layouts": [
                {k: v for k, v in layout.items() if k != "coords"}
                for layout in ds["layouts"]
            ],
        }
        for ds_key, ds in data.items()
    }
    with open(out / "fig5_umap_layouts_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"  [saved] {out / 'fig5_umap_layouts_summary.json'}")


def fig_global_triptych(out: Path, data: Dict[str, Dict[str, Any]], ds_key: str = "global_snp") -> None:
    ds = data[ds_key]
    x0, x1, y0, y1 = _panel_limits(ds["layouts"])
    fig, axes = plt.subplots(1, 3, figsize=(170 * MM, 58 * MM))

    for ax, layout in zip(axes, ds["layouts"]):
        metric_text = (
            f"J_edge={layout['J_edge_pca_vs_umap']:.3f}\n"
            f"J_nbr={layout['J_nbr_pca_vs_umap']:.3f}"
        )
        title = LAYOUT_LABELS[layout["name"]]
        _plot_layout(
            ax,
            layout["coords"],
            ds["partition"],
            title,
            metric_text,
            (x0, x1),
            (y0, y1),
            point_size=1.0,
        )

    fig.tight_layout(w_pad=1.4)
    _save(fig, out, "fig5_umap_layouts_global")


def fig_all_panels(out: Path, data: Dict[str, Dict[str, Any]]) -> None:
    fig, axes = plt.subplots(4, 3, figsize=(170 * MM, 185 * MM))

    for row_idx, ds_key in enumerate(DATASET_ORDER):
        ds = data[ds_key]
        x0, x1, y0, y1 = _panel_limits(ds["layouts"])
        point_size = 0.8 if ds["n"] > 1000 else 4.0

        for col_idx, layout in enumerate(ds["layouts"]):
            ax = axes[row_idx, col_idx]
            metric_text = (
                f"J_edge={layout['J_edge_pca_vs_umap']:.3f} | "
                f"J_nbr={layout['J_nbr_pca_vs_umap']:.3f}"
            )
            title = (
                LAYOUT_LABELS[layout["name"]]
                if row_idx == 0
                else ""
            )
            _plot_layout(
                ax,
                layout["coords"],
                ds["partition"],
                title,
                metric_text,
                (x0, x1),
                (y0, y1),
                point_size=point_size,
            )
            if col_idx == 0:
                ax.set_ylabel(f"{ds['label']}\n(n={ds['n']:,})", fontsize=8)

    fig.suptitle(
        "ADR-0014 real UMAP comparisons: baseline / local-stress / global-stress",
        fontsize=10,
        y=0.995,
    )
    fig.text(
        0.5, 0.003,
        "Colors come from a fixed 2-way partition in PCA-30D; UMAP only changes layout.",
        ha="center", va="bottom", fontsize=7, color="#555555",
    )
    fig.tight_layout(rect=(0.03, 0.02, 1.0, 0.98), h_pad=1.0, w_pad=0.8)
    _save(fig, out, "fig5_umap_layouts_all")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate real ADR-0014 UMAP comparison figures."
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=ROOT / "docs" / "poster" / "figures",
        help="Output directory for PNG/PDF figures.",
    )
    parser.add_argument(
        "--only-global",
        action="store_true",
        help="Generate only the representative Global SNP triptych.",
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    _style()

    data = {"global_snp": _compute_layouts("global_snp")}
    _write_summary(args.outdir, data)
    fig_global_triptych(args.outdir, data)

    if args.only_global:
        return

    for ds_key in DATASET_ORDER[1:]:
        data[ds_key] = _compute_layouts(ds_key)

    _write_summary(args.outdir, data)
    fig_all_panels(args.outdir, data)


if __name__ == "__main__":
    main()
