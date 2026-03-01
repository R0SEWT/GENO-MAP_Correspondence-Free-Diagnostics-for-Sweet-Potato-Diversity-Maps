#!/usr/bin/env python3
"""Generate paper-ready figures from validation JSON artefacts.

All figures are saved as 300 DPI PNG.  Each figure also writes a companion
``*_data.json`` with the underlying numbers for reproducibility.

Figures produced
----------------
  fig_trustworthiness_comparison  — bar chart per dataset / method
  fig_stability_heatmap           — Jaccard inter-seed heat-map
  fig_imputation_sensitivity      — delta trustworthiness + Jaccard
  fig_embedding_scatter           — 2D UMAP scatter (PCA baseline)
  fig_pca_variance                — cumulative explained variance
  fig_missingness_structured      — sorted heatmap (reuses EDA data)

Usage::

    python scripts/plot_results.py --experiment-dir experiments/global_snp/v2-baseline
    python scripts/plot_results.py --experiment-dir experiments/ --recursive
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

DPI = 300
BBOX = "tight"
PALETTE = {
    "PCA": "#4C72B0",
    "Autoencoder": "#DD8452",
    "Transformer": "#55A868",
    "baseline": "#4C72B0",
    "alt": "#DD8452",
}


def _save(fig: plt.Figure, path: Path, data: Any = None) -> None:
    """Save figure and optional companion data JSON."""
    fig.savefig(path, dpi=DPI, bbox_inches=BBOX)
    plt.close(fig)
    print(f"  [saved] {path}")
    if data is not None:
        data_path = path.with_suffix(".json")
        with open(data_path, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)


# ---------------------------------------------------------------------------
# Figure generators
# ---------------------------------------------------------------------------


def fig_trustworthiness_comparison(
    validations: Dict[str, Dict],
    out_dir: Path,
) -> None:
    """Bar chart: trustworthiness per dataset (grouped by method)."""
    datasets = sorted(validations.keys())
    means = [validations[d]["trustworthiness"]["mean"] for d in datasets]
    stds = [validations[d]["trustworthiness"]["std"] for d in datasets]

    fig, ax = plt.subplots(figsize=(max(6, len(datasets) * 1.5), 5))
    x = np.arange(len(datasets))
    bars = ax.bar(x, means, yerr=stds, capsize=4, color=PALETTE["PCA"],
                  edgecolor="white", linewidth=0.5, label="PCA baseline")

    ax.set_xticks(x)
    ax.set_xticklabels([d.replace("_", "\n") for d in datasets], fontsize=9)
    ax.set_ylabel("Trustworthiness")
    ax.set_title("Trustworthiness por dataset (PCA → UMAP → kNN)")
    ax.set_ylim(0.8, 1.0)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    _save(fig, out_dir / "fig_trustworthiness_comparison.png",
          data={"datasets": datasets, "means": means, "stds": stds})


def fig_stability_heatmap(
    validations: Dict[str, Dict],
    out_dir: Path,
) -> None:
    """Heatmap of Jaccard stability metrics per dataset."""
    datasets = sorted(validations.keys())
    metrics = ["pca_neighbour_jaccard_mean", "umap_neighbour_jaccard_mean", "edge_jaccard_mean"]
    labels = ["PCA neighbours", "UMAP neighbours", "Edge overlap"]

    matrix = []
    for d in datasets:
        stab = validations[d].get("stability")
        row = []
        for m in metrics:
            row.append(stab.get(m, float("nan")) if stab else float("nan"))
        matrix.append(row)
    matrix = np.array(matrix)

    fig, ax = plt.subplots(figsize=(max(6, len(labels) * 2), max(4, len(datasets) * 0.8)))
    im = ax.imshow(matrix, cmap="YlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_yticks(range(len(datasets)))
    ax.set_yticklabels([d.replace("_", "\n") for d in datasets], fontsize=9)
    for i in range(len(datasets)):
        for j in range(len(labels)):
            v = matrix[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.3f}", ha="center", va="center", fontsize=8,
                        color="white" if v > 0.7 else "black")
    fig.colorbar(im, ax=ax, label="Jaccard")
    ax.set_title("Multi-seed stability (Jaccard)")
    fig.tight_layout()
    _save(fig, out_dir / "fig_stability_heatmap.png",
          data={"datasets": datasets, "metrics": labels,
                "matrix": matrix.tolist()})


def fig_imputation_sensitivity(
    validations: Dict[str, Dict],
    out_dir: Path,
) -> None:
    """Grouped bar: imputation sensitivity (neighbour Jaccard + trust delta)."""
    datasets = sorted(validations.keys())
    nbr_jaccards, trust_deltas = [], []
    for d in datasets:
        sens = validations[d].get("imputation_sensitivity", {})
        comp = sens.get("most_frequent_vs_median", {})
        nbr_jaccards.append(comp.get("neighbour_jaccard", float("nan")))
        trust_deltas.append(comp.get("trustworthiness_delta", float("nan")))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    x = np.arange(len(datasets))

    axes[0].bar(x, nbr_jaccards, color=PALETTE["baseline"], edgecolor="white")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([d.replace("_", "\n") for d in datasets], fontsize=8)
    axes[0].set_ylabel("Jaccard (neighbours)")
    axes[0].set_title("Neighbourhood overlap: mode vs median")
    axes[0].set_ylim(0, 1)
    axes[0].grid(axis="y", alpha=0.3)

    colors = ["#55A868" if d >= 0 else "#C44E52" for d in trust_deltas]
    axes[1].bar(x, trust_deltas, color=colors, edgecolor="white")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([d.replace("_", "\n") for d in datasets], fontsize=8)
    axes[1].set_ylabel("Δ Trustworthiness")
    axes[1].set_title("Trustworthiness: median − mode")
    axes[1].axhline(0, color="grey", linewidth=0.5)
    axes[1].grid(axis="y", alpha=0.3)

    fig.suptitle("Imputation sensitivity (ADR-0003)", fontsize=12)
    fig.tight_layout()
    _save(fig, out_dir / "fig_imputation_sensitivity.png",
          data={"datasets": datasets, "neighbour_jaccard": nbr_jaccards,
                "trust_delta": trust_deltas})


def fig_embedding_scatter(
    nodes_file: Path,
    out_dir: Path,
    dataset_name: str = "",
) -> None:
    """2D scatter from nodes JSON (UMAP embedding)."""
    with open(nodes_file, "r") as f:
        nodes = json.load(f)

    xs = [n["embedding"][0] for n in nodes if n["embedding"]]
    ys = [n["embedding"][1] for n in nodes if n["embedding"]]

    # Try to colour by Institution if available
    institutions = [n.get("meta", {}).get("Institution", "unknown") for n in nodes if n["embedding"]]
    unique_inst = sorted(set(institutions))
    cmap = plt.cm.get_cmap("tab20", max(len(unique_inst), 1))
    inst_to_idx = {inst: i for i, inst in enumerate(unique_inst)}
    colors = [cmap(inst_to_idx[inst]) for inst in institutions]

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(xs, ys, c=colors, s=6, alpha=0.6, edgecolors="none")
    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")
    title = f"Embedding scatter — {dataset_name}" if dataset_name else "Embedding scatter"
    ax.set_title(title)
    ax.set_aspect("equal")
    ax.grid(alpha=0.15)

    # Legend if few institutions
    if 1 < len(unique_inst) <= 15:
        from matplotlib.lines import Line2D
        handles = [Line2D([0], [0], marker="o", linestyle="", color=cmap(inst_to_idx[inst]),
                          markersize=5, label=inst[:30])
                   for inst in unique_inst]
        ax.legend(handles=handles, fontsize=7, loc="upper right", ncol=1)

    fig.tight_layout()
    _save(fig, out_dir / f"fig_embedding_scatter_{dataset_name or 'combined'}.png",
          data={"n_nodes": len(xs), "institutions": unique_inst})


def fig_pca_variance(
    validations: Dict[str, Dict],
    out_dir: Path,
) -> None:
    """Cumulative PCA variance plot per dataset."""
    datasets = sorted(validations.keys())
    fig, ax = plt.subplots(figsize=(8, 5))

    for d in datasets:
        var = validations[d].get("pca_loadings", {}).get("variance_explained", [])
        if var:
            cumvar = np.cumsum(var)
            ax.plot(range(1, len(cumvar) + 1), cumvar, marker="o", markersize=3,
                    label=d.replace("_", " "), linewidth=1.5)

    ax.set_xlabel("Number of PCs")
    ax.set_ylabel("Cumulative explained variance")
    ax.set_title("PCA variance explained")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_ylim(0, 1.05)
    fig.tight_layout()

    data_out: Dict[str, Any] = {}
    for d in datasets:
        var = validations[d].get("pca_loadings", {}).get("variance_explained", [])
        data_out[d] = var
    _save(fig, out_dir / "fig_pca_variance.png", data=data_out)


# ---------------------------------------------------------------------------
# Collector: find all validation JSONs
# ---------------------------------------------------------------------------


def collect_validations(experiment_dir: Path, recursive: bool = False) -> Dict[str, Dict]:
    """Find ``*_validation.json`` files and index by dataset name."""
    pattern = "**/*_validation.json" if recursive else "*_validation.json"
    validations: Dict[str, Dict] = {}
    for f in sorted(experiment_dir.glob(pattern)):
        with open(f, "r") as fp:
            data = json.load(fp)
        # Infer dataset name from directory or filename
        parts = f.stem.replace("_validation", "").split("_", 1)
        ds_name = parts[-1] if len(parts) > 1 else f.parent.name
        # Prefer longer / more specific name
        if ds_name in validations:
            ds_name = f"{f.parent.name}_{ds_name}"
        validations[ds_name] = data
    return validations


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate paper-ready figures from validation JSONs.")
    parser.add_argument("--experiment-dir", required=True, type=Path,
                        help="Directory with *_validation.json files")
    parser.add_argument("--recursive", action="store_true",
                        help="Search subdirectories recursively")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="Output directory for figures (default: experiment-dir/figures)")
    parser.add_argument("--nodes", type=Path, default=None,
                        help="Specific nodes JSON for scatter plot")
    parser.add_argument("--dataset-name", default="",
                        help="Dataset label for scatter plot")
    args = parser.parse_args()

    out_dir = args.out_dir or args.experiment_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    validations = collect_validations(args.experiment_dir, args.recursive)
    if not validations:
        print(f"[warn] No *_validation.json found in {args.experiment_dir}")
        return

    print(f"Found {len(validations)} validation(s): {list(validations.keys())}")
    print(f"Output → {out_dir}\n")

    # Generate figures
    fig_trustworthiness_comparison(validations, out_dir)
    fig_stability_heatmap(validations, out_dir)
    fig_imputation_sensitivity(validations, out_dir)
    fig_pca_variance(validations, out_dir)

    # Scatter if nodes file provided
    if args.nodes and args.nodes.exists():
        fig_embedding_scatter(args.nodes, out_dir, args.dataset_name)
    else:
        # Look for node files in experiment dir
        node_files = list(args.experiment_dir.glob("**/*_nodes.json"))
        for nf in node_files[:3]:  # max 3 scatter plots
            ds = nf.stem.replace("_nodes", "").split("_", 1)[-1]
            fig_embedding_scatter(nf, out_dir, ds)

    # Also copy figures to docs/figures for paper
    docs_fig_dir = Path(__file__).resolve().parent.parent / "docs" / "figures"
    docs_fig_dir.mkdir(parents=True, exist_ok=True)
    for png in out_dir.glob("fig_*.png"):
        import shutil
        dest = docs_fig_dir / png.name
        shutil.copy2(png, dest)
        print(f"  [copy] {png.name} → docs/figures/")

    print(f"\nDone. {len(list(out_dir.glob('fig_*.png')))} figures generated.")


if __name__ == "__main__":
    main()
