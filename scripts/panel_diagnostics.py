#!/usr/bin/env python3
"""Módulo 1 — Panel QA & Geometry Diagnostics.

Correspondence-free characterization of each genotyping panel.
No cross-panel alignment; all metrics are intra-dataset.

Metrics:
  - Effective rank (participation ratio of PCA eigenvalues)
  - PC1 / PC2 dominance (variance explained)
  - kNN distance statistics (mean, std, CV, skewness)
  - Reciprocal kNN rate (fraction of mutual edges)
  - Giant component ratio (connectivity)
  - Automatic QA flags

Usage::

    python scripts/panel_diagnostics.py                    # all datasets
    python scripts/panel_diagnostics.py --dataset global_snp
    python scripts/panel_diagnostics.py --k 20 --out-dir experiments/diagnostics
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd
from scipy import stats as sp_stats
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).resolve().parent))
from load_dart import filter_missingness, load_genotypes  # noqa: E402

# ╭──────────────────────────────────────────────────────────────╮
# │  Dataset registry                                           │
# ╰──────────────────────────────────────────────────────────────╯

DATASETS = {
    "global_snp": {
        "path": "data/10.21223P30BVZYY_Genetic_diversity/SNP_Genotypes.csv",
        "metric": "cosine",
        "label": "Global SNP",
    },
    "global_silico": {
        "path": "data/10.21223P30BVZYY_Genetic_diversity/SilicoDArT_Genotypes.csv",
        "metric": "jaccard",
        "label": "Global SilicoDArT",
    },
    "lowdensity_snp": {
        "path": "data/10.21223P3UBDJ44_LowDensity/01_Report_DSp25-515_SNPs_Filtered_by _reads.csv",
        "metric": "cosine",
        "label": "LowDensity SNP",
    },
    "lowdensity_silico": {
        "path": "data/10.21223P3UBDJ44_LowDensity/02_Report_DSp25-515_Silico-DArT_Filtered_by_reads.csv",
        "metric": "jaccard",
        "label": "LowDensity SilicoDArT",
    },
}


# ╭──────────────────────────────────────────────────────────────╮
# │  Metric functions                                           │
# ╰──────────────────────────────────────────────────────────────╯

def effective_rank(eigenvalues: np.ndarray) -> float:
    """Participation ratio: (Σλ)² / Σλ², measures effective dimensionality."""
    ev = eigenvalues[eigenvalues > 0]
    return float((ev.sum() ** 2) / (ev ** 2).sum())


def pca_spectrum(X: np.ndarray, n_components: int = 50) -> Dict[str, Any]:
    """Compute PCA spectrum and related diagnostics."""
    nc = min(n_components, X.shape[1], X.shape[0])
    pca = PCA(n_components=nc, random_state=42)
    pca.fit(X)

    ev = pca.explained_variance_
    evr = pca.explained_variance_ratio_

    return {
        "n_components_fit": int(nc),
        "eigenvalues": ev.tolist(),
        "variance_ratios": evr.tolist(),
        "pc1_var_pct": round(float(evr[0]) * 100, 4),
        "pc2_var_pct": round(float(evr[1]) * 100, 4) if len(evr) > 1 else 0.0,
        "pc1_pc2_ratio": round(float(evr[0] / evr[1]), 2) if len(evr) > 1 and evr[1] > 0 else float("inf"),
        "cumvar_5": round(float(evr[:5].sum()) * 100, 2),
        "cumvar_10": round(float(evr[:min(10, nc)].sum()) * 100, 2),
        "cumvar_30": round(float(evr[:min(30, nc)].sum()) * 100, 2),
        "effective_rank": round(effective_rank(ev), 2),
    }


def knn_distance_stats(
    X: np.ndarray, k: int = 15, metric: str = "cosine"
) -> Dict[str, Any]:
    """kNN distance distribution statistics."""
    k_eff = min(k, X.shape[0] - 1)
    nbrs = NearestNeighbors(n_neighbors=k_eff + 1, metric=metric)
    nbrs.fit(X)
    dists, inds = nbrs.kneighbors(X)

    # Exclude self (first column)
    knn_dists = dists[:, 1:].ravel()

    return {
        "k": k_eff,
        "metric": metric,
        "mean": round(float(np.mean(knn_dists)), 6),
        "std": round(float(np.std(knn_dists)), 6),
        "cv": round(float(np.std(knn_dists) / np.mean(knn_dists)), 4) if np.mean(knn_dists) > 0 else 0,
        "median": round(float(np.median(knn_dists)), 6),
        "p5": round(float(np.percentile(knn_dists, 5)), 6),
        "p95": round(float(np.percentile(knn_dists, 95)), 6),
        "skewness": round(float(sp_stats.skew(knn_dists)), 4),
        "kurtosis": round(float(sp_stats.kurtosis(knn_dists)), 4),
    }


def reciprocal_knn_rate(X: np.ndarray, k: int = 15, metric: str = "cosine") -> Dict[str, Any]:
    """Fraction of kNN edges that are mutual (i→j AND j→i)."""
    k_eff = min(k, X.shape[0] - 1)
    nbrs = NearestNeighbors(n_neighbors=k_eff + 1, metric=metric)
    nbrs.fit(X)
    _, inds = nbrs.kneighbors(X)

    # Build directed edge set (excluding self)
    directed = set()
    for i in range(inds.shape[0]):
        for j in inds[i, 1:]:
            directed.add((i, j))

    mutual = sum(1 for (i, j) in directed if (j, i) in directed)
    total = len(directed)
    rate = mutual / total if total > 0 else 1.0

    return {
        "k": k_eff,
        "total_directed_edges": total,
        "mutual_edges": mutual,
        "reciprocal_rate": round(rate, 4),
    }


def component_analysis(X: np.ndarray, k: int = 15, metric: str = "cosine") -> Dict[str, Any]:
    """Analyze connected components of the kNN graph."""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components

    n = X.shape[0]
    k_eff = min(k, n - 1)
    nbrs = NearestNeighbors(n_neighbors=k_eff + 1, metric=metric)
    nbrs.fit(X)
    dists, inds = nbrs.kneighbors(X)

    # Build symmetric adjacency
    rows, cols = [], []
    for i in range(n):
        for j in inds[i, 1:]:
            rows.extend([i, j])
            cols.extend([j, i])

    adj = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    n_comp, labels = connected_components(adj, directed=False)

    # Giant component
    comp_sizes = np.bincount(labels)
    giant = int(comp_sizes.max())

    return {
        "n_components": int(n_comp),
        "giant_component_size": giant,
        "giant_component_ratio": round(giant / n, 4),
        "component_sizes": sorted(comp_sizes.tolist(), reverse=True)[:10],
    }


def degree_stats(X: np.ndarray, k: int = 15, metric: str = "cosine") -> Dict[str, Any]:
    """Degree distribution of undirected kNN graph."""
    n = X.shape[0]
    k_eff = min(k, n - 1)
    nbrs = NearestNeighbors(n_neighbors=k_eff + 1, metric=metric)
    nbrs.fit(X)
    _, inds = nbrs.kneighbors(X)

    # Count undirected degree (union of in + out neighbors)
    neighbors_sets = [set() for _ in range(n)]
    for i in range(n):
        for j in inds[i, 1:]:
            neighbors_sets[i].add(j)
            neighbors_sets[j].add(i)

    degrees = np.array([len(s) for s in neighbors_sets])

    return {
        "mean_degree": round(float(np.mean(degrees)), 2),
        "std_degree": round(float(np.std(degrees)), 2),
        "min_degree": int(degrees.min()),
        "max_degree": int(degrees.max()),
        "cv_degree": round(float(np.std(degrees) / np.mean(degrees)), 4) if np.mean(degrees) > 0 else 0,
    }


def missingness_stats(X_raw: pd.DataFrame) -> Dict[str, Any]:
    """Compute missingness statistics before imputation."""
    total_cells = X_raw.size
    missing_mask = X_raw.isna() | (X_raw == "-") | (X_raw == "")
    n_missing = int(missing_mask.sum().sum())
    rate = n_missing / total_cells if total_cells > 0 else 0

    sample_miss = missing_mask.mean(axis=1)
    marker_miss = missing_mask.mean(axis=0)

    return {
        "total_cells": total_cells,
        "missing_cells": n_missing,
        "missing_rate": round(rate, 4),
        "sample_miss_mean": round(float(sample_miss.mean()), 4),
        "sample_miss_max": round(float(sample_miss.max()), 4),
        "marker_miss_mean": round(float(marker_miss.mean()), 4),
        "marker_miss_max": round(float(marker_miss.max()), 4),
    }


# ╭──────────────────────────────────────────────────────────────╮
# │  QA Flags                                                   │
# ╰──────────────────────────────────────────────────────────────╯

def compute_flags(diag: Dict[str, Any]) -> List[str]:
    """Generate automatic QA flags from diagnostics."""
    flags = []

    # 1D-dominant: PC1 explains >50% of variance
    if diag["pca"]["pc1_var_pct"] > 50:
        flags.append(f"1D-DOMINANT (PC1={diag['pca']['pc1_var_pct']:.1f}%)")

    # Low effective rank
    if diag["pca"]["effective_rank"] < 5:
        flags.append(f"LOW-RANK (eff_rank={diag['pca']['effective_rank']:.1f})")

    # Collapsed neighborhood: very low CV of kNN distances
    if diag["knn_distances"]["cv"] < 0.05:
        flags.append(f"COLLAPSED-NEIGHBORHOOD (cv={diag['knn_distances']['cv']:.3f})")

    # High missingness
    if diag["missingness"]["missing_rate"] > 0.10:
        flags.append(f"HIGH-MISSINGNESS ({diag['missingness']['missing_rate']*100:.1f}%)")

    # Disconnected graph
    if diag["components"]["n_components"] > 1:
        flags.append(f"DISCONNECTED ({diag['components']['n_components']} components)")

    # Low reciprocity
    if diag["reciprocal_knn"]["reciprocal_rate"] < 0.4:
        flags.append(f"LOW-RECIPROCITY ({diag['reciprocal_knn']['reciprocal_rate']:.2f})")

    # Very few samples relative to markers
    ratio = diag["n_samples"] / diag["n_markers"]
    if ratio < 0.02:
        flags.append(f"EXTREME-WIDE (n/p={ratio:.4f})")

    return flags if flags else ["OK"]


# ╭──────────────────────────────────────────────────────────────╮
# │  Main                                                       │
# ╰──────────────────────────────────────────────────────────────╯

def diagnose_panel(
    ds_key: str, k: int = 15, out_dir: Path | None = None
) -> Dict[str, Any]:
    """Run full diagnostics for one dataset."""
    ds = DATASETS[ds_key]
    root = Path(__file__).resolve().parent.parent
    print(f"\n{'='*60}")
    print(f"  PANEL DIAGNOSTICS: {ds['label']} ({ds_key})")
    print(f"{'='*60}")

    t0 = time.time()

    # Load raw data
    X_raw, sample_ids, _ = load_genotypes(root / ds["path"])
    miss = missingness_stats(X_raw)

    # Filter + impute
    X_raw, _ = filter_missingness(X_raw, sample_thresh=0.50, marker_thresh=0.50)
    imp = SimpleImputer(strategy="most_frequent")
    X = imp.fit_transform(X_raw)

    n_samples, n_markers = X.shape
    print(f"  Shape: {n_samples} × {n_markers}")
    print(f"  Missing rate: {miss['missing_rate']*100:.2f}%")

    # PCA spectrum
    print("  Computing PCA spectrum...")
    pca_diag = pca_spectrum(X)
    print(f"    PC1={pca_diag['pc1_var_pct']:.2f}%  PC2={pca_diag['pc2_var_pct']:.2f}%  "
          f"eff_rank={pca_diag['effective_rank']:.1f}")

    # Use PCA 30D for kNN (same as pipeline)
    # Always use cosine on PCA space (Jaccard is meaningless for continuous PCs)
    knn_metric = "cosine"
    nc_knn = min(30, pca_diag["n_components_fit"])
    pca_knn = PCA(n_components=nc_knn, random_state=42)
    X_pca30 = pca_knn.fit_transform(X)

    # kNN diagnostics
    print(f"  Computing kNN diagnostics (k={k}, metric={knn_metric} on PCA space)...")
    knn_dist = knn_distance_stats(X_pca30, k=k, metric=knn_metric)
    print(f"    mean_dist={knn_dist['mean']:.6f}  cv={knn_dist['cv']:.4f}  skew={knn_dist['skewness']:.4f}")

    recip = reciprocal_knn_rate(X_pca30, k=k, metric=knn_metric)
    print(f"    reciprocal_rate={recip['reciprocal_rate']:.4f}")

    comp = component_analysis(X_pca30, k=k, metric=knn_metric)
    print(f"    components={comp['n_components']}  giant_ratio={comp['giant_component_ratio']:.4f}")

    deg = degree_stats(X_pca30, k=k, metric=knn_metric)
    print(f"    mean_degree={deg['mean_degree']:.1f}  std={deg['std_degree']:.1f}")

    elapsed = time.time() - t0

    result = {
        "dataset": ds_key,
        "label": ds["label"],
        "n_samples": n_samples,
        "n_markers": n_markers,
        "ratio_n_p": round(n_samples / n_markers, 4),
        "missingness": miss,
        "pca": pca_diag,
        "knn_distances": knn_dist,
        "reciprocal_knn": recip,
        "components": comp,
        "degree_stats": deg,
        "elapsed_s": round(elapsed, 1),
    }

    # QA flags
    flags = compute_flags(result)
    result["flags"] = flags
    print(f"  Flags: {flags}")
    print(f"  Done in {elapsed:.1f}s")

    # Save
    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)
        with open(out_dir / f"diagnostics_{ds_key}.json", "w") as f:
            # Remove large eigenvalue arrays from saved version
            save_result = {k: v for k, v in result.items()}
            save_result["pca"] = {k: v for k, v in result["pca"].items()
                                  if k not in ("eigenvalues", "variance_ratios")}
            json.dump(save_result, f, indent=2, ensure_ascii=False)

    return result


def main():
    parser = argparse.ArgumentParser(description="Panel QA & Geometry Diagnostics")
    parser.add_argument("--dataset", default=None,
                        choices=list(DATASETS.keys()),
                        help="Single dataset (default: all)")
    parser.add_argument("--k", type=int, default=15, help="kNN neighbors")
    parser.add_argument("--out-dir", type=Path,
                        default=Path("experiments/diagnostics"))
    args = parser.parse_args()

    datasets = [args.dataset] if args.dataset else list(DATASETS.keys())
    all_results = []

    for ds_key in datasets:
        result = diagnose_panel(ds_key, k=args.k, out_dir=args.out_dir)
        all_results.append(result)

    # Print summary table
    print(f"\n{'='*80}")
    print("  PANEL DIAGNOSTICS SUMMARY")
    print(f"{'='*80}")
    print(f"{'Dataset':25s} {'n':>6s} {'p':>7s} {'n/p':>6s} {'miss%':>6s} "
          f"{'PC1%':>6s} {'eff_rk':>7s} {'recip':>6s} {'comp':>5s} {'flags'}")
    for r in all_results:
        flags_short = "; ".join(r["flags"])
        print(f"{r['label']:25s} {r['n_samples']:6d} {r['n_markers']:7d} "
              f"{r['ratio_n_p']:6.3f} {r['missingness']['missing_rate']*100:5.1f}% "
              f"{r['pca']['pc1_var_pct']:5.1f}% {r['pca']['effective_rank']:7.1f} "
              f"{r['reciprocal_knn']['reciprocal_rate']:6.3f} "
              f"{r['components']['n_components']:5d} {flags_short}")

    # Save combined
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    combined_path = out_dir / "diagnostics_all.json"
    combined = []
    for r in all_results:
        save_r = {k: v for k, v in r.items()}
        save_r["pca"] = {k: v for k, v in r["pca"].items()
                         if k not in ("eigenvalues", "variance_ratios")}
        combined.append(save_r)
    with open(combined_path, "w") as f:
        json.dump(combined, f, indent=2, ensure_ascii=False)
    print(f"\n[saved] {combined_path}")


if __name__ == "__main__":
    main()
