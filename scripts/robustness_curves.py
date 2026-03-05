#!/usr/bin/env python3
"""Módulo 2 — Robustness Curves.

Correspondence-free robustness analysis: measures sensitivity of
graph topology and PCA geometry to controlled perturbations.

Perturbation axes:
  1. Marker subsampling  — random columns at {5,10,20,50,80}%
  2. Missing injection   — +{0,5,10,20}% MCAR on top of existing missingness
  3. Imputation strategy — mode vs median

Metrics (all intra-dataset):
  - Jaccard of kNN neighbors  (mean per-node set overlap)
  - Jaccard of kNN edge sets  (global graph overlap)
  - kNN distance drift         (shift in distance distribution)
  - PCA subspace stability     (principal angles between subspaces)

Usage::

    python scripts/robustness_curves.py                           # all datasets
    python scripts/robustness_curves.py --dataset global_snp --seeds 42
    python scripts/robustness_curves.py --profile quick           # fewer combos
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd
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

PROFILES = {
    "quick": {
        "marker_fractions": [0.20, 0.50, 0.80],
        "missing_extra": [0.0, 0.10],
        "seeds": [42],
    },
    "full": {
        "marker_fractions": [0.05, 0.10, 0.20, 0.50, 0.80],
        "missing_extra": [0.0, 0.05, 0.10, 0.20],
        "seeds": [42, 52, 62],
    },
}


# ╭──────────────────────────────────────────────────────────────╮
# │  Core pipeline: impute → PCA → kNN                         │
# ╰──────────────────────────────────────────────────────────────╯

def _impute(X_raw: np.ndarray, strategy: str = "most_frequent") -> np.ndarray:
    """Impute NaN / -9 using sklearn."""
    imp = SimpleImputer(strategy=strategy, missing_values=np.nan)
    return imp.fit_transform(X_raw)


def _pipeline(
    X: np.ndarray,
    k: int,
    metric: str = "cosine",
    n_pca: int = 50,
    seed: int = 42,
) -> Dict[str, Any]:
    """PCA → kNN on first 30 PCs.  Always uses cosine on PCA space."""
    nc = min(n_pca, X.shape[0], X.shape[1])
    pca = PCA(n_components=nc, random_state=seed)
    X_pca = pca.fit_transform(X)

    knn_dim = min(30, nc)
    feats = X_pca[:, :knn_dim]

    # Always cosine on PCA space (Jaccard meaningless for continuous features)
    k_eff = min(k, X.shape[0] - 1)
    nbrs = NearestNeighbors(n_neighbors=k_eff + 1, metric="cosine")
    nbrs.fit(feats)
    dists, inds = nbrs.kneighbors(feats)

    return {
        "pca": pca,
        "X_pca": X_pca,
        "features": feats,
        "knn_dists": dists[:, 1:],   # exclude self
        "knn_inds": inds[:, 1:],
        "k_eff": k_eff,
    }


# ╭──────────────────────────────────────────────────────────────╮
# │  Metric functions                                           │
# ╰──────────────────────────────────────────────────────────────╯

def jaccard_neighbors(inds_a: np.ndarray, inds_b: np.ndarray) -> float:
    """Mean per-node Jaccard of kNN neighbor sets (by index)."""
    scores = []
    for a_row, b_row in zip(inds_a, inds_b):
        sa, sb = set(a_row), set(b_row)
        if not sa and not sb:
            scores.append(1.0)
        elif not sa or not sb:
            scores.append(0.0)
        else:
            scores.append(len(sa & sb) / len(sa | sb))
    return float(np.mean(scores))


def jaccard_edges(inds_a: np.ndarray, inds_b: np.ndarray) -> float:
    """Jaccard of global directed edge sets."""
    def _edge_set(inds: np.ndarray) -> set:
        edges = set()
        for i, row in enumerate(inds):
            for j in row:
                edges.add((i, j))
        return edges
    sa, sb = _edge_set(inds_a), _edge_set(inds_b)
    if not sa and not sb:
        return 1.0
    return len(sa & sb) / len(sa | sb)


def knn_distance_drift(dists_ref: np.ndarray, dists_pert: np.ndarray) -> Dict[str, float]:
    """Distribution shift of kNN distances (KS statistic + mean shift)."""
    from scipy.stats import ks_2samp

    flat_ref = dists_ref.ravel()
    flat_pert = dists_pert.ravel()

    ks_stat, ks_p = ks_2samp(flat_ref, flat_pert)
    mean_ref = float(np.mean(flat_ref))
    mean_pert = float(np.mean(flat_pert))
    rel_shift = (mean_pert - mean_ref) / mean_ref if mean_ref > 0 else 0

    return {
        "ks_statistic": round(float(ks_stat), 6),
        "ks_pvalue": float(ks_p),
        "mean_ref": round(mean_ref, 6),
        "mean_pert": round(mean_pert, 6),
        "relative_shift": round(float(rel_shift), 6),
    }


def pca_subspace_similarity(
    X_pca_ref: np.ndarray, X_pca_pert: np.ndarray, n_comp: int = 10
) -> Dict[str, Any]:
    """Subspace similarity via principal angles between PCA score matrices.

    Uses the sample-space PCA scores (n × nc) so that comparisons work
    even when the marker counts differ (e.g., marker subsampling).

    Principal angles are computed from SVD of Q1^T Q2 where Q1, Q2 are
    orthonormal bases obtained from QR of the score matrices.
    """
    nc = min(n_comp, X_pca_ref.shape[1], X_pca_pert.shape[1])
    # QR for orthonormal basis of each subspace in sample space
    Q1, _ = np.linalg.qr(X_pca_ref[:, :nc])    # (n, nc)
    Q2, _ = np.linalg.qr(X_pca_pert[:, :nc])

    M = Q1.T @ Q2
    s = np.linalg.svd(M, compute_uv=False)
    s = np.clip(s, -1, 1)
    angles_deg = np.degrees(np.arccos(s))

    return {
        "n_components": int(nc),
        "cosines": [round(float(c), 6) for c in s],
        "angles_deg": [round(float(a), 4) for a in angles_deg],
        "mean_angle_deg": round(float(np.mean(angles_deg)), 4),
        "max_angle_deg": round(float(np.max(angles_deg)), 4),
        "subspace_similarity": round(float(np.mean(s)), 6),
    }


# ╭──────────────────────────────────────────────────────────────╮
# │  Perturbation functions                                     │
# ╰──────────────────────────────────────────────────────────────╯

def subsample_markers(X: np.ndarray, frac: float, rng: np.random.Generator) -> np.ndarray:
    """Keep a random fraction of columns (markers)."""
    n_keep = max(1, int(X.shape[1] * frac))
    cols = rng.choice(X.shape[1], size=n_keep, replace=False)
    return X[:, np.sort(cols)]


def inject_missing(
    X: np.ndarray, extra_rate: float, rng: np.random.Generator
) -> np.ndarray:
    """Inject additional MCAR missingness on a numeric matrix (NaN)."""
    if extra_rate <= 0:
        return X.copy()
    X_out = X.copy().astype(float)
    mask = rng.random(X_out.shape) < extra_rate
    X_out[mask] = np.nan
    return X_out


# ╭──────────────────────────────────────────────────────────────╮
# │  Experiment runner                                          │
# ╰──────────────────────────────────────────────────────────────╯

def run_marker_subsampling(
    X_imputed: np.ndarray,
    ref_result: Dict,
    fractions: List[float],
    seeds: List[int],
    k: int,
    metric: str,
) -> List[Dict[str, Any]]:
    """Run marker subsampling curve: compare perturbed topology to reference."""
    rows = []
    for frac in fractions:
        for seed in seeds:
            rng = np.random.default_rng(seed)
            X_sub = subsample_markers(X_imputed, frac, rng)
            pert = _pipeline(X_sub, k, metric, seed=seed)

            jn = jaccard_neighbors(ref_result["knn_inds"], pert["knn_inds"])
            je = jaccard_edges(ref_result["knn_inds"], pert["knn_inds"])
            drift = knn_distance_drift(ref_result["knn_dists"], pert["knn_dists"])
            pa = pca_subspace_similarity(ref_result["X_pca"], pert["X_pca"])

            rows.append({
                "perturbation": "marker_subsampling",
                "fraction": frac,
                "seed": seed,
                "jaccard_neighbors": round(jn, 4),
                "jaccard_edges": round(je, 4),
                "distance_ks": drift["ks_statistic"],
                "distance_rel_shift": drift["relative_shift"],
                "subspace_similarity": pa["subspace_similarity"],
                "mean_angle_deg": pa["mean_angle_deg"],
            })
    return rows


def run_missing_injection(
    X_raw_numeric: np.ndarray,
    ref_result: Dict,
    extra_rates: List[float],
    seeds: List[int],
    k: int,
    metric: str,
    imputation: str = "most_frequent",
) -> List[Dict[str, Any]]:
    """Run missing injection curve: add MCAR then re-impute and compare."""
    rows = []
    for rate in extra_rates:
        for seed in seeds:
            rng = np.random.default_rng(seed)
            X_miss = inject_missing(X_raw_numeric, rate, rng)
            X_imp = _impute(X_miss, strategy=imputation)
            pert = _pipeline(X_imp, k, metric, seed=seed)

            jn = jaccard_neighbors(ref_result["knn_inds"], pert["knn_inds"])
            je = jaccard_edges(ref_result["knn_inds"], pert["knn_inds"])
            drift = knn_distance_drift(ref_result["knn_dists"], pert["knn_dists"])
            pa = pca_subspace_similarity(ref_result["X_pca"], pert["X_pca"])

            rows.append({
                "perturbation": "missing_injection",
                "extra_rate": rate,
                "imputation": imputation,
                "seed": seed,
                "jaccard_neighbors": round(jn, 4),
                "jaccard_edges": round(je, 4),
                "distance_ks": drift["ks_statistic"],
                "distance_rel_shift": drift["relative_shift"],
                "subspace_similarity": pa["subspace_similarity"],
                "mean_angle_deg": pa["mean_angle_deg"],
            })
    return rows


def run_imputation_comparison(
    X_raw_numeric: np.ndarray,
    k: int,
    metric: str,
    seeds: List[int],
) -> List[Dict[str, Any]]:
    """Compare mode vs median imputation on same data."""
    rows = []
    for seed in seeds:
        X_mode = _impute(X_raw_numeric, "most_frequent")
        X_median = _impute(X_raw_numeric, "median")

        ref = _pipeline(X_mode, k, metric, seed=seed)
        pert = _pipeline(X_median, k, metric, seed=seed)

        jn = jaccard_neighbors(ref["knn_inds"], pert["knn_inds"])
        je = jaccard_edges(ref["knn_inds"], pert["knn_inds"])
        drift = knn_distance_drift(ref["knn_dists"], pert["knn_dists"])
        pa = pca_subspace_similarity(ref["X_pca"], pert["X_pca"])

        rows.append({
            "perturbation": "imputation_comparison",
            "strategy_ref": "mode",
            "strategy_alt": "median",
            "seed": seed,
            "jaccard_neighbors": round(jn, 4),
            "jaccard_edges": round(je, 4),
            "distance_ks": drift["ks_statistic"],
            "distance_rel_shift": drift["relative_shift"],
            "subspace_similarity": pa["subspace_similarity"],
            "mean_angle_deg": pa["mean_angle_deg"],
        })
    return rows


# ╭──────────────────────────────────────────────────────────────╮
# │  Dataset-level entry point                                  │
# ╰──────────────────────────────────────────────────────────────╯

def run_dataset(
    ds_key: str,
    marker_fractions: List[float],
    missing_extra: List[float],
    seeds: List[int],
    k: int = 15,
    out_dir: Path | None = None,
) -> Dict[str, Any]:
    """Full robustness analysis for one dataset."""
    ds = DATASETS[ds_key]
    root = Path(__file__).resolve().parent.parent
    print(f"\n{'='*60}")
    print(f"  ROBUSTNESS CURVES: {ds['label']} ({ds_key})")
    print(f"{'='*60}")

    t0 = time.time()

    # Load & filter
    X_raw, sample_ids, _ = load_genotypes(root / ds["path"])
    X_raw, _ = filter_missingness(X_raw, sample_thresh=0.50, marker_thresh=0.50)
    n, p = X_raw.shape
    print(f"  Shape: {n} × {p}")

    # Create numeric version (for missing injection)
    X_numeric = X_raw.apply(pd.to_numeric, errors="coerce").values.astype(float)

    # Reference pipeline (mode imputation, seed=42)
    X_ref_imp = _impute(X_numeric, "most_frequent")
    ref = _pipeline(X_ref_imp, k, ds["metric"], seed=42)
    print(f"  Reference pipeline built (PCA {ref['pca'].n_components_}D, k={ref['k_eff']})")

    # 1. Marker subsampling curves
    print(f"  [1/3] Marker subsampling ({len(marker_fractions)} × {len(seeds)} seeds)...")
    marker_rows = run_marker_subsampling(
        X_ref_imp, ref, marker_fractions, seeds, k, ds["metric"]
    )

    # 2. Missing injection curves
    print(f"  [2/3] Missing injection ({len(missing_extra)} × {len(seeds)} seeds)...")
    missing_rows = run_missing_injection(
        X_numeric, ref, missing_extra, seeds, k, ds["metric"]
    )

    # 3. Imputation comparison
    print(f"  [3/3] Imputation comparison ({len(seeds)} seeds)...")
    imp_rows = run_imputation_comparison(
        X_numeric, k, ds["metric"], seeds
    )

    elapsed = time.time() - t0

    result = {
        "dataset": ds_key,
        "label": ds["label"],
        "n_samples": n,
        "n_markers": p,
        "k": k,
        "metric": ds["metric"],
        "marker_subsampling": marker_rows,
        "missing_injection": missing_rows,
        "imputation_comparison": imp_rows,
        "elapsed_s": round(elapsed, 1),
    }

    # Summary
    def _avg(rows, col):
        vals = [r[col] for r in rows]
        return round(float(np.mean(vals)), 4) if vals else None

    print(f"\n  Marker subsampling (avg jaccard_neighbors per fraction):")
    for frac in marker_fractions:
        subset = [r for r in marker_rows if r["fraction"] == frac]
        avg_jn = _avg(subset, "jaccard_neighbors")
        avg_ss = _avg(subset, "subspace_similarity")
        print(f"    {frac*100:5.0f}%: J_nbr={avg_jn}  subspace_sim={avg_ss}")

    print(f"  Missing injection (avg jaccard_neighbors per rate):")
    for rate in missing_extra:
        subset = [r for r in missing_rows if r["extra_rate"] == rate]
        avg_jn = _avg(subset, "jaccard_neighbors")
        avg_ss = _avg(subset, "subspace_similarity")
        print(f"    +{rate*100:4.0f}%: J_nbr={avg_jn}  subspace_sim={avg_ss}")

    print(f"  Imputation (mode vs median): J_nbr={_avg(imp_rows, 'jaccard_neighbors')}  "
          f"subspace_sim={_avg(imp_rows, 'subspace_similarity')}")
    print(f"  Done in {elapsed:.1f}s")

    # Save
    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)
        with open(out_dir / f"robustness_{ds_key}.json", "w") as f:
            json.dump(result, f, indent=2, ensure_ascii=False)
        print(f"  [saved] {out_dir / f'robustness_{ds_key}.json'}")

    return result


# ╭──────────────────────────────────────────────────────────────╮
# │  CLI                                                        │
# ╰──────────────────────────────────────────────────────────────╯

def main():
    parser = argparse.ArgumentParser(description="Módulo 2 — Robustness Curves")
    parser.add_argument("--dataset", default=None, choices=list(DATASETS.keys()),
                        help="Single dataset (default: all)")
    parser.add_argument("--profile", default="full", choices=list(PROFILES.keys()))
    parser.add_argument("--k", type=int, default=15)
    parser.add_argument("--seeds", nargs="+", type=int, default=None,
                        help="Override seeds from profile")
    parser.add_argument("--out-dir", type=Path,
                        default=Path("experiments/robustness"))
    args = parser.parse_args()

    profile = PROFILES[args.profile]
    seeds = args.seeds or profile["seeds"]
    marker_fracs = profile["marker_fractions"]
    missing_extra = profile["missing_extra"]

    datasets = [args.dataset] if args.dataset else list(DATASETS.keys())
    all_results = []

    for ds_key in datasets:
        result = run_dataset(
            ds_key, marker_fracs, missing_extra, seeds,
            k=args.k, out_dir=args.out_dir,
        )
        all_results.append(result)

    # Save combined
    combined_path = args.out_dir / "robustness_all.json"
    with open(combined_path, "w") as f:
        json.dump(all_results, f, indent=2, ensure_ascii=False)
    print(f"\n[saved] {combined_path}")


if __name__ == "__main__":
    main()
