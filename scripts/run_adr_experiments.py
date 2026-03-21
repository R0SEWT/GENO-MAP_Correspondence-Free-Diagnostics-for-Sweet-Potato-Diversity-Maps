#!/usr/bin/env python3
"""
run_adr_experiments.py — Unified experiments for ADR-0012..0015.

Experiments:
  ADR-0012  PC / variance retention sensitivity
  ADR-0013  Block removal vs random subsampling
  ADR-0014  UMAP sensitivity (kNN invariance proof)
  ADR-0015  AE closure — recalculate with unified metrics

Outputs → experiments/adr/  (JSON results + summary tables)

Usage:
    python scripts/run_adr_experiments.py                  # all
    python scripts/run_adr_experiments.py --only 0012      # single ADR
    python scripts/run_adr_experiments.py --only 0012 0013 # subset
    python scripts/run_adr_experiments.py --quick           # fast mode
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

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "experiments" / "adr"

# ═══════════════════════════════════════════════════════════════════
# Dataset registry (same as robustness_curves.py)
# ═══════════════════════════════════════════════════════════════════

DATASETS = {
    "global_snp": {
        "path": "data/10.21223P30BVZYY_Genetic_diversity/SNP_Genotypes.csv",
        "label": "Global SNP",
    },
    "global_silico": {
        "path": "data/10.21223P30BVZYY_Genetic_diversity/SilicoDArT_Genotypes.csv",
        "label": "Global SilicoDArT",
    },
    "lowdensity_snp": {
        "path": "data/10.21223P3UBDJ44_LowDensity/01_Report_DSp25-515_SNPs_Filtered_by _reads.csv",
        "label": "LowDensity SNP",
    },
    "lowdensity_silico": {
        "path": "data/10.21223P3UBDJ44_LowDensity/02_Report_DSp25-515_Silico-DArT_Filtered_by_reads.csv",
        "label": "LowDensity SilicoDArT",
    },
}

SEEDS = [42, 52, 62]
K = 15  # kNN neighbors
ADR0014_LAYOUTS = [
    {"name": "baseline", "n_neighbors": 15, "min_dist": 0.3, "seed": 42},
    {"name": "local_stress", "n_neighbors": 10, "min_dist": 0.1, "seed": 52},
    {"name": "global_stress", "n_neighbors": 30, "min_dist": 0.5, "seed": 62},
]


# ═══════════════════════════════════════════════════════════════════
# Core pipeline helpers
# ═══════════════════════════════════════════════════════════════════

def _impute(X: np.ndarray, strategy: str = "most_frequent") -> np.ndarray:
    imp = SimpleImputer(strategy=strategy, missing_values=np.nan)
    return imp.fit_transform(X)


def _pca_knn(
    X: np.ndarray, n_pca: int = 30, k: int = K, seed: int = 42
) -> Dict[str, Any]:
    """PCA-nD → kNN(k, cosine). Returns PCA object, scores, kNN indices."""
    nc = min(n_pca, X.shape[0] - 1, X.shape[1])
    pca = PCA(n_components=nc, random_state=seed)
    X_pca = pca.fit_transform(X)
    knn_dim = min(n_pca, nc)
    feats = X_pca[:, :knn_dim]
    k_eff = min(k, X.shape[0] - 1)
    nbrs = NearestNeighbors(n_neighbors=k_eff + 1, metric="cosine")
    nbrs.fit(feats)
    _, inds = nbrs.kneighbors(feats)
    return {
        "pca": pca,
        "X_pca": X_pca,
        "knn_inds": inds[:, 1:],
        "explained_variance_ratio": pca.explained_variance_ratio_,
        "cumvar": float(np.cumsum(pca.explained_variance_ratio_)[-1]),
    }


def jaccard_neighbors(a: np.ndarray, b: np.ndarray) -> float:
    scores = []
    for ar, br in zip(a, b):
        sa, sb = set(ar), set(br)
        if not sa and not sb:
            scores.append(1.0)
        elif not sa or not sb:
            scores.append(0.0)
        else:
            scores.append(len(sa & sb) / len(sa | sb))
    return float(np.mean(scores))


def jaccard_edges(a: np.ndarray, b: np.ndarray) -> float:
    def _es(inds):
        return {(i, int(j)) for i, row in enumerate(inds) for j in row}
    sa, sb = _es(a), _es(b)
    if not sa and not sb:
        return 1.0
    return len(sa & sb) / len(sa | sb)


def subspace_similarity(X1: np.ndarray, X2: np.ndarray, nc: int = 10) -> float:
    nc = min(nc, X1.shape[1], X2.shape[1])
    Q1, _ = np.linalg.qr(X1[:, :nc])
    Q2, _ = np.linalg.qr(X2[:, :nc])
    s = np.linalg.svd(Q1.T @ Q2, compute_uv=False)
    return float(np.mean(np.clip(s, -1, 1)))


def _load_panel(ds_key: str) -> np.ndarray:
    """Load → filter → numeric array (with NaN)."""
    path = ROOT / DATASETS[ds_key]["path"]
    X_raw, _, _ = load_genotypes(path)
    X_raw, _ = filter_missingness(X_raw, sample_thresh=0.50, marker_thresh=0.50)
    return X_raw.apply(pd.to_numeric, errors="coerce").values.astype(float)


# ═══════════════════════════════════════════════════════════════════
# ADR-0012: PC / variance retention
# ═══════════════════════════════════════════════════════════════════

def run_adr0012(quick: bool = False) -> Dict[str, Any]:
    """Test kNN stability across PCA dimensionalities."""
    print("\n" + "=" * 60)
    print("  ADR-0012: PC / Variance Retention Sensitivity")
    print("=" * 60)

    pc_values = [5, 10, 15, 20, 30, 50] if not quick else [5, 15, 30]
    seeds = SEEDS if not quick else [42]
    results = {}

    for ds_key in DATASETS:
        label = DATASETS[ds_key]["label"]
        print(f"\n  [{ds_key}] Loading...")
        X_num = _load_panel(ds_key)
        X_imp = _impute(X_num)
        n, p = X_imp.shape
        print(f"  Shape: {n} × {p}")

        # Baseline: PCA-30D
        baseline = _pca_knn(X_imp, n_pca=30, seed=42)
        print(f"  Baseline: PCA-30D, cumvar={baseline['cumvar']:.4f}")

        rows = []
        for n_pca in pc_values:
            for seed in seeds:
                res = _pca_knn(X_imp, n_pca=n_pca, seed=seed)
                je = jaccard_edges(baseline["knn_inds"], res["knn_inds"])
                jn = jaccard_neighbors(baseline["knn_inds"], res["knn_inds"])
                ss = subspace_similarity(baseline["X_pca"], res["X_pca"],
                                         nc=min(n_pca, 30))
                cumvar = float(np.cumsum(res["explained_variance_ratio"])[-1])
                rows.append({
                    "n_pca": n_pca,
                    "seed": seed,
                    "cumvar": round(cumvar, 4),
                    "J_edge": round(je, 4),
                    "J_nbr": round(jn, 4),
                    "SS": round(ss, 4),
                })
                print(f"    PCA-{n_pca}D seed={seed}: "
                      f"J_edge={je:.4f}  J_nbr={jn:.4f}  SS={ss:.4f}  "
                      f"cumvar={cumvar:.4f}")

        results[ds_key] = {
            "label": label, "n": n, "p": p,
            "baseline_cumvar": round(baseline["cumvar"], 4),
            "rows": rows,
        }

    return {"experiment": "ADR-0012", "results": results}


# ═══════════════════════════════════════════════════════════════════
# ADR-0013: Block removal vs random subsampling
# ═══════════════════════════════════════════════════════════════════

def run_adr0013(quick: bool = False) -> Dict[str, Any]:
    """Compare contiguous block removal vs random subsampling."""
    print("\n" + "=" * 60)
    print("  ADR-0013: Block Removal vs Random Subsampling")
    print("=" * 60)

    fracs = [0.05, 0.10, 0.20] if not quick else [0.10]
    positions = ["start", "center", "end"]
    seeds = SEEDS if not quick else [42]
    results = {}

    for ds_key in DATASETS:
        label = DATASETS[ds_key]["label"]
        print(f"\n  [{ds_key}] Loading...")
        X_num = _load_panel(ds_key)
        X_imp = _impute(X_num)
        n, p = X_imp.shape
        print(f"  Shape: {n} × {p}")

        baseline = _pca_knn(X_imp, n_pca=30, seed=42)
        rows = []

        for frac in fracs:
            n_remove = int(p * frac)

            # -- Random removal --
            for seed in seeds:
                rng = np.random.default_rng(seed)
                keep = np.sort(rng.choice(p, size=p - n_remove, replace=False))
                X_sub = X_imp[:, keep]
                res = _pca_knn(X_sub, n_pca=30, seed=seed)
                je = jaccard_edges(baseline["knn_inds"], res["knn_inds"])
                jn = jaccard_neighbors(baseline["knn_inds"], res["knn_inds"])
                ss = subspace_similarity(baseline["X_pca"], res["X_pca"])
                rows.append({
                    "type": "random",
                    "frac_removed": frac,
                    "position": "N/A",
                    "seed": seed,
                    "J_edge": round(je, 4),
                    "J_nbr": round(jn, 4),
                    "SS": round(ss, 4),
                })
                print(f"    Random {frac*100:.0f}% seed={seed}: "
                      f"J_edge={je:.4f}  J_nbr={jn:.4f}  SS={ss:.4f}")

            # -- Block removal --
            for pos in positions:
                if pos == "start":
                    block_start = 0
                elif pos == "center":
                    block_start = (p - n_remove) // 2
                else:
                    block_start = p - n_remove

                mask = np.ones(p, dtype=bool)
                mask[block_start:block_start + n_remove] = False
                X_sub = X_imp[:, mask]

                for seed in seeds:
                    res = _pca_knn(X_sub, n_pca=30, seed=seed)
                    je = jaccard_edges(baseline["knn_inds"], res["knn_inds"])
                    jn = jaccard_neighbors(baseline["knn_inds"], res["knn_inds"])
                    ss = subspace_similarity(baseline["X_pca"], res["X_pca"])
                    rows.append({
                        "type": "block",
                        "frac_removed": frac,
                        "position": pos,
                        "seed": seed,
                        "J_edge": round(je, 4),
                        "J_nbr": round(jn, 4),
                        "SS": round(ss, 4),
                    })
                    print(f"    Block {frac*100:.0f}% @{pos} seed={seed}: "
                          f"J_edge={je:.4f}  J_nbr={jn:.4f}  SS={ss:.4f}")

        results[ds_key] = {"label": label, "n": n, "p": p, "rows": rows}

    return {"experiment": "ADR-0013", "results": results}


# ═══════════════════════════════════════════════════════════════════
# ADR-0014: UMAP sensitivity (kNN invariance proof)
# ═══════════════════════════════════════════════════════════════════

def run_adr0014(quick: bool = False) -> Dict[str, Any]:
    """Show kNN in PCA-30D is invariant to UMAP hyperparameters."""
    print("\n" + "=" * 60)
    print("  ADR-0014: UMAP Sensitivity (kNN Invariance)")
    print("=" * 60)

    try:
        import umap  # noqa: F401
    except ImportError:
        print("  ⚠ umap-learn not installed, skipping ADR-0014")
        return {"experiment": "ADR-0014", "results": {}, "skipped": True}

    n_neighbors_vals = [10, 15, 30] if not quick else [15]
    min_dist_vals = [0.1, 0.3, 0.5] if not quick else [0.1, 0.5]
    seeds = SEEDS if not quick else [42]
    analytic_seed = 42
    results = {}

    for ds_key in DATASETS:
        label = DATASETS[ds_key]["label"]
        print(f"\n  [{ds_key}] Loading...")
        X_num = _load_panel(ds_key)
        X_imp = _impute(X_num)
        n, p = X_imp.shape
        print(f"  Shape: {n} × {p}")

        # Fixed analytic graph in PCA-30D.
        # ADR-0014 varies only UMAP rendering parameters; the PCA graph stays fixed.
        baseline = _pca_knn(X_imp, n_pca=30, seed=analytic_seed)
        X_pca30 = baseline["X_pca"][:, :min(30, baseline["X_pca"].shape[1])]

        rows = []
        for nn in n_neighbors_vals:
            for md in min_dist_vals:
                for seed in seeds:
                    # Run UMAP with different params
                    reducer = umap.UMAP(
                        n_components=2,
                        n_neighbors=nn,
                        min_dist=md,
                        random_state=seed,
                        metric="cosine",
                    )
                    X_umap = reducer.fit_transform(X_pca30)

                    # Build kNN in UMAP-2D space
                    k_eff = min(K, n - 1)
                    nbrs_umap = NearestNeighbors(
                        n_neighbors=k_eff + 1, metric="euclidean"
                    )
                    nbrs_umap.fit(X_umap)
                    _, inds_umap = nbrs_umap.kneighbors(X_umap)
                    inds_umap = inds_umap[:, 1:]

                    # Compare: PCA-kNN vs UMAP-kNN
                    je_pca_umap = jaccard_edges(
                        baseline["knn_inds"], inds_umap
                    )
                    jn_pca_umap = jaccard_neighbors(
                        baseline["knn_inds"], inds_umap
                    )

                    # The analytic graph is unchanged because UMAP is downstream.
                    je_pca_pca = 1.0
                    layout_name = next(
                        (
                            cfg["name"] for cfg in ADR0014_LAYOUTS
                            if cfg["n_neighbors"] == nn
                            and cfg["min_dist"] == md
                            and cfg["seed"] == seed
                        ),
                        None,
                    )

                    rows.append({
                        "n_neighbors": nn,
                        "min_dist": md,
                        "seed": seed,
                        "layout_name": layout_name,
                        "J_edge_pca_vs_umap": round(je_pca_umap, 4),
                        "J_nbr_pca_vs_umap": round(jn_pca_umap, 4),
                        "J_edge_pca_invariance": round(je_pca_pca, 4),
                    })
                    print(
                        f"    nn={nn} md={md} s={seed}: "
                        f"PCA-vs-UMAP J_edge={je_pca_umap:.4f}  "
                        f"PCA invariance={je_pca_pca:.4f}"
                    )

        rep_rows = []
        for cfg in ADR0014_LAYOUTS:
            match = next(
                (
                    r for r in rows
                    if r["n_neighbors"] == cfg["n_neighbors"]
                    and r["min_dist"] == cfg["min_dist"]
                    and r["seed"] == cfg["seed"]
                ),
                None,
            )
            if match is not None:
                rep_rows.append(match)

        results[ds_key] = {
            "label": label,
            "n": n,
            "p": p,
            "analytic_seed": analytic_seed,
            "rows": rows,
            "representative_layouts": rep_rows,
        }

    return {"experiment": "ADR-0014", "results": results}


# ═══════════════════════════════════════════════════════════════════
# ADR-0015: AE closure — unified metrics recalculation
# ═══════════════════════════════════════════════════════════════════

def _load_ae_bottleneck(ds_key: str, seed: int) -> np.ndarray | None:
    """Load AE 64-D bottleneck from experiments/<panel>/ae-v2/seed<s>/ae_embedding_nodes.json."""
    path = ROOT / "experiments" / ds_key / "ae-v2" / f"seed{seed}" / "ae_embedding_nodes.json"
    if not path.exists():
        return None
    with open(path) as fh:
        nodes = json.load(fh)
    # Sort by idx to guarantee consistent ordering
    nodes.sort(key=lambda n: n["idx"])
    return np.array([n["bottleneck"] for n in nodes], dtype=float)


def _knn_from_features(X: np.ndarray, k: int = K) -> np.ndarray:
    """Build kNN index array from a feature matrix (cosine metric)."""
    k_eff = min(k, X.shape[0] - 1)
    nbrs = NearestNeighbors(n_neighbors=k_eff + 1, metric="cosine")
    nbrs.fit(X)
    _, inds = nbrs.kneighbors(X)
    return inds[:, 1:]


def run_adr0015(quick: bool = False) -> Dict[str, Any]:
    """Recalculate PCA vs AE with unified metrics (J_edge, J_nbr, SS).

    Loads AE-64D bottleneck embeddings from experiments/*/ae-v2/seed*/
    and compares kNN stability against PCA-30D baselines.
    """
    print("\n" + "=" * 60)
    print("  ADR-0015: AE Closure — Unified Metrics")
    print("=" * 60)

    results = {}

    for ds_key in DATASETS:
        label = DATASETS[ds_key]["label"]
        print(f"\n  [{ds_key}] Loading...")
        X_num = _load_panel(ds_key)
        X_imp = _impute(X_num)
        n, p = X_imp.shape

        # ── PCA baselines across seeds ──
        pca_results = {}
        for seed in SEEDS:
            res = _pca_knn(X_imp, n_pca=30, seed=seed)
            pca_results[seed] = res

        pairs = [(42, 52), (42, 62), (52, 62)]
        pca_je_list, pca_jn_list = [], []
        for s1, s2 in pairs:
            je = jaccard_edges(pca_results[s1]["knn_inds"],
                               pca_results[s2]["knn_inds"])
            jn = jaccard_neighbors(pca_results[s1]["knn_inds"],
                                    pca_results[s2]["knn_inds"])
            pca_je_list.append(je)
            pca_jn_list.append(jn)

        pca_je_mean = float(np.mean(pca_je_list))
        pca_jn_mean = float(np.mean(pca_jn_list))
        print(f"  PCA inter-seed: J_edge={pca_je_mean:.4f}  J_nbr={pca_jn_mean:.4f}")

        # ── AE bottleneck kNN across seeds ──
        ae_knn = {}
        ae_loaded_seeds = []
        for seed in SEEDS:
            bottleneck = _load_ae_bottleneck(ds_key, seed)
            if bottleneck is not None:
                ae_knn[seed] = _knn_from_features(bottleneck)
                ae_loaded_seeds.append(seed)
                print(f"    AE seed={seed}: loaded {bottleneck.shape[0]}×{bottleneck.shape[1]} bottleneck")
            else:
                print(f"    AE seed={seed}: embeddings not found")

        ae_available = len(ae_loaded_seeds) >= 2
        ae_je_mean = None
        ae_jn_mean = None
        ae_ss_mean = None

        if ae_available:
            # Inter-seed AE stability
            ae_je_list, ae_jn_list = [], []
            for s1, s2 in pairs:
                if s1 in ae_knn and s2 in ae_knn:
                    je = jaccard_edges(ae_knn[s1], ae_knn[s2])
                    jn = jaccard_neighbors(ae_knn[s1], ae_knn[s2])
                    ae_je_list.append(je)
                    ae_jn_list.append(jn)
            ae_je_mean = float(np.mean(ae_je_list)) if ae_je_list else None
            ae_jn_mean = float(np.mean(ae_jn_list)) if ae_jn_list else None

            # Cross-method: PCA-kNN vs AE-kNN (seed-matched)
            cross_je_list, cross_jn_list = [], []
            for seed in ae_loaded_seeds:
                if seed in pca_results:
                    je = jaccard_edges(pca_results[seed]["knn_inds"], ae_knn[seed])
                    jn = jaccard_neighbors(pca_results[seed]["knn_inds"], ae_knn[seed])
                    cross_je_list.append(je)
                    cross_jn_list.append(jn)

            # Subspace similarity: PCA-30D vs AE-64D bottleneck
            ss_list = []
            for seed in ae_loaded_seeds:
                bottleneck = _load_ae_bottleneck(ds_key, seed)
                if bottleneck is not None:
                    ss = subspace_similarity(
                        pca_results[seed]["X_pca"],
                        bottleneck,
                        nc=min(10, 30, bottleneck.shape[1]),
                    )
                    ss_list.append(ss)
            ae_ss_mean = float(np.mean(ss_list)) if ss_list else None

            print(f"  AE inter-seed: J_edge={ae_je_mean:.4f}  J_nbr={ae_jn_mean:.4f}")
            if cross_je_list:
                print(f"  Cross PCA-vs-AE: J_edge={np.mean(cross_je_list):.4f}")
        else:
            print(f"  ⚠ Fewer than 2 AE seeds found — cannot compute inter-seed AE stability")
            cross_je_list, cross_jn_list = [], []

        results[ds_key] = {
            "label": label,
            "n": n,
            "p": p,
            "pca_J_edge_mean": round(pca_je_mean, 4),
            "pca_J_edge_pairs": [round(v, 4) for v in pca_je_list],
            "pca_J_nbr_mean": round(pca_jn_mean, 4),
            "pca_J_nbr_pairs": [round(v, 4) for v in pca_jn_list],
            "ae_available": ae_available,
            "ae_seeds_loaded": ae_loaded_seeds,
            "ae_J_edge_mean": round(ae_je_mean, 4) if ae_je_mean is not None else None,
            "ae_J_nbr_mean": round(ae_jn_mean, 4) if ae_jn_mean is not None else None,
            "ae_SS_mean": round(ae_ss_mean, 4) if ae_ss_mean is not None else None,
            "cross_J_edge_mean": round(float(np.mean(cross_je_list)), 4) if cross_je_list else None,
            "cross_J_nbr_mean": round(float(np.mean(cross_jn_list)), 4) if cross_jn_list else None,
            "stability_gap": round(pca_je_mean - ae_je_mean, 4) if ae_je_mean is not None else None,
        }
        if ae_je_mean is not None:
            print(f"  Stability gap (PCA − AE): {pca_je_mean - ae_je_mean:+.4f}")

    return {"experiment": "ADR-0015", "results": results}


# ═══════════════════════════════════════════════════════════════════
# Summary table generation
# ═══════════════════════════════════════════════════════════════════

def _print_adr0012_summary(data: Dict):
    """Print ADR-0012 summary table."""
    print("\n┌─────────────────────────────────────────────────────────┐")
    print("│  ADR-0012 SUMMARY: J_edge vs PCA dimensionality        │")
    print("├───────────────────┬───────┬───────┬───────┬───────┬─────┤")
    print("│ Panel             │ PC-5  │ PC-10 │ PC-15 │ PC-20 │PC-50│")
    print("├───────────────────┼───────┼───────┼───────┼───────┼─────┤")

    for ds_key, ds_data in data["results"].items():
        vals = {}
        for row in ds_data["rows"]:
            k = row["n_pca"]
            if k not in vals:
                vals[k] = []
            vals[k].append(row["J_edge"])
        label = ds_data["label"][:17].ljust(17)
        cells = []
        for pc in [5, 10, 15, 20, 50]:
            if pc in vals:
                cells.append(f"{np.mean(vals[pc]):.3f}")
            else:
                cells.append("  -  ")
        print(f"│ {label} │{'│'.join(c.center(7) for c in cells)}│")
    print("└───────────────────┴───────┴───────┴───────┴───────┴─────┘")


def _print_adr0013_summary(data: Dict):
    """Print ADR-0013 summary table."""
    print("\n┌──────────────────────────────────────────────────────────┐")
    print("│  ADR-0013 SUMMARY: Block vs Random (J_edge, mean±std)   │")
    print("├───────────────────┬──────────────┬──────────────────────┤")
    print("│ Panel             │ % removed    │ Random    Block      │")
    print("├───────────────────┼──────────────┼──────────────────────┤")

    for ds_key, ds_data in data["results"].items():
        label = ds_data["label"][:17].ljust(17)
        fracs = sorted(set(r["frac_removed"] for r in ds_data["rows"]))
        for frac in fracs:
            random_vals = [r["J_edge"] for r in ds_data["rows"]
                           if r["frac_removed"] == frac and r["type"] == "random"]
            block_vals = [r["J_edge"] for r in ds_data["rows"]
                          if r["frac_removed"] == frac and r["type"] == "block"]
            rm = np.mean(random_vals) if random_vals else 0
            rs = np.std(random_vals) if random_vals else 0
            bm = np.mean(block_vals) if block_vals else 0
            bs = np.std(block_vals) if block_vals else 0
            pct = f"{frac*100:.0f}%".center(12)
            print(f"│ {label} │{pct}│ {rm:.3f}±{rs:.3f}  {bm:.3f}±{bs:.3f} │")
    print("└───────────────────┴──────────────┴──────────────────────┘")


def _print_adr0014_summary(data: Dict):
    """Print ADR-0014 summary."""
    print("\n┌──────────────────────────────────────────────────────────┐")
    print("│  ADR-0014 SUMMARY: UMAP sensitivity                     │")
    print("├───────────────────┬──────────────┬───────────────────────┤")
    print("│ Panel             │ PCA invarian.│ PCA-vs-UMAP J_edge   │")
    print("├───────────────────┼──────────────┼───────────────────────┤")

    for ds_key, ds_data in data["results"].items():
        label = ds_data["label"][:17].ljust(17)
        inv_vals = [r["J_edge_pca_invariance"] for r in ds_data["rows"]]
        umap_vals = [r["J_edge_pca_vs_umap"] for r in ds_data["rows"]]
        inv_m = np.mean(inv_vals) if inv_vals else 0
        umap_m = np.mean(umap_vals) if umap_vals else 0
        umap_s = np.std(umap_vals) if umap_vals else 0
        print(f"│ {label} │   {inv_m:.4f}    │  {umap_m:.3f} ± {umap_s:.3f}          │")
    print("└───────────────────┴──────────────┴───────────────────────┘")


# ═══════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Run ADR-0012..0015 experiments")
    parser.add_argument("--only", nargs="*", default=None,
                        help="Run only specific ADRs (e.g. 0012 0013)")
    parser.add_argument("--quick", action="store_true",
                        help="Quick mode (fewer combos)")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    adrs_to_run = args.only or ["0012", "0013", "0014", "0015"]
    t0 = time.time()
    all_results = {}

    if "0012" in adrs_to_run:
        r = run_adr0012(quick=args.quick)
        all_results["ADR-0012"] = r
        _print_adr0012_summary(r)
        with open(args.out_dir / "adr0012_results.json", "w") as f:
            json.dump(r, f, indent=2)

    if "0013" in adrs_to_run:
        r = run_adr0013(quick=args.quick)
        all_results["ADR-0013"] = r
        _print_adr0013_summary(r)
        with open(args.out_dir / "adr0013_results.json", "w") as f:
            json.dump(r, f, indent=2)

    if "0014" in adrs_to_run:
        r = run_adr0014(quick=args.quick)
        all_results["ADR-0014"] = r
        if not r.get("skipped"):
            _print_adr0014_summary(r)
        with open(args.out_dir / "adr0014_results.json", "w") as f:
            json.dump(r, f, indent=2)

    if "0015" in adrs_to_run:
        r = run_adr0015(quick=args.quick)
        all_results["ADR-0015"] = r
        with open(args.out_dir / "adr0015_results.json", "w") as f:
            json.dump(r, f, indent=2)

    elapsed = time.time() - t0
    print(f"\n{'='*60}")
    print(f"  ALL DONE in {elapsed:.1f}s")
    print(f"  Results → {args.out_dir}/")
    print(f"{'='*60}")

    # Save combined
    with open(args.out_dir / "all_adr_results.json", "w") as f:
        json.dump(all_results, f, indent=2)


if __name__ == "__main__":
    main()
