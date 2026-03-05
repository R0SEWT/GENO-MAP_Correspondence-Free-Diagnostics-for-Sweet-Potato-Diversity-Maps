#!/usr/bin/env python3
"""Validate embedding quality, stability, imputation sensitivity, and PCA loadings.

Produces a ``*_validation.json`` file with all metrics, suitable for
``plot_results.py`` and ``short_paper_v2``.

Metrics implemented (matching paper v1 §2.4):
  1. Trustworthiness  (sklearn)
  2. Multi-seed stability  (Jaccard of kNN neighbourhoods + global edge overlap)
  3. Imputation sensitivity  (mode vs median, delta trustworthiness + Jaccard)
  4. PCA interpretability  (top-10 marker loadings, cumulative variance)
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import warnings

import numpy as np
import pandas as pd
from sklearn.manifold import trustworthiness

# Suppress noisy UMAP / sklearn warnings
warnings.filterwarnings("ignore", message=".*n_jobs.*overridden.*")
warnings.filterwarnings("ignore", message=".*Graph is not fully connected.*")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from load_dart import filter_missingness, load_genotypes  # noqa: E402
from gpu_utils import (  # noqa: E402
    detect_device,
    print_device_info,
    smart_impute,
    smart_knn,
    smart_pca,
)

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _jaccard_neighbours(
    ids_a: List[List[str]], ids_b: List[List[str]]
) -> float:
    """Mean per-node Jaccard of kNN neighbour sets."""
    scores = []
    for a, b in zip(ids_a, ids_b):
        sa, sb = set(a), set(b)
        if not sa and not sb:
            scores.append(1.0)
        elif not sa or not sb:
            scores.append(0.0)
        else:
            scores.append(len(sa & sb) / len(sa | sb))
    return float(np.mean(scores))


def _jaccard_edges(
    edges_a: List[Dict[str, Any]], edges_b: List[Dict[str, Any]]
) -> float:
    """Jaccard overlap of global edge sets."""
    def _edge_set(edges: List[Dict[str, Any]]) -> set:
        return {(e["source"], e["target"]) for e in edges}
    sa, sb = _edge_set(edges_a), _edge_set(edges_b)
    if not sa and not sb:
        return 1.0
    return len(sa & sb) / len(sa | sb)


def _build_knn(features: np.ndarray, ids: List[str], k: int, metric: str, device: str = "cpu"):
    """Return (edges_list, neighbour_ids_per_node).  GPU-accelerated when device=cuda."""
    return smart_knn(features, ids, k, metric, device)


def _run_pipeline(
    X_raw: pd.DataFrame,
    imputation: str,
    seed: int,
    k: int,
    metric: str,
    n_pca: int = 50,
    device: str = "cpu",
) -> Dict[str, Any]:
    """Full pipeline: impute → PCA → UMAP → kNN.  Returns artefacts."""
    X_imp = smart_impute(X_raw, strategy=imputation, device=device)

    n_comp = min(n_pca, X_imp.shape[0], X_imp.shape[1])
    X_pca, pca = smart_pca(X_imp, n_comp, seed, device)

    # UMAP (CPU — no cuML dependency)
    try:
        import umap  # type: ignore
        reducer = umap.UMAP(n_components=2, random_state=seed)
        embedding_2d = reducer.fit_transform(X_pca)
    except Exception:
        embedding_2d = X_pca[:, :2] if X_pca.shape[1] >= 2 else np.column_stack([X_pca[:, 0], np.zeros(len(X_pca))])

    ids = [str(s) for s in X_raw.index]
    graph_feat = X_pca[:, :min(30, X_pca.shape[1])]
    edges, nbrs = smart_knn(graph_feat, ids, k, metric, device)

    return {
        "embedding_2d": embedding_2d,
        "X_pca": X_pca,
        "pca": pca,
        "edges": edges,
        "neighbours": nbrs,
        "ids": ids,
        "X_imputed": X_imp,
    }


# ---------------------------------------------------------------------------
# metric functions
# ---------------------------------------------------------------------------


def compute_trustworthiness(
    X_high: np.ndarray, X_low: np.ndarray, k: int = 15
) -> float:
    """Trustworthiness score (Venna & Kaski, 2006)."""
    k_eff = min(k, X_high.shape[0] - 2)
    if k_eff < 1:
        return float("nan")
    return float(trustworthiness(X_high, X_low, n_neighbors=k_eff))


def compute_stability(
    X_raw: pd.DataFrame,
    seeds: List[int],
    k: int,
    metric: str,
    imputation: str = "most_frequent",
    device: str = "cpu",
) -> Dict[str, Any]:
    """Multi-seed stability: Jaccard of PCA-space kNN and UMAP-space kNN."""
    runs = []
    for s in seeds:
        runs.append(_run_pipeline(X_raw, imputation, s, k, metric, device=device))

    # Pairwise Jaccard of neighbour sets in PCA-space
    pca_jaccards = []
    umap_jaccards = []
    edge_jaccards = []
    for i in range(len(runs)):
        for j in range(i + 1, len(runs)):
            pca_jaccards.append(_jaccard_neighbours(runs[i]["neighbours"], runs[j]["neighbours"]))
            # UMAP-space kNN
            k_eff_u = min(k, len(runs[i]["ids"]) - 1)
            if k_eff_u >= 1:
                _, nbrs_i = smart_knn(runs[i]["embedding_2d"], runs[i]["ids"], k, metric, device)
                _, nbrs_j = smart_knn(runs[j]["embedding_2d"], runs[j]["ids"], k, metric, device)
                umap_jaccards.append(_jaccard_neighbours(nbrs_i, nbrs_j))
            edge_jaccards.append(_jaccard_edges(runs[i]["edges"], runs[j]["edges"]))

    return {
        "seeds": seeds,
        "pca_neighbour_jaccard_mean": float(np.mean(pca_jaccards)) if pca_jaccards else None,
        "pca_neighbour_jaccard_std": float(np.std(pca_jaccards)) if pca_jaccards else None,
        "umap_neighbour_jaccard_mean": float(np.mean(umap_jaccards)) if umap_jaccards else None,
        "umap_neighbour_jaccard_std": float(np.std(umap_jaccards)) if umap_jaccards else None,
        "edge_jaccard_mean": float(np.mean(edge_jaccards)) if edge_jaccards else None,
        "edge_jaccard_std": float(np.std(edge_jaccards)) if edge_jaccards else None,
    }


def compute_imputation_sensitivity(
    X_raw: pd.DataFrame,
    seed: int,
    k: int,
    metric: str,
    strategies: List[str] = ("most_frequent", "median"),
    device: str = "cpu",
) -> Dict[str, Any]:
    """Compare graph topology under different imputation strategies."""
    results: Dict[str, Any] = {"strategies": list(strategies)}
    runs = {}
    for strat in strategies:
        runs[strat] = _run_pipeline(X_raw, strat, seed, k, metric, device=device)

    baseline = strategies[0]
    for strat in strategies[1:]:
        nbr_j = _jaccard_neighbours(runs[baseline]["neighbours"], runs[strat]["neighbours"])
        edge_j = _jaccard_edges(runs[baseline]["edges"], runs[strat]["edges"])
        trust_base = compute_trustworthiness(runs[baseline]["X_imputed"], runs[baseline]["embedding_2d"], k)
        trust_alt = compute_trustworthiness(runs[strat]["X_imputed"], runs[strat]["embedding_2d"], k)
        results[f"{baseline}_vs_{strat}"] = {
            "neighbour_jaccard": nbr_j,
            "edge_jaccard": edge_j,
            "trustworthiness_baseline": trust_base,
            "trustworthiness_alt": trust_alt,
            "trustworthiness_delta": trust_alt - trust_base,
        }
    return results


def compute_pca_loadings(
    pca_model: PCA, marker_names: List[str], top_n: int = 10
) -> Dict[str, Any]:
    """Extract top-N marker loadings for PC1 and PC2, plus cumulative variance."""
    result: Dict[str, Any] = {}
    var_explained = pca_model.explained_variance_ratio_.tolist()
    result["variance_explained"] = var_explained
    result["cumulative_variance_5pc"] = float(sum(var_explained[:5]))

    components = pca_model.components_
    for pc_idx, pc_name in enumerate(["PC1", "PC2"]):
        if pc_idx >= components.shape[0]:
            break
        loadings = components[pc_idx]
        top_idx = np.argsort(np.abs(loadings))[::-1][:top_n]
        result[f"{pc_name}_top_markers"] = [
            {"marker": marker_names[i] if i < len(marker_names) else str(i),
             "loading": float(loadings[i])}
            for i in top_idx
        ]
    return result


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate embeddings and kNN graphs.")
    parser.add_argument("--input", required=True, type=Path, help="Genotype CSV")
    parser.add_argument("--out-prefix", required=True, type=Path, help="Output prefix")
    parser.add_argument("--metric", default="euclidean")
    parser.add_argument("--neighbors", type=int, default=15)
    parser.add_argument("--seeds", default="42,52,62",
                        help="Comma-separated seeds for stability")
    parser.add_argument("--sample-thresh", type=float, default=0.50)
    parser.add_argument("--marker-thresh", type=float, default=0.50)
    parser.add_argument("--max-samples", type=int, default=0,
                        help="Cap samples (0 = no limit)")
    parser.add_argument("--max-markers", type=int, default=0,
                        help="Cap markers (0 = no limit)")
    parser.add_argument("--imputation", default="most_frequent")
    parser.add_argument("--skip-stability", action="store_true",
                        help="Skip multi-seed stability (saves time)")
    parser.add_argument("--skip-sensitivity", action="store_true",
                        help="Skip imputation sensitivity (saves time)")
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"],
                        help="Compute device (auto detects CUDA)")
    args = parser.parse_args()

    seeds = [int(s) for s in args.seeds.split(",")]
    out_prefix = args.out_prefix
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    device = detect_device(args.device)
    print_device_info(device)

    t0 = time.time()
    print(f"[validate] Loading {args.input} ...")

    # Load and filter
    X_raw, sample_ids, _ = load_genotypes(args.input)
    if args.max_samples and X_raw.shape[0] > args.max_samples:
        X_raw = X_raw.sample(n=args.max_samples, random_state=seeds[0])
    if args.max_markers and X_raw.shape[1] > args.max_markers:
        X_raw = X_raw.sample(n=args.max_markers, axis=1, random_state=seeds[0])
    X_raw, filter_report = filter_missingness(
        X_raw, sample_thresh=args.sample_thresh, marker_thresh=args.marker_thresh,
    )
    sample_ids = list(X_raw.index.astype(str))
    marker_names = list(X_raw.columns.astype(str))
    print(f"[validate] Shape after filter: {X_raw.shape}")

    validation: Dict[str, Any] = {
        "input": str(args.input),
        "shape": list(X_raw.shape),
        "metric": args.metric,
        "neighbors": args.neighbors,
        "seeds": seeds,
        "imputation": args.imputation,
        "filter": filter_report,
    }

    # 1. Trustworthiness (primary seed)
    print("[validate] Computing trustworthiness ...")
    run0 = _run_pipeline(X_raw, args.imputation, seeds[0], args.neighbors, args.metric, device=device)
    trust_vals = []
    for s in seeds:
        r = _run_pipeline(X_raw, args.imputation, s, args.neighbors, args.metric, device=device)
        t = compute_trustworthiness(r["X_imputed"], r["embedding_2d"], args.neighbors)
        trust_vals.append(t)
    validation["trustworthiness"] = {
        "mean": float(np.mean(trust_vals)),
        "std": float(np.std(trust_vals)),
        "per_seed": dict(zip(seeds, trust_vals)),
    }
    print(f"    trust = {validation['trustworthiness']['mean']:.4f} ± {validation['trustworthiness']['std']:.4f}")

    # 2. Multi-seed stability
    if not args.skip_stability and len(seeds) > 1:
        print("[validate] Computing stability ...")
        validation["stability"] = compute_stability(
            X_raw, seeds, args.neighbors, args.metric, args.imputation, device=device,
        )
        print(f"    PCA-Jaccard = {validation['stability']['pca_neighbour_jaccard_mean']:.4f}")
    else:
        validation["stability"] = None

    # 3. Imputation sensitivity
    if not args.skip_sensitivity:
        print("[validate] Computing imputation sensitivity ...")
        validation["imputation_sensitivity"] = compute_imputation_sensitivity(
            X_raw, seeds[0], args.neighbors, args.metric, device=device,
        )
        for key, val in validation["imputation_sensitivity"].items():
            if isinstance(val, dict) and "neighbour_jaccard" in val:
                print(f"    {key}: neighbour_jaccard={val['neighbour_jaccard']:.4f}")
    else:
        validation["imputation_sensitivity"] = None

    # 4. PCA loadings
    print("[validate] Extracting PCA loadings ...")
    validation["pca_loadings"] = compute_pca_loadings(run0["pca"], marker_names)
    print(f"    cumulative variance (PC1-5) = {validation['pca_loadings']['cumulative_variance_5pc']:.4f}")

    elapsed = time.time() - t0
    validation["elapsed_seconds"] = round(elapsed, 1)

    # Save
    out_file = out_prefix.parent / f"{out_prefix.name}_validation.json"
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(validation, f, ensure_ascii=False, indent=2)
    print(f"[validate] Saved → {out_file}  ({elapsed:.1f}s)")


if __name__ == "__main__":
    main()
