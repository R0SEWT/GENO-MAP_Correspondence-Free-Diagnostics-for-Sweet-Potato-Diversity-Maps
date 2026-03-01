#!/usr/bin/env python3
"""Utility to turn DArT/DArTSeq genotype tables into embeddings and kNN graphs.

Supported formats:
- sample_columns: first column is row label, columns 2..n are samples (e.g., SNP_Genotypes.csv).
- marker_metrics: leading metadata columns, then sample columns start where header begins with a digit
  (e.g., Report_DSp25-515_* files).

Genotype loading is delegated to ``load_dart.py`` (ADR-0002/0003/0004).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import warnings

import numpy as np
import pandas as pd

# Suppress noisy UMAP / sklearn warnings
warnings.filterwarnings("ignore", message=".*n_jobs.*overridden.*")
warnings.filterwarnings("ignore", message=".*Graph is not fully connected.*")

# --- centralised loader (ADR-0002 / 0003 / 0004) ---
sys.path.insert(0, str(Path(__file__).resolve().parent))
from load_dart import (  # noqa: E402
    check_duplicate_guard,
    filter_missingness,
    load_genotypes,
)
from gpu_utils import (  # noqa: E402
    detect_device,
    print_device_info,
    smart_impute,
    smart_knn,
    smart_pca,
)


def subsample_matrix(
    X: pd.DataFrame, sample_ids: List[str], sample_meta: Dict[str, Dict[str, str]], max_samples: int, max_markers: int, seed: int
) -> Tuple[pd.DataFrame, List[str], Dict[str, Dict[str, str]]]:
    if max_markers and X.shape[1] > max_markers:
        X = X.sample(n=max_markers, axis=1, random_state=seed)
    if max_samples and X.shape[0] > max_samples:
        chosen = X.sample(n=max_samples, axis=0, random_state=seed)
        chosen_ids = chosen.index.tolist()
        X = chosen
        sample_ids = chosen_ids
        sample_meta = {sid: sample_meta.get(sid, {}) for sid in sample_ids}
    return X, sample_ids, sample_meta


def build_embedding(X: np.ndarray, seed: int, device: str = "cpu") -> Tuple[np.ndarray, np.ndarray]:
    max_components = min(50, X.shape[1], X.shape[0])
    if max_components < 2:
        emb = np.column_stack([X[:, 0], np.zeros_like(X[:, 0])])
        return emb, emb

    X_pca, _ = smart_pca(X, max_components, seed, device)
    try:
        import umap  # type: ignore

        reducer = umap.UMAP(n_components=2, random_state=seed)
        embedding = reducer.fit_transform(X_pca)
    except Exception:
        embedding = X_pca[:, :2]
    return embedding, X_pca


def build_graph_features(pca_features: np.ndarray) -> np.ndarray:
    keep = min(30, pca_features.shape[1])
    return pca_features[:, :keep]


def knn_edges(features: np.ndarray, sample_ids: List[str], neighbors: int, metric: str, device: str = "cpu") -> List[Dict[str, object]]:
    edges, _ = smart_knn(features, sample_ids, neighbors, metric, device)
    return edges


def main() -> None:
    parser = argparse.ArgumentParser(description="Build embeddings and kNN graph from genotype tables.")
    parser.add_argument("--input", required=True, help="Path to genotype CSV")
    parser.add_argument("--format", choices=["auto", "sample_columns", "marker_metrics"], default="auto")
    parser.add_argument("--metric", default="euclidean", help="Distance metric for kNN (e.g., euclidean, cosine, jaccard)")
    parser.add_argument("--neighbors", type=int, default=15, help="Number of neighbors for the graph")
    parser.add_argument("--max-samples", type=int, default=500, help="Optional cap on samples for quick runs")
    parser.add_argument("--max-markers", type=int, default=5000, help="Optional cap on markers for quick runs")
    parser.add_argument(
        "--limit-rows",
        type=int,
        default=None,
        help="Limit rows when reading (nrows). Default uses max-markers+10 for speed; set 0 to read all.",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out-prefix", required=True, help="Output prefix (e.g., data/outputs/global_snp)")
    parser.add_argument("--sample-thresh", type=float, default=0.50,
                        help="Drop samples with missingness above this fraction (ADR-0003)")
    parser.add_argument("--marker-thresh", type=float, default=0.50,
                        help="Drop markers with missingness above this fraction (ADR-0003)")
    parser.add_argument("--imputation", default="most_frequent",
                        choices=["most_frequent", "median", "mean"],
                        help="Imputation strategy after filtering")
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"],
                        help="Compute device (auto detects CUDA)")
    args = parser.parse_args()

    path = Path(args.input)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    device = detect_device(args.device)
    print_device_info(device)

    # ADR-0002: warn on duplicate
    if check_duplicate_guard(path):
        print(f"[WARN] {path.name} is the known Wild/LowDensity SNP duplicate (ADR-0002).")
        print("       Consider using the LowDensity copy instead.")

    # --- Load via centralised loader (ADR-0004: BOM, sentinel, NaN) ---
    limit_rows = args.limit_rows
    if limit_rows is None and args.max_markers:
        limit_rows = args.max_markers + 10
    if limit_rows == 0:
        limit_rows = None

    X_df, sample_ids, sample_meta = load_genotypes(
        path, max_rows=limit_rows, fmt=args.format if args.format != "auto" else "auto",
    )

    X_df, sample_ids, sample_meta = subsample_matrix(
        X_df, sample_ids, sample_meta, args.max_samples, args.max_markers, args.seed
    )

    # --- ADR-0003: filter missingness before imputation ---
    missing_rate_raw = float(X_df.isna().sum().sum() / max(X_df.shape[0] * X_df.shape[1], 1))
    X_df, filter_report = filter_missingness(
        X_df,
        sample_thresh=args.sample_thresh,
        marker_thresh=args.marker_thresh,
    )
    # Update sample_ids after potential sample drops
    sample_ids = list(X_df.index.astype(str))
    sample_meta = {sid: sample_meta.get(sid, {}) for sid in sample_ids}
    missing_rate = float(X_df.isna().sum().sum() / max(X_df.shape[0] * X_df.shape[1], 1))

    if filter_report["samples_dropped"] or filter_report["markers_dropped"]:
        print(f"[filter] Dropped {filter_report['samples_dropped']} samples, "
              f"{filter_report['markers_dropped']} markers (ADR-0003)")
        print(f"[filter] Shape: {filter_report['shape_before']} → {filter_report['shape_after']}")

    X_imputed = smart_impute(X_df, strategy=args.imputation, device=device)

    embedding, pca_features = build_embedding(X_imputed, args.seed, device=device)
    graph_features = build_graph_features(pca_features)
    edges = knn_edges(graph_features, sample_ids, args.neighbors, args.metric, device=device)

    nodes = []
    for idx, sid in enumerate(sample_ids):
        meta = sample_meta.get(str(sid), {})
        nodes.append(
            {
                "id": str(sid),
                "idx": idx,
                "embedding": embedding[idx].tolist() if embedding is not None else None,
                "meta": meta,
            }
        )

    stats = {
        "format": args.format,
        "samples": len(sample_ids),
        "markers": int(X_df.shape[1]),
        "missing_rate_raw": missing_rate_raw,
        "missing_rate_after_filter": missing_rate,
        "metric": args.metric,
        "neighbors": args.neighbors,
        "input": str(path),
        "limit_rows": limit_rows,
        "imputation": args.imputation,
        "filter": filter_report,
        "seed": args.seed,
    }

    with (out_prefix.parent / f"{out_prefix.name}_nodes.json").open("w", encoding="utf-8") as f:
        json.dump(nodes, f, ensure_ascii=True, indent=2)
    with (out_prefix.parent / f"{out_prefix.name}_edges.json").open("w", encoding="utf-8") as f:
        json.dump(edges, f, ensure_ascii=True, indent=2)
    with (out_prefix.parent / f"{out_prefix.name}_stats.json").open("w", encoding="utf-8") as f:
        json.dump(stats, f, ensure_ascii=True, indent=2)

    print(f"[done] nodes: {len(nodes)}, edges: {len(edges)}, stats: {stats}")


if __name__ == "__main__":
    main()
