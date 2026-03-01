#!/usr/bin/env python3
"""Utility to turn DArT/DArTSeq genotype tables into embeddings and kNN graphs.

Supported formats:
- sample_columns: first column is row label, columns 2..n are samples (e.g., SNP_Genotypes.csv).
- marker_metrics: leading metadata columns, then sample columns start where header begins with a digit
  (e.g., Report_DSp25-515_* files).
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.neighbors import NearestNeighbors


def guess_separator(path: Path) -> str:
    line = path.open(encoding="utf-8-sig").readline()
    return ";" if line.count(";") > line.count(",") else ","


def detect_format(header: List[str]) -> str:
    if header and header[0] == "Sample_code":
        return "sample_columns"
    if "AlleleID" in header:
        return "marker_metrics"
    return "sample_columns"


def load_sample_columns(
    path: Path, sep: str, nrows: int | None
) -> Tuple[pd.DataFrame, List[str], Dict[str, Dict[str, str]]]:
    df = pd.read_csv(path, sep=sep, low_memory=False, nrows=nrows)
    sample_ids = [str(c) for c in df.columns[1:]]
    first_col = df.columns[0]
    labels = df[first_col].astype(str)
    marker_mask = labels.str.match(r"^\\d") | labels.str.contains("|", regex=False)
    meta_rows = df.loc[~marker_mask]
    marker_rows = df.loc[marker_mask]

    sample_meta: Dict[str, Dict[str, str]] = {sid: {} for sid in sample_ids}
    for _, row in meta_rows.iterrows():
        key = str(row.iloc[0])
        for sid, val in zip(sample_ids, row.iloc[1:]):
            if pd.notna(val) and str(val).strip():
                sample_meta[sid][key] = str(val)

    geno = marker_rows.set_index(first_col)
    geno = geno.replace({"-": np.nan, "": np.nan})
    geno = geno.apply(pd.to_numeric, errors="coerce")
    X = geno.T  # samples x markers
    return X, sample_ids, sample_meta


def load_marker_metrics(
    path: Path, sep: str, nrows: int | None
) -> Tuple[pd.DataFrame, List[str], Dict[str, Dict[str, str]]]:
    df = pd.read_csv(path, sep=sep, low_memory=False, nrows=nrows)
    header = [str(c) for c in df.columns]
    sample_start = next(i for i, name in enumerate(header) if name[0].isdigit())
    sample_cols = header[sample_start:]
    marker_id_col = header[0]
    geno = df[sample_cols]
    geno = geno.replace({"-": np.nan, "": np.nan})
    geno = geno.apply(pd.to_numeric, errors="coerce")
    geno.index = df[marker_id_col].astype(str)
    X = geno.T  # samples x markers
    sample_meta: Dict[str, Dict[str, str]] = {str(s): {} for s in sample_cols}
    return X, sample_cols, sample_meta


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


def build_embedding(X: np.ndarray, seed: int) -> Tuple[np.ndarray, np.ndarray]:
    max_components = min(50, X.shape[1], X.shape[0])
    if max_components < 2:
        # Degenerate case: only one feature or sample; pad with zeros.
        emb = np.column_stack([X[:, 0], np.zeros_like(X[:, 0])])
        return emb, emb

    pca = PCA(n_components=max_components, random_state=seed)
    X_pca = pca.fit_transform(X)
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


def knn_edges(features: np.ndarray, sample_ids: List[str], neighbors: int, metric: str) -> List[Dict[str, object]]:
    k = min(neighbors, len(sample_ids) - 1)
    if k < 1:
        return []
    nn = NearestNeighbors(n_neighbors=k + 1, metric=metric)
    nn.fit(features)
    distances, indices = nn.kneighbors(features)
    edges: List[Dict[str, object]] = []
    for i, (row_idx, row_dist) in enumerate(zip(indices, distances)):
        for j, dist in zip(row_idx[1:], row_dist[1:]):  # skip self
            edges.append({"source": sample_ids[i], "target": sample_ids[j], "distance": float(dist)})
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
    args = parser.parse_args()

    path = Path(args.input)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    sep = guess_separator(path)
    header_line = path.open(encoding="utf-8-sig").readline().strip()
    header = re.split(sep, header_line)
    fmt = detect_format(header) if args.format == "auto" else args.format

    limit_rows = args.limit_rows
    if limit_rows is None and args.max_markers:
        limit_rows = args.max_markers + 10
    if limit_rows == 0:
        limit_rows = None

    if fmt == "sample_columns":
        X_df, sample_ids, sample_meta = load_sample_columns(path, sep, limit_rows)
    else:
        X_df, sample_ids, sample_meta = load_marker_metrics(path, sep, limit_rows)

    X_df, sample_ids, sample_meta = subsample_matrix(
        X_df, sample_ids, sample_meta, args.max_samples, args.max_markers, args.seed
    )
    missing_rate = float(pd.isna(X_df).sum().sum() / (X_df.shape[0] * X_df.shape[1]))

    imputer = SimpleImputer(strategy="most_frequent")
    X_imputed = imputer.fit_transform(X_df)

    embedding, pca_features = build_embedding(X_imputed, args.seed)
    graph_features = build_graph_features(pca_features)
    edges = knn_edges(graph_features, sample_ids, args.neighbors, args.metric)

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
        "format": fmt,
        "samples": len(sample_ids),
        "markers": int(X_df.shape[1]),
        "missing_rate": missing_rate,
        "metric": args.metric,
        "neighbors": args.neighbors,
        "input": str(path),
        "separator": sep,
        "limit_rows": limit_rows,
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
