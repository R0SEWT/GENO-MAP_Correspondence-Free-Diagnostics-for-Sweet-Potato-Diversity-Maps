#!/usr/bin/env python3
"""Ensemble bottleneck embeddings across seeds.

Averages bottleneck vectors from multiple trained autoencoders to produce a
single stable embedding — mitigating the stochastic instability observed in
ae-v1 (ADR-0006). Then re-runs UMAP + kNN on the averaged bottleneck.

Usage::

    python scripts/ensemble_embeddings.py \
        --seeds-dir experiments/global_snp/ae-v2 \
        --out-dir experiments/global_snp/ae-v2/ensemble

    # Or specify seed dirs explicitly:
    python scripts/ensemble_embeddings.py \
        --seed-dirs experiments/global_snp/ae-v2/seed42 \
                    experiments/global_snp/ae-v2/seed52 \
                    experiments/global_snp/ae-v2/seed62 \
        --out-dir experiments/global_snp/ae-v2/ensemble
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np


def load_seed_bottlenecks(seed_dir: Path) -> tuple[list[str], np.ndarray]:
    """Load sample IDs and bottleneck vectors from one seed's nodes JSON."""
    nodes_path = seed_dir / "ae_embedding_nodes.json"
    if not nodes_path.exists():
        raise FileNotFoundError(f"No nodes file: {nodes_path}")

    with open(nodes_path) as f:
        nodes = json.load(f)

    ids = [n["id"] for n in nodes]
    bottlenecks = np.array([n["bottleneck"] for n in nodes], dtype=np.float32)
    return ids, bottlenecks


def ensemble_bottlenecks(
    seed_dirs: list[Path],
) -> tuple[list[str], np.ndarray, dict[str, Any]]:
    """Average bottleneck embeddings across seeds.

    Returns (sample_ids, averaged_bottleneck, stats).
    """
    all_ids: list[list[str]] = []
    all_bn: list[np.ndarray] = []

    for sd in seed_dirs:
        ids, bn = load_seed_bottlenecks(sd)
        all_ids.append(ids)
        all_bn.append(bn)

    # Verify same samples in same order
    ref_ids = all_ids[0]
    for i, ids in enumerate(all_ids[1:], 1):
        if ids != ref_ids:
            # Try reorder
            id_to_idx = {sid: j for j, sid in enumerate(ids)}
            reorder = [id_to_idx[sid] for sid in ref_ids]
            all_bn[i] = all_bn[i][reorder]
            print(f"  [ensemble] Re-ordered seed {i} to match seed 0")

    stacked = np.stack(all_bn, axis=0)  # (n_seeds, n_samples, bottleneck_dim)
    averaged = stacked.mean(axis=0)

    # Compute per-sample std across seeds (measure of consensus)
    per_sample_std = stacked.std(axis=0).mean(axis=1)  # (n_samples,)

    stats = {
        "n_seeds": len(seed_dirs),
        "n_samples": averaged.shape[0],
        "bottleneck_dim": averaged.shape[1],
        "mean_cross_seed_std": float(per_sample_std.mean()),
        "max_cross_seed_std": float(per_sample_std.max()),
        "seed_dirs": [str(d) for d in seed_dirs],
    }

    return ref_ids, averaged, stats


def export_ensemble(
    sample_ids: list[str],
    embeddings: np.ndarray,
    out_dir: Path,
    stats: dict[str, Any],
    k: int = 15,
) -> None:
    """Export ensemble embedding as nodes/edges/stats JSON."""
    from sklearn.neighbors import NearestNeighbors

    # UMAP 2D on averaged bottleneck
    try:
        import umap  # type: ignore
        reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=15,
                            min_dist=0.1, metric="cosine")
        emb_2d = reducer.fit_transform(embeddings)
    except Exception:
        from sklearn.decomposition import PCA
        emb_2d = PCA(n_components=2, random_state=42).fit_transform(embeddings)

    nodes = []
    for idx, sid in enumerate(sample_ids):
        nodes.append({
            "id": str(sid),
            "idx": idx,
            "embedding": emb_2d[idx].tolist(),
            "bottleneck": embeddings[idx].tolist(),
        })

    # kNN
    k_actual = min(k, len(sample_ids) - 1)
    edges: list = []
    if k_actual >= 1:
        nn = NearestNeighbors(n_neighbors=k_actual + 1, metric="cosine")
        nn.fit(embeddings)
        dists, inds = nn.kneighbors(embeddings)
        for i in range(len(sample_ids)):
            for j, d in zip(inds[i][1:], dists[i][1:]):
                edges.append({
                    "source": sample_ids[i],
                    "target": sample_ids[j],
                    "distance": float(d),
                })

    stats["method"] = "autoencoder_ensemble"
    stats["edges"] = len(edges)

    out_dir.mkdir(parents=True, exist_ok=True)
    for suffix, data in [("_nodes.json", nodes), ("_edges.json", edges), ("_stats.json", stats)]:
        fpath = out_dir / f"ae_embedding{suffix}"
        with open(fpath, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
        print(f"  [saved] {fpath}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Ensemble AE embeddings across seeds.")
    parser.add_argument("--seeds-dir", type=Path, default=None,
                        help="Parent dir containing seed* subdirs")
    parser.add_argument("--seed-dirs", type=Path, nargs="+", default=None,
                        help="Explicit list of seed directories")
    parser.add_argument("--out-dir", type=Path, required=True,
                        help="Output directory for ensemble embedding")
    parser.add_argument("-k", type=int, default=15,
                        help="kNN neighbors for graph construction")
    args = parser.parse_args()

    if args.seed_dirs:
        seed_dirs = args.seed_dirs
    elif args.seeds_dir:
        seed_dirs = sorted(args.seeds_dir.glob("seed*"))
    else:
        print("ERROR: Provide --seeds-dir or --seed-dirs")
        sys.exit(1)

    if len(seed_dirs) < 2:
        print(f"ERROR: Need ≥2 seed dirs, found {len(seed_dirs)}")
        sys.exit(1)

    print(f"[ensemble] Averaging {len(seed_dirs)} seeds: {[d.name for d in seed_dirs]}")

    sample_ids, averaged, stats = ensemble_bottlenecks(seed_dirs)
    print(f"[ensemble] Averaged bottleneck: {averaged.shape}")
    print(f"[ensemble] Cross-seed std: mean={stats['mean_cross_seed_std']:.4f}, "
          f"max={stats['max_cross_seed_std']:.4f}")

    export_ensemble(sample_ids, averaged, args.out_dir, stats, k=args.k)
    print(f"\n[ensemble] Done → {args.out_dir}")


if __name__ == "__main__":
    main()
