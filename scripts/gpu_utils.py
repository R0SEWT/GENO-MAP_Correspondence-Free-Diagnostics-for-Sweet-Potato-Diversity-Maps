#!/usr/bin/env python3
"""GPU-accelerated utilities for the OXOR genomic pipeline.

Provides drop-in replacements for sklearn PCA, kNN, and imputation
using PyTorch CUDA.  Falls back transparently to CPU / sklearn when
CUDA is unavailable.

Usage::

    from gpu_utils import detect_device, print_device_info
    from gpu_utils import smart_pca, smart_knn, smart_impute

    device = detect_device("auto")   # → "cuda" or "cpu"
    print_device_info(device)

    X_imp = smart_impute(X_df, "most_frequent", device)
    X_pca, pca_obj = smart_pca(X_imp, 50, seed=42, device=device)
    edges, nbrs = smart_knn(X_pca[:, :30], ids, k=20, metric="cosine", device=device)
"""
from __future__ import annotations

import warnings
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Conditional PyTorch import
# ---------------------------------------------------------------------------

_HAS_TORCH = False
try:
    import torch
    import torch.nn.functional as F

    _HAS_TORCH = True
except ImportError:
    pass


# ---------------------------------------------------------------------------
# Device helpers
# ---------------------------------------------------------------------------


def detect_device(preference: str = "auto") -> str:
    """Resolve device preference to ``'cuda'`` or ``'cpu'``."""
    if preference == "auto":
        if _HAS_TORCH and torch.cuda.is_available():
            return "cuda"
        return "cpu"
    if preference == "cuda":
        if not _HAS_TORCH or not torch.cuda.is_available():
            warnings.warn("CUDA requested but not available — falling back to CPU")
            return "cpu"
    return preference


def print_device_info(device: str) -> None:
    """Print GPU model and VRAM when running on CUDA."""
    if device.startswith("cuda") and _HAS_TORCH and torch.cuda.is_available():
        name = torch.cuda.get_device_name(0)
        vram_gb = torch.cuda.get_device_properties(0).total_memory / 1e9
        print(f"[gpu] {name}  ({vram_gb:.1f} GB VRAM)")
    else:
        print(f"[device] {device}")


# ---------------------------------------------------------------------------
# PCA  (torch.pca_lowrank  →  sklearn-compatible interface)
# ---------------------------------------------------------------------------


class TorchPCA:
    """GPU-accelerated PCA via ``torch.pca_lowrank``.

    Exposes ``components_`` and ``explained_variance_ratio_`` so that
    downstream code (e.g. ``compute_pca_loadings``) works unchanged.
    """

    def __init__(self, n_components: int, random_state: int = 42):
        self.n_components = n_components
        self.random_state = random_state
        self.components_: np.ndarray | None = None
        self.explained_variance_ratio_: np.ndarray | None = None

    def fit_transform(self, X: np.ndarray, device: str = "cuda") -> np.ndarray:
        torch.manual_seed(self.random_state)
        X_t = torch.tensor(X, dtype=torch.float32, device=device)
        mean = X_t.mean(dim=0)
        X_c = X_t - mean

        q = min(self.n_components, X_c.shape[0] - 1, X_c.shape[1])
        U, S, V = torch.pca_lowrank(X_c, q=q, center=False, niter=4)

        X_pca = (U[:, :q] * S[:q]).cpu().numpy()

        total_ss = (X_c**2).sum().item()
        self.explained_variance_ratio_ = (S[:q] ** 2 / total_ss).cpu().numpy()
        self.components_ = V[:, :q].T.cpu().numpy()  # (q, n_features)

        return X_pca


def smart_pca(
    X: np.ndarray,
    n_components: int,
    seed: int,
    device: str,
) -> Tuple[np.ndarray, Any]:
    """PCA with auto-backend.  Returns ``(X_pca, pca_object)``."""
    n_comp = min(n_components, X.shape[0] - 1, X.shape[1])
    if n_comp < 1:
        return X[:, :1].copy(), None

    if device.startswith("cuda"):
        pca = TorchPCA(n_comp, random_state=seed)
        return pca.fit_transform(X, device), pca

    from sklearn.decomposition import PCA

    pca = PCA(n_components=n_comp, random_state=seed)
    return pca.fit_transform(X), pca


# ---------------------------------------------------------------------------
# kNN  (torch.cdist / cosine)
# ---------------------------------------------------------------------------


def smart_knn(
    features: np.ndarray,
    sample_ids: List[str],
    k: int,
    metric: str,
    device: str,
) -> Tuple[List[Dict[str, Any]], List[List[str]]]:
    """kNN graph.  Returns ``(edges, neighbour_ids_per_node)``."""
    k_eff = min(k, len(sample_ids) - 1)
    if k_eff < 1:
        return [], [[] for _ in sample_ids]

    if device.startswith("cuda"):
        X = torch.tensor(features, dtype=torch.float32, device=device)
        if metric == "cosine":
            X_n = F.normalize(X, dim=1)
            dist = 1.0 - (X_n @ X_n.T)
        else:
            dist = torch.cdist(X, X)
        dist.fill_diagonal_(float("inf"))
        top_d, top_i = dist.topk(k_eff, largest=False, dim=1)
        dists_np = top_d.cpu().numpy()
        inds_np = top_i.cpu().numpy()
    else:
        from sklearn.neighbors import NearestNeighbors

        nn = NearestNeighbors(n_neighbors=k_eff + 1, metric=metric)
        nn.fit(features)
        d_full, i_full = nn.kneighbors(features)
        dists_np = d_full[:, 1:]  # skip self
        inds_np = i_full[:, 1:]

    edges: List[Dict[str, Any]] = []
    neighbours: List[List[str]] = []
    for i in range(len(sample_ids)):
        nbrs = [sample_ids[int(j)] for j in inds_np[i]]
        neighbours.append(nbrs)
        for j_idx, d in zip(inds_np[i], dists_np[i]):
            edges.append(
                {
                    "source": sample_ids[i],
                    "target": sample_ids[int(j_idx)],
                    "distance": float(d),
                }
            )
    return edges, neighbours


# ---------------------------------------------------------------------------
# Imputation
# ---------------------------------------------------------------------------


def smart_impute(
    X: np.ndarray | pd.DataFrame,
    strategy: str = "most_frequent",
    device: str = "cpu",
) -> np.ndarray:
    """Impute NaN values.  GPU-accelerated for CUDA; sklearn for CPU."""
    X_np = (
        X.values.astype(np.float32)
        if isinstance(X, pd.DataFrame)
        else np.asarray(X, dtype=np.float32)
    )

    if device.startswith("cuda"):
        X_t = torch.tensor(X_np, dtype=torch.float32, device=device)
        nan_mask = torch.isnan(X_t)
        if not nan_mask.any():
            return X_t.cpu().numpy()

        if strategy == "mean":
            fill = torch.nanmean(X_t, dim=0)
        elif strategy == "median":
            fill = torch.nanmedian(X_t, dim=0).values
        else:  # most_frequent — optimised for genotype values {0, 1, 2}
            X_safe = X_t.clone()
            X_safe[nan_mask] = -1.0
            counts = torch.stack(
                [(X_safe == float(v)).sum(dim=0) for v in range(3)]
            )
            fill = counts.argmax(dim=0).float()
            # Columns where every value is NaN → fill with 0
            all_nan = nan_mask.all(dim=0)
            fill[all_nan] = 0.0

        X_t = torch.where(nan_mask, fill.unsqueeze(0), X_t)
        return X_t.cpu().numpy()

    from sklearn.impute import SimpleImputer

    return SimpleImputer(strategy=strategy).fit_transform(X_np)
