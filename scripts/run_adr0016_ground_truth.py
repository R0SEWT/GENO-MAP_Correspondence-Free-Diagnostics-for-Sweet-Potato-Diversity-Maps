#!/usr/bin/env python3
"""
run_adr0016_ground_truth.py — ADR-0016: Ground-truth validation.

Global SNP and Global SilicoDArT share 5 930 / 5 970 CIP accession
codes.  This lets us compare two independently-built kNN graphs
**with** sample-level correspondence (ground truth) and verify
that the correspondence-free metrics (J_edge, J_nbr, SS) yield
conclusions consistent with correspondence-based metrics.

The experiment proceeds as follows:

  1.  Load Global SNP and Global SilicoDArT.
  2.  Restrict both to the intersection of shared CIP accessions
      (ground truth alignment).
  3.  Build PCA→kNN graphs for each panel independently.
  4.  Compute **correspondence-free** metrics:
        – J_edge  (global edge overlap by node index)
        – J_nbr   (per-node kNN neighbour Jaccard by index)
        → These treat nodes as if alignment is unknown.
  5.  Compute **correspondence-based** metrics (using the shared IDs):
        – J_edge_gt  (edges compared after re-indexing to shared ID order)
        – J_nbr_gt   (per-accession neighbour overlap in accession space)
  6.  Report both sets of metrics and the gap.

  Additionally:
  7.  Same exercise for LD SNP ↔ LD SilicoDArT (476 shared IDs).

Usage:
    python scripts/run_adr0016_ground_truth.py
"""
from __future__ import annotations

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

K = 15
N_PCA = 30
SEEDS = [42, 52, 62]


# ═══════════════════════════════════════════════════════════════════
# Dataset pairs with shared ID namespaces
# ═══════════════════════════════════════════════════════════════════

PAIRS = {
    "global": {
        "panel_a": {
            "key": "global_snp",
            "path": "data/10.21223P30BVZYY_Genetic_diversity/SNP_Genotypes.csv",
            "label": "Global SNP",
        },
        "panel_b": {
            "key": "global_silico",
            "path": "data/10.21223P30BVZYY_Genetic_diversity/SilicoDArT_Genotypes.csv",
            "label": "Global SilicoDArT",
        },
    },
    "lowdensity": {
        "panel_a": {
            "key": "lowdensity_snp",
            "path": "data/10.21223P3UBDJ44_LowDensity/01_Report_DSp25-515_SNPs_Filtered_by _reads.csv",
            "label": "LowDensity SNP",
        },
        "panel_b": {
            "key": "lowdensity_silico",
            "path": "data/10.21223P3UBDJ44_LowDensity/02_Report_DSp25-515_Silico-DArT_Filtered_by_reads.csv",
            "label": "LowDensity SilicoDArT",
        },
    },
}


# ═══════════════════════════════════════════════════════════════════
# Pipeline helpers
# ═══════════════════════════════════════════════════════════════════

def _impute(X: np.ndarray) -> np.ndarray:
    imp = SimpleImputer(strategy="most_frequent", missing_values=np.nan)
    return imp.fit_transform(X)


def _pca_knn(
    X: np.ndarray, n_pca: int = N_PCA, k: int = K, seed: int = 42,
) -> Dict[str, Any]:
    """PCA-nD → kNN(k, cosine). Returns PCA scores + kNN indices."""
    nc = min(n_pca, X.shape[0] - 1, X.shape[1])
    pca = PCA(n_components=nc, random_state=seed)
    X_pca = pca.fit_transform(X)
    feats = X_pca[:, : min(n_pca, nc)]
    k_eff = min(k, X.shape[0] - 1)
    nbrs = NearestNeighbors(n_neighbors=k_eff + 1, metric="cosine")
    nbrs.fit(feats)
    _, inds = nbrs.kneighbors(feats)
    return {
        "X_pca": X_pca,
        "knn_inds": inds[:, 1:],   # exclude self
        "cumvar": float(np.cumsum(pca.explained_variance_ratio_)[-1]),
    }


# ═══════════════════════════════════════════════════════════════════
# Correspondence-free metrics (index-based, no alignment)
# ═══════════════════════════════════════════════════════════════════

def jaccard_neighbors(a: np.ndarray, b: np.ndarray) -> float:
    """Mean per-node Jaccard of kNN neighbour index sets."""
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
    """Jaccard of global directed edge sets (index-based)."""
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


# ═══════════════════════════════════════════════════════════════════
# Correspondence-based metrics (accession-aligned)
# ═══════════════════════════════════════════════════════════════════

def jaccard_neighbors_gt(
    knn_a: np.ndarray, knn_b: np.ndarray,
    ids_a: List[str], ids_b: List[str],
    shared_ids: List[str],
) -> float:
    """Per-accession Jaccard of kNN neighbours in accession space.

    For each shared accession, look up its kNN in graph A (by index),
    map those indices back to accession IDs, do the same for graph B,
    then compute Jaccard on the accession-ID sets.
    """
    idx_a = {sid: i for i, sid in enumerate(ids_a)}
    idx_b = {sid: i for i, sid in enumerate(ids_b)}

    scores = []
    for sid in shared_ids:
        ia, ib = idx_a[sid], idx_b[sid]
        # kNN indices → accession IDs
        nbrs_a = {ids_a[j] for j in knn_a[ia]}
        nbrs_b = {ids_b[j] for j in knn_b[ib]}
        # Keep only neighbours that are also shared (so we compare on same universe)
        nbrs_a_shared = nbrs_a & set(shared_ids)
        nbrs_b_shared = nbrs_b & set(shared_ids)
        if not nbrs_a_shared and not nbrs_b_shared:
            scores.append(1.0)
        elif not nbrs_a_shared or not nbrs_b_shared:
            scores.append(0.0)
        else:
            scores.append(
                len(nbrs_a_shared & nbrs_b_shared)
                / len(nbrs_a_shared | nbrs_b_shared)
            )
    return float(np.mean(scores))


def jaccard_edges_gt(
    knn_a: np.ndarray, knn_b: np.ndarray,
    ids_a: List[str], ids_b: List[str],
    shared_ids: List[str],
) -> float:
    """Jaccard of global directed edge sets in accession space.

    Edges (i, j) are converted to (accession_i, accession_j), restricted
    to the shared accession set, then compared.
    """
    shared_set = set(shared_ids)

    def _edge_set(knn, ids):
        edges = set()
        for i, row in enumerate(knn):
            src = ids[i]
            if src not in shared_set:
                continue
            for j in row:
                tgt = ids[j]
                if tgt in shared_set:
                    edges.add((src, tgt))
        return edges

    ea, eb = _edge_set(knn_a, ids_a), _edge_set(knn_b, ids_b)
    if not ea and not eb:
        return 1.0
    return len(ea & eb) / len(ea | eb)


# ═══════════════════════════════════════════════════════════════════
# Loading with accession resolution
# ═══════════════════════════════════════════════════════════════════

def _load_with_accessions(ds: Dict) -> Tuple[pd.DataFrame, List[str]]:
    """Load a panel and return (DataFrame samples×markers, accession_ids).

    For sample_columns format, the header row has numeric sample codes
    and row 3 ('Accession') has the CIP accession names.
    For marker_metrics format, sample IDs are already the column headers.
    """
    path = ROOT / ds["path"]
    X_raw, sample_ids, meta = load_genotypes(path)

    # Try to get accession from metadata
    accessions = []
    if meta:
        for sid in sample_ids:
            acc = meta.get(sid, {}).get("Accession", sid)
            accessions.append(str(acc).strip())
    else:
        accessions = [str(s).strip() for s in sample_ids]

    # Filter missingness
    X_filt, info = filter_missingness(X_raw, sample_thresh=0.50,
                                       marker_thresh=0.50)
    # Map accessions through surviving samples
    surviving_samples = list(X_filt.index)
    acc_map = dict(zip(sample_ids, accessions))
    final_accessions = [acc_map.get(s, s) for s in surviving_samples]

    X_num = X_filt.apply(pd.to_numeric, errors="coerce").values.astype(float)
    return X_num, final_accessions, surviving_samples


# ═══════════════════════════════════════════════════════════════════
# Main experiment
# ═══════════════════════════════════════════════════════════════════

def run_pair(pair_name: str, pair_cfg: Dict) -> Dict[str, Any]:
    """Run ground-truth validation for one pair of panels."""
    print(f"\n{'=' * 60}")
    print(f"  ADR-0016: Ground-truth validation — {pair_name}")
    print(f"{'=' * 60}")

    # Load both panels with accession IDs
    print(f"\n  Loading {pair_cfg['panel_a']['label']}...")
    X_a, acc_a, sids_a = _load_with_accessions(pair_cfg["panel_a"])
    print(f"    Shape: {X_a.shape}, accessions: {len(set(acc_a))}")

    print(f"  Loading {pair_cfg['panel_b']['label']}...")
    X_b, acc_b, sids_b = _load_with_accessions(pair_cfg["panel_b"])
    print(f"    Shape: {X_b.shape}, accessions: {len(set(acc_b))}")

    # Find shared accessions
    shared = sorted(set(acc_a) & set(acc_b))
    print(f"\n  Shared accessions: {len(shared)}")
    print(f"  Panel A total: {len(acc_a)}, Panel B total: {len(acc_b)}")

    if len(shared) < 10:
        print("  ⚠ Too few shared accessions, skipping pair.")
        return {"pair": pair_name, "error": "too_few_shared", "n_shared": len(shared)}

    # Handle duplicate accessions: keep first occurrence only
    def _deduplicate(X, accessions):
        seen = set()
        keep = []
        for i, acc in enumerate(accessions):
            if acc not in seen and acc in set(shared):
                seen.add(acc)
                keep.append(i)
        return X[keep], [accessions[i] for i in keep]

    X_a_shared, acc_a_shared = _deduplicate(X_a, acc_a)
    X_b_shared, acc_b_shared = _deduplicate(X_b, acc_b)

    # Sort both by accession to ensure *identical* ordering
    order_a = np.argsort(acc_a_shared)
    order_b = np.argsort(acc_b_shared)
    X_a_aligned = X_a_shared[order_a]
    X_b_aligned = X_b_shared[order_b]
    acc_aligned = [acc_a_shared[i] for i in order_a]
    acc_b_check = [acc_b_shared[i] for i in order_b]
    assert acc_aligned == acc_b_check, "Accession alignment mismatch!"

    n_aligned = len(acc_aligned)
    print(f"  Aligned samples (deduplicated): {n_aligned}")

    # Impute
    X_a_imp = _impute(X_a_aligned)
    X_b_imp = _impute(X_b_aligned)

    # Run for multiple seeds
    seed_results = []
    for seed in SEEDS:
        print(f"\n  --- seed={seed} ---")

        # Build independent graphs
        res_a = _pca_knn(X_a_imp, seed=seed)
        res_b = _pca_knn(X_b_imp, seed=seed)

        knn_a = res_a["knn_inds"]
        knn_b = res_b["knn_inds"]

        # CORRESPONDENCE-FREE metrics (index-based, as if no alignment known)
        # Since both matrices have the same row order (aligned by accession),
        # the index-based comparison is actually a correspondence-based
        # comparison — the indices refer to the same accession.
        #
        # To get a TRUE correspondence-free comparison, we need to compare
        # the graphs *without* alignment. We do this by randomly permuting
        # the rows of one graph.

        # Metric 1: Correspondence-based (aligned — ground truth)
        j_edge_gt = jaccard_edges(knn_a, knn_b)
        j_nbr_gt = jaccard_neighbors(knn_a, knn_b)
        ss_gt = subspace_similarity(res_a["X_pca"], res_b["X_pca"])

        print(f"    GT (aligned):      J_edge={j_edge_gt:.4f}  "
              f"J_nbr={j_nbr_gt:.4f}  SS={ss_gt:.4f}")

        # Metric 2: Correspondence-free (random permutation baseline)
        rng = np.random.default_rng(seed)
        perm = rng.permutation(n_aligned)
        knn_b_perm = knn_b[perm]
        # Remap indices through permutation inverse
        inv_perm = np.argsort(perm)
        knn_b_remapped = inv_perm[knn_b_perm]

        j_edge_cf = jaccard_edges(knn_a, knn_b_remapped)
        j_nbr_cf = jaccard_neighbors(knn_a, knn_b_remapped)
        ss_cf = subspace_similarity(
            res_a["X_pca"], res_b["X_pca"][perm]
        )

        print(f"    CF (shuffled):     J_edge={j_edge_cf:.4f}  "
              f"J_nbr={j_nbr_cf:.4f}  SS={ss_cf:.4f}")

        # Metric 3: Accession-space GT (full panels, not just aligned subset)
        # Use the accession-based functions on FULL panels
        res_a_full = _pca_knn(_impute(X_a), seed=seed)
        res_b_full = _pca_knn(_impute(X_b), seed=seed)

        j_edge_acc = jaccard_edges_gt(
            res_a_full["knn_inds"], res_b_full["knn_inds"],
            acc_a, acc_b, shared,
        )
        j_nbr_acc = jaccard_neighbors_gt(
            res_a_full["knn_inds"], res_b_full["knn_inds"],
            acc_a, acc_b, shared,
        )
        ss_full = subspace_similarity(
            res_a_full["X_pca"][:min(len(acc_a), len(acc_b))],
            res_b_full["X_pca"][:min(len(acc_a), len(acc_b))],
        )

        print(f"    ACC (full→shared): J_edge={j_edge_acc:.4f}  "
              f"J_nbr={j_nbr_acc:.4f}  SS={ss_full:.4f}")

        # Metric 4: Intra-panel stability (same panel, different seeds)
        # as a calibration reference
        res_a_s2 = _pca_knn(X_a_imp, seed=seed + 100)
        j_edge_intra = jaccard_edges(knn_a, res_a_s2["knn_inds"])
        j_nbr_intra = jaccard_neighbors(knn_a, res_a_s2["knn_inds"])

        print(f"    INTRA (panel A):   J_edge={j_edge_intra:.4f}  "
              f"J_nbr={j_nbr_intra:.4f}")

        seed_results.append({
            "seed": seed,
            "n_aligned": n_aligned,
            # Ground truth (aligned by accession)
            "J_edge_gt": round(j_edge_gt, 4),
            "J_nbr_gt": round(j_nbr_gt, 4),
            "SS_gt": round(ss_gt, 4),
            # Shuffled (correspondence-free null)
            "J_edge_cf_shuffled": round(j_edge_cf, 4),
            "J_nbr_cf_shuffled": round(j_nbr_cf, 4),
            "SS_cf_shuffled": round(ss_cf, 4),
            # Accession-space on full panels
            "J_edge_acc": round(j_edge_acc, 4),
            "J_nbr_acc": round(j_nbr_acc, 4),
            "SS_full": round(ss_full, 4),
            # Intra-panel calibration
            "J_edge_intra": round(j_edge_intra, 4),
            "J_nbr_intra": round(j_nbr_intra, 4),
            # PCA variance
            "cumvar_a": round(res_a["cumvar"], 4),
            "cumvar_b": round(res_b["cumvar"], 4),
        })

    return {
        "pair": pair_name,
        "panel_a": pair_cfg["panel_a"]["label"],
        "panel_b": pair_cfg["panel_b"]["label"],
        "n_a": len(acc_a),
        "n_b": len(acc_b),
        "n_shared": len(shared),
        "n_aligned": n_aligned,
        "seeds": seed_results,
    }


# ═══════════════════════════════════════════════════════════════════
# Summary table
# ═══════════════════════════════════════════════════════════════════

def _print_summary(results: Dict):
    print("\n" + "=" * 72)
    print("  ADR-0016 SUMMARY")
    print("=" * 72)

    for pair_name, pair_res in results.items():
        if "error" in pair_res:
            print(f"\n  {pair_name}: SKIPPED ({pair_res['error']})")
            continue

        print(f"\n  ── {pair_name}: {pair_res['panel_a']}  ↔  {pair_res['panel_b']}")
        print(f"     shared accessions: {pair_res['n_shared']} "
              f"  aligned: {pair_res['n_aligned']}")
        print()
        print(f"     {'metric':<20} {'GT (aligned)':>14} {'CF (shuffled)':>14} "
              f"{'ACC (full)':>14} {'INTRA (ref)':>14}")
        print(f"     {'─' * 20} {'─' * 14} {'─' * 14} {'─' * 14} {'─' * 14}")

        # Average across seeds
        seeds = pair_res["seeds"]
        for metric_key, label in [
            ("J_edge", "J_edge"),
            ("J_nbr", "J_nbr"),
            ("SS", "SS"),
        ]:
            gt = np.mean([s[f"{metric_key}_gt"] for s in seeds])
            cf = np.mean([s.get(f"{metric_key}_cf_shuffled", float("nan")) for s in seeds])
            acc = np.mean([s.get(f"{metric_key}_acc", float("nan")) for s in seeds])
            intra = np.mean([s.get(f"{metric_key}_intra", float("nan")) for s in seeds])

            gt_s = f"{gt:.4f}"
            cf_s = f"{cf:.4f}" if not np.isnan(cf) else "—"
            acc_s = f"{acc:.4f}" if not np.isnan(acc) else "—"
            intra_s = f"{intra:.4f}" if not np.isnan(intra) else "—"

            print(f"     {label:<20} {gt_s:>14} {cf_s:>14} {acc_s:>14} {intra_s:>14}")

    print()
    print("  Interpretation:")
    print("    GT (aligned):  graphs compared with known sample correspondence")
    print("    CF (shuffled): graphs compared after random row permutation (null)")
    print("    ACC (full):    full panels, edges mapped to accession space")
    print("    INTRA (ref):   same panel, PCA with different seed (upper bound)")
    print()
    print("  If GT >> CF → correspondence reveals shared structure")
    print("  If GT ≈ INTRA → cross-technology graphs are as stable as intra-panel")
    print("  If GT << INTRA → SNP & SilicoDArT capture different structure")


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════

def main():
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    results = {}
    for pair_name, pair_cfg in PAIRS.items():
        results[pair_name] = run_pair(pair_name, pair_cfg)

    _print_summary(results)

    out_path = OUT_DIR / "adr0016_ground_truth.json"
    with open(out_path, "w") as f:
        json.dump(
            {"experiment": "ADR-0016", "description": "Ground-truth validation",
             "results": results},
            f, indent=2,
        )
    print(f"\n  Results saved to {out_path}")
    print(f"  Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
