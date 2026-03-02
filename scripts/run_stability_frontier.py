#!/usr/bin/env python3
"""Level A — Representation Stability Frontier.

Answers: **When do learned embeddings (AE) surpass PCA in genomics?**

Procedure:
  1. Load the largest dataset (Global SNP, n=5970)
  2. Subsample at different n: {100, 200, 500, 1000, 2000, 3000, 5000, 5970}
  3. For each n × seed: run PCA baseline AND AE v1, measure trust + stability
  4. Plot the crossover frontier n* where AE > PCA

Usage::

    python scripts/run_stability_frontier.py --profile quick   # 3 sizes, 1 seed
    python scripts/run_stability_frontier.py --profile full    # all sizes, 3 seeds

Output: experiments/frontier/results.json + figures
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from load_dart import filter_missingness, load_genotypes  # noqa: E402

# ╭──────────────────────────────────────────────────────────────╮
# │  Profiles                                                   │
# ╰──────────────────────────────────────────────────────────────╯

PROFILES = {
    "quick": {
        "sample_sizes": [200, 1000, 5970],
        "seeds": [42],
        "ae_epochs": 50,
        "ae_patience": 15,
    },
    "local": {
        "sample_sizes": [100, 200, 500, 1000, 2000, 3000, 5000, 5970],
        "seeds": [42, 52, 62],
        "ae_epochs": 200,
        "ae_patience": 20,
    },
    "full": {
        "sample_sizes": [50, 100, 200, 500, 1000, 1500, 2000, 3000, 4000, 5000, 5970],
        "seeds": [42, 52, 62],
        "ae_epochs": 200,
        "ae_patience": 20,
    },
}

# ╭──────────────────────────────────────────────────────────────╮
# │  Data sources                                               │
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
}


# ╭──────────────────────────────────────────────────────────────╮
# │  PCA baseline computation                                   │
# ╰──────────────────────────────────────────────────────────────╯

def run_pca_condition(X_high: np.ndarray, seed: int, k: int = 15):
    """Run PCA → kNN and return trust + edges for one condition."""
    from sklearn.decomposition import PCA
    from sklearn.manifold import trustworthiness
    from sklearn.neighbors import NearestNeighbors

    n_comp = min(50, X_high.shape[1], X_high.shape[0])
    pca = PCA(n_components=n_comp, random_state=seed)
    X_pca = pca.fit_transform(X_high)

    # Trustworthiness: high-D → PCA space
    k_eff = min(k, X_high.shape[0] - 1)
    trust = trustworthiness(X_high, X_pca[:, :30], n_neighbors=k_eff)

    # kNN edges on PCA 30D
    knn_dim = min(30, X_pca.shape[1])
    nbrs = NearestNeighbors(n_neighbors=k_eff + 1, metric="cosine")
    nbrs.fit(X_pca[:, :knn_dim])
    _, inds = nbrs.kneighbors(X_pca[:, :knn_dim])
    edges = set()
    for i in range(inds.shape[0]):
        for j in inds[i, 1:]:
            edges.add((i, j))

    return trust, edges, X_pca


def run_ae_condition(
    X_high: np.ndarray,
    seed: int,
    epochs: int = 200,
    patience: int = 20,
    k: int = 15,
):
    """Train AE and return trust + edges for one condition."""
    from train_autoencoder import train_autoencoder
    from sklearn.manifold import trustworthiness
    from sklearn.neighbors import NearestNeighbors

    result = train_autoencoder(
        X_high,
        bottleneck=64,
        hidden=512,
        n_blocks=2,
        dropout=0.2,
        noise_frac=0.15,
        sparsity_lambda=1e-4,
        lr=1e-3,
        batch_size=256,
        epochs=epochs,
        patience=patience,
        seed=seed,
    )

    embeddings = result["embeddings"]  # (n, 64)

    k_eff = min(k, X_high.shape[0] - 1)
    trust = trustworthiness(X_high, embeddings, n_neighbors=k_eff)

    nbrs = NearestNeighbors(n_neighbors=k_eff + 1, metric="cosine")
    nbrs.fit(embeddings)
    _, inds = nbrs.kneighbors(embeddings)
    edges = set()
    for i in range(inds.shape[0]):
        for j in inds[i, 1:]:
            edges.add((i, j))

    return trust, edges, result


def edge_jaccard(e1: set, e2: set) -> float:
    inter = len(e1 & e2)
    union = len(e1 | e2)
    return inter / union if union > 0 else 1.0


# ╭──────────────────────────────────────────────────────────────╮
# │  Merge helper                                              │
# ╰──────────────────────────────────────────────────────────────╯

def _merge_results(path: Path, new_results: List[Dict]) -> List[Dict]:
    """Merge new results into existing file, replacing entries with same n_samples."""
    existing = []
    if path.exists():
        with open(path) as f:
            existing = json.load(f)
    by_n = {r["n_samples"]: r for r in existing}
    for r in new_results:
        by_n[r["n_samples"]] = r
    return sorted(by_n.values(), key=lambda r: r["n_samples"])


# ╭──────────────────────────────────────────────────────────────╮
# │  Main experiment loop                                       │
# ╰──────────────────────────────────────────────────────────────╯

def run_frontier(
    ds_key: str,
    profile: dict,
    out_dir: Path,
) -> List[Dict[str, Any]]:
    """Run the full frontier experiment for one dataset."""
    ds = DATASETS[ds_key]
    print(f"\n{'='*60}")
    print(f"  FRONTIER: {ds['label']} ({ds_key})")
    print(f"{'='*60}")

    # Load data once
    root = Path(__file__).resolve().parent.parent
    X_raw, sample_ids, _ = load_genotypes(root / ds["path"])
    X_raw, _ = filter_missingness(X_raw, sample_thresh=0.50, marker_thresh=0.50)

    from sklearn.impute import SimpleImputer
    imp = SimpleImputer(strategy="most_frequent")
    X_full = imp.fit_transform(X_raw)
    sample_ids = list(X_raw.index.astype(str))
    n_total = X_full.shape[0]
    n_markers = X_full.shape[1]

    print(f"  Full dataset: {n_total} samples × {n_markers} markers")

    sample_sizes = [s for s in profile["sample_sizes"] if s <= n_total]
    seeds = profile["seeds"]
    ae_epochs = profile["ae_epochs"]
    ae_patience = profile["ae_patience"]

    results = []
    seed_pairs = [(s1, s2) for i, s1 in enumerate(seeds) for s2 in seeds[i+1:]]

    for n_sub in sample_sizes:
        print(f"\n  ── n = {n_sub} ──")
        t0 = time.time()

        # Deterministic subsample (always same rows for given n)
        rng = np.random.RandomState(0)  # fixed for reproducibility
        if n_sub >= n_total:
            idx_sub = np.arange(n_total)
        else:
            idx_sub = np.sort(rng.choice(n_total, size=n_sub, replace=False))
        X_sub = X_full[idx_sub]

        pca_trusts = {}
        pca_edges = {}
        ae_trusts = {}
        ae_edges = {}
        ae_val_losses = {}
        ae_epochs_trained = {}

        for seed in seeds:
            # PCA
            trust_pca, edges_pca, _ = run_pca_condition(X_sub, seed)
            pca_trusts[seed] = trust_pca
            pca_edges[seed] = edges_pca
            print(f"    PCA seed={seed}: trust={trust_pca:.4f}")

            # AE
            trust_ae, edges_ae, ae_result = run_ae_condition(
                X_sub, seed, epochs=ae_epochs, patience=ae_patience
            )
            ae_trusts[seed] = trust_ae
            ae_edges[seed] = edges_ae
            ae_val_losses[seed] = ae_result["best_val_loss"]
            ae_epochs_trained[seed] = ae_result["best_epoch"]
            print(f"    AE  seed={seed}: trust={trust_ae:.4f}  "
                  f"val_loss={ae_result['best_val_loss']:.4f}  "
                  f"epochs={ae_result['best_epoch']}")

        # Compute stability (edge Jaccard between seed pairs)
        pca_stab = [edge_jaccard(pca_edges[s1], pca_edges[s2]) for s1, s2 in seed_pairs]
        ae_stab = [edge_jaccard(ae_edges[s1], ae_edges[s2]) for s1, s2 in seed_pairs]

        elapsed = time.time() - t0
        pca_trust_mean = np.mean(list(pca_trusts.values()))
        ae_trust_mean = np.mean(list(ae_trusts.values()))
        pca_stab_mean = np.mean(pca_stab) if pca_stab else float("nan")
        ae_stab_mean = np.mean(ae_stab) if ae_stab else float("nan")

        result = {
            "dataset": ds_key,
            "n_samples": n_sub,
            "n_markers": n_markers,
            "ratio_n_p": round(n_sub / n_markers, 4),
            "pca_trust_mean": round(pca_trust_mean, 5),
            "pca_trust_std": round(float(np.std(list(pca_trusts.values()))), 5),
            "ae_trust_mean": round(ae_trust_mean, 5),
            "ae_trust_std": round(float(np.std(list(ae_trusts.values()))), 5),
            "trust_delta": round(ae_trust_mean - pca_trust_mean, 5),
            "ae_wins_trust": bool(ae_trust_mean > pca_trust_mean),
            "pca_stab_mean": round(pca_stab_mean, 5) if not np.isnan(pca_stab_mean) else None,
            "ae_stab_mean": round(ae_stab_mean, 5) if not np.isnan(ae_stab_mean) else None,
            "stab_delta": round(ae_stab_mean - pca_stab_mean, 5) if (not np.isnan(ae_stab_mean) and not np.isnan(pca_stab_mean)) else None,
            "ae_val_loss_mean": round(float(np.mean(list(ae_val_losses.values()))), 5),
            "ae_epochs_trained": list(ae_epochs_trained.values()),
            "seeds": seeds,
            "elapsed_s": round(elapsed, 1),
        }
        results.append(result)

        win = "AE ✓" if result["ae_wins_trust"] else "PCA ✓"
        print(f"    → Trust: PCA={pca_trust_mean:.4f}  AE={ae_trust_mean:.4f}  "
              f"Δ={result['trust_delta']:+.4f}  [{win}]")
        if result["pca_stab_mean"] is not None:
            print(f"    → Stab:  PCA={pca_stab_mean:.4f}  AE={ae_stab_mean:.4f}  "
                  f"Δ={result['stab_delta']:+.4f}")
        print(f"    → Time: {elapsed:.0f}s")

        # Save incremental results (merge with any existing)
        out_dir.mkdir(parents=True, exist_ok=True)
        results_path = out_dir / f"results_{ds_key}.json"
        merged = _merge_results(results_path, results)
        with open(results_path, "w") as f:
            json.dump(merged, f, indent=2, ensure_ascii=False)

    return results


def plot_frontier(results: List[Dict], out_dir: Path, ds_label: str):
    """Generate frontier plots."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [warn] matplotlib not available, skipping plots")
        return

    df = pd.DataFrame(results)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # ── Panel 1: Trust vs n ──
    ax = axes[0]
    ax.plot(df["n_samples"], df["pca_trust_mean"], "o-", color="#4A90D9",
            label="PCA baseline", markersize=8, linewidth=2)
    ax.fill_between(df["n_samples"],
                     df["pca_trust_mean"] - df["pca_trust_std"],
                     df["pca_trust_mean"] + df["pca_trust_std"],
                     alpha=0.2, color="#4A90D9")
    ax.plot(df["n_samples"], df["ae_trust_mean"], "s-", color="#E57373",
            label="AE v1 (bn64)", markersize=8, linewidth=2)
    ax.fill_between(df["n_samples"],
                     df["ae_trust_mean"] - df["ae_trust_std"],
                     df["ae_trust_mean"] + df["ae_trust_std"],
                     alpha=0.2, color="#E57373")
    ax.set_xlabel("Number of samples (n)")
    ax.set_ylabel("Trustworthiness (k=15)")
    ax.set_title(f"Trust vs Sample Size — {ds_label}")
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_xscale("log")

    # ── Panel 2: Stability vs n ──
    ax = axes[1]
    stab_pca = df.dropna(subset=["pca_stab_mean"])
    stab_ae = df.dropna(subset=["ae_stab_mean"])
    if len(stab_pca) > 0:
        ax.plot(stab_pca["n_samples"], stab_pca["pca_stab_mean"], "o-",
                color="#4A90D9", label="PCA baseline", markersize=8, linewidth=2)
    if len(stab_ae) > 0:
        ax.plot(stab_ae["n_samples"], stab_ae["ae_stab_mean"], "s-",
                color="#E57373", label="AE v1 (bn64)", markersize=8, linewidth=2)
    ax.set_xlabel("Number of samples (n)")
    ax.set_ylabel("Edge Jaccard (seed-to-seed)")
    ax.set_title(f"Stability vs Sample Size — {ds_label}")
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_xscale("log")

    # ── Panel 3: Trust delta (AE - PCA) ──
    ax = axes[2]
    colors = ["#2E7D32" if d > 0 else "#C62828" for d in df["trust_delta"]]
    ax.bar(range(len(df)), df["trust_delta"], color=colors, alpha=0.8)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels([str(n) for n in df["n_samples"]], rotation=45, ha="right")
    ax.set_xlabel("Number of samples (n)")
    ax.set_ylabel("Δ Trust (AE − PCA)")
    ax.set_title(f"AE Advantage — {ds_label}")
    ax.axhline(0, color="black", linewidth=1)
    ax.grid(axis="y", alpha=0.3)
    for i, (delta, n) in enumerate(zip(df["trust_delta"], df["n_samples"])):
        ax.text(i, delta + (0.002 if delta >= 0 else -0.005),
                f"{delta:+.3f}", ha="center", va="bottom" if delta >= 0 else "top",
                fontsize=8, fontweight="bold")

    fig.suptitle(f"Representation Stability Frontier — {ds_label}",
                 fontsize=14, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig_path = out_dir / f"frontier_{results[0]['dataset']}.png"
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [saved] {fig_path}")


# ╭──────────────────────────────────────────────────────────────╮
# │  Main                                                       │
# ╰──────────────────────────────────────────────────────────────╯

def main():
    parser = argparse.ArgumentParser(
        description="Level A — Representation Stability Frontier"
    )
    parser.add_argument("--profile", default="local",
                        choices=list(PROFILES.keys()))
    parser.add_argument("--dataset", default="global_snp",
                        choices=list(DATASETS.keys()),
                        help="Dataset to use (default: global_snp)")
    parser.add_argument("--sizes", type=str, default=None,
                        help="Comma-separated sample sizes (overrides profile)")
    parser.add_argument("--out-dir", type=Path,
                        default=Path("experiments/frontier"))
    args = parser.parse_args()

    profile = dict(PROFILES[args.profile])  # copy
    if args.sizes:
        profile["sample_sizes"] = sorted(int(x) for x in args.sizes.split(","))
    print(f"[frontier] Profile: {args.profile}")
    print(f"[frontier] Sizes: {profile['sample_sizes']}")
    print(f"[frontier] Seeds: {profile['seeds']}")
    print(f"[frontier] AE: {profile['ae_epochs']} epochs, patience={profile['ae_patience']}")

    t0 = time.time()
    results = run_frontier(args.dataset, profile, args.out_dir)
    total = time.time() - t0

    # Load merged results for final summary & plot
    merged_path = args.out_dir / f"results_{args.dataset}.json"
    with open(merged_path) as f:
        all_results = json.load(f)

    print(f"\n{'='*60}")
    print(f"  FRONTIER COMPLETE — {total:.0f}s total")
    print(f"{'='*60}")

    # Summary
    print(f"\n{'n':>6s}  {'PCA trust':>10s}  {'AE trust':>10s}  {'Δ':>8s}  {'Winner':>8s}")
    for r in all_results:
        win = "AE" if r["ae_wins_trust"] else "PCA"
        print(f"{r['n_samples']:6d}  {r['pca_trust_mean']:10.4f}  "
              f"{r['ae_trust_mean']:10.4f}  {r['trust_delta']:+8.4f}  {win:>8s}")

    # Find crossover
    crossovers = []
    for i in range(1, len(all_results)):
        if all_results[i-1]["trust_delta"] <= 0 and all_results[i]["trust_delta"] > 0:
            crossovers.append((all_results[i-1]["n_samples"], all_results[i]["n_samples"]))
    if crossovers:
        lo, hi = crossovers[0]
        print(f"\n  ★ Crossover n* ∈ [{lo}, {hi}]")
    elif all(r["trust_delta"] > 0 for r in all_results):
        print(f"\n  ★ AE always wins (even at n={all_results[0]['n_samples']})")
    else:
        print(f"\n  ★ PCA always wins (up to n={all_results[-1]['n_samples']})")

    # Plot (using all merged data)
    plot_frontier(all_results, args.out_dir, DATASETS[args.dataset]["label"])

    print(f"\n[frontier] Results saved to {args.out_dir}")


if __name__ == "__main__":
    main()
