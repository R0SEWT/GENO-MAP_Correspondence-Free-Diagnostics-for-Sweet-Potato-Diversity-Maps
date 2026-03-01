#!/usr/bin/env python3
"""Level 1 — Autoencoder Genómico.

Trains a denoising autoencoder on DArT genotype data and exports the
bottleneck embeddings in the same JSON format as ``build_embeddings.py``
for downstream kNN graph construction and validation.

Architecture (per docs/roadmap/scale.md):
  Input(p) → Dense → ResBlock → ... → Bottleneck(d) → ... → Dense → Output(p)

Usage::

    # Quick dev run
    python scripts/train_autoencoder.py \\
        --input data/10.21223P30BVZYY_Genetic_diversity/SNP_Genotypes.csv \\
        --out-dir experiments/autoencoder/dev-test \\
        --epochs 10 --batch-size 256

    # Full run (4060 Ti)
    python scripts/train_autoencoder.py \\
        --input data/10.21223P30BVZYY_Genetic_diversity/SNP_Genotypes.csv \\
        --out-dir experiments/autoencoder/v2-ae \\
        --epochs 200 --patience 20 --bottleneck 64
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from load_dart import filter_missingness, load_genotypes  # noqa: E402

# PyTorch imports (fail early with helpful message)
try:
    import torch
    import torch.nn as nn
    from torch.utils.data import DataLoader, TensorDataset
except ImportError:
    print("ERROR: PyTorch not installed.  Run: pip install torch")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------


class ResBlock(nn.Module):
    """Residual MLP block with dropout."""

    def __init__(self, dim: int, dropout: float = 0.2):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(dim, dim),
            nn.BatchNorm1d(dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(dim, dim),
            nn.BatchNorm1d(dim),
        )
        self.act = nn.GELU()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.act(x + self.net(x))


class GenoAutoencoder(nn.Module):
    """Denoising autoencoder for genotype data.

    Parameters
    ----------
    input_dim : int
        Number of markers.
    bottleneck : int
        Embedding dimension.
    hidden : int
        Width of hidden layers.
    n_blocks : int
        Number of residual blocks in encoder and decoder.
    dropout : float
        Dropout rate.
    """

    def __init__(
        self,
        input_dim: int,
        bottleneck: int = 64,
        hidden: int = 512,
        n_blocks: int = 2,
        dropout: float = 0.2,
    ):
        super().__init__()
        self.input_dim = input_dim
        self.bottleneck_dim = bottleneck

        # Encoder
        enc_layers: list = [nn.Linear(input_dim, hidden), nn.GELU(), nn.Dropout(dropout)]
        for _ in range(n_blocks):
            enc_layers.append(ResBlock(hidden, dropout))
        enc_layers.append(nn.Linear(hidden, bottleneck))
        self.encoder = nn.Sequential(*enc_layers)

        # Decoder
        dec_layers: list = [nn.Linear(bottleneck, hidden), nn.GELU(), nn.Dropout(dropout)]
        for _ in range(n_blocks):
            dec_layers.append(ResBlock(hidden, dropout))
        dec_layers.append(nn.Linear(hidden, input_dim))
        self.decoder = nn.Sequential(*dec_layers)

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        return self.decoder(z)

    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        z = self.encode(x)
        x_hat = self.decode(z)
        return x_hat, z


# ---------------------------------------------------------------------------
# Training
# ---------------------------------------------------------------------------


def train_autoencoder(
    X: np.ndarray,
    *,
    bottleneck: int = 64,
    hidden: int = 512,
    n_blocks: int = 2,
    dropout: float = 0.2,
    noise_frac: float = 0.15,
    sparsity_lambda: float = 1e-4,
    lr: float = 1e-3,
    batch_size: int = 256,
    epochs: int = 200,
    patience: int = 20,
    val_frac: float = 0.1,
    seed: int = 42,
    device: str = "auto",
) -> Dict[str, Any]:
    """Train and return model + history + embeddings.

    Returns dict with keys:
      model, embeddings, history, best_epoch, device, elapsed
    """
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"[ae] Device: {device}")

    torch.manual_seed(seed)
    np.random.seed(seed)

    n_total = X.shape[0]
    n_val = max(1, int(n_total * val_frac))
    n_train = n_total - n_val
    perm = np.random.permutation(n_total)
    idx_train, idx_val = perm[:n_train], perm[n_train:]

    X_tensor = torch.tensor(X, dtype=torch.float32)
    X_train = X_tensor[idx_train]
    X_val = X_tensor[idx_val]

    train_ds = TensorDataset(X_train)
    train_dl = DataLoader(train_ds, batch_size=batch_size, shuffle=True, drop_last=False)

    model = GenoAutoencoder(
        input_dim=X.shape[1],
        bottleneck=bottleneck,
        hidden=hidden,
        n_blocks=n_blocks,
        dropout=dropout,
    ).to(device)

    optimiser = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimiser, mode="min", factor=0.5, patience=max(5, patience // 3),
    )
    mse_loss = nn.MSELoss()

    history: Dict[str, List[float]] = {"train_loss": [], "val_loss": [], "lr": []}
    best_val = float("inf")
    best_epoch = 0
    best_state = None

    t0 = time.time()
    for epoch in range(1, epochs + 1):
        # --- Train ---
        model.train()
        epoch_loss = 0.0
        for (batch_x,) in train_dl:
            batch_x = batch_x.to(device)
            # Masking noise (denoising objective)
            mask = (torch.rand_like(batch_x) < noise_frac).float()
            x_noisy = batch_x * (1 - mask)

            x_hat, z = model(x_noisy)
            recon = mse_loss(x_hat, batch_x)
            sparsity = z.abs().mean() * sparsity_lambda
            loss = recon + sparsity

            optimiser.zero_grad()
            loss.backward()
            optimiser.step()
            epoch_loss += loss.item() * batch_x.size(0)

        epoch_loss /= n_train

        # --- Val ---
        model.eval()
        with torch.no_grad():
            x_val_dev = X_val.to(device)
            x_hat_val, z_val = model(x_val_dev)
            val_loss = mse_loss(x_hat_val, x_val_dev).item()

        scheduler.step(val_loss)
        current_lr = optimiser.param_groups[0]["lr"]

        history["train_loss"].append(epoch_loss)
        history["val_loss"].append(val_loss)
        history["lr"].append(current_lr)

        if val_loss < best_val:
            best_val = val_loss
            best_epoch = epoch
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}

        if epoch % 10 == 0 or epoch == 1:
            print(f"  epoch {epoch:4d}  train={epoch_loss:.6f}  val={val_loss:.6f}  lr={current_lr:.2e}")

        # Early stopping
        if epoch - best_epoch >= patience:
            print(f"  Early stop at epoch {epoch} (best={best_epoch})")
            break

    elapsed = time.time() - t0

    # Load best weights
    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()

    # Extract embeddings for ALL samples
    with torch.no_grad():
        embeddings = model.encode(X_tensor.to(device)).cpu().numpy()

    return {
        "model": model,
        "embeddings": embeddings,
        "history": history,
        "best_epoch": best_epoch,
        "best_val_loss": best_val,
        "device": device,
        "elapsed_seconds": round(elapsed, 1),
    }


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------


def export_embeddings(
    embeddings: np.ndarray,
    sample_ids: List[str],
    sample_meta: Dict[str, Dict[str, str]],
    out_prefix: Path,
    stats_extra: Dict[str, Any] | None = None,
) -> None:
    """Write nodes/edges/stats JSON compatible with build_embeddings.py output."""
    from sklearn.decomposition import PCA as SkPCA
    from sklearn.neighbors import NearestNeighbors

    # 2D projection for visualisation
    n_comp = min(2, embeddings.shape[1])
    if embeddings.shape[1] > 2:
        pca_2d = SkPCA(n_components=2, random_state=42)
        emb_2d = pca_2d.fit_transform(embeddings)
    else:
        emb_2d = embeddings[:, :2]

    # Try UMAP on the bottleneck
    try:
        import umap  # type: ignore
        reducer = umap.UMAP(n_components=2, random_state=42)
        emb_2d = reducer.fit_transform(embeddings)
    except Exception:
        pass

    nodes = []
    for idx, sid in enumerate(sample_ids):
        nodes.append({
            "id": str(sid),
            "idx": idx,
            "embedding": emb_2d[idx].tolist(),
            "bottleneck": embeddings[idx].tolist(),
            "meta": sample_meta.get(str(sid), {}),
        })

    # kNN on bottleneck
    k = min(15, len(sample_ids) - 1)
    edges: list = []
    if k >= 1:
        nn = NearestNeighbors(n_neighbors=k + 1, metric="cosine")
        nn.fit(embeddings)
        dists, inds = nn.kneighbors(embeddings)
        for i in range(len(sample_ids)):
            for j, d in zip(inds[i][1:], dists[i][1:]):
                edges.append({
                    "source": sample_ids[i],
                    "target": sample_ids[j],
                    "distance": float(d),
                })

    stats: Dict[str, Any] = {
        "method": "autoencoder",
        "samples": len(sample_ids),
        "bottleneck_dim": embeddings.shape[1],
        "edges": len(edges),
    }
    if stats_extra:
        stats.update(stats_extra)

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix, data in [("_nodes.json", nodes), ("_edges.json", edges), ("_stats.json", stats)]:
        fpath = out_prefix.parent / f"{out_prefix.name}{suffix}"
        with open(fpath, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
        print(f"  [saved] {fpath}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description="Train genotype autoencoder (Level 1).")
    parser.add_argument("--input", required=True, type=Path, help="Genotype CSV")
    parser.add_argument("--out-dir", required=True, type=Path, help="Output directory")
    parser.add_argument("--bottleneck", type=int, default=64, help="Embedding dim")
    parser.add_argument("--hidden", type=int, default=512, help="Hidden layer width")
    parser.add_argument("--n-blocks", type=int, default=2, help="Residual blocks")
    parser.add_argument("--dropout", type=float, default=0.2)
    parser.add_argument("--noise-frac", type=float, default=0.15,
                        help="Fraction of inputs masked for denoising")
    parser.add_argument("--sparsity-lambda", type=float, default=1e-4)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--patience", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    parser.add_argument("--sample-thresh", type=float, default=0.50)
    parser.add_argument("--marker-thresh", type=float, default=0.50)
    parser.add_argument("--imputation", default="most_frequent",
                        choices=["most_frequent", "median", "mean"])
    parser.add_argument("--max-samples", type=int, default=0, help="0=no limit")
    parser.add_argument("--max-markers", type=int, default=0, help="0=no limit")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[ae] Input:      {args.input}")
    print(f"[ae] Output:     {args.out_dir}")
    print(f"[ae] Bottleneck: {args.bottleneck}  Hidden: {args.hidden}  Blocks: {args.n_blocks}")

    # --- Load data ---
    print("[ae] Loading genotypes ...")
    X_raw, sample_ids, sample_meta = load_genotypes(args.input)

    if args.max_samples and X_raw.shape[0] > args.max_samples:
        X_raw = X_raw.sample(n=args.max_samples, random_state=args.seed)
        sample_ids = list(X_raw.index.astype(str))
        sample_meta = {s: sample_meta.get(s, {}) for s in sample_ids}
    if args.max_markers and X_raw.shape[1] > args.max_markers:
        X_raw = X_raw.sample(n=args.max_markers, axis=1, random_state=args.seed)

    X_raw, filter_report = filter_missingness(
        X_raw, sample_thresh=args.sample_thresh, marker_thresh=args.marker_thresh,
    )
    sample_ids = list(X_raw.index.astype(str))
    sample_meta = {s: sample_meta.get(s, {}) for s in sample_ids}
    marker_names = list(X_raw.columns.astype(str))
    print(f"[ae] Shape after filter: {X_raw.shape}")
    print(f"[ae] Filter report: {filter_report['samples_dropped']} samples, {filter_report['markers_dropped']} markers dropped")

    # --- Impute ---
    from sklearn.impute import SimpleImputer
    imp = SimpleImputer(strategy=args.imputation)
    X_imp = imp.fit_transform(X_raw)
    print(f"[ae] Imputed ({args.imputation}), shape: {X_imp.shape}")

    # --- Train ---
    result = train_autoencoder(
        X_imp,
        bottleneck=args.bottleneck,
        hidden=args.hidden,
        n_blocks=args.n_blocks,
        dropout=args.dropout,
        noise_frac=args.noise_frac,
        sparsity_lambda=args.sparsity_lambda,
        lr=args.lr,
        batch_size=args.batch_size,
        epochs=args.epochs,
        patience=args.patience,
        seed=args.seed,
        device=args.device,
    )

    print(f"\n[ae] Training complete: best_epoch={result['best_epoch']}, "
          f"best_val={result['best_val_loss']:.6f}, elapsed={result['elapsed_seconds']}s")

    # --- Save model checkpoint ---
    ckpt_path = args.out_dir / "model_checkpoint.pt"
    torch.save({
        "model_state_dict": result["model"].state_dict(),
        "bottleneck": args.bottleneck,
        "hidden": args.hidden,
        "n_blocks": args.n_blocks,
        "dropout": args.dropout,
        "input_dim": X_imp.shape[1],
        "best_epoch": result["best_epoch"],
        "best_val_loss": result["best_val_loss"],
        "seed": args.seed,
    }, ckpt_path)
    print(f"  [saved] {ckpt_path}")

    # --- Save training history ---
    history_path = args.out_dir / "training_history.json"
    with open(history_path, "w") as f:
        json.dump(result["history"], f, indent=2)
    print(f"  [saved] {history_path}")

    # --- Training curves figure ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        epochs_range = range(1, len(result["history"]["train_loss"]) + 1)

        axes[0].plot(epochs_range, result["history"]["train_loss"], label="train")
        axes[0].plot(epochs_range, result["history"]["val_loss"], label="val")
        axes[0].axvline(result["best_epoch"], color="grey", linestyle="--", alpha=0.5, label=f"best={result['best_epoch']}")
        axes[0].set_xlabel("Epoch")
        axes[0].set_ylabel("Loss")
        axes[0].set_title("Training curves")
        axes[0].legend()
        axes[0].grid(alpha=0.3)

        axes[1].plot(epochs_range, result["history"]["lr"])
        axes[1].set_xlabel("Epoch")
        axes[1].set_ylabel("Learning rate")
        axes[1].set_title("LR schedule")
        axes[1].grid(alpha=0.3)
        axes[1].set_yscale("log")

        fig.tight_layout()
        fig_path = args.out_dir / "fig_training_curves.png"
        fig.savefig(fig_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {fig_path}")
    except Exception as e:
        print(f"  [warn] Could not generate training curves figure: {e}")

    # --- Export embeddings (same format as build_embeddings.py) ---
    out_prefix = args.out_dir / "ae_embedding"
    stats_extra = {
        "input": str(args.input),
        "bottleneck": args.bottleneck,
        "hidden": args.hidden,
        "n_blocks": args.n_blocks,
        "epochs_trained": result["best_epoch"],
        "best_val_loss": result["best_val_loss"],
        "elapsed_seconds": result["elapsed_seconds"],
        "imputation": args.imputation,
        "filter": filter_report,
        "markers": len(marker_names),
        "seed": args.seed,
    }
    export_embeddings(result["embeddings"], sample_ids, sample_meta, out_prefix, stats_extra)

    print(f"\n[ae] All done.  Embeddings: {out_prefix}_nodes.json")


if __name__ == "__main__":
    main()
