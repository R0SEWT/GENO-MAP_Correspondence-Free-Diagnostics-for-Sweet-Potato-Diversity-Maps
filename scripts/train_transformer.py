#!/usr/bin/env python3
"""Level 2 — Masked Genotype Transformer (MGT).

Learns contextual genotype representations via masked prediction of loci,
analogous to BERT's MLM but operating on genotype vectors.

Architecture:
  Value embedding: {0, 1, 2, MASK} → d_model
  Sinusoidal positional encoding (marker index, parameter-free)
  [CLS] token prepended
  TransformerEncoder(n_layers, n_heads)
  Prediction head: d_model → n_classes (per masked position)
  Projection: d_model → bottleneck_dim (CLS representation)

Training:
  - Random marker sub-sampling (context_len per iteration)
  - 15% random masking → CrossEntropy on masked positions
  - Marker sub-sampling acts as data augmentation

Inference:
  - Monte Carlo CLS: average [CLS] over K random subsets
  - UMAP 2D on CLS embeddings
  - kNN graph on CLS embeddings

Usage::

    python scripts/train_transformer.py \\
        --input data/10.21223P30BVZYY_Genetic_diversity/SNP_Genotypes.csv \\
        --out-dir experiments/global_snp/mgt-v1/seed42 \\
        --epochs 50 --context-len 512

    # With fewer params for small datasets
    python scripts/train_transformer.py \\
        --input data/.../SNPs.csv --out-dir experiments/.../seed42 \\
        --d-model 64 --n-heads 4 --n-layers 2 --d-ff 128 --context-len 256
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from load_dart import filter_missingness, load_genotypes  # noqa: E402

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    from torch.utils.data import DataLoader, Dataset
except ImportError:
    print("ERROR: PyTorch not installed.  Run: pip install torch")
    sys.exit(1)


# ╭──────────────────────────────────────────────────────────────╮
# │  Model                                                      │
# ╰──────────────────────────────────────────────────────────────╯

# Token vocabulary: 0=geno_0, 1=geno_1, 2=geno_2, 3=MASK, 4=CLS
MASK_TOKEN = 3
CLS_TOKEN = 4
VOCAB_SIZE = 5  # 0,1,2,MASK,CLS


class SinusoidalPositionalEncoding(nn.Module):
    """Fixed sinusoidal encoding for marker positions (parameter-free).

    Supports arbitrary marker indices (up to ``max_pos``), so random
    sub-sampling of markers preserves positional identity.
    """

    def __init__(self, d_model: int, max_pos: int = 100_000):
        super().__init__()
        pe = torch.zeros(max_pos, d_model)
        pos = torch.arange(0, max_pos, dtype=torch.float32).unsqueeze(1)
        div = torch.exp(
            torch.arange(0, d_model, 2, dtype=torch.float32)
            * (-math.log(10000.0) / d_model)
        )
        pe[:, 0::2] = torch.sin(pos * div)
        pe[:, 1::2] = torch.cos(pos * div)
        self.register_buffer("pe", pe)  # (max_pos, d_model)

    def forward(self, indices: torch.Tensor) -> torch.Tensor:
        """indices: (B, L) of marker positions → (B, L, d_model)"""
        return self.pe[indices]


class GenoTransformer(nn.Module):
    """Masked Genotype Transformer.

    Parameters
    ----------
    max_markers : int
        Maximum marker index (for positional encoding buffer).
    d_model : int
        Transformer hidden dimension.
    n_heads : int
        Number of attention heads.
    n_layers : int
        Number of TransformerEncoder layers.
    d_ff : int
        Feed-forward inner dimension.
    dropout : float
        Dropout rate.
    n_values : int
        Number of genotype classes (3 for SNP: 0/1/2, 2 for SilicoDArT).
    bottleneck_dim : int
        Dimension of the projected CLS representation.
    """

    def __init__(
        self,
        max_markers: int = 100_000,
        d_model: int = 128,
        n_heads: int = 4,
        n_layers: int = 4,
        d_ff: int = 256,
        dropout: float = 0.1,
        n_values: int = 3,
        bottleneck_dim: int = 64,
    ):
        super().__init__()
        self.d_model = d_model
        self.bottleneck_dim = bottleneck_dim
        self.n_values = n_values

        # Token embedding (0, 1, 2, MASK, CLS)
        self.value_embedding = nn.Embedding(VOCAB_SIZE, d_model)

        # Positional encoding (sinusoidal, parameter-free)
        self.pos_encoding = SinusoidalPositionalEncoding(d_model, max_markers + 1)

        # Transformer encoder
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=n_heads,
            dim_feedforward=d_ff,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,  # Pre-LN for better training stability
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=n_layers)

        # Prediction head: marker token → genotype class
        self.pred_head = nn.Sequential(
            nn.Linear(d_model, d_model),
            nn.GELU(),
            nn.LayerNorm(d_model),
            nn.Linear(d_model, n_values),
        )

        # CLS projection to bottleneck
        self.bottleneck_proj = nn.Sequential(
            nn.Linear(d_model, bottleneck_dim),
            nn.LayerNorm(bottleneck_dim),
        )

        self._init_weights()

    def _init_weights(self):
        """Xavier uniform initialization."""
        for p in self.parameters():
            if p.dim() > 1:
                nn.init.xavier_uniform_(p)

    def forward(
        self,
        tokens: torch.Tensor,
        marker_indices: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Parameters
        ----------
        tokens : (B, L)
            Genotype values or MASK token (CLS prepended by caller).
        marker_indices : (B, L)
            Absolute marker positions (CLS position = max_markers).

        Returns
        -------
        cls_emb : (B, bottleneck_dim)
            Projected CLS embedding.
        logits : (B, L-1, n_values)
            Genotype predictions for marker tokens (excluding CLS at pos 0).
        """
        x = self.value_embedding(tokens) + self.pos_encoding(marker_indices)
        x = self.encoder(x)

        cls_out = x[:, 0, :]          # CLS token
        marker_out = x[:, 1:, :]      # marker tokens

        cls_emb = self.bottleneck_proj(cls_out)
        logits = self.pred_head(marker_out)

        return cls_emb, logits

    def count_parameters(self) -> int:
        return sum(p.numel() for p in self.parameters() if p.requires_grad)


# ╭──────────────────────────────────────────────────────────────╮
# │  Dataset                                                    │
# ╰──────────────────────────────────────────────────────────────╯


class GenotypeMLMDataset(Dataset):
    """Genotype dataset with random marker sub-sampling + masking.

    Each ``__getitem__`` call:
      1. Randomly selects ``context_len`` markers.
      2. Masks ~``mask_ratio`` of them.
      3. Prepends a [CLS] token.
      4. Returns (tokens, marker_indices, targets, mask_flags).
    """

    def __init__(
        self,
        X: np.ndarray,
        context_len: int = 512,
        mask_ratio: float = 0.15,
    ):
        self.X = torch.tensor(X, dtype=torch.long)
        self.n_samples, self.n_markers = self.X.shape
        self.context_len = min(context_len, self.n_markers)
        self.mask_ratio = mask_ratio
        # CLS gets a fixed position index = n_markers (beyond real markers)
        self.cls_pos = self.n_markers

    def __len__(self) -> int:
        return self.n_samples

    def __getitem__(self, idx: int):
        genotype = self.X[idx]

        # 1) Random marker sub-sample (sorted for stable positional order)
        perm = torch.randperm(self.n_markers)[: self.context_len]
        perm, _ = perm.sort()
        values = genotype[perm].clone()           # (L,)
        targets = genotype[perm].clone()           # (L,)

        # 2) Random masking
        n_mask = max(1, int(self.context_len * self.mask_ratio))
        mask_perm = torch.randperm(self.context_len)[:n_mask]
        mask_flags = torch.zeros(self.context_len, dtype=torch.bool)
        mask_flags[mask_perm] = True
        values[mask_perm] = MASK_TOKEN

        # 3) Prepend CLS
        cls_tok = torch.tensor([CLS_TOKEN], dtype=torch.long)
        tokens = torch.cat([cls_tok, values])                   # (L+1,)
        cls_idx = torch.tensor([self.cls_pos], dtype=torch.long)
        marker_indices = torch.cat([cls_idx, perm])             # (L+1,)

        return tokens, marker_indices, targets, mask_flags


class GenotypeInferenceDataset(Dataset):
    """Fixed subset of markers for deterministic inference (no masking)."""

    def __init__(self, X: np.ndarray, marker_subset: np.ndarray):
        self.X = torch.tensor(X, dtype=torch.long)
        self.marker_subset = torch.tensor(marker_subset, dtype=torch.long)
        self.n_samples = self.X.shape[0]
        self.cls_pos = self.X.shape[1]  # beyond real markers

    def __len__(self) -> int:
        return self.n_samples

    def __getitem__(self, idx: int):
        values = self.X[idx, self.marker_subset]
        cls_tok = torch.tensor([CLS_TOKEN], dtype=torch.long)
        tokens = torch.cat([cls_tok, values])
        cls_idx = torch.tensor([self.cls_pos], dtype=torch.long)
        marker_indices = torch.cat([cls_idx, self.marker_subset])
        return tokens, marker_indices


# ╭──────────────────────────────────────────────────────────────╮
# │  Training                                                   │
# ╰──────────────────────────────────────────────────────────────╯


def train_transformer(
    X: np.ndarray,
    *,
    context_len: int = 512,
    mask_ratio: float = 0.15,
    d_model: int = 128,
    n_heads: int = 4,
    n_layers: int = 4,
    d_ff: int = 256,
    dropout: float = 0.1,
    bottleneck_dim: int = 64,
    n_values: int = 3,
    lr: float = 5e-4,
    weight_decay: float = 1e-4,
    batch_size: int = 64,
    epochs: int = 50,
    patience: int = 15,
    val_frac: float = 0.1,
    n_mc_subsets: int = 20,
    seed: int = 42,
    device: str = "auto",
) -> Dict[str, Any]:
    """Train Masked Genotype Transformer and extract CLS embeddings.

    Returns dict with keys:
      model, embeddings, history, best_epoch, best_val_loss,
      best_val_acc, device, elapsed_seconds
    """
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"[mgt] Device: {device}")

    torch.manual_seed(seed)
    np.random.seed(seed)

    n_total, n_markers = X.shape

    # ── Val split ──
    n_val = max(1, int(n_total * val_frac))
    n_train = n_total - n_val
    perm = np.random.permutation(n_total)
    idx_train, idx_val = perm[:n_train], perm[n_train:]

    train_ds = GenotypeMLMDataset(X[idx_train], context_len, mask_ratio)
    val_ds = GenotypeMLMDataset(X[idx_val], context_len, mask_ratio)
    train_dl = DataLoader(train_ds, batch_size=batch_size, shuffle=True, drop_last=False)
    val_dl = DataLoader(val_ds, batch_size=batch_size, shuffle=False)

    # ── Model ──
    model = GenoTransformer(
        max_markers=n_markers,
        d_model=d_model,
        n_heads=n_heads,
        n_layers=n_layers,
        d_ff=d_ff,
        dropout=dropout,
        n_values=n_values,
        bottleneck_dim=bottleneck_dim,
    ).to(device)

    n_params = model.count_parameters()
    print(f"[mgt] Parameters: {n_params:,}  ({n_params/1e6:.2f}M)")
    print(f"[mgt] Context: {context_len} markers  Mask: {mask_ratio:.0%}")
    print(f"[mgt] Arch: d={d_model} heads={n_heads} layers={n_layers} ff={d_ff} bn={bottleneck_dim}")

    # ── Optimiser + scheduler ──
    optimiser = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    # Linear warmup for 5 epochs, then cosine decay
    warmup_epochs = min(5, epochs // 4)

    def lr_lambda(epoch: int) -> float:
        if epoch < warmup_epochs:
            return (epoch + 1) / warmup_epochs
        progress = (epoch - warmup_epochs) / max(1, epochs - warmup_epochs)
        return 0.5 * (1 + math.cos(math.pi * progress))

    scheduler = torch.optim.lr_scheduler.LambdaLR(optimiser, lr_lambda)

    ce_loss = nn.CrossEntropyLoss()

    history: Dict[str, List[float]] = {
        "train_loss": [], "val_loss": [],
        "train_acc": [], "val_acc": [],
        "lr": [],
    }
    best_val = float("inf")
    best_val_acc = 0.0
    best_epoch = 0
    best_state = None

    t0 = time.time()
    for epoch in range(1, epochs + 1):
        # ── Train ──
        model.train()
        epoch_loss, epoch_correct, epoch_masked = 0.0, 0, 0
        for tokens, marker_indices, targets, mask_flags in train_dl:
            tokens = tokens.to(device)
            marker_indices = marker_indices.to(device)
            targets = targets.to(device)
            mask_flags = mask_flags.to(device)

            cls_emb, logits = model(tokens, marker_indices)
            # logits: (B, L, n_values), targets: (B, L)
            # Only compute loss on masked positions
            masked_logits = logits[mask_flags]       # (N_masked, n_values)
            masked_targets = targets[mask_flags]     # (N_masked,)
            loss = ce_loss(masked_logits, masked_targets)

            optimiser.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimiser.step()

            epoch_loss += loss.item() * masked_targets.size(0)
            epoch_correct += (masked_logits.argmax(-1) == masked_targets).sum().item()
            epoch_masked += masked_targets.size(0)

        train_loss = epoch_loss / max(1, epoch_masked)
        train_acc = epoch_correct / max(1, epoch_masked)

        # ── Val ──
        model.eval()
        val_loss_sum, val_correct, val_masked = 0.0, 0, 0
        with torch.no_grad():
            for tokens, marker_indices, targets, mask_flags in val_dl:
                tokens = tokens.to(device)
                marker_indices = marker_indices.to(device)
                targets = targets.to(device)
                mask_flags = mask_flags.to(device)

                _, logits = model(tokens, marker_indices)
                masked_logits = logits[mask_flags]
                masked_targets = targets[mask_flags]
                loss = ce_loss(masked_logits, masked_targets)

                val_loss_sum += loss.item() * masked_targets.size(0)
                val_correct += (masked_logits.argmax(-1) == masked_targets).sum().item()
                val_masked += masked_targets.size(0)

        val_loss = val_loss_sum / max(1, val_masked)
        val_acc = val_correct / max(1, val_masked)

        scheduler.step()
        current_lr = optimiser.param_groups[0]["lr"]

        history["train_loss"].append(train_loss)
        history["val_loss"].append(val_loss)
        history["train_acc"].append(train_acc)
        history["val_acc"].append(val_acc)
        history["lr"].append(current_lr)

        if val_loss < best_val:
            best_val = val_loss
            best_val_acc = val_acc
            best_epoch = epoch
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}

        if epoch % 5 == 0 or epoch == 1:
            print(
                f"  epoch {epoch:3d}  "
                f"loss={train_loss:.4f}/{val_loss:.4f}  "
                f"acc={train_acc:.3f}/{val_acc:.3f}  "
                f"lr={current_lr:.2e}"
            )

        # Early stopping
        if epoch - best_epoch >= patience:
            print(f"  Early stop at epoch {epoch} (best={best_epoch})")
            break

    elapsed = time.time() - t0

    # ── Load best weights ──
    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()

    # ── Extract CLS embeddings via Monte Carlo sub-sampling ──
    print(f"[mgt] Extracting CLS embeddings (MC={n_mc_subsets} subsets) ...")
    X_all = torch.tensor(X, dtype=torch.long)
    all_cls = []

    for mc in range(n_mc_subsets):
        mc_seed = seed + mc * 1000
        rng = np.random.RandomState(mc_seed)
        subset = np.sort(rng.choice(n_markers, size=min(context_len, n_markers), replace=False))

        inf_ds = GenotypeInferenceDataset(X, subset)
        inf_dl = DataLoader(inf_ds, batch_size=batch_size, shuffle=False)

        cls_parts = []
        with torch.no_grad():
            for tokens, marker_indices in inf_dl:
                tokens = tokens.to(device)
                marker_indices = marker_indices.to(device)
                cls_emb, _ = model(tokens, marker_indices)
                cls_parts.append(cls_emb.cpu())

        all_cls.append(torch.cat(cls_parts, dim=0))  # (N, bn_dim)

    # Average CLS across MC subsets → stable representation
    cls_stack = torch.stack(all_cls, dim=0)  # (MC, N, bn_dim)
    embeddings = cls_stack.mean(dim=0).numpy()  # (N, bn_dim)
    cls_std = cls_stack.std(dim=0).mean().item()  # cross-MC variation

    print(f"[mgt] CLS embeddings shape: {embeddings.shape}  cross-MC std: {cls_std:.4f}")

    return {
        "model": model,
        "embeddings": embeddings,
        "history": history,
        "best_epoch": best_epoch,
        "best_val_loss": best_val,
        "best_val_acc": best_val_acc,
        "cls_mc_std": cls_std,
        "device": device,
        "elapsed_seconds": round(elapsed, 1),
        "n_params": n_params,
    }


# ╭──────────────────────────────────────────────────────────────╮
# │  Export (same JSON format as train_autoencoder.py)           │
# ╰──────────────────────────────────────────────────────────────╯


def export_embeddings(
    embeddings: np.ndarray,
    sample_ids: List[str],
    sample_meta: Dict[str, Dict[str, str]],
    out_prefix: Path,
    stats_extra: Dict[str, Any] | None = None,
) -> None:
    """Write nodes/edges/stats JSON compatible with AE / baseline output."""
    from sklearn.neighbors import NearestNeighbors

    # 2D via UMAP on CLS embeddings
    try:
        import umap

        reducer = umap.UMAP(n_components=2, random_state=42)
        emb_2d = reducer.fit_transform(embeddings)
    except Exception:
        from sklearn.decomposition import PCA as SkPCA

        emb_2d = SkPCA(n_components=2, random_state=42).fit_transform(embeddings)

    nodes = []
    for idx, sid in enumerate(sample_ids):
        nodes.append(
            {
                "id": str(sid),
                "idx": idx,
                "embedding": emb_2d[idx].tolist(),
                "bottleneck": embeddings[idx].tolist(),
                "meta": sample_meta.get(str(sid), {}),
            }
        )

    # kNN on CLS embeddings
    k = min(15, len(sample_ids) - 1)
    edges: list = []
    if k >= 1:
        nbrs = NearestNeighbors(n_neighbors=k + 1, metric="cosine")
        nbrs.fit(embeddings)
        dists, inds = nbrs.kneighbors(embeddings)
        for i in range(len(sample_ids)):
            for j, d in zip(inds[i][1:], dists[i][1:]):
                edges.append(
                    {
                        "source": sample_ids[i],
                        "target": sample_ids[j],
                        "distance": float(d),
                    }
                )

    stats: Dict[str, Any] = {
        "method": "masked_genotype_transformer",
        "samples": len(sample_ids),
        "bottleneck_dim": embeddings.shape[1],
        "edges": len(edges),
    }
    if stats_extra:
        stats.update(stats_extra)

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix, data in [
        ("_nodes.json", nodes),
        ("_edges.json", edges),
        ("_stats.json", stats),
    ]:
        fpath = out_prefix.parent / f"{out_prefix.name}{suffix}"
        with open(fpath, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
        print(f"  [saved] {fpath}")


# ╭──────────────────────────────────────────────────────────────╮
# │  Main                                                       │
# ╰──────────────────────────────────────────────────────────────╯


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Train Masked Genotype Transformer (Level 2)."
    )
    parser.add_argument("--input", required=True, type=Path, help="Genotype CSV")
    parser.add_argument("--out-dir", required=True, type=Path, help="Output directory")

    # Architecture
    parser.add_argument("--d-model", type=int, default=128)
    parser.add_argument("--n-heads", type=int, default=4)
    parser.add_argument("--n-layers", type=int, default=4)
    parser.add_argument("--d-ff", type=int, default=256)
    parser.add_argument("--bottleneck", type=int, default=64)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--n-values", type=int, default=3,
                        help="Genotype classes (3 for SNP, 2 for SilicoDArT)")

    # Training
    parser.add_argument("--context-len", type=int, default=512,
                        help="Markers sub-sampled per iteration")
    parser.add_argument("--mask-ratio", type=float, default=0.15)
    parser.add_argument("--lr", type=float, default=5e-4)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--patience", type=int, default=15)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])

    # Inference
    parser.add_argument("--n-mc-subsets", type=int, default=20,
                        help="Monte Carlo subsets for CLS extraction")

    # Data
    parser.add_argument("--sample-thresh", type=float, default=0.50)
    parser.add_argument("--marker-thresh", type=float, default=0.50)
    parser.add_argument("--imputation", default="most_frequent",
                        choices=["most_frequent", "median", "mean"])
    parser.add_argument("--max-samples", type=int, default=0, help="0=no limit")
    parser.add_argument("--max-markers", type=int, default=0, help="0=no limit")

    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[mgt] Input:      {args.input}")
    print(f"[mgt] Output:     {args.out_dir}")
    print(f"[mgt] Arch:       d={args.d_model} h={args.n_heads} L={args.n_layers} ff={args.d_ff}")

    # ── Load data ──
    print("[mgt] Loading genotypes ...")
    X_raw, sample_ids, sample_meta = load_genotypes(args.input)

    if args.max_samples and X_raw.shape[0] > args.max_samples:
        X_raw = X_raw.sample(n=args.max_samples, random_state=args.seed)
        sample_ids = list(X_raw.index.astype(str))
        sample_meta = {s: sample_meta.get(s, {}) for s in sample_ids}
    if args.max_markers and X_raw.shape[1] > args.max_markers:
        X_raw = X_raw.sample(n=args.max_markers, axis=1, random_state=args.seed)

    X_raw, filter_report = filter_missingness(
        X_raw,
        sample_thresh=args.sample_thresh,
        marker_thresh=args.marker_thresh,
    )
    sample_ids = list(X_raw.index.astype(str))
    sample_meta = {s: sample_meta.get(s, {}) for s in sample_ids}
    marker_names = list(X_raw.columns.astype(str))
    print(f"[mgt] Shape after filter: {X_raw.shape}")

    # ── Impute ──
    from sklearn.impute import SimpleImputer

    imp = SimpleImputer(strategy=args.imputation)
    X_imp = imp.fit_transform(X_raw)
    # Round to nearest int and clip to valid range
    X_int = np.clip(np.round(X_imp).astype(int), 0, args.n_values - 1)
    print(f"[mgt] Imputed ({args.imputation}), tokenised to {{0..{args.n_values-1}}}")
    unique_vals = np.unique(X_int)
    print(f"[mgt] Unique token values: {unique_vals.tolist()}")

    # ── Train ──
    result = train_transformer(
        X_int,
        context_len=args.context_len,
        mask_ratio=args.mask_ratio,
        d_model=args.d_model,
        n_heads=args.n_heads,
        n_layers=args.n_layers,
        d_ff=args.d_ff,
        dropout=args.dropout,
        bottleneck_dim=args.bottleneck,
        n_values=args.n_values,
        lr=args.lr,
        weight_decay=args.weight_decay,
        batch_size=args.batch_size,
        epochs=args.epochs,
        patience=args.patience,
        n_mc_subsets=args.n_mc_subsets,
        seed=args.seed,
        device=args.device,
    )

    print(
        f"\n[mgt] Training complete: best_epoch={result['best_epoch']}, "
        f"val_loss={result['best_val_loss']:.4f}, "
        f"val_acc={result['best_val_acc']:.3f}, "
        f"elapsed={result['elapsed_seconds']}s"
    )

    # ── Save checkpoint ──
    ckpt_path = args.out_dir / "model_checkpoint.pt"
    torch.save(
        {
            "model_state_dict": result["model"].state_dict(),
            "d_model": args.d_model,
            "n_heads": args.n_heads,
            "n_layers": args.n_layers,
            "d_ff": args.d_ff,
            "bottleneck": args.bottleneck,
            "dropout": args.dropout,
            "n_values": args.n_values,
            "max_markers": X_int.shape[1],
            "best_epoch": result["best_epoch"],
            "best_val_loss": result["best_val_loss"],
            "best_val_acc": result["best_val_acc"],
            "seed": args.seed,
        },
        ckpt_path,
    )
    print(f"  [saved] {ckpt_path}")

    # ── Save training history ──
    history_path = args.out_dir / "training_history.json"
    with open(history_path, "w") as f:
        json.dump(result["history"], f, indent=2)
    print(f"  [saved] {history_path}")

    # ── Training curves figure ──
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        ep_range = range(1, len(result["history"]["train_loss"]) + 1)

        # Loss
        axes[0].plot(ep_range, result["history"]["train_loss"], label="train")
        axes[0].plot(ep_range, result["history"]["val_loss"], label="val")
        axes[0].axvline(
            result["best_epoch"], color="grey", ls="--", alpha=0.5,
            label=f"best={result['best_epoch']}",
        )
        axes[0].set_xlabel("Epoch")
        axes[0].set_ylabel("CE Loss")
        axes[0].set_title("Training Loss")
        axes[0].legend()
        axes[0].grid(alpha=0.3)

        # Accuracy
        axes[1].plot(ep_range, result["history"]["train_acc"], label="train")
        axes[1].plot(ep_range, result["history"]["val_acc"], label="val")
        axes[1].set_xlabel("Epoch")
        axes[1].set_ylabel("Accuracy")
        axes[1].set_title("MLM Accuracy")
        axes[1].legend()
        axes[1].grid(alpha=0.3)

        # LR
        axes[2].plot(ep_range, result["history"]["lr"])
        axes[2].set_xlabel("Epoch")
        axes[2].set_ylabel("Learning Rate")
        axes[2].set_title("LR Schedule")
        axes[2].set_yscale("log")
        axes[2].grid(alpha=0.3)

        fig.tight_layout()
        fig_path = args.out_dir / "fig_training_curves.png"
        fig.savefig(fig_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {fig_path}")
    except Exception as e:
        print(f"  [warn] Could not generate training curves figure: {e}")

    # ── Export embeddings ──
    out_prefix = args.out_dir / "ae_embedding"  # same name for notebook compat
    stats_extra = {
        "input": str(args.input),
        "arch": f"d{args.d_model}_h{args.n_heads}_L{args.n_layers}_ff{args.d_ff}",
        "bottleneck": args.bottleneck,
        "context_len": args.context_len,
        "mask_ratio": args.mask_ratio,
        "n_values": args.n_values,
        "n_params": result["n_params"],
        "epochs_trained": result["best_epoch"],
        "best_val_loss": result["best_val_loss"],
        "best_val_acc": result["best_val_acc"],
        "cls_mc_std": result["cls_mc_std"],
        "n_mc_subsets": args.n_mc_subsets,
        "elapsed_seconds": result["elapsed_seconds"],
        "imputation": args.imputation,
        "filter": filter_report,
        "markers": len(marker_names),
        "seed": args.seed,
    }
    export_embeddings(
        result["embeddings"], sample_ids, sample_meta, out_prefix, stats_extra
    )

    print(f"\n[mgt] All done.  Embeddings: {out_prefix}_nodes.json")


if __name__ == "__main__":
    main()
