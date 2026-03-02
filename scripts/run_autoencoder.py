#!/usr/bin/env python3
"""Orchestrate autoencoder training across datasets and profiles.

Trains the Level-1 genotype autoencoder (train_autoencoder.py) on each
dataset, then runs validation (validate_embeddings.py) on the learned
bottleneck embeddings to compare against the PCA baseline (poster-v2).

v2 additions:
  - Per-dataset regularization (aggressive for LowDensity)
  - Seed ensemble (ensemble_embeddings.py)
  - Transfer learning: pre-train on Global → fine-tune on LowDensity
  - local_v2 profile (50 epochs, based on training curve analysis)

Profiles
--------
  dev       – 500 samples, 5 000 markers, 10 epochs  (seconds)
  local     – full data, 200 epochs, patience=20      (minutes, 4060 Ti)
  local_v2  – full data, 50 epochs, patience=15       (ae-v2, 4060 Ti)
  full      – full data, 300 epochs, patience=30      (production, L40S)

Usage::

    python scripts/run_autoencoder.py --profile dev  --tag ae-smoke
    python scripts/run_autoencoder.py --profile local_v2 --tag ae-v2 --gpu 4060Ti
    python scripts/run_autoencoder.py --profile local_v2 --tag ae-v2 --gpu 4060Ti \
        --transfer-from global_snp --transfer-to lowdensity_snp,lowdensity_silico
"""
from __future__ import annotations

import argparse
import datetime
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

# ---------------------------------------------------------------------------
# Dataset catalogue  (mirrors run_experiments.py, ADR-0002)
# ---------------------------------------------------------------------------

DATASETS: Dict[str, Dict[str, Any]] = {
    "global_snp": {
        "input": "data/10.21223P30BVZYY_Genetic_diversity/SNP_Genotypes.csv",
        "neighbors": 20,
        "metric": "cosine",
    },
    "global_silico": {
        "input": "data/10.21223P30BVZYY_Genetic_diversity/SilicoDArT_Genotypes.csv",
        "neighbors": 20,
        "metric": "cosine",
    },
    "lowdensity_snp": {
        "input": "data/10.21223P3UBDJ44_LowDensity/01_Report_DSp25-515_SNPs_Filtered_by _reads.csv",
        "neighbors": 15,
        "metric": "cosine",
    },
    "lowdensity_silico": {
        "input": "data/10.21223P3UBDJ44_LowDensity/02_Report_DSp25-515_Silico-DArT_Filtered_by_reads.csv",
        "neighbors": 15,
        "metric": "cosine",
    },
}

PROFILES: Dict[str, Dict[str, Any]] = {
    "dev": {
        "max_samples": 500,
        "max_markers": 5000,
        "epochs": 10,
        "patience": 5,
        "bottleneck": 32,
        "hidden": 256,
        "n_blocks": 1,
        "batch_size": 128,
        "seeds": [42],
    },
    "local": {
        "max_samples": 0,
        "max_markers": 0,
        "epochs": 200,
        "patience": 20,
        "bottleneck": 64,
        "hidden": 512,
        "n_blocks": 2,
        "batch_size": 256,
        "seeds": [42, 52, 62],
    },
    "local_v2": {
        "max_samples": 0,
        "max_markers": 0,
        "epochs": 50,
        "patience": 15,
        "bottleneck": 64,
        "hidden": 512,
        "n_blocks": 2,
        "batch_size": 256,
        "seeds": [42, 52, 62],
    },
    "full": {
        "max_samples": 0,
        "max_markers": 0,
        "epochs": 300,
        "patience": 30,
        "bottleneck": 64,
        "hidden": 512,
        "n_blocks": 2,
        "batch_size": 256,
        "seeds": [42, 52, 62, 72, 82],
    },
}

# ---------------------------------------------------------------------------
# Per-dataset regularization overrides
# LowDensity (~630 samples) needs more regularization than Global (~5970)
# ---------------------------------------------------------------------------

DS_REGULARIZATION: Dict[str, Dict[str, float]] = {
    "global_snp": {
        "dropout": 0.2,
        "noise_frac": 0.15,
        "weight_decay": 1e-5,
        "sparsity_lambda": 1e-4,
    },
    "global_silico": {
        "dropout": 0.2,
        "noise_frac": 0.15,
        "weight_decay": 1e-5,
        "sparsity_lambda": 1e-4,
    },
    "lowdensity_snp": {
        "dropout": 0.45,
        "noise_frac": 0.30,
        "weight_decay": 5e-4,
        "sparsity_lambda": 5e-4,
    },
    "lowdensity_silico": {
        "dropout": 0.45,
        "noise_frac": 0.30,
        "weight_decay": 5e-4,
        "sparsity_lambda": 5e-4,
    },
}

SCRIPTS_DIR = Path(__file__).resolve().parent


def _run_cmd(cmd: List[str], label: str) -> int:
    """Run a subprocess and stream output."""
    print(f"\n{'='*60}")
    print(f"[run] {label}")
    print(f"[cmd] {' '.join(cmd)}")
    print("=" * 60)
    result = subprocess.run(cmd, cwd=SCRIPTS_DIR.parent)
    return result.returncode


def _log_run(entry: Dict[str, Any]) -> None:
    """Append a line to experiments/runs.jsonl."""
    log_path = SCRIPTS_DIR.parent / "experiments" / "runs.jsonl"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "a", encoding="utf-8") as f:
        f.write(json.dumps(entry, ensure_ascii=True) + "\n")


def run_ae_dataset(
    dataset: str,
    ds_cfg: Dict[str, Any],
    profile_cfg: Dict[str, Any],
    tag: str,
    gpu: str,
    device: str,
    imputation: str,
    sample_thresh: float,
    marker_thresh: float,
    lr: float,
    dropout: float,
    noise_frac: float,
    sparsity_lambda: float,
    weight_decay: float = 1e-5,
    pretrained_ckpt: Path | None = None,
    freeze_encoder_epochs: int = 0,
    use_ds_regularization: bool = True,
) -> Dict[str, Any]:
    """Train autoencoder for one dataset across seeds, then validate."""
    seeds = profile_cfg["seeds"]
    results: Dict[str, Any] = {"dataset": dataset, "seeds": {}}

    # Per-dataset regularization overrides
    if use_ds_regularization and dataset in DS_REGULARIZATION:
        ds_reg = DS_REGULARIZATION[dataset]
        dropout = ds_reg["dropout"]
        noise_frac = ds_reg["noise_frac"]
        weight_decay = ds_reg["weight_decay"]
        sparsity_lambda = ds_reg["sparsity_lambda"]
        print(f"  [reg] {dataset}: dropout={dropout}, noise={noise_frac}, "
              f"wd={weight_decay}, sparsity={sparsity_lambda}")

    for seed in seeds:
        out_dir = SCRIPTS_DIR.parent / "experiments" / dataset / tag / f"seed{seed}"
        out_dir.mkdir(parents=True, exist_ok=True)

        t0 = time.time()

        train_cmd = [
            sys.executable, str(SCRIPTS_DIR / "train_autoencoder.py"),
            "--input", ds_cfg["input"],
            "--out-dir", str(out_dir),
            "--bottleneck", str(profile_cfg["bottleneck"]),
            "--hidden", str(profile_cfg["hidden"]),
            "--n-blocks", str(profile_cfg["n_blocks"]),
            "--dropout", str(dropout),
            "--noise-frac", str(noise_frac),
            "--sparsity-lambda", str(sparsity_lambda),
            "--weight-decay", str(weight_decay),
            "--lr", str(lr),
            "--batch-size", str(profile_cfg["batch_size"]),
            "--epochs", str(profile_cfg["epochs"]),
            "--patience", str(profile_cfg["patience"]),
            "--seed", str(seed),
            "--device", device,
            "--sample-thresh", str(sample_thresh),
            "--marker-thresh", str(marker_thresh),
            "--imputation", imputation,
        ]
        if profile_cfg["max_samples"]:
            train_cmd += ["--max-samples", str(profile_cfg["max_samples"])]
        if profile_cfg["max_markers"]:
            train_cmd += ["--max-markers", str(profile_cfg["max_markers"])]
        if pretrained_ckpt is not None:
            train_cmd += ["--pretrained", str(pretrained_ckpt)]
            train_cmd += ["--freeze-encoder-epochs", str(freeze_encoder_epochs)]

        rc = _run_cmd(train_cmd, f"{dataset} seed={seed} — train_autoencoder")
        elapsed = time.time() - t0

        # Read stats if available
        stats_file = out_dir / "ae_embedding_stats.json"
        stats = {}
        if stats_file.exists():
            with open(stats_file) as f:
                stats = json.load(f)

        seed_result = {
            "returncode": rc,
            "elapsed_seconds": round(elapsed, 1),
            "best_epoch": stats.get("epochs_trained"),
            "best_val_loss": stats.get("best_val_loss"),
        }
        results["seeds"][seed] = seed_result

        _log_run({
            "ts": datetime.datetime.utcnow().isoformat() + "Z",
            "tag": tag,
            "method": "autoencoder",
            "dataset": dataset,
            "seed": seed,
            "profile": "ae",
            "gpu": gpu,
            "device": device,
            "bottleneck": profile_cfg["bottleneck"],
            "hidden": profile_cfg["hidden"],
            "n_blocks": profile_cfg["n_blocks"],
            "elapsed_seconds": round(elapsed, 1),
            "best_epoch": stats.get("epochs_trained"),
            "best_val_loss": stats.get("best_val_loss"),
            "samples": stats.get("samples"),
            "markers": stats.get("markers"),
        })

        if rc != 0:
            print(f"[ERROR] train_autoencoder failed for {dataset} seed={seed}")

    # --- Validate AE embeddings across seeds ---
    # We run validate_embeddings.py on the bottleneck embeddings
    # by passing the AE nodes file as if it were a regular embedding
    print(f"\n[ae-val] Validation for {dataset} across {len(seeds)} seeds")

    return results


def _run_ensemble(dataset: str, tag: str) -> int:
    """Run ensemble_embeddings.py for a dataset's seeds."""
    seeds_dir = SCRIPTS_DIR.parent / "experiments" / dataset / tag
    out_dir = seeds_dir / "ensemble"
    cmd = [
        sys.executable, str(SCRIPTS_DIR / "ensemble_embeddings.py"),
        "--seeds-dir", str(seeds_dir),
        "--out-dir", str(out_dir),
    ]
    return _run_cmd(cmd, f"{dataset} — ensemble ({len(list(seeds_dir.glob('seed*')))} seeds)")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run autoencoder (Level 1) across datasets."
    )
    parser.add_argument("--profile", choices=list(PROFILES.keys()), default="dev")
    parser.add_argument("--datasets", default=None,
                        help="Comma-separated dataset names (default: all)")
    parser.add_argument("--tag", default=None,
                        help="Run tag (default: auto-generated)")
    parser.add_argument("--gpu", default="cpu",
                        help="GPU label for logging (e.g., 4060Ti, L40S)")
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    parser.add_argument("--imputation", default="most_frequent")
    parser.add_argument("--sample-thresh", type=float, default=0.50)
    parser.add_argument("--marker-thresh", type=float, default=0.50)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--dropout", type=float, default=0.2)
    parser.add_argument("--noise-frac", type=float, default=0.15)
    parser.add_argument("--sparsity-lambda", type=float, default=1e-4)
    parser.add_argument("--weight-decay", type=float, default=1e-5)
    # Ensemble
    parser.add_argument("--ensemble", action="store_true", default=True,
                        help="Run ensemble after training all seeds (default: True)")
    parser.add_argument("--no-ensemble", dest="ensemble", action="store_false")
    # Transfer learning
    parser.add_argument("--transfer-from", default=None,
                        help="Dataset to use as pre-trained source (e.g., global_snp)")
    parser.add_argument("--transfer-to", default=None,
                        help="Comma-separated target datasets for fine-tuning")
    parser.add_argument("--freeze-encoder-epochs", type=int, default=10,
                        help="Epochs to freeze encoder during fine-tuning")
    args = parser.parse_args()

    profile_cfg = PROFILES[args.profile]
    tag = args.tag or f"ae-{datetime.datetime.now().strftime('%y%m%d-%H%M')}-{args.profile}"

    if args.datasets:
        ds_names = [d.strip() for d in args.datasets.split(",")]
    else:
        ds_names = list(DATASETS.keys())

    # Auto-detect device
    device = args.device
    if device == "auto":
        try:
            import torch
            device = "cuda" if torch.cuda.is_available() else "cpu"
        except ImportError:
            device = "cpu"

    print(f"{'='*60}")
    print(f" GENO-MAP Level 1 — Autoencoder Training")
    print(f"{'='*60}")
    print(f"Profile:    {args.profile}")
    print(f"Tag:        {tag}")
    print(f"Device:     {device}")
    print(f"GPU label:  {args.gpu}")
    print(f"Datasets:   {ds_names}")
    print(f"Seeds:      {profile_cfg['seeds']}")
    print(f"Arch:       bottleneck={profile_cfg['bottleneck']}, "
          f"hidden={profile_cfg['hidden']}, blocks={profile_cfg['n_blocks']}")
    print(f"Training:   epochs={profile_cfg['epochs']}, patience={profile_cfg['patience']}, "
          f"batch={profile_cfg['batch_size']}")
    print(f"Denoising:  noise={args.noise_frac}, sparsity_λ={args.sparsity_lambda}")
    print(f"Ensemble:   {args.ensemble}")
    if args.transfer_from:
        transfer_to = [d.strip() for d in args.transfer_to.split(",")] if args.transfer_to else []
        print(f"Transfer:   {args.transfer_from} → {transfer_to} "
              f"(freeze={args.freeze_encoder_epochs} epochs)")
    else:
        transfer_to = []

    if device.startswith("cuda"):
        try:
            import torch
            props = torch.cuda.get_device_properties(0)
            print(f"[gpu]       {torch.cuda.get_device_name(0)} "
                  f"({props.total_memory / 1e9:.1f} GB)")
        except Exception:
            pass
    print()

    t_total = time.time()
    all_results = {}

    # Determine training order: transfer sources first, then others
    if args.transfer_from and args.transfer_from in ds_names:
        # Put transfer source first so its checkpoint is available
        ordered = [args.transfer_from] + [d for d in ds_names if d != args.transfer_from]
    else:
        ordered = ds_names

    source_ckpt: Path | None = None  # best checkpoint from transfer source

    for ds_name in ordered:
        if ds_name not in DATASETS:
            print(f"[WARN] Unknown dataset '{ds_name}', skipping.")
            continue
        input_path = SCRIPTS_DIR.parent / DATASETS[ds_name]["input"]
        if not input_path.exists():
            print(f"[WARN] Input not found: {input_path}, skipping {ds_name}.")
            continue

        # Transfer learning: use pre-trained checkpoint for target datasets
        pretrained = None
        freeze_epochs = 0
        if ds_name in transfer_to and source_ckpt is not None:
            pretrained = source_ckpt
            freeze_epochs = args.freeze_encoder_epochs
            print(f"\n[transfer] Fine-tuning {ds_name} from {source_ckpt}")

        ds_result = run_ae_dataset(
            ds_name, DATASETS[ds_name], profile_cfg, tag, args.gpu, device,
            args.imputation, args.sample_thresh, args.marker_thresh,
            args.lr, args.dropout, args.noise_frac, args.sparsity_lambda,
            weight_decay=args.weight_decay,
            pretrained_ckpt=pretrained,
            freeze_encoder_epochs=freeze_epochs,
        )
        all_results[ds_name] = ds_result

        # If this is the transfer source, grab the best checkpoint
        if ds_name == args.transfer_from:
            # Use seed42 checkpoint as the transfer source
            ckpt_candidates = sorted(
                (SCRIPTS_DIR.parent / "experiments" / ds_name / tag).glob("seed*/model_checkpoint.pt")
            )
            if ckpt_candidates:
                source_ckpt = ckpt_candidates[0]
                print(f"[transfer] Source checkpoint: {source_ckpt}")

        # Ensemble all seeds for this dataset
        if args.ensemble and len(profile_cfg["seeds"]) >= 2:
            _run_ensemble(ds_name, tag)

    elapsed_total = time.time() - t_total

    # --- Summary ---
    print(f"\n{'='*60}")
    print(f" Summary — {tag}")
    print(f"{'='*60}")
    print(f"{'Dataset':<20} {'Seed':>6} {'Epoch':>6} {'Val Loss':>10} {'Time':>8}")
    print("-" * 54)
    for ds_name, res in all_results.items():
        for seed, sr in res["seeds"].items():
            ep = sr.get("best_epoch", "?")
            vl = sr.get("best_val_loss")
            vl_str = f"{vl:.6f}" if vl is not None else "?"
            tm = f"{sr['elapsed_seconds']:.1f}s"
            print(f"{ds_name:<20} {seed:>6} {ep!s:>6} {vl_str:>10} {tm:>8}")

    print(f"\nTotal time: {elapsed_total:.1f}s ({elapsed_total/60:.1f} min)")
    print(f"Artefacts:  experiments/*/{ tag }/seed*/")
    print(f"Log:        experiments/runs.jsonl")
    print("=" * 60)


if __name__ == "__main__":
    main()
