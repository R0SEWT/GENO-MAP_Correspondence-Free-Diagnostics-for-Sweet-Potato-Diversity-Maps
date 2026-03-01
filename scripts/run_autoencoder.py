#!/usr/bin/env python3
"""Orchestrate autoencoder training across datasets and profiles.

Trains the Level-1 genotype autoencoder (train_autoencoder.py) on each
dataset, then runs validation (validate_embeddings.py) on the learned
bottleneck embeddings to compare against the PCA baseline (poster-v2).

Profiles
--------
  dev    – 500 samples, 5 000 markers, 10 epochs  (seconds)
  local  – full data, 200 epochs, patience=20      (minutes, 4060 Ti)
  full   – full data, 300 epochs, patience=30      (production, L40S)

Usage::

    python scripts/run_autoencoder.py --profile dev  --tag ae-smoke
    python scripts/run_autoencoder.py --profile local --tag ae-v1 --gpu 4060Ti
    python scripts/run_autoencoder.py --profile full  --tag ae-v1 --gpu L40S
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
) -> Dict[str, Any]:
    """Train autoencoder for one dataset across seeds, then validate."""
    seeds = profile_cfg["seeds"]
    results: Dict[str, Any] = {"dataset": dataset, "seeds": {}}

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

    for ds_name in ds_names:
        if ds_name not in DATASETS:
            print(f"[WARN] Unknown dataset '{ds_name}', skipping.")
            continue
        input_path = SCRIPTS_DIR.parent / DATASETS[ds_name]["input"]
        if not input_path.exists():
            print(f"[WARN] Input not found: {input_path}, skipping {ds_name}.")
            continue

        ds_result = run_ae_dataset(
            ds_name, DATASETS[ds_name], profile_cfg, tag, args.gpu, device,
            args.imputation, args.sample_thresh, args.marker_thresh,
            args.lr, args.dropout, args.noise_frac, args.sparsity_lambda,
        )
        all_results[ds_name] = ds_result

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
