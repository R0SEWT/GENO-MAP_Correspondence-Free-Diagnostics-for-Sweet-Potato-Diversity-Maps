#!/usr/bin/env python3
"""Orchestrate Masked Genotype Transformer training across datasets.

Runs train_transformer.py on each dataset with multiple seeds,
following the same pattern as run_autoencoder.py.

Profiles
--------
  dev       – 500 samples, 5 000 markers, 10 epochs  (seconds)
  local     – full data, 50 epochs, patience=15       (4060 Ti)
  full      – full data, 100 epochs, patience=25      (L40S)

Usage::

    python scripts/run_transformer.py --profile dev   --tag mgt-smoke
    python scripts/run_transformer.py --profile local --tag mgt-v1 --gpu 4060Ti
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

# ── Dataset catalogue ──

DATASETS: Dict[str, Dict[str, Any]] = {
    "global_snp": {
        "input": "data/10.21223P30BVZYY_Genetic_diversity/SNP_Genotypes.csv",
        "n_values": 3,
    },
    "global_silico": {
        "input": "data/10.21223P30BVZYY_Genetic_diversity/SilicoDArT_Genotypes.csv",
        "n_values": 2,
    },
    "lowdensity_snp": {
        "input": "data/10.21223P3UBDJ44_LowDensity/01_Report_DSp25-515_SNPs_Filtered_by _reads.csv",
        "n_values": 3,
    },
    "lowdensity_silico": {
        "input": "data/10.21223P3UBDJ44_LowDensity/02_Report_DSp25-515_Silico-DArT_Filtered_by_reads.csv",
        "n_values": 2,
    },
}

PROFILES: Dict[str, Dict[str, Any]] = {
    "dev": {
        "max_samples": 500,
        "max_markers": 5000,
        "epochs": 10,
        "patience": 5,
        "d_model": 64,
        "n_heads": 4,
        "n_layers": 2,
        "d_ff": 128,
        "bottleneck": 32,
        "context_len": 256,
        "batch_size": 64,
        "n_mc_subsets": 5,
        "seeds": [42],
    },
    "local": {
        "max_samples": 0,
        "max_markers": 0,
        "epochs": 50,
        "patience": 15,
        "d_model": 128,
        "n_heads": 4,
        "n_layers": 4,
        "d_ff": 256,
        "bottleneck": 64,
        "context_len": 512,
        "batch_size": 64,
        "n_mc_subsets": 20,
        "seeds": [42, 52, 62],
    },
    "full": {
        "max_samples": 0,
        "max_markers": 0,
        "epochs": 100,
        "patience": 25,
        "d_model": 256,
        "n_heads": 8,
        "n_layers": 6,
        "d_ff": 512,
        "bottleneck": 64,
        "context_len": 1024,
        "batch_size": 64,
        "n_mc_subsets": 30,
        "seeds": [42, 52, 62],
    },
}

# Per-dataset overrides for small datasets
DS_OVERRIDES: Dict[str, Dict[str, Any]] = {
    "lowdensity_snp": {"dropout": 0.2},
    "lowdensity_silico": {"dropout": 0.2},
}

ROOT = Path(__file__).resolve().parent.parent


def _run_single(
    ds_name: str,
    ds_conf: Dict[str, Any],
    profile: Dict[str, Any],
    tag: str,
    seed: int,
    extra_args: List[str],
) -> Dict[str, Any]:
    """Run train_transformer.py for one dataset+seed."""
    out_dir = ROOT / "experiments" / ds_name / tag / f"seed{seed}"
    out_dir.mkdir(parents=True, exist_ok=True)

    overrides = DS_OVERRIDES.get(ds_name, {})

    cmd = [
        sys.executable,
        str(ROOT / "scripts" / "train_transformer.py"),
        "--input", str(ROOT / ds_conf["input"]),
        "--out-dir", str(out_dir),
        "--d-model", str(profile["d_model"]),
        "--n-heads", str(profile["n_heads"]),
        "--n-layers", str(profile["n_layers"]),
        "--d-ff", str(profile["d_ff"]),
        "--bottleneck", str(profile["bottleneck"]),
        "--context-len", str(profile["context_len"]),
        "--batch-size", str(profile["batch_size"]),
        "--epochs", str(profile["epochs"]),
        "--patience", str(profile["patience"]),
        "--seed", str(seed),
        "--n-values", str(ds_conf["n_values"]),
        "--n-mc-subsets", str(profile["n_mc_subsets"]),
        "--dropout", str(overrides.get("dropout", 0.1)),
    ]
    if profile.get("max_samples"):
        cmd += ["--max-samples", str(profile["max_samples"])]
    if profile.get("max_markers"):
        cmd += ["--max-markers", str(profile["max_markers"])]
    cmd += extra_args

    print(f"\n{'='*70}")
    print(f"[run] {ds_name}  seed={seed}  tag={tag}")
    print(f"[run] {' '.join(cmd[-10:])}")
    print(f"{'='*70}\n")

    t0 = time.time()
    result = subprocess.run(cmd, cwd=str(ROOT))
    elapsed = round(time.time() - t0, 1)

    return {
        "dataset": ds_name,
        "seed": seed,
        "tag": tag,
        "out_dir": str(out_dir),
        "returncode": result.returncode,
        "elapsed_seconds": elapsed,
    }


def _log_run(run_info: Dict[str, Any]) -> None:
    """Append run metadata to experiments/runs.jsonl."""
    run_info["timestamp"] = datetime.datetime.now().isoformat()
    log_path = ROOT / "experiments" / "runs.jsonl"
    with open(log_path, "a") as f:
        f.write(json.dumps(run_info, ensure_ascii=False) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Orchestrate MGT training across datasets."
    )
    parser.add_argument(
        "--profile",
        required=True,
        choices=list(PROFILES.keys()),
        help="Training profile",
    )
    parser.add_argument("--tag", required=True, help="Experiment tag")
    parser.add_argument("--gpu", default="", help="GPU label for logging")
    parser.add_argument(
        "--datasets",
        nargs="*",
        default=None,
        help="Subset of datasets (default: all)",
    )
    args, extra = parser.parse_known_args()

    profile = PROFILES[args.profile]
    ds_names = args.datasets or list(DATASETS.keys())

    print(f"[run] Profile: {args.profile}  Tag: {args.tag}  GPU: {args.gpu or 'auto'}")
    print(f"[run] Datasets: {ds_names}")
    print(f"[run] Seeds: {profile['seeds']}")
    print(f"[run] Arch: d={profile['d_model']} h={profile['n_heads']} L={profile['n_layers']} "
          f"ff={profile['d_ff']} ctx={profile['context_len']} bn={profile['bottleneck']}")

    total_t0 = time.time()
    runs: List[Dict[str, Any]] = []

    for ds_name in ds_names:
        if ds_name not in DATASETS:
            print(f"[warn] Unknown dataset: {ds_name}, skipping")
            continue
        ds_conf = DATASETS[ds_name]

        for seed in profile["seeds"]:
            run_info = _run_single(ds_name, ds_conf, profile, args.tag, seed, extra)
            run_info["gpu"] = args.gpu
            run_info["profile"] = args.profile
            _log_run(run_info)
            runs.append(run_info)

    total_elapsed = round(time.time() - total_t0, 1)

    # ── Summary ──
    print(f"\n{'='*70}")
    print(f"[run] COMPLETE — {len(runs)} runs in {total_elapsed}s")
    print(f"{'='*70}")
    for r in runs:
        status = "OK" if r["returncode"] == 0 else f"FAIL({r['returncode']})"
        print(f"  {r['dataset']:25s} seed={r['seed']}  {status:6s}  {r['elapsed_seconds']:6.1f}s")


if __name__ == "__main__":
    main()
