#!/usr/bin/env python3
"""Orchestrate embedding + validation runs across datasets and profiles.

Profiles
--------
  dev    – 500 samples, 5 000 markers, 1 seed  (seconds)
  local  – full data, 3 seeds   (minutes, 4060 Ti)
  full   – full data, 5 seeds   (production, L40S)

Usage::

    python scripts/run_experiments.py --profile dev  --tag smoke-test
    python scripts/run_experiments.py --profile local --tag v2-baseline
    python scripts/run_experiments.py --profile full  --tag poster-v2 --gpu L40S
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
# Dataset catalogue  (ADR-0002: excludes Wild SNP)
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
    # wild_snp EXCLUDED — identical to lowdensity_snp (ADR-0002)
}

PROFILES: Dict[str, Dict[str, Any]] = {
    "dev": {
        "max_samples": 500,
        "max_markers": 5000,
        "seeds": [42],
        "skip_stability": True,
        "skip_sensitivity": True,
    },
    "local": {
        "max_samples": 0,
        "max_markers": 0,
        "seeds": [42, 52, 62],
        "skip_stability": False,
        "skip_sensitivity": False,
    },
    "full": {
        "max_samples": 0,
        "max_markers": 0,
        "seeds": [42, 52, 62, 72, 82],
        "skip_stability": False,
        "skip_sensitivity": False,
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


def _log_run(
    tag: str,
    dataset: str,
    seed: int,
    profile: str,
    gpu: str,
    out_prefix: Path,
    elapsed: float,
    stats_file: Path | None,
) -> None:
    """Append a line to experiments/runs.jsonl."""
    run_entry: Dict[str, Any] = {
        "ts": datetime.datetime.utcnow().isoformat() + "Z",
        "tag": tag,
        "dataset": dataset,
        "seed": seed,
        "profile": profile,
        "gpu": gpu,
        "out_prefix": str(out_prefix),
        "elapsed_seconds": round(elapsed, 1),
    }
    if stats_file and stats_file.exists():
        with open(stats_file, "r") as f:
            stats = json.load(f)
        run_entry["samples"] = stats.get("samples")
        run_entry["markers"] = stats.get("markers")
        run_entry["missing_rate_raw"] = stats.get("missing_rate_raw")
        run_entry["missing_rate_after_filter"] = stats.get("missing_rate_after_filter")

    log_path = SCRIPTS_DIR.parent / "experiments" / "runs.jsonl"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "a", encoding="utf-8") as f:
        f.write(json.dumps(run_entry, ensure_ascii=True) + "\n")


def run_dataset(
    dataset: str,
    ds_cfg: Dict[str, Any],
    profile_cfg: Dict[str, Any],
    tag: str,
    gpu: str,
    sample_thresh: float,
    marker_thresh: float,
    imputation: str,
    device: str = "cpu",
) -> None:
    """Run build_embeddings + validate_embeddings for one dataset."""
    seeds = profile_cfg["seeds"]
    max_samples = profile_cfg["max_samples"]
    max_markers = profile_cfg["max_markers"]

    for seed in seeds:
        out_dir = SCRIPTS_DIR.parent / "experiments" / dataset / tag
        out_prefix = out_dir / f"seed{seed}_{dataset}"
        out_dir.mkdir(parents=True, exist_ok=True)

        t0 = time.time()

        # --- build embeddings ---
        build_cmd = [
            sys.executable, str(SCRIPTS_DIR / "build_embeddings.py"),
            "--input", ds_cfg["input"],
            "--metric", ds_cfg["metric"],
            "--neighbors", str(ds_cfg["neighbors"]),
            "--max-samples", str(max_samples),
            "--max-markers", str(max_markers),
            "--limit-rows", "0",
            "--seed", str(seed),
            "--sample-thresh", str(sample_thresh),
            "--marker-thresh", str(marker_thresh),
            "--imputation", imputation,
            "--device", device,
            "--out-prefix", str(out_prefix),
        ]
        rc = _run_cmd(build_cmd, f"{dataset} seed={seed} — build_embeddings")
        if rc != 0:
            print(f"[ERROR] build_embeddings failed for {dataset} seed={seed}")
            continue

        elapsed_build = time.time() - t0
        stats_file = out_dir / f"{out_prefix.name}_stats.json"
        _log_run(tag, dataset, seed, profile_cfg.__class__.__name__ if hasattr(profile_cfg, '__class__') else "custom",
                 gpu, out_prefix, elapsed_build, stats_file)

    # --- validate (uses all seeds together) ---
    val_prefix = out_dir / f"all_{dataset}"
    val_cmd = [
        sys.executable, str(SCRIPTS_DIR / "validate_embeddings.py"),
        "--input", ds_cfg["input"],
        "--metric", ds_cfg["metric"],
        "--neighbors", str(ds_cfg["neighbors"]),
        "--seeds", ",".join(str(s) for s in seeds),
        "--sample-thresh", str(sample_thresh),
        "--marker-thresh", str(marker_thresh),
        "--imputation", imputation,
        "--device", device,
        "--out-prefix", str(val_prefix),
    ]
    if max_samples:
        val_cmd += ["--max-samples", str(max_samples)]
    if max_markers:
        val_cmd += ["--max-markers", str(max_markers)]
    if profile_cfg.get("skip_stability"):
        val_cmd.append("--skip-stability")
    if profile_cfg.get("skip_sensitivity"):
        val_cmd.append("--skip-sensitivity")

    _run_cmd(val_cmd, f"{dataset} — validate_embeddings")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run experiments across datasets.")
    parser.add_argument("--profile", choices=list(PROFILES.keys()), default="dev")
    parser.add_argument("--datasets", default=None,
                        help="Comma-separated dataset names (default: all)")
    parser.add_argument("--tag", default=None,
                        help="Run tag (default: auto-generated)")
    parser.add_argument("--gpu", default="cpu",
                        help="GPU label for logging (e.g., 4060Ti, L40S)")
    parser.add_argument("--sample-thresh", type=float, default=0.50)
    parser.add_argument("--marker-thresh", type=float, default=0.50)
    parser.add_argument("--imputation", default="most_frequent")
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"],
                        help="Compute device (auto detects CUDA)")
    args = parser.parse_args()

    profile_cfg = PROFILES[args.profile]
    tag = args.tag or f"{datetime.datetime.now().strftime('%y%m%d-%H%M')}-{args.profile}"

    if args.datasets:
        ds_names = [d.strip() for d in args.datasets.split(",")]
    else:
        ds_names = list(DATASETS.keys())

    # Auto-detect compute device
    device = args.device
    if device == "auto":
        try:
            import torch
            device = "cuda" if torch.cuda.is_available() else "cpu"
        except ImportError:
            device = "cpu"

    print(f"Profile: {args.profile}  |  Tag: {tag}  |  GPU: {args.gpu}  |  Device: {device}")
    print(f"Datasets: {ds_names}")
    print(f"Seeds: {profile_cfg['seeds']}")
    if device.startswith("cuda"):
        try:
            import torch
            print(f"[gpu] {torch.cuda.get_device_name(0)}  "
                  f"({torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB)")
        except Exception:
            pass
    print()

    t_total = time.time()
    for ds_name in ds_names:
        if ds_name not in DATASETS:
            print(f"[WARN] Unknown dataset '{ds_name}', skipping.")
            continue
        input_path = SCRIPTS_DIR.parent / DATASETS[ds_name]["input"]
        if not input_path.exists():
            print(f"[WARN] Input not found: {input_path}, skipping {ds_name}.")
            continue
        run_dataset(
            ds_name, DATASETS[ds_name], profile_cfg, tag, args.gpu,
            args.sample_thresh, args.marker_thresh, args.imputation,
            device=device,
        )

    elapsed_total = time.time() - t_total
    print(f"\n{'='*60}")
    print(f"All done.  Total: {elapsed_total:.1f}s")
    print(f"Artefacts: experiments/*/{ tag }/")
    print(f"Log: experiments/runs.jsonl")
    print("=" * 60)


if __name__ == "__main__":
    main()
