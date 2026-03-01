#!/usr/bin/env bash
# Run embedding + kNN graph for sweetpotato wild relatives SNP dataset.
set -euo pipefail

DATASET="wild_snp"
INPUT="data/10.21223P33VYY8C_Wild/01_Report_DSp25-515_SNPs_Filtered_by _reads.csv"
METRIC="${METRIC:-cosine}"
NEIGHBORS="${NEIGHBORS:-15}"
MAX_SAMPLES="${MAX_SAMPLES:-0}"
MAX_MARKERS="${MAX_MARKERS:-0}"
LIMIT_ROWS="${LIMIT_ROWS:-0}"
SEED="${SEED:-42}"
RUN_TAG="${RUN_TAG:-$(date +%y%m%d-%H%M)}"
OUT_PREFIX="experiments/${DATASET}/${RUN_TAG}-${DATASET}"
mkdir -p "$(dirname "$OUT_PREFIX")"

python scripts/build_embeddings.py \
  --input "$INPUT" \
  --metric "$METRIC" \
  --neighbors "$NEIGHBORS" \
  --max-samples "$MAX_SAMPLES" \
  --max-markers "$MAX_MARKERS" \
  --limit-rows "$LIMIT_ROWS" \
  --seed "$SEED" \
  --out-prefix "$OUT_PREFIX"

python - <<PY
import json, datetime, pathlib, os
stats = json.load(open(f"{os.environ['OUT_PREFIX']}_stats.json"))
run = {
    "ts": datetime.datetime.utcnow().isoformat() + "Z",
    "dataset": os.environ["DATASET"],
    "input": os.environ["INPUT"],
    "format": stats.get("format"),
    "metric": os.environ["METRIC"],
    "neighbors": int(os.environ["NEIGHBORS"]),
    "max_samples": int(os.environ["MAX_SAMPLES"]),
    "max_markers": int(os.environ["MAX_MARKERS"]),
    "limit_rows": int(os.environ["LIMIT_ROWS"]),
    "seed": int(os.environ["SEED"]),
    "samples": stats["samples"],
    "markers": stats["markers"],
    "missing_rate": stats["missing_rate"],
    "outputs": os.environ["OUT_PREFIX"],
    "notes": os.environ.get("NOTES", "")
}
path = pathlib.Path("experiments/runs.jsonl")
path.parent.mkdir(parents=True, exist_ok=True)
with path.open("a", encoding="utf-8") as f:
    f.write(json.dumps(run, ensure_ascii=True) + "\n")
print("logged:", run)
PY
