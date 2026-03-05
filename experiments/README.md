# Tracking de experimentos

Recomendación: un único archivo JSONL (`experiments/runs.jsonl`) con una línea por corrida. Claves mínimas sugeridas:
- `ts`: ISO timestamp
- `dataset`: nombre corto (p.ej. global_snp, lowdensity_snp, wild_silico)
- `input`: ruta del CSV
- `format`: sample_columns | marker_metrics
- `metric`: euclidean | cosine | jaccard
- `neighbors`, `max_samples`, `max_markers`, `limit_rows`, `seed`
- `samples`, `markers`, `missing_rate` (del `_stats.json`)
- `outputs`: prefix o rutas de `*_nodes.json`/`*_edges.json`
- `notes`: decisiones o hallazgos

Snippet para registrar (editando los valores que correspondan):
```bash
python - <<'PY'
import json, datetime, pathlib
run = {
    "ts": datetime.datetime.utcnow().isoformat() + "Z",
    "dataset": "global_snp",
    "input": "data/10.21223P30BVZYY_Genetic_diversity/SNP_Genotypes.csv",
    "format": "sample_columns",
    "metric": "cosine",
    "neighbors": 15,
    "max_samples": 200,
    "max_markers": 5000,
    "limit_rows": 5010,
    "seed": 42,
    "samples": 200,
    "markers": 5000,
    "missing_rate": 0.06,
    "outputs": "data/outputs/global_snp_preview",
    "notes": "Preview para layout póster"
}
path = pathlib.Path("experiments/runs.jsonl")
path.parent.mkdir(parents=True, exist_ok=True)
with path.open("a", encoding="utf-8") as f:
    f.write(json.dumps(run, ensure_ascii=True) + "\n")
print("logged to", path)
PY
```

Sugerencias:
- Usa un `dataset` consistente para facilitar filtros.
- Si repites corridas con distinto `seed`, agrega `run_id` o `tag`.
- No guardes datos brutos aquí (ya están en `data/`); solo metadata de corridas.
