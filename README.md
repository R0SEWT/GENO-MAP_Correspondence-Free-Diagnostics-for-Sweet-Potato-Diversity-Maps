# Oxor – Diversidad genómica andina

Pipeline ligero para construir embeddings 2D y grafos kNN a partir de matrices DArT/DArTSeq (SNP y SilicoDArT) y generar insumos para el póster (nubes de diversidad, grafos de similitud, “recetas” de cruza).

## Setup rápido (uv o pip)
- Requiere Python 3.10+.
- Con uv: `uv venv && source .venv/bin/activate && uv sync`
- Con pip: `python -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt`

## Datos locales (no versionados)
- `data/10.21223P30BVZYY_Genetic_diversity/` (sweetpotato global SNP/SilicoDArT)
- `data/10.21223P3UBDJ44_LowDensity/` (sweetpotato low-density DArTSeq)
- `data/10.21223P33VYY8C_Wild/` (sweetpotato wild relatives)
- Referencias y DOIs en `data/metadata.json`. `data/` está en `.gitignore`.

## Scripts clave
- `scripts/build_embeddings.py`: lee matrices anchas, imputa faltantes, reduce dimensión (UMAP si está instalado, si no PCA), y construye grafo kNN. Salidas: `<prefix>_nodes.json`, `<prefix>_edges.json`, `<prefix>_stats.json`.
  - Preview rápido:  
    `python scripts/build_embeddings.py --input data/10.21223P30BVZYY_Genetic_diversity/SNP_Genotypes.csv --max-samples 200 --max-markers 5000 --neighbors 15 --metric cosine --out-prefix data/outputs/global_snp_preview`
  - Full run (CPU) leer todo:  
    `python scripts/build_embeddings.py --input data/10.21223P30BVZYY_Genetic_diversity/SNP_Genotypes.csv --max-samples 0 --max-markers 0 --limit-rows 0 --neighbors 20 --metric cosine --out-prefix data/outputs/global_snp_full`
  - Para SilicoDArT (binario) usa `--metric jaccard`.
  - Campos útiles en `_stats.json`: samples, markers, missing_rate, metric, neighbors, limit_rows.
- Drivers listos (`scripts/run_*.sh`) que ejecutan y registran en `experiments/runs.jsonl`:
  - `run_global_snp.sh`, `run_global_silico.sh`, `run_lowdensity_snp.sh`, `run_lowdensity_silico.sh`, `run_wild_snp.sh`
  - Overrides via env: `METRIC`, `NEIGHBORS`, `MAX_SAMPLES`, `MAX_MARKERS`, `LIMIT_ROWS`, `RUN_TAG`, `SEED`, `NOTES`.

## Qué graficar con las salidas
- `*_nodes.json`: lista de nodos con `id`, `embedding: [x, y]` y `meta` (ej. institución, accession si viene en la tabla). Úsalo para nubes 2D (UMAP/PCA) coloreadas por origen.
- `*_edges.json`: lista de aristas kNN (`source`, `target`, `distance`) para grafos de similitud. Útil en layouts ForceAtlas2/Fruchterman o para resaltar comunidades.
- `*_stats.json`: conteos (samples, markers), métrica y missing_rate para anotar en la figura o pie de gráfico.
- Visual flow: renderiza puntos del embedding; opcionalmente superpone aristas cercanas para ver grupos densos. Colorea por `meta` (p.ej. accession/origen) o por dataset. 

## Tracking de experimentos
- Guarda cada corrida en `experiments/runs.jsonl` (JSONL) con claves mínimas: `ts`, `dataset`, `input`, `metric`, `neighbors`, `limit_rows`, `max_samples`, `max_markers`, `missing_rate`, `outputs`, `notes`.
- Hay guía en `experiments/README.md` y puedes usar el snippet de ejemplo allí para registrar.

## Escalar a GPU (A100 o GPU pequeña)
- El script corre en CPU; si tienes RAPIDS, instala `cuml` y usa UMAP GPU (drop-in del paquete umap-learn). Si no, mantén PCA+UMAP CPU y solo aumenta `--max-markers`/`--max-samples`.
- Estrategia: prototipo con `--max-markers` reducido en GPU pequeña; en la A100 correr con `--limit-rows 0` y más vecinos. Ajusta `--metric` según tipo de marcador.

## Notebooks
- `notebooks/data_catalog.ipynb`: catálogo de datos reproducible → genera `data/DATA_CATALOG.md`
- `notebooks/revision_datasets.ipynb`: exploración y chequeo de los CSV descargados.

## Documentación
- `docs/agents.md` — contexto del proyecto para agentes IA (Copilot, Claude, etc.)
- `docs/addr/` — Architecture Decision Records (ADRs)
- `docs/memoria/` — trabajo de pre-curaduría (paper GENO-MAP)
- `docs/roadmap/` — ideas y planes futuros (escalado GPU, transformers, etc.)
