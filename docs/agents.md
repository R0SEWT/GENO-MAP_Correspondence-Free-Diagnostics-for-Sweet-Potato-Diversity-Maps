# 🤖 Contexto del Proyecto para Agentes IA

> Este archivo proporciona contexto estructurado para que agentes de IA
> (GitHub Copilot, Claude, etc.) trabajen eficazmente en este repositorio.

---

## Proyecto

**Nombre:** Oxor / GENO-MAP  
**Propósito:** Exploración visual de diversidad genómica en colecciones de germoplasma de cultivos andinos (camote, yacón, ulluco, oca, mashua).  
**Idioma:** Documentación, comentarios y variables en **español**. Código y nombres de función en inglés.

---

## Datos

- **Fuente:** [CIP Dataverse](https://data.cipotato.org) — Centro Internacional de la Papa
- **Formatos:** Matrices DArT/DArTSeq (SNP y SilicoDArT) en CSV
  - `sample_columns`: primera columna = etiqueta de fila, columnas restantes = muestras
  - `marker_metrics`: columnas iniciales de metadatos de marcador, seguidas de columnas numéricas
- **Catálogo:** `data/metadata.json` es la fuente de verdad para los 11 datasets curados
- **Universo:** `data/cipotato_datasets_latest100.json` — 100 datasets más recientes del CIP
- **Datos locales:** `data/` está en `.gitignore` — no versionar CSVs ni outputs

### Descarga de datos

La API del CIP Dataverse (`/api/access/datafile/{id}`) puede devolver HTTP 403 incluso
para datasets públicos. Esto es una restricción de red/tamaño, **no de acceso**.
Descargar manualmente desde https://data.cipotato.org usando el DOI.
Ver `docs/addr/0001-descarga-cip-dataverse.md`.

---

## Pipeline

```
CSV DArT/DArTSeq
  → Detección automática de formato y separador
  → Imputación (moda por defecto, mediana como alternativa)
  → PCA (min(50, p, n) componentes)
  → UMAP (2D sobre PCA)
  → Grafo kNN (k=15–20, coseno para SNP, Jaccard para SilicoDArT)
  → Artefactos JSON (nodes, edges, stats)
```

**Métricas de validación:** Trustworthiness, estabilidad inter-semilla (Jaccard),
sensibilidad a la imputación, varianza explicada PCA.

---

## Archivos Clave

| Archivo | Descripción |
|---------|-------------|
| `scripts/build_embeddings.py` | Pipeline principal: lectura → imputación → PCA → UMAP → kNN |
| `notebooks/data_catalog.ipynb` | Catálogo de datos reproducible → genera `data/DATA_CATALOG.md` |
| `data/metadata.json` | Catálogo curado: 11 datasets, 5 cultivos |
| `data/DATA_CATALOG.md` | Markdown generado con inventario completo |
| `scripts/run_*.sh` | Drivers de ejecución por dataset |

---

## Convenciones de Código

- **Python** ≥ 3.10, dependencias mínimas (`pandas`, `numpy`, `scikit-learn`, `umap-learn`, `matplotlib`)
- **Notebooks:** patrón reproducible — re-ejecutar regenera artefactos
- **Paths:** relativos desde la raíz del proyecto o desde `notebooks/` con `../data`
- **Métricas DArT:** usar nombres unificados (`CallRate`, `PIC`, `Reproducibility`, `AvgReadDepth`)
- **Seeds:** reproducibilidad con `seed=42` (variantes: 52, 62 para estabilidad)

---

## Documentación

| Carpeta | Contenido |
|---------|-----------|
| `docs/agents.md` | Este archivo — contexto para agentes IA |
| `docs/addr/` | Architecture Decision Records (ADRs) |
| `docs/memoria/` | Trabajo de pre-curaduría: paper GENO-MAP |
| `docs/roadmap/` | Ideas y planes futuros (ej: escalado GPU) |
