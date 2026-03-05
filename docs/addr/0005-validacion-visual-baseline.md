# ADR-0005: Validación Visual del Baseline poster-v2

**Estado:** Aceptada  
**Fecha:** 2026-03-01

## Contexto

El baseline **poster-v2** (PCA→UMAP→kNN, 5 seeds × 4 datasets) produjo métricas
numéricas satisfactorias (trustworthiness 0.91–0.98, PCA-Jaccard ≥ 0.89). Sin embargo,
las métricas agregadas pueden ocultar artefactos locales: clusters espurios, hubs
dominantes, outliers técnicos o distribuciones degeneradas. Se ejecutó un diagnóstico
visual completo en `notebooks/visual_diagnostics.ipynb` con 7 figuras + tabla resumen.

## Hallazgos

### 1. Estructura real en UMAP (Fig 01)

Los 4 datasets muestran clusters visuales genuinos en el espacio UMAP seed=42.
LowDensity presenta los grupos más separados (3–5 clusters claros), consistente
con menor número de muestras y mayor diversidad genética relativa. No se observan
anillos, líneas rectas u otros artefactos típicos de embeddings degenerados.

### 2. Estabilidad entre seeds (Fig 02)

Probada con 5 seeds (42, 52, 62, 72, 82) para `global_snp` y `lowdensity_snp`.
Los clusters mantienen su topología relativa; solo rotan/reflejan entre seeds.
Esto explica por qué UMAP-Jaccard (~0.22) es bajo sin que haya un problema real:
UMAP reordena vecindarios locales entre ejecuciones pero preserva la macro-estructura.

### 3. Densidad KDE (Fig 03)

Distribuciones multimodales en todos los datasets, confirmando que los clusters
no son artefactos de la proyección. Global presenta distribución más suave;
LowDensity muestra picos más marcados.

### 4. PCA scree (Fig 04)

| Dataset | Var. PC1 | Var. PC5 acum. | Var. PC50 acum. |
|---------|----------|----------------|-----------------|
| global_snp | 9.5% | 15.9% | 33.5% |
| global_silico | 10.2% | 17.1% | 35.4% |
| lowdensity_snp | 8.8% | 21.7% | 48.9% |
| lowdensity_silico | 9.5% | 22.2% | 50.6% |

La varianza está difusa en muchos componentes (típico de datos genómicos).
**Nota:** En poster-v1, `lowdensity_silico` mostraba PC1=99.99%. Esto fue
corregido al verificar el filtro de missingness; en poster-v2 reproduce valores
normales (~50% en PC50).

### 5. Distribución de grado kNN (Fig 05)

| Dataset | Grado medio | p99 | Hubs (>p99×1.1) |
|---------|-------------|-----|------------------|
| global_snp | 20.2 | 48 | 35 |
| global_silico | 20.1 | 46 | 26 |
| lowdensity_snp | 15.1 | 33 | 1 |
| lowdensity_silico | 15.1 | 34 | 4 |

Distribución razonablemente homogénea. Los hubs en global (26–35) representan
<0.6% de nodos. No hay nodos con grado patológicamente alto.

### 6. Outliers 3σ del centroide (Fig 06)

| Dataset | Outliers | % |
|---------|----------|---|
| global_snp | 140 | 2.3% |
| global_silico | 27 | 0.5% |
| lowdensity_snp | 21 | 3.3% |
| lowdensity_silico | 0 | 0.0% |

**global_snp** tiene 140 outliers que corresponden a muestras en clusters
periféricos (no aislados aleatoriamente). Son probablemente grupos genéticos
legítimos pero divergentes. **lowdensity_snp** muestra 21 outliers concentrados
en un solo grupo (IDs `908625126005_*`), posiblemente un batch técnico. Se
conservan porque podrían ser biológicamente relevantes.

### 7. Distancia de aristas kNN (Fig 07)

| Dataset | Dist. media (coseno) | Mediana |
|---------|---------------------|---------|
| global_snp | 0.048 | 0.025 |
| global_silico | 0.049 | 0.026 |
| lowdensity_snp | 0.133 | 0.089 |
| lowdensity_silico | 0.141 | 0.102 |

Distribución right-skewed (cola larga) en todos los datasets. LowDensity tiene
distancias ~3× mayores, consistente con menor densidad de muestras y mayor
diversidad inter-cluster. No se observan modos artificiales ni gaps.

## Decisión

1. **El baseline poster-v2 se acepta** como línea base limpia para comparación
   con modelos aprendidos (Nivel 1: autoencoder, Nivel 2: masked genotype).

2. **No se excluyen outliers.** Los 140 de global_snp y los 21 de lowdensity_snp
   son grupos coherentes, no errores de carga. Se documentan para revisión futura
   con metadatos biológicos.

3. **UMAP-Jaccard bajo (~0.22) no es un defecto:** la inestabilidad de UMAP solo
   afecta al embedding 2D de visualización, no al grafo kNN (construido sobre PCA,
   Jaccard ≥ 0.89).

4. **Proceder a Level 1** (autoencoder genómico) para generar embeddings aprendidos
   y comparar estabilidad/trustworthiness vs baseline PCA.

## Consecuencias

- Las figuras se archivan en `docs/figures/diagnostics/` (01–07 + tabla).
- El notebook `visual_diagnostics.ipynb` queda como artefacto reproducible.
- Los métricas de outliers/hubs servirán como benchmark para detectar mejoras
  o regresiones en los embeddings del autoencoder.
