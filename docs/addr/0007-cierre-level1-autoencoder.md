# ADR-0007: Cierre del Level 1 — Autoencoder Genómico

**Estado:** Aceptada  
**Fecha:** 2026-03-01  
**Supersede:** Extiende ADR-0006

## Contexto

ADR-0006 reportó los resultados iniciales del autoencoder Level 1 (`ae-v1`, 200
epochs). El AE superó al baseline PCA en trustworthiness para datasets grandes
(Global, +1 pp) pero empeoró significativamente en datasets pequeños (LowDensity,
−8 a −10 pp). La estabilidad inter-seed fue muy inferior al baseline.

Se propusieron tres mejoras:
1. Regularización más agresiva para LowDensity
2. Ensemble de bottleneck entre seeds
3. Transfer learning (Global → LowDensity)

Este ADR documenta los resultados de `ae-v2` (50 epochs) con las tres mejoras
implementadas, y cierra la exploración del Level 1.

## Diseño experimental ae-v2

### Perfil `local_v2`

| Parámetro | Valor |
|-----------|-------|
| Epochs máx. | 50 (basado en análisis de curvas ae-v1: convergencia < 50 ep) |
| Patience | 15 |
| Bottleneck | 64D |
| Hidden | 512 |
| ResBlocks | 2 |
| Seeds | 42, 52, 62 |

### Regularización per-dataset

| Parámetro | Global (5970 muestras) | LowDensity (630 muestras) |
|-----------|------------------------|---------------------------|
| Dropout | 0.20 | **0.45** |
| Noise mask | 0.15 | **0.30** |
| Weight decay | 1e-5 | **5e-4** |
| Sparsity λ | 1e-4 | **5e-4** |

### Ensemble

Promedio aritmético de vectores bottleneck (64D) entre los 3 seeds, seguido de
UMAP 2D + kNN sobre el bottleneck promediado. Script: `ensemble_embeddings.py`.

### Transfer learning

1. Pre-entrenar AE en Global SNP (5970 muestras, 20k markers)
2. Cargar checkpoint `seed42` como base
3. Fine-tune en LowDensity (630 muestras, 38k–62k markers)
4. Freeze encoder por 10 epochs, luego unfreeze con lr reducido (×0.1)
5. Las capas con dimensiones incompatibles (primera/última Dense) se re-inicializan

## Resultados

### Trustworthiness (k=15, bottleneck 64D)

| Dataset | Baseline PCA | ae-v1 (200 ep) | ae-v2 reg (50 ep) | ae-v2 ensemble | ae-v2 transfer | ae-v2 TL ens |
|---------|-------------|----------------|--------------------|--------------------|-----------------|-----------|
| Global SNP | 0.976 | **0.988** | 0.982 | 0.985 | — | — |
| Global SilicoDArT | 0.979 | **0.989** | 0.977 | 0.980 | — | — |
| LowDensity SNP | **0.933** | 0.836 | 0.844 | 0.784 | 0.658 | 0.767 |
| LowDensity SilicoDArT | **0.911** | 0.829 | 0.853 | 0.770 | 0.772 | 0.823 |

### Estabilidad kNN (Edge Jaccard)

| Dataset | Baseline PCA | ae-v1 s2s | ae-v2 s2s | ae-v2 s↔ens |
|---------|-------------|-----------|-----------|-------------|
| Global SNP | 0.889 | 0.524 | 0.438 | 0.568 |
| Global SilicoDArT | 0.885 | 0.498 | 0.333 | 0.466 |
| LowDensity SNP | 0.909 | 0.163 | 0.173 | 0.146 |
| LowDensity SilicoDArT | 0.908 | 0.207 | 0.233 | 0.178 |

### Val loss ae-v2

| Dataset | Seed 42 | Seed 52 | Seed 62 | Media |
|---------|---------|---------|---------|-------|
| Global SNP | 0.377 | 0.380 | 0.382 | 0.380 (vs ae-v1 0.274) |
| Global SilicoDArT | 0.115 | 0.118 | 0.118 | 0.117 (vs ae-v1 0.085) |
| LowDensity SNP | 0.068 | 0.068 | 0.077 | 0.071 (vs ae-v1 0.067) |
| LowDensity SilicoDArT | 0.054 | 0.052 | 0.053 | 0.053 (vs ae-v1 0.049) |

## Análisis por mejora

### 1. Regularización agresiva → Mejora marginal (+2 pp)

- LowDensity SilicoDArT: 0.853 vs 0.829 de ae-v1 (**+2.4 pp**)
- LowDensity SNP: 0.844 vs 0.836 de ae-v1 (**+0.8 pp**)
- La regularización reduce el colapso pero no lo elimina
- Con 630 muestras y ~1.2M parámetros, el problema es fundamentalmente de volumen

### 2. Ensemble de bottleneck → Perjudica la trustworthiness

- Global SNP: 0.985 (ensemble) vs 0.988 (ae-v1 individual) → **−0.3 pp**
- LowDensity SNP: 0.784 (ensemble) vs 0.844 (ae-v2 reg) → **−6 pp**
- **Diagnóstico**: Promediar bottleneck suaviza las representaciones, borrando la
  estructura local fina que cada seed captura individualmente. El promedio converge
  a una representación "consenso" que pierde la discriminación de vecindarios.
- La estabilidad seed↔ensemble sube para Global (0.57 vs 0.44 s2s) pero no compensa
  la pérdida en trustworthiness.

### 3. Transfer learning → Perjudica significativamente

- LowDensity SNP transfer: **0.658** vs 0.844 regularizado → **−18.6 pp**
- LowDensity SilicoDArT transfer: **0.772** vs 0.853 regularizado → **−8.1 pp**
- **Causa raíz**: Los markers son diferentes entre Global SNP (20k loci) y
  LowDensity SNP (62k loci). Solo las capas internas (ResBlocks, bottleneck→hidden)
  se transfieren; las capas input/output se re-inicializan. La representación
  interna pre-entrenada en 20k markers no transfiere bien a 62k markers.
- El freeze de encoder por 10 epochs empeora aún más: fuerza representaciones
  internas irrelevantes para los nuevos datos.

## Decisión

### Cierre del Level 1 Autoencoder

1. **ae-v1 (200 epochs)** queda como la mejor configuración AE Level 1:
   - Mejor trust para Global: 0.988–0.989 (+1 pp vs baseline)
   - Trust para LowDensity: 0.829–0.836 (−8 a −10 pp vs baseline)
   - Tiempo de entrenamiento: ~79 min (4 datasets × 3 seeds)

2. **Las tres mejoras de ae-v2 no justifican su adopción**:
   - Regularización: mejora marginal, no resuelve el problema fundamental
   - Ensemble: daña la calidad del embedding
   - Transfer learning: incompatible entre datasets con markers distintos

3. **Para el short paper**:
   - **Sección de resultados**: reportar ae-v1 como exploración Level 1, destacar
     que supera PCA solo en datasets grandes
   - **Sección de discusión**: documentar que el ratio muestras/parámetros limita
     la aplicabilidad de autoencoders en datasets genómicos pequeños
   - **Tabla de paper**: incluir baseline PCA + ae-v1 bottleneck como las dos
     condiciones principales
   - **Figuras para paper**: Fig 06 (trust comparativo) y Fig 07 (estabilidad)

4. **El baseline PCA→UMAP→kNN se mantiene como pipeline canónico** del proyecto.

5. **Para Level 2 (Masked Genotype Transformer)**: considerar:
   - Modelos más pequeños para datasets <1000 muestras
   - Arquitecturas que compartan tokenización entre datasets
   - Evaluación previa de si el cuello de botella es muestras o markers

## Lecciones aprendidas

1. **Regularización no sustituye datos**: Con n=630 y p=60k, ninguna cantidad de
   dropout/noise/weight-decay compensa la falta de muestras.
2. **El ensemble aritmético daña estructura local**: Para embeddings de calidad,
   mejor seleccionar el seed con mejor val_loss que promediar.
3. **Transfer learning requiere alineamiento de features**: Si los markers difieren
   entre datasets, transferir pesos internos del encoder es contraproducente.
4. **50 epochs suficientes para evaluación**: Las curvas de ae-v1 muestran que la
   convergencia ocurre antes de epoch 50 en datasets pequeños y ~100 en grandes.
   Para búsqueda de hiperparámetros, 50 epochs es un proxy razonable.
5. **PCA es sorprendentemente competitivo**: La trustworthiness de PCA 2D (0.91–0.98)
   es difícil de superar con modelos más complejos, especialmente con pocas muestras.

## Artefactos

| Artefacto | Ubicación |
|-----------|-----------|
| ADR ae-v1 | `docs/addr/0006-comparacion-ae-vs-baseline.md` |
| Notebook completo | `notebooks/ae_vs_baseline.ipynb` |
| Fig 06: Trust multi-method | `docs/figures/comparison/06_trust_v2_comparison.png` |
| Fig 07: Estabilidad multi-method | `docs/figures/comparison/07_stability_v2_comparison.png` |
| Script entrenamiento | `scripts/train_autoencoder.py` |
| Script ensemble | `scripts/ensemble_embeddings.py` |
| Script orquestador | `scripts/run_autoencoder.py` |
| Runs ae-v1 | `experiments/*/ae-v1/seed{42,52,62}/` |
| Runs ae-v2 | `experiments/*/ae-v2/{seed*,ensemble*}/` |
| Log de runs | `experiments/runs.jsonl` |
