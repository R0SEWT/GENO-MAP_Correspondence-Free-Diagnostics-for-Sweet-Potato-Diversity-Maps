# ADR-0006: Comparación Autoencoder Level 1 vs Baseline PCA

**Estado:** Aceptada  
**Fecha:** 2026-03-01

## Contexto

El **baseline poster-v2** (PCA→UMAP→kNN) fue validado visual y numéricamente
(ADR-0005). El siguiente paso del roadmap (Level 1) propone reemplazar PCA con un
**autoencoder denoising con bloques residuales** para obtener embeddings no-lineales
de mayor calidad.

Se entrenó `ae-v1` con la siguiente configuración:

| Parámetro       | Valor                     |
|-----------------|---------------------------|
| Arquitectura    | ResBlock × 2              |
| Bottleneck      | 64 dimensiones            |
| Hidden          | 512 neuronas              |
| Denoising       | 15% masking aleatorio     |
| Loss            | MSE + sparsity (λ=1e-4)   |
| Optimizer       | AdamW, lr=1e-3            |
| Scheduler       | ReduceLROnPlateau         |
| Early stopping  | patience=20               |
| Epochs máx.     | 200 (perfil `local`)      |
| Seeds           | 42, 52, 62                |
| GPU             | RTX 4060 Ti 16 GB         |

Tiempo total de entrenamiento: **78.7 minutos** (4 datasets × 3 seeds).

## Resultados

### Trustworthiness (k=15, media ± std sobre 3 seeds)

| Dataset              | Baseline PCA 2D  | AE Bottleneck 64D | AE UMAP 2D | Δ (AE bn − PCA) |
|----------------------|-------------------|--------------------|-------------|------------------|
| Global SNP           | 0.976 ± 0.002    | **0.988 ± 0.001**  | 0.962       | **+0.012**       |
| Global SilicoDArT    | 0.979 ± 0.001    | **0.989 ± 0.001**  | 0.973       | **+0.010**       |
| LowDensity SNP       | **0.933 ± 0.002**| 0.836 ± 0.024      | 0.793       | **−0.097**       |
| LowDensity SilicoDArT| **0.911 ± 0.005**| 0.829 ± 0.022      | 0.795       | **−0.082**       |

### Estabilidad inter-seed (Jaccard edge overlap, k=15)

| Dataset              | Baseline PCA | AE Bottleneck |
|----------------------|--------------|---------------|
| Global SNP           | 0.888        | 0.524         |
| Global SilicoDArT    | 0.885        | 0.498         |
| LowDensity SNP       | 0.909        | 0.163         |
| LowDensity SilicoDArT| 0.908       | 0.208         |

### Curvas de entrenamiento

- Convergencia suave en los 4 datasets, sin overfitting visible.
- LowDensity converge en ~50–95 epochs (early stop); Global usa ~198 epochs.
- Val loss final: global_snp ≈ 0.274, global_silico ≈ 0.085, lowdensity_snp ≈ 0.067,
  lowdensity_silico ≈ 0.049.

### Distribución de distancias kNN

- AE produce distancias inter-vecinos más compactas (media coseno menor).
- LowDensity AE muestra distancias extremadamente pequeñas (μ ≈ 0.001–0.006),
  síntoma de colapso parcial del espacio de embeddings.

## Diagnóstico

### ¿Por qué AE gana en Global pero pierde en LowDensity?

1. **Ratio muestras/parámetros**: Global tiene ~5970 muestras vs LowDensity ~630.
   El autoencoder tiene ~1.2M parámetros — con 630 muestras el ratio es ~1:1900,
   insuficiente para aprender representaciones generalizables.
2. **Colapso parcial**: Las distancias de LowDensity AE muestran concentración
   excesiva cerca de cero, indicando que el bottleneck colapsa muestras distintas
   a representaciones casi idénticas.
3. **Estabilidad baja**: Cada seed inicializa pesos diferentes → representaciones
   diferentes. PCA es determinista; el AE no.

### Implicaciones para el roadmap

- Para datasets **grandes** (≥5000 muestras), el AE bottleneck supera al baseline
  en trustworthiness (+1 pp).
- Para datasets **pequeños** (<1000 muestras), PCA sigue siendo superior.
- La **estabilidad** del AE es un problema estructural que requiere soluciones
  arquitectónicas (ver Decisión).

## Decisión

1. **No reemplazar el baseline con ae-v1 de forma universal.**
   El AE Level 1 es un avance marginal para Global (+1 pp trust) pero una
   regresión significativa para LowDensity (−8 a −10 pp).

2. **Aceptar ae-v1 como "prueba de concepto validada"** para confirmar que
   autoencoders denoising pueden producir embeddings competitivos sobre datasets
   genómicos grandes.

3. **Próximos pasos concretos** antes de Level 2 (Transformer):
   - **Regularización LowDensity**: probar dropout más agresivo, reducir bottleneck
     a 32D, o usar weight decay mayor para evitar colapso.
   - **Ensemble / media de embeddings**: promediar bottleneck de múltiples seeds
     para mejorar estabilidad.
   - **Transfer learning**: pre-entrenar en Global, fine-tune en LowDensity.
   - **Data augmentation**: generar variantes sintéticas (bootstrap genómico).

4. **Conservar el baseline PCA como referencia canónica** en el paper hasta que
   un modelo Level ≥ 2 lo supere en todos los datasets y métricas.

## Artefactos

| Artefacto | Ubicación |
|-----------|-----------|
| Notebook de comparación | `notebooks/ae_vs_baseline.ipynb` |
| Fig 01: Trustworthiness | `docs/figures/comparison/01_trustworthiness_comparison.png` |
| Fig 02: Estabilidad | `docs/figures/comparison/02_stability_comparison.png` |
| Fig 03: UMAP side-by-side | `docs/figures/comparison/03_umap_side_by_side.png` |
| Fig 04: Distribución distancias | `docs/figures/comparison/04_distance_comparison.png` |
| Fig 05: Curvas de entrenamiento | `docs/figures/comparison/05_training_curves.png` |
| Runs AE ae-v1 | `experiments/*/ae-v1/seed{42,52,62}/` |
| Orchestrador AE | `scripts/run_autoencoder.py` |
