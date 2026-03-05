# ADR-0008: Level A — Frontera de Estabilidad Representacional

**Estado:** Aceptada  
**Fecha:** 2026-03-01  
**Contexto:** Nivel A del roadmap revisado (scale.md)

## Pregunta

¿A partir de qué tamaño muestral $n^*$ las representaciones aprendidas (AE)
superan a PCA en trustworthiness? ¿El cruce es limpio o difuso?

## Diseño experimental

| Parámetro | Valor |
|-----------|-------|
| Dataset | Global SNP (n=5970, p=20069) |
| Tamaños de muestra | 14 puntos: {100, 200, 500, 1000, 2000, 3000, 3500, 4000, 4100, 4200, 4300, 4500, 5000, 5970} |
| Seeds por condición | 3 (42, 52, 62) |
| Método baseline | PCA (50 componentes) → trustworthiness k=15 en primeros 30 PCs |
| Método aprendido | AE v1 (bn64, hidden=512, 2 bloques, 200 epochs, patience=20) |
| Métricas | Trustworthiness (k=15), Edge Jaccard (estabilidad inter-seed) |
| Submuestreo | Determinístico (RandomState(0)) para reproducibilidad |

Total: 14 tamaños × 3 seeds × 2 métodos = 84 condiciones.  
Tiempo total: ~25 minutos (GPU RTX 4060 Ti).

## Resultados

### Trust vs tamaño muestral

| n | n/p | PCA Trust | AE Trust | Δ Trust | Winner |
|---|-----|-----------|----------|---------|--------|
| 100 | 0.005 | 0.939 | 0.886 | −0.053 | PCA |
| 200 | 0.010 | 0.932 | 0.870 | −0.062 | PCA |
| 500 | 0.025 | 0.962 | 0.875 | −0.087 | PCA |
| 1000 | 0.050 | 0.971 | 0.961 | −0.010 | PCA |
| 2000 | 0.100 | 0.978 | 0.927 | −0.052 | PCA |
| 3000 | 0.150 | 0.982 | 0.981 | −0.001 | PCA |
| 3500 | 0.174 | 0.983 | 0.983 | −0.000 | PCA |
| 4000 | 0.199 | 0.984 | 0.978 | −0.006 | PCA |
| 4100 | 0.204 | 0.984 | 0.984 | +0.000 | **AE** |
| 4200 | 0.209 | 0.984 | 0.984 | −0.000 | PCA |
| 4300 | 0.214 | 0.984 | 0.983 | −0.001 | PCA |
| 4500 | 0.224 | 0.985 | 0.986 | +0.002 | **AE** |
| 5000 | 0.249 | 0.985 | 0.987 | +0.002 | **AE** |
| 5970 | 0.298 | 0.986 | 0.988 | +0.001 | **AE** |

### Stability (Edge Jaccard inter-seed)

| n | PCA Stab | AE Stab | Δ Stab |
|---|----------|---------|--------|
| 100 | 0.805 | 0.394 | −0.411 |
| 500 | 0.897 | 0.208 | −0.689 |
| 1000 | 0.985 | 0.458 | −0.527 |
| 3000 | 0.987 | 0.499 | −0.488 |
| 5000 | 0.988 | 0.531 | −0.456 |
| 5970 | 0.987 | 0.524 | −0.463 |

## Hallazgos clave

### 1. Crossover n* ∈ [4000, 4500] (n/p ≈ 0.20–0.22)

El AE supera marginalmente a PCA en trustworthiness cuando n/p > 0.20.
Pero la ventaja máxima es Δ = +0.002 — estadísticamente insignificante.

### 2. La zona de transición es difusa, no abrupta

En el rango n ∈ [3000, 5970], |Δtrust| < 0.005. A n=4100 el AE gana por
+0.0001, a n=4200 pierde por −0.0004. Los métodos son esencialmente
indistinguibles en trust para n > 3000.

### 3. La estabilidad es estructuralmente inferior — siempre

AE stability ≈ 0.20–0.53 vs PCA ≈ 0.80–0.99 en **todos** los tamaños
de muestra. Incluso donde el AE "gana" en trust (n=5000), su estabilidad
es 0.53 vs PCA 0.99 (-0.46).

Este gap de estabilidad **no se cierra con más datos**. Es una propiedad
intrínseca de la inicialización aleatoria + dropout + scheduling del AE.

### 4. PCA escala mejor que AE en estabilidad

PCA stability mejora monotónicamente: 0.80 → 0.99 conforme n crece.
AE stability oscila sin tendencia clara: 0.39 → 0.21 → 0.46 → 0.53.

## Implicaciones

1. **Para n < 3000** (ratio n/p < 0.15): PCA domina en trust y estabilidad.
   No hay justificación para usar AE.

2. **Para n > 4000** (ratio n/p > 0.20): AE iguala trust de PCA pero con
   estabilidad inaceptable (~0.5 vs ~1.0). Tampoco hay justificación.

3. **El cuello de botella del AE no es n/p sino estabilidad**: Resolver
   la inestabilidad inter-seed requeriría:
   - Ensemble determinístico (pero ADR-0007 mostró que ensemble baja trust)
   - Redes más pequeñas (pero reducir capacidad baja trust)
   - Weight averaging o EMA (no probado)

4. **Regla práctica para genómica de DArT**: Si n/p < 0.20 (que incluye
   el 99% de los datasets CIP existentes), PCA es óptimo.

## Decisión

- **Level A cerrado**: La frontera de representación está mapeada.
- El crossover trust ocurre en n/p ≈ 0.20, pero no es útil porque
  la estabilidad no cruza nunca.
- **Proceder a Level B** (Marker Subsampling Robustness).

## Artefactos

- `scripts/run_stability_frontier.py` — Script del experimento
- `experiments/frontier/results_global_snp.json` — Datos completos (14 × 3 seeds)
- `experiments/frontier/frontier_global_snp.png` — Figura generada por script
- `notebooks/ae_vs_baseline.ipynb` — Sección 13, Fig 10
