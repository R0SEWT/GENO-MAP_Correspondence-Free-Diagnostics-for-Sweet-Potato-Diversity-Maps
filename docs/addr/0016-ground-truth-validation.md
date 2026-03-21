# ADR-0016: Verificación Empírica de Namespaces Disjuntos y Validación con Ground Truth

**Estado:** Aceptada  
**Fecha:** 2026-03-07  
**Contexto:** El pipeline GENO-MAP compara mapas de diversidad mediante métricas correspondence-free (J_edge, J_nbr, SS). Todo el enfoque se justifica bajo la premisa de que los paneles usan sistemas de ID incompatibles. Se necesita verificar empíricamente esta premisa y, donde sea posible, validar las métricas contra ground truth.

## Problema

> 1. ¿Los namespaces de identificadores son realmente disjuntos entre colecciones?
> 2. Cuando existe correspondencia conocida (ground truth), ¿las métricas correspondence-based confirman que los grafos capturan estructura compartida?

## Hallazgo 1: Verificación de namespaces

| Comparación | IDs compartidos | Namespace |
|---|---|---|
| Global SNP ↔ Global SilicoDArT | 5 930 / 5 970 (99.3%) | CIP accession (`CIP 400363`) |
| LD SNP ↔ LD SilicoDArT | 476 / 635 (74.9%) | barcode\_plate\_well (`908625126003_A_7`) |
| **Global ↔ LowDensity** | **0 / ∞** | **Completamente disjuntos** |

**Conclusión**: Dentro de cada colección, SNP y SilicoDArT comparten namespace (son las mismas muestras genotipadas con dos tecnologías). Entre colecciones (Global vs LD), la intersección es vacía. El enfoque correspondence-free es **necesario** para comparar entre colecciones y **opcional** (pero útil como baseline) dentro de la misma colección.

## Hallazgo 2: Validación con Ground Truth

Se explotó el namespace compartido intra-colección para validar las métricas. Se compararon los grafos PCA-30D→kNN(k=15, cosine) de SNP vs SilicoDArT alineando muestras por su accession CIP (ground truth).

### Global: SNP ↔ SilicoDArT (n = 5 930 alineadas)

| Métrica | GT (alineado) | CF (shuffled) | INTRA (calibración) |
|---|---|---|---|
| J_edge | **0.5933** | 0.0013 | 0.9175 |
| J_nbr | **0.6220** | 0.0013 | 0.9232 |
| SS | **0.9859** | 0.0334 | — |

### LowDensity: SNP ↔ SilicoDArT (n = 630 alineadas)

| Métrica | GT (alineado) | CF (shuffled) | INTRA (calibración) |
|---|---|---|---|
| J_edge | **0.6674** | 0.0115 | 0.9243 |
| J_nbr | **0.7026** | 0.0119 | 0.9291 |
| SS | **0.9343** | 0.1011 | — |

### Interpretación

| Comparación | Interpretación |
|---|---|
| GT >> CF (×450 en J_edge) | La correspondencia revela estructura real; no es artefacto de indexación |
| SS_GT ≈ 0.99 (Global), 0.93 (LD) | Los subespacios PCA de SNP y SilicoDArT son casi idénticos — ambas tecnologías capturan la misma estructura de diversidad |
| GT < INTRA (0.59 vs 0.92) | Los grafos SNP y SilicoDArT no son idénticos (tecnologías distintas generan variación), pero capturan la misma estructura subyacente (SS ≈ 1) |
| LD J_edge_GT > Global J_edge_GT | El panel más pequeño (630 vs 5 930) produce grafos más estables entre tecnologías |

## Decisión

1. **Se confirma** que el enfoque correspondence-free es necesario para comparaciones inter-colección (Global ↔ LD).
2. **Se documenta** que las métricas de estabilidad del grafo son capaces de detectar estructura real, validadas contra ground truth con un factor ×450 sobre la null (shuffled).
3. **No se modifican** los ADRs de metodología existentes (0009-0015): ninguno requiere correspondencia muestral.
4. **Se agrega** como oportunidad de trabajo futuro: usar el ground truth intra-colección como benchmark formal de las métricas.

## Artefactos

| Entregable | Ubicación |
|---|---|
| Script del experimento | `scripts/run_adr0016_ground_truth.py` |
| Resultados JSON | `experiments/adr/adr0016_ground_truth.json` |
| Esta ADR | `docs/addr/0016-ground-truth-validation.md` |

## Semillas y parámetros

- Seeds: 42, 52, 62 (consistente con ADRs anteriores)
- PCA: 30 componentes
- kNN: k=15, métrica cosine
- Imputation: mode (most_frequent)
- Filtrado: sample_thresh=0.50, marker_thresh=0.50
