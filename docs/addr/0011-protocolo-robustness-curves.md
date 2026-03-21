# ADR-0011: Protocolo de Robustness Curves (Prioridad 2)

**Estado:** Aceptada  
**Fecha:** 2026-03-06  
**Contexto:** Las robustness curves son la contribución principal del paper (OE2). Este ADR formaliza qué perturbaciones se aplican, qué métricas se reportan, y qué protocolos de ejecución se siguen.

## Problema

Las curvas de robustez ya están implementadas (Tablas 6–8), pero el protocolo no está documentado de forma unificada. Necesitamos un estándar que cubra las tres perturbaciones base y sea extensible a block removal (ADR-0013).

## Decisión

### Perturbaciones base (imprescindibles)

| Perturbación | Niveles | Descripción |
|-------------|---------|-------------|
| **Marker subsampling** | {5, 10, 20, 50, 80}% retención | Eliminación aleatoria uniforme de columnas (marcadores) |
| **MCAR injection** | +{0, 5, 10, 20}% | Inyección de NaN en posiciones aleatorias sobre datos numéricos, seguida de re-imputación |
| **Imputation comparison** | mode vs. median | Mismos datos, distinta estrategia de imputación |

### Métricas reportadas (por perturbación)

Para cada condición, comparando contra baseline (pipeline completo sin perturbación):

| Métrica | Símbolo | Definición (ver ADR-0010) |
|---------|---------|--------------------------|
| Neighbor stability | $J_{\text{nbr}}$ | Jaccard medio por nodo de kNN |
| Graph stability | $J_{\text{edge}}$ | Jaccard global de aristas |
| Subspace stability | SS | Cosenos de ángulos principales |

### Protocolo de ejecución

1. **Baseline:** Pipeline completo (todos los marcadores, imputación mode, seed 42).
2. **Perturbación:** Aplicar el cambio (subsample / inject / impute).
3. **Reconstrucción:** PCA-30D → kNN (k=15, coseno) sobre datos perturbados.
4. **Comparación:** Calcular $J_{\text{nbr}}$, $J_{\text{edge}}$, SS entre baseline y perturbado.
5. **Seeds:** 3 seeds (42, 52, 62) para el subsampleo/injection random; reportar media ± std.
6. **Espacio de comparación:** Siempre PCA-30D (no UMAP).

### Qué NO se hace (explícitamente)

- No se reportan métricas en espacio UMAP (la capa estocástica se evalúa por separado, ADR-0014).
- No se prueba MNAR (reconocido como limitación en §7 del paper).
- No se varía $k$ ni la métrica de distancia dentro de las curvas de robustez.

## Estado actual

Las Tablas 6–8 del short paper ya reportan estos resultados con el protocolo descrito. Este ADR formaliza y congela el diseño para que sea referenciable.

## Consecuencias

- Las curvas de robustez son el corazón de OE2 y sostienen el claim principal del paper.
- El póster (bloque 3: "Robustness Under Perturbation") visualiza exactamente estos datos.
- Block removal (ADR-0013) extiende este protocolo sin modificarlo.
