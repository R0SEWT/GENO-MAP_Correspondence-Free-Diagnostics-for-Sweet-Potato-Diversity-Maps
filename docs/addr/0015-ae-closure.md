# ADR-0015: Cierre del Benchmark AE — No Más Variantes

**Estado:** Aceptada  
**Fecha:** 2026-03-06  
**Contexto:** Se entrenó un autoencoder (AE) con arquitectura 128–64–32–16–32–64–128, dropout 0.3, 50 épocas. Los resultados muestran que PCA ≥ AE en trustworthiness y continuity para todos los paneles. Se necesita cerrar formalmente el benchmark para evitar la tentación de seguir probando variantes.

## Problema

> El AE se diseñó como ablación ("¿una representación no lineal mejora el mapa?"). La respuesta experimental es NO para estos datos y esta tarea. ¿Se sigue explorando o se cierra?

## Decisión

**Se cierra.** El AE pasa de "candidato alternativo" a "ablación negativa documentada". No se probarán más arquitecturas, no se añadirá variational AE, no se probará transformer.

## Justificación

1. **PCA ≥ AE en todas las métricas canónicas** ($J_{\text{edge}}$, $J_{\text{nbr}}$, SS, Trustworthiness, Continuity) para los 4 paneles.
2. **Reproducibilidad**: PCA es determinista (salvo signo de eigenvectors, que no afecta distancias). AE depende de inicialización, learning rate schedule, batch size.
3. **Parsimonia**: el pipeline PCA→kNN→UMAP no tiene hiperparámetros de entrenamiento. El AE añade ≥6 hiperparámetros sin mejora.
4. **Scope**: el objetivo (ADR-0009) es evaluar robustez del mapa, no optimizar la representación.

## Entregables finales del benchmark AE

| Entregable | Estado |
|-----------|--------|
| Tabla comparativa PCA vs AE (Trustworthiness, $J_{\text{edge}}$, $J_{\text{nbr}}$, SS) por panel | ✅ Consolidada con métricas unificadas para cierre metodológico |
| Figura fig6 (PCA vs AE, escala 0–1) | ✅ En el póster (bloque 5) |
| Frase de resultado en short paper | ✅ Integrada en §4.6 |
| Nota de cierre en esta ADR | Este documento |

### Consolidación realizada

Con las métricas unificadas de ADR-0010:

1. Se consolidó la comparación PCA vs AE usando $J_{\text{edge}}$ como métrica primaria de reproducibilidad del grafo.
2. Se confirmó la brecha estable a favor de PCA en los cuatro paneles.
3. Se verificó consistencia con la figura comparativa del póster.
4. El benchmark queda cerrado como ablación negativa documentada.

## Lo que NO se hace

- No se prueba Variational AE (VAE).
- No se prueba transformer / attention-based encoder.
- No se hace grid search de arquitecturas AE.
- No se optimiza el AE para mejorar trustworthiness (eso sería cambiar el objetivo).
- No se añade el AE como rama permanente del pipeline.

## Relación con el póster

- **Bloque 5** ya tiene fig6 (PCA vs AE, escala unificada 0–1). No requiere cambios.
- **Bloque 6** (Stability Regime) muestra solo PCA. No se añade AE aquí.
- Si el recálculo con métricas unificadas confirma los valores actuales, cero cambios al póster.

## Prioridad

**Media** (Fase B). Resultado ya consolidado y suficiente para cerrar la rama AE dentro del scope del paper.

## Consecuencias

- El AE queda documentado como evidencia de que PCA es suficiente, no como limitación.
- Futuros trabajos pueden explorar representaciones no lineales, pero este paper no lo hará.
- La sección de limitaciones (§6) puede mencionar: "only one AE architecture was tested; deeper exploration of non-linear embeddings is left for future work."

## Estado actual

- El cierre del benchmark AE ya está reflejado en el short paper y en el póster beta.
- Resultado operativo consolidado: PCA mantiene una ventaja clara de estabilidad inter-seed frente a AE en los 4 paneles.
- No se abrirán nuevas variantes AE dentro de este proyecto.
