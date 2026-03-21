# ADR-0013: Experimento de Remoción por Bloques vs. Submuestreo Aleatorio

**Estado:** Aceptada  
**Fecha:** 2026-03-06  
**Contexto:** ADR-0011 define perturbaciones de submuestreo aleatorio de marcadores. Pero un reviewer podría objetar que la robustez observada se debe a redundancia genómica: marcadores adyacentes son casi idénticos, así que eliminar uno al azar no destruye información. Se necesita un stress-test más agresivo.

## Problema

> ¿La topología del mapa resiste la eliminación de *bloques contiguos* de marcadores, o solo resiste la eliminación aleatoria (que preserva la cobertura genómica por redundancia)?

## Diseño experimental

### Variables

| Factor | Niveles |
|--------|---------|
| Tipo de eliminación | Aleatorio, Bloque contiguo |
| % eliminado | 5%, 10%, 20% |
| Posición del bloque | Inicio, Centro, Final del panel |
| Semilla (aleatorio) | 42, 52, 62 |

### Nota sobre "contiguos"

No disponemos de coordenadas cromosómicas. Usamos el **orden de columnas del archivo DArT** como proxy de posición genómica (los marcadores se reportan agrupados por cromosoma/scaffold).

### Protocolo

Para cada panel × % × tipo:

1. Partir de la matriz numérica ya imputada del panel.
2. Eliminar marcadores según el esquema.
3. PCA-30D sobre los marcadores restantes.
4. kNN (k=15, coseno).
5. Comparar contra baseline (todos los marcadores): $J_{\text{nbr}}$, $J_{\text{edge}}$, SS.

### Hipótesis esperada

- $J_{\text{edge}}^{\text{bloque}} < J_{\text{edge}}^{\text{aleatorio}}$ a igualdad de %.
- Pero si $J_{\text{edge}}^{\text{bloque}}$ sigue >0.85 al 10%, el claim de robustez se fortalece bajo el argumento más duro.

## Entregables esperados

- Tabla comparativa: $J_{\text{edge}}$ (aleatorio vs. bloque) por panel × %.
- Gráfico: barras pareadas (aleatorio/bloque) por nivel de eliminación.
- Frase de resultado para el póster (bloque 3): "Even contiguous block removal (10% of markers) preserves >X% of kNN edges."

## Lo que NO se hace

- No se simulan posiciones cromosómicas reales.
- No se hace GWAS ni linkage analysis.
- No se eliminan bloques >20% (sería mutilación, no robustez).

## Prioridad

**Alta** (Fase B). Es el complemento directo de ADR-0011 y la evidencia más fuerte contra el argumento de "robustez por redundancia".

## Relación con el póster

Resultado directo para bloque 3 (robustness curves). Si el resultado es limpio, una línea de texto adicional puede mencionarlo. No requiere nueva figura.

## Consecuencias

- Si bloques de 10% destruyen la topología, se documenta como limitación honesta (§6 del paper).
- Si resiste, es el claim más fuerte del trabajo.

## Estado actual

- Experimento ejecutado en modo completo para los 4 paneles, con 5%, 10% y 20% de remoción.
- Resultado observado: a 10% de remoción, los bloques contiguos preservan $J_{\text{edge}} > 0.85$ en todos los paneles y se mantienen cerca del submuestreo aleatorio.
- El hallazgo ya se incorporó al póster beta como `fig8_block_removal` y como evidencia adicional en robustez.
