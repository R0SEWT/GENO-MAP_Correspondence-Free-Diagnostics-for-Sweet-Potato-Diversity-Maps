# ADR-0012: Sensibilidad al Número de PCs / Varianza Retenida

**Estado:** Aceptada  
**Fecha:** 2026-03-06  
**Contexto:** El pipeline usa PCA-30D por convención, reteniendo >85% de varianza en todos los paneles. Se necesita un experimento que responda: ¿cuánta varianza retenida es suficiente para preservar la topología útil del mapa?

## Problema

El número de PCs (30) se eligió pragmáticamente. No hay evidencia experimental de que 30 sea el punto óptimo, ni que el 85% de varianza sea el umbral correcto. La pregunta real es:

> ¿A partir de cuántos PCs (o qué % de varianza) la topología kNN ya no cambia significativamente?

## Diseño experimental

### Opción A — Número fijo de PCs

$$k \in \{5, 10, 15, 20, 30, 50\}$$

Para cada $k$:
1. PCA-$k$D sobre datos imputados.
2. kNN (k=15, coseno) sobre los primeros $k$ PCs.
3. Comparar contra baseline (PCA-30D): $J_{\text{nbr}}$, $J_{\text{edge}}$, SS.

### Opción B — Varianza acumulada

Determinar el número mínimo de componentes para alcanzar cada umbral:

$$\text{varianza acumulada} \in \{50\%, 65\%, 75\%, 85\%, 90\%\}$$

Y comparar contra baseline.

### Métrica clave

$$R(k) = J_{\text{edge}}(\text{kNN}_k, \text{kNN}_{\text{baseline}})$$

Donde $\text{kNN}_k$ es el grafo construido en PCA-$k$D y $\text{kNN}_{\text{baseline}}$ es el grafo en PCA-30D.

Complementar con SS entre subespacios.

### Paneles

Los 4 datasets: global_snp, global_silico, lowdensity_snp, lowdensity_silico.

### Seeds

3 seeds (42, 52, 62); mismos datos, variando solo la dimensionalidad PCA.

## Entregables esperados

- Tabla: $R(k)$ por panel × número de PCs.
- Curva: $J_{\text{edge}}$ vs. PCs (o vs. % varianza).
- Heurística: "a partir de $k = X$ (correspondiente a ~Y% varianza), la topología se estabiliza ($J_{\text{edge}} > 0.95$)".

## Prioridad

**Alta** (Fase B del plan operativo). Conecta directamente con el claim de robustez y puede simplificar el pipeline si se demuestra que menos de 30 PCs son suficientes.

## Relación con el póster

Si el resultado es interesante, podría mencionarse en el bloque 3 (robustness) como una línea adicional de evidencia. No requiere un bloque nuevo.

## Consecuencias

- Si $k=15$ es suficiente, el pipeline se simplifica y se documenta.
- Si el umbral es >30, se ajusta y se documenta por qué.
- El resultado alimenta §5.4 del short paper (discusión de robustez).

## Estado actual

- Experimento ejecutado en los 4 paneles con seeds 42, 52 y 62.
- Resultado observado: $J_{\text{edge}}$ alcanza su máximo en PCA-30D; 5--20 PCs pierden topología y 50 PCs introduce degradación por dimensiones ruidosas.
- El hallazgo ya se incorporó al póster beta como `fig7_pc_sensitivity`.
