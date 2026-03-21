# ADR-0010: Unificación de Métricas de Estabilidad

**Estado:** Aceptada  
**Fecha:** 2026-03-06  
**Contexto:** Las Tablas 3 y 10 del short paper reportan estabilidad con métricas similares pero no idénticas, generando ambigüedad. Se necesita una definición operativa única para todo el paper y el póster.

## Problema

Actualmente se mezclan:
- "Jaccard vecinos (PCA)" — Jaccard medio por nodo de conjuntos kNN.
- "Jaccard aristas (PCA)" — Jaccard de conjuntos globales de aristas.
- "Edge Jaccard (seed-to-seed)" — usado en la tabla AE sin clarificar si es lo mismo.

El lector no puede comparar directamente Tabla 3 (baseline) con Tabla 10 (AE) porque no queda claro si se usa el mismo cálculo.

## Decisión

Se adoptan **tres métricas de estabilidad**, cada una con definición operativa explícita y nombre canónico:

### 1. Graph Stability ($J_{\text{edge}}$)

$$J_{\text{edge}}(G_1, G_2) = \frac{|E_1 \cap E_2|}{|E_1 \cup E_2|}$$

Jaccard sobre el conjunto global de aristas del grafo kNN construido en el **espacio analítico** (PCA-30D). Se computan entre pares de seeds; se reporta la media de los $\binom{3}{2} = 3$ pares.

**Uso:** Medida primaria de reproducibilidad del grafo. Aparece en Tablas 3 y 10, y en el póster (bloque 5, 6).

### 2. Neighbor Stability ($J_{\text{nbr}}$)

$$J_{\text{nbr}} = \frac{1}{n} \sum_{i=1}^{n} \frac{|N_k^{(1)}(i) \cap N_k^{(2)}(i)|}{|N_k^{(1)}(i) \cup N_k^{(2)}(i)|}$$

Jaccard medio por nodo de los conjuntos de $k$ vecinos más cercanos. Más sensible a cambios locales que $J_{\text{edge}}$.

**Uso:** Curvas de robustez (Tablas 6–8), donde mide sensibilidad local de vecindarios bajo perturbaciones.

### 3. Subspace Stability (SS)

$$\text{SS} = \frac{1}{d} \sum_{j=1}^{d} \cos(\theta_j)$$

Promedio de cosenos de los ángulos principales entre dos subespacios PCA de $d=30$ dimensiones, calculados vía descomposición QR → SVD de las matrices de scores.

**Uso:** Mide si la geometría global del subespacio PCA cambia. Aparece en curvas de robustez y en el póster (bloque 3, inset SS).

## Protocolo de cálculo

- **Seeds:** 42, 52, 62 (fijas para todo el proyecto).
- **k:** 15 (fijo).
- **Espacio:** PCA-30D con métrica coseno.
- **Pares:** $\{(42,52), (42,62), (52,62)\}$; se reporta media ± std.
- **Función canonical:** `compute_stability(knn_1, knn_2) → (J_edge, J_nbr)` y `subspace_similarity(scores_1, scores_2) → SS`.

## Consecuencias

- Todas las tablas del short paper se regeneran con estas tres métricas.
- El póster ya usa $J_{\text{edge}}$ y SS de manera consistente (no requiere cambios).
- Cualquier comparación futura (AE, UMAP, block removal) usa estas mismas funciones.
- Se elimina la ambigüedad entre "Jaccard vecinos" y "Edge Jaccard".
