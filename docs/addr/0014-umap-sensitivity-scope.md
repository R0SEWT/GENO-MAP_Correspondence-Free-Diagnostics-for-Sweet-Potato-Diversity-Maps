# ADR-0014: Alcance del Análisis de Sensibilidad UMAP

**Estado:** Aceptada  
**Fecha:** 2026-03-06  

## Contexto

En el pipeline GENO-MAP, UMAP se utiliza únicamente como **capa de visualización bidimensional** del espacio analítico definido por PCA-30D y el grafo kNN correspondiente.

El análisis estructural del pipeline (diagnósticos geométricos, estabilidad topológica y curvas de robustez) se realiza exclusivamente en el **espacio PCA-30D**, que constituye el espacio analítico del método.

Sin embargo, un reviewer podría cuestionar si el layout UMAP introduce artefactos visuales o dependencias críticas de hiperparámetros que puedan afectar la interpretación del mapa.

Por tanto, se considera necesario demostrar explícitamente que:

- los diagnósticos estructurales del pipeline no dependen de UMAP, y
- UMAP actúa únicamente como **rendering perceptual** del grafo analítico.

## Problema

> ¿Cambios razonables en los hiperparámetros de UMAP alteran la interpretación estructural del mapa, o únicamente su apariencia visual?

## Posición

UMAP se considera **una capa de rendering perceptual**, no una representación analítica.

Por tanto:

- el grafo kNN se construye **exclusivamente en PCA-30D**,
- los diagnósticos estructurales se calculan **antes de UMAP**,
- UMAP no participa en ninguna decisión analítica del pipeline.

El objetivo del análisis de sensibilidad no es encontrar la “mejor” configuración de UMAP, sino demostrar que **las conclusiones analíticas del pipeline son invariantes al layout visual**.

## Diseño experimental

Se adopta un análisis en dos capas:

1. **Validación ejecutada**: barrido de configuraciones representativas para demostrar invariancia analítica.
2. **Comunicación del resultado**: selección de tres layouts representativos para mostrar que el mapa puede cambiar perceptualmente sin alterar la capa analítica.

### 1. Validación ejecutada

Se ejecutó el siguiente grid sobre los 4 paneles:

| Parámetro | Valores |
|-----------|---------|
| `n_neighbors` | 10, 15, 30 |
| `min_dist` | 0.1, 0.3, 0.5 |
| `random_state` | 42, 52, 62 |
| `metric` | cosine |

Total: 27 configuraciones por panel.

Este grid no se interpreta como búsqueda de hiperparámetros, sino como **stress test metodológico** para verificar que la capa analítica permanece invariante.

### 2. Layouts representativos para comunicación

De ese espacio de configuraciones, se pueden seleccionar tres layouts visualmente distintos para mostrar el argumento de forma simple.

#### Baseline layout

Configuración estándar utilizada en el pipeline.

| Parámetro | Valor |
|-----------|-------|
| `n_neighbors` | 15 |
| `min_dist` | 0.3 |
| `metric` | cosine |
| `random_state` | 42 |

#### Local-stress layout

Favorece estructuras locales más compactas y mayor separación visual.

| Parámetro | Valor |
|-----------|-------|
| `n_neighbors` | 10 |
| `min_dist` | 0.1 |
| `metric` | cosine |
| `random_state` | 52 |

#### Global-stress layout

Favorece una estructura más global y suavizada.

| Parámetro | Valor |
|-----------|-------|
| `n_neighbors` | 30 |
| `min_dist` | 0.5 |
| `metric` | cosine |
| `random_state` | 62 |

Estas configuraciones inducen **layouts perceptualmente distintos**, permitiendo evaluar si la interpretación del mapa depende de la proyección visual.

## Métricas evaluadas

### 1. Invariancia analítica (criterio principal)

Se verifica que el grafo analítico permanece invariante:

$$
J_{\text{edge}}^{\text{PCA}} = \text{constante}
$$

donde el grafo kNN se construye en PCA-30D.

Esto confirma que:

- los diagnósticos QA,
- las curvas de robustez,
- y la estructura topológica

no dependen del layout UMAP.

### 2. Divergencia visual del layout

Se comparan visualmente layouts UMAP representativos para confirmar que:

- el mapa puede cambiar perceptualmente,
- pero la estructura analítica permanece intacta.

Las visualizaciones pueden presentarse como:

```text
Baseline | Local-stress | Global-stress
```

con los mismos datos y escala visual.

### 3. Comparación opcional de grafos en UMAP-2D

Como análisis secundario, se calcula:

$$
J_{\text{edge}}(\text{kNN}_{\text{PCA}}, \text{kNN}_{\text{UMAP}})
$$

Se espera que el solapamiento sea moderado o bajo debido a la reducción de dimensionalidad, lo que refuerza que **UMAP-2D no debe utilizarse como espacio analítico**.

## Lo que NO se hace

Para evitar confusión metodológica, se establece explícitamente que:

- no se optimizan hiperparámetros de UMAP,
- no se selecciona un “mejor UMAP”,
- no se utiliza UMAP para construir grafos ni clustering,
- no se comparan layouts mediante métricas geométricas (Procrustes, etc.),
- UMAP no se integra en ninguna etapa de decisión analítica.

## Prioridad

**Media-baja**

Este análisis sirve principalmente como defensa metodológica frente a posibles críticas de reviewers, ya que el resultado esperado es que el grafo analítico en PCA-30D permanezca invariante.

## Relación con el póster

Este experimento no introduce un claim nuevo.

Sirve para reforzar la separación conceptual entre:

```text
análisis estructural (PCA + kNN)
visualización perceptual (UMAP)
```

En el póster puede mostrarse como evidencia de soporte o resumirse en una frase metodológica; no requiere convertirse en figura principal.

## Consecuencias

Si el experimento confirma la hipótesis (resultado esperado):

- se demuestra que UMAP es una **capa de rendering independiente del análisis**,
- se refuerza el argumento metodológico del pipeline,
- se evita la interpretación errónea de UMAP como espacio analítico.

En caso de que layouts extremos produzcan interpretaciones visuales ambiguas, se documentará como **artefacto perceptual** y se recomendará complementar el análisis con visualizaciones alternativas (p. ej., heatmaps de distancia).

## Resultado observado

La estructura analítica del pipeline:

```text
PCA-30D → kNN graph
```

permanece invariante bajo cambios en los hiperparámetros de UMAP. En la ejecución completa sobre 27 configuraciones por panel, el grafo analítico en PCA-30D mantuvo $J_{\text{edge}} = 1.0$, mientras que el kNN inducido en UMAP-2D compartió solo una fracción de las aristas del espacio analítico.

Conclusión: **GENO-MAP separa correctamente análisis estructural y visualización perceptual**.
