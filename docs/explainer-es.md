# TLDR

Este trabajo muestra que:

> Es posible validar rigurosamente mapas de diversidad genómica incluso cuando los datasets no comparten identificadores, y que en matrices ultra-anchas el embedding PCA es más estable que modelos profundos.


# La idea central de la investigación

Este trabajo intenta resolver un problema muy concreto en **genómica agrícola**.

Hoy en día podemos medir **decenas de miles de marcadores genéticos** para cada planta.
Eso genera matrices enormes donde cada fila es una variedad y cada columna un marcador del genoma.

Con esos datos queremos responder preguntas como:

* ¿Qué tan diversa es una colección de plantas?
* ¿Qué variedades son genéticamente similares?
* ¿Existen grupos o subpoblaciones genéticas?
* ¿Hay accesiones duplicadas en el banco de germoplasma?

Para explorarlo visualmente se suelen construir **mapas de diversidad genética**.

Pero aparece un problema práctico importante.

---

# El problema real

Los datos genómicos suelen provenir de **paneles de genotipado distintos**.

Cada panel puede usar **un sistema de identificadores diferente**.

Por ejemplo, un dataset puede identificar una planta como:

```
CIP_12345
```

y otro como:

```
Plate_7_Well_B03
```

Sin una tabla externa que conecte ambos sistemas de IDs, **no hay forma de saber si representan la misma planta**.

Esto rompe muchos métodos que intentan comparar datasets entre sí.

---

# La idea de GENO-MAP

En lugar de intentar alinear datasets con identificadores incompatibles, el proyecto propone otra estrategia:

> **evaluar cada panel de datos por separado sin requerir correspondencia entre muestras.**

A este enfoque lo llamamos:

**correspondence-free validation**.

La idea es analizar la **geometría del dataset** y responder preguntas como:

* ¿El panel induce una estructura genética real?
* ¿Los resultados dependen demasiado de decisiones del pipeline?
* ¿Las relaciones entre muestras son estables?

---

# Qué datos analizamos

Usamos datasets de **camote (*Ipomoea batatas*)** provenientes del banco de germoplasma del CIP.

Cada dataset tiene la forma:

```
muestras × marcadores
```

Por ejemplo:

| dataset           | muestras | marcadores |
| ----------------- | -------- | ---------- |
| Global SNP        | 5970     | 20069      |
| Global Silico     | 5970     | 57715      |
| LowDensity SNP    | 630      | 62732      |
| LowDensity Silico | 635      | 38272      |

Esto produce matrices **ultra-anchas**.

En algunos casos tenemos algo como:

```
630 muestras
62 000 marcadores
```

Es decir:

[
n \ll p
]

Este régimen es complicado porque muchos métodos de machine learning se vuelven **inestables** en este escenario.

---

# Qué hace el pipeline GENO-MAP

El pipeline transforma estas matrices gigantes en algo que los genetistas puedan explorar.

Produce tres cosas principales:

1. **mapas 2D de diversidad genética**
2. **grafos de similitud entre plantas**
3. **diagnósticos automáticos de calidad del dataset**

El flujo del método es:

```
genotype matrix
      ↓
imputation (missing values)
      ↓
PCA (espacio analítico)
      ↓
UMAP (visualización)
      ↓
k-NN graph
```

Un detalle importante:

* **PCA se usa para calcular relaciones reales entre muestras**
* **UMAP se usa solo para visualizar**

Esto es importante porque UMAP es estocástico.

---

# Qué es exactamente el mapa de diversidad

Imagina que cada planta es un punto en un espacio de **20 000 dimensiones**.

Cada dimensión corresponde a un marcador genético.

Eso es imposible de visualizar directamente.

Entonces hacemos una reducción de dimensionalidad:

```
20000 dimensiones
↓
PCA
↓
30 dimensiones
↓
UMAP
↓
2 dimensiones
```

En el mapa final:

* puntos cercanos → plantas genéticamente similares
* puntos lejanos → plantas genéticamente distintas

---

# Cómo evaluamos si el mapa es confiable

El paper introduce dos herramientas.

---

## 1. Diagnósticos geométricos del panel

Se analizan propiedades del embedding, por ejemplo:

* dimensionalidad efectiva
* dominancia de PC1
* reciprocidad de vecinos
* componentes desconectadas

Esto permite detectar problemas como:

```
HIGH-MISSINGNESS
EXTREME-WIDE
DISCONNECTED
```

que indican datasets potencialmente problemáticos.

---

## 2. Curvas de robustez

Aquí viene una de las ideas más interesantes.

Degradamos artificialmente el dataset para ver qué tan estable es la estructura.

Por ejemplo:

* eliminar 95 % de los marcadores
* agregar datos faltantes
* cambiar el método de imputación

Luego medimos cuánto cambia el grafo de similitud.

Las métricas principales son:

**Jaccard de vecinos**

[
J = \frac{|N_i \cap N_i'|}{|N_i \cup N_i'|}
]

y **similaridad del subespacio PCA**.

---

# Resultado importante

Incluso cuando eliminamos **95 % de los marcadores**, el subespacio PCA casi no cambia:

[
SS \ge 0.91
]

Esto significa que la estructura genética global **no depende de marcadores específicos**.

La señal está distribuida por todo el genoma.

En otras palabras:

> **la señal genómica es redundante y robusta.**

---

# Segundo resultado importante

Los vecindarios del grafo se degradan **de forma gradual**.

No ocurre un colapso repentino.

Por ejemplo:

```
5% markers   → J ≈ 0.43
20% markers  → J ≈ 0.61
50% markers  → J ≈ 0.75
80% markers  → J ≈ 0.84
```

Esto permite definir **umbrales de calidad** para datasets genómicos.

---

# El experimento con autoencoders

También probamos un enfoque moderno:
un **autoencoder neuronal** para aprender el embedding.

La idea era que un modelo no lineal podría capturar mejor la estructura genética.

Pero ocurrió algo interesante.

---

# Resultado inesperado

El autoencoder mejora ligeramente una métrica llamada **trustworthiness**, pero pierde mucha estabilidad.

Por ejemplo:

```
PCA stability ≈ 0.89
AE stability  ≈ 0.52
```

En datasets pequeños el problema es aún peor.

Esto significa que:

* cada entrenamiento produce un embedding distinto
* el grafo cambia entre ejecuciones

Por eso la conclusión es:

> En matrices genómicas **n ≪ p**, PCA sigue siendo el método más confiable.

---

# El insight conceptual del trabajo

El aporte principal no es solo el pipeline.

Es el **marco de validación**.

El trabajo muestra que es posible evaluar mapas de diversidad genética usando:

* diagnósticos geométricos
* estabilidad bajo perturbaciones

sin necesidad de alinear datasets con identificadores incompatibles.

---
