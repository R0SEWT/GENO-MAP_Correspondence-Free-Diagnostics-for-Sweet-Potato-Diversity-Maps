# TL;DR

Este trabajo muestra que:

> Es posible validar rigurosamente mapas de diversidad genómica incluso cuando los datasets no comparten identificadores, analizando directamente la **geometría del dataset** y su **estabilidad bajo perturbaciones**.

Además encontramos que en matrices genómicas **ultra-anchas** ((n \ll p)), el embedding **PCA es sustancialmente más estable que modelos profundos**.

---

# El problema

En genómica agrícola moderna es común trabajar con matrices donde:

```
filas   → muestras (plantas)
columnas → marcadores genéticos
```

Un dataset típico puede verse así:

```
630 muestras
62 000 marcadores
```

Esto produce el régimen estadístico:

$$
n \ll p
$$

donde el número de variables supera ampliamente al número de muestras.

Con estos datos queremos estudiar **diversidad genética**, por ejemplo:

* identificar subpoblaciones
* detectar duplicados en bancos de germoplasma
* explorar relaciones genéticas entre accesiones

Para hacerlo, los genetistas suelen construir **mapas de diversidad genética** mediante reducción de dimensionalidad.

---

# El obstáculo práctico

En la práctica los datasets provienen de **paneles de genotipado distintos**.

Cada panel puede usar un sistema diferente de identificadores de muestras.

Por ejemplo:

dataset A

```
CIP_12345
```

dataset B

```
Plate_7_Well_B03
```

Sin una tabla externa que conecte ambos sistemas de IDs, **no sabemos si representan la misma planta**.

Esto rompe muchos enfoques de validación que requieren comparar embeddings entre datasets.

---

# La idea central de GENO-MAP

El proyecto propone cambiar completamente la perspectiva.

En lugar de intentar alinear datasets con identificadores incompatibles, evaluamos **cada panel de manera independiente**.

La pregunta pasa a ser:

> ¿La estructura geométrica inducida por el dataset es consistente y estable?

A este enfoque lo llamamos:

**correspondence-free validation**

La validación no depende de saber qué muestra corresponde a cuál en otro dataset.

Depende únicamente de propiedades **intrínsecas de la geometría del embedding**.

---

# Geometría del dataset

Un dataset genómico puede interpretarse como una nube de puntos en un espacio de alta dimensión.

Si tenemos 20 000 marcadores, cada planta vive en:

$$
\mathbb{R}^{20000}
$$

El objetivo de la reducción de dimensionalidad es encontrar una representación más compacta que preserve relaciones genéticas.

En GENO-MAP usamos una jerarquía:

```
genotype matrix
      ↓
imputation
      ↓
PCA (espacio analítico)
      ↓
UMAP (visualización)
      ↓
kNN graph
```

Aquí hay una decisión importante:

* **PCA define el espacio analítico**
* **UMAP solo se usa para visualización**

Todas las relaciones entre muestras se calculan en el espacio PCA.

---

# Qué produce el pipeline

El sistema genera tres artefactos principales:

### 1. Mapas 2D de diversidad

Representaciones visuales donde:

* puntos cercanos → plantas genéticamente similares
* puntos lejanos → plantas genéticamente distintas

---

### 2. Grafos de similitud genética

Se construye un grafo **k-nearest neighbors** donde cada nodo es una planta y las aristas conectan vecinos genéticos.

Esto permite estudiar estructura poblacional o detectar duplicados.

---

### 3. Diagnósticos automáticos del dataset

El pipeline calcula métricas que describen la **geometría del embedding**.

---

# Diagnósticos geométricos

Se analizan varias propiedades estructurales del espacio PCA.

---

## Dimensionalidad efectiva

Mide cuántas dimensiones contienen señal genética real.

Se estima a partir del espectro de autovalores del PCA.

Un rank efectivo bajo indica que la diversidad genética está concentrada en pocos ejes.

---

## Dominancia del primer componente

Se mide la fracción de varianza explicada por PC1.

Valores muy altos pueden indicar:

* estructura poblacional fuerte
* o artefactos técnicos.

---

## Reciprocidad de vecinos

Se construye un grafo k-NN y se mide qué proporción de relaciones son recíprocas.

Alta reciprocidad indica que la geometría del embedding es coherente.

Reciprocidad baja puede indicar ruido o distancias poco informativas.

---

## Componentes desconectadas

Se analiza el número de componentes conectadas del grafo.

Muchas componentes pueden indicar:

* datasets fragmentados
* missingness extremo
* errores de genotipado.

---

A partir de estos diagnósticos el sistema genera **flags automáticos** como:

```
HIGH-MISSINGNESS
EXTREME-WIDE
DISCONNECTED
```

que alertan sobre posibles problemas en el dataset.

---

# Validación mediante perturbaciones

La segunda parte del framework evalúa **estabilidad estructural**.

La idea es simple:

> si la estructura genética es real, debería sobrevivir perturbaciones razonables del dataset.

Introducimos perturbaciones controladas como:

* subsampling de marcadores
* inyección de missing values
* cambios en el método de imputación.

Luego medimos cuánto cambia la estructura del embedding.

---

# Métricas de estabilidad

Dos métricas principales.

---

## Estabilidad del subespacio PCA

Se mide la similitud entre subespacios PCA antes y después de la perturbación.

Valores cercanos a 1 indican que la estructura global del dataset se conserva.

---

## Estabilidad de vecindarios

Para cada muestra se comparan los vecinos antes y después de la perturbación usando índice de Jaccard:

$$
J = \frac{|N_i \cap N_i'|}{|N_i \cup N_i'|}
$$

Esto mide si las relaciones locales entre plantas se mantienen.

---

# Resultado importante

Incluso cuando eliminamos **95 % de los marcadores**, el subespacio PCA se mantiene muy estable:

$$
SS \ge 0.91
$$

Esto sugiere que la señal genética está **distribuida a lo largo del genoma**, en lugar de depender de marcadores específicos.

---

# Estabilidad del grafo

Los vecindarios del grafo se degradan **de forma gradual**.

Por ejemplo:

```
5% markers   → J ≈ 0.43
20% markers  → J ≈ 0.61
50% markers  → J ≈ 0.75
80% markers  → J ≈ 0.84
```

No se observa un colapso abrupto de la estructura.

Esto permite definir **umbrales de calidad para datasets genómicos**.

---

# Experimento con autoencoders

También probamos embeddings aprendidos mediante **autoencoders neuronales**.

La hipótesis era que modelos no lineales podrían capturar mejor la estructura genética.

---

# Resultado

Los autoencoders mejoran ligeramente una métrica llamada **trustworthiness**, pero pierden estabilidad.

Por ejemplo:

```
PCA stability ≈ 0.89
AE stability  ≈ 0.52
```

Además cada entrenamiento produce embeddings diferentes.

Esto genera grafos de similitud inconsistentes.

---

# Interpretación

La razón está en el régimen estadístico del problema.

En datasets donde:

$$
n \ll p
$$

PCA tiene ventajas estructurales:

* solución analítica
* embedding determinístico
* comportamiento reproducible.

Los autoencoders, en cambio, dependen de optimización no convexa y pueden converger a representaciones distintas.

---

# Insight conceptual

El aporte principal del trabajo no es solo el pipeline.

Es el **marco de validación**.

Mostramos que es posible evaluar mapas de diversidad genética usando:

* diagnósticos geométricos
* estabilidad bajo perturbaciones

sin necesidad de alinear datasets con identificadores incompatibles.

En otras palabras:

> La validez de un mapa de diversidad genética depende de que la **geometría del dataset sea estable**, no de conocer la correspondencia entre muestras.
