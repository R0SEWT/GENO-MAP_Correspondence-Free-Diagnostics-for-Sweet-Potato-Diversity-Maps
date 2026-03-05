## 1) Qué NO puedes afirmar / hacer (y por qué)

Con 0 IDs en común, no puedes:

* comparar embeddings “punto a punto”
* medir Jaccard de vecinos entre paneles para los mismos nodos
* decir “este cluster en LowDensity corresponde a este cluster en Global”

Eso sería *hallucinated alignment*.

---

## 2) Qué SÍ puedes hacer: “comparación sin correspondencia” (invariantes)

La idea es: aunque no haya matching, sí puedes comparar **propiedades del espacio**, no de individuos.

### A. Invariantes espectrales / geométricos (muy fuertes)

En PCA-space (y también en distancias):

* **effective rank / participation ratio**
* **espectro de autovalores** (cómo cae la varianza)
* **PC1 dominance** (ya lo viste en lowdensity_silico)
* **curva de distancias kNN** (densidad local)
* **conectividad del grafo**: #componentes, giant component ratio
* **distribución de grados** (aunque kNN es casi regular, hay efectos por ties/reciprocidad)
* **reciprocal kNN rate** (proporción de aristas mutuas)

Estas métricas responden:

> “¿El dataset/panel induce una geometría rica o degenerada?”

Y eso es curaduría + QA real.

### B. Robustez bajo intervenciones (tu mejor arma)

Sin IDs compartidos, tu validación “comparativa” debe ser:

* subsampling de marcadores (5/10/20/50%)
* inyección controlada de missing (+5/+10/+20%)
* variación de imputación (modo/mediana/AE-imputer)
* variación de k (15/20/30)
* seeds UMAP (solo visual)

Comparas **curvas**, no instancias.

Ejemplo de claim defendible:

> “LowDensity panels are more sensitive to missingness and marker subsampling than Global panels, as measured by topology drift and effective rank.”

Eso es ciencia útil y honesta.

### C. “Overlap biológico plausible” pasa a ser *hipótesis*, no análisis

Lo dices tal cual lo pusiste:

* plausible
* no verificable sin manifest
* deliverable estándar pendiente

Y listo.

---

## 3) Roadmap orientado correctamente (con namespaces disjuntos)

### Módulo 1 — Panel QA & Geometry Diagnostics (obligatorio)

**Objetivo:** caracterizar y detectar degeneraciones por panel.

**Outputs (paper/poster)**

* tabla por dataset: effective rank, PC1/PC2 ratio, kNN distance stats, component ratio
* flags automáticos: “1D-dominant”, “collapsed neighborhood”, “high missingness sensitivity”

Esto convierte tu proyecto en **herramienta de control de calidad**.

---

### Módulo 2 — Robustness Curves (tu “comparación entre paneles” real)

**Objetivo:** medir sensibilidad estructural a decisiones razonables.

**Experimentos**

* markers subsampling: 5/10/20/50/80%
* missing injection: +0/+5/+10/+20% (MCAR)
* imputación: mode vs median vs learned-imputer (si lo metes)

**Métricas**

* Jaccard vecinos/aristas (intra-dataset)
* drift de distribución de distancias kNN
* estabilidad del subespacio PCA (principal angles)

**Outputs**

* 2–3 figuras tipo “curvas” por dataset/panel

Esto reemplaza “alignment por ID” de manera sólida.

---

### Módulo 3 — Decision Support (core-set selection) por panel

**Objetivo:** producir subsets representativos *dentro de cada panel*.

Aunque no puedas cruzar Global↔LowDensity, sí puedes entregar:

* “core set” para Global
* “core set” para LowDensity
* “diversity coverage vs k”

Eso es aplicable de inmediato en CIP, sin mapping.

---

### Módulo 4 — Manifest acquisition as an explicit deliverable (future work con valor)

Esto debe aparecer como:

* **Limitación estructural**: IDs disjuntos
* **Acción concreta**: obtener manifest DArT (plate+well → CIP accession)
* **Resultado esperado**: permitir cross-panel matching, batch correction, transfer alignment


---

## 4) Cómo lo vendes sin que parezca “nos faltó data”

Clave: no lo vendas como carencia, sino como **diseño robusto ante falta de correspondencia**.

> *Because genotyping panels use disjoint identifier namespaces (CIP accession IDs vs DArT plate/well codes) and no manifest was available, we focus on correspondence-free validation: geometry diagnostics and robustness curves under controlled perturbations.*

]

## 5) Bonus: si quieres una “aproximación de alignment” sin manifest (con advertencia)

Solo si necesitas un *hint* exploratorio (no claim), podrías hacer:

* comparar distribuciones de metadatos (si LowDensity tiene alguno)
* buscar coincidencias por atributos no-ID (origen, institución, país, etc.)

Pero dado que LowDensity parece no tener mapeo a CIP accession, **no vale la pena**; mejor mantenerlo limpio.

---

## 6) Idas de implementacion evaluar viabilidad

1. Implementar `panel_diagnostics.py`:

   * effective rank
   * PC1 dominance
   * kNN distance stats
   * reciprocal kNN rate
   * component ratio

2. Implementar `robustness_curves.py`:

   * marker subsampling + missing injection

3. Escribir en el short paper:

   * “No manifest, no alignment; we do correspondence-free validation.”

