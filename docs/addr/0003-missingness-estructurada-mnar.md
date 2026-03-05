# ADR-0003: Missingness estructurada — no es MCAR

**Estado:** Aceptada  
**Fecha:** 2026-03-01

## Contexto

El análisis exploratorio de missingness (`notebooks/eda_variables.ipynb`, Sección 4)
reveló que los datos faltantes en los genotipos DArT/DArTSeq **no son aleatorios**.

### Evidencia cuantitativa

| Dataset | Missing global | Mediana por marcador | Mediana por muestra | Muestras >50% |
|---------|---------------|---------------------|--------------------|--------------:|
| Global SNP | 6.5% | 4.5% | 5.7% | 0 |
| Global SilicoDArT | 2.3% | 2.2% | 2.1% | 0 |
| LowDensity SNP | **18.7%** | **15.6%** | **17.3%** | **4** |
| LowDensity SilicoDArT | 8.2% | 7.7% | 8.2% | 1 |

### Evidencia visual

Los heatmaps de missingness (celdas 10–11) muestran:

1. **Heatmap sin ordenar (LowDensity SNP):** Columnas verticales densas en muestras
   ~100–150, indicando un grupo de muestras con missingness coordinada.

2. **Heatmap ordenado por missingness:** Patrón de **"escalera"** en la esquina
   inferior-derecha — un subconjunto de marcadores (>40% de los subsampleados) falta
   consistentemente en un subconjunto de muestras. Esto es la firma clásica de
   **batch effect** o **placa de laboratorio**.

3. **Distribución bimodal:** La missingness por marcador en LowDensity SNP muestra
   un pico en ~0% y una cola larga hasta 50%, no la distribución unimodal esperada
   bajo MCAR.

### Outliers detectados (IQR method)

- **Global SNP:** 475 muestras outlier (>10.5% miss) de 5,970
- **LowDensity SNP:** 23 muestras outlier (>29.7% miss) de 635, 4 con >50%
- **LowDensity SilicoDArT:** 15 muestras outlier (>13.2% miss) de 635

## Decisión

1. **Clasificar la missingness como MAR/MNAR**, no MCAR. Los datos faltantes son
   informativos: correlacionan con batch de secuenciación, no con el genotipo real.

2. **No usar imputación simple** (media, moda, kNN sin ajuste por batch). Métodos
   que asumen MCAR producirían sesgos.

3. **Umbrales de filtrado propuestos** (a validar en análisis posteriores):
   - Descartar **muestras** con >50% missing (4 en LowDensity SNP)
   - Descartar **marcadores** con >50% missing
   - Considerar umbral más estricto (>20%) para análisis de estructura poblacional

4. **Usar el patrón de missingness como feature** en el pipeline:
   - El vector de missingness por muestra puede servir como covariable de batch
   - Considerar análisis de co-missingness para identificar placas de laboratorio
   - El modelo de imputación (si se usa) debe ser batch-aware

5. **Estrategia de imputación recomendada** (para trabajo futuro):
   - Masked genotype modeling (como en el roadmap `docs/roadmap/scale.md`)
   - Imputación condicional al batch inferido del patrón de missing
   - Nunca imputar >30% de los valores de un marcador

## Consecuencias

- **Pro:** Evita sesgos por imputación naive que distorsionaría estructura genética
- **Pro:** El patrón de missing aporta información sobre diseño experimental
- **Pro:** Los umbrales protegen contra marcadores/muestras de baja calidad
- **Con:** Reducción del dataset efectivo tras filtrado (~5–10% de muestras)
- **Con:** Mayor complejidad en el pipeline de preprocesamiento
- **Riesgo:** Los umbrales propuestos son heurísticos — deben validarse con
  análisis de sensibilidad (¿cambian los clusters al variar el umbral?)
