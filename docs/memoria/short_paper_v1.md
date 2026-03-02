# GENO-MAP: Visual Mapping of Genomic Diversity in Germplasm Collections

## Abstract

Exploring genomic diversity in large germplasm collections is challenged by the high dimensionality of genotyping data. In this poster we present GENO-MAP, a reproducible approach to visually map and explore DArT and DArTSeq matrices using two-dimensional representations and k-nearest neighbor (kNN) similarity graphs. GENO-MAP combines explicit imputation, hierarchical dimensionality reduction (PCA → UMAP), and kNN graph construction on deterministic spaces, producing structured outputs designed for interactive visualization. We additionally explore a denoising autoencoder as a nonlinear alternative to PCA for the reduction step, finding that the learned bottleneck improves local structure preservation on large datasets ($n \approx 6\,000$) but degrades on small ones ($n < 700$) due to overparameterisation. The approach is illustrated on large-scale sweet potato (*Ipomoea batatas*) datasets, showing stable behavior under reasonable variations of analysis parameters. Rather than inferring causal biological relationships, GENO-MAP aims to facilitate inspection, comparison, and curation of genomic collections, supporting human decision-making in AI-assisted agricultural research. 

**Keywords:** genomic diversity, exploratory analysis, kNN graphs, PCA, UMAP, autoencoder, germplasm.

---

## 1. Introducción 

Los bancos de germoplasma de cultivos andinos albergan miles de accesiones genotipadas con tecnologías DArT y DArTSeq, generando matrices anchas (muestras × marcadores) con decenas de miles de loci. Explorar visualmente la estructura de diversidad en estas matrices es un paso esencial para (i) identificar grupos genéticos, (ii) detectar redundancias o brechas en la colección, y (iii) priorizar combinaciones de cruza que maximicen la diversidad genética. Sin embargo, las dimensiones de estos datos (> 50 000 marcadores) dificultan su visualización directa. Los métodos de reducción de dimensionalidad como UMAP y PCA permiten proyectar estos datos a 2D preservando la estructura local, mientras que los grafos kNN capturan relaciones de similitud entre muestras. En este trabajo presentamos GENO-MAP, un pipeline de código abierto que automatiza la cadena completa: lectura de matrices DArT/DArTSeq → imputación → reducción PCA+UMAP → grafo kNN → artefactos JSON para visualización, con un marco de validación integrado que evalúa calidad del embedding, estabilidad, sensibilidad a la imputación e interpretabilidad de marcadores. --- ## 2. Metodología ### 2.1 Datos Se utilizaron tres fuentes de datos de camote (*Ipomoea batatas*) provenientes del repositorio CIP Dataverse, con cinco datasets derivados: 

| Dataset | Fuente DOI | Tipo | Muestras | Marcadores | Missing rate | 
|---------|-----------|------|----------|------------|-------------| 
| global_snp | 10.21223/P3/S2IMOS | SNP | 5 970 | 20 101 | 5.52% | | global_silico | 10.21223/P3/S2IMOS | SilicoDArT | 5 970 | 57 715 | 3.21% | | lowdensity_snp | 10.21223/P3/UBDJ44 | SNP | 635 | 62 736 | 18.17% | | lowdensity_silico | 10.21223/P3/UBDJ44 | SilicoDArT | 649 | 38 272 | 4.85% | | wild_snp | 10.21223/P3/3VYY8C | SNP (parientes silvestres) | 635 | 62 736 | 18.17% | 


Los formatos de entrada se detectan automáticamente: sample_columns (primera columna = etiqueta de fila, columnas = muestras) y marker_metrics (columnas iniciales de metadatos, seguidas de columnas numéricas de muestra). ### 2.2 Pipeline de procesamiento El pipeline se implementa en Python 3.10+ con dependencias mínimas: NumPy, pandas, scikit-learn, umap-learn y matplotlib. El flujo consta de cuatro etapas: 
1. **Lectura y parseo**: Detección automática de separador (, / ;) y formato. Se separan filas de metadatos (no numéricas) de filas de marcadores. La matriz se transpone a la orientación muestras × marcadores. 
2. **Imputación**: Los valores faltantes (codificados como -, vacío, o NaN) se imputan mediante la estrategia de moda (most_frequent) por defecto con sklearn.impute.SimpleImputer. Se evalúa alternativamente la mediana como sensibilidad. 
3. **Reducción de dimensionalidad**: - **PCA**: Se calcula con $\min(50, p, n)$ componentes para obtener una representación intermedia. - **UMAP**: Se aplica sobre las componentes PCA con 2 dimensiones de salida ($n_{\text{components}} = 2$). Si UMAP no está disponible, se retienen las primeras 2 PCs. 
4. **Grafo kNN**: Se construye sobre las primeras 30 componentes PCA (o el embedding UMAP, según configuración). Se usa sklearn.neighbors.NearestNeighbors con $k \in \{15, 20\}$ y métrica coseno (SNP) o Jaccard (SilicoDArT). Cada arista se almacena como (source, target, distance). 

### 2.3 Salidas 

Para cada dataset se generan tres archivos JSON: - *_nodes.json: Lista de nodos con id, embedding: [x, y] y meta (metadatos asociados). - *_edges.json: Lista de aristas kNN con source, target y distance. - *_stats.json: Estadísticas del dataset (muestras, marcadores, missing rate, métrica, vecinos). ### 2.4 Marco de validación El módulo de validación (validate_embeddings.py) evalúa cuatro dimensiones: 

1. **Calidad del embedding**: Trustworthiness (Venna & Kaski, 2006) que mide la preservación de vecindarios locales entre el espacio original y el embedding 2D.

 2. **Estabilidad**: Se ejecutan múltiples semillas (42, 52, 62) y se compara: - Solapamiento de vecinos por nodo (Jaccard medio) en el espacio de grafo. - Solapamiento global de aristas (Jaccard del conjunto de aristas). - Solapamiento de vecinos en el espacio del embedding UMAP. 
 
 3. **Sensibilidad a la imputación**: Se compara la topología del grafo entre imputación por moda y por mediana, midiendo Jaccard de vecindarios y aristas respecto a la línea base. 
 
 4. **Interpretabilidad de marcadores**: Se reportan las cargas (loadings) PCA de los 10 marcadores con mayor contribución a PC1 y PC2, junto con la varianza explicada acumulada. 

### 2.5 Autoencoder como alternativa a PCA (Level 1)

Como alternativa a la reducción lineal PCA, se implementa un autoencoder denoising (DAE) que aprende una representación no lineal comprimida.

**Arquitectura** (`GenoAutoencoder`):
- Encoder: Dense(input → 512) → ResBlock × 2 → Dense(512 → 64)
- Decoder: Dense(64 → 512) → ResBlock × 2 → Dense(512 → input)
- Cada ResBlock: LayerNorm → Dense → GELU → Dropout → Dense + skip connection
- Bottleneck de 64 dimensiones

**Entrenamiento**:
- Denoising: 15% de inputs enmascarados aleatoriamente
- Pérdida: MSE reconstrucción + penalización de dispersión ($\lambda$ sparsity sobre activación del bottleneck)
- Optimizador: AdamW, scheduler ReduceLROnPlateau
- Early stopping con paciencia configurable
- 3 seeds (42, 52, 62) para evaluación de estabilidad

**Variantes evaluadas** (tag `ae-v2`):
- Regularización agresiva para datasets pequeños: dropout 0.45, noise mask 0.30, weight decay 5×10⁻⁴
- Ensemble: promedio de vectores bottleneck entre seeds, UMAP sobre el bottleneck promediado
- Transfer learning: pre-entrenamiento en dataset grande (Global, $n=5\,970$) → fine-tuning en dataset pequeño (LowDensity, $n=630$) con freeze de encoder por 10 epochs

El grafo kNN y el embedding UMAP se construyen sobre el bottleneck de 64D en lugar de las 30 primeras PCs.

 
  ## 3. Configuración experimental 
  ### 3.1 Corridas de producción
  
   Se ejecutó el perfil full (sin límite de muestras ni marcadores, limit_rows=0) sobre los cinco datasets, con semilla fija (seed=42), paralelización de kNN (knn_jobs=-1) y generación de figuras PNG. El tag de referencia es poster-260206-2135. 
   
   ### 3.2 Validación 
   
   La validación se ejecutó con tres semillas (42, 52, 62), espacio de grafo PCA (graph-space=pca), campo de metadato Institution, e imputaciones comparadas (most_frequent, median). Se computaron 3 pares de semillas para estabilidad. ### 3.3 Corridas autoencoder

Se ejecutaron dos rondas de entrenamiento sobre cuatro datasets (excluyendo wild_snp por ser duplicado de lowdensity_snp):
- **ae-v1**: 200 epochs, patience 20, configuración base (dropout 0.20, noise 0.15).
- **ae-v2**: 50 epochs, patience 15, regularización agresiva para LowDensity (dropout 0.45, noise 0.30), ensemble de 3 seeds, y transfer learning (Global SNP → LowDensity SNP/SilicoDArT con freeze de encoder 10 epochs).

Ambas rondas: 3 seeds, GPU NVIDIA RTX 4060 Ti (16.7 GB), PyTorch 2.5.1 + CUDA 12.1.

### 3.4 Infraestructura - **Hardware**: CPU (sin GPU). - **Software**: Python 3.10, scikit-learn 1.3+, umap-learn 0.5.5, pandas 2.1+, NumPy 1.26+. - **Orquestación**: CLI unificado (run_experiments.py) con tracking automático en experiments/runs.jsonl. --- 
   
   ## 4. Resultados 
   
   ### 4.1 Calidad del embedding (Trustworthiness) | Dataset | Trustworthiness (media ± std) | Rango | |---------|-------------------------------|-------| | global_snp | 0.9835 ± 0.0002 | [0.9833, 0.9838] | | global_silico | 0.9864 ± 0.0001 | [0.9863, 0.9864] | | lowdensity_snp | 0.9436 ± 0.0031 | [0.9413, 0.9479] | | lowdensity_silico | 0.9214 ± 0.0011 | [0.9203, 0.9230] | | wild_snp | 0.9436 ± 0.0031 | [0.9413, 0.9479] | Todos los datasets superan 0.92, lo que indica una preservación excelente de las relaciones de vecindario local en la proyección 2D. ### 4.2 Estabilidad del grafo kNN | Dataset | Jaccard vecinos (PCA) | Jaccard aristas (PCA) | Jaccard vecinos (UMAP) | |---------|-----------------------|-----------------------|------------------------| | global_snp | 0.9887 ± 0.0021 | 0.9882 ± 0.0022 | 0.6521 ± 0.0026 | | global_silico | 1.0000 ± 0.0000 | 1.0000 ± 0.0000 | 0.6686 ± 0.0021 | | lowdensity_snp | 0.9894 ± 0.0012 | 0.9887 ± 0.0013 | 0.6543 ± 0.0022 | | lowdensity_silico | 1.0000 ± 0.0000 | 1.0000 ± 0.0000 | 0.7115 ± 0.0135 | | wild_snp | 0.9894 ± 0.0012 | 0.9887 ± 0.0013 | 0.6543 ± 0.0022 | Los grafos construidos sobre PCA son altamente estables (Jaccard ≥ 0.989), mientras que los vecindarios UMAP muestran mayor variabilidad (~0.65–0.71), consistente con la naturaleza estocástica de UMAP. ### 4.3 Sensibilidad a la imputación | Dataset | Estrategia | TrusVtworthiness | Jaccard vecinos vs. baseline | Jaccard aristas vs. baseline | |---------|-----------|-----------------|------------------------------|------------------------------| | global_snp | most_frequent (base) | 0.9835 | — | — | | global_snp | median | 0.9832 | 0.9154 | 0.9105 | | lowdensity_snp | most_frequent (base) | 0.9479 | — | — | | lowdensity_snp | median | 0.9444 | 0.9918 | 0.9912 | | global_silico | most_frequent (base) | 0.9865 | — | — | | global_silico | median | 0.9861 | 1.0000 | 1.0000 | Los datasets SilicoDArT (binarios) son insensibles al cambio de imputación. Los SNP muestran leves diferencias (Jaccard > 0.91 para global_snp, > 0.99 para lowdensity_snp), indicando robustez general. ### 4.4 Interpretabilidad (varianza explicada por PCA) | Dataset | PC1 (%) | PC2 (%) | PC1–5 acumulado (%) | |---------|---------|---------|---------------------| | global_snp | 9.49 | 2.38 | 15.76 | | global_silico | 10.20 | 2.72 | 17.09 | | lowdensity_snp | 8.84 | 5.26 | 21.58 | | lowdensity_silico | 99.99 | < 0.01 | ~100.00 | | wild_snp | 8.84 | 5.26 | 21.58 | El dataset lowdensity_silico presenta PC1 con varianza casi total (99.99%), lo que sugiere una estructura unidimensional dominante en los marcadores SilicoDArT de baja densidad.

### 4.5 Autoencoder vs. Baseline PCA

#### Trustworthiness (k = 15, bottleneck 64D vs. PCA 30D)

| Dataset | Baseline PCA | AE v1 (200 ep) | AE v2 reg (50 ep) | AE v2 ensemble | AE v2 transfer |
|---------|-------------|----------------|--------------------|--------------------|-----------------|
| global_snp | 0.976 | **0.988** | 0.982 | 0.985 | — |
| global_silico | 0.979 | **0.989** | 0.977 | 0.980 | — |
| lowdensity_snp | **0.933** | 0.836 | 0.844 | 0.784 | 0.658 |
| lowdensity_silico | **0.911** | 0.829 | 0.853 | 0.770 | 0.772 |

El AE mejora la trustworthiness en +1 pp para los datasets Global ($n \approx 6\,000$), pero la degrada en −8 a −10 pp para LowDensity ($n \approx 630$).

#### Estabilidad kNN (Edge Jaccard, seed-to-seed)

| Dataset | Baseline PCA | AE v1 | AE v2 reg | AE v2 seed↔ens |
|---------|-------------|-------|-----------|-------------|
| global_snp | 0.889 | 0.524 | 0.438 | 0.568 |
| global_silico | 0.885 | 0.498 | 0.333 | 0.466 |
| lowdensity_snp | 0.909 | 0.163 | 0.173 | 0.146 |
| lowdensity_silico | 0.908 | 0.207 | 0.233 | 0.178 |

La estabilidad del AE es sustancialmente inferior al baseline PCA en todos los datasets. Para LowDensity, la estabilidad del AE desciende a ~0.17 (vs. 0.91 del baseline), indicando que los grafos kNN derivados del AE no son reproducibles entre seeds.

 --- ## 5. Discusión ### Preservación de estructura local Los valores de trustworthiness (0.92–0.99) confirman que el pipeline PCA→UMAP preserva bien la estructura de vecindarios locales. Los datasets más grandes (global_snp, global_silico, ~6 000 muestras) obtienen los valores más altos, beneficiándose de la densidad muestral para la estimación de vecindarios. ### Estabilidad: PCA vs. UMAP La separación clara entre la estabilidad en PCA (~1.0) y en UMAP (~0.65) ilustra un trade-off conocido: PCA produce representaciones deterministas y estables, mientras UMAP introduce estocasticidad que afecta los vecindarios del embedding pero no necesariamente la calidad perceptual. La construcción del grafo kNN sobre componentes PCA (y no sobre el embedding UMAP) garantiza reproducibilidad del grafo. ### Elección de métrica El uso de coseno para SNP y Jaccard para SilicoDArT refleja la naturaleza de los datos: los SNP son cuantitativos (0/1/2 o dosificación alélica) donde la coseno captura similitud de perfil, mientras los SilicoDArT son presencia/ausencia donde Jaccard es la distancia natural. ### Sensibilidad a la imputación La baja sensibilidad observada (Jaccard > 0.91 en todos los casos) indica que el pipeline es robusto frente a la estrategia de imputación, un resultado deseable dado que los datos genómicos frecuentemente presentan tasas de missing del 5–18%. ### Estructura del dataset lowdensity_silico La concentración de varianza en PC1 (99.99%) para lowdensity_silico merece atención: podría reflejar un batch effect o una señal biológica genuinamente dominante (e.g., separación entre grupos taxonómicos). Esto requiere investigación adicional con metadatos de accesión.

### Autoencoder: beneficio condicionado al tamaño muestral

El autoencoder denoising de 64D mejora la preservación local (trustworthiness +1 pp) en los datasets Global ($n \approx 6\,000$), donde la relación muestras/parámetros (~1.2 M parámetros) es favorable. Sin embargo, en LowDensity ($n \approx 630$), la sobreparametrización (~1 900× más parámetros que muestras) produce un colapso parcial del bottleneck: las representaciones pierden discriminación local y el trust desciende 8–10 pp respecto al baseline PCA.

Las tres estrategias de mitigación evaluadas (\S 4.5) no resuelven el problema fundamental:
- **Regularización agresiva**: mejora marginal (+2 pp), insuficiente para cerrar la brecha.
- **Ensemble de bottleneck**: el promedio aritmético suaviza la estructura local que cada seed captura individualmente, empeorando el trust (−6 pp vs. semilla individual).
- **Transfer learning**: los datasets emplean paneles de markers distintos (20k vs. 62k loci), por lo que solo las capas internas se transfieren; las capas de entrada/salida se reinicializan, anulando el beneficio esperado (trust −18 pp en el peor caso).

Estos resultados sugieren que, para colecciones genómicas con $n < 1\,000$, la reducción lineal PCA sigue siendo preferible a autoencoders densos. Arquitecturas más eficientes en datos (e.g., transformers con masked prediction) podrían superar esta barrera (cf. \S 8, trabajo futuro).

 --- ## 6. Conclusiones 1. **GENO-MAP** proporciona un pipeline reproducible, minimalista y extensible para convertir matrices de genotipado DArT/DArTSeq en representaciones visualizables (embeddings 2D) y navegables (grafos kNN).
 2. La calidad de los embeddings, medida por trustworthiness (> 0.92), valida que las proyecciones 2D preservan fielmente la estructura local de la diversidad genómica. 3. La estabilidad del grafo kNN en el espacio PCA (Jaccard ≥ 0.989) asegura que las relaciones de similitud son reproducibles y no dependen de la semilla de UMAP. 4. El pipeline es robusto frente a estrategias de imputación alternativas, lo que simplifica las decisiones metodológicas para el usuario. 5. Un autoencoder denoising de 64D mejora la trustworthiness en +1 pp para datasets grandes ($n \approx 6\,000$), pero su estabilidad inter-seed y desempeño en datasets pequeños ($n < 700$) son inferiores al baseline PCA, confirmando que la reducción lineal sigue siendo más adecuada para colecciones con pocas muestras. 6. Las salidas JSON permiten una integración directa con herramientas de visualización interactiva (D3, Sigma.js, Gephi) para exploración exploratoria por curadores de germoplasma. --- ## 7. Limitaciones 1. **Imputación simplista**: Se utiliza imputación por moda/mediana sin considerar la estructura poblacional (e.g., imputación basada en LD o métodos como BEAGLE). Esto podría distorsionar relaciones en datasets con alta tasa de missing. 2. **Metadatos limitados**: Los datasets disponibles contienen metadatos escasos (solo Institution en algunos casos, con una única etiqueta), lo que impide evaluar la separación biológica real (silhouette, ARI). La misma etiqueta Institution = CIP/NA domina en la mayoría de las muestras. 3. **Validación restringida al embedding**: No se evalúa rendimiento downstream (e.g., clasificación de accesiones, predicción de rasgos fenotípicos). 4. **Escalabilidad**: Aunque funcional en CPU para ~6 000 muestras × 57 000 marcadores, el pipeline no escala eficientemente a datasets significativamente mayores sin GPU (RAPIDS/cuML). 5. **Métricas de grafos**: No se calculan métricas topológicas avanzadas del grafo kNN (modularidad, centralidad, detección de comunidades formal) que podrían revelar subestructura. 6. **Un solo cultivo**: Los resultados se limitan a camote; la generalización a otros cultivos andinos (yacón, ulluco, oca, mashua) referenciados en los metadatos del proyecto no se ha probado debido a restricciones de acceso a los datos. 7. **Autoencoder limitado a Level 1**: La exploración del autoencoder se restringió a arquitecturas densas (ResBlock) sin probar variantes convolucionales, attention-based, o modelos pre-entrenados. La mejora observada (+1 pp) en datasets grandes no justifica la complejidad adicional frente a PCA. --- ## 8. Trabajo futuro 1. **Representaciones no lineales avanzadas**: El autoencoder Level 1 mostró potencial limitado. Arquitecturas como Masked Genotype Transformers (Level 2), aprendizaje contrastivo (Level 3) y redes sobre grafos genómicos (Level 4) podrían superar la barrera de muestras pequeñas al incorporar inductive biases más adecuados para datos genómicos. 2. **Imputación aprendida**: Los autoencoders denoising reconstruyen la entrada completa a partir de versiones corruptas; esta capacidad podría explotarse directamente para imputación genotípica, reemplazando métodos simples (moda/mediana) por imputación basada en la estructura aprendida (Level 5). 2. **Detección de comunidades**: Aplicar algoritmos de clustering en grafos (Louvain, Leiden) sobre el grafo kNN para identificar formalmente subpoblaciones y comparar con clasificaciones taxonómicas tradicionales. 3. **Integración multi-cultivo**: Extender el pipeline a yacón, ulluco, oca y mashua cuando los datasets restringidos sean accesibles, habilitando análisis comparativos de diversidad entre cultivos. 4. **Visualización interactiva web**: Desarrollar un front-end (e.g., React + D3/Sigma.js) que consuma directamente los JSON de nodos y aristas, permitiendo a los curadores filtrar, buscar y explorar la diversidad genómica de forma interactiva. 5. **Validación biológica**: Correlacionar los clusters del embedding con datos fenotípicos (rendimiento, resistencia a plagas, contenido de betacaroteno en camote) para establecer relevancia agronómica de los agrupamientos. 6. **Escalabilidad GPU**: Integrar RAPIDS/cuML para UMAP y kNN acelerados por GPU, habilitando el procesamiento de colecciones con > 50 000 muestras. 7. **Recomendación de cruzas**: Explotar la estructura del grafo kNN para generar "recetas de cruza" que maximicen la diversidad genética entre progenitores seleccionados, optimizando heterosis y cobertura de alelos. 8. **Benchmarking de reducción de dimensionalidad**: Comparar UMAP con alternativas recientes (TriMap, PaCMAP, t-SNE parametrizado) en términos de trustworthiness, estabilidad y separación de grupos biológicos. --- ## Referencias - Venna, J., & Kaski, S. (2006). Local multidimensional scaling. *Neural Networks*, 19(6–7), 889–899. - McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. *arXiv:1802.03426*. - Jaccoud, D., et al. (2001). Diversity Arrays: a solid state technology for sequence information independent genotyping. *Nucleic Acids Research*, 29(4), e25. - Kilian, A., et al. (2012). Diversity Arrays Technology: A Generic Genome Profiling Technology on Open Platforms. *Methods in Molecular Biology*, 888, 67–89. - CIP Dataverse. International Potato Center genetic datasets. https://data.cipotato.org/ --- *Pipeline disponible en: [GENO-MAP – GitHub repository]* *Datos fuente: CIP Dataverse (DOIs referenciados en Tabla 1)*