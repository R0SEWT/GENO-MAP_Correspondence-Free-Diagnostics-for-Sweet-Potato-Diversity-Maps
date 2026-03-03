  # GENO-MAP: Correspondence-Free Diagnostics and Robustness Curves for Sweet Potato DArT/DArTSeq Diversity Maps

  ## Abstract

  Exploring genomic diversity in germplasm collections is challenged by ultra-wide genotype matrices and, in practice, by the lack of shared identifiers across independently generated genotyping panels. We present **GENO-MAP**, a reproducible pipeline that converts sweet potato (*Ipomoea batatas*) DArT/DArTSeq matrices into 2D visual maps and k-nearest neighbor (kNN) graphs for exploration and curation. GENO-MAP applies explicit imputation, PCA as a fixed analytical space for graph construction, UMAP as a visualization layer, and exports structured JSON artifacts for interactive use.

  Because the studied panels use disjoint identifier namespaces (CIP accession codes vs DArT plate/well coordinates), cross-panel alignment is infeasible without an external manifest. We therefore introduce a **correspondence-free validation** framework combining (i) per-panel geometry diagnostics (effective rank, PC dominance, kNN reciprocity, connected components, QA flags) and (ii) robustness curves under controlled perturbations (marker subsampling, MCAR missing injection, and imputation strategy changes). Across four panels (630–5,970 samples × 20,069–62,732 markers), PCA subspace structure remains stable under severe marker dropout (subspace similarity ≥ 0.91 at 5% markers), while local kNN neighborhoods degrade smoothly with perturbation and remain largely insensitive to mode-vs-median imputation (Jaccard ≥ 0.91). A denoiVsing autoencoder benchmark yields at most marginal trustworthiness gains at higher sample-to-feature ratios, but exhibits substantially lower topological stability than PCA, supporting PCA as the most reliable baseline under a fixed graph-construction protocol in the studied n≪p regimes.

  **Keywords:** genomic diversity, correspondence-free validation, kNN graphs, PCA, robustness curves, panel diagnostics, germplasm, sweet potato.

  ---

  ## 1. Introducción

  Los bancos de germoplasma de cultivos andinos albergan miles de accesiones genotipadas con tecnologías DArT y DArTSeq, generando matrices anchas (muestras × marcadores) con decenas de miles de loci. Explorar visualmente la estructura de diversidad en estas matrices es esencial para (i) identificar grupos genéticos, (ii) detectar redundancias o brechas en la colección, y (iii) priorizar combinaciones de cruza que maximicen la diversidad genética. Sin embargo, las dimensiones de estos datos (> 50 000 marcadores) dificultan su visualización directa, y una complicación adicional surge en la práctica: los distintos paneles de genotipado emplean *espacios de identificadores disjuntos* (accesiones CIP vs. coordenadas plate/well DArT), lo que imposibilita el alignment cruzado sin un manifest externo. En consecuencia, este trabajo prioriza **robustez intra-cultivo** y comparación entre paneles **sin correspondencia**, en lugar de generalización multi-cultivo o alineamiento inter-panel espurio.

  En este trabajo presentamos GENO-MAP, un pipeline de código abierto que automatiza la cadena completa: lectura de matrices DArT/DArTSeq → imputación → reducción PCA+UMAP → grafo kNN → artefactos JSON para visualización. Más allá del pipeline de producción, desarrollamos un marco de validación **correspondence-free** que opera rigurosamente dentro de la limitación de IDs disjuntos, caracterizando cada panel individualmente mediante diagnósticos geométricos y evaluando la robustez de la topología del grafo bajo perturbaciones controladas. Como extensión, benchmarkeamos un autoencoder denoising como alternativa no lineal a PCA para la etapa de reducción, incluyendo un análisis de frontera de estabilidad representacional que identifica el crossover muestral $n^*$ donde el autoencoder iguala al baseline.

  ---

  ## 2. Datos

  Se utilizaron tres fuentes de datos de camote (*Ipomoea batatas*) provenientes del repositorio CIP Dataverse. Tras detectar que `wild_snp` es un duplicado del archivo `lowdensity_snp` bajo hashing de matriz y coincidencia de dimensiones, se excluye del análisis (ADR-0002), reteniendo cuatro datasets.

  | Dataset | DOI | Tipo | $n$ | $p$ | $n/p$ | Missing |
  |---------|-----|------|-----|-----|-------|---------|
  | global_snp | 10.21223/P30BVZYY | SNP | 5 970 | 20 069 | 0.297 | 5.52% |
  | global_silico | 10.21223/P30BVZYY | SilicoDArT | 5 970 | 57 715 | 0.103 | 3.21% |
  | lowdensity_snp | 10.21223/P3UBDJ44 | SNP | 630 | 62 732 | 0.010 | 18.17% |
  | lowdensity_silico | 10.21223/P3UBDJ44 | SilicoDArT | 635 | 38 272 | 0.017 | 4.65% |

  **Tabla 1.** Datasets utilizados. Los paneles Global y LowDensity emplean espacios de identificadores disjuntos (accesiones CIP vs. coordenadas plate/well DArT). El ratio $n/p$ varía de 0.010 a 0.297, reflejando la naturaleza ultra-ancha de las matrices genómicas.

  Los formatos de entrada se detectan automáticamente: *sample_columns* (primera columna = etiqueta de fila, columnas = muestras) y *marker_metrics* (columnas iniciales de metadatos, seguidas de columnas numéricas). La matriz se transpone a orientación muestras × marcadores.

  ---

  ## 3. Metodología

  ### 3.1 Pipeline de procesamiento

  El pipeline se implementa en Python 3.11.14 (vía uv) con dependencias mínimas: NumPy, pandas, scikit-learn, umap-learn 0.5.11 y matplotlib. El flujo consta de cuatro etapas:

  1. **Lectura y parseo**: Detección automática de separador (, / ;) y formato. Se separan filas de metadatos (no numéricas) de filas de marcadores.
  2. **Imputación**: Los valores faltantes (codificados como `-`, vacío, o NaN) se imputan mediante moda (`most_frequent`) por defecto con `sklearn.impute.SimpleImputer`. Se evalúa alternativamente la mediana como sensibilidad.
  3. **Reducción de dimensionalidad**:
    - **PCA**: $\min(50, p, n)$ componentes para representación intermedia. Las primeras 30 PCs se usan como features deterministas para el grafo.
    - **UMAP**: 2 dimensiones de salida sobre las PCs, exclusivamente para visualización.
  4. **Grafo kNN**: Sobre las primeras 30 PCs con `sklearn.neighbors.NearestNeighbors`, $k = 15$ y métrica coseno. Las aristas se almacenan como tuplas (source, target, distance). Para comparabilidad *cross-panel* en el marco correspondence-free, se fija una única métrica (coseno) en el espacio PCA-30D para todos los paneles; el efecto de métricas alternativas (p.ej., Jaccard en binarios) queda como trabajo futuro.

  ### 3.2 Salidas

  Para cada dataset se generan tres archivos JSON: `*_nodes.json` (nodos con id, embedding $[x, y]$ y metadatos), `*_edges.json` (aristas kNN), y `*_stats.json` (estadísticas del dataset).

  ### 3.3 Marco de validación (correspondence-free)

  Dado que los paneles emplean espacios de identificadores disjuntos, toda validación opera intra-panel. Se evalúan cuatro dimensiones:

  1. **Calidad del embedding**: Trustworthiness (Venna & Kaski, 2006) — preservación de vecindarios locales entre espacio original y embedding 2D.
  2. **Estabilidad**: Múltiples semillas (42, 52, 62); solapamiento kNN (Jaccard medio por nodo y global de aristas) en espacio PCA y UMAP.
  3. **Sensibilidad a la imputación**: Topología del grafo entre moda y mediana — Jaccard de vecindarios y aristas vs. línea base.
  4. **Interpretabilidad de marcadores**: Loadings PCA de los 10 marcadores con mayor contribución a PC1/PC2, varianza explicada acumulada.

  ### 3.4 Panel QA & Geometry Diagnostics (Módulo 1)

  Caracterización automática de la geometría de cada panel:

  - **Rango efectivo** (participation ratio): $\text{eff\_rank} = (\sum \lambda_i)^2 / \sum \lambda_i^2$, mide dimensionalidad intrínseca.
  - **Dominancia PC1**: Ratio varianza PC1/PC2 y porcentaje de varianza de PC1.
  - **Estadísticas de distancias kNN**: Media, CV, asimetría y curtosis de la distribución de distancias en PCA-30D (métrica coseno).
  - **Tasa de reciprocidad kNN**: Fracción de aristas mutuas ($i \to j$ AND $j \to i$).
  - **Componentes conexas**: Número y tamaño del giant component, reportadas sobre el grafo **no dirigido** obtenido al simetrizar aristas kNN (union), para detectar outliers y subgrupos aislados.
  - **Flags automáticos**: `1D-DOMINANT` (PC1 > 50%), `COLLAPSED-NEIGHBORHOOD` (CV < 0.05), `HIGH-MISSINGNESS` (> 10%), `EXTREME-WIDE` ($n/p < 0.02$), `DISCONNECTED` (> 1 componente), `LOW-RECIPROCITY` (< 0.40).

  ### 3.5 Robustness Curves (Módulo 2)

  Análisis de sensibilidad estructural bajo perturbaciones controladas:

  - **Marker subsampling**: Retención aleatoria al {5, 10, 20, 50, 80}% de columnas.
  - **Missing injection**: Inyección de +{0, 5, 10, 20}% MCAR sobre los datos numéricos, seguida de re-imputación.
  - **Imputation comparison**: Moda vs. mediana sobre los mismos datos.

  Métricas de comparación (referencia = pipeline completo sin perturbación):

  - **Jaccard de vecinos** ($J_{\text{nbr}}$): Solapamiento medio por nodo de los conjuntos kNN.
  - **Drift de distancias kNN**: Estadístico KS de dos muestras entre distribuciones de distancias.
  - **Similitud del subespacio PCA**: Ángulos principales entre subspace scores (QR → SVD), promediar cosenos.

  Cada combinación se ejecuta con 3 seeds (42, 52, 62); se reportan medias ± std.

  ### 3.6 Autoencoder como alternativa a PCA

  Como alternativa no lineal, se implementa un autoencoder denoising (DAE):

  **Arquitectura** (`GenoAutoencoder`): Encoder Dense(input → 512) → ResBlock ×2 → Dense(512 → 64); Decoder simétrico. Cada ResBlock: LayerNorm → Dense → GELU → Dropout → Dense + skip connection. Bottleneck de 64 dimensiones.

  **Entrenamiento**: Denoising (15% mask), pérdida MSE + sparsity penalty, AdamW, ReduceLROnPlateau, early stopping. 3 seeds para evaluación de estabilidad. GPU NVIDIA RTX 4060 Ti (16.7 GB), PyTorch 2.5.1+cu121.

  **Variantes evaluadas**:

  - **ae-v1**: 200 epochs, patience 20, configuración base (dropout 0.20, noise 0.15).
  - **ae-v2**: 50 epochs, patience 15, con tres estrategias de mitigación para datasets pequeños:
    - Regularización agresiva (dropout 0.45, noise 0.30, weight decay 5×10⁻⁴).
    - Ensemble de bottleneck (promedio aritmético de 3 seeds).
    - Transfer learning (Global → LowDensity, freeze encoder 10 epochs).

  **Frontera de estabilidad** (Level A): Subsampleo de Global SNP a 14 tamaños muestrales (50–5 970) × 3 seeds × 2 métodos, midiendo trustworthiness y estabilidad kNN para identificar el crossover $n^*$.

  ### 3.7 Infraestructura

  - **Hardware**: NVIDIA GeForce RTX 4060 Ti (16.7 GB VRAM). GPU utilizada para entrenamiento AE y kNN acelerado.
  - **Software**: Python 3.11.14 (vía uv), PyTorch 2.5.1+cu121, scikit-learn 1.3+, umap-learn 0.5.11, pandas 2.1+, NumPy 1.26+.
  - **Orquestación**: CLI unificado (`run_experiments.py`) con tracking en `experiments/runs.jsonl`.

  ---

  ## 4. Resultados

  ### 4.1 Calidad del embedding (Trustworthiness)

  | Dataset | Trust (media ± std) | Rango |
  |---------|---------------------|-------|
  | global_snp | 0.9835 ± 0.0002 | [0.9833, 0.9838] |
  | global_silico | 0.9864 ± 0.0001 | [0.9863, 0.9864] |
  | lowdensity_snp | 0.9436 ± 0.0031 | [0.9413, 0.9479] |
  | lowdensity_silico | 0.9214 ± 0.0011 | [0.9203, 0.9230] |

  **Tabla 2.** Trustworthiness del pipeline PCA→UMAP (k=15, 3 seeds). Todos los datasets superan 0.92.

  ### 4.2 Estabilidad del grafo kNN

  | Dataset | Jaccard vecinos (PCA) | Jaccard aristas (PCA) | Jaccard vecinos (UMAP) |
  |---------|-----------------------|-----------------------|------------------------|
  | global_snp | 0.989 ± 0.002 | 0.988 ± 0.002 | 0.652 ± 0.003 |
  | global_silico | 1.000 ± 0.000 | 1.000 ± 0.000 | 0.669 ± 0.002 |
  | lowdensity_snp | 0.989 ± 0.001 | 0.989 ± 0.001 | 0.654 ± 0.002 |
  | lowdensity_silico | 1.000 ± 0.000 | 1.000 ± 0.000 | 0.712 ± 0.014 |

  **Tabla 3.** Estabilidad inter-seed del grafo kNN cuando el grafo se construye en el espacio PCA-30D con la misma configuración y datos (k=15, coseno). La variabilidad proviene principalmente de la capa UMAP (visualización), no del espacio PCA.

  ### 4.3 Sensibilidad a la imputación

  | Dataset | Estrategia | Trust | J vecinos vs. base | J aristas vs. base |
  |---------|-----------|-------|---------------------|---------------------|
  | global_snp | most_frequent | 0.984 | — | — |
  | global_snp | median | 0.983 | 0.915 | 0.911 |
  | global_silico | most_frequent | 0.987 | — | — |
  | global_silico | median | 0.986 | 1.000 | 1.000 |
  | lowdensity_snp | most_frequent | 0.948 | — | — |
  | lowdensity_snp | median | 0.944 | 0.992 | 0.991 |

  **Tabla 4.** Sensibilidad a imputación. SilicoDArT binario es insensible; SNP muestra leves diferencias (Jaccard > 0.91).

  ### 4.4 Panel QA & Geometry Diagnostics

  | Panel | $n$ | $p$ | $n/p$ | Miss% | PC1% | Eff. Rank | Recip. | Comp. | Flags |
  |-------|-----|-----|-------|-------|------|-----------|--------|-------|-------|
  | Global SNP | 5 970 | 20 069 | 0.297 | 5.5 | 9.49 | 10.4 | 0.648 | 9 | DISCONNECTED |
  | Global Silico | 5 970 | 57 715 | 0.103 | 3.2 | 10.20 | 10.0 | 0.647 | 5 | DISCONNECTED |
  | LD SNP | 630 | 62 732 | 0.010 | 18.2 | 8.87 | 16.2 | 0.699 | 2 | HIGH-MISS; DISCONN; EXTREME-WIDE |
  | LD Silico | 635 | 38 272 | 0.017 | 4.7 | 9.46 | 16.1 | 0.686 | 1 | EXTREME-WIDE |

  **Tabla 5.** Diagnósticos geométricos por panel (k=15, coseno en PCA-30D). Los cuatro paneles presentan rango efectivo moderado (10–16), sin colapso 1D (PC1 < 11%), y reciprocidad kNN saludable (0.65–0.70). Los paneles LowDensity son flaggueados como EXTREME-WIDE ($n/p \leq 0.017$). Los small components (2–9) en los datasets SNP son compatibles con outliers o estratificación fuerte; se reportan como señal de QA y no como subpoblaciones biológicas confirmadas en ausencia de metadatos informativos.

  ### 4.5 Robustness Curves

  #### 4.5.1 Marker subsampling

  | Markers | Global SNP | Global Silico | LD SNP | LD Silico |
  |---------|-----------|---------------|--------|-----------|
  | 5% | $J_{\text{nbr}}$ = 0.43, SS = 0.92 | 0.56, 0.94 | 0.59, 0.97 | 0.53, 0.93 |
  | 10% | 0.51, 0.95 | 0.65, 0.98 | 0.66, 0.99 | 0.61, 0.97 |
  | 20% | 0.61, 0.98 | 0.73, 0.99 | 0.73, 0.99 | 0.68, 0.99 |
  | 50% | 0.75, 1.00 | 0.84, 1.00 | 0.83, 1.00 | 0.80, 1.00 |
  | 80% | 0.84, 1.00 | 0.90, 1.00 | 0.89, 1.00 | 0.88, 1.00 |

  **Tabla 6.** Jaccard de vecinos ($J_{\text{nbr}}$) y similitud del subespacio PCA (SS) bajo marker subsampling (3 seeds). La SS se mantiene > 0.91 incluso al 5%, indicando que la geometría principal es robusta al dropout de marcadores. $J_{\text{nbr}}$ degrada gracefully: al 20%, se preserva ~61–73% de los vecindarios locales.

  #### 4.5.2 Missing injection (MCAR)

  | Extra miss. | Global SNP | Global Silico | LD SNP | LD Silico |
  |-------------|-----------|---------------|--------|-----------|
  | +0% | 0.99, 1.00 | 0.99, 1.00 | 0.99, 1.00 | 0.99, 1.00 |
  | +5% | 0.88, 1.00 | 0.92, 1.00 | 0.89, 1.00 | 0.89, 1.00 |
  | +10% | 0.83, 1.00 | 0.90, 1.00 | 0.87, 1.00 | 0.85, 1.00 |
  | +20% | 0.77, 1.00 | 0.86, 1.00 | 0.83, 1.00 | 0.81, 1.00 |

  **Tabla 7.** Jaccard de vecinos y similitud subespacio bajo inyección de missing MCAR (3 seeds). El subespacio PCA es prácticamente invariante (SS > 0.998 en todos los casos). La topología local es más sensible: +20% MCAR reduce $J_{\text{nbr}}$ a 0.77–0.86, con Global SNP siendo el más afectado (probablemente por la mayor varianza en patrones de imputación).

  #### 4.5.3 Imputation comparison (mode vs. median)

  | Panel | $J_{\text{nbr}}$ | SS |
  |-------|-------------------|-----|
  | Global SNP | 0.907 | 1.000 |
  | Global Silico | 1.000 | 1.000 |
  | LD SNP | 0.993 | 1.000 |
  | LD Silico | 0.995 | 1.000 |

  **Tabla 8.** Efecto de cambiar moda por mediana. SilicoDArT es completamente insensible (datos binarios → mediana ≈ moda). SNP muestra ligero efecto en Global SNP (Jaccard 0.91), donde la mayor variabilidad alélica (0/1/2) permite que moda y mediana difieran.

  ### 4.6 Autoencoder vs. Baseline PCA

  #### Trustworthiness (k=15, bottleneck 64D vs. PCA 30D)

  | Dataset | PCA | AE v1 (200ep) | AE v2 reg | AE v2 ens. | AE v2 transfer |
  |---------|-----|---------------|-----------|------------|----------------|
  | global_snp | 0.976 | **0.988** | 0.982 | 0.985 | — |
  | global_silico | 0.979 | **0.989** | 0.977 | 0.980 | — |
  | lowdensity_snp | **0.933** | 0.836 | 0.844 | 0.784 | 0.658 |
  | lowdensity_silico | **0.911** | 0.829 | 0.853 | 0.770 | 0.772 |

  **Tabla 9.** El AE mejora trust en +1 pp para Global ($n \approx 6\,000$) pero la degrada en −8 a −10 pp para LowDensity ($n \approx 630$).

  #### Estabilidad kNN del autoencoder (Edge Jaccard, seed-to-seed)

  | Dataset | PCA | AE v1 | AE v2 reg | AE v2 ens. |
  |---------|-----|-------|-----------|------------|
  | global_snp | 0.889 | 0.524 | 0.438 | 0.568 |
  | global_silico | 0.885 | 0.498 | 0.333 | 0.466 |
  | lowdensity_snp | 0.909 | 0.163 | 0.173 | 0.146 |
  | lowdensity_silico | 0.908 | 0.207 | 0.233 | 0.178 |

  **Tabla 10.** Estabilidad del grafo derivado de embeddings aprendidos (AE) frente a variación de inicialización/entrenamiento. Para PCA, se reporta el mismo protocolo de comparación aplicado en el experimento de embedding (seeds y pipeline completos), lo que incluye fuentes adicionales de variación (p.ej., normalización/entrenamiento/serialización) respecto al baseline de Tabla 3. Para LowDensity, Jaccard ~0.17 (vs. 0.91 PCA), indicando grafos no reproducibles entre seeds.

  > *Nota:* para evitar ambigüedad, en la versión camera-ready se unificará la definición operativa de estabilidad en un único protocolo.

  #### Frontera de estabilidad (Level A)

  El crossover muestral se sitúa en $n^* \in [4\,000, 4\,500]$ ($n/p \approx 0.20$–0.22 para Global SNP con $p = 20\,069$), donde AE alcanza trust igual o ligeramente superior a PCA ($\Delta \leq +0.002$). Sin embargo, la estabilidad kNN del AE nunca converge a la de PCA: incluso a $n = 5\,970$ (máximo disponible), el gap de estabilidad es ~0.46 Jaccard units (AE: 0.524 vs. PCA: 0.889).

  ---

  ## 5. Discusión

  ### 5.1 Preservación de estructura local

  Los valores de trustworthiness (0.92–0.99) confirman que el pipeline PCA→UMAP preserva bien la estructura de vecindarios locales. Los datasets más grandes obtienen los valores más altos, beneficiándose de la densidad muestral.

  ### 5.2 Estabilidad: PCA vs. UMAP

  La separación entre estabilidad PCA (~1.0) y UMAP (~0.65) ilustra un trade-off conocido: PCA produce representaciones deterministas, mientras UMAP introduce estocasticidad que afecta los vecindarios del embedding pero no la calidad perceptual. La decisión de construir el grafo kNN sobre PCs — y no sobre el embedding UMAP — garantiza reproducibilidad.

  ### 5.3 Geometría de los paneles

  Los cuatro paneles muestran rango efectivo moderado (10–16 dimensiones con varianza significativa), sin colapso unidimensional (PC1 < 11%). La reciprocidad kNN de 0.65–0.70 indica vecindarios saludables donde la mayoría de relaciones son bidireccionales. Las pequeñas componentes desconexas en los paneles SNP (2–9 componentes) son compatibles con outliers o estratificación fuerte; se reportan como señal de QA y no como subpoblaciones biológicas confirmadas en ausencia de metadatos informativos.

  Los flags automáticos `EXTREME-WIDE` para los paneles LowDensity ($n/p \leq 0.017$) alertan sobre la sobreparametrización implícita: con ~60× más marcadores que muestras, la estimación de relaciones de similitud depende fuertemente de la etapa de reducción dimensional.

  ### 5.4 Robustez bajo perturbaciones

  **Marker subsampling.** La similitud del subespacio PCA permanece por encima de 0.91 incluso al retener solo 5% de los marcadores (~1 000–3 000 loci). Esto confirma que la estructura genómica principal está distribuida a lo largo del genoma y no depende de un subconjunto específico de marcadores. Los paneles con más marcadores (SilicoDArT: 38–57k) toleran mejor el dropout ($J_{\text{nbr}}$ = 0.56 vs. 0.43 al 5%), probablemente porque la redundancia informativa es mayor.

  **Missing injection.** El subespacio PCA es prácticamente invariante a +20% MCAR (SS > 0.998), pero la topología kNN local es más sensible ($J_{\text{nbr}}$ = 0.77–0.86). La degradación es monótona y predecible, lo que permite establecer umbrales de calidad para colecciones con alta tasa de missing.

  **Imputation.** El cambio moda↔mediana produce un efecto menor, con Jaccard > 0.91 en todos los paneles. Los datos SilicoDArT (binarios) son completamente insensibles (Jaccard = 1.0). Esto confirma que la elección de imputación simple no es una decisión crítica para la topología del grafo.

  ### 5.5 Autoencoder: beneficio condicionado y frontera de viabilidad

  El autoencoder denoising de 64D mejora marginalmente la trustworthiness en los datasets Global (+1 pp), pero con un costo severo en estabilidad inter-seed. En LowDensity ($n \approx 630$), la sobreparametrización (~1 900× más parámetros que muestras) produce un colapso parcial del bottleneck, y el trust desciende 8–10 pp.

  Las tres estrategias de mitigación evaluadas no resuelven el problema fundamental:

  1. **Regularización agresiva**: mejora marginal (+2 pp), insuficiente para cerrar la brecha. La regularización no puede sustituir datos.
  2. **Ensemble de bottleneck**: el promedio aritmético suaviza la estructura local ($-6$ pp vs. semilla individual), neutralizando la ventaja del AE.
  3. **Transfer learning**: los paneles emplean conjuntos de marcadores distintos (20k vs. 62k loci), así que solo las capas internas se transfieren; las capas de entrada/salida se reinicializan, produciendo el peor resultado ($-18$ pp).

  La frontera de estabilidad (Fig. 10) muestra que en nuestros paneles de camote, cuando $n/p < 0.20$ el AE no muestra ventajas robustas y pierde estabilidad; en este régimen, PCA resulta el baseline más confiable. El crossover a $n^* \approx 4\,000$–$4\,500$ ofrece una ventaja negligible ($\Delta \leq 0.002$) que no justifica la complejidad adicional ni la pérdida de estabilidad.

  ### 5.6 Limitación de IDs disjuntos y framing

  Los paneles de genotipado operan con namespaces disjuntos: Global emplea accesiones CIP mientras LowDensity usa coordenadas plate/well DArT. Sin un manifest externo que mapee plate/well → CIP accession, cualquier claim de alignment cruzado sería *hallucinated alignment*. El marco correspondence-free adoptado en este trabajo (diagnósticos geométricos + curvas de robustez) proporciona validación rigurosa sin requirir correspondencia directa, respondiendo preguntas como:

  > *¿El panel induce una geometría rica o degenerada?*
  > *¿La topología del grafo es sensible a decisiones metodológicas razonables?*
  > *¿LowDensity es más frágil que Global ante missingness y marker dropout?*

  La respuesta a la tercera pregunta es matizada: LowDensity no es sistemáticamente más frágil. Su $J_{\text{nbr}}$ bajo marker subsampling es comparable o superior al de Global SNP, probablemente porque el mayor número de marcadores (62k vs. 20k) compensa el menor $n$.

  ---

  ## 6. Conclusiones

  1. **GENO-MAP** proporciona un pipeline reproducible y minimalista para convertir matrices DArT/DArTSeq en embeddings 2D y grafos kNN, con trustworthiness > 0.92 y estabilidad kNN (Jaccard ≥ 0.989) en espacio PCA.

  2. El marco **correspondence-free** permite validación rigurosa de paneles con IDs disjuntos, sin incurrir en alignment espurio. Los diagnósticos geométricos (rango efectivo, reciprocidad, flags QA) y las curvas de robustez (marker subsampling, missing injection, imputation) operan intra-panel y son directamente comunicables a curadores de germoplasma.

  3. La **robustez del pipeline** es notable: la similitud del subespacio PCA se mantiene > 0.91 incluso al retener solo 5% de los marcadores, y la topología kNN es insensible a la estrategia de imputación (Jaccard > 0.91 en todos los casos).

  4. El autoencoder denoising muestra mejoras de trustworthiness pequeñas y condicionadas al tamaño muestral, pero mantiene un déficit sustancial de estabilidad topológica frente a PCA. En los regímenes $n \ll p$ observados en estos paneles de camote, PCA permanece como el baseline más estable y operacionalmente recomendable.

  5. Los **flags automáticos** de QA (EXTREME-WIDE, DISCONNECTED, HIGH-MISSINGNESS) constituyen una herramienta inmediata de control de calidad para bancos de germoplasma, alertando sobre condiciones que podrían comprometer interpretaciones downstream.

  6. Las salidas JSON permiten integración directa con herramientas de visualización interactiva (D3, Sigma.js, Gephi) para exploración por curadores.

  ---

  ## 7. Limitaciones

  1. **Imputación simplista**: Moda/mediana sin estructura poblacional (e.g., BEAGLE, LD-based). Podría distorsionar relaciones en datasets con alta tasa de missing, aunque los análisis de robustez (§4.5.2) muestran que el impacto es gradual y predecible.

  2. **Metadatos limitados**: Solo Institution disponible en algunos datasets, con etiqueta única (CIP/NA), impidiendo validación biológica formal (silhouette, ARI).

  3. **Validación restringida al embedding**: No se evalúa rendimiento downstream (clasificación, predicción fenotípica).

  4. **Métricas topológicas**: No se calculan modularidad, centralidad, ni detección formal de comunidades sobre el grafo kNN.

  5. **Un solo cultivo**: Resultados limitados a camote (*Ipomoea batatas*); generalización a otros cultivos andinos no probada. El alcance del estudio es **intra-cultivo**: robustez y QA comparativos entre paneles de camote; no se pretende generalización cross-crop.

  6. **Autoencoder limitado a arquitecturas densas**: No se exploraron variantes convolucionales, attention-based, ni modelos pre-entrenados genómicos.

  7. **IDs disjuntos entre paneles**: Los paneles Global y LowDensity emplean namespaces de identificadores no solapados (accesiones CIP vs. plate/well DArT). Sin un manifest DArT (plate/well → CIP accession), la correspondencia cruzada es imposible. Todas las comparaciones y análisis de robustez se realizan intra-panel (*correspondence-free validation*).

  8. **Perturbaciones limitadas**: Las curvas de robustez asumen MCAR (Missing Completely At Random). El missing real en datos genómicos puede ser MNAR (e.g., GC-bias, longitud de fragmento), lo que podría producir degradaciones no uniformes no capturadas por nuestro protocolo.

  ---

  ## 8. Trabajo futuro

  ### Módulo 3 — Decision Support (core-set selection)

  Selección de subconjuntos representativos *dentro de cada panel* para optimizar la cobertura de diversidad genómica. Incluye selección greedy de core-sets en el espacio PCA-30D, curva de diversidad vs. tamaño $k$, y exportación de IDs seleccionados para uso directo en bancos de germoplasma.

  ### Módulo 4 — Manifest acquisition

  La correspondencia entre paneles (Global ↔ LowDensity) requiere un manifest DArT que mapee plate/well a accesiones CIP. Se documenta como limitación estructural explícita y acción inmediata: obtener dicho manifest para habilitar batch correction, transfer alignment, y validación cruzada.

  ### Líneas complementarias

  1. **Detección de comunidades**: Louvain/Leiden sobre kNN para identificar subpoblaciones.
  2. **Integración multi-cultivo**: Extensión a yacón, ulluco, oca y mashua.
  3. **Visualización interactiva web**: Front-end React + D3/Sigma.js.
  4. **Validación biológica**: Correlación con datos fenotípicos.
  5. **Escalabilidad GPU**: RAPIDS/cuML para > 50 000 muestras.
  6. **Recomendación de cruzas**: Optimización de heterosis via grafos kNN.
  7. **Perturbaciones MNAR**: Modelar missing no aleatorio (GC-bias, locus-dependent) para curvas de robustez más realistas.

  ---

  ## Figuras

  - **Fig. 1–9**: Embeddings 2D, grafos kNN, loadings PCA (ver notebook §1–9).
  - **Fig. 10**: Frontera de estabilidad representacional — crossover $n^* \approx 4\,000$–$4\,500$, gap de estabilidad ~0.46 (notebook §13).
  - **Fig. 11**: Diagnósticos geométricos por panel — varianza acumulada PCA, rango efectivo, ratio PC1/PC2 (notebook §14).
  - **Fig. 12**: Curvas de marker subsampling — Jaccard neighbors, similitud subespacio, KS drift (notebook §15).
  - **Fig. 13**: Curvas de missing injection MCAR — topología local vs. global (notebook §15).
  - **Fig. 14**: Resumen de robustez — imputation sensitivity + heatmap de Jaccard por fracción (notebook §15).

  ---

  ## Referencias

  - Venna, J., & Kaski, S. (2006). Local multidimensional scaling. *Neural Networks*, 19(6–7), 889–899.
  - McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. *arXiv:1802.03426*.
  - Jaccoud, D., et al. (2001). Diversity Arrays: a solid state technology for sequence information independent genotyping. *Nucleic Acids Research*, 29(4), e25.
  - Kilian, A., et al. (2012). Diversity Arrays Technology: A Generic Genome Profiling Technology on Open Platforms. *Methods in Molecular Biology*, 888, 67–89.
  - Roy, J., et al. (2005). Participation ratio and minimum embedding dimension. *Physical Review E*, 72(2), 026106.
  - Franco, J., et al. (2006). A sampling strategy for conserving genetic diversity when forming core subsets. *Crop Science*, 46(2), 854–861.
  - CIP Dataverse. International Potato Center genetic datasets. https://data.cipotato.org/

  ---

  *Pipeline disponible en: [GENO-MAP – GitHub repository]*
  *Datos fuente: CIP Dataverse (DOIs referenciados en Tabla 1)*
  *Fecha de compilación: marzo 2026*
