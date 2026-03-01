# 🔬 EDA: Variables y Missingness — OXOR / GENO-MAP

> **Generado automáticamente** por `notebooks/eda_variables.ipynb`  
> **Fecha:** 2026-03-01 01:41

---

## TL;DR

- **Archivos analizados:** 4 genotipos + 3 métricas
- ⚠️ **Duplicado confirmado:** LowDensity SNP y Wild SNP son el MISMO archivo (MD5 idéntico)
- **Global SNP:** 1,998×5,970, missingness global = 6.5%
- **Global SilicoDArT:** 1,998×5,970, missingness global = 2.3%
- **LowDensity SNP:** 2,000×635, missingness global = 18.7%
- **LowDensity SilicoDArT:** 2,000×635, missingness global = 8.2%
- 🔍 **Missingness alta (>10%):** LowDensity SNP — revisar heatmaps para patrones de parche

---

## 1. Diccionario de Variables (Data Dictionary xlsx)


### Hoja: Hoja1 (236,610 filas → 4 definiciones únicas)


| Variable name   | Abbreviation   | Context of use   | Description                                | Unit   | Category 1   | Category 2   | Category 3      |
|:----------------|:---------------|:-----------------|:-------------------------------------------|:-------|:-------------|:-------------|:----------------|
| Sample_code     | N              | sample code      | Numeric identifier of samples              | number | nan          | nan          | nan             |
| Country_origin  | ORIGCNTRY      | Passport data    | Country of origin of sweetpotato germplasm | text   | nan          | nan          | nan             |
| Accession       | ACCSNUMBER     | Passport data    | Accession identifier                       | text   | nan          | nan          | nan             |
| 7529063         | DARTGENOTYPE   | Data             | Presence/absence of a SilicoDArT fragment. | number | 1 =presence  | 0 =absence   | 9 =missing data |


## 2. Inventario de Columnas por Archivo


| Archivo                   |   Total cols |   ID |   Secuencia |   Genómica |   Calidad |   Muestra |   Otro | Sep   |
|:--------------------------|-------------:|-----:|------------:|-----------:|----------:|----------:|-------:|:------|
| Global SNP geno           |         5971 |    0 |           0 |          0 |         0 |      5970 |      1 | ','   |
| Global SNP metrics        |           25 |    3 |           2 |         10 |        10 |         0 |      0 | ','   |
| Global SilicoDArT geno    |         5971 |    0 |           0 |          0 |         0 |      5970 |      1 | ','   |
| Global SilicoDArT metrics |           14 |    1 |           2 |          4 |         7 |         0 |      0 | ','   |
| LowDensity SNP            |          667 |    2 |           1 |         10 |        11 |       635 |      8 | ','   |
| LowDensity SilicoDArT     |          650 |    1 |           1 |          0 |         9 |       635 |      4 | ','   |
| Wild SNP                  |          667 |    2 |           1 |         10 |        11 |       635 |      8 | ','   |
| Wild SilicoDArT metrics   |           15 |    1 |           2 |          5 |         7 |         0 |      0 | ';'   |


### Detalle de columnas no-muestra



**Global SNP geno** (1 metadatos, 5970 muestras)


| Columna | Categoría |
|---------|-----------|

| ﻿Sample_code | Otro |




**Global SNP metrics** (25 metadatos, 0 muestras)


| Columna | Categoría |
|---------|-----------|

| ﻿AlleleID | ID |

| CloneID | ID |

| AlleleSequence | Secuencia |

| TrimmedSequence | Secuencia |

| Chrom_SweetPotato_unkn_v4 | Genómica |

| ChromPos_SweetPotato_unkn_v4 | Genómica |

| AlnCnt_SweetPotato_unkn_v4 | Genómica |

| AlnEvalue_SweetPotato_unkn_v4 | Genómica |

| SNP | Genómica |

| SnpPosition | Genómica |

| CallRate | Calidad |

| OneRatioRef | Calidad |

| OneRatioSnp | Genómica |

| FreqHomRef | Calidad |

| FreqHomSnp | Genómica |

| FreqHets | Calidad |

| PICRef | Calidad |

| PICSnp | Genómica |

| AvgPIC | Calidad |

| AvgCountRef | Calidad |

| AvgCountSnp | Genómica |

| RepAvg | Calidad |

| rdepth | Calidad |

| maf | Calidad |

| locusName | ID |




**Global SilicoDArT geno** (1 metadatos, 5970 muestras)


| Columna | Categoría |
|---------|-----------|

| Sample_code | Otro |




**Global SilicoDArT metrics** (14 metadatos, 0 muestras)


| Columna | Categoría |
|---------|-----------|

| ﻿CloneID | ID |

| AlleleSequence | Secuencia |

| TrimmedSequence | Secuencia |

| Chrom_SweetPotato_unkn_v4 | Genómica |

| ChromPos_SweetPotato_unkn_v4 | Genómica |

| AlnCnt_SweetPotato_unkn_v4 | Genómica |

| AlnEvalue_SweetPotato_unkn_v4 | Genómica |

| CallRate | Calidad |

| OneRatio | Calidad |

| PIC | Calidad |

| AvgReadDepth | Calidad |

| StDevReadDepth | Calidad |

| Qpmr | Calidad |

| Reproducibility | Calidad |




**LowDensity SNP** (32 metadatos, 635 muestras)


| Columna | Categoría |
|---------|-----------|

| AlleleID | ID |

| CloneID | ID |

| ClusterTempIndex | Otro |

| AlleleSequence | Secuencia |

| ClusterConsensusSequence | Otro |

| ClusterSize | Otro |

| AlleleSeqDist | Otro |

| SNP | Genómica |

| SnpPosition | Genómica |

| CallRate | Calidad |

| OneRatioRef | Calidad |

| OneRatioSnp | Genómica |

| FreqHomRef | Calidad |

| FreqHomSnp | Genómica |

| FreqHets | Calidad |

| PICRef | Calidad |

| PICSnp | Genómica |

| AvgPIC | Calidad |

| AvgCountRef | Calidad |

| AvgCountSnp | Genómica |

| RatioAvgCountRefAvgCountSnp | Genómica |

| FreqHetsMinusFreqMinHom | Calidad |

| AlleleCountsCorrelation | Otro |

| aggregateTagsTotal | Otro |

| DerivedCorrMinusSeedCorr | Otro |

| RepRef | Otro |

| RepSNP | Genómica |

| RepAvg | Calidad |

| PicRepRef | Calidad |

| PicRepSNP | Genómica |

| TotalPicRepRefTest | Calidad |

| TotalPicRepSnpTest | Genómica |




**LowDensity SilicoDArT** (15 metadatos, 635 muestras)


| Columna | Categoría |
|---------|-----------|

| CloneID | ID |

| ClusterTempIndex | Otro |

| AlleleSequence | Secuencia |

| ClusterConsensusSequence | Otro |

| ClusterSize | Otro |

| CallRate | Calidad |

| OneRatio | Calidad |

| PIC | Calidad |

| AvgReadDepth | Calidad |

| StDevReadDepth | Calidad |

| Qpmr | Calidad |

| aggregateTagsTotal | Otro |

| Reproducibility | Calidad |

| PicRep | Calidad |

| TotalPicRepTest | Calidad |




**Wild SNP** (32 metadatos, 635 muestras)


| Columna | Categoría |
|---------|-----------|

| AlleleID | ID |

| CloneID | ID |

| ClusterTempIndex | Otro |

| AlleleSequence | Secuencia |

| ClusterConsensusSequence | Otro |

| ClusterSize | Otro |

| AlleleSeqDist | Otro |

| SNP | Genómica |

| SnpPosition | Genómica |

| CallRate | Calidad |

| OneRatioRef | Calidad |

| OneRatioSnp | Genómica |

| FreqHomRef | Calidad |

| FreqHomSnp | Genómica |

| FreqHets | Calidad |

| PICRef | Calidad |

| PICSnp | Genómica |

| AvgPIC | Calidad |

| AvgCountRef | Calidad |

| AvgCountSnp | Genómica |

| RatioAvgCountRefAvgCountSnp | Genómica |

| FreqHetsMinusFreqMinHom | Calidad |

| AlleleCountsCorrelation | Otro |

| aggregateTagsTotal | Otro |

| DerivedCorrMinusSeedCorr | Otro |

| RepRef | Otro |

| RepSNP | Genómica |

| RepAvg | Calidad |

| PicRepRef | Calidad |

| PicRepSNP | Genómica |

| TotalPicRepRefTest | Calidad |

| TotalPicRepSnpTest | Genómica |




**Wild SilicoDArT metrics** (15 metadatos, 0 muestras)


| Columna | Categoría |
|---------|-----------|

| ﻿CloneID | ID |

| AlleleSequence | Secuencia |

| TrimmedSequence | Secuencia |

| Chrom_SweetPotato_NSP323_v3 | Genómica |

| ChromPosTag_SweetPotato_NSP323_v3 | Genómica |

| AlnCnt_SweetPotato_NSP323_v3 | Genómica |

| AlnEvalue_SweetPotato_NSP323_v3 | Genómica |

| Strand_SweetPotato_NSP323_v3 | Genómica |

| CallRate | Calidad |

| OneRatio | Calidad |

| PIC | Calidad |

| AvgReadDepth | Calidad |

| StDevReadDepth | Calidad |

| Qpmr | Calidad |

| Reproducibility | Calidad |



### Tabla cruzada de columnas de calidad


| Columna                 | Global SNP geno   | Global SNP metrics   | Global SilicoDArT geno   | Global SilicoDArT metrics   | LowDensity SNP   | LowDensity SilicoDArT   | Wild SNP   | Wild SilicoDArT metrics   |
|:------------------------|:------------------|:---------------------|:-------------------------|:----------------------------|:-----------------|:------------------------|:-----------|:--------------------------|
| AvgCountRef             | —                 | ✅                   | —                        | —                           | ✅               | —                       | ✅         | —                         |
| AvgPIC                  | —                 | ✅                   | —                        | —                           | ✅               | —                       | ✅         | —                         |
| AvgReadDepth            | —                 | —                    | —                        | ✅                          | —                | ✅                      | —          | ✅                        |
| CallRate                | —                 | ✅                   | —                        | ✅                          | ✅               | ✅                      | ✅         | ✅                        |
| FreqHets                | —                 | ✅                   | —                        | —                           | ✅               | —                       | ✅         | —                         |
| FreqHetsMinusFreqMinHom | —                 | —                    | —                        | —                           | ✅               | —                       | ✅         | —                         |
| FreqHomRef              | —                 | ✅                   | —                        | —                           | ✅               | —                       | ✅         | —                         |
| OneRatio                | —                 | —                    | —                        | ✅                          | —                | ✅                      | —          | ✅                        |
| OneRatioRef             | —                 | ✅                   | —                        | —                           | ✅               | —                       | ✅         | —                         |
| PIC                     | —                 | —                    | —                        | ✅                          | —                | ✅                      | —          | ✅                        |
| PICRef                  | —                 | ✅                   | —                        | —                           | ✅               | —                       | ✅         | —                         |
| PicRep                  | —                 | —                    | —                        | —                           | —                | ✅                      | —          | —                         |
| PicRepRef               | —                 | —                    | —                        | —                           | ✅               | —                       | ✅         | —                         |
| Qpmr                    | —                 | —                    | —                        | ✅                          | —                | ✅                      | —          | ✅                        |
| RepAvg                  | —                 | ✅                   | —                        | —                           | ✅               | —                       | ✅         | —                         |
| Reproducibility         | —                 | —                    | —                        | ✅                          | —                | ✅                      | —          | ✅                        |
| StDevReadDepth          | —                 | —                    | —                        | ✅                          | —                | ✅                      | —          | ✅                        |
| TotalPicRepRefTest      | —                 | —                    | —                        | —                           | ✅               | —                       | ✅         | —                         |
| TotalPicRepTest         | —                 | —                    | —                        | —                           | —                | ✅                      | —          | —                         |
| maf                     | —                 | ✅                   | —                        | —                           | —                | —                       | —          | —                         |
| rdepth                  | —                 | ✅                   | —                        | —                           | —                | —                       | —          | —                         |


## 3. Verificación de Duplicados


- **LowDensity SNP:** `01_Report_DSp25-515_SNPs_Filtered_by _reads.csv` — 100,932,615 bytes — MD5: `ae014cb2a8a17cc215f73b545c571f8b`

- **Wild SNP:** `01_Report_DSp25-515_SNPs_Filtered_by _reads.csv` — 100,932,615 bytes — MD5: `ae014cb2a8a17cc215f73b545c571f8b`

- **Resultado:** ⚠️ IDÉNTICOS — mismo archivo en ambos datasets


## 4. Missingness Profunda


### Correlación CallRate ↔ Missingness observada


- **Global SNP:** CallRate mediana = 0.9781, 
missingness esperada mediana = 2.2%, 
marcadores con >20% miss = 551

- **Global SilicoDArT:** CallRate mediana = 0.9769, 
missingness esperada mediana = 2.3%, 
marcadores con >20% miss = 0

- **Wild SilicoDArT:** CallRate mediana = 0.9980, 
missingness esperada mediana = 0.2%, 
marcadores con >20% miss = 0



### Missingness por cromosoma


- **Global SNP:** 20,101 marcadores, 3,419 no mapeados (17%)

- **Global SilicoDArT:** 57,715 marcadores, 17,840 no mapeados (31%)

- **Wild SilicoDArT:** 236,607 marcadores, 143,437 no mapeados (61%)



### Outliers de missingness por muestra


- **Global SNP:** 5970 muestras, 475 outliers (>10.5% miss), 
max = 34.4%

- **Global SilicoDArT:** 5970 muestras, 237 outliers (>5.8% miss), 
max = 18.0%

- **LowDensity SNP:** 635 muestras, 23 outliers (>29.7% miss), 
max = 59.2%

- **LowDensity SilicoDArT:** 635 muestras, 15 outliers (>13.2% miss), 
max = 50.9%



## 5. Valores Centinela y Codificación


**Global SNP:** 
3 valores únicos. 
✅ Solo valores esperados


**Global SilicoDArT:** 
2 valores únicos. 
✅ Solo valores esperados


**LowDensity SNP:** 
2 valores únicos. 
✅ Solo valores esperados


**LowDensity SilicoDArT:** 
2 valores únicos. 
✅ Solo valores esperados


### Valores centinela en columnas genómicas


- **Global SNP:** 3,419 marcadores no mapeados (17.0%). 
ChromPos=0 y AlnEvalue=999 co-ocurren con Chrom=* al 100%

- **Global SilicoDArT:** 17,840 marcadores no mapeados (30.9%). 
ChromPos=0 y AlnEvalue=999 co-ocurren con Chrom=* al 100%

- **Wild SilicoDArT:** 143,437 marcadores no mapeados (60.6%). 
ChromPos=0 y AlnEvalue=999 co-ocurren con Chrom=* al 100%


