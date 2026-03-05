# 🧬 DATA_CATALOG — OXOR / GENO-MAP

> **Generado automáticamente** por `notebooks/data_catalog.ipynb`  
> **Fecha:** 2026-03-01 00:55  
> **Fuente de datos:** [CIP Dataverse](https://data.cipotato.org)

---

## TL;DR

- **Universo:** 100 datasets escaneados del CIP Dataverse (últimos 100 publicados)
- **Catálogo curado:** 11 datasets de genotipado de 5 cultivos andinos
- **Descargados localmente:** 3 datasets (1210 MB total)
- **Restringidos:** 7 datasets (publicados pero requieren permisos CIP)
- **Accesibles sin descargar:** 1
- **Genomas de referencia:** SweetPotato_unkn_v4, SweetPotato_NSP323_v3
- ⚠️ **Múltiples genomas de referencia**: integración cruzada requiere re-mapeo

---

## 1. Universo CIP Dataverse
- **Total datasets escaneados:** 100
- **Relacionados con genotipado:** 20
- **Seleccionados para este proyecto:** 11

### Datasets de genotipado detectados

| doi_code   | title                                                                                                                                                                                                          | fecha      | en_catalogo   |
|:-----------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-----------|:--------------|
| nan        | Dataset for: At least 20 advanced clones from LBHTxLTVR population with resistance to late blight, virus, heat and drought tolerance, high tuber yield and quality for french fries and/or chips and earliness | 2026-01-28 | False         |
| FFORMJ     | Dataset for: Advanced clones of group B3 cycle 3, population B with late blight resistance in Oxapampa, 2014                                                                                                   | 2026-01-28 | False         |
| OTYRIV     | Dataset for: Advanced clones of group B3 cycle 3, population B with late blight resistance in Oxapampa, 2012                                                                                                   | 2026-01-28 | False         |
| LCGD4P     | Dataset for: Advanced clones of group B3 cycle 3, population B with late blight resistance in Oxapampa, 2013                                                                                                   | 2026-01-28 | False         |
| 3VYY8C     | Replication data for: Genetic Diversity Analysis of a large collection of sweetpotato wild relatives, Ipomoea spp.                                                                                             | 2026-01-28 | True          |
| UBDJ44     | Dataset for: Genotyping of sweetpotato genebank accessions using low-density DArTSeq markers                                                                                                                   | 2026-01-28 | True          |
| SFXXDC     | Dataset for: Advanced clones of  B3C0, group B3, population B in Oxapampa-Peru                                                                                                                                 | 2026-01-27 | False         |
| TNNRYW     | Dataset for: Advanced clones of  B3C0, group B3, population B in Comas-Peru                                                                                                                                    | 2026-01-27 | False         |
| JW6UOC     | Dataset for: High density genotyping of cultivated yacon genebank accessions using SNP markers from DArTseq                                                                                                    | 2026-01-27 | True          |
| 5KHIK4     | Dataset for: High density genotyping of cultivated yacon genebank accessions using SilicoDArT markers                                                                                                          | 2026-01-27 | True          |
| CD01KI     | Dataset for: High density genotyping of cultivated ullucus genebank accessions using SNP markers from DArTseq                                                                                                  | 2026-01-27 | True          |
| IXXW5N     | Dataset for: High density genotyping of cultivated ullucus genebank accessions using SilicoDArT markers                                                                                                        | 2026-01-27 | True          |
| FYYLLO     | Dataset for: Genotypic Dataset of In Situ Native Ulluco Varieties from the Andenes de Cuyocuyo Agrobiodiversity Zone using SilicoDArT markers                                                                  | 2026-01-27 | True          |
| M613LS     | Dataset for: Genotypic Dataset of In Situ Native Oca Varieties from the Andenes de Cuyocuyo Agrobiodiversity Zone using SilicoDArT markers                                                                     | 2026-01-27 | True          |
| TFOJNZ     | Dataset for: Genotypic Dataset of In Situ Native Mashua Varieties from the Andenes de Cuyocuyo Agrobiodiversity Zone using SilicoDArT markers                                                                  | 2026-01-27 | True          |
| EKOYZ1     | Dataset for: Heritability for yield components in LBHT potato clones, under warm conditions                                                                                                                    | 2026-01-27 | False         |
| D3XEEH     | Dataset for: Genotyping of sweetpotato genebank accessions using amplicon-based targeted markers (DArTag)                                                                                                      | 2025-12-31 | True          |
| RUTCZH     | Dataset for: Characterization database of 36 genotypes of native potatoes of the andigenum group from the Cumbal Indigenous Reservation - 2025                                                                 | 2025-12-31 | False         |
| HVXVLB     | Dataset for: Characterization database of 64 genotypes of native potatoes of the Phurejas group from the Cumbal Indigenous Reservation - 2025                                                                  | 2025-12-17 | False         |
| S2IMOS     | Replication data for: Revealing the genetic diversity in the global ex situ sweetpotato collection with next generation sequencing technologies.                                                               | 2025-12-17 | True          |

## 2. Catálogo Curado (metadata.json)

**11 datasets** de **5 cultivos**: Sweetpotato, Yacon, Ulluco, Oca, Mashua

| Status | Cantidad |
|--------|----------|
| 🔒 Restringido | 7 |
| ✅ Descargado | 3 |
| 🌐 Accesible (no descargado) | 1 |

| Cultivo     | Tipo                             | DOI                                | Status                       |
|:------------|:---------------------------------|:-----------------------------------|:-----------------------------|
| Sweetpotato | low-density DArTSeq              | https://doi.org/10.21223/P3/UBDJ44 | ✅ Descargado                |
| Sweetpotato | DArTag 3.1K amplicones           | https://doi.org/10.21223/P3/D3XEEH | 🔒 Restringido               |
| Sweetpotato | diversidad global NGS            | https://doi.org/10.21223/P3/S2IMOS | ✅ Descargado                |
| Sweetpotato | wild relatives                   | https://doi.org/10.21223/P3/3VYY8C | ✅ Descargado                |
| Yacon       | SNP                              | https://doi.org/10.21223/P3/JW6UOC | 🌐 Accesible (no descargado) |
| Yacon       | SilicoDArT                       | https://doi.org/10.21223/P3/5KHIK4 | 🔒 Restringido               |
| Ulluco      | 300 accesiones SNP               | https://doi.org/10.21223/P3/CD01KI | 🔒 Restringido               |
| Ulluco      | 300 accesiones SilicoDArT        | https://doi.org/10.21223/P3/IXXW5N | 🔒 Restringido               |
| Ulluco      | 68 accesiones in situ SilicoDArT | https://doi.org/10.21223/P3/FYYLLO | 🔒 Restringido               |
| Oca         | 56 accesiones in situ SilicoDArT | https://doi.org/10.21223/P3/M613LS | 🔒 Restringido               |
| Mashua      | 28 accesiones in situ SilicoDArT | https://doi.org/10.21223/P3/TFOJNZ | 🔒 Restringido               |

## 3. Inventario de Archivos Locales

**9 archivos**, **1210 MB** (1.2 GB)

| Cultivo     | Dataset               | Archivo                                               |   Tamaño (MB) |
|:------------|:----------------------|:------------------------------------------------------|--------------:|
| Sweetpotato | low-density DArTSeq   | 01_Report_DSp25-515_SNPs_Filtered_by _reads.csv       |          96.3 |
| Sweetpotato | low-density DArTSeq   | 02_Report_DSp25-515_Silico-DArT_Filtered_by_reads.csv |          54   |
| Sweetpotato | diversidad global NGS | SNP_Genotypes.csv                                     |         229.5 |
| Sweetpotato | diversidad global NGS | SNP_metrics.csv                                       |           6.8 |
| Sweetpotato | diversidad global NGS | SilicoDArT_Genotypes.csv                              |         657.8 |
| Sweetpotato | diversidad global NGS | SilicoDArT_metrics.csv                                |          12.3 |
| Sweetpotato | wild relatives        | 01_Report_DSp25-515_SNPs_Filtered_by _reads.csv       |          96.3 |
| Sweetpotato | wild relatives        | 02_SilicoDArT_metrics.csv                             |          49.5 |
| Sweetpotato | wild relatives        | 03_SilicoDArT_Data dictionary.xlsx                    |           7.7 |

## 4. Dimensiones y Formato

| Archivo                                               | Dataset               | Filas   | Columnas   | Delimitador   | Formato                   | Muestras (est.)   |   Tamaño (MB) |
|:------------------------------------------------------|:----------------------|:--------|:-----------|:--------------|:--------------------------|:------------------|--------------:|
| 01_Report_DSp25-515_SNPs_Filtered_by _reads.csv       | low-density DArTSeq   | 62,736  | 667        | ','           | marker_metrics + muestras | 635               |          96.3 |
| 02_Report_DSp25-515_Silico-DArT_Filtered_by_reads.csv | low-density DArTSeq   | 38,272  | 650        | ','           | marker_metrics + muestras | 635               |          54   |
| SNP_Genotypes.csv                                     | diversidad global NGS | 20,103  | 5,971      | ','           | marker_metrics + muestras | 5970              |         229.5 |
| SNP_metrics.csv                                       | diversidad global NGS | 20,101  | 25         | ','           | marker_metrics            | —                 |           6.8 |
| SilicoDArT_Genotypes.csv                              | diversidad global NGS | 57,717  | 5,971      | ','           | marker_metrics + muestras | 5970              |         657.8 |
| SilicoDArT_metrics.csv                                | diversidad global NGS | 57,715  | 14         | ','           | marker_metrics            | —                 |          12.3 |
| 01_Report_DSp25-515_SNPs_Filtered_by _reads.csv       | wild relatives        | 62,736  | 667        | ','           | marker_metrics + muestras | 635               |          96.3 |
| 02_SilicoDArT_metrics.csv                             | wild relatives        | 236,607 | 15         | ';'           | marker_metrics            | —                 |          49.5 |

## 5. Genomas de Referencia

- **SweetPotato_unkn_v4**: diversidad global NGS
- **SweetPotato_NSP323_v3**: wild relatives

> ⚠️ Múltiples genomas de referencia: integración cruzada requiere re-mapeo.

## 6. Métricas de Calidad por Marcador

| Dataset           | N marcadores   |   CallRate (media) |   CallRate (mediana) |   CallRate (min) |   PIC (media) |   PIC (mediana) |   PIC (min) |   Reproducibility (media) |   Reproducibility (mediana) |   Reproducibility (min) |   AvgReadDepth (media) |   AvgReadDepth (mediana) |   AvgReadDepth (min) |
|:------------------|:---------------|-------------------:|---------------------:|-----------------:|--------------:|----------------:|------------:|--------------------------:|----------------------------:|------------------------:|-----------------------:|-------------------------:|---------------------:|
| Global SNP        | 20,101         |             0.9448 |               0.9781 |           0.2136 |        0.2297 |          0.2398 |      0.0004 |                    0.9761 |                      0.9844 |                  0.8605 |                29.9346 |                  21.4    |               2.5    |
| Global SilicoDArT | 57,715         |             0.9679 |               0.9769 |           0.806  |        0.3091 |          0.3446 |      0      |                    0.9902 |                      0.9945 |                  0.95   |                18.9498 |                  14.3524 |               5      |
| Wild SilicoDArT   | 236,607        |             0.9926 |               0.998  |           0.8112 |        0.0735 |          0.0252 |      0.002  |                    0.9993 |                      1      |                  0.95   |                39.2332 |                  32.6786 |               5.0061 |

## 7. Missingness en Genotipos (estimación por muestreo)

| Dataset                                                                    |   N muestras (cols) |   Filas leídas |   Missingness (%) |
|:---------------------------------------------------------------------------|--------------------:|---------------:|------------------:|
| Sweetpotato low-density DArTSeq — 01_Report_DSp25-515_SNPs_Filtered_by _re |                 635 |            500 |              17.6 |
| Sweetpotato low-density DArTSeq — 02_Report_DSp25-515_Silico-DArT_Filtered |                 635 |            500 |               9.1 |
| Sweetpotato diversidad global NGS — SNP_Genotypes.csv                      |                5970 |            500 |               8   |
| Sweetpotato diversidad global NGS — SilicoDArT_Genotypes.csv               |                5970 |            500 |               2.4 |
| Sweetpotato wild relatives — 01_Report_DSp25-515_SNPs_Filtered_by _re      |                 635 |            500 |              17.6 |

## 8. Tabla Resumen Consolidada

| Cultivo     | Tipo                             | DOI    | Status         | Archivos   | Tamaño (MB)   | Genoma Ref            |
|:------------|:---------------------------------|:-------|:---------------|:-----------|:--------------|:----------------------|
| Sweetpotato | low-density DArTSeq              | UBDJ44 | ✅ Descargado  | 2          | 150.2         | —                     |
| Sweetpotato | DArTag 3.1K amplicones           | D3XEEH | 🔒 Restringido | —          | —             | —                     |
| Sweetpotato | diversidad global NGS            | S2IMOS | ✅ Descargado  | 4          | 906.4         | SweetPotato_unkn_v4   |
| Sweetpotato | wild relatives                   | 3VYY8C | ✅ Descargado  | 3          | 153.4         | SweetPotato_NSP323_v3 |
| Yacon       | SNP                              | JW6UOC | 🌐 Accesible   | —          | —             | —                     |
| Yacon       | SilicoDArT                       | 5KHIK4 | 🔒 Restringido | —          | —             | —                     |
| Ulluco      | 300 accesiones SNP               | CD01KI | 🔒 Restringido | —          | —             | —                     |
| Ulluco      | 300 accesiones SilicoDArT        | IXXW5N | 🔒 Restringido | —          | —             | —                     |
| Ulluco      | 68 accesiones in situ SilicoDArT | FYYLLO | 🔒 Restringido | —          | —             | —                     |
| Oca         | 56 accesiones in situ SilicoDArT | M613LS | 🔒 Restringido | —          | —             | —                     |
| Mashua      | 28 accesiones in situ SilicoDArT | TFOJNZ | 🔒 Restringido | —          | —             | —                     |

## 9. Nota sobre Acceso a Datasets

### Datasets restringidos (requieren permisos CIP)

- **Sweetpotato** — DArTag 3.1K amplicones (`D3XEEH`)
- **Yacon** — SilicoDArT (`5KHIK4`)
- **Ulluco** — 300 accesiones SNP (`CD01KI`)
- **Ulluco** — 300 accesiones SilicoDArT (`IXXW5N`)
- **Ulluco** — 68 accesiones in situ SilicoDArT (`FYYLLO`)
- **Oca** — 56 accesiones in situ SilicoDArT (`M613LS`)
- **Mashua** — 28 accesiones in situ SilicoDArT (`TFOJNZ`)

### Datasets accesibles (descarga manual)

- **Yacon** — SNP (`JW6UOC`)

> **⚠️ La API del CIP Dataverse** puede devolver HTTP 403 para archivos grandes incluso en
> datasets públicos. Esto es una restricción de red/tamaño, no de acceso. Descargar manualmente
> desde https://data.cipotato.org. Ver `docs/addr/0001-descarga-cip-dataverse.md`.

---

## 10. Estructura de Carpetas

```
data/
├── cipotato_datasets_latest100.json    # Universo: 100 datasets CIP Dataverse
├── metadata.json                        # Catálogo curado: 11 datasets seleccionados
├── DATA_CATALOG.md                      # Este archivo (generado por notebook)
├── 10.21223P3UBDJ44_LowDensity/
│   ├── 01_Report_DSp25-515_SNPs_Filtered_by _reads.csv (96 MB)
│   ├── 02_Report_DSp25-515_Silico-DArT_Filtered_by_reads.csv (54 MB)
├── 10.21223P30BVZYY_Genetic_diversity/
│   ├── SNP_Genotypes.csv (230 MB)
│   ├── SNP_metrics.csv (7 MB)
│   ├── SilicoDArT_Genotypes.csv (658 MB)
│   ├── SilicoDArT_metrics.csv (12 MB)
├── 10.21223P33VYY8C_Wild/
│   ├── 01_Report_DSp25-515_SNPs_Filtered_by _reads.csv (96 MB)
│   ├── 02_SilicoDArT_metrics.csv (50 MB)
│   ├── 03_SilicoDArT_Data dictionary.xlsx (8 MB)
└── outputs/
    └── *.json, *.png                    # Artefactos generados
```
