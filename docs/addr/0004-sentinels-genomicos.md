# ADR-0004: Valores centinela en columnas genómicas

**Estado:** Aceptada  
**Fecha:** 2026-03-01

## Contexto

Los archivos DArT usan valores centinela (sentinel values) para codificar condiciones
especiales. El EDA (`notebooks/eda_variables.ipynb`, Sección 5) identificó dos sistemas
de codificación que coexisten:

### 1. Sentinels en columnas genómicas (archivos de métricas)

Cuando un marcador **no mapea** al genoma de referencia, DArT usa:

| Columna | Valor sentinel | Significado |
|---------|---------------|-------------|
| `Chrom_*` | `*` | Cromosoma desconocido |
| `ChromPos_*` / `ChromPosTag_*` | `0` | Posición desconocida |
| `AlnEvalue_*` | `999` | Sin alineamiento |
| `AlnCnt_*` | `0` | Cero alineamientos |

**Co-ocurrencia verificada al 100%:** En los 3 archivos de métricas analizados,
`Chrom=*` siempre implica `ChromPos=0` y `AlnEvalue=999`.

| Dataset | Marcadores no mapeados | % del total |
|---------|----------------------|-------------|
| Global SNP | 3,419 | 17% |
| Global SilicoDArT | 17,840 | 31% |
| Wild SilicoDArT | 143,437 | **61%** |

La diferencia en % de no mapeados refleja dos genomas de referencia distintos:
- Global: `SweetPotato_unkn_v4`
- Wild: `SweetPotato_NSP323_v3` (menos completo → más marcadores sin mapear)

### 2. Sentinels en columnas de genotipo

| Fuente | Missing codificado como | Encontrado en |
|--------|------------------------|---------------|
| Archivos CSV | `"-"` (guión) o celda vacía | Todos los CSV |
| Data dictionary xlsx | `9` | Definición formal SilicoDArT |

La auditoría de valores en genotipos (celda 15) confirmó que los CSV solo contienen
`{0, 1}` (SilicoDArT) o `{0, 1, 2}` (SNP) + NaN. **No se encontró `9`** en ningún
CSV — el sentinel `9` del xlsx no se usa en la práctica.

## Decisión

1. **Para marcadores no mapeados (`Chrom=*`):**
   - Mantenerlos en el dataset, no filtrar automáticamente
   - Crear columna derivada booleana `is_mapped` para facilitar filtrado selectivo
   - En análisis que requieran posición genómica (ej: LD decay, GWAS), excluirlos
   - En análisis de diversidad genética (PCA, distancias), incluirlos

2. **Para valores numéricos sentinel (`ChromPos=0`, `AlnEvalue=999`):**
   - Reemplazar por `NaN` en el pipeline de preprocesamiento
   - Nunca tratarlos como valores numéricos reales (ej: no promediar `AlnEvalue`
     incluyendo 999s)

3. **Para missing en genotipos:**
   - Estandarizar a `NaN` en pandas desde la carga: `"-"`, `""`, `" "` → `NaN`
   - Ignorar el sentinel `9` del xlsx (no está en los datos reales)
   - Documentar que `load_genotypes_subsample()` del notebook ya hace esta
     transformación

4. **Nomenclatura de cromosomas:**
   - Global usa `chr1`–`chr15` (lowercase, sin padding)
   - Wild usa `Chr01`–`Chr15` (capitalized, zero-padded)
   - Normalizar a formato lowercase sin padding (`chr1`) en el pipeline

## Consecuencias

- **Pro:** Pipeline robusto ante valores que parecen numéricos pero son sentinels
- **Pro:** Consistencia entre datasets con distintas convenciones DArT
- **Pro:** Los marcadores no mapeados (hasta 61%) se preservan para análisis de
  diversidad donde la posición genómica no es necesaria
- **Con:** Requiere paso explícito de normalización al cargar cada archivo
- **Mitigación:** Centralizar la lógica de carga en una función reutilizable
  (candidata: `scripts/load_dart.py` o similar)
