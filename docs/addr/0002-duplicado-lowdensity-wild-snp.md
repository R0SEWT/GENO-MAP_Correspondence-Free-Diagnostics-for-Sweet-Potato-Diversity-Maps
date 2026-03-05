# ADR-0002: Duplicado LowDensity SNP ↔ Wild SNP

**Estado:** Aceptada  
**Fecha:** 2026-03-01

## Contexto

Los datasets **LowDensity** (`10.21223/P3UBDJ44`) y **Wild** (`10.21223/P33VYY8C`)
fueron publicados por separado en el CIP Dataverse, pero ambos contienen un archivo
con nombre idéntico:

```
01_Report_DSp25-515_SNPs_Filtered_by _reads.csv
```

La orden DArT `DSp25-515` generó este reporte SNP que fue empaquetado en ambos
datasets. Tanto el tamaño (100,932,615 bytes), las dimensiones (62,736 × 667) como
el hash MD5 son idénticos:

| Propiedad | LowDensity | Wild |
|-----------|-----------|------|
| Tamaño | 100,932,615 bytes | 100,932,615 bytes |
| MD5 | `ae014cb2a8a17cc215f73b545c571f8b` | `ae014cb2a8a17cc215f73b545c571f8b` |
| Dimensiones | 62,736 × 667 | 62,736 × 667 |
| Resultado | **IDÉNTICOS** | **IDÉNTICOS** |

**Evidencia:** Verificación MD5 en `notebooks/eda_variables.ipynb`, celda 6 (Sección 3).

## Decisión

1. **Usar solo UNA copia** del archivo SNP `DSp25-515` en todo análisis downstream.
   Se conserva la del dataset **LowDensity** como referencia canónica, por ser el
   dataset que incluye también el SilicoDArT correspondiente.

2. **No combinar ambas copias** en análisis que integren múltiples datasets. Hacerlo
   duplicaría 62,736 marcadores × 635 muestras, inflando artificialmente el tamaño
   de la muestra.

3. **Documentar en `metadata.json`** la relación: agregar un campo `"shared_order"`
   o nota equivalente que vincule ambos datasets a la orden DArT `DSp25-515`.

4. El dataset **Wild** aporta dos archivos únicos que no están en LowDensity:
   - `02_SilicoDArT_metrics.csv` (métricas SilicoDArT, 236,607 marcadores)
   - `03_SilicoDArT_Data dictionary.xlsx` (diccionario de variables)

   Estos sí se usan normalmente desde el dataset Wild.

## Consecuencias

- **Pro:** Evita doble conteo de muestras y marcadores en análisis combinados
- **Pro:** Clarifica la procedencia real: una sola orden DArT, dos publicaciones
- **Con:** Requiere lógica explícita para excluir Wild SNP en pipelines automáticos
- **Mitigación:** El notebook `eda_variables.ipynb` ya implementa esta exclusión
  (`if is_duplicate: GENO_FILES` omite Wild SNP)
