# Architecture Decision Records (ADRs)

Registro de decisiones técnicas del proyecto Oxor / GENO-MAP. TLDR 

## Formato

Cada ADR sigue esta estructura:

```markdown
# ADR-NNNN: Título descriptivo

**Estado:** Aceptada | Propuesta | Rechazada | Reemplazada  
**Fecha:** YYYY-MM-DD

## Contexto
¿Qué problema o situación motivó esta decisión?

## Decisión
¿Qué se decidió hacer?

## Consecuencias
¿Qué implica esta decisión? Pros, contras, riesgos.
```

## Convención de nombres

`NNNN-titulo-kebab-case.md` — ejemplo: `0001-descarga-cip-dataverse.md`

## ADRs existentes

| # | Título | Estado |
|---|--------|--------|
| 0001 | [Descarga CIP Dataverse](0001-descarga-cip-dataverse.md) | Aceptada |
| 0002 | [Duplicado LowDensity ↔ Wild SNP](0002-duplicado-lowdensity-wild-snp.md) | Aceptada |
| 0003 | [Missingness Estructurada (MNAR)](0003-missingness-estructurada-mnar.md) | Aceptada |
| 0004 | [Valores Centinela Genómicos](0004-sentinels-genomicos.md) | Aceptada |
| 0005 | [Validación Visual Baseline](0005-validacion-visual-baseline.md) | Aceptada |
| 0006 | [Comparación AE Level 1 vs Baseline](0006-comparacion-ae-vs-baseline.md) | Aceptada |
| 0007 | [Cierre Level 1 — Autoencoder Genómico](0007-cierre-level1-autoencoder.md) | Aceptada |
| 0008 | [Level A — Frontera de Estabilidad Representacional](0008-level-a-frontera-estabilidad.md) | Aceptada |
| 0009 | [Congelación de Objetivos y Scope](0009-congelacion-objetivos-scope.md) | Propuesta |
| 0010 | [Unificación de Métricas de Estabilidad](0010-unificacion-metricas-estabilidad.md) | Propuesta |
| 0011 | [Protocolo de Robustness Curves](0011-protocolo-robustness-curves.md) | Propuesta |
| 0012 | [Sensibilidad al Número de PCs / Varianza](0012-sensibilidad-pcs-varianza.md) | Aceptada |
| 0013 | [Block Removal vs. Submuestreo Aleatorio](0013-block-removal-experiment.md) | Aceptada |
| 0014 | [Alcance Sensibilidad UMAP](0014-umap-sensitivity-scope.md) | Aceptada |
| 0015 | [Cierre Benchmark AE](0015-ae-closure.md) | Aceptada |
