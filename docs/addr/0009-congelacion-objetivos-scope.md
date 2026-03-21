# ADR-0009: Congelación de Objetivo General y Específicos

**Estado:** Aceptada  
**Fecha:** 2026-03-06  
**Contexto:** El scope del proyecto se ha ido ampliando con cada experimento nuevo (AE, transformers, frontera de estabilidad). Es necesario congelar los objetivos para evitar scope-creep y garantizar que la narrativa del paper y el póster sean coherentes.

## Problema

Sin un objetivo general y específicos congelados, cada experimento nuevo desplaza el foco del proyecto. El póster (camera-ready) y el short paper deben contar la misma historia.

## Decisión


### Objetivo general

Evaluar si los mapas de diversidad genómica obtenidos con PCA→UMAP+kNN son estructuralmente robustos en paneles de camote con IDs disjuntos, bajo perturbaciones razonables del preprocesamiento y del muestreo de marcadores.

### Objetivos específicos

1. **OE1 — Geometría y calidad estructural:** Caracterizar la geometría y calidad estructural de cada panel sin correspondencia entre IDs, mediante diagnósticos automáticos (rango efectivo, reciprocidad kNN, componentes conexas, flags de QA).

2. **OE2 — Robustez del mapa:** Cuantificar la robustez del mapa frente a subsampling de marcadores, inyección de missingness (MCAR), y estrategia de imputación (moda vs. mediana), reportando estabilidad del subespacio PCA, topología kNN local y drift de distancias.

3. **OE3 — Alternativas no lineales:** Evaluar si alternativas no lineales (autoencoder denoising) ofrecen ventajas reales sobre PCA en régimen $n \ll p$, considerando tanto trustworthiness como estabilidad inter-seed.

## Hipótesis

H1 — Robustez estructural:
La geometría de diversidad genómica obtenida mediante PCA es estable bajo perturbaciones razonables del conjunto de marcadores y de los patrones de missingness.

H2 — Baseline representacional:
En matrices genómicas ultra-anchas ($n \ll p$), las representaciones basadas en PCA presentan mayor estabilidad topológica que representaciones aprendidas mediante autoencoders.

### Fuera de scope (explícitamente)

- Detección de comunidades / clustering formal.
- Core-set selection.
- Transformers / modelos pre-entrenados genómicos.
- Validación biológica con fenotipo.
- Integración multi-cultivo.
- Perturbaciones MNAR (reconocido como limitación, no como experimento).

## Consecuencias

- Todo experimento nuevo se evalúa contra estos tres OE. Si no responde a ninguno, no entra en el paper.
- El póster se estructura como: OE1 → bloques 1, 4; OE2 → bloques 2, 3; OE3 → bloques 5, 6, 7.
- El short paper se actualiza para reflejar estos objetivos en §1 y §6 (conclusiones).
- Trabajo futuro (Módulos 3, 4) se mantienen como extensiones documentadas pero no como scope activo.
