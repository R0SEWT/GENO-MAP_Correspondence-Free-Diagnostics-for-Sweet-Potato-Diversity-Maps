# ADR-0001: Descarga de archivos del CIP Dataverse

**Estado:** Aceptada  
**Fecha:** 2026-03-01

## Contexto

El CIP Dataverse (https://data.cipotato.org) publica datasets de genotipado como públicos
(DOI accesible, metadatos visibles vía API). Sin embargo, al intentar descargar archivos
programáticamente vía la API REST (`/api/access/datafile/{id}`), el servidor devuelve
**HTTP 403 Forbidden**.

Esto ocurre incluso cuando:
- El dataset no tiene flag `restricted` en sus metadatos API
- El DOI es público y accesible desde cualquier navegador
- Los metadatos del dataset (nombres de archivo, tamaños) se obtienen sin problema

Los archivos de genotipado DArT/DArTSeq son grandes (1–7 MB por archivo individual,
pero las matrices completas pueden superar los 200 MB).

## Decisión

1. **No marcar datasets como `restricted` en `metadata.json`** cuando la API devuelve 403.
   El 403 es una restricción de red/API, no de acceso al dataset.

2. **Descargar manualmente** desde la interfaz web del CIP Dataverse usando el DOI directo.
   Los archivos son descargables vía navegador sin autenticación.

3. **El notebook `data_catalog.ipynb`** intenta la descarga vía API como primer paso.
   Si falla con 403, informa al usuario sin modificar `metadata.json`.

4. Los datasets genuinamente restringidos (que requieren permisos del CIP) se marcan
   explícitamente como `"restricted": true` en `metadata.json` — esta información viene
   de la documentación del Dataverse, no de errores HTTP.

## Consecuencias

- **Pro:** Se evita confundir restricciones de red con restricciones de acceso
- **Pro:** `metadata.json` refleja el estado real de acceso, no errores transitorios
- **Con:** Requiere paso manual de descarga para datasets nuevos
- **Mitigación:** El notebook lista URLs directas y muestra instrucciones claras
