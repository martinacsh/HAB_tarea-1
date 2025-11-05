# Proyecto de Análisis Funcional y Propagación en Redes

# Análisis funcional de **COX4I2**, **ND1** y **ATP6**

## Introducción

Este proyecto implementa un análisis funcional centrado en tres genes mitocondriales clave de la fosforilación oxidativa (OXPHOS):

- **COX4I2** — subunidad reguladora del **Complejo IV** (citocromo c oxidasa), ajusta la eficiencia respiratoria y la respuesta a oxígeno.  
- **ND1 (MT-ND1)** — subunidad esencial del **Complejo I** (NADH deshidrogenasa), punto de entrada de electrones a la cadena.  
- **ATP6 (MT-ATP6)** — subunidad del canal **F₀** de la **ATP sintasa** (Complejo V), necesaria para la síntesis de ATP.

**Qué haremos:** recuperar anotación gene-céntrica (GO/KEGG/Reactome) y realizar enriquecimiento funcional (**ORA**) para identificar procesos y vías sobre-representados en la lista de genes.

---

## Bases teóricas (resumen)

- **Gene Ontology (GO):** vocabulario controlado para describir **Procesos Biológicos (BP)**, **Funciones Moleculares (MF)** y **Componentes Celulares (CC)**.  
- **KEGG / Reactome:** bases de vías (*pathways*) que agrupan interacciones y reacciones moleculares.  
- **Análisis de sobre-representación (ORA):** contrasta si ciertos términos/vías aparecen más de lo esperado por azar en nuestra lista.  
  - Prueba **hipergeométrica**.  
  - Corrección por múltiples pruebas (**FDR**).

---

## Metodología

### 1) Entrada

- `data/genes_input.txt` con símbolos de genes (no se modifica).  
- Aunque nos enfocamos en **COX4I2**, **ND1** y **ATP6**, el flujo acepta **cualquier lista**.

### 2) Mapeo y anotación gene-céntrica

- Conversión de símbolos a **Entrez/Ensembl/UniProt** usando **MyGene.info**.  
- Recuperación, por gen, de términos **GO (BP/MF/CC)** y vías **KEGG/Reactome**.

### 3) Enriquecimiento (ORA)

- Motor **g:Profiler** contra **GO:BP/MF/CC**, **KEGG** y **Reactome**.  
- Cálculo de **p-valores** y **FDR (Benjamini–Hochberg)**.  
- Priorización de términos/vías más significativos.

### 4) Exportación y visualización

**CSV**
- `results/id_mapping.csv` — equivalencias de IDs por gen.  
- `results/functional_annotations.csv` — anotaciones GO/KEGG/Reactome por gen.  
- `results/enrichment_table.csv` — resultados de ORA.

**Figuras**
- `results/plots/*.png` — barplots con **−log10(FDR)** por base de datos.

**Resumen**
- `results/resumen.csv` — mapeos, anotaciones y *top* de enriquecimiento.

---

## Justificación de métodos (resumen)

- **Anotación gene-céntrica (MyGene.info):** unifica **Entrez/Ensembl/UniProt** y recupera **GO/KEGG/Reactome** por gen, asegurando **trazabilidad** y **desambiguación** (crítica en genes mitocondriales). Aporta **contexto** incluso con listas cortas.  
- **Inferencia funcional (ORA con g:Profiler):** adecuada para listas **no rankeadas** y **pequeñas**; prueba **hipergeométrica** y control **FDR (BH)**. **GSEA** no se usa por defecto (requiere ranking global y pierde potencia en listas pequeñas), pero se contempla si se dispone de un **universo rankeado**.  
- **Visualización:** barplots de **−log10(FDR)** por base (**GO: BP/MF/CC**, **KEGG**, **Reactome**) para comunicar significancia; la **tabla completa** de enriquecimiento (incluyendo genes de intersección) se mantiene en **CSV** para análisis detallado.  
- **Reproducibilidad/robustez:** CLI con `--all` (mapeo+anotación+ORA+gráficos+resumen), **valores por defecto seguros**, normalización de organismo (**MyGene→"human"**, **g:Profiler→"hsapiens"**), opción `--verbose`, **dependencias fijadas** y **reintentos** ante fallos puntuales del servicio de enriquecimiento.  
- **Extensiones y limitaciones:** posible **propagación en redes** (p. ej., **Random Walk with Restart** sobre PPI) para priorizar vecinos y repetir ORA capturando módulos. Con **listas muy pequeñas**, la **potencia** es limitada; se compensa con la **evidencia estable** de la anotación gene-céntrica.


---

### 5) Robustez y reproducibilidad

- Normalización de organismo automática: **human** (MyGene) / **hsapiens** (g:Profiler).  
- Reintento si g:Profiler falla (p. ej., **error 500**) y continuidad del *pipeline*.  
- Dependencias fijadas en `requirements.txt`.

---

## Resultados esperados (para COX4I2, ND1, ATP6)

- **GO:BP:** respiración mitocondrial, transporte de electrones, síntesis de ATP, mantenimiento del potencial de membrana, respuesta a hipoxia (por el papel regulador de **COX4I2**).  
- **KEGG/Reactome:** **Oxidative phosphorylation (OXPHOS)** y módulos de **Complejo I/IV/V**.

> **Nota:** con listas pequeñas, el poder estadístico del **ORA** es limitado; la anotación gene-céntrica aporta contexto clave.


