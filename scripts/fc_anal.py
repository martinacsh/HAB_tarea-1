#!/usr/bin/env python3
"""
func_an.py — Análisis funcional de genes (COX4I2, ND1, ATP6)

Descripción general
-------------------
Script de línea de comandos (CLI) para:
  1) Leer una lista de genes (símbolos) desde un archivo de texto.
  2) Mapear identificadores (Entrez, Ensembl, UniProt) y recuperar anotaciones básicas (GO, KEGG, Reactome) usando MyGene.info.
  3) (Opcional) Ejecutar un análisis de sobre-representación (ORA) mediante g:Profiler
     sobre GO (BP/MF/CC), KEGG y Reactome con corrección FDR (Benjamini–Hochberg).
  4) Exportar resultados a CSV y, si se solicita, generar gráficos sencillos.

Bases de datos y librerías
--------------------------
- MyGene.info (vía librería `mygene`): servicio para convertir IDs y recuperar
  anotaciones gene-céntricas (GO/KEGG/Reactome por gen).
- g:Profiler (librería `gprofiler-official`): enriquecimiento funcional de conjuntos
  de genes contra GO, KEGG y Reactome. Se usa el organismo "hsapiens" por defecto.
- pandas/numpy: manipulación de datos y tablas.
- matplotlib: generación de gráficos (barplots) de términos enriquecidos.
- statsmodels: NO es estrictamente necesario si usas g:Profiler (que ya aplica FDR),
  pero útil si extiendes el script a métodos propios.

Notas metodológicas
-------------------
- ORA (Over-Representation Analysis): prueba de hipergeométrica que contrasta la
  proporción de "hits" de tu lista dentro de un conjunto (p.ej. un término GO)
  frente a lo esperado por azar, aplicando corrección por múltiples pruebas (FDR).
- Tamaño de lista pequeño: con 3 genes (COX4I2, ND1, ATP6) el poder estadístico es
  limitado. Integra el resultado con anotación gene-céntrica.

Uso
----
python scripts/func_an.py \
  --input data/genes_input.txt \
  --outdir results/ \
  --organism hsapiens \
  --run-enrichment \
  --plot

Argumentos mínimos: --input y --outdir/--output. No modifica el archivo de entrada.

Salida (por defecto)
--------------------
- <outdir>/id_mapping.csv               — mapeo de IDs y nombres canónicos
- <outdir>/functional_annotations.csv   — anotaciones por gen (GO/KEGG/Reactome)
- <outdir>/enrichment_table.csv         — (si --run-enrichment) tabla ORA g:Profiler
- <outdir>/plots/*.png                  — (si --plot) barplots de términos top

Requisitos (requirements.txt sugerido)
--------------------------------------
- pandas, numpy, mygene, gprofiler-official, matplotlib

"""
from __future__ import annotations
import argparse
import sys
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional

import pandas as pd

# Importes perezosos (se cargan dentro de funciones si están disponibles)
# mygene y gprofiler-official no son parte de la stdlib


# ----------------------------- Utilidades básicas -----------------------------

def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="[%(asctime)s] %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def read_gene_list(path: Path) -> List[str]:
    """
    Lee símbolos génicos desde un archivo de texto.
    - Acepta formatos con uno por línea o separados por comas/espacios.
    - Ignora líneas vacías y comentarios que comiencen con '#'.
    - Devuelve una lista única (sin duplicados), en mayúsculas (símbolos estándar).
    """
    raw = []
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        # Permite coma/espacio como separador adicional
        parts = [p for chunk in line.split(',') for p in chunk.split()] if ',' in line or ' ' in line else [line]
        raw.extend(parts)
    # normaliza
    genes = [g.strip().upper() for g in raw if g.strip()]
    # dedup conservando orden
    seen = set()
    unique = []
    for g in genes:
        if g not in seen:
            unique.append(g)
            seen.add(g)
    if not unique:
        raise ValueError("La lista de genes está vacía tras el preprocesado.")
    logging.info("Genes de entrada (%d): %s", len(unique), ", ".join(unique))
    return unique


# ----------------------------- Mapeo con MyGene -------------------------------

def map_and_annotate_with_mygene(genes: List[str], organism: str = "hsapiens") -> Dict[str, pd.DataFrame]:
    """
    Usa MyGene.info para:
      - Mapear símbolos a IDs (Entrez/Ensembl/UniProt) y nombre canónico
      - Recuperar anotaciones GO (BP/MF/CC) y pathways (KEGG/Reactome) por gen

    Devuelve un diccionario con dos DataFrames:
      - mapping: una fila por símbolo de entrada (con columnas de IDs y nombre)
      - annotations: filas (gene, source, term_id, term_name)
    """
    try:
        import mygene  # type: ignore
    except ImportError as e:
        raise ImportError("Se requiere la librería 'mygene'. Añádela a requirements.txt e instálala.") from e

    mg = mygene.MyGeneInfo()

    # Campos a solicitar (ver documentación de mygene)
    fields = [
        "entrezgene",
        "ensembl.gene",
        "uniprot.Swiss-Prot",
        "name",
        "go.BP",
        "go.MF",
        "go.CC",
        "pathway.kegg",
        "pathway.reactome",
    ]

    logging.info("Consultando MyGene.info para mapeo y anotación…")
    qres = mg.querymany(
        genes,
        scopes="symbol",
        fields=",".join(fields),
        species=organism,
        as_dataframe=True,
        returnall=False,
        verbose=False,
    )

    # qres es un DataFrame indexado por símbolo consultado (si as_dataframe=True)
    # Normalizamos a un DF 'mapping' limpio y otro DF 'annotations'
    mapping_rows = []
    annot_rows = []

    def _ensure_list(x: Any) -> List[Any]:
        if x is None or (isinstance(x, float) and pd.isna(x)):
            return []
        if isinstance(x, list):
            return x
        return [x]

    # Iteración por símbolo de entrada
    for query_symbol, row in qres.iterrows():
        # Algunos símbolos pueden no resolverse; protegemos accesos
        entrez = row.get("entrezgene", None)
        ensembl = None
        if isinstance(row.get("ensembl.gene", None), list):
            ensembl = row.get("ensembl.gene")[0]
        else:
            ensembl = row.get("ensembl.gene", None)
        uniprot = None
        if isinstance(row.get("uniprot.Swiss-Prot", None), list):
            uniprot = row.get("uniprot.Swiss-Prot")[0]
        else:
            uniprot = row.get("uniprot.Swiss-Prot", None)
        name = row.get("name", None)

        mapping_rows.append({
            "query_symbol": query_symbol,
            "symbol": row.get("symbol", query_symbol),
            "name": name,
            "entrezgene": entrez,
            "ensembl_gene": ensembl,
            "uniprot_swissprot": uniprot,
        })

        # GO terms (cada categoría puede venir como lista de dicts con id/name)
        for cat in ("BP", "MF", "CC"):
            go_list = row.get(f"go.{cat}")
            for term in _ensure_list(go_list):
                term_id = term.get("id") if isinstance(term, dict) else None
                term_name = term.get("term") if isinstance(term, dict) else str(term)
                if term_id or term_name:
                    annot_rows.append({
                        "gene": row.get("symbol", query_symbol),
                        "source": f"GO:{cat}",
                        "term_id": term_id,
                        "term_name": term_name,
                    })
        # KEGG
        for term in _ensure_list(row.get("pathway.kegg")):
            term_id = term.get("id") if isinstance(term, dict) else None
            term_name = term.get("name") if isinstance(term, dict) else str(term)
            if term_id or term_name:
                annot_rows.append({
                    "gene": row.get("symbol", query_symbol),
                    "source": "KEGG",
                    "term_id": term_id,
                    "term_name": term_name,
                })
        # Reactome
        for term in _ensure_list(row.get("pathway.reactome")):
            term_id = term.get("id") if isinstance(term, dict) else None
            term_name = term.get("name") if isinstance(term, dict) else str(term)
            if term_id or term_name:
                annot_rows.append({
                    "gene": row.get("symbol", query_symbol),
                    "source": "REAC",
                    "term_id": term_id,
                    "term_name": term_name,
                })

    mapping_df = pd.DataFrame(mapping_rows)
    annotations_df = pd.DataFrame(annot_rows)

    return {"mapping": mapping_df, "annotations": annotations_df}


# ------------------------- Enriquecimiento con g:Profiler ---------------------

def run_gprofiler_enrichment(genes: List[str], organism: str = "hsapiens",
                              sources: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Ejecuta ORA usando g:Profiler y devuelve un DataFrame con los resultados.

    Parámetros
    ----------
    genes: lista de símbolos (p.ej., ["COX4I2", "ND1", "ATP6"]) 
    organism: código de organismo para g:Profiler (por defecto 'hsapiens')
    sources: lista de fuentes g:Profiler (ej.: ['GO:BP','GO:MF','GO:CC','KEGG','REAC'])

    Devuelve
    --------
    DataFrame con columnas típicas de g:Profiler, incluyendo p_value y p_value_adjusted.
    """
    if sources is None:
        sources = ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"]

    try:
        from gprofiler import GProfiler  # type: ignore
    except ImportError as e:
        raise ImportError("Se requiere 'gprofiler-official'. Añádelo a requirements.txt e instálalo.") from e

    logging.info("Ejecutando enriquecimiento con g:Profiler…")
    gp = GProfiler(return_dataframe=True)
    # g:Profiler espera símbolos válidos y organismo en formato del servicio (hsapiens)
    res_df = gp.profile(
        organism=organism,
        query=genes,
        sources=sources,
        user_threshold=1.0,  # dejamos que devuelva todo, filtraremos nosotros por FDR si se desea
    )
    # Normalización de nombres de columnas por conveniencia
    rename_map = {
        "p_value": "pvalue",
        "p_value_adjusted": "fdr",
        "source": "db",
        "term_id": "term_id",
        "name": "term_name",
        "term_size": "term_size",
        "query_size": "query_size",
        "intersection_size": "overlap",
        "intersections": "intersection_genes",
    }
    res_df = res_df.rename(columns=rename_map)

    # Orden sugerido: por FDR y luego por p-value
    sort_cols = [c for c in ("fdr", "pvalue") if c in res_df.columns]
    if sort_cols:
        res_df = res_df.sort_values(by=sort_cols, ascending=True)

    return res_df


# ---------------------------------- Plots -------------------------------------

def plot_top_terms(res_df: pd.DataFrame, outdir: Path, topn: int = 10) -> None:
    """Genera barplots de −log10(FDR) para las distintas fuentes (db)."""
    import math
    import matplotlib.pyplot as plt

    if res_df.empty:
        logging.warning("No hay resultados de enriquecimiento para graficar.")
        return

    out_plots = outdir / "plots"
    out_plots.mkdir(parents=True, exist_ok=True)

    # Filtra términos con FDR disponible
    df = res_df.copy()
    if "fdr" not in df.columns:
        logging.warning("La columna 'fdr' no está presente; los gráficos usarán p-value.")
        df["fdr"] = df.get("pvalue", 1.0)

    df["neglog10_fdr"] = df["fdr"].apply(lambda x: -math.log10(x) if x > 0 else 0)

    # Por cada base de datos, guardamos un barplot con top N términos
    for db in sorted(df["db"].unique()):
        sub = df[df["db"] == db].head(topn)
        if sub.empty:
            continue
        plt.figure(figsize=(10, max(3, 0.4 * len(sub))))
        plt.barh(sub["term_name"].astype(str), sub["neglog10_fdr"].astype(float))
        plt.xlabel("-log10(FDR)")
        plt.ylabel("Término")
        plt.title(f"Top {len(sub)} términos — {db}")
        plt.gca().invert_yaxis()
        fname = out_plots / f"{db.lower()}_barplot.png"
        plt.tight_layout()
        plt.savefig(fname, dpi=150)
        plt.close()
        logging.info("Guardado plot: %s", fname)


# --------------------------------- CLI/Main -----------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Análisis funcional (anotación + enriquecimiento opcional)"
    )
    parser.add_argument("--input", required=True, help="Ruta al archivo de genes de entrada (símbolos)")
    parser.add_argument("--outdir", required=False, default="results", help="Directorio de salida (por defecto: results)")
    parser.add_argument("--output", required=False, default=None, help="Ruta a un CSV principal de salida (opcional)")
    parser.add_argument("--organism", required=False, default="hsapiens", help="Organismo para MyGene/g:Profiler (por defecto: hsapiens)")
    parser.add_argument("--run-enrichment", action="store_true", help="Ejecutar g:Profiler ORA")
    parser.add_argument("--plot", action="store_true", help="Generar barplots con resultados de enriquecimiento")
    parser.add_argument("--topn", type=int, default=10, help="N de términos a graficar por base de datos (default=10)")
    parser.add_argument("--verbose", action="store_true", help="Modo detallado de logging")

    args = parser.parse_args(argv)
    setup_logging(args.verbose)

    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        genes = read_gene_list(input_path)
    except Exception as e:
        logging.error("Error leyendo genes desde %s: %s", input_path, e)
        return 2

    # 1) Mapeo y anotación gene-céntrica con MyGene.info
    try:
        mg_res = map_and_annotate_with_mygene(genes, organism=args.organism)
    except ImportError as e:
        logging.error(str(e))
        return 3
    except Exception as e:
        logging.exception("Fallo en mapeo/anotación con MyGene.info: %s", e)
        return 4

    mapping_df = mg_res["mapping"]
    ann_df = mg_res["annotations"]

    # Guardar salidas base
    mapping_csv = outdir / "id_mapping.csv"
    annotations_csv = outdir / "functional_annotations.csv"
    mapping_df.to_csv(mapping_csv, index=False)
    ann_df.to_csv(annotations_csv, index=False)
    logging.info("Guardado: %s (%d filas)", mapping_csv, len(mapping_df))
    logging.info("Guardado: %s (%d filas)", annotations_csv, len(ann_df))

    enrichment_csv = None
    enrich_df = pd.DataFrame()

    # 2) Enriquecimiento opcional con g:Profiler
    if args.run_enrichment:
        try:
            enrich_df = run_gprofiler_enrichment(genes, organism=args.organism)
            enrichment_csv = outdir / "enrichment_table.csv"
            enrich_df.to_csv(enrichment_csv, index=False)
            logging.info("Guardado: %s (%d filas)", enrichment_csv, len(enrich_df))
        except ImportError as e:
            logging.error(str(e))
            return 5
        except Exception as e:
            logging.exception("Fallo en enriquecimiento con g:Profiler: %s", e)
            return 6

    # 3) Gráficos (si procede)
    if args.plot and not enrich_df.empty:
        try:
            plot_top_terms(enrich_df, outdir, topn=args.topn)
        except Exception as e:
            logging.exception("Fallo generando gráficos: %s", e)
            # No abortamos el script si fallan los plots

    # 4) CSV principal (si el usuario quiere un único archivo de salida)
    if args.output is not None:
        # intentamos componer un resumen compacto
        summary_path = Path(args.output)
        try:
            # Tomamos top 20 términos por FDR si existen
            if not enrich_df.empty and "fdr" in enrich_df.columns:
                top_terms = enrich_df.sort_values("fdr").head(20)[
                    ["db", "term_id", "term_name", "fdr", "overlap", "term_size", "intersection_genes"]
                ].copy()
            else:
                top_terms = pd.DataFrame(columns=["db", "term_id", "term_name", "fdr", "overlap", "term_size", "intersection_genes"])

            # Guardamos un CSV "multi-hoja" concatenando con separadores
            with open(summary_path, "w", encoding="utf-8") as fh:
                fh.write("# id_mapping\n")
                mapping_df.to_csv(fh, index=False)
                fh.write("\n# functional_annotations\n")
                ann_df.to_csv(fh, index=False)
                fh.write("\n# enrichment_top20\n")
                top_terms.to_csv(fh, index=False)
            logging.info("Guardado resumen: %s", summary_path)
        except Exception as e:
            logging.exception("No se pudo escribir el CSV de salida principal (--output): %s", e)
            return 7

    logging.info("Análisis funcional completado correctamente.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
