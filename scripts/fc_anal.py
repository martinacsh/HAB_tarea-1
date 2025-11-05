#!/usr/bin/env python3
"""
func_an.py — Análisis funcional de genes (anotación + enriquecimiento opcional)

Objetivo
--------
- Ejecutar TODO con **un solo comando** si quieres: `python scripts/func_an.py --all`
- Por defecto usa `data/genes_input.txt`, escribe en `results/`, y genera `results/resumen.csv`.

Qué hace
--------
1) Lee símbolos génicos desde `--input` (por defecto: `data/genes_input.txt`).
2) Mapea IDs y anota GO/KEGG/Reactome por gen usando **MyGene.info** (`mygene`).
3) (Opcional) Ejecuta enriquecimiento ORA con **g:Profiler** (`gprofiler-official`) para GO (BP/MF/CC), KEGG y Reactome, con FDR.
4) Exporta CSV y, si se solicita, genera barplots (`matplotlib`).

Notas importantes
-----------------
- **Especie/organismo**: MyGene y g:Profiler usan códigos distintos. Este script **normaliza**:
  - MyGene ← `human` (o 9606)
  - g:Profiler ← `hsapiens`
- Si g:Profiler falla (p.ej., error 500), el script **no se detiene** y continúa sin ORA.

Dependencias
------------
- pandas, numpy, mygene, gprofiler-official (para ORA), matplotlib

Uso rápido
---------
- Todo en uno: `python scripts/func_an.py --all`
- Personalizado:
  `python scripts/func_an.py --input data/genes_input.txt --outdir results --organism hsapiens --run-enrichment --plot --output results/resumen.csv`
"""
from __future__ import annotations
import argparse
import sys
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional

from time import sleep
import pandas as pd

# -----------------------------------------------------------------------------
# Logging y utilidades
# -----------------------------------------------------------------------------

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
    - Acepta uno por línea o separados por comas/espacios.
    - Ignora líneas vacías y las que comienzan por '#'.
    - Devuelve lista única (sin duplicados), en mayúsculas.
    """
    if not path.exists():
        raise FileNotFoundError(f"No existe el archivo de entrada: {path}")

    raw: List[str] = []
    text = path.read_text(encoding="utf-8", errors="ignore")
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = [p for chunk in line.split(',') for p in chunk.split()] if (',' in line or ' ' in line) else [line]
        raw.extend(parts)
    genes = [g.strip().upper() for g in raw if g.strip()]

    seen = set()
    unique: List[str] = []
    for g in genes:
        if g not in seen:
            unique.append(g)
            seen.add(g)
    if not unique:
        raise ValueError("La lista de genes está vacía tras el preprocesado.")
    logging.info("Genes de entrada (%d): %s", len(unique), ", ".join(unique))
    return unique


# -----------------------------------------------------------------------------
# Normalización de organismos para MyGene (human) y g:Profiler (hsapiens)
# -----------------------------------------------------------------------------

def _normalize_mygene_species(organism: str) -> str:
    """MyGene acepta 'human' o taxid '9606', NO 'hsapiens'."""
    org = (organism or '').strip().lower()
    mapping = {
        'hsapiens': 'human', 'homo_sapiens': 'human', 'homo sapiens': 'human',
        'human': 'human', 'hsa': 'human', '9606': 'human',
        'mmusculus': 'mouse', 'mus_musculus': 'mouse', 'mus musculus': 'mouse',
        'mouse': 'mouse', 'mmu': 'mouse', '10090': 'mouse',
        'rnorvegicus': 'rat', 'rattus_norvegicus': 'rat',
        'rat': 'rat', 'rnor': 'rat', '10116': 'rat',
    }
    return mapping.get(org, 'human')


def _normalize_gprofiler_organism(organism: str) -> str:
    """g:Profiler espera 'hsapiens' para humano."""
    o = (organism or '').strip().lower()
    if o in {'human', 'homo_sapiens', 'homo sapiens', 'hsa', '9606', 'hsapiens'}:
        return 'hsapiens'
    return o


# -----------------------------------------------------------------------------
# Mapeo y anotación con MyGene.info
# -----------------------------------------------------------------------------

def map_and_annotate_with_mygene(genes: List[str], organism: str = "hsapiens") -> Dict[str, pd.DataFrame]:
    """
    Usa MyGene.info para:
      - Mapear símbolos a IDs (Entrez/Ensembl/UniProt) y nombre canónico
      - Recuperar anotaciones GO (BP/MF/CC) y pathways (KEGG/Reactome) por gen

    Devuelve un dict con:
      - mapping: DataFrame una fila por símbolo (IDs, nombres)
      - annotations: DataFrame (gene, source, term_id, term_name)
    """
    try:
        import mygene  # type: ignore
    except ImportError as e:
        raise ImportError("Se requiere la librería 'mygene'. Añádela a requirements.txt e instálala.") from e

    mg = mygene.MyGeneInfo()
    mygene_species = _normalize_mygene_species(organism)

    fields = [
        "entrezgene",
        "ensembl.gene",
        "uniprot.Swiss-Prot",
        "name",
        "symbol",
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
        species=mygene_species,
        as_dataframe=True,
        returnall=False,
        verbose=False,
    )

    mapping_rows: List[Dict[str, Any]] = []
    annot_rows: List[Dict[str, Any]] = []

    def _ensure_list(x: Any) -> List[Any]:
        if x is None or (isinstance(x, float) and pd.isna(x)):
            return []
        return x if isinstance(x, list) else [x]

    for query_symbol, row in qres.iterrows():
        # IDs
        entrez = row.get("entrezgene")
        ensembl_val = row.get("ensembl.gene")
        ensembl = ensembl_val[0] if isinstance(ensembl_val, list) else ensembl_val
        uniprot_val = row.get("uniprot.Swiss-Prot")
        uniprot = uniprot_val[0] if isinstance(uniprot_val, list) else uniprot_val
        symbol = row.get("symbol", query_symbol)
        name = row.get("name")

        mapping_rows.append({
            "query_symbol": query_symbol,
            "symbol": symbol,
            "name": name,
            "entrezgene": entrez,
            "ensembl_gene": ensembl,
            "uniprot_swissprot": uniprot,
        })

        # GO terms
        for cat in ("BP", "MF", "CC"):
            go_list = row.get(f"go.{cat}")
            for term in _ensure_list(go_list):
                term_id = term.get("id") if isinstance(term, dict) else None
                term_name = term.get("term") if isinstance(term, dict) else (str(term) if term is not None else None)
                if term_id or term_name:
                    annot_rows.append({
                        "gene": symbol,
                        "source": f"GO:{cat}",
                        "term_id": term_id,
                        "term_name": term_name,
                    })
        # KEGG
        for term in _ensure_list(row.get("pathway.kegg")):
            term_id = term.get("id") if isinstance(term, dict) else None
            term_name = term.get("name") if isinstance(term, dict) else (str(term) if term is not None else None)
            if term_id or term_name:
                annot_rows.append({
                    "gene": symbol,
                    "source": "KEGG",
                    "term_id": term_id,
                    "term_name": term_name,
                })
        # Reactome
        for term in _ensure_list(row.get("pathway.reactome")):
            term_id = term.get("id") if isinstance(term, dict) else None
            term_name = term.get("name") if isinstance(term, dict) else (str(term) if term is not None else None)
            if term_id or term_name:
                annot_rows.append({
                    "gene": symbol,
                    "source": "REAC",
                    "term_id": term_id,
                    "term_name": term_name,
                })

    mapping_df = pd.DataFrame(mapping_rows)
    annotations_df = pd.DataFrame(annot_rows)
    return {"mapping": mapping_df, "annotations": annotations_df}


# -----------------------------------------------------------------------------
# Enriquecimiento con g:Profiler (con normalización y reintento)
# -----------------------------------------------------------------------------

def run_gprofiler_enrichment(genes: List[str], organism: str = "hsapiens",
                              sources: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Ejecuta ORA usando g:Profiler y devuelve un DataFrame con p-value y FDR.
    """
    if sources is None:
        sources = ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"]

    try:
        from gprofiler import GProfiler  # type: ignore
    except ImportError as e:
        raise ImportError("Se requiere 'gprofiler-official'. Añádelo a requirements.txt e instálalo.") from e

    gp = GProfiler(return_dataframe=True)
    gp_org = _normalize_gprofiler_organism(organism)

    # Reintento ligero por errores temporales (p.ej. 500)
    last_err: Optional[BaseException] = None
    for attempt in range(2):
        try:
            res_df = gp.profile(
                organism=gp_org,
                query=genes,
                sources=sources,
                user_threshold=1.0,
            )
            break
        except AssertionError as e:
            last_err = e
            logging.warning("g:Profiler falló (intento %d/2): %s", attempt + 1, e)
            sleep(1.0)
    else:
        # No se pudo recuperar; propaga el último error
        raise last_err  # type: ignore[misc]

    # Normalización de columnas
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

    sort_cols = [c for c in ("fdr", "pvalue") if c in res_df.columns]
    if sort_cols:
        res_df = res_df.sort_values(by=sort_cols, ascending=True)

    return res_df


# -----------------------------------------------------------------------------
# Gráficos
# -----------------------------------------------------------------------------

def plot_top_terms(res_df: pd.DataFrame, outdir: Path, topn: int = 10) -> None:
    import math
    import matplotlib.pyplot as plt

    if res_df.empty:
        logging.warning("No hay resultados de enriquecimiento para graficar.")
        return

    out_plots = outdir / "plots"
    out_plots.mkdir(parents=True, exist_ok=True)

    df = res_df.copy()
    if "fdr" not in df.columns:
        logging.warning("La columna 'fdr' no está; se usará 'pvalue'.")
        df["fdr"] = df.get("pvalue", 1.0)

    df["neglog10_fdr"] = df["fdr"].apply(lambda x: -math.log10(x) if x > 0 else 0)

    for db in sorted(df["db"].dropna().unique()):
        sub = df[df["db"] == db].head(topn)
        if sub.empty:
            continue
        # Evita NaN en nombres
        labels = sub["term_name"].fillna(sub["term_id"]).astype(str)
        values = sub["neglog10_fdr"].astype(float)

        plt.figure(figsize=(10, max(3, 0.45 * len(sub))))
        plt.barh(labels, values)
        plt.xlabel("-log10(FDR)")
        plt.ylabel("Término")
        plt.title(f"Top {len(sub)} términos — {db}")
        plt.gca().invert_yaxis()
        fname = out_plots / f"{db.lower()}_barplot.png"
        plt.tight_layout()
        plt.savefig(fname, dpi=150)
        plt.close()
        logging.info("Guardado plot: %s", fname)


# -----------------------------------------------------------------------------
# CLI / main
# -----------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Anotación gene-céntrica + enriquecimiento (ORA) opcional"
    )
    parser.add_argument("--input", required=False, default="data/genes_input.txt",
                        help="Archivo de genes (por defecto: data/genes_input.txt)")
    parser.add_argument("--outdir", required=False, default="results",
                        help="Directorio de salida (por defecto: results)")
    parser.add_argument("--output", required=False, default="results/resumen.csv",
                        help="CSV resumen (por defecto: results/resumen.csv)")
    parser.add_argument("--organism", required=False, default="hsapiens",
                        help="Organismo (se adapta a cada servicio). Por defecto: hsapiens")
    parser.add_argument("--run-enrichment", action="store_true",
                        help="Ejecutar g:Profiler ORA")
    parser.add_argument("--plot", action="store_true",
                        help="Generar barplots de términos enriquecidos")
    parser.add_argument("--topn", type=int, default=10, help="N de términos a graficar por base de datos (10)")
    parser.add_argument("--verbose", action="store_true", help="Logging detallado")
    parser.add_argument("--all", action="store_true",
                        help="Ejecuta todo: ORA + plots y genera resumen (usa defaults)")

    args = parser.parse_args(argv)

    # Modo "todo en uno"
    if args.all:
        args.run_enrichment = True
        args.plot = True
        # Mantiene defaults de input/outdir/output definidos arriba

    setup_logging(args.verbose)

    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) Lectura de genes
    try:
        genes = read_gene_list(input_path)
    except Exception as e:
        logging.error("Error leyendo genes desde %s: %s", input_path, e)
        return 2

    # 2) Mapeo y anotación con MyGene
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

    # Guardados base
    mapping_csv = outdir / "id_mapping.csv"
    annotations_csv = outdir / "functional_annotations.csv"
    try:
        mapping_df.to_csv(mapping_csv, index=False)
        ann_df.to_csv(annotations_csv, index=False)
        logging.info("Guardado: %s (%d filas)", mapping_csv, len(mapping_df))
        logging.info("Guardado: %s (%d filas)", annotations_csv, len(ann_df))
    except Exception as e:
        logging.exception("No se pudieron escribir las tablas base: %s", e)
        return 5

    # 3) Enriquecimiento (robusto)
    enrich_df = pd.DataFrame()
    if args.run_enrichment:
        try:
            enrich_df = run_gprofiler_enrichment(genes, organism=args.organism)
            enrichment_csv = outdir / "enrichment_table.csv"
            enrich_df.to_csv(enrichment_csv, index=False)
            logging.info("Guardado: %s (%d filas)", enrichment_csv, len(enrich_df))
        except ImportError as e:
            logging.warning("%s. Continuo sin enriquecimiento.", e)
        except Exception as e:
            logging.warning("Enriquecimiento falló: %s. Continuo sin ORA.", e)

    # 4) Gráficos (si procede)
    if args.plot and not enrich_df.empty:
        try:
            plot_top_terms(enrich_df, outdir, topn=args.topn)
        except Exception as e:
            logging.warning("Fallo generando gráficos: %s", e)

    # 5) CSV resumen único
    if args.output:
        summary_path = Path(args.output)
        try:
            # Top 20 términos por FDR si existen
            if not enrich_df.empty and "fdr" in enrich_df.columns:
                cols = [c for c in ["db", "term_id", "term_name", "fdr", "overlap", "term_size", "intersection_genes"] if c in enrich_df.columns]
                top_terms = enrich_df.sort_values("fdr").head(20)[cols].copy()
            else:
                top_terms = pd.DataFrame(columns=["db", "term_id", "term_name", "fdr", "overlap", "term_size", "intersection_genes"])

            with open(summary_path, "w", encoding="utf-8") as fh:
                fh.write("# id_mapping\n")
                mapping_df.to_csv(fh, index=False)
                fh.write("\n# functional_annotations\n")
                ann_df.to_csv(fh, index=False)
                fh.write("\n# enrichment_top20\n")
                top_terms.to_csv(fh, index=False)
            logging.info("Guardado resumen: %s", summary_path)
        except Exception as e:
            logging.warning("No se pudo escribir el CSV resumen (%s): %s", summary_path, e)

    logging.info("Análisis funcional completado.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
