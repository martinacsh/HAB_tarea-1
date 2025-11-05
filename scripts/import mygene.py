import mygene

# Para obtener la información de cada gen
mg = mygene.MyGeneInfo()

# Lista de genes a analizar
genes = ["COX4I2", "ND1", "ATP6"]

# Consultar información funcional
results = mg.querymany(
    genes,
    scopes="symbol",
    fields="symbol,name,entrezgene,go.BP,go.MF,go.CC,pathway.kegg",
    species="human"
)

# Mostrar resultados resumidos
for r in results:
    print(f"\nGene: {r.get('symbol')}")
    print(f"Name: {r.get('name')}")
    print(f"Entrez ID: {r.get('entrezgene')}")
    if "pathway" in r:
        print(f"KEGG Pathways: {r['pathway'].get('kegg', 'N/A')}")
    if "go" in r:
        print("GO Terms:")
        for category, terms in r["go"].items():
            print(f"  {category}: {[t['term'] for t in terms]}")
