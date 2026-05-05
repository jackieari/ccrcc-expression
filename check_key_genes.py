import pandas as pd

genes_of_interest = ["VHL", "CA9", "PBRM1", "SETD2", "BAP1"]

df = pd.read_csv("degs_significant_annotated.csv")

print(f"{'Gene':<10} {'log2FoldChange':>16} {'padj':>14} {'Found'}")
print("-" * 50)

for gene in genes_of_interest:
    match = df[df["gene_symbol"] == gene]
    if not match.empty:
        row = match.iloc[0]
        print(f"{gene:<10} {row['log2FoldChange']:>16.4f} {row['padj']:>14.2e}  yes")
    else:
        print(f"{gene:<10} {'N/A':>16} {'N/A':>14}  not in significant DEGs")
