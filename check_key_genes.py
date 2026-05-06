"""
Report DESeq2 statistics for a curated set of KIRC-relevant genes.

Genes include known clear cell RCC tumor suppressors (VHL, PBRM1, SETD2, BAP1)
and hypoxia/angiogenesis response markers (CA9, EGLN3, HILPDA, VEGFA, UMOD).

Input:
  - deseq2_results_annotated.csv : annotated full DESeq2 results

Output:
  - Formatted summary table printed to stdout (no file written)
"""

import pandas as pd

# Canonical clear cell RCC driver genes and hypoxia markers from the literature
genes_of_interest = ["VHL", "CA9", "PBRM1", "SETD2", "BAP1", "EGLN3", "HILPDA", "VEGFA", "UMOD"]

df = pd.read_csv("deseq2_results_annotated.csv")

print("Key Genes — TCGA-KIRC Tumor vs Normal (DESeq2)")
print(f"{'Gene':<10} {'log2FoldChange':>16} {'padj':>14} {'Significant'}")
print("-" * 56)

for gene in genes_of_interest:
    match = df[df["gene_symbol"] == gene]
    if not match.empty:
        row = match.iloc[0]
        padj = row["padj"]
        # Match the same thresholds used in filter_degs.py: padj < 0.05 and >=2-fold
        sig = "yes" if pd.notna(padj) and padj < 0.05 and abs(row["log2FoldChange"]) >= 1 else "no"
        padj_str = f"{padj:.2e}" if pd.notna(padj) else "NA"
        print(f"{gene:<10} {row['log2FoldChange']:>16.4f} {padj_str:>14}  {sig}")
    else:
        print(f"{gene:<10} {'N/A':>16} {'N/A':>14}  not found")
