"""
Filter DESeq2 results to statistically and biologically significant DEGs.

Input:
  - deseq2_results.csv : full DESeq2 output (gene, log2FoldChange, pvalue, padj)

Output:
  - degs_significant.csv : rows passing padj < 0.05 AND |log2FC| >= 1 (2-fold change),
                           sorted by adjusted p-value ascending
"""

import pandas as pd

df = pd.read_csv("deseq2_results.csv")

# Require both statistical (padj < 0.05) and biological (>=2-fold) significance
# to reduce low-effect false positives common in large RNA-seq datasets
sig = df[
    (df["padj"] < 0.05) &
    (df["log2FoldChange"].abs() >= 1)
].copy()

sig_up = sig[sig["log2FoldChange"] >= 1]
sig_down = sig[sig["log2FoldChange"] <= -1]

sig.sort_values("padj").to_csv("degs_significant.csv", index=False)

print(f"Total significant DEGs: {len(sig)}")
print(f"  Upregulated in tumor:   {len(sig_up)}")
print(f"  Downregulated in tumor: {len(sig_down)}")
print("Saved to degs_significant.csv")
