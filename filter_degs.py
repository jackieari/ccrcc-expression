import pandas as pd

df = pd.read_csv("deseq2_results.csv")

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
