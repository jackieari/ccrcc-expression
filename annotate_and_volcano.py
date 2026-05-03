"""
Annotate Ensembl gene IDs with gene symbols and create a volcano plot.

Inputs:
  - deseq2_results.csv
  - degs_significant.csv

Outputs:
  - gene_annotations.tsv
  - deseq2_results_annotated.csv
  - degs_significant_annotated.csv
  - volcano_plot.png
  - volcano_plot.pdf
"""

from pathlib import Path
import urllib.parse
import urllib.request

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


RESULTS_FILE = Path("deseq2_results.csv")
SIG_FILE = Path("degs_significant.csv")
ANNOTATION_FILE = Path("gene_annotations.tsv")
ANNOTATED_RESULTS_FILE = Path("deseq2_results_annotated.csv")
ANNOTATED_SIG_FILE = Path("degs_significant_annotated.csv")
VOLCANO_PNG = Path("volcano_plot.png")
VOLCANO_PDF = Path("volcano_plot.pdf")

PADJ_CUTOFF = 0.05
LOG2FC_CUTOFF = 1.0


def strip_version(gene_id: str) -> str:
    """Turn ENSG000001.12 into ENSG000001 so it can match Ensembl annotations."""
    return str(gene_id).split(".", 1)[0]


def download_ensembl_annotations() -> None:
    """Download human gene ID to gene-symbol annotations from Ensembl BioMart."""
    query = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1"
       count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Attribute name="ensembl_gene_id" />
    <Attribute name="external_gene_name" />
    <Attribute name="gene_biotype" />
    <Attribute name="description" />
  </Dataset>
</Query>
"""
    data = urllib.parse.urlencode({"query": query}).encode()
    request = urllib.request.Request(
        "https://www.ensembl.org/biomart/martservice",
        data=data,
        headers={"User-Agent": "ccrcc-expression-volcano/1.0"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        ANNOTATION_FILE.write_text(response.read().decode("utf-8"))


def load_annotations() -> pd.DataFrame:
    if not ANNOTATION_FILE.exists():
        print(f"Downloading Ensembl annotations to {ANNOTATION_FILE}...")
        download_ensembl_annotations()
    else:
        print(f"Using existing {ANNOTATION_FILE}")

    annotations = pd.read_csv(ANNOTATION_FILE, sep="\t")
    annotations = annotations.rename(
        columns={
            "Gene stable ID": "gene_base",
            "Gene name": "gene_symbol",
            "Gene type": "gene_biotype",
            "Gene description": "gene_description",
        }
    )
    annotations["gene_symbol"] = annotations["gene_symbol"].fillna(annotations["gene_base"])
    return annotations[["gene_base", "gene_symbol", "gene_biotype", "gene_description"]]


def annotate_results(path: Path, annotations: pd.DataFrame) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["gene_base"] = df["gene"].map(strip_version)
    df = df.merge(annotations, on="gene_base", how="left")
    df["gene_symbol"] = df["gene_symbol"].fillna(df["gene_base"])

    columns = ["gene", "gene_symbol", "gene_biotype", "gene_description"]
    columns += [col for col in df.columns if col not in columns + ["gene_base"]]
    return df[columns]


def add_volcano_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    positive_padj = df.loc[df["padj"] > 0, "padj"]
    padj_floor = positive_padj.min() / 10

    df["padj_for_plot"] = df["padj"].clip(lower=padj_floor)
    df["neg_log10_padj"] = -np.log10(df["padj_for_plot"])

    df["volcano_group"] = "Not significant"
    df.loc[
        (df["padj"] < PADJ_CUTOFF) & (df["log2FoldChange"] >= LOG2FC_CUTOFF),
        "volcano_group",
    ] = "Up in tumor"
    df.loc[
        (df["padj"] < PADJ_CUTOFF) & (df["log2FoldChange"] <= -LOG2FC_CUTOFF),
        "volcano_group",
    ] = "Down in tumor"
    return df


def make_volcano_plot(df: pd.DataFrame) -> None:
    df = add_volcano_columns(df.dropna(subset=["log2FoldChange", "padj"]))

    colors = {
        "Not significant": "#b8b8b8",
        "Up in tumor": "#d73027",
        "Down in tumor": "#4575b4",
    }
    sizes = {
        "Not significant": 8,
        "Up in tumor": 14,
        "Down in tumor": 14,
    }

    fig, ax = plt.subplots(figsize=(11, 8), dpi=180)

    for group in ["Not significant", "Down in tumor", "Up in tumor"]:
        subset = df[df["volcano_group"] == group]
        ax.scatter(
            subset["log2FoldChange"],
            subset["neg_log10_padj"],
            s=sizes[group],
            c=colors[group],
            alpha=0.35 if group == "Not significant" else 0.75,
            linewidths=0,
            rasterized=True,
            label=f"{group} ({len(subset):,})",
        )

    ax.axvline(-LOG2FC_CUTOFF, color="#555555", linestyle="--", linewidth=1)
    ax.axvline(LOG2FC_CUTOFF, color="#555555", linestyle="--", linewidth=1)
    ax.axhline(-np.log10(PADJ_CUTOFF), color="#555555", linestyle="--", linewidth=1)

    label_df = df[df["volcano_group"].isin(["Up in tumor", "Down in tumor"])].copy()
    label_df["label_score"] = label_df["neg_log10_padj"] + label_df["log2FoldChange"].abs()
    for _, row in label_df.nlargest(12, "label_score").iterrows():
        ax.annotate(
            row["gene_symbol"],
            (row["log2FoldChange"], row["neg_log10_padj"]),
            xytext=(4 if row["log2FoldChange"] >= 0 else -4, 4),
            textcoords="offset points",
            ha="left" if row["log2FoldChange"] >= 0 else "right",
            fontsize=8,
            fontweight="bold",
        )

    ax.set_title("TCGA-KIRC Tumor vs Normal Volcano Plot", fontsize=16, weight="bold")
    ax.set_xlabel("log2 fold change: tumor vs normal", fontsize=12)
    ax.set_ylabel("-log10 adjusted p-value", fontsize=12)
    ax.legend(frameon=False, loc="upper right")
    ax.grid(True, color="#eeeeee", linewidth=0.8)
    ax.set_axisbelow(True)
    fig.tight_layout()

    fig.savefig(VOLCANO_PNG, bbox_inches="tight")
    fig.savefig(VOLCANO_PDF, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    annotations = load_annotations()
    print(f"Loaded {len(annotations):,} gene annotations")

    results = annotate_results(RESULTS_FILE, annotations)
    sig = annotate_results(SIG_FILE, annotations)

    results.to_csv(ANNOTATED_RESULTS_FILE, index=False)
    sig.to_csv(ANNOTATED_SIG_FILE, index=False)
    print(f"Wrote {ANNOTATED_RESULTS_FILE}")
    print(f"Wrote {ANNOTATED_SIG_FILE}")

    make_volcano_plot(results)
    print(f"Wrote {VOLCANO_PNG}")
    print(f"Wrote {VOLCANO_PDF}")


if __name__ == "__main__":
    main()
