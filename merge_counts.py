"""
Merges TCGA-KIRC STAR count TSV files into a single gene x sample matrix.

Each sample lives in its own UUID folder under DATA_DIR:
    <DATA_DIR>/<uuid>/<uuid>.rna_seq.augmented_star_gene_counts.tsv

Output columns are the folder UUIDs (use a manifest to map to TCGA barcodes).
Rows are Ensembl gene IDs; the summary rows (N_unmapped, etc.) are dropped.
The count type written to the matrix is controlled by COUNT_COLUMN.
"""

import glob
import os
import sys

import pandas as pd
from tqdm import tqdm  # optional – remove if not installed

# ── configuration ─────────────────────────────────────────────────────────────
DATA_DIR = r"C:\Users\jacki\data"
OUT_FILE = r"C:\Users\jacki\OneDrive\Desktop\bioinfo\tcga_kirc_counts.csv"

# Choose one: unstranded | stranded_first | stranded_second
# For TCGA RNA-seq the library is typically unstranded
COUNT_COLUMN = "unstranded"
# ──────────────────────────────────────────────────────────────────────────────

SUMMARY_PREFIXES = ("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")


def load_counts(tsv_path: str, sample_id: str) -> pd.Series:
    """Read one sample's STAR counts TSV, drop summary rows, return as a named Series."""
    df = pd.read_csv(
        tsv_path,
        sep="\t",
        comment="#",
        usecols=["gene_id", COUNT_COLUMN],
        dtype={"gene_id": str, COUNT_COLUMN: "Int64"},
    )
    # Drop STAR summary rows
    df = df[~df["gene_id"].str.startswith(SUMMARY_PREFIXES)]
    df = df.set_index("gene_id")[COUNT_COLUMN]
    df.name = sample_id
    return df


def main():
    """Discover all sample TSVs under DATA_DIR, load them, and write a genes×samples CSV."""
    pattern = os.path.join(DATA_DIR, "*", "*.rna_seq.augmented_star_gene_counts.tsv")
    tsv_files = sorted(glob.glob(pattern))

    if not tsv_files:
        sys.exit(f"No TSV files found with pattern:\n  {pattern}")

    print(f"Found {len(tsv_files)} files. Building matrix …")

    series_list = []
    try:
        iterator = tqdm(tsv_files, unit="file")
    except NameError:
        iterator = tsv_files  # tqdm not installed

    for path in iterator:
        sample_id = os.path.basename(os.path.dirname(path))  # UUID folder name
        series_list.append(load_counts(path, sample_id))

    print("Concatenating …")
    matrix = pd.concat(series_list, axis=1)  # genes x samples
    matrix.index.name = "gene_id"

    print(f"Matrix shape: {matrix.shape[0]} genes × {matrix.shape[1]} samples")
    print(f"Writing to {OUT_FILE} …")
    matrix.to_csv(OUT_FILE)
    print("Done.")


if __name__ == "__main__":
    main()
