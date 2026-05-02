"""
PyDESeq2 differential expression: TCGA-KIRC tumor vs normal.

Inputs:
  - tcga_kirc_counts.csv   : genes x samples (columns = file UUIDs)
  - metadata.repository.2026-05-02.json : GDC metadata with file_id and entity_submitter_id

Outputs:
  - deseq2_results.csv : gene, log2FoldChange, pvalue, padj
"""

import json
import re
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

COUNTS_FILE = "tcga_kirc_counts.csv"
METADATA_FILE = "metadata.repository.2026-05-02.json"
OUTPUT_FILE = "deseq2_results.csv"

# ---------------------------------------------------------------------------
# 1. Parse metadata: file_id -> "tumor" or "normal"
# ---------------------------------------------------------------------------
with open(METADATA_FILE, "r") as fh:
    metadata = json.load(fh)

records = []
for entry in metadata:
    file_id = entry.get("file_id")
    entities = entry.get("associated_entities", [])
    if not file_id or not entities:
        continue
    sample_id = entities[0].get("entity_submitter_id", "")
    # TCGA sample type code is the 4th hyphen-delimited field (e.g. "01B", "11A")
    parts = sample_id.split("-")
    if len(parts) < 4:
        continue
    type_code = parts[3]          # e.g. "01B" or "11A"
    numeric_code = re.match(r"(\d+)", type_code)
    if not numeric_code:
        continue
    code = int(numeric_code.group(1))
    if 1 <= code <= 9:
        label = "tumor"
    elif 10 <= code <= 19:
        label = "normal"
    else:
        continue                  # skip metastatic / other types
    records.append({"file_id": file_id, "sample_id": sample_id, "condition": label})

sample_labels = pd.DataFrame(records).set_index("file_id")
print(f"Metadata parsed: {len(sample_labels)} samples "
      f"({(sample_labels.condition == 'tumor').sum()} tumor, "
      f"{(sample_labels.condition == 'normal').sum()} normal)")

# ---------------------------------------------------------------------------
# 2. Load counts and filter to samples present in metadata
# ---------------------------------------------------------------------------
counts = pd.read_csv(COUNTS_FILE, index_col=0)   # genes x samples
print(f"Counts matrix loaded: {counts.shape[0]} genes x {counts.shape[1]} samples")

shared = counts.columns.intersection(sample_labels.index)
counts = counts[shared]
sample_labels = sample_labels.loc[shared]
print(f"After filtering to shared samples: {counts.shape[1]} samples retained")

# PyDESeq2 expects samples as rows, genes as columns (transposed)
counts_t = counts.T.astype(int)

# Align condition labels
conditions = sample_labels.loc[counts_t.index, "condition"]

# ---------------------------------------------------------------------------
# 3. Run PyDESeq2
# ---------------------------------------------------------------------------
dds = DeseqDataSet(
    counts=counts_t,
    metadata=conditions.to_frame(name="condition"),
    design_factors="condition",
    ref_level=["condition", "normal"],   # normal is the reference
    n_cpus=4,
    quiet=False,
)
dds.deseq2()

stat_res = DeseqStats(dds, contrast=["condition", "tumor", "normal"], quiet=False)
stat_res.summary()

# ---------------------------------------------------------------------------
# 4. Export results
# ---------------------------------------------------------------------------
results = stat_res.results_df.reset_index().rename(columns={"gene_id": "gene"})
out_cols = ["gene", "log2FoldChange", "pvalue", "padj"]
results[out_cols].to_csv(OUTPUT_FILE, index=False)
print(f"Results written to {OUTPUT_FILE} ({len(results)} genes)")
