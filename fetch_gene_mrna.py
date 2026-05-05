"""
Fetch RefSeq mRNA sequences for a list of human genes from NCBI and write
them to a FASTA file with headers: >NM_XXXXXX Human [gene name]

Strategy (each step is tried only if the previous finds no NM_):
  1. Direct nuccore search (fast, covers most genes).
  2. Gene DB -> elink -> batch-fetch all linked RefSeq RNA IDs (handles genes
     like PBRM1 where XM_ records crowd out NM_ in the first search page).
  3. Allow NR_ accessions as a last resort (lncRNAs like PVT1).
"""

from Bio import Entrez, SeqIO
import time
import sys

Entrez.email = "dinoparadise45@gmail.com"

GENES = [
    "CA9", "EGLN3", "HILPDA", "VEGFA", "UMOD",
    "MFSD4A", "PVT1", "NDUFA4L2", "SPAG4",
    "VHL", "PBRM1", "SETD2", "BAP1",
]

OUTPUT_FILE = "human_gene_mrna_sequences.fasta"
DELAY = 0.4  # seconds between requests (NCBI: 3 req/sec without API key)


# ---------------------------------------------------------------------------
# NCBI helper functions
# ---------------------------------------------------------------------------

def esearch(db: str, term: str, retmax: int = 50) -> list[str]:
    handle = Entrez.esearch(db=db, term=term, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    time.sleep(DELAY)
    return record["IdList"]


def accessions_for_ids(ids: list[str]) -> list[str]:
    """Return accession strings (e.g. NM_001216.3) for a list of GI/UID IDs."""
    if not ids:
        return []
    handle = Entrez.efetch(
        db="nuccore", id=",".join(ids), rettype="acc", retmode="text"
    )
    text = handle.read().strip()
    handle.close()
    time.sleep(DELAY)
    return text.splitlines()


def pick_best(accessions: list[str], prefix: str) -> str | None:
    """Return the lowest-versioned accession with the given prefix, or None."""
    candidates = [a for a in accessions if a.startswith(prefix)]
    if not candidates:
        return None
    candidates.sort(key=lambda a: int(a.split(".")[-1]) if "." in a else 0)
    return candidates[0]


def all_linked_accessions(gene_uid: str) -> list[str]:
    """
    Use Gene->nuccore elink to collect every RefSeq RNA accession for a gene,
    fetching in batches of 50 to stay within URL-length limits.
    """
    handle = Entrez.elink(
        dbfrom="gene", db="nuccore", id=gene_uid, linkname="gene_nuccore_refseqrna"
    )
    lr = Entrez.read(handle)
    handle.close()
    time.sleep(DELAY)

    ids = [
        str(link["Id"])
        for block in lr
        for ls in block.get("LinkSetDb", [])
        for link in ls["Link"]
    ]

    all_accs: list[str] = []
    for i in range(0, len(ids), 50):
        all_accs.extend(accessions_for_ids(ids[i : i + 50]))

    return all_accs


# ---------------------------------------------------------------------------
# Per-gene pipeline
# ---------------------------------------------------------------------------

def get_accession(gene: str) -> str | None:
    """
    Try up to three strategies to find the best RefSeq accession for a gene.
    Returns an accession string like 'NM_001216.3', or None.
    """
    # --- Strategy 1: direct nuccore search (fast) ---------------------------
    ids = esearch(
        "nuccore",
        f"{gene}[Gene Name] AND Homo sapiens[Organism] AND RefSeq[Filter] AND mRNA[Filter]",
        retmax=50,
    )
    if ids:
        accs = accessions_for_ids(ids)
        nm = pick_best(accs, "NM_")
        if nm:
            return nm

    # --- Strategy 2: Gene DB -> elink -> batch fetch (catches PBRM1-style) --
    gene_ids = esearch(
        "gene",
        f"{gene}[Gene Name] AND Homo sapiens[Organism]",
        retmax=5,
    )
    if gene_ids:
        accs = all_linked_accessions(gene_ids[0])
        nm = pick_best(accs, "NM_")
        if nm:
            return nm

    # --- Strategy 3: NR_ fallback for lncRNAs (e.g. PVT1) ------------------
    ids = esearch(
        "nuccore",
        f"{gene}[Gene Name] AND Homo sapiens[Organism] AND RefSeq[Filter]",
        retmax=50,
    )
    if ids:
        accs = accessions_for_ids(ids)
        nr = pick_best(accs, "NR_")
        if nr:
            return nr

    return None


def fetch_sequence(accession: str) -> SeqIO.SeqRecord:
    handle = Entrez.efetch(
        db="nuccore", id=accession, rettype="fasta", retmode="text"
    )
    record = next(SeqIO.parse(handle, "fasta"))
    handle.close()
    time.sleep(DELAY)
    return record


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    results: list[SeqIO.SeqRecord] = []
    failed: list[str] = []

    print(f"Fetching RefSeq mRNA sequences for {len(GENES)} human genes...")
    print("-" * 65)

    for gene in GENES:
        print(f"  {gene:<12}", end="", flush=True)
        try:
            accession = get_accession(gene)
            if accession is None:
                print("NOT FOUND")
                failed.append(gene)
                continue

            record = fetch_sequence(accession)
            base = accession.split(".")[0]  # strip .version suffix
            record.id = base
            record.name = base
            record.description = f"Human {gene}"
            results.append(record)

            tag = "NR" if base.startswith("NR_") else "NM"
            print(f"-> {base}  ({len(record.seq):,} nt)  [{tag}]")

        except Exception as exc:
            print(f"ERROR: {exc}")
            failed.append(gene)
            time.sleep(1)

    print("-" * 65)
    print(f"\nFetched {len(results)}/{len(GENES)} sequences.")

    if results:
        with open(OUTPUT_FILE, "w") as fh:
            SeqIO.write(results, fh, "fasta")
        print(f"Saved  -> {OUTPUT_FILE}")

    if failed:
        print(f"Failed: {', '.join(failed)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
