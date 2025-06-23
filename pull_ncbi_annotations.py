import sys
import csv
import pandas as pd
from Bio import Entrez
from tqdm import tqdm

Entrez.email = "your.email@example.com"  # Replace with your email

STANDARD_RANKS = [
    "superkingdom", "kingdom", "phylum", "class", "order",
    "family", "subfamily", "genus"
]

def get_taxonomy_with_fallback(accession):
    try:
        handle = Entrez.elink(dbfrom="nuccore", db="taxonomy", id=accession, linkname="nuccore_taxonomy")
        result = Entrez.read(handle)
        handle.close()

        links = result[0].get("LinkSetDb", [])
        if not links:
            print(f"[!] No taxonomy found for {accession}")
            return None

        taxid = links[0]["Link"][0]["Id"]

        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        record = Entrez.read(handle)[0]
        handle.close()

        lineage_ex = record.get("LineageEx", [])
        lineage_dict = {
            entry["Rank"]: entry["ScientificName"]
            for entry in lineage_ex
            if entry["Rank"] != "no rank"
        }

        taxonomy = {
            "accession": accession,
            "taxid": taxid,
            "description": record.get("ScientificName", "")
        }

        for rank in STANDARD_RANKS:
            taxonomy[rank] = lineage_dict.get(rank, "unclassified")

        return taxonomy

    except Exception as e:
        print(f"[ERROR] Failed for {accession}: {e}")
        return None


def read_accessions_from_csv(csv_path, accession_column="accession"):
    accessions = []
    with open(csv_path, newline="", encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            acc = row.get(accession_column)
            if acc:
                accessions.append(acc.strip())
    return accessions


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python taxonomy_from_csv.py path/to/accessions.csv [accession_column]")
        sys.exit(1)

    csv_path = sys.argv[1]
    accession_column = sys.argv[2] if len(sys.argv) > 2 else "accession"

    accessions = read_accessions_from_csv(csv_path, accession_column)
    print(f"Read {len(accessions)} accessions from {csv_path}")

    taxonomy_data = []
    for acc in tqdm(accessions):
        tax = get_taxonomy_with_fallback(acc)
        if tax:
            taxonomy_data.append(tax)

    if taxonomy_data:
        df = pd.DataFrame(taxonomy_data)
        # Reorder columns nicely
        cols = ["accession", "description"] + STANDARD_RANKS + ["taxid"]
        df = df[cols]
        output_path = csv_path.rsplit(".", 1)[0] + "_taxonomy_output.csv"
        df.to_csv(output_path, index=False)
        print(f"Saved taxonomy info to {output_path}")
    else:
        print("No taxonomy data fetched.")
