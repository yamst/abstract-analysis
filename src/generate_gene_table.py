import re
import os
import pandas as pd

# Output directory
_OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output")

df = pd.read_csv(os.path.join(_OUTPUT_DIR, "model_organisms_dataset.csv"))

def first_sentence(text):
    if not isinstance(text, str):
        return ""
    lines = text.strip().splitlines()
    clean_lines = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if re.match(r'^(PMID|DOI|PMCID|©|\d{4}\s)', line):
            continue
        line = re.sub(r'^[A-Z][A-Z\s/]+:\s*', '', line)
        if line:
            clean_lines.append(line)
    full = " ".join(clean_lines)
    m = re.search(r'.{30,}?[.!?](?=\s|$)', full)
    if m:
        return m.group(0).strip()
    return full[:200].strip()

def parse_genes(cell):
    if not isinstance(cell, str) or cell.strip() in ("", "nan"):
        return []
    return [g.strip() for g in cell.split(",") if g.strip() and g.strip() != "nan"]

rows = []
for _, row in df.iterrows():
    brain_genes = parse_genes(row.get("brain_genes_found", ""))
    gut_genes   = parse_genes(row.get("gut_genes_found", ""))
    if not brain_genes and not gut_genes:
        continue  # abstract did not pass filters

    disorder = row["disorder"]
    abstract = row["abstract"]
    subject  = first_sentence(abstract)

    for gene in brain_genes:
        rows.append({"subject": subject, "tissue": "Brain", "disorder": disorder, "gene": gene})
    for gene in gut_genes:
        rows.append({"subject": subject, "tissue": "Gut", "disorder": disorder, "gene": gene})

result = pd.DataFrame(rows).drop_duplicates()
result = result.sort_values(["disorder", "tissue", "gene", "subject"]).reset_index(drop=True)

result.to_csv(os.path.join(_OUTPUT_DIR, "gene_tissue_table.csv"), index=False)
print(f"Saved gene_tissue_table.csv — {len(result):,} rows, {result['gene'].nunique()} unique genes")
try:
    print(result.head(10).to_string(index=False))
except UnicodeEncodeError:
    print("(first 10 rows contain non-ASCII characters; see CSV for full output)")
