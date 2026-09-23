import pandas as pd

SEP = "=" * 60

files = [
    ("abstracts_brain_genes.csv",  "brain_genes",  "abstracts_brain_genes.txt"),
    ("abstracts_gut_genes.csv",    "gut_genes",    "abstracts_gut_genes.txt"),
    ("abstracts_shared_genes.csv", "shared_genes", "abstracts_shared_genes.txt"),
]

for csv_path, gene_col, out_path in files:
    df = pd.read_csv(csv_path)
    lines = []
    for i, row in df.iterrows():
        genes = row[gene_col] if pd.notna(row[gene_col]) else ""
        abstract = row["abstract"] if pd.notna(row["abstract"]) else ""
        lines.append(f"[{i}, {row['disorder']}, {genes}, [{abstract}]]")
        lines.append(SEP)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"Saved {out_path}  ({len(df):,} entries)")
