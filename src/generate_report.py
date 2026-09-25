import re
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib_venn import venn2
from collections import defaultdict

# Output directory
_OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output")

df = pd.read_csv(os.path.join(_OUTPUT_DIR, "model_organisms_dataset.csv"), low_memory=False)

DISORDER_ORDER = [
    "Schizophrenia", "Bipolar Disorder", "Autism", "ADHD",
    "Major Depression", "Anxiety", "PTSD", "OCD",
    "Alzheimer's", "Parkinson's",
]


def parse_genes(cell):
    if not isinstance(cell, str) or cell.strip() in ("", "nan"):
        return []
    return [g.strip() for g in cell.split(",") if g.strip() and g.strip() != "nan"]

def first_sentence(text):
    if not isinstance(text, str):
        return ""
    lines = text.strip().splitlines()
    clean = []
    for line in lines:
        line = line.strip()
        if not line or re.match(r"^(PMID|DOI|PMCID|©|\d{4}\s)", line):
            continue
        line = re.sub(r"^[A-Z][A-Z\s/]+:\s*", "", line)
        if line:
            clean.append(line)
    full = " ".join(clean)
    m = re.search(r".{30,}?[.!?](?=\s|$)", full)
    return m.group(0).strip() if m else full[:200].strip()


brain_genes_by_disorder = defaultdict(set)
gut_genes_by_disorder   = defaultdict(set)
abstract_counts = defaultdict(lambda: defaultdict(int))

from collections import Counter
brain_paper_counts = defaultdict(lambda: defaultdict(int))  
gut_paper_counts   = defaultdict(lambda: defaultdict(int))
gene_reg_votes     = defaultdict(lambda: defaultdict(list)) 

for _, row in df.iterrows():
    brain_genes = parse_genes(row.get("brain_genes_found", ""))
    gut_genes   = parse_genes(row.get("gut_genes_found", ""))
    if not brain_genes and not gut_genes:
        continue

    disorder = row["disorder"]

    # Parse regulation dict stored as string
    reg_raw = row.get("gene_regulation", "")
    reg_map = {}
    if isinstance(reg_raw, str) and reg_raw.startswith("{"):
        try:
            import ast
            reg_map = ast.literal_eval(reg_raw)
        except Exception:
            pass

    if brain_genes:
        abstract_counts[disorder]["Brain"] += 1
        for gene in brain_genes:
            brain_genes_by_disorder[disorder].add(gene)
            brain_paper_counts[gene][disorder] += 1
            if gene in reg_map:
                gene_reg_votes[gene][disorder].append(reg_map[gene])

    if gut_genes:
        abstract_counts[disorder]["Gut"] += 1
        for gene in gut_genes:
            gut_genes_by_disorder[disorder].add(gene)
            gut_paper_counts[gene][disorder] += 1
            if gene in reg_map:
                gene_reg_votes[gene][disorder].append(reg_map[gene])


rows = []
all_genes = set(brain_paper_counts.keys()) | set(gut_paper_counts.keys())
for gene in sorted(all_genes):
    disorders_involved = set(brain_paper_counts[gene]) | set(gut_paper_counts[gene])
    for disorder in DISORDER_ORDER:
        if disorder not in disorders_involved:
            continue
        n_brain = brain_paper_counts[gene][disorder]
        n_gut   = gut_paper_counts[gene][disorder]
        if n_brain > 0 and n_gut > 0:
            tissue = "Both"
        elif n_brain > 0:
            tissue = "Brain"
        else:
            tissue = "Gut"
        votes = gene_reg_votes[gene][disorder]
        up    = votes.count("Upregulated")
        down  = votes.count("Downregulated")
        if up > 0 and down > 0:   reg = "Mixed"
        elif up > 0:               reg = "Upregulated"
        elif down > 0:             reg = "Downregulated"
        else:                      reg = "No directional language"

        rows.append({
            "gene":         gene,
            "disorder":     disorder,
            "tissue":       tissue,
            "brain_papers": n_brain,
            "gut_papers":   n_gut,
            "regulation":   reg,
        })

gene_table = pd.DataFrame(rows)
gene_table = gene_table.sort_values(["gene", "disorder"]).reset_index(drop=True)

def row_tissue(r):
    if r["brain_papers"] > 0 and r["gut_papers"] > 0:
        return "Both"
    elif r["brain_papers"] > 0:
        return "Brain"
    else:
        return "Gut"
gene_table["overall_tissue"] = gene_table.apply(row_tissue, axis=1)

gene_table.to_csv(os.path.join(_OUTPUT_DIR, "consolidated_gene_table.csv"), index=False)
print(f"Saved consolidated_gene_table.csv — {len(gene_table):,} rows, "
      f"{gene_table['gene'].nunique()} unique genes")

#summary stat
summary_rows = []
all_brain = set()
all_gut   = set()
for disorder in DISORDER_ORDER:
    bg = brain_genes_by_disorder[disorder]
    gg = gut_genes_by_disorder[disorder]
    shared = bg & gg
    all_brain |= bg
    all_gut   |= gg
    summary_rows.append({
        "disorder":         disorder,
        "brain_abstracts":  abstract_counts[disorder]["Brain"],
        "gut_abstracts":    abstract_counts[disorder]["Gut"],
        "total_abstracts":  abstract_counts[disorder]["Brain"] + abstract_counts[disorder]["Gut"],
        "brain_genes":      len(bg),
        "gut_genes":        len(gg),
        "shared_genes":     len(shared),
        "brain_gene_list":  ", ".join(sorted(bg)),
        "gut_gene_list":    ", ".join(sorted(gg)),
        "shared_gene_list": ", ".join(sorted(shared)),
    })

summary = pd.DataFrame(summary_rows)
summary.to_csv(os.path.join(_OUTPUT_DIR, "summary_stats.csv"), index=False)
print(f"Saved summary_stats.csv")
print(summary[["disorder","brain_abstracts","gut_abstracts","brain_genes","gut_genes","shared_genes"]].to_string(index=False))

legend_rows = []
for _, s in summary.iterrows():
    legend_rows.append({
        "Disorder":               s["disorder"],
        "Brain abstracts":        s["brain_abstracts"],
        "Gut abstracts":          s["gut_abstracts"],
        "Total abstracts":        s["total_abstracts"],
        "Unique brain genes":     s["brain_genes"],
        "Unique gut genes":       s["gut_genes"],
        "Genes in both (shared)": s["shared_genes"],
    })
# Totals row
legend_rows.append({
    "Disorder":               "TOTAL (all disorders)",
    "Brain abstracts":        sum(r["Brain abstracts"] for r in legend_rows),
    "Gut abstracts":          sum(r["Gut abstracts"] for r in legend_rows),
    "Total abstracts":        sum(r["Total abstracts"] for r in legend_rows),
    "Unique brain genes":     gene_table[gene_table["tissue"].isin(["Brain","Both"])]["gene"].nunique(),
    "Unique gut genes":       gene_table[gene_table["tissue"].isin(["Gut","Both"])]["gene"].nunique(),
    "Genes in both (shared)": gene_table[gene_table["tissue"] == "Both"]["gene"].nunique(),
})
legend_df = pd.DataFrame(legend_rows)

with pd.ExcelWriter(os.path.join(_OUTPUT_DIR, "gene_report.xlsx"), engine="openpyxl") as writer:
    gene_table.to_excel(writer, sheet_name="Genes", index=False)
    legend_df.to_excel(writer, sheet_name="Summary", index=False)
print("Saved gene_report.xlsx")


n_disorders = len(DISORDER_ORDER)
ncols = 4
nrows = (n_disorders + 1 + ncols - 1) // ncols  # +1 for overall

fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4.5, nrows * 4))
axes = axes.flatten()

BRAIN_COLOR = "#4C9BE8"
GUT_COLOR   = "#E87B4C"

def draw_venn(ax, brain_set, gut_set, title):
    shared = brain_set & gut_set
    brain_only = brain_set - gut_set
    gut_only   = gut_set - brain_set

    if not brain_set and not gut_set:
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes, fontsize=10, color="gray")
        ax.set_title(title, fontsize=11, fontweight="bold", pad=8)
        ax.axis("off")
        return

    # Draw even if one side is empty
    v = venn2(
        subsets=(len(brain_only), len(gut_only), len(shared)),
        set_labels=("Brain", "Gut"),
        ax=ax,
        set_colors=(BRAIN_COLOR, GUT_COLOR),
        alpha=0.6,
    )
    # Label shared genes if few enough
    if v.get_label_by_id("11") and shared:
        genes_str = "\n".join(sorted(shared)) if len(shared) <= 8 else f"{len(shared)} genes"
        v.get_label_by_id("11").set_text(genes_str)
        v.get_label_by_id("11").set_fontsize(7)

    ax.set_title(title, fontsize=11, fontweight="bold", pad=8)

# Per disorder venns
for i, disorder in enumerate(DISORDER_ORDER):
    draw_venn(axes[i],
              brain_genes_by_disorder[disorder],
              gut_genes_by_disorder[disorder],
              disorder)

# Overall venn
draw_venn(axes[n_disorders], all_brain, all_gut, "Overall (All Disorders)")

# Hide unused axes
for j in range(n_disorders + 1, len(axes)):
    axes[j].axis("off")

fig.suptitle("Brain vs Gut Gene Overlap by Disorder\n(genes expressed in disorder models, 2000–2026)",
             fontsize=14, fontweight="bold", y=1.01)
plt.tight_layout()
fig.savefig(os.path.join(_OUTPUT_DIR, "venn_diagrams.png"), dpi=150, bbox_inches="tight")
print("Saved venn_diagrams.png")
