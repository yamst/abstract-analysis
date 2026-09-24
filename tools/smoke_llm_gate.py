"""Smoke: run the LLM gate over real regex-passing abstracts (no PubMed fetch).

Reads passing_abstracts.csv (abstracts the regex gate already admitted), samples
N, and runs them through llm_qualify.qualify_batch_sync. Prints each Judgment
alongside the regex-detected genes, so you can eyeball whether the LLM agrees with
the regex gate or tightens it (flags over-inclusions the regex admitted).

Needs the LLM server running (llama-cpp-python / llama-server / LM Studio) at
LLM_BASE_URL (default http://127.0.0.1:8080).
"""

import sys
import argparse
import pandas as pd

from llm_qualify import qualify_batch_sync

# Windows console defaults to cp1252 and crashes on Greek letters (e.g. the β
# in "Aβ"); force UTF-8 so LLM reasons print cleanly without PYTHONUTF8=1.
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


def parse_genes(cell):
    out = []
    if isinstance(cell, str):
        out = [g.strip() for g in cell.split(",") if g.strip() and g.strip() != "nan"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=12, help="number of abstracts to sample")
    ap.add_argument("--csv", default="passing_abstracts.csv")
    args = ap.parse_args()

    df = pd.read_csv(args.csv, low_memory=False)
    n = min(args.n, len(df))
    sample = df.sample(n=n, random_state=1) if n < len(df) else df

    items = []
    for i, (_, row) in enumerate(sample.iterrows()):
        body = str(row.get("abstract", "")).strip()
        mentioned = sorted(set(
            parse_genes(row.get("brain_genes_found"))
            + parse_genes(row.get("gut_genes_found"))
            + parse_genes(row.get("shared_genes_found"))
        ))
        items.append((f"smoke-{i}", body, mentioned))

    print(f"Qualifying {len(items)} abstracts …\n")
    judgments = qualify_batch_sync(items)

    agree = 0
    for i, (_, row) in enumerate(sample.iterrows()):
        iid = f"smoke-{i}"
        j = judgments.get(iid)
        if j is None:
            print(f"[{row['disorder']}] (no judgment)\n")
            continue
        if j.qualifies:
            agree += 1
        snip = " ".join(str(row.get("abstract", "")).split())[:120]
        genes = ", ".join(items[i][2]) or "(none)"
        print(f"[{row['disorder']}] qualifies={j.qualifies}  tissue={j.tissue_hint}")
        print(f"  regex genes: {genes}")
        print(f"  llm reason  : {j.reason}")
        print(f"  abstract    : {snip}…\n")

    print("=" * 50)
    print(f"LLM agrees with regex on {agree}/{len(items)} "
          f"(regex admitted all {len(items)}; LLM rejected {len(items)-agree}).")


if __name__ == "__main__":
    main()
