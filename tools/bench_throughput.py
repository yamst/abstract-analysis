"""Measure LLM-gate throughput on the local server and extrapolate to full scale.

Times a fresh batch of N uncached abstracts (ids prefixed bench- so they never hit
the smoke cache), prints abstracts/min, and projects the wall-clock for 50k/100k/
250k abstracts — the numbers that decide whether a full scan is feasible on this
GPU and whether Qwen2.5-3B is the right size.
"""
import time
import sys
import argparse

import pandas as pd

from llm_qualify import qualify_batch_sync, LLM_CONCURRENCY, LLM_BASE_URL

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


def parse_genes(cell):
    if not isinstance(cell, str):
        return []
    return [g.strip() for g in cell.split(",") if g.strip() and g.strip() != "nan"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=20)
    ap.add_argument("--csv", default="passing_abstracts.csv")
    args = ap.parse_args()

    df = pd.read_csv(args.csv, low_memory=False)
    n = min(args.n, len(df))
    # random_state=42 (different from smoke's 1) + bench- prefix => uncached
    sample = df.sample(n=n, random_state=42) if n < len(df) else df
    items = []
    for i, (_, row) in enumerate(sample.iterrows()):
        body = str(row.get("abstract", "")).strip()
        mentioned = sorted(set(
            parse_genes(row.get("brain_genes_found"))
            + parse_genes(row.get("gut_genes_found"))
            + parse_genes(row.get("shared_genes_found"))
        ))
        items.append((f"bench-{i}", body, mentioned))

    print(f"base_url={LLM_BASE_URL}  concurrency={LLM_CONCURRENCY}  n={len(items)}")
    t0 = time.time()
    judgments = qualify_batch_sync(items)
    dt = time.time() - t0

    rate = len(items) / dt
    nq = sum(1 for v in judgments.values() if v.qualifies)
    errs = sum(1 for v in judgments.values() if "error" in v.reason or v.reason in ("empty_response", "json_parse_failed"))
    print(f"elapsed={dt:.1f}s  rate={rate:.2f} abstracts/s = {rate*60:.0f} abstracts/min")
    print(f"qualifies={nq}/{len(items)}  errors={errs}")
    print("extrapolation (single-GPU, this concurrency):")
    for total in (50000, 100000, 250000):
        hrs = total / rate / 3600
        print(f"  {total:>7,} abstracts -> {total/rate/60:,.0f} min = {hrs:.1f} h")


if __name__ == "__main__":
    main()
