"""Eval harness: old regex qualification gate vs. new LLM gate vs. hand labels.

Proves the accuracy claim BEFORE committing to the multi-day full-scale run.

Workflow:
  1. Run scanner.py once with the regex gate (USE_LLM_QUALIFIER unset) to produce
     model_organisms_dataset.csv.
  2. python eval_qualify.py --sample 200
        -> writes gold_labels.csv with `label` blank (stratified: half regex-pass,
           half regex-fail, balanced across the 10 disorders).
  3. Hand-fill the `label` column (TRUE / FALSE) per the inclusion criteria in the
     plan (esp. the differential disease-vs-control bar).
  4. Start the LLM server (llama-server), then:
     python eval_qualify.py
        -> confusion matrix + precision/recall/F1 + Cohen's kappa for both gates
           vs. labels, plus the disagreement list for manual review.

Self-contained: the old gate's verdict is reused from model_organisms_dataset.csv
(`gene_context` column) so scanner.py is NOT imported (its module-level pipeline
would otherwise run). Only llm_qualify is imported, to run the new gate live.
"""

import sys
import argparse
from collections import Counter

import pandas as pd

try:
    from llm_qualify import qualify_batch_sync
except ImportError as e:
    print(f"llm_qualify import failed ({e}); install httpx and run the LLM server first.")
    sys.exit(1)


def parse_genes(cell):
    if not isinstance(cell, str) or cell.strip() in ("", "nan"):
        return []
    return [g.strip() for g in cell.split(",") if g.strip() and g.strip() != "nan"]


def to_bool(v):
    if isinstance(v, bool):
        return v
    if isinstance(v, str):
        return v.strip().lower() in ("true", "1", "yes", "t")
    return bool(v)


def _strat_sample(df_sub, n_per_disorder):
    if df_sub.empty:
        return df_sub.iloc[:0]
    parts = []
    for _d, g in df_sub.groupby("disorder"):
        parts.append(g.sample(min(len(g), n_per_disorder), random_state=1))
    return pd.concat(parts)


def sample_gold(n, src="model_organisms_dataset.csv", out="gold_labels.csv"):
    try:
        df = pd.read_csv(src, low_memory=False)
    except FileNotFoundError:
        sys.exit(f"{src} not found — run scanner.py once (regex gate) first.")
    df = df.dropna(subset=["abstract"]).copy()
    for c in ("disorder", "gene_context", "genes_found"):
        if c not in df.columns:
            sys.exit(f"{src} missing required column `{c}`.")
    df["gc"] = df["gene_context"].apply(to_bool)
    n_dis = max(1, df["disorder"].nunique())
    half = max(1, n // 2)
    per = max(1, half // n_dis)
    sub = pd.concat([_strat_sample(df[df["gc"]], per),
                     _strat_sample(df[~df["gc"]], per)])
    sub = sub.sample(frac=1, random_state=1).head(n).reset_index(drop=True)
    out_df = sub[["disorder", "abstract", "gene_context", "genes_found"]].copy()
    out_df.columns = ["disorder", "body", "regex_gate", "genes_hint"]
    out_df["label"] = ""   # user fills TRUE / FALSE
    out_df.to_csv(out, index=False)
    npos = (out_df["regex_gate"].astype(str).str.lower() == "true").sum()
    print(f"Wrote {out} — {len(out_df)} rows ({npos} regex-pass, "
          f"{len(out_df)-npos} regex-fail). Fill the `label` column, then re-run without --sample.")


def cohen_kappa(y1, y2):
    n = len(y1)
    if n == 0:
        return float("nan")
    po = sum(a == b for a, b in zip(y1, y2)) / n
    c1, c2 = Counter(y1), Counter(y2)
    pe = sum((c1[v] / n) * (c2[v] / n) for v in set(y1) | set(y2))
    if pe >= 1:
        return float("nan")
    return (po - pe) / (1 - pe)


def metrics(preds, labels):
    tp = sum(p and l for p, l in zip(preds, labels))
    fp = sum(p and not l for p, l in zip(preds, labels))
    tn = sum(not p and not l for p, l in zip(preds, labels))
    fn = sum(not p and l for p, l in zip(preds, labels))
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    acc = (tp + tn) / len(preds) if preds else 0.0
    return {"TP": tp, "FP": fp, "TN": tn, "FN": fn,
            "precision": prec, "recall": rec, "f1": f1, "accuracy": acc}


def _print_metrics(name, preds, labels):
    m = metrics(preds, labels)
    k = cohen_kappa(preds, labels)
    print(f"\n{name} vs. labels:")
    print(f"  TP={m['TP']} FP={m['FP']} TN={m['TN']} FN={m['FN']}")
    print(f"  precision={m['precision']:.3f}  recall={m['recall']:.3f}  "
          f"f1={m['f1']:.3f}  accuracy={m['accuracy']:.3f}  kappa={k:.3f}")
    return m, k


def run_eval(gold="gold_labels.csv"):
    try:
        df = pd.read_csv(gold, dtype=str)
    except FileNotFoundError:
        sys.exit(f"{gold} not found — run `eval_qualify.py --sample N` first.")
    for c in ("disorder", "body", "regex_gate", "genes_hint", "label"):
        if c not in df.columns:
            sys.exit(f"{gold} missing column `{c}`.")
    df = df.dropna(subset=["label"]).copy()
    df = df[df["label"].str.strip() != ""]
    if df.empty:
        sys.exit("No labeled rows — fill the `label` column (TRUE/FALSE) and re-run.")
    df = df.reset_index(drop=True)
    labels = df["label"].apply(to_bool).tolist()
    regex_preds = df["regex_gate"].apply(to_bool).tolist()

    items = [(str(i), row["body"], parse_genes(row["genes_hint"]))
             for i, row in df.iterrows()]
    print(f"Running LLM gate on {len(items)} labeled abstracts …")
    judgments = qualify_batch_sync(items)
    llm_preds = [(str(i) in judgments) and judgments[str(i)].qualifies
                 for i in range(len(df))]
    df["llm_gate"] = llm_preds
    df["llm_reason"] = [judgments.get(str(i)).reason if judgments.get(str(i)) else ""
                       for i in range(len(df))]

    print("\n" + "=" * 60)
    print(f"EVAL — {len(df)} labeled abstracts")
    print("=" * 60)
    rm, rk = _print_metrics("REGEX gate ", regex_preds, labels)
    lm, lk = _print_metrics("LLM gate   ", llm_preds, labels)

    print("\nDecision:")
    if lk >= rk and lm["precision"] >= rm["precision"]:
        print(f"  SHIP the LLM gate (kappa {lk:.3f} >= {rk:.3f}, "
              f"precision {lm['precision']:.3f} >= {rm['precision']:.3f}).")
    else:
        print(f"  DO NOT ship yet — iterate the rubric/few-shots in llm_qualify.py "
              f"(bump PROMPT_VERSION) and re-run.")

    print("\nDisagreements (regex != LLM, or either != label):")
    shown = 0
    for i, row in df.iterrows():
        if regex_preds[i] != llm_preds[i] or llm_preds[i] != labels[i] or regex_preds[i] != labels[i]:
            shown += 1
            snip = " ".join(str(row["body"]).split())[:160]
            print(f"\n  [{row['disorder']}] regex={regex_preds[i]} llm={llm_preds[i]} "
                  f"label={labels[i]}")
            print(f"    llm_reason: {row['llm_reason']}")
            print(f"    body: {snip}…")
            if shown >= 40:
                print("  …(truncated, see gold_labels.csv for the rest)")
                break


def main():
    ap = argparse.ArgumentParser(description="Eval regex-vs-LLM qualification gate.")
    ap.add_argument("--sample", type=int, metavar="N",
                    help="Write gold_labels.csv with N stratified rows (label blank) "
                         "drawn from model_organisms_dataset.csv.")
    ap.add_argument("--gold", default="gold_labels.csv",
                    help="Path to the labeled gold CSV (default: gold_labels.csv).")
    args = ap.parse_args()
    if args.sample:
        sample_gold(args.sample)
    else:
        run_eval(args.gold)


if __name__ == "__main__":
    main()
