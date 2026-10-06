# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project does

A text-mining pipeline that fetches PubMed abstracts for 10 psychiatric and neurodegenerative disorders (Schizophrenia, Bipolar, Autism, ADHD, Major Depression, Anxiety, PTSD, OCD, Alzheimer's, Parkinson's), then identifies genes mentioned in a genuine gene-expression context and classifies each by tissue (brain vs. gut). Output is a consolidated, disorder-organized gene table plus Venn diagrams of brain/gut gene overlap.

## Running the pipeline

Each script is a top-level script with module-level execution code — there is no `if __name__ == "__main__"` guard, so `python <script>.py` runs the whole pipeline. Run in order:

```bash
python scanner.py              # fetch PubMed + filter + classify → model_organisms_dataset.csv, *.txt, passing_abstracts.csv, per-disorder *.png
python generate_report.py     # reads model_organisms_dataset.csv → consolidated_gene_table.csv, summary_stats.csv, gene_report.xlsx, venn_diagrams.png
python generate_gene_table.py # reads model_organisms_dataset.csv → gene_tissue_table.csv
```

[scanner_debug.py](scanner_debug.py) is a fast smoke-test variant (single disorder, 20 abstracts, year 2018, fixed candidate gene lists) for iterating on regex/filter tuning without a full PubMed run.

Dependencies (no `requirements.txt` — install manually): `biopython`, `pandas`, `matplotlib`, `matplotlib-venn`, `openpyxl`. The optional LLM qualification gate additionally needs `httpx`.

[scanner.py](scanner.py) requires network access to PubMed and genenames.org. On first run it downloads the HGNC gene list to [hgnc_full.txt](hgnc_full.txt) and caches it. Set the `NCBI_API_KEY` environment variable to raise the Entrez rate limit (sleep drops from 0.34 s to 0.11 s between fetches). The Entrez contact email is hardcoded as `dasfour@outlook.com`, and SSL verification is globally disabled (`ssl._create_default_https_context = ssl._create_unverified_context`). Full-run knobs are `DEBUG_YEAR_START` (default 2000) and `DEBUG_PER_YEAR` (default `None` = no cap) near the bottom of the fetch module.

## Optional: local-LLM qualification gate (hybrid)

The contextual qualification decision inside `detect_gene_context` (whether an abstract reports *differential* gene expression in brain/gut tissue of a real disease model, vs a control) can instead be made by a small local LLM. This is a **hybrid**: the LLM replaces ONLY that gate; HGNC gene extraction, per-gene tissue classification, regulation, and multi-disorder attribution stay as the existing regex (the LLM would be worse there — it hallucinates gene symbols).

- **`USE_LLM_QUALIFIER=1`** on [scanner.py](scanner.py) enables the LLM gate; unset (default) runs the original regex `detect_gene_context`. On import failure (no `httpx` / server down) scanner prints a warning and falls back to the regex gate.
- **`USE_LLM_GATE=1`** on [scanner_debug.py](scanner_debug.py) batch-qualifies the 20 smoke abstracts and prints each `Judgment`, then exits; default-off still runs the verbose sentence trace.
- [llm_qualify.py](llm_qualify.py) is the client: async batched `httpx` calls to a local llama.cpp `llama-server` (OpenAI-compatible `/v1/chat/completions`, `response_format=json_object`, `temperature=0`, `max_tokens=96`). It writes a resume-safe SQLite cache `llm_cache.sqlite` keyed by `sha256(MODEL_ID|PROMPT_VERSION|item_id)`; the gate is disorder-agnostic so a PMID's judgment is reused across all 10 disorder fetches. Bump `PROMPT_VERSION` to re-judge after editing the rubric/few-shots.
- Env knobs (all optional): `LLM_MODEL` (default `qwen2.5-3b-instruct-q4_k_m`), `LLM_BASE_URL` (default `http://127.0.0.1:8080`), `LLM_CONCURRENCY` (default 8), `LLM_CACHE_PATH`, `LLM_PROMPT_VERSION`. Start the server e.g. `llama-server -m Qwen2.5-3B-Instruct-Q4_K_M.gguf -c 4096 -ngl 99 -np 8 --cont-batching --host 127.0.0.1 --port 8080` (`-np` is the throughput knob).
- [eval_qualify.py](eval_qualify.py) proves accuracy before the multi-hour full run: `python eval_qualify.py --sample 200` writes a stratified `gold_labels.csv` (drawn from `model_organisms_dataset.csv` — half regex-pass, half regex-fail, balanced across disorders) with the `label` column blank; hand-fill it, start the server, then `python eval_qualify.py` prints confusion matrix + precision/recall/F1 + Cohen's kappa for both gates vs labels plus the disagreement list. Ship rule: LLM gate only if `kappa_new >= kappa_old` AND `precision_new >= precision_old`.

## Architecture

### Data flow

[scanner.py](scanner.py) is the core (~2100 lines). It fetches abstracts, filters them, detects genes, and writes the primary intermediate `model_organisms_dataset.csv` — one row per abstract×disorder with columns for tissue, gene_context, brain/gut/shared gene lists, and the regulation dict. The downstream scripts ([generate_report.py](generate_report.py), [generate_gene_table.py](generate_gene_table.py)) are pure post-processing of that CSV; they re-run nothing live.

`model_organisms_dataset.csv` is not committed — only downstream outputs ([passing_abstracts.csv](passing_abstracts.csv), [consolidated_gene_table.csv](consolidated_gene_table.csv), [gene_report.xlsx](gene_report.xlsx), [venn_diagrams.png](venn_diagrams.png)) and the [hgnc_full.txt](hgnc_full.txt) cache are present. [reformat.py](reformat.py) reads `*_genes.csv` files that [scanner.py](scanner.py) does not produce (scanner emits `*_genes_context.txt`), so it is currently stale/unused.

### scanner.py internal structure (single file, top-to-bottom)

1. **Config & data loading** (~lines 1–153): HGNC loader (`_load_hgnc`) validates gene symbols and resolves aliases → approved symbol; `Entrez` setup; the `DISORDERS` dict mapping each disorder to a PubMed/MeSH query; `RODENT_SPECIES`/`SPECIES_KEYWORDS`/`BRAIN_KEYWORDS`/`GUT_KEYWORDS` lists.

2. **Filter regexes** (~156–968): a large set of compiled patterns plus `_abstract_*` predicate functions deciding whether an abstract qualifies. Two kinds:
   - **Exclusion filters** reject abstracts whose gene changes stem from a confound: therapeutic drugs (`_abstract_uses_therapeutic_drug`), model-inducer compounds (`_abstract_uses_model_inducer`), non-drug interventions (`_abstract_uses_nondrug_intervention`), epigenetics/miRNA (`_abstract_uses_epigenetics`), behavioral-primary studies, genetic-association/GWAS, computational/reuse-of-public-data, fluid-only measurements, blood-derived cell lines, and LTD-vs-Major-Depression disambiguation.
   - **Positive requirements** an abstract must also meet: a real disease model (`_has_disease_model` — patient tissue or animal model) and an explicit expression-measurement signal (`_EXPRESSION_MEASUREMENT_RE`) found in the results section (`_get_results_section`).

3. **Gene detection** (~980–1694): two paths in `detect_mentioned_genes`:
   - Alias / common-name path via `_TIER2_ALIASES` (e.g. "occludin" → OCLN, "tau" → MAPT) — context-independent.
   - Symbol path: any ALL-CAPS token (regex `[A-Z][A-Z0-9]{1,8}`) that is an HGNC-approved symbol or alias, appears in a sentence with biological vocabulary (`_BIO_CONTEXT_RE`), and is not in the `_NON_GENE_SYMBOLS` blocklist. A `BANNED_ACRONYM_CONTEXT` map handles ambiguous tokens (e.g. `SI`, `TH`) by requiring disqualifying nearby words.

   `detect_tissue_for_gene` classifies each gene per-sentence with results-section priority (background/intro mentions are ignored when the results section places the gene in a tissue); `detect_gene_regulation` labels Up/Down/Mixed from nearby regulation verbs; `detect_gene_context` applies proximity rules (1–3 words, with allowed `_LINKING_WORDS` in between) between a regulation term and "gene"/a detected gene symbol — this is the gate that marks an abstract as `gene_context=True`.

4. **Fetch + main loop** (~1696–1960): `fetch_abstracts` does chronological year-by-year PubMed `esearch`/`efetch`, with retries and PMID dedup. The main loop runs two passes: Pass 1 builds per-disorder brain/gut gene pools and tissue/gene-context tallies; Pass 2 annotates each row (brain/gut/shared gene lists, regulation dict) using the completed pools. **Multi-disorder attribution**: a passing abstract is duplicated into every disorder it mentions (`detect_disorders_in_abstract`), so the same abstract can appear under multiple disorders. Under `USE_LLM_QUALIFIER=1`, Pass 1 is split into three phases — **A** (cheap precompute: tissue, genes, tallies; no LLM), **B** (one batched `qualify_batch_sync` call per disorder via [llm_qualify.py](llm_qualify.py)), **C** (apply each `Judgment.qualifies` as `has_gene_context`, then the unchanged tissue tally / per-gene block / tissue-gene downgrade safety net / multi-disorder attribution) — preserving the exact tuple shape Pass 2 consumes. `detect_gene_context` is retained (not deleted) as the regex fallback and the eval baseline.

5. **Outputs & plots** (~1960–end): CSV/text outputs, per-disorder brain-vs-gut Venn PNGs, and console reports of disorder pairs sharing the most genes.

### Tuning the filters

Most adjustments are made by editing the regex pattern lists or the `_abstract_*` predicate functions in [scanner.py](scanner.py), not the downstream scripts. When a false positive/negative appears, identify which filter admits or rejects it, then smoke-test against [scanner_debug.py](scanner_debug.py) before a full run. Note `scanner_debug.py` uses a *fixed* `BRAIN_GENES`/`GUT_GENES` candidate list and an `ALL_CANDIDATE_GENES` approach, whereas production [scanner.py](scanner.py) uses *dynamic* HGNC-validated detection — the two can diverge, so confirm final filter behavior in `scanner.py` itself.
