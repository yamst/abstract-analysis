# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project does

A text-mining pipeline that fetches PubMed abstracts for 10 psychiatric and neurodegenerative disorders (Schizophrenia, Bipolar, Autism, ADHD, Major Depression, Anxiety, PTSD, OCD, Alzheimer's, Parkinson's), then identifies genes mentioned in a genuine gene-expression context and classifies each by tissue (brain vs. gut). Output is a consolidated, disorder-organized gene table plus Venn diagrams of brain/gut gene overlap.

## Running the pipeline

Scripts are in `src/` and run as top-level modules with module-level execution code (no `if __name__ == "__main__"` guard). Run in order:

```bash
python src/scanner.py              # fetch PubMed + filter + classify → model_organisms_dataset.csv, *.txt, passing_abstracts.csv, per-disorder *.png
python src/generate_report.py      # reads model_organisms_dataset.csv → consolidated_gene_table.csv, summary_stats.csv, gene_report.xlsx, venn_diagrams.png
python src/generate_gene_table.py  # reads model_organisms_dataset.csv → gene_tissue_table.csv
```

All outputs go to `output/`.

`src/scanner_debug.py` is a fast smoke-test (single disorder, 20 abstracts, year 2018, fixed candidate gene lists) for iterating on regex/filter tuning without a full PubMed run.

### Dependencies

Install from `requirements.txt`:

```bash
pip install -r requirements.txt
```

Packages: `biopython`, `pandas`, `matplotlib`, `matplotlib-venn`, `openpyxl`, `httpx`.

### Environment variables

| Variable | Default | Description |
|----------|---------|-------------|
| `NCBI_API_KEY` | - | NCBI Entrez API key (raises rate limit: 0.11s vs 0.34s between fetches) |
| `DEBUG_YEAR_START` | 2000 | Start year for PubMed search |
| `DEBUG_PER_YEAR` | None | Cap abstracts per year (for testing) |
| `USE_LLM_QUALIFIER` | 0 | Enable LLM gate (1=on) |
| `LLM_MODEL` | qwen2.5:3b | Ollama model name |
| `LLM_BASE_URL` | http://127.0.0.1:11434 | Ollama server URL |

The Entrez contact email is hardcoded as `dasfour@outlook.com`, and SSL verification is globally disabled.

## Optional: LLM qualification gate (hybrid)

The qualification decision inside `detect_gene_context` (whether an abstract reports *differential* gene expression in brain/gut tissue of a real disease model, vs a control) can instead be made by a local LLM. This is a **hybrid**: the LLM replaces ONLY that gate; HGNC gene extraction, per-gene tissue classification, regulation, and multi-disorder attribution stay as regex (the LLM would be worse there — it hallucinates gene symbols).

### Setup (one-time)

```bash
python setup/install_ollama.py
```

This downloads Ollama Windows build with CUDA DLLs and the Qwen2.5-3B-Instruct Q4_K_M model (~2 GB).

### Run with LLM gate

```bash
# Terminal 1: Start the server
ollama/ollama.exe serve

# Terminal 2: Run scanner
set USE_LLM_QUALIFIER=1
set LLM_BASE_URL=http://127.0.0.1:11434
set LLM_MODEL=qwen2.5:3b-instruct-q4_k_m
python src/scanner.py
```

### How it works

- `src/llm_qualify.py` is the client: async batched calls to Ollama's OpenAI-compatible endpoint (`response_format=json_object`, `temperature=0`, `max_tokens=96`).
- SQLite cache (`llm_cache.sqlite`) keyed by `sha256(MODEL_ID|PROMPT_VERSION|item_id)`; a PMID's judgment is reused across all 10 disorder fetches.
- `src/eval_qualify.py` proves accuracy before a full run: `python src/eval_qualify.py --sample 200` creates a stratified `gold_labels.csv` for hand-labeling; running again prints confusion matrix + precision/recall/F1 + Cohen's kappa for both gates. Ship rule: LLM gate only if `kappa_new >= kappa_old` AND `precision_new >= precision_old`.

## Architecture

### Data flow

`src/scanner.py` is the core (~2100 lines). It fetches abstracts, filters them, detects genes, and writes `output/model_organisms_dataset.csv` — one row per abstract×disorder with columns for tissue, gene_context, brain/gut/shared gene lists, and regulation dict. Downstream scripts (`src/generate_report.py`, `src/generate_gene_table.py`) are pure post-processing of that CSV; they re-run nothing live.

`output/model_organisms_dataset.csv` is not committed — only downstream outputs (`passing_abstracts.csv`, `consolidated_gene_table.csv`, `gene_report.xlsx`, `venn_diagrams.png`) and `hgnc_full.txt` cache are present. `src/reformat.py` reads `*_genes.csv` files that scanner does not produce (scanner emits `*_genes_context.txt`), so it is currently stale/unused.

### scanner.py internal structure (single file, top-to-bottom)

1. **Config & data loading** (~lines 1–153): HGNC loader (`_load_hgnc`) validates gene symbols and resolves aliases → approved symbol; `Entrez` setup; the `DISORDERS` dict mapping each disorder to a PubMed/MeSH query; `RODENT_SPECIES`/`SPECIES_KEYWORDS`/`BRAIN_KEYWORDS`/`GUT_KEYWORDS` lists.

2. **Filter regexes** (~156–968): compiled patterns plus `_abstract_*` predicate functions. Two kinds:
   - **Exclusion filters** reject abstracts whose gene changes stem from a confound: therapeutic drugs, model-inducer compounds, non-drug interventions, epigenetics/miRNA, behavioral-primary studies, genetic-association/GWAS, computational/reuse-of-public-data, fluid-only measurements, blood-derived cell lines, and LTD-vs-Major-Depression disambiguation.
   - **Positive requirements**: a real disease model (`_has_disease_model` — patient tissue or animal model) and an explicit expression-measurement signal (`_EXPRESSION_MEASUREMENT_RE`) found in the results section.

3. **Gene detection** (~980–1694): two paths in `detect_mentioned_genes`:
   - Alias/common-name path via `_TIER2_ALIASES` (e.g. "occludin" → OCLN, "tau" → MAPT) — context-independent.
   - Symbol path: any ALL-CAPS token that is an HGNC-approved symbol or alias, appears in a sentence with biological vocabulary, and is not in the `_NON_GENE_SYMBOLS` blocklist. `BANNED_ACRONYM_CONTEXT` handles ambiguous tokens by requiring disqualifying nearby words.

   `detect_tissue_for_gene` classifies each gene per-sentence with results-section priority; `detect_gene_regulation` labels Up/Down/Mixed; `detect_gene_context` applies proximity rules between a regulation term and "gene"/a detected gene symbol — this is the gate marking an abstract as `gene_context=True`.

4. **Fetch + main loop** (~1696–1960): `fetch_abstracts` does year-by-year PubMed `esearch`/`efetch` with retries and PMID dedup. The main loop runs two passes: Pass 1 builds per-disorder brain/gut gene pools and tissue/gene-context tallies; Pass 2 annotates each row using completed pools. **Multi-disorder attribution**: a passing abstract is duplicated into every disorder it mentions. Under `USE_LLM_QUALIFIER=1`, Pass 1 splits into phases — **A** (precompute: tissue, genes, tallies), **B** (batched LLM calls), **C** (apply judgments, then unchanged downstream logic).

5. **Outputs & plots** (~1960–end): CSV/text outputs to `output/`, per-disorder Venn PNGs, console reports.

### Tuning the filters

Edit regex pattern lists or `_abstract_*` predicates in `src/scanner.py`. Smoke-test with `src/scanner_debug.py` before a full run. Note: `scanner_debug.py` uses fixed candidate gene lists; production uses dynamic HGNC-validated detection. Confirm final behavior in `scanner.py` itself.

## Known issues

See [KNOWN_ISSUES.md](KNOWN_ISSUES.md) for current limitations. Notably: multi-disorder attribution is too permissive (abstracts attributed to ALL mentioned disorders, not just studied ones).

## Tools directory

`tools/` contains diagnostic and benchmark scripts:
- `tools/smoke_llm_gate.py` — quick LLM gate test
- `tools/bench_throughput.py` — benchmark LLM throughput
- `tools/diag_*.py` — diagnostics for context, DLL issues, imports
