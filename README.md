# Abstract Analysis Pipeline

![PubMed](https://img.shields.io/badge/PubMed-326599?logo=pubmed&logoColor=white)
![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=white)
![Windows](https://img.shields.io/badge/Windows-0078D6?logo=data:image/svg%2Bxml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCA0NDggNTEyIj48cGF0aCBmaWxsPSIjZmZmIiBkPSJNMCAxNzZoMjAwdjE2MEgwVjE3NnptMjQ4IDBoMjAwdjE2MEgyNDhWMTc2ek0wIDM1MmgyMDB2MTYwSDBWMzUyem0yNDggMGgyMDB2MTYwSDI0OFYzNTJ6Ii8+PC9zdmc+&logoColor=white)
![Ollama](https://img.shields.io/badge/Ollama-000000?logo=ollama&logoColor=white)
![SQLite](https://img.shields.io/badge/SQLite-003B57?logo=sqlite&logoColor=white)
![Pandas](https://img.shields.io/badge/Pandas-150458?logo=pandas&logoColor=white)

A text-mining pipeline that fetches PubMed abstracts for 10 psychiatric and neurodegenerative disorders, identifies genes mentioned in a gene-expression context, and classifies each by tissue (brain vs. gut).

## 📋 Disorders Covered

Schizophrenia, Bipolar, Autism, ADHD, Major Depression, Anxiety, PTSD, OCD, Alzheimer's, Parkinson's

## 🚀 Quick Start

### 📦 Prerequisites

Install Python dependencies:

```bash
pip install -r requirements.txt
```

### 🔑 Set NCBI API Key (Recommended)

Set your NCBI API key to increase rate limits:

```bash
set NCBI_API_KEY=your_key_here
```

### ▶️ Run the Pipeline

```bash
python src/scanner.py              # Fetch PubMed + filter + classify
python src/generate_report.py      # Generate summary and Venn diagrams
python src/generate_gene_table.py # Generate consolidated gene table
```

All outputs are saved to `output/`.

### 🤖 Optional: LLM Qualification Gate

The LLM qualification gate replaces the regex-based qualification decision with a local LLM judgment. This is a **hybrid**: the LLM only decides qualification; gene extraction, tissue classification, and regulation detection remain regex-based.

#### Setup (one-time after fresh clone)

```bash
python setup/install_ollama.py
```

This downloads everything needed:
- Ollama Windows build with CUDA DLLs
- Qwen2.5-3B-Instruct Q4_K_M model (~2 GB)

#### Run with LLM Gate

```bash
# Terminal 1: Start the server
ollama\ollama.exe serve

# Terminal 2: Run scanner
set USE_LLM_QUALIFIER=1
set LLM_BASE_URL=http://127.0.0.1:11434
set LLM_MODEL=qwen2.5-3b-instruct-q4_k_m
python src/scanner.py
```

## 📁 Output Files

All outputs are in `output/`:

### scanner.py
| File | Description |
|------|-------------|
| `model_organisms_dataset.csv` | Primary intermediate: one row per abstract×disorder with all metadata |
| `passing_abstracts.csv` | Abstracts that passed all filters (gene_context=True), with year extracted |
| `abstracts_brain_genes_context.txt` | Full abstract text for abstracts with brain genes, formatted for review |
| `abstracts_gut_genes_context.txt` | Full abstract text for abstracts with gut genes, formatted for review |
| `abstracts_shared_genes_context.txt` | Full abstract text for abstracts with shared brain/gut genes |
| `shared_genes_summary.txt` | Count per disorder of shared-gene abstracts |

### generate_report.py
| File | Description |
|------|-------------|
| `consolidated_gene_table.csv` | One row per gene×disorder with tissue, paper counts, regulation consensus |
| `summary_stats.csv` | Per-disorder counts: abstracts, unique genes, shared genes, gene lists |
| `gene_report.xlsx` | Excel workbook: "Genes" sheet (full table) + "Summary" sheet |
| `venn_diagrams.png` | 11-panel figure: per-disorder + overall brain/gut gene overlap |

### generate_gene_table.py
| File | Description |
|------|-------------|
| `gene_tissue_table.csv` | Flat table: subject sentence, tissue, disorder, gene (one row per gene mention) |

## 🔧 Configuration

Environment variables:

| Variable | Default | Description |
|----------|---------|-------------|
| `NCBI_API_KEY` | - | NCBI Entrez API key (faster fetches) |
| `DEBUG_YEAR_START` | 2000 | Start year for PubMed search |
| `DEBUG_PER_YEAR` | None | Cap abstracts per year (for testing) |
| `USE_LLM_QUALIFIER` | 0 | Enable LLM gate (1=on) |
| `LLM_MODEL` | qwen2.5:3b | Ollama model name |

## 📂 Project Structure

```
abstract-analysis/
├── README.md
├── CLAUDE.md
├── requirements.txt
├── .gitignore
├── Modelfile               # Ollama model configuration
├── hgnc_full.txt           # HGNC gene list (auto-downloaded)
├── models/                 # Model weights (gitignored)
│   └── qwen2.5-3b-instruct-q4_k_m.gguf
├── ollama/                 # Bundled Ollama binaries (tracked)
│   ├── ollama.exe
│   └── lib/ollama/...
├── src/
│   ├── scanner.py          # Core pipeline (~2100 lines)
│   ├── scanner_debug.py    # Fast smoke test (20 abstracts)
│   ├── llm_qualify.py      # LLM client with SQLite cache
│   ├── eval_qualify.py     # Accuracy evaluation vs hand labels
│   ├── generate_report.py  # Post-processing
│   ├── generate_gene_table.py
│   └── reformat.py
├── setup/
│   └── install_ollama.py   # Downloads Ollama from GitHub
├── tools/                  # Diagnostic and benchmark scripts
└── output/                 # Generated files (gitignored)
```

## 🧪 Testing

Quick smoke test (20 abstracts, single disorder):

```bash
python src/scanner_debug.py
```

## 🤝 Credits

![Claude](https://img.shields.io/badge/Claude_Code-claude.com/claude--code-D97757?logo=claude&logoColor=white)
![Anthropic](https://img.shields.io/badge/Anthropic-anthropic.com-19C37D?logo=anthropic&logoColor=white)

Built with 🧡 using Claude Code by Anthropic

---

🦙 LLM inference via [Ollama](https://ollama.com/)

🔀 [OpenRouter](https://openrouter.ai/) ready for cloud inference

---

## 📝 Note for Future Improvements

> Read KNOWN_ISSUES.txt

