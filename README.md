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
pip install biopython pandas matplotlib matplotlib-venn openpyxl httpx
```

### 🔑 Set NCBI API Key (Recommended)

Set your NCBI API key to increase rate limits:

```bash
set NCBI_API_KEY=your_key_here
```

### ▶️ Run the Pipeline

```bash
python scanner.py              # Fetch PubMed + filter + classify
python generate_report.py      # Generate summary and Venn diagrams
python generate_gene_table.py # Generate consolidated gene table
```

### 🤖 Optional: LLM Qualification Gate

To use the local LLM qualification gate (hybrid approach):

```bash
python install_llm.py          # One-click setup: downloads Ollama + model
set USE_LLM_QUALIFIER=1        # Enable LLM gate
python scanner.py
```

## 📁 Output Files

| File | Description |
|------|-------------|
| `model_organisms_dataset.csv` | Primary output: one row per abstract×disorder |
| `passing_abstracts.csv` | Abstracts that passed all filters |
| `consolidated_gene_table.csv` | Gene counts by disorder |
| `gene_report.xlsx` | Excel workbook with all results |
| `venn_diagrams.png` | Brain/gut gene overlap per disorder |

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
├── scanner.py           # Core pipeline (~2100 lines)
├── scanner_debug.py     # Fast smoke test (20 abstracts)
├── llm_qualify.py       # LLM client with SQLite cache
├── eval_qualify.py      # Accuracy evaluation vs hand labels
├── install_llm.py       # One-click LLM setup for Windows
├── generate_report.py   # Post-processing
├── generate_gene_table.py
├── Modelfile            # Ollama model configuration
└── tools/               # Diagnostic and benchmark scripts
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

> Possibly good idea to add direct disorder check.
