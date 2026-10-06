"""One-click setup for LLM qualification gate.

Downloads everything needed for a fresh clone:
1. Ollama Windows build (with CUDA DLLs)
2. Qwen2.5-3B-Instruct Q4_K_M model
3. Creates model in Ollama from Modelfile

Run once on a new machine:
    python setup/install_ollama.py

Then start the server:
    ollama/ollama.exe serve
"""
import json
import os
import shutil
import ssl
import subprocess
import sys
import urllib.request
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(HERE)
REPO = "ollama/ollama"
ZIP_PATH = os.path.join(HERE, "ollama-windows.zip")
DEST_DIR = os.path.join(PROJECT_ROOT, "ollama")
MODELS_DIR = os.path.join(PROJECT_ROOT, "models")
MODEL_NAME = "qwen2.5-3b-instruct-q4_k_m"
GGUF_URL = "https://huggingface.co/Qwen/Qwen2.5-3B-Instruct-GGUF/resolve/main/qwen2.5-3b-instruct-q4_k_m.gguf"


def _ctx():
    try:
        return ssl.create_default_context()
    except Exception:
        return ssl.create_unverified_context()


def get_latest():
    url = f"https://api.github.com/repos/{REPO}/releases/latest"
    req = urllib.request.Request(url, headers={"User-Agent": "claude-code"})
    with urllib.request.urlopen(req, context=_ctx(), timeout=30) as r:
        return json.load(r)


def ollama_exe():
    return os.path.join(DEST_DIR, "ollama.exe")


def download_file(url, dest_path, desc="file"):
    """Download a file with progress."""
    print(f"downloading {desc}...", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "claude-code"})
    with urllib.request.urlopen(req, context=_ctx(), timeout=300) as r, open(dest_path, "wb") as f:
        shutil.copyfileobj(r, f)
    print(f"downloaded -> {dest_path}  ({os.path.getsize(dest_path):,} bytes)", flush=True)


def main():
    print("=" * 60, flush=True)
    print("LLM Setup for Abstract Analysis Pipeline", flush=True)
    print("=" * 60, flush=True)

    # 1. Download Ollama if not present
    if os.path.isfile(ollama_exe()):
        print("[OK] Ollama already installed", flush=True)
    else:
        print("\n[1/3] Downloading Ollama...", flush=True)
        rel = get_latest()
        print("latest tag:", rel.get("tag_name"), flush=True)
        asset = None
        for a in rel.get("assets", []):
            if a["name"].lower() == "ollama-windows-amd64.zip":
                asset = a
                break
        if not asset:
            print("ERROR: ollama-windows-amd64.zip not found in release", file=sys.stderr)
            sys.exit(2)
        download_file(asset["browser_download_url"], ZIP_PATH, f"Ollama ({asset['size']/1e6:.0f} MB)")
        print("extracting...", flush=True)
        os.makedirs(DEST_DIR, exist_ok=True)
        with zipfile.ZipFile(ZIP_PATH) as z:
            z.extractall(DEST_DIR)
        print(f"extracted to {DEST_DIR}", flush=True)

    # 2. Download model GGUF if not present
    os.makedirs(MODELS_DIR, exist_ok=True)
    gguf_path = os.path.join(MODELS_DIR, f"{MODEL_NAME}.gguf")

    if os.path.isfile(gguf_path):
        print(f"[OK] Model already downloaded: {gguf_path}", flush=True)
    else:
        print("\n[2/3] Downloading Qwen2.5-3B-Instruct Q4_K_M (2 GB)...", flush=True)
        download_file(GGUF_URL, gguf_path, "Qwen2.5-3B-Instruct Q4_K_M (~2 GB)")

    # 3. Create model in Ollama
    print("\n[3/3] Creating model in Ollama...", flush=True)
    modelfile = os.path.join(PROJECT_ROOT, "Modelfile")
    result = subprocess.run(
        [ollama_exe(), "create", MODEL_NAME, "-f", modelfile],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"[WARN] Failed to create model: {result.stderr}", flush=True)
        print("You may need to run: ollama\\ollama.exe create qwen2.5-3b-instruct-q4_k_m -f Modelfile", flush=True)
    else:
        print(f"[OK] Model {MODEL_NAME} created", flush=True)

    print("\n" + "=" * 60, flush=True)
    print("Setup complete!", flush=True)
    print("=" * 60, flush=True)
    print("\nTo start the server:", flush=True)
    print("    ollama\\ollama.exe serve", flush=True)
    print("\nThen in a new terminal:", flush=True)
    print("    set USE_LLM_QUALIFIER=1", flush=True)
    print("    set LLM_BASE_URL=http://127.0.0.1:11434", flush=True)
    print(f"    set LLM_MODEL={MODEL_NAME}", flush=True)
    print("    python src/scanner.py", flush=True)


if __name__ == "__main__":
    main()
