"""Download the portable Ollama Windows build from GitHub releases and extract
ollama.exe locally. Ollama ships llama.cpp with runtime CPU dispatch (picks the
AVX2 path on the i9-9980HK, so no AVX-512 crash) and bundles its own CUDA
runtime, so it works on this machine where the abetlen cu124 wheel did not.
"""
import json
import os
import shutil
import ssl
import sys
import urllib.request
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = "ollama/ollama"
ZIP_PATH = os.path.join(HERE, "ollama-windows.zip")
DEST_DIR = os.path.join(HERE, "ollama")


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


def main():
    rel = get_latest()
    print("latest tag:", rel.get("tag_name"), flush=True)
    asset = None
    for a in rel.get("assets", []):
        if a["name"].lower() == "ollama-windows-amd64.zip":
            asset = a
            break
    if not asset:
        print("windows assets:", flush=True)
        for a in rel.get("assets", []):
            if "windows" in a["name"].lower():
                print("  ", a["name"], a["browser_download_url"], flush=True)
        print("ERROR: ollama-windows-amd64.zip not found", file=sys.stderr)
        sys.exit(2)
    size_mb = asset.get("size", 0) / 1e6
    print(f"downloading {asset['browser_download_url']}  ({size_mb:.0f} MB) ...", flush=True)
    req = urllib.request.Request(asset["browser_download_url"], headers={"User-Agent": "claude-code"})
    with urllib.request.urlopen(req, context=_ctx(), timeout=180) as r, open(ZIP_PATH, "wb") as f:
        shutil.copyfileobj(r, f)
    print(f"downloaded -> {ZIP_PATH}  {os.path.getsize(ZIP_PATH)} bytes", flush=True)
    os.makedirs(DEST_DIR, exist_ok=True)
    with zipfile.ZipFile(ZIP_PATH) as z:
        z.extractall(DEST_DIR)
    print(f"extracted to {DEST_DIR}", flush=True)
    exe = None
    for root, _dirs, files in os.walk(DEST_DIR):
        for fn in files:
            p = os.path.join(root, fn)
            print("  ", p, flush=True)
            if fn.lower() == "ollama.exe":
                exe = p
    print("OLLAMA_EXE=" + (exe or "NOT FOUND"), flush=True)


if __name__ == "__main__":
    main()
