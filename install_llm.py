"""One-click LLM setup for Windows. Downloads Ollama and pulls the model.

Run once on a new machine:
    python install_llm.py

Then run scanner with LLM gate:
    set USE_LLM_QUALIFIER=1
    python scanner.py
"""
import subprocess
import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OLLAMA_DIR = os.path.join(HERE, "ollama")

# Default model (matches llm_qualify.py default)
DEFAULT_MODEL = os.environ.get("LLM_MODEL", "qwen2.5:3b")


def ollama_exe():
    """Path to ollama.exe if downloaded, else 'ollama' from PATH."""
    exe = os.path.join(OLLAMA_DIR, "ollama.exe")
    return exe if os.path.isfile(exe) else "ollama"


def is_ollama_installed():
    """Check if ollama command is available."""
    try:
        subprocess.run([ollama_exe(), "--version"], capture_output=True, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def download_ollama():
    """Download portable Ollama Windows build if not present."""
    if is_ollama_installed():
        print("[OK] Ollama already installed")
        return

    print("[1/2] Downloading Ollama...")
    try:
        # Try running install_ollama.py if it exists
        import install_ollama as installer
        installer.main()
        print("[OK] Ollama downloaded")
    except ImportError:
        print("[WARN] install_ollama.py not found, please install Ollama manually from:")
        print("       https://github.com/ollama/ollama/releases")
        sys.exit(1)


def pull_model(model=DEFAULT_MODEL):
    """Pull the model from Ollama registry."""
    print(f"[2/2] Pulling model {model} (this may take a few minutes)...")
    try:
        subprocess.run([ollama_exe(), "pull", model], check=True)
        print(f"[OK] Model {model} ready")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to pull model: {e}")
        sys.exit(1)


def main():
    print("=" * 50)
    print("LLM Setup for Abstract Analysis Pipeline")
    print("=" * 50)

    download_ollama()
    pull_model()

    print("\n" + "=" * 50)
    print("Setup complete!")
    print("=" * 50)
    print("\nTo run with LLM qualification gate:")
    print(f"    set USE_LLM_QUALIFIER=1")
    print(f"    python scanner.py")
    print("\nOr with PowerShell:")
    print(f"    $env:USE_LLM_QUALIFIER=1")
    print(f"    python scanner.py")


if __name__ == "__main__":
    main()
