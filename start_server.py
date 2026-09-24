"""Launch llama-cpp-python's OpenAI server with CUDA runtime DLLs preloaded.

Run instead of `python -m llama_cpp.server`:
    python start_server.py --model models/qwen2.5-3b-instruct-q4_k_m.gguf \
        --n_ctx 4096 --n_gpu_layers -1 --n_parallel 1 \
        --host 127.0.0.1 --port 8080
"""
import llm_dll_setup  # noqa: F401  (registers + preloads CUDA DLLs before import)
import runpy

runpy.run_module("llama_cpp.server", run_name="__main__", alter_sys=True)
