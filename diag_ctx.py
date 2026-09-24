"""Bisect the STATUS_ILLEGAL_INSTRUCTION (0xc000001d) crash.

The server got through model load (all layers on GPU) and crashed in
llama_init_from_model (context/graph creation). 0xc000001d = illegal CPU
instruction. Test: does a tiny CPU-only context init survive? If yes, the GPU
offload path is the culprit; if no, the wheel's ggml CPU code needs an
instruction this CPU lacks (-> need a different wheel).
"""
import llm_dll_setup  # noqa: F401
import traceback
import llama_cpp

MP = "models/qwen2.5-3b-instruct-q4_k_m.gguf"
print("llama_cpp", getattr(llama_cpp, "__version__", "?"), flush=True)


def try_init(label, **kw):
    print(f"\n=== {label}: {kw} ===", flush=True)
    try:
        m = llama_cpp.Llama(model_path=MP, verbose=False, **kw)
        print(f"  ctx OK", flush=True)
        out = m.create_completion("The gene is", max_tokens=4, echo=False)
        print(f"  gen OK -> {out['choices'][0]['text']!r}", flush=True)
        del m
        return True
    except Exception as e:
        print(f"  FAIL {label}: {type(e).__name__}: {str(e)[:200]}", flush=True)
        traceback.print_exc()
        return False


try_init("CPU n_ctx=1024", n_ctx=1024, n_gpu_layers=0, n_batch=512)
try_init("GPU n_ctx=1024", n_ctx=1024, n_gpu_layers=-1, n_batch=512)
print("\nDONE", flush=True)
