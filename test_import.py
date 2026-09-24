"""Verify llama_cpp imports after preloading the CUDA runtime DLLs."""
import llm_dll_setup  # noqa: F401  (registers + preloads CUDA DLLs)
import llama_cpp
print("IMPORT_OK", getattr(llama_cpp, "__version__", "?"))
