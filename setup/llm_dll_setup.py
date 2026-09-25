"""Register + preload the CUDA runtime DLLs so llama_cpp imports on Windows.

The abetlen cu124 CUDA wheel ships ggml-cuda.dll but NOT the cuBLAS/cuDART/NVRTC
runtime it links against. We pull those via the nvidia-*-cu12 pip wheels. But
llama-cpp-python's own `ctypes.CDLL(lib_path, winmode=...)` uses a restricted
search mode that ignores os.add_dll_directory — so we PRELOAD the CUDA runtime
DLLs into the process with the default winmode first. Once a DLL is loaded into
the process, it's found regardless of the search mode used by later loads.

Import this module (or call setup()) before importing llama_cpp.
"""
import os
import site
import ctypes

_DIRS = []
for _sp in site.getsitepackages():
    _nv = os.path.join(_sp, "nvidia")
    if os.path.isdir(_nv):
        for _sub in os.listdir(_nv):
            _b = os.path.join(_nv, _sub, "bin")
            if os.path.isdir(_b):
                os.add_dll_directory(_b)
                _DIRS.append(_b)
    _lib = os.path.join(_sp, "llama_cpp", "lib")
    if os.path.isdir(_lib):
        os.add_dll_directory(_lib)
        _DIRS.append(_lib)


def _find(name):
    for d in _DIRS:
        p = os.path.join(d, name)
        if os.path.isfile(p):
            return p
    return None


def preload_cuda():
    """Load the CUDA runtime DLLs in dependency order before llama_cpp imports."""
    for name in ("cudart64_12.dll", "nvrtc-builtins64_129.dll",
                 "cublasLt64_12.dll", "cublas64_12.dll", "nvrtc64_120_0.dll"):
        p = _find(name)
        if p:
            try:
                ctypes.CDLL(p)
            except OSError as e:
                print(f"[warn] preload {name} failed: {e}")


preload_cuda()
