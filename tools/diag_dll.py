"""Diagnose the llama_cpp DLL load failure: register all bin/ dirs, then try
loading each DLL in dependency order, and check whether the VC++ runtime that
ggml-cuda.dll needs is actually on the system."""
import os
import site
import ctypes

dirs = []
for sp in site.getsitepackages():
    nv = os.path.join(sp, "nvidia")
    if os.path.isdir(nv):
        for sub in os.listdir(nv):
            b = os.path.join(nv, sub, "bin")
            if os.path.isdir(b):
                os.add_dll_directory(b)
                dirs.append(b)
    lib = os.path.join(sp, "llama_cpp", "lib")
    if os.path.isdir(lib):
        os.add_dll_directory(lib)
        dirs.append(lib)
print("registered dirs:", dirs)

print("\nVC++ runtime presence:")
for name in ("MSVCP140.dll", "VCRUNTIME140.dll", "VCRUNTIME140_1.dll"):
    found = []
    for d in (r"C:\Windows\System32", r"C:\Windows\SysWOW64", *site.getsitepackages()):
        p = os.path.join(d, name)
        if os.path.isfile(p):
            found.append(p)
    print(f"  {name:22s} -> {found or 'NOT FOUND'}")


def find(name):
    for d in dirs:
        p = os.path.join(d, name)
        if os.path.isfile(p):
            return p
    return name


print("\nLoad each DLL in dependency order:")
for name in ["cudart64_12.dll", "cublasLt64_12.dll", "cublas64_12.dll",
            "nvrtc64_120_0.dll", "ggml-base.dll", "ggml.dll",
            "ggml-cuda.dll", "llama.dll"]:
    try:
        ctypes.CDLL(find(name))
        print(f"  OK    {name}")
    except OSError as e:
        print(f"  FAIL  {name}  ->  {str(e).splitlines()[-1]}")
