import os
import sys

# Set writable directories for Vercel / AWS Lambda serverless environments
if os.environ.get("VERCEL") == "1" or "AWS_LAMBDA_FUNCTION_NAME" in os.environ:
    os.environ["HOME"] = "/tmp"
    os.environ["XDG_CACHE_HOME"] = "/tmp/.cache"

import math
import traceback
from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse


if not hasattr(math, 'comb'):
    def math_comb(n, k):
        if k < 0 or k > n:
            return 0
        if k == 0 or k == n:
            return 1
        k = min(k, n - k)
        numerator = 1
        denominator = 1
        for i in range(1, k + 1):
            numerator *= n - i + 1
            denominator *= i
        return numerator // denominator
    math.comb = math_comb

# Add the parent directory to Python path so we can import backend.app
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(BASE_DIR)
if PARENT_DIR not in sys.path:
    sys.path.insert(0, PARENT_DIR)

app = FastAPI()

@app.get("/api/debug")
def debug_endpoint():
    import platform
    info = {
        "sys.path": sys.path,
        "os.getcwd": os.getcwd(),
        "glibc_version": platform.libc_ver(),
        "environ": {k: v for k, v in os.environ.items() if "KEY" not in k and "PASSWORD" not in k and "SECRET" not in k},
    }
    
    # Search for libexpat in system library directories
    expat_files = []
    for lib_dir in ["/lib64", "/usr/lib64", "/lib", "/usr/lib"]:
        if os.path.exists(lib_dir):
            try:
                for f in os.listdir(lib_dir):
                    if "libexpat" in f:
                        expat_files.append(os.path.join(lib_dir, f))
            except Exception as e:
                expat_files.append(f"error reading {lib_dir}: {e}")
    info["libexpat_files"] = expat_files

    try:
        import cobra
        info["cobra_import"] = "SUCCESS"
    except Exception as e:
        info["cobra_import"] = f"FAILED: {str(e)}"
        info["cobra_import_trace"] = traceback.format_exc()
        
    try:
        from backend.app import app as real_app
        info["backend_app_import"] = "SUCCESS"
    except Exception as e:
        info["backend_app_import"] = f"FAILED: {str(e)}"
        info["backend_app_import_trace"] = traceback.format_exc()
        
    return info

@app.api_route("/{path:path}", methods=["GET", "POST", "PUT", "DELETE", "OPTIONS", "HEAD", "PATCH"])
async def catch_all(request: Request, path: str):
    try:
        from backend.app import app as real_app
        return await real_app(request.scope, request._receive, request._send)
    except Exception as e:
        tb = traceback.format_exc()
        return JSONResponse(
            status_code=500,
            content={
                "error": "Failed to load real backend app",
                "detail": str(e),
                "traceback": tb
            }
        )

