from fastapi import FastAPI, HTTPException, Header, Response, Request
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from contextlib import asynccontextmanager
import logging
import sys
import os
import json
import mimetypes
from fastapi.responses import JSONResponse, FileResponse as _FileResp, Response as _Resp
from starlette.exceptions import HTTPException as StarletteHTTPException

if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    BACKEND_DIR = os.path.join(sys._MEIPASS, "backend")
    PARENT_DIR = sys._MEIPASS
else:
    BACKEND_DIR = os.path.dirname(os.path.abspath(__file__))
    PARENT_DIR = os.path.dirname(BACKEND_DIR)

if BACKEND_DIR not in sys.path:
    sys.path.insert(0, BACKEND_DIR)
if PARENT_DIR not in sys.path:
    sys.path.insert(0, PARENT_DIR)

try:
    from security import PUBLIC_DEPLOYMENT
except ImportError:
    from backend.security import PUBLIC_DEPLOYMENT

import run_server
from model_loader import load_model_if_needed
from services.reference_data import (
    load_all_reference_data,
    check_essentiality,
    check_abasy_role,
    ESSENTIAL_GENES,
    PRODORIC_PWMS,
    BRENDA_KCAT_MAPPINGS,
    STRING_INTERACTIONS,
    ABASY_ROLES,
    RHEA_MAPPINGS,
    CHEBI_MAPPINGS,
    COG_ANNOTATIONS
)

try:
    from db_manager import get_db_manager
    _DB_MANAGER_AVAILABLE = True
except ImportError:
    get_db_manager = None
    _DB_MANAGER_AVAILABLE = False

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("app")


@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("Initializing FBA simulator service...")
    if PUBLIC_DEPLOYMENT:
        logger.info("Deferring metabolic model load in the public serverless deployment.")
    else:
        try:
            load_model_if_needed()
        except Exception as e:
            logger.warning(f"Initial model load failed (will retry on demand): {str(e)}")
    
    # Initialize run_server mappings and caches
    try:
        run_server.load_gene_mappings()
        run_server.load_organism_kegg_links()
        logger.info("Successfully loaded gene mappings and KEGG links from run_server.")
    except Exception as e:
        logger.warning(f"Failed to load run_server mappings/caches: {str(e)}")

    # Load reference JSON databases
    try:
        load_all_reference_data(PARENT_DIR)
    except Exception as e:
        logger.error(f"Error loading reference datasets: {e}")

    yield
    if _DB_MANAGER_AVAILABLE and get_db_manager():
        get_db_manager().close()


def _application_version() -> str:
    try:
        with open(os.path.join(PARENT_DIR, "web", "version.json"), encoding="utf-8") as stream:
            return str(json.load(stream).get("version", "0.0.0"))
    except (OSError, ValueError, TypeError):
        return "0.0.0"


APP_VERSION = _application_version()
app = FastAPI(
    title="Cgl Regulation FBA Simulator API",
    version=APP_VERSION,
    lifespan=lifespan,
    docs_url=None if PUBLIC_DEPLOYMENT else "/docs",
    redoc_url=None if PUBLIC_DEPLOYMENT else "/redoc",
    openapi_url=None if PUBLIC_DEPLOYMENT else "/openapi.json",
)

_configured_origins = [
    origin.strip()
    for origin in os.environ.get(
        "CGL_ALLOWED_ORIGINS",
        "https://cgl-regulation.vercel.app,http://127.0.0.1:8000,http://localhost:8000",
    ).split(",")
    if origin.strip()
]
app.add_middleware(
    CORSMiddleware,
    allow_origins=_configured_origins,
    allow_credentials=False,
    allow_methods=["GET", "POST", "OPTIONS"],
    allow_headers=[
        "Accept", "Content-Type", "If-None-Match",
        "X-AI-API-Key", "X-Gemini-API-Key", "X-AI-Provider",
        "X-AI-Model", "X-AI-Base-URL",
    ],
)

# Enable Gzip compression to speed up transfer of large datasets
app.add_middleware(GZipMiddleware, minimum_size=1000)


@app.middleware("http")
async def security_headers(request: Request, call_next):
    response = await call_next(request)
    response.headers["X-Content-Type-Options"] = "nosniff"
    response.headers["X-Frame-Options"] = "DENY"
    response.headers["Referrer-Policy"] = "strict-origin-when-cross-origin"
    response.headers["Permissions-Policy"] = "camera=(), microphone=(), geolocation=()"
    response.headers["Content-Security-Policy"] = (
        "default-src 'self'; object-src 'none'; base-uri 'self'; frame-ancestors 'none'; "
        "script-src 'self' 'unsafe-inline' https://cdn.jsdelivr.net; "
        "style-src 'self' 'unsafe-inline' https://fonts.googleapis.com https://cdnjs.cloudflare.com; "
        "font-src 'self' data: https://fonts.gstatic.com https://cdnjs.cloudflare.com; "
        "img-src 'self' data: blob: https:; connect-src 'self' https: http://127.0.0.1:* http://localhost:*; "
        "worker-src 'self' blob:; form-action 'self'"
    )
    return response


@app.exception_handler(StarletteHTTPException)
async def sanitized_http_exception(request: Request, exc: StarletteHTTPException):
    if exc.status_code >= 500:
        logger.error("Server error on %s: %s", request.url.path, exc.detail)
        return JSONResponse(status_code=exc.status_code, content={"detail": "Internal server error"})
    return JSONResponse(status_code=exc.status_code, content={"detail": exc.detail}, headers=exc.headers)


@app.exception_handler(Exception)
async def sanitized_unhandled_exception(request: Request, exc: Exception):
    logger.exception("Unhandled error on %s", request.url.path)
    return JSONResponse(status_code=500, content={"detail": "Internal server error"})


# ── Register Routers ─────────────────────────────────────────────────────────
try:
    from routers.system import router as system_router
    from routers.gene import router as gene_router
    from routers.simulation import router as simulation_router
    from routers.ai import router as ai_router
    from routers.network import router as network_router
except ImportError:
    from backend.routers.system import router as system_router
    from backend.routers.gene import router as gene_router
    from backend.routers.simulation import router as simulation_router
    from backend.routers.ai import router as ai_router
    from backend.routers.network import router as network_router

app.include_router(system_router)
app.include_router(gene_router)
app.include_router(simulation_router)
app.include_router(ai_router)
app.include_router(network_router)


# ── Static file serving with Cache-Control headers ────────────────────────────
web_dir = os.path.join(PARENT_DIR, "web")
data_dir = os.path.join(PARENT_DIR, "data", "reference")


def _etag(path: str) -> str:
    stat = os.stat(path)
    return f'"{stat.st_size}-{int(stat.st_mtime)}"'


def _file_response(path: str, cache_seconds: int, request: Request = None) -> _Resp:
    if not os.path.isfile(path):
        raise HTTPException(status_code=404)
    etag = _etag(path)
    if request:
        client_etag = request.headers.get("If-None-Match", "")
        if client_etag == etag:
            return _Resp(status_code=304)
    mt, _ = mimetypes.guess_type(path)
    resp = _FileResp(path, media_type=mt or "application/octet-stream")
    resp.headers["Cache-Control"] = f"public, max-age={cache_seconds}"
    resp.headers["ETag"] = etag
    resp.headers["Vary"] = "Accept-Encoding"
    return resp


@app.get("/data/{path:path}")
async def serve_data(path: str, request: Request):
    """Serve data reference files with 1-hour browser caching."""
    full = os.path.realpath(os.path.join(data_dir, path))
    base = os.path.realpath(data_dir)
    if os.path.commonpath((base, full)) != base:
        raise HTTPException(status_code=403)
    return _file_response(full, cache_seconds=3600, request=request)


@app.get("/{path:path}")
async def serve_static(path: str, request: Request):
    """Serve web static files. JS/CSS/images get 7-day cache, HTML gets no-cache."""
    if not path or path == "/":
        path = "index.html"
    if path == "api" or path.startswith("api/"):
        raise HTTPException(status_code=404, detail="API endpoint not found")
    full = os.path.realpath(os.path.join(web_dir, path))
    base = os.path.realpath(web_dir)
    if os.path.commonpath((base, full)) != base:
        raise HTTPException(status_code=403)
    if not os.path.isfile(full):
        full = os.path.join(web_dir, "index.html")
    if full.endswith(".html"):
        return _file_response(full, cache_seconds=0, request=request)
    if os.path.basename(full) in ("version.json", "version_local.json"):
        return _file_response(full, cache_seconds=0, request=request)
    if full.endswith((".js", ".css", ".ico", ".png", ".svg", ".woff", ".woff2")):
        return _file_response(full, cache_seconds=604800, request=request)
    return _file_response(full, cache_seconds=3600, request=request)
