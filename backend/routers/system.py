from fastapi import APIRouter
import os
import sys
import json
import logging
from pathlib import Path

try:
    from db_manager import get_db_manager
    _DB_MANAGER_AVAILABLE = True
except ImportError:
    try:
        from backend.db_manager import get_db_manager
        _DB_MANAGER_AVAILABLE = True
    except ImportError:
        get_db_manager = None
        _DB_MANAGER_AVAILABLE = False

try:
    from security import PUBLIC_DEPLOYMENT
except ImportError:
    try:
        from backend.security import PUBLIC_DEPLOYMENT
    except ImportError:
        PUBLIC_DEPLOYMENT = False

try:
    from services.provenance import build_provenance
except ImportError:
    from backend.services.provenance import build_provenance

router = APIRouter(tags=["System"])
logger = logging.getLogger("app.routers.system")


def get_app_version() -> str:
    try:
        if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
            parent_dir = sys._MEIPASS
        else:
            backend_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            parent_dir = os.path.dirname(backend_dir)
        with open(os.path.join(parent_dir, "web", "version.json"), encoding="utf-8") as stream:
            return str(json.load(stream).get("version", "0.0.0"))
    except (OSError, ValueError, TypeError):
        return "0.0.0"


APP_VERSION = get_app_version()


@router.get("/api/health")
def health():
    db_path = getattr(get_db_manager(), "_db_path", "") if _DB_MANAGER_AVAILABLE else ""
    return {
        "status": "ok",
        "app": "cgl-regulation",
        "version": APP_VERSION,
        "database": "available" if db_path and os.path.isfile(db_path) else "unavailable",
        "deployment": "public" if PUBLIC_DEPLOYMENT else "local",
    }


@router.get("/api/provenance")
def provenance():
    backend_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    root = os.path.dirname(backend_dir)
    manager = get_db_manager() if _DB_MANAGER_AVAILABLE else None
    return build_provenance(
        manager=manager,
        root=Path(root),
        version=APP_VERSION,
        deployment="public" if PUBLIC_DEPLOYMENT else "local",
    )


@router.get("/api/check-update")
def check_update():
    return {
        "version": APP_VERSION,
        "up_to_date": True,
        "message": "Using latest release."
    }
