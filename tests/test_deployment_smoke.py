import json
import os
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient


ROOT = Path(__file__).resolve().parents[1]
BACKEND = ROOT / "backend"
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))


@pytest.fixture(scope="module")
def client():
    from backend.app import app

    with TestClient(app) as test_client:
        yield test_client


def test_health_identifies_application_and_database(client):
    response = client.get("/api/health")
    assert response.status_code == 200
    payload = response.json()
    assert payload["app"] == "cgl-regulation"
    assert payload["status"] == "ok"
    assert payload["database"] == "available"


def test_security_headers_and_unknown_api(client):
    response = client.get("/api/health")
    assert response.headers["x-content-type-options"] == "nosniff"
    assert response.headers["x-frame-options"] == "DENY"
    assert "frame-ancestors 'none'" in response.headers["content-security-policy"]

    missing = client.get("/api/debug")
    assert missing.status_code == 404
    assert "sys.path" not in missing.text
    assert "traceback" not in missing.text.lower()


def test_cors_rejects_untrusted_origin(client):
    response = client.options(
        "/api/health",
        headers={
            "Origin": "https://evil.example",
            "Access-Control-Request-Method": "GET",
        },
    )
    assert response.status_code == 400
    assert "access-control-allow-origin" not in response.headers


def test_public_database_contains_new_workflows():
    database = ROOT / "data" / "deploy" / "cgl_regulation_public.db"
    assert database.is_file()
    assert database.stat().st_size < 75 * 1024 * 1024
    connection = sqlite3.connect(f"file:{database.as_posix()}?mode=ro", uri=True)
    try:
        assert connection.execute("SELECT COUNT(*) FROM condition_analysis_runs").fetchone()[0] > 0
        assert connection.execute("SELECT COUNT(*) FROM v_intervention_target_priorities").fetchone()[0] > 0
    finally:
        connection.close()


def test_serverless_entrypoint_has_no_debug_route():
    env = os.environ.copy()
    env.update({
        "VERCEL": "1",
        "CGL_DATABASE_PATH": str(ROOT / "data" / "deploy" / "cgl_regulation_public.db"),
    })
    code = """
import json
from fastapi.testclient import TestClient
from api.index import app
paths = sorted(route.path for route in app.routes)
with TestClient(app) as client:
    statuses = {
        path: client.get(path).status_code
        for path in (
            "/api/health",
            "/api/condition-regulation/runs",
            "/api/intervention-targets?limit=1",
        )
    }
print(json.dumps({"paths": paths, "statuses": statuses}))
"""
    result = subprocess.run(
        [sys.executable, "-c", code], cwd=ROOT, env=env,
        capture_output=True, text=True, timeout=45, check=True,
    )
    output = json.loads(result.stdout.strip().splitlines()[-1])
    paths = output["paths"]
    assert "/api/debug" not in paths
    assert "/docs" not in paths
    assert "/openapi.json" not in paths
    assert "/api/condition-regulation/runs" in paths
    assert set(output["statuses"].values()) == {200}
