import os, sys, pytest
from fastapi.testclient import TestClient

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "backend"))
from app import app
from ai_handlers import build_grounded_omics_context, handle_ai_engineering_command

client = TestClient(app)

def test_build_grounded_omics_context():
    context = build_grounded_omics_context("cg0350") # glxR
    assert context is not None
    assert "Hard Empirical Data" in context
    assert "cg0350" in context

def test_handle_ai_engineering_command_dummy():
    res = handle_ai_engineering_command(
        command="crispri",
        gene="cg0350",
        provider="openai",
        api_key="DummyKeyForTest"
    )
    assert res is not None
    assert res["command"] == "crispri"
    assert res["gene"] == "cg0350"
    assert "grounded_facts" in res

def test_ai_engineering_command_api_endpoint():
    resp = client.post("/api/ai/engineering_command", json={
        "command": "bottleneck",
        "gene": "cg0350",
        "provider": "openai",
        "api_key": "DummyKeyForTest"
    })
    assert resp.status_code == 200
    data = resp.json()
    assert data["command"] == "bottleneck"
    assert data["gene"] == "cg0350"
    assert "Hard Empirical Data" in data["grounded_facts"]
