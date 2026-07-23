import os
import sys
import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "backend"))
from app import app
from services.graph_rag_service import get_graph_rag_service

client = TestClient(app)


def test_graph_rag_service_instance():
    svc = get_graph_rag_service()
    assert svc is not None


def test_graph_rag_causal_paths_query():
    svc = get_graph_rag_service()
    paths = svc.get_causal_paths("cg0350", "cg0351", max_depth=3)
    assert isinstance(paths, list)


def test_graph_rag_reasoning_summary():
    svc = get_graph_rag_service()
    res = svc.query_graph_rag_reasoning(source="cg0350", target="cg0351", max_depth=2)
    assert isinstance(res, dict)
    assert "source" in res
    assert "graph_rag_context" in res
    assert "causal_paths" in res


def test_api_graph_rag_reasoning_endpoint():
    resp = client.get("/api/ai/graph_rag/reasoning?source=cg0350&target=cg0351")
    assert resp.status_code == 200
    data = resp.json()
    assert data["source"] == "cg0350"
    assert data["target"] == "cg0351"
    assert "graph_rag_context" in data
