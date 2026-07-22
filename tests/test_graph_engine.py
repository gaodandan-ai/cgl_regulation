import os
import sys
import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "backend"))
from app import app
from graph_engine import get_graph_engine

client = TestClient(app)

def test_graph_engine_instance():
    ge = get_graph_engine()
    assert ge is not None
    assert ge.graph.number_of_nodes() > 0

def test_graph_motifs_ffl():
    ge = get_graph_engine()
    motifs = ge.detect_motifs("ffl", limit=10)
    assert isinstance(motifs, list)
    assert len(motifs) > 0

def test_api_graph_cascade_endpoint():
    resp = client.get("/api/graph/cascade?source=sigH&target=cg0350")
    assert resp.status_code == 200
    data = resp.json()
    assert "source" in data
    assert "paths" in data

def test_api_graph_motifs_endpoint():
    resp = client.get("/api/graph/motifs?type=ffl&limit=5")
    assert resp.status_code == 200
    data = resp.json()
    assert data["motif_type"] == "ffl"
    assert data["count"] > 0
    assert len(data["items"]) > 0
