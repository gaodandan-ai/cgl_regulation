import os, sys, pytest
from fastapi.testclient import TestClient

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "backend"))
from app import app
from db_manager import get_db_manager

client = TestClient(app)

def test_genomic_track_db_query():
    db = get_db_manager()
    data = db.get_genomic_track_data("cg0350", window_bp=10000)
    assert data is not None
    assert data["query_locus"] == "cg0350"
    assert "target" in data
    assert "genes" in data
    assert len(data["genes"]) > 0

def test_genomic_track_api_endpoint():
    resp = client.get("/api/genomic_tracks/cg0350?window_bp=10000")
    assert resp.status_code == 200
    payload = resp.json()
    assert payload["query_locus"] == "cg0350"
    assert len(payload["genes"]) > 0
    assert "window" in payload
