import os
import sys
import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "backend"))
from app import app
from db_manager import get_db_manager

client = TestClient(app)

def test_collectf_db_lookup():
    db = get_db_manager()
    sites = db.get_collectf_tfbs("GlxR")
    assert isinstance(sites, list)
    assert len(sites) > 0
    assert "sequence" in sites[0]
    assert "source_tag" in sites[0]

def test_collectf_deduplication_tags():
    db = get_db_manager()
    sites = db.get_collectf_tfbs()
    source_tags = [s["source_tag"] for s in sites]
    assert any("RegPrecise + CollecTF" in tag for tag in source_tags)
    assert any("CollecTF (Novel Validated)" in tag for tag in source_tags)

def test_api_collectf_tfbs_endpoint():
    resp = client.get("/api/tfbs/collectf?locus=cg0350")
    assert resp.status_code == 200
    data = resp.json()
    assert data["count"] > 0
    assert len(data["sites"]) > 0
