import os
import sys
import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "backend"))
from app import app
from db_manager import get_db_manager

client = TestClient(app)

def test_imodulon_db_lookup():
    db = get_db_manager()
    imodulons = db.get_imodulons_for_gene("cg0369")
    assert isinstance(imodulons, list)
    assert len(imodulons) > 0
    assert "imodulon_id" in imodulons[0]

def test_rf_edge_scores_db_lookup():
    db = get_db_manager()
    scores = db.get_rf_edge_scores("cg0012", min_confidence=0.3)
    assert isinstance(scores, list)
    assert len(scores) > 0
    assert "predicted_confidence" in scores[0]

def test_tf_hierarchy_rankings_db_lookup():
    db = get_db_manager()
    rankings = db.get_tf_hierarchy_rankings()
    assert isinstance(rankings, list)
    assert len(rankings) > 0
    assert "tier" in rankings[0]

def test_rewired_edges_db_lookup():
    db = get_db_manager()
    edges = db.get_rewired_edges("ssuR")
    assert isinstance(edges, list)
    assert len(edges) > 0

def test_api_imodulon_gene_endpoint():
    resp = client.get("/api/imodulon/gene/cg0369")
    assert resp.status_code == 200
    data = resp.json()
    assert data["count"] > 0

def test_api_rf_scores_endpoint():
    resp = client.get("/api/network/rf-scores?locus=cg0012&min_confidence=0.5")
    assert resp.status_code == 200
    data = resp.json()
    assert data["count"] > 0

def test_api_tf_hierarchy_rankings_endpoint():
    resp = client.get("/api/tf/hierarchy-rankings")
    assert resp.status_code == 200
    data = resp.json()
    assert data["count"] > 0

def test_api_network_rewired_endpoint():
    resp = client.get("/api/network/rewired?locus=ssuR")
    assert resp.status_code == 200
    data = resp.json()
    assert data["count"] > 0
