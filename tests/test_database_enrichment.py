import os
import sys
import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "backend"))
from app import app
from db_manager import get_db_manager

client = TestClient(app)

def test_gene_coordinates_db():
    db = get_db_manager()
    coords = db.get_gene_coordinates("cg0001")
    assert coords is not None
    assert coords["locus_tag"] == "cg0001"
    assert "start_pos" in coords
    assert "strand" in coords

def test_tf_effector_info_db():
    db = get_db_manager()
    info = db.get_tf_effector_info("cg0350") # glxR
    assert info is not None
    assert info["tf_name"] == "glxR"
    assert "CRP/FNR" in info["tf_family"]
    assert "cAMP" in info["effector_molecule"]

def test_extended_edges_db():
    db = get_db_manager()
    edges = db.get_extended_edges("cg0350", mode="all")
    assert isinstance(edges, list)
    assert len(edges) > 0

    srna_edges = db.get_extended_edges("scgl257.1", mode="all", edge_type="srna_mrna")
    assert isinstance(srna_edges, list)
    assert len(srna_edges) > 0

def test_gene_coordinates_api_endpoint():
    resp = client.get("/api/gene/coordinates/cg0001")
    assert resp.status_code == 200
    data = resp.json()
    assert data["locus_tag"] == "cg0001"

def test_tf_effectors_api_endpoint():
    resp = client.get("/api/tf/effectors/glxR")
    assert resp.status_code == 200
    data = resp.json()
    assert "CRP/FNR" in data["tf_family"]

def test_extended_network_api_endpoint():
    resp = client.get("/api/network/extended?locus=cg0350&mode=strong")
    assert resp.status_code == 200
    data = resp.json()
    assert data["count"] > 0
