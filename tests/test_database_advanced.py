import os
import sys
import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "backend"))
from app import app
from db_manager import get_db_manager

client = TestClient(app)

def test_full_gene_profile_view():
    db = get_db_manager()
    profile = db.get_full_gene_profile("cg0350") # glxR
    assert profile is not None
    assert profile["cg_locus"] == "cg0350"
    assert profile["gene_name"] == "glxR"
    assert "CRP/FNR" in profile["tf_family"]
    assert "cAMP" in profile["effector_molecule"]
    assert "Global Regulator" in str(profile["abasy_role"])

def test_genomic_neighborhood_spatial_query():
    db = get_db_manager()
    genes = db.get_genomic_neighborhood("cg0001", window_bp=10000)
    assert isinstance(genes, list)
    assert len(genes) > 0
    # Center gene cg0001 should be present
    loci = [g["locus_tag"] for g in genes]
    assert "cg0001" in loci

def test_allosteric_feedback_loops_view():
    db = get_db_manager()
    loops = db.get_allosteric_feedback_loops("cAMP")
    assert isinstance(loops, list)
    assert len(loops) > 0
    assert loops[0]["tf_locus"] == "cg0350"

def test_srna_competition_view():
    db = get_db_manager()
    targets = db.get_srna_target_competition("scgl257.1")
    assert isinstance(targets, list)
    assert len(targets) > 0

def test_api_full_gene_profile():
    resp = client.get("/api/gene/profile/cg0350")
    assert resp.status_code == 200
    data = resp.json()
    assert data["cg_locus"] == "cg0350"

def test_api_genomic_neighborhood():
    resp = client.get("/api/gene/neighborhood/cg0001?window_bp=15000")
    assert resp.status_code == 200
    data = resp.json()
    assert data["count"] > 0
    assert len(data["genes"]) > 0

def test_api_allosteric_feedback():
    resp = client.get("/api/network/allosteric-feedback?query=cAMP")
    assert resp.status_code == 200
    data = resp.json()
    assert data["count"] > 0

def test_api_srna_competition():
    resp = client.get("/api/network/srna-competition?srna_id=scgl257.1")
    assert resp.status_code == 200
    data = resp.json()
    assert data["count"] > 0
