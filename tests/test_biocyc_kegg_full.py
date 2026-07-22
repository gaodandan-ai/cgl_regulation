import os
import sys
import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "backend"))
from app import app
from db_manager import get_db_manager

client = TestClient(app)

def test_full_pathways_db_lookup():
    db = get_db_manager()
    pathways = db.get_pathways_for_gene("cg0350")
    assert isinstance(pathways, list)
    assert len(pathways) > 0
    assert "pathway_id" in pathways[0]
    assert "pathway_name" in pathways[0]

def test_genes_in_pathway_db_lookup():
    db = get_db_manager()
    genes = db.get_genes_in_pathway("cgl00010")
    assert isinstance(genes, list)
    assert len(genes) > 0

def test_api_pathways_for_gene_endpoint():
    resp = client.get("/api/pathway/gene/cg0350")
    assert resp.status_code == 200
    data = resp.json()
    assert data["count"] > 0
    assert len(data["pathways"]) > 0

def test_api_genes_in_pathway_endpoint():
    resp = client.get("/api/pathway/info/cgl00010")
    assert resp.status_code == 200
    data = resp.json()
    assert data["count"] > 0
    assert len(data["genes"]) > 0
