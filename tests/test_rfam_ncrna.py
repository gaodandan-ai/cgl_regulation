import os
import sys
import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "backend"))
from app import app
from db_manager import get_db_manager

client = TestClient(app)

def test_rfam_ncrnas_db_lookup():
    db = get_db_manager()
    ncrnas = db.get_ncrnas()
    assert isinstance(ncrnas, list)
    assert len(ncrnas) > 0
    assert "rfam_acc" in ncrnas[0]
    assert "rna_type" in ncrnas[0]

def test_srna_targets_db_lookup():
    db = get_db_manager()
    targets = db.get_srna_targets("scgl257.1")
    assert isinstance(targets, list)
    assert len(targets) > 0
    assert "binding_energy_kcal" in targets[0]
    assert "target_region_type" in targets[0]

def test_api_ncrna_list_endpoint():
    resp = client.get("/api/ncrna/list?rna_type=Riboswitch")
    assert resp.status_code == 200
    data = resp.json()
    assert data["count"] > 0
    assert len(data["ncrnas"]) > 0

def test_api_srna_targets_endpoint():
    resp = client.get("/api/ncrna/targets?locus=scgl257.1")
    assert resp.status_code == 200
    data = resp.json()
    assert data["count"] > 0
    assert len(data["targets"]) > 0
