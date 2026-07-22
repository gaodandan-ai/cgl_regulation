import os
import sys
import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "backend"))
from app import app
from db_manager import get_db_manager

client = TestClient(app)

def test_ica_condition_activities_db_lookup():
    db = get_db_manager()
    activities = db.get_condition_specific_regulons("Glucose")
    assert isinstance(activities, list)
    assert len(activities) > 0
    assert "activity_score" in activities[0]
    assert "imodulon_id" in activities[0]

def test_imodulon_regulon_overlaps_db_lookup():
    db = get_db_manager()
    overlaps = db.get_imodulon_regulon_overlap()
    assert isinstance(overlaps, list)
    assert len(overlaps) > 0
    assert "f1_score" in overlaps[0]

def test_api_imodulon_condition_endpoint():
    resp = client.get("/api/imodulon/condition?condition=Glucose")
    assert resp.status_code == 200
    data = resp.json()
    assert data["count"] > 0
    assert len(data["activities"]) > 0

def test_api_imodulon_overlap_endpoint():
    resp = client.get("/api/imodulon/overlap")
    assert resp.status_code == 200
    data = resp.json()
    assert data["count"] > 0
    assert len(data["overlaps"]) > 0
