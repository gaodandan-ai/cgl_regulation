import os
import sys
import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "backend"))
from app import app
from db_manager import get_db_manager

client = TestClient(app)


def test_chipseq_tables_exist_and_populated():
    db = get_db_manager()
    conn = db.get_connection()
    assert conn is not None, "Database connection failed"

    cursor = conn.cursor()
    cursor.execute("SELECT COUNT(*) FROM chipseq_peaks;")
    peak_count = cursor.fetchone()[0]
    assert peak_count > 0, "chipseq_peaks table should contain peak records"

    cursor.execute("SELECT COUNT(*) FROM chipseq_regulations;")
    reg_count = cursor.fetchone()[0]
    assert reg_count > 0, "chipseq_regulations table should contain regulation records"


def test_spatial_confidence_tiers():
    db = get_db_manager()
    conn = db.get_connection()
    cursor = conn.cursor()
    cursor.execute("SELECT DISTINCT spatial_confidence FROM chipseq_peaks;")
    tiers = {r[0] for r in cursor.fetchall()}
    assert "PROMOTER_DIRECT" in tiers, "Should classify PROMOTER_DIRECT peaks"
    assert "GENE_BODY_INTERNAL" in tiers or "INTERGENIC_PROMOTER" in tiers


def test_genomic_tracks_includes_chipseq_peaks():
    db = get_db_manager()
    data = db.get_genomic_track_data("cg0350", window_bp=20000)
    assert data is not None
    assert "peaks" in data
    peaks = data["peaks"]
    assert len(peaks) >= 0


def test_chipseq_peaks_api_endpoint():
    resp = client.get("/api/chipseq_peaks/cg0350")
    assert resp.status_code == 200
    payload = resp.json()
    assert payload["query"] == "cg0350"
    assert "as_target_count" in payload
    assert "as_tf_count" in payload


def test_cross_validation_functional_categories():
    db = get_db_manager()
    conn = db.get_connection()
    cursor = conn.cursor()
    cursor.execute("SELECT DISTINCT functional_category FROM regulations WHERE functional_category IS NOT NULL;")
    cats = {r[0] for r in cursor.fetchall()}
    assert "DIRECT_FUNCTIONAL_REGULATION" in cats or "COMPUTATIONAL_PREDICTION_UNVERIFIED" in cats

    cursor.execute("SELECT COUNT(*) FROM v_chipseq_functional_regulations;")
    view_count = cursor.fetchone()[0]
    assert view_count > 0
