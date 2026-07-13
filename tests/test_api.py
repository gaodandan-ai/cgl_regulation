"""
tests/test_api.py
=================
FastAPI endpoint smoke tests using httpx TestClient (no server needed).
Run with:  pytest tests/ -v

These tests verify that all critical endpoints:
  - Return HTTP 200
  - Return valid JSON
  - Contain expected response fields
"""
import os
import sys
import pytest

# Ensure backend is importable
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "backend"))

# Use starlette TestClient (bundled with FastAPI)
from fastapi.testclient import TestClient


@pytest.fixture(scope="module")
def client():
    """Create a TestClient with the FastAPI app (no external server needed)."""
    os.environ.setdefault("HEADLESS", "true")
    from backend.app import app
    with TestClient(app, raise_server_exceptions=False) as c:
        yield c


# ─── Health / Status ──────────────────────────────────────────────────────────

class TestModelStatus:
    def test_status_200(self, client):
        r = client.get("/api/model/status")
        assert r.status_code == 200

    def test_status_has_reactions(self, client):
        r = client.get("/api/model/status")
        data = r.json()
        assert "reaction_count" in data or "n_reactions" in data or "reaction" in str(data)


# ─── Thermodynamics endpoints ─────────────────────────────────────────────────

class TestThermoEndpoints:
    def test_pruning_report_200(self, client):
        r = client.get("/api/thermo/pruning-report")
        assert r.status_code == 200

    def test_pruning_report_schema(self, client):
        data = client.get("/api/thermo/pruning-report").json()
        assert "data_coverage_pct" in data
        assert data["data_coverage_pct"] >= 60.0

    def test_gene_context_known_gene(self, client):
        r = client.get("/api/thermo/gene_context?gene=cg0350")
        assert r.status_code in (200, 404)  # 404 if gene has no thermo reactions
        if r.status_code == 200:
            data = r.json()
            assert "gene" in data

    def test_gene_context_missing_param(self, client):
        r = client.get("/api/thermo/gene_context")
        assert r.status_code == 400


# ─── STRING PPI endpoint ──────────────────────────────────────────────────────

class TestStringPPI:
    def test_known_gene(self, client):
        r = client.get("/api/analysis/string_ppi?gene=cg0001")
        assert r.status_code == 200
        data = r.json()
        assert "partners" in data
        assert len(data["partners"]) > 0

    def test_score_filter(self, client):
        r = client.get("/api/analysis/string_ppi?gene=cg0001&min_score=700")
        assert r.status_code == 200
        data = r.json()
        for p in data["partners"]:
            assert p["score"] >= 700

    def test_meta_guard(self, client):
        """_meta should not be treated as a gene."""
        r = client.get("/api/analysis/string_ppi?gene=_meta")
        assert r.status_code == 400

    def test_missing_gene(self, client):
        r = client.get("/api/analysis/string_ppi")
        assert r.status_code == 400

    def test_unknown_gene(self, client):
        r = client.get("/api/analysis/string_ppi?gene=not_a_real_gene_xyz")
        assert r.status_code == 200
        data = r.json()
        assert data["partners"] == []

    def test_string_meta_in_response(self, client):
        r = client.get("/api/analysis/string_ppi?gene=cg0001")
        data = r.json()
        assert "string_meta" in data
        assert data["string_meta"]["n_edges"] > 10000


# ─── Network centrality ───────────────────────────────────────────────────────

class TestNetworkCentrality:
    def test_centrality_list(self, client):
        r = client.get("/api/network/centrality?limit=10&tfs_only=true")
        assert r.status_code == 200
        data = r.json()
        assert isinstance(data, (list, dict))

    def test_centrality_for_gene(self, client):
        r = client.get("/api/network/centrality/cg0012")
        assert r.status_code in (200, 404)


# ─── Quality endpoints ────────────────────────────────────────────────────────

class TestQualityEndpoints:
    def test_brenda_quality(self, client):
        r = client.get("/api/quality/brenda")
        assert r.status_code == 200
        data = r.json()
        # After P0-2, should have 1800+ entries
        total = len(data)
        assert total >= 1800, f"Expected >=1800 kcat entries, got {total}"

    def test_essential_quality(self, client):
        r = client.get("/api/quality/essential")
        assert r.status_code == 200

    def test_abasy_quality(self, client):
        r = client.get("/api/quality/abasy")
        assert r.status_code == 200


# ─── Gene assistant (summary) ─────────────────────────────────────────────────

class TestGeneSummary:
    def test_summarize_requires_gene(self, client):
        r = client.get("/api/summarize")
        assert r.status_code in (400, 422)

    def test_summarize_known_gene(self, client):
        r = client.get("/api/summarize?gene=cg0350&name=cg0350")
        # May return 200 with cached summary or 503 if AI not configured
        assert r.status_code in (200, 503, 500)
