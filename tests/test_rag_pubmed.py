import os
import sys
import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "backend"))
from app import app
from rag_service import RAGService

client = TestClient(app)


def test_rag_service_fetch_europe_pmc():
    svc = RAGService()
    articles = svc.fetch_europe_pmc_abstracts("cg0350 glxR", max_results=3)
    assert isinstance(articles, list)
    if articles:
        art = articles[0]
        assert "title" in art
        assert "abstract" in art
        assert "pmid" in art


def test_rag_service_ingest_articles(tmp_path):
    svc = RAGService()
    dummy_articles = [{
        "pmid": "99999999",
        "title": "Test Corynebacterium glutamicum Regulation Article",
        "authors": "Test Author et al.",
        "journal": "Test Journal",
        "year": "2026",
        "doi": "10.1000/test.doi",
        "abstract": "GlxR regulates multiple genes in Corynebacterium glutamicum."
    }]
    count = svc.ingest_fetched_articles(dummy_articles, tag="test")
    assert count == 1


def test_api_fetch_pubmed_endpoint():
    resp = client.post("/api/ai/literature/fetch_pubmed", json={
        "query": "cg0350",
        "max_results": 2,
        "auto_ingest": False
    })
    assert resp.status_code == 200
    data = resp.json()
    assert "query" in data
    assert "articles" in data
    assert isinstance(data["articles"], list)
