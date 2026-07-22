import os
import sys
import time
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "backend"))
from db_manager import get_db_manager

def test_database_connection():
    db = get_db_manager()
    conn = db.get_connection()
    assert conn is not None

def test_canonical_locus_map():
    db = get_db_manager()
    # Test CG locus
    res_cg = db.get_canonical_locus("cg0001")
    assert res_cg is not None
    assert res_cg["canonical_cg"] == "cg0001"

    # Test Cgl locus
    res_cgl = db.get_canonical_locus("Cgl0419")
    assert res_cgl is not None
    assert res_cgl["canonical_cg"] == "cg0499"

def test_essential_genes():
    db = get_db_manager()
    res = db.get_essential_gene("cg0001")
    assert res is not None
    assert "essentiality" in res

def test_abasy_roles():
    db = get_db_manager()
    role = db.get_abasy_role("cg0350") # glxR
    assert role is not None
    assert "Global Regulator" in role

def test_string_interactions():
    db = get_db_manager()
    interactions = db.get_string_interactions("cg0001", min_score=400)
    assert isinstance(interactions, list)
    assert len(interactions) > 0

def test_fts5_literature_search():
    db = get_db_manager()
    results = db.search_literature_fts("sigH")
    assert isinstance(results, list)
    assert len(results) > 0

def test_query_speed_benchmark():
    db = get_db_manager()
    t0 = time.time()
    for _ in range(100):
        db.get_canonical_locus("cg0001")
        db.get_essential_gene("cg0350")
        db.get_abasy_role("cg0350")
    t1 = time.time()
    elapsed = t1 - t0
    print(f"100 DB queries took: {elapsed*1000:.2f} ms")
    assert elapsed < 1.0 # Must finish 100 queries in under 1 second!
