"""scripts/add_centrality_endpoint.py — injects network centrality endpoints into app.py"""

ENDPOINTS_CODE = r'''
# ── Network Centrality Endpoints ───────────────────────────────────────────────
_CENTRALITY_DATA = None

def _load_centrality():
    global _CENTRALITY_DATA
    if _CENTRALITY_DATA is not None:
        return _CENTRALITY_DATA
    path = os.path.join(os.path.dirname(BACKEND_DIR), "data", "reference", "network_centrality.json")
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        _CENTRALITY_DATA = json.load(f)
    return _CENTRALITY_DATA

@app.get("/api/network/centrality")
def get_network_centrality(limit: int = 30, tfs_only: bool = True):
    """Return network centrality metrics for TFs in the regulatory network."""
    data = _load_centrality()
    if data is None:
        raise HTTPException(status_code=503, detail="Centrality data not available. Run scripts/network_centrality.py first.")
    nodes = data.get("nodes", {})
    result = [v for v in nodes.values() if (not tfs_only or v.get("is_tf"))]
    result.sort(key=lambda x: x.get("importance", 0), reverse=True)
    return {
        "_meta": data.get("_meta", {}),
        "top_tfs": result[:limit],
        "total_tfs": sum(1 for v in nodes.values() if v.get("is_tf")),
        "total_nodes": len(nodes),
    }

@app.get("/api/network/centrality/{locus}")
def get_centrality_for_gene(locus: str):
    """Return centrality metrics for a specific gene locus tag."""
    data = _load_centrality()
    if data is None:
        raise HTTPException(status_code=503, detail="Centrality data not available.")
    nodes = data.get("nodes", {})
    locus_lower = locus.strip().lower()
    entry = nodes.get(locus_lower) or nodes.get(locus)
    if entry is None:
        for k, v in nodes.items():
            if k.lower() == locus_lower:
                entry = v
                break
    if entry is None:
        raise HTTPException(status_code=404, detail="Gene not found in centrality data.")
    return entry

'''

MARKER = '    except Exception as e:\n        raise HTTPException(status_code=500, detail=f"Failed to retrieve pruning report: {str(e)}")\n\n@app.get("/api/analysis/string_ppi")'
REPLACEMENT = '    except Exception as e:\n        raise HTTPException(status_code=500, detail=f"Failed to retrieve pruning report: {str(e)}")' + ENDPOINTS_CODE + '@app.get("/api/analysis/string_ppi")'

with open("backend/app.py", "r", encoding="utf-8") as f:
    content = f.read()

if MARKER in content:
    content = content.replace(MARKER, REPLACEMENT, 1)
    with open("backend/app.py", "w", encoding="utf-8") as f:
        f.write(content)
    print("SUCCESS: centrality endpoints added to app.py")
else:
    print("ERROR: marker not found")
    # Debug
    idx = content.find("Failed to retrieve pruning report")
    print(f"  Found 'Failed to retrieve pruning report' at char {idx}")
    print(f"  Context: {repr(content[idx-50:idx+100])}")
