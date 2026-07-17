#!/usr/bin/env python3
"""
run_server.py — Compatibility shim + server entry point
========================================================
Business logic has been refactored into focused sub-modules under backend/:
  - backend/gene_utils.py       : gene identifier utilities
  - backend/kegg_client.py      : KEGG REST cache & pathway loading
  - backend/metabolic_mapper.py : metabolic model parsing
  - backend/bio_handlers.py     : query handler functions
  - backend/sequence_tools.py   : sequence alignment & homolog search

This file re-exports all symbols that app.py depends on so that the
existing `import run_server` calls continue to work without modification.
"""
import os
import sys
import webbrowser
import threading
import time
import json

PORT    = int(os.environ.get("PORT", 8000))
if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    ROOT_DIR = sys._MEIPASS
else:
    ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

# ── Ensure backend/ is importable when run_server.py is executed directly ────
_backend_dir = os.path.join(ROOT_DIR, "backend")
if _backend_dir not in sys.path:
    sys.path.insert(0, _backend_dir)

# ── Re-export from sub-modules (app.py uses `run_server.<symbol>`) ────────────
from backend.gene_utils import (
    get_absolute_path,
    load_gene_mappings,
    normalize_gene_locus,
    expand_gene_aliases,
    split_mapping_values,
    first_row_value,
    infer_model_name_from_file,
    safe_float,
    extract_genes_from_gpr_rule,
    reaction_equation_from_metabolites,
    CG_TO_CGL, CGL_TO_CG, GENE_NAMES, NAME_TO_CG,
)
from backend.kegg_client import (
    load_kegg_cache,
    save_kegg_cache,
    load_kegg_pathway_names,
    load_organism_kegg_links,
    find_matching_kegg_pathways,
    get_gene_pathways_and_go,
    KEGG_PATHWAY_NAMES, KEGG_CACHE_HIT, KEGG_CACHE_LOADED,
    KEGG_CACHE_DIR, KEGG_CACHE_FILE,
    GENE_TO_PATHWAYS, PATHWAY_TO_GENES,
    ORGANISM_PATHWAYS_LOADED, PATHWAY_NAMES_MUTEX,
    GENE_PATHWAYS_CACHE, PATHWAY_REGULATION_CACHE,
)
from backend.metabolic_mapper import (
    find_ecgl1_root,
    load_ecgl1_metadata,
    parse_ecgl1_json_mappings,
    parse_sbml_gene_reaction_mappings,
    load_metabolic_model_mappings,
    METABOLIC_MODEL_DIR, METABOLIC_MODEL_CACHE,
)
from backend.bio_handlers import (
    hypergeom_sf,
    evidence_weight,
    calculate_tf_pathway_impact,
    get_regulatory_targets_for_tf,
    handle_regulon_enrichment,
    handle_go_enrichment,
    handle_metabolic_impact,
    handle_metabolic_pathways,
    handle_imodulon_reactions,
    handle_imodulon_simulation,
    handle_tf_simulation,
    handle_pathway_regulation,
)
from backend.sequence_tools import (
    run_needleman_wunsch,
    handle_homolog_alignment,
)

# RAGService instance has been moved to backend/ai_handlers.py

# Species abbreviations map from CoryneRegNet7 prefixes to user-friendly names
SPECIES_MAP = {
    "B_s": "B. subtilis",
    "E_c": "E. coli",
    "M_t": "M. tuberculosis",
    "C_g": "C. glutamicum",
    "C_a": "C. aurimucosum",
    "C_c": "C. callunae",
    "C_d": "C. diphtheriae",
    "C_e": "C. efficiens",
    "C_f": "C. falsenii",
    "C_h": "C. halotolerans",
    "C_i": "C. imitans",
    "C_j": "C. jeikeium",
    "C_k": "C. kroppenstedtii",
    "C_l": "C. lipophiloflavum",
    "C_m": "C. minutissimum",
    "C_p": "C. pseudotuberculosis",
    "C_r": "C. resistens",
    "C_s": "C. striatum",
    "C_t": "C. tuberculostearicum",
    "C_u": "C. urealyticum",
    "C_v": "C. viteruminis",
    "C_x": "C. xerosis",
    "[_f": "B. flavum",
}



def open_browser():
    # Wait 1 second to make sure the server has started
    time.sleep(1.0)
    url = f"http://localhost:{PORT}/index.html"
    print(f"Opening network explorer at: {url}")
    webbrowser.open(url)

if __name__ == "__main__":
    import uvicorn
    
    # Start browser in a background thread if not in headless mode
    if os.environ.get("HEADLESS", "false").lower() != "true":
        browser_thread = threading.Thread(target=open_browser)
        browser_thread.daemon = True
        browser_thread.start()
        
    print(f"Local Server successfully starting on port {PORT} using FastAPI & Uvicorn...")
    try:
        from backend.app import app
        uvicorn.run(app, host="0.0.0.0", port=PORT, reload=False)
    except KeyboardInterrupt:
        print("\nStopping local server. Goodbye!")
        sys.exit(0)
    except Exception as e:
        print(f"Error starting server: {e}")
        sys.exit(1)
