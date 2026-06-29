from fastapi import FastAPI, HTTPException, Header
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
import logging
import sys
import os

# Add backend directory and parent directory to sys.path
BACKEND_DIR = os.path.dirname(os.path.abspath(__file__))
if BACKEND_DIR not in sys.path:
    sys.path.insert(0, BACKEND_DIR)

PARENT_DIR = os.path.dirname(BACKEND_DIR)
if PARENT_DIR not in sys.path:
    sys.path.insert(0, PARENT_DIR)

import run_server

from model_loader import get_model_status, load_model_if_needed
from simulation import run_baseline_simulation, run_gene_knockout, run_gene_set_knockout, run_tf_perturbation, run_fva_analysis, run_dynamic_rfba, run_ecfba_simulation, run_mfa_comparison
from schemas import (
    ModelStatusResponse,
    ReactionSearchResponse,
    BaselineSimulationResponse,
    GeneKnockoutRequest,
    GeneKnockoutResponse,
    GeneSetKnockoutRequest,
    GeneSetKnockoutResponse,
    TFPerturbationRequest,
    TFPerturbationResponse,
    GlutamateCandidatesResponse,
    GlutamateCandidateSchema,
    FVARequest,
    FVAResponse,
    RFBARequest,
    RFBAResponse,
    ECFBARequest,
    ECFBAResponse,
    MFAComparisonResponse
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("app")

app = FastAPI(title="Cgl Regulation FBA Simulator API", version="0.1.0")

# Enable CORS for frontend integration across ports
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.on_event("startup")
def startup_event():
    logger.info("Initializing FBA simulator service...")
    try:
        load_model_if_needed()
    except Exception as e:
        logger.warning(f"Initial model load failed (will retry on demand): {str(e)}")
    
    # Initialize run_server mappings and caches
    try:
        run_server.load_gene_mappings()
        run_server.load_organism_kegg_links()
        logger.info("Successfully loaded gene mappings and KEGG links from run_server.")
    except Exception as e:
        logger.warning(f"Failed to load run_server mappings/caches: {str(e)}")

@app.get("/api/model/status", response_model=ModelStatusResponse)
def model_status():
    status = get_model_status()
    return status

def classify_glutamate_reaction(rxn):
    rxn_id_lower = rxn.id.lower()
    rxn_name_lower = rxn.name.lower()
    rxn_formula = rxn.reaction
    rxn_formula_lower = rxn_formula.lower()
    
    # Identify if extracellular glutamate is involved
    # Extracellular glutamate is usually 'glu__L_e' or similar
    has_extracellular = any('glu__L_e' in met.id or 'glu_e' in met.id.lower() for met in rxn.metabolites)
    has_intracellular = any('glu__L_c' in met.id or 'glu_c' in met.id.lower() for met in rxn.metabolites)
    
    is_exchange = rxn_id_lower.startswith('ex_') or '_ex' in rxn_id_lower
    
    if is_exchange:
        return "exchange", "high", "Reaction ID suggests exchange and equation represents extracellular L-glutamate boundary flux."
    elif 'export' in rxn_name_lower or 'export' in rxn_id_lower or 'secretion' in rxn_name_lower:
        return "export", "high", "Reaction name or equation explicitly suggests extracellular glutamate secretion or export."
    elif 'transport' in rxn_name_lower or (has_extracellular and has_intracellular):
        return "transport", "medium", "Reaction represents transport of L-glutamate across cellular compartments."
    elif has_intracellular and not has_extracellular:
        if 'synth' in rxn_name_lower or 'dehydrogenase' in rxn_name_lower or 'transaminase' in rxn_name_lower:
            return "biosynthesis", "medium", "Intracellular enzymatic reaction converting reactants to L-glutamate."
        elif 'decarboxylase' in rxn_name_lower or 'kinase' in rxn_name_lower or 'synthase' in rxn_name_lower:
            return "consumption", "medium", "Intracellular reaction consuming L-glutamate."
        else:
            return "uncertain", "low", "Intracellular glutamate conversion reaction of uncertain direction."
    else:
        return "uncertain", "low", "Glutamate-associated reaction of uncertain category or compartment."

@app.get("/api/model/reactions/glutamate-candidates", response_model=GlutamateCandidatesResponse)
def get_glutamate_candidates():
    candidates = []
    warnings = []
    
    try:
        model = load_model_if_needed()
    except Exception as e:
        logger.error(f"Failed to load model for candidates list: {str(e)}")
        return {"candidates": [], "warnings": [f"Model offline or missing: {str(e)}"]}
        
    for rxn in model.reactions:
        rxn_id_lower = rxn.id.lower()
        rxn_name_lower = rxn.name.lower()
        
        is_glu_related = (
            'glu' in rxn_id_lower or
            'glutamate' in rxn_name_lower
        )
        
        if not is_glu_related:
            for met in rxn.metabolites:
                if 'glu__l' in met.id.lower() or 'glutamate' in (met.name or "").lower():
                    is_glu_related = True
                    break
                    
        if is_glu_related:
            classification, confidence, reason = classify_glutamate_reaction(rxn)
            candidates.append({
                "reactionId": rxn.id,
                "name": rxn.name,
                "equation": rxn.reaction,
                "lowerBound": float(rxn.lower_bound),
                "upperBound": float(rxn.upper_bound),
                "classification": classification,
                "confidence": confidence,
                "reason": reason
            })
            
    # Check if any exchange/export candidate was found
    has_export_or_exchange = any(c["classification"] in ("exchange", "export") for c in candidates)
    if not has_export_or_exchange:
        warnings.append("No high-confidence glutamate export or exchange reaction was identified in the loaded model. Please select a transport or uncertain candidate for tracking manually.")
        
    return {"candidates": candidates, "warnings": warnings}


@app.get("/api/model/reactions/search", response_model=ReactionSearchResponse)
def search_reactions(q: str = ""):
    query_clean = q.strip().lower()
    matches = []
    
    if not query_clean:
        return {"query": q, "matches": []}
        
    try:
        model = load_model_if_needed()
    except Exception as e:
        logger.warning(f"Reaction search fallback (model offline/missing): {str(e)}")
        return {"query": q, "matches": []}
        
    for rxn in model.reactions:
        rxn_id_lower = rxn.id.lower()
        rxn_name_lower = rxn.name.lower()
        rxn_formula_lower = rxn.reaction.lower()
        
        # Check if query matches ID, name, or formula
        matches_rxn = (
            query_clean in rxn_id_lower or
            query_clean in rxn_name_lower or
            query_clean in rxn_formula_lower
        )
        
        # Check if query matches any metabolites
        if not matches_rxn:
            for met in rxn.metabolites:
                if query_clean in met.id.lower() or query_clean in (met.name or "").lower():
                    matches_rxn = True
                    break
                    
        if matches_rxn:
            matches.append({
                "reactionId": rxn.id,
                "name": rxn.name,
                "equation": rxn.reaction,
                "lowerBound": float(rxn.lower_bound),
                "upperBound": float(rxn.upper_bound),
                "metabolites": [met.id for met in rxn.metabolites]
            })
            if len(matches) >= 100:  # Cap matches
                break
                
    return {"query": q, "matches": matches}

@app.post("/api/simulation/baseline", response_model=BaselineSimulationResponse)
def baseline_simulation():
    try:
        model = load_model_if_needed()
        result = run_baseline_simulation(model)
        return result
    except Exception as e:
        logger.error(f"Baseline FBA failed: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/simulation/gene-knockout", response_model=GeneKnockoutResponse)
def gene_knockout(req: GeneKnockoutRequest):
    try:
        model = load_model_if_needed()
        result = run_gene_knockout(model, req.geneId, req.objective, req.trackReactionIds, req.method)
        return result
    except Exception as e:
        logger.error(f"Gene knockout simulation failed: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/simulation/gene-set-knockout", response_model=GeneSetKnockoutResponse)
def gene_set_knockout(req: GeneSetKnockoutRequest):
    try:
        model = load_model_if_needed()
        result = run_gene_set_knockout(model, req.geneIds, req.objective, req.trackReactionIds, req.method)
        return result
    except Exception as e:
        logger.error(f"Gene set knockout simulation failed: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/simulation/tf-perturbation", response_model=TFPerturbationResponse)
def tf_perturbation(req: TFPerturbationRequest):
    if req.mode != "knockout":
        raise HTTPException(status_code=400, detail="Only 'knockout' perturbation mode is currently supported in v0.1.")
    try:
        model = load_model_if_needed()
        result = run_tf_perturbation(model, req.tfId, req.targetGeneIds, req.objective, req.trackReactionIds, req.method)
        return result
    except Exception as e:
        logger.error(f"TF target perturbation simulation failed: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/simulation/flux-variability", response_model=FVAResponse)
def flux_variability(req: FVARequest):
    try:
        model = load_model_if_needed()
        
        # Resolve knockout list depending on mode
        knockout_genes = []
        if req.mode == "gene-knockout":
            if req.geneId:
                knockout_genes.append(req.geneId)
        elif req.mode == "tf-perturbation":
            if req.targetGeneIds:
                knockout_genes.extend(req.targetGeneIds)
                
        status, fva_ranges, warnings = run_fva_analysis(
            model,
            knockout_genes,
            req.objective,
            req.trackReactionIds,
            req.fractionOfOptimum
        )
        
        return {
            "status": status,
            "fractionOfOptimum": req.fractionOfOptimum,
            "fvaRanges": fva_ranges,
            "warnings": warnings
        }
    except Exception as e:
        logger.error(f"FVA simulation endpoint failed: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/simulation/rfba", response_model=RFBAResponse)
def dynamic_rfba(req: RFBARequest):
    try:
        model = load_model_if_needed()
        result = run_dynamic_rfba(
            model,
            req.tfPerturbations,
            req.initialGlucose,
            req.initialBiomass,
            req.timeSteps
        )
        return result
    except Exception as e:
        logger.error(f"Dynamic rFBA simulation failed: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/simulation/ecfba", response_model=ECFBAResponse)
def ecfba_simulation(req: ECFBARequest):
    try:
        root_dir = os.path.dirname(os.path.dirname(__file__))
        json_model_path = os.path.join(root_dir, "data", "reference", "model", "ecCGL1-main", "ecCGL1-main", "model", "iCW773_irr_enz_constraint.json")
        result = run_ecfba_simulation(
            json_model_path,
            req.proteinPoolLimit,
            req.enzymePerturbations,
            req.targetProduct,
            req.temperature,
            req.calibrateTimepoint
        )
        if isinstance(result, dict):
            if "pool_limit" in result and "poolLimit" not in result:
                result["poolLimit"] = result["pool_limit"]
            if "pool_usage" in result and "poolUsage" not in result:
                result["poolUsage"] = result["pool_usage"]
        return result
    except Exception as e:
        logger.error(f"ec-FBA simulation failed: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

# ── Heat stress analysis API endpoints ────────────────────────────────────────
# Data files are excluded from GitHub via .gitignore (unpublished).
# Locally (where data files exist) these endpoints return full results.
# On GitHub/deployed environments the file is absent → natural 404.
# ──────────────────────────────────────────────────────────────────────────────

@app.get("/api/analysis/rna-seq")
def get_rna_seq_analysis():
    file_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "reference", "rna_seq_analysis_results.json")
    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="Heat stress RNA-Seq analysis data is not publicly available yet. It will be released upon publication.")
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        return data
    except Exception as e:
        logger.error(f"Failed to read RNA-Seq analysis file: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Failed to read analysis results: {str(e)}")

@app.get("/api/analysis/dynamic-grn")
def get_dynamic_grn():
    data = get_rna_seq_analysis()
    return data.get("dynamic_grn", {})

@app.get("/api/analysis/causal-grn")
def get_causal_grn():
    data = get_rna_seq_analysis()
    return data.get("causal_grn", [])

@app.get("/api/analysis/metabolic-coupling")
def get_metabolic_coupling():
    data = get_rna_seq_analysis()
    return data.get("metabolic_coupling", {})

@app.get("/api/analysis/tf-motif-enrichment")
def get_tf_motif_enrichment():
    data = get_rna_seq_analysis()
    return data.get("motif_enrichment", {})

@app.get("/api/summarize")
def summarize(
    gene: str = "",
    name: str = "",
    x_ai_api_key: str = Header(None, alias="X-AI-API-Key"),
    x_gemini_api_key: str = Header(None, alias="X-Gemini-API-Key"),
    x_ai_provider: str = Header("google", alias="X-AI-Provider"),
    x_ai_model: str = Header("", alias="X-AI-Model"),
    x_ai_base_url: str = Header("", alias="X-AI-Base-URL"),
):
    api_key = x_ai_api_key or x_gemini_api_key or ""
    try:
        handler_instance = run_server.CustomHTTPRequestHandler.__new__(run_server.CustomHTTPRequestHandler)
        result = handler_instance.perform_summarize(gene, name, api_key, x_ai_provider, x_ai_model, x_ai_base_url)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/pathway")
def pathway(
    pathway: str = "",
    x_ai_api_key: str = Header(None, alias="X-AI-API-Key"),
    x_gemini_api_key: str = Header(None, alias="X-Gemini-API-Key"),
    x_ai_provider: str = Header("google", alias="X-AI-Provider"),
    x_ai_model: str = Header("", alias="X-AI-Model"),
    x_ai_base_url: str = Header("", alias="X-AI-Base-URL"),
):
    api_key = x_ai_api_key or x_gemini_api_key or ""
    try:
        handler_instance = run_server.CustomHTTPRequestHandler.__new__(run_server.CustomHTTPRequestHandler)
        result = handler_instance.perform_pathway_analysis(pathway, api_key, x_ai_provider, x_ai_model, x_ai_base_url)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/gene_assistant")
def gene_assistant(
    query: str = "",
    x_ai_api_key: str = Header(None, alias="X-AI-API-Key"),
    x_gemini_api_key: str = Header(None, alias="X-Gemini-API-Key"),
    x_ai_provider: str = Header("google", alias="X-AI-Provider"),
    x_ai_model: str = Header("", alias="X-AI-Model"),
    x_ai_base_url: str = Header("", alias="X-AI-Base-URL"),
):
    api_key = x_ai_api_key or x_gemini_api_key or ""
    try:
        handler_instance = run_server.CustomHTTPRequestHandler.__new__(run_server.CustomHTTPRequestHandler)
        result = handler_instance.perform_gene_analysis(query, api_key, x_ai_provider, x_ai_model, x_ai_base_url)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/protein_domain")
def protein_domain(
    gene: str = "",
    x_ai_api_key: str = Header(None, alias="X-AI-API-Key"),
    x_gemini_api_key: str = Header(None, alias="X-Gemini-API-Key"),
    x_ai_provider: str = Header("google", alias="X-AI-Provider"),
    x_ai_model: str = Header("", alias="X-AI-Model"),
    x_ai_base_url: str = Header("", alias="X-AI-Base-URL"),
):
    api_key = x_ai_api_key or x_gemini_api_key or ""
    try:
        handler_instance = run_server.CustomHTTPRequestHandler.__new__(run_server.CustomHTTPRequestHandler)
        result = handler_instance.perform_protein_domain_analysis(gene, api_key, x_ai_provider, x_ai_model, x_ai_base_url)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/binding_site")
def binding_site(
    gene: str = "",
    x_ai_api_key: str = Header(None, alias="X-AI-API-Key"),
    x_gemini_api_key: str = Header(None, alias="X-Gemini-API-Key"),
    x_ai_provider: str = Header("google", alias="X-AI-Provider"),
    x_ai_model: str = Header("", alias="X-AI-Model"),
    x_ai_base_url: str = Header("", alias="X-AI-Base-URL"),
):
    api_key = x_ai_api_key or x_gemini_api_key or ""
    try:
        handler_instance = run_server.CustomHTTPRequestHandler.__new__(run_server.CustomHTTPRequestHandler)
        result = handler_instance.perform_binding_site_analysis(gene, api_key, x_ai_provider, x_ai_model, x_ai_base_url)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/predict_motif")
def predict_motif(tf: str = ""):
    try:
        handler_instance = run_server.CustomHTTPRequestHandler.__new__(run_server.CustomHTTPRequestHandler)
        result = handler_instance.perform_motif_prediction(tf)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/predict_binding_affinity")
def predict_binding_affinity(tf: str = "", sequence: str = "", temperature: float = 30.0):
    if not tf or not sequence:
        raise HTTPException(status_code=400, detail="Missing tf or sequence parameter")
    try:
        handler_instance = run_server.CustomHTTPRequestHandler.__new__(run_server.CustomHTTPRequestHandler)
        motif_res = handler_instance.perform_motif_prediction(tf)
        if "error" in motif_res:
            raise HTTPException(status_code=400, detail=motif_res["error"])
        
        pwm = motif_res.get("pwm")
        if not pwm:
            raise HTTPException(status_code=400, detail="Could not resolve PWM motif matrix for the TF")
            
        from backend.thermodynamics import scan_sequence_for_affinity
        affinity_res = scan_sequence_for_affinity(pwm, sequence, temperature)
        if "error" in affinity_res:
            raise HTTPException(status_code=400, detail=affinity_res["error"])
            
        return {
            "tf": tf,
            "tf_name": motif_res.get("tf_name", tf),
            "consensus": motif_res.get("consensus", ""),
            "targets_count": motif_res.get("targets_count", 0),
            **affinity_res
        }
    except HTTPException as he:
        raise he
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/kegg_pathways")
def kegg_pathways(cg: str = "", cgl: str = ""):
    try:
        result = run_server.get_gene_pathways_and_go(cg, cgl)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/pathway_regulation")
def pathway_regulation(pathway: str = "", query: str = ""):
    target = pathway or query
    try:
        result = run_server.handle_pathway_regulation(target)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/metabolic_impact")
def metabolic_impact(gene: str = "", query: str = ""):
    target = gene or query
    try:
        result = run_server.handle_metabolic_impact(target)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/metabolic_pathways")
def metabolic_pathways(pathway: str = "", query: str = ""):
    target = pathway or query
    try:
        result = run_server.handle_metabolic_pathways(target)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/test_ai")
def test_ai(
    x_ai_api_key: str = Header(None, alias="X-AI-API-Key"),
    x_gemini_api_key: str = Header(None, alias="X-Gemini-API-Key"),
    x_ai_provider: str = Header("google", alias="X-AI-Provider"),
    x_ai_model: str = Header("", alias="X-AI-Model"),
    x_ai_base_url: str = Header("", alias="X-AI-Base-URL"),
):
    api_key = x_ai_api_key or x_gemini_api_key or ""
    try:
        handler_instance = run_server.CustomHTTPRequestHandler.__new__(run_server.CustomHTTPRequestHandler)
        prompt = "Hello! Please return a single word: Success."
        response = handler_instance.call_llm_api(prompt, x_ai_provider, api_key, x_ai_model, x_ai_base_url)
        return {"status": "success", "message": f"连接成功！AI 响应：{response}"}
    except Exception as e:
        return {"status": "error", "message": str(e)}

@app.get("/api/regulon_enrichment")
def regulon_enrichment(tf: str = ""):
    try:
        result = run_server.handle_regulon_enrichment(tf)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/homolog_alignment")
def homolog_alignment(gene_name: str = "", accession: str = ""):
    try:
        result = run_server.handle_homolog_alignment(gene_name, accession)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/list_organisms")
def list_organisms():
    try:
        organisms = []
        folder = os.path.join(os.getcwd(), 'data', 'reference', 'AllOrganismsFiles')
        if os.path.exists(folder):
            for filename in os.listdir(folder):
                if filename.endswith('_regulations.csv'):
                    org_id = filename[:-16]
                    if not org_id:
                        continue
                    name = org_id
                    parts = org_id.split('_', 2)
                    if len(parts) >= 2:
                        key = f"{parts[0]}_{parts[1]}"
                        rest = parts[2] if len(parts) > 2 else ""
                        if key in run_server.SPECIES_MAP:
                            clean_rest = rest.replace('_', ' ').strip()
                            name = f"{run_server.SPECIES_MAP[key]} {clean_rest}".strip()
                        else:
                            name = org_id.replace('_', ' ')
                    else:
                        name = org_id.replace('_', ' ')
                    rna_file = f"{org_id}_rna_regulation.csv"
                    has_rna = os.path.exists(os.path.join(folder, rna_file))
                    organisms.append({
                        "id": org_id,
                        "name": name,
                        "has_rna": has_rna
                    })
        
        def sort_key(x):
            is_default = (x['id'] == 'C_g_DSM_20300_=_ATCC_13032')
            return (not is_default, x['name'])
        organisms.sort(key=sort_key)
        return organisms
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/simulation/mfa_comparison", response_model=MFAComparisonResponse)
async def mfa_comparison_endpoint():
    """
    Compare FBA-simulated fluxes against published 13C-MFA literature values
    for C. glutamicum ATCC 13032 wild-type (aerobic glucose minimal medium).
    Reference: Cheng et al. 2017, Becker & Wittmann 2011.
    """
    try:
        import run_server
        model_path = None
        try:
            model_path = run_server.get_model_json_path()
        except AttributeError:
            pass
            
        if not model_path:
            root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            model_path = os.path.join(root_dir, "data", "reference", "model", "ecCGL1-main", "ecCGL1-main", "model", "iCW773_irr_enz_constraint.json")
            if not os.path.exists(model_path):
                model_path = os.path.join(root_dir, "data", "reference", "model", "ecCGL1", "model", "iCW773_irr_enz_constraint.json")

        if not model_path or not os.path.exists(model_path):
            raise HTTPException(status_code=503, detail=f"No model found at path: {model_path}")
            
        result = run_mfa_comparison(model_path)
        return result
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# Mount static files
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
web_dir = os.path.join(ROOT_DIR, "web")
data_dir = os.path.join(ROOT_DIR, "data", "reference")

if os.path.exists(data_dir):
    app.mount("/data", StaticFiles(directory=data_dir), name="data")

if os.path.exists(web_dir):
    app.mount("/", StaticFiles(directory=web_dir, html=True), name="static")


