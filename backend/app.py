from fastapi import FastAPI, HTTPException, Header
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
import logging
import sys
import os
import json
import math

if not hasattr(math, 'comb'):
    def math_comb(n, k):
        if k < 0 or k > n:
            return 0
        if k == 0 or k == n:
            return 1
        k = min(k, n - k)
        numerator = 1
        denominator = 1
        for i in range(1, k + 1):
            numerator *= n - i + 1
            denominator *= i
        return numerator // denominator
    math.comb = math_comb

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
    MFAComparisonResponse,
    PathwayReactionsRequest
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("app")

app = FastAPI(title="Cgl Regulation FBA Simulator API", version="0.3.0")

# Enable CORS for frontend integration across ports
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

ESSENTIAL_GENES = {}
PRODORIC_PWMS = {}
BRENDA_KCAT_MAPPINGS = {}
STRING_INTERACTIONS = {}
ABASY_ROLES = {}

def check_essentiality(gene_id: str):
    """
    Check if a gene locus tag (e.g. cg0001) or its aliases are classified as essential.
    """
    if not gene_id or not ESSENTIAL_GENES:
        return None
    
    # Try direct lookup
    g_lower = gene_id.strip().lower()
    if g_lower in ESSENTIAL_GENES:
        return ESSENTIAL_GENES[g_lower]
        
    # Try resolving aliases
    aliases = run_server.expand_gene_aliases(g_lower)
    for alias in aliases:
        a_lower = alias.lower()
        if a_lower in ESSENTIAL_GENES:
            return ESSENTIAL_GENES[a_lower]
            
    return None

def check_abasy_role(gene_id: str):
    """
    Check if a gene locus tag has an Abasy role classification.
    """
    if not gene_id or not ABASY_ROLES:
        return None
        
    g_lower = gene_id.strip().lower()
    if g_lower in ABASY_ROLES:
        return ABASY_ROLES[g_lower]
        
    # Resolve aliases
    aliases = run_server.expand_gene_aliases(g_lower)
    for alias in aliases:
        a_lower = alias.lower()
        if a_lower in ABASY_ROLES:
            return ABASY_ROLES[a_lower]
            
    return None

@app.on_event("startup")
def startup_event():
    global ESSENTIAL_GENES, PRODORIC_PWMS, BRENDA_KCAT_MAPPINGS, STRING_INTERACTIONS, ABASY_ROLES
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

    # Load essential genes
    try:
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        eg_path = os.path.join(root_dir, "data", "reference", "essential_genes.json")
        if os.path.exists(eg_path):
            with open(eg_path, "r", encoding="utf-8") as f:
                ESSENTIAL_GENES = json.load(f)
            logger.info(f"Loaded {len(ESSENTIAL_GENES)} essential genes.")
        else:
            logger.warning(f"essential_genes.json not found at {eg_path}")
    except Exception as e:
        logger.error(f"Error loading essential_genes.json: {e}")

    # Load PRODORIC motifs
    try:
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        pp_path = os.path.join(root_dir, "data", "reference", "prodoric_pwms.json")
        if os.path.exists(pp_path):
            with open(pp_path, "r", encoding="utf-8") as f:
                PRODORIC_PWMS = json.load(f)
            logger.info(f"Loaded {len(PRODORIC_PWMS)} PRODORIC motifs.")
        else:
            logger.warning(f"prodoric_pwms.json not found at {pp_path}")
    except Exception as e:
        logger.error(f"Error loading prodoric_pwms.json: {e}")

    # Load BRENDA kcat mappings
    try:
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        bk_path = os.path.join(root_dir, "data", "reference", "brenda_kcat_mappings.json")
        if os.path.exists(bk_path):
            with open(bk_path, "r", encoding="utf-8") as f:
                BRENDA_KCAT_MAPPINGS = json.load(f)
            logger.info(f"Loaded {len(BRENDA_KCAT_MAPPINGS)} BRENDA kcat mappings.")
        else:
            logger.warning(f"brenda_kcat_mappings.json not found at {bk_path}")
    except Exception as e:
        logger.error(f"Error loading brenda_kcat_mappings.json: {e}")

    # Load STRING interactions
    try:
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        si_path = os.path.join(root_dir, "data", "reference", "string_interactions.json")
        if os.path.exists(si_path):
            with open(si_path, "r", encoding="utf-8") as f:
                STRING_INTERACTIONS = json.load(f)
            logger.info(f"Loaded {len(STRING_INTERACTIONS)} STRING interactions.")
        else:
            logger.warning(f"string_interactions.json not found at {si_path}")
    except Exception as e:
        logger.error(f"Error loading string_interactions.json: {e}")

    # Load Abasy roles
    try:
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        ab_path = os.path.join(root_dir, "data", "reference", "abasy_roles.json")
        if os.path.exists(ab_path):
            with open(ab_path, "r", encoding="utf-8") as f:
                ABASY_ROLES = json.load(f)
            logger.info(f"Loaded {len(ABASY_ROLES)} Abasy roles.")
        else:
            logger.warning(f"abasy_roles.json not found at {ab_path}")
    except Exception as e:
        logger.error(f"Error loading abasy_roles.json: {e}")

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
            
            # Check reaction essential genes
            essential_genes_in_rxn = []
            is_rxn_essential = False
            
            # Check reaction global regulators (Abasy Atlas)
            global_regulators_in_rxn = []
            has_global_reg = False
            
            for g in rxn.genes:
                eg_info = check_essentiality(g.id)
                if eg_info:
                    essential_genes_in_rxn.append(f"{g.id} ({eg_info.get('gene', '')})")
                    is_rxn_essential = True
                    
                ab_info = check_abasy_role(g.id)
                if ab_info:
                    role = ab_info.get("role", "")
                    if role in ("Global Regulator", "Basal Machinery"):
                        global_regulators_in_rxn.append(f"{g.id} ({role})")
                        has_global_reg = True
                        
            candidates.append({
                "reactionId": rxn.id,
                "name": rxn.name,
                "equation": rxn.reaction,
                "lowerBound": float(rxn.lower_bound),
                "upperBound": float(rxn.upper_bound),
                "classification": classification,
                "confidence": confidence,
                "reason": reason,
                "isEssential": is_rxn_essential,
                "essentialGenes": essential_genes_in_rxn,
                "hasGlobalRegulator": has_global_reg,
                "globalRegulators": global_regulators_in_rxn
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

@app.post("/api/model/pathway/reactions")
def get_pathway_reactions(req: PathwayReactionsRequest):
    try:
        model = load_model_if_needed()
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    
    results = []
    for rxn_id in req.reactionIds:
        # Try original ID first, then strip R_ prefix (pathway data uses R_ but model may not)
        lookup_id = rxn_id
        if rxn_id not in model.reactions and rxn_id.startswith("R_"):
            lookup_id = rxn_id[2:]
        if lookup_id in model.reactions:
            rxn = model.reactions.get_by_id(lookup_id)
            reactants = [m.id for m in rxn.reactants]
            products = [m.id for m in rxn.products]
            results.append({
                "reactionId": rxn_id,  # Return original ID to match frontend expectations
                "name": rxn.name,
                "equation": rxn.reaction,
                "reactants": reactants,
                "products": products
            })
    return results

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
            req.calibrateTimepoint,
            brenda_kcat_mappings=BRENDA_KCAT_MAPPINGS
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
        tf_lower = tf.strip().lower()
        
        # Check PRODORIC_PWMS first
        if tf_lower in PRODORIC_PWMS:
            pwm_data = PRODORIC_PWMS[tf_lower]
            return {
                "tf": tf_lower,
                "tf_name": pwm_data.get("tf_name", tf),
                "consensus": pwm_data.get("consensus", ""),
                "pwm": pwm_data.get("pwm"),
                "nsites": pwm_data.get("targets_count", 0),
                "source": "PRODORIC (Local DB)"
            }
        else:
            # Check by TF name (e.g. "glxr" -> "cg0350")
            for k, v in PRODORIC_PWMS.items():
                if v.get("tf_name", "").lower() == tf_lower:
                    return {
                        "tf": k,
                        "tf_name": v.get("tf_name", tf),
                        "consensus": v.get("consensus", ""),
                        "pwm": v.get("pwm"),
                        "nsites": v.get("targets_count", 0),
                        "source": "PRODORIC (Local DB)"
                    }
                    
        handler_instance = run_server.CustomHTTPRequestHandler.__new__(run_server.CustomHTTPRequestHandler)
        result = handler_instance.perform_motif_prediction(tf)
        return result
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/predict_binding_affinity")
def predict_binding_affinity(tf: str = "", sequence: str = "", temperature: float = 30.0):
    if not tf or not sequence:
        raise HTTPException(status_code=400, detail="Missing tf or sequence parameter")
    try:
        tf_lower = tf.strip().lower()
        pwm = None
        tf_name = tf
        consensus = ""
        targets_count = 0
        
        # Check PRODORIC_PWMS first
        if tf_lower in PRODORIC_PWMS:
            pwm_data = PRODORIC_PWMS[tf_lower]
            pwm = pwm_data.get("pwm")
            tf_name = pwm_data.get("tf_name", tf)
            consensus = pwm_data.get("consensus", "")
            targets_count = pwm_data.get("targets_count", 0)
        else:
            # Check by TF name (e.g. "glxr" -> "cg0350")
            for k, v in PRODORIC_PWMS.items():
                if v.get("tf_name", "").lower() == tf_lower:
                    pwm = v.get("pwm")
                    tf_name = v.get("tf_name", tf)
                    consensus = v.get("consensus", "")
                    targets_count = v.get("targets_count", 0)
                    break
        
        if not pwm:
            # Fallback to dynamic de novo prediction flow
            handler_instance = run_server.CustomHTTPRequestHandler.__new__(run_server.CustomHTTPRequestHandler)
            motif_res = handler_instance.perform_motif_prediction(tf)
            if "error" in motif_res:
                raise HTTPException(status_code=400, detail=motif_res["error"])
            
            pwm = motif_res.get("pwm")
            tf_name = motif_res.get("tf_name", tf)
            consensus = motif_res.get("consensus", "")
            targets_count = motif_res.get("targets_count", 0)
            
        if not pwm:
            raise HTTPException(status_code=400, detail="Could not resolve PWM motif matrix for the TF")
            
        from backend.thermodynamics import scan_sequence_for_affinity
        affinity_res = scan_sequence_for_affinity(pwm, sequence, temperature)
        if "error" in affinity_res:
            raise HTTPException(status_code=400, detail=affinity_res["error"])
            
        return {
            "tf": tf,
            "tf_name": tf_name,
            "consensus": consensus,
            "targets_count": targets_count,
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


# ── Thermodynamic Gene Context API ────────────────────────────────────────────
_THERMO_DATA_CACHE = None

def _load_thermo_gene_data():
    """Load thermo data, cached after first call."""
    global _THERMO_DATA_CACHE
    if _THERMO_DATA_CACHE is not None:
        return _THERMO_DATA_CACHE
    path = os.path.join(os.path.dirname(BACKEND_DIR), "data", "reference", "thermo_dgr_data.json")
    if not os.path.exists(path):
        return {}
    with open(path, "r", encoding="utf-8") as f:
        _THERMO_DATA_CACHE = json.load(f)
    return _THERMO_DATA_CACHE


@app.get("/api/thermo/gene_context")
def get_thermo_gene_context(gene: str = ""):
    """
    Return thermodynamic status for all reactions associated with a given gene.
    Used to annotate ecFBA results with thermodynamic context.
    
    Returns:
      - reactions: list of {id, name, direction_locked, dgr0, dgr_min, dgr_max, confidence}
      - thermo_support_level: 'strong' | 'moderate' | 'weak' | 'none'
      - n_locked: number of this gene's reactions that are direction-locked
      - ko_thermo_confidence: 0.0-1.0 score
    """
    if not gene:
        raise HTTPException(status_code=400, detail="Missing gene parameter")

    thermo_data = _load_thermo_gene_data()
    thermo_rxns = thermo_data.get("reactions", {})

    # Get the loaded model to find gene→reaction associations
    from model_loader import load_model_if_needed
    model = load_model_if_needed()
    if model is None:
        raise HTTPException(status_code=503, detail="Model not loaded")

    # Expand aliases
    gene_lower = gene.strip().lower()
    aliases = set()
    try:
        aliases = run_server.expand_gene_aliases(gene_lower)
    except Exception:
        aliases = {gene_lower, gene.strip()}

    # Find reactions associated with this gene
    gene_reactions = []
    for rxn in model.reactions:
        rxn_gene_ids = {g.id.lower() for g in rxn.genes}
        if aliases & rxn_gene_ids:
            gene_reactions.append(rxn.id)

    if not gene_reactions:
        return {
            "gene": gene,
            "gene_reactions": [],
            "thermo_annotated": [],
            "thermo_support_level": "none",
            "n_locked": 0,
            "ko_thermo_confidence": 0.0,
            "message": "No reactions found for this gene in the metabolic model"
        }

    # Annotate each reaction with thermo data
    annotated = []
    n_locked = 0
    n_confirmed = 0
    confidence_sum = 0.0
    conf_weights = {"HIGH": 1.0, "MED": 0.6, "LOW": 0.3}

    for rxn_id in gene_reactions:
        rxn = model.reactions.get_by_id(rxn_id)
        entry = thermo_rxns.get(rxn_id, {})

        direction_locked = entry.get("direction_locked", "none")
        dgr0 = entry.get("dgr_prime_0")
        dgr_min = entry.get("dgr_prime_min")
        dgr_max = entry.get("dgr_prime_max")
        conf = entry.get("confidence", "NONE")
        conf_weight = conf_weights.get(conf, 0.0)

        is_locked = direction_locked in ("forward", "reverse")
        if is_locked:
            n_locked += 1
            confidence_sum += conf_weight
        elif dgr0 is not None:
            n_confirmed += 1
            confidence_sum += conf_weight * 0.5

        annotated.append({
            "reaction_id": rxn_id,
            "reaction_name": rxn.name or rxn_id,
            "direction_locked": direction_locked,
            "current_lb": rxn.lower_bound,
            "current_ub": rxn.upper_bound,
            "has_thermo_data": dgr0 is not None,
            "dgr_prime_0": dgr0,
            "dgr_prime_min": dgr_min,
            "dgr_prime_max": dgr_max,
            "confidence": conf,
            "thermo_label": (
                f"ΔG'°={dgr0:.1f} kJ/mol" if dgr0 is not None else "No ΔG' data"
            ),
        })

    # Compute overall thermo support level
    n_rxns = len(gene_reactions)
    lock_fraction = n_locked / n_rxns if n_rxns > 0 else 0
    data_fraction = (n_locked + n_confirmed) / n_rxns if n_rxns > 0 else 0
    ko_confidence = min(confidence_sum / max(n_rxns, 1), 1.0)

    if lock_fraction >= 0.5:
        support_level = "strong"
    elif lock_fraction > 0 or data_fraction >= 0.5:
        support_level = "moderate"
    elif data_fraction > 0:
        support_level = "weak"
    else:
        support_level = "none"

    return {
        "gene": gene,
        "gene_reactions": gene_reactions,
        "thermo_annotated": sorted(annotated, key=lambda x: x["direction_locked"] != "none", reverse=True),
        "thermo_support_level": support_level,
        "n_locked": n_locked,
        "n_confirmed": n_confirmed,
        "n_no_data": n_rxns - n_locked - n_confirmed,
        "total_reactions": n_rxns,
        "lock_fraction": round(lock_fraction, 3),
        "ko_thermo_confidence": round(ko_confidence, 3),
    }

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
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/imodulon/reactions")
def imodulon_reactions(imodulon: str = ""):
    try:
        result = run_server.handle_imodulon_reactions(imodulon)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/imodulon/simulation")
def imodulon_simulation(imodulon: str = ""):
    try:
        result = run_server.handle_imodulon_simulation(imodulon)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/engineering/simulation")
def engineering_simulation(tf: str = ""):
    try:
        result = run_server.handle_tf_simulation(tf)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/quality/essential")
def get_essential_genes():
    """Return the database of C. glutamicum essential genes."""
    return ESSENTIAL_GENES

@app.get("/api/quality/brenda")
def get_brenda_mappings():
    """Return the database of C. glutamicum BRENDA kcat mappings."""
    return BRENDA_KCAT_MAPPINGS

@app.get("/api/quality/abasy")
def get_abasy_roles():
    """Return the database of C. glutamicum Abasy roles."""
    return ABASY_ROLES

@app.get("/api/thermo/pruning-report")
def get_thermo_pruning_report():
    """
    Return the thermodynamic directionality pruning report for the loaded model.
    Shows how many reactions had their bounds tightened based on ΔrG' feasibility analysis.
    """
    try:
        from thermo_pruner import get_pruning_report
        report = get_pruning_report()
        return report
    except ImportError:
        return {
            "enabled": False,
            "message": "Thermodynamic pruning module not available."
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to retrieve pruning report: {str(e)}")
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

@app.get("/api/analysis/string_ppi")
def get_string_ppi(gene: str = ""):
    """Return STRING protein-protein interaction partners for a gene locus tag."""
    if not gene:
        raise HTTPException(status_code=400, detail="Missing gene parameter")
    gene_lower = gene.strip().lower()
    
    # Try direct lookup
    partners = STRING_INTERACTIONS.get(gene_lower)
    
    # If not found, try mapping aliases
    if not partners:
        aliases = run_server.expand_gene_aliases(gene_lower)
        for alias in aliases:
            a_lower = alias.lower()
            if a_lower in STRING_INTERACTIONS:
                partners = STRING_INTERACTIONS[a_lower]
                break
                
    if not partners:
        return {"gene": gene, "partners": []}
        
    return {"gene": gene, "partners": partners}

@app.get("/api/quality/icgb21fr")
def quality_icgb21fr():
    """Compute regulatory gene coverage statistics against the iCGB21FR model."""
    import cobra, csv as csv_mod
    try:
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        model_path = os.path.join(root_dir, "data", "reference", "model", "iCGB21FR.xml")
        if not os.path.exists(model_path):
            raise HTTPException(status_code=404, detail=f"iCGB21FR.xml not found at {model_path}")

        model = cobra.io.read_sbml_model(model_path)

        # Build gene-to-reactions and gene-to-subsystems maps (with alias expansion)
        gene_to_rxns: dict = {}
        gene_to_paths: dict = {}
        for rxn in model.reactions:
            subsystem = (rxn.subsystem or "").strip()
            for gene in rxn.genes:
                g_id = gene.id.strip().lower()
                all_ids = run_server.expand_gene_aliases(g_id)
                all_ids.add(g_id)
                for aid in all_ids:
                    gene_to_rxns.setdefault(aid.lower(), set()).add(rxn.id)
                    if subsystem:
                        gene_to_paths.setdefault(aid.lower(), set()).add(subsystem)

        # Load regulatory genes from regulations.csv
        reg_path = os.path.join(root_dir, "data", "reference", "regulations.csv")
        reg_genes_raw: set = set()
        with open(reg_path, "r", encoding="utf-8") as csvf:
            reader = csv_mod.DictReader(csvf)
            for row in reader:
                for field in ("TF_locusTag", "TF_altLocusTag", "TG_locusTag", "TG_altLocusTag", "TF", "Target", "tf_locus", "target_locus"):
                    val = row.get(field, "").strip().lower()
                    if val:
                        reg_genes_raw.add(val)

        # Expand aliases for all regulatory genes
        unique_reg_genes: set = set()
        for rg in reg_genes_raw:
            unique_reg_genes.add(rg)
            for alias in run_server.expand_gene_aliases(rg):
                unique_reg_genes.add(alias.lower())

        # Compute coverage
        mapped_rxn_genes = 0
        mapped_path_genes = 0
        unique_rxns: set = set()
        unique_paths: set = set()
        unmapped = []

        for gene in unique_reg_genes:
            rxns = gene_to_rxns.get(gene, set())
            paths = gene_to_paths.get(gene, set())
            if rxns:
                mapped_rxn_genes += 1
                unique_rxns.update(rxns)
            if paths:
                mapped_path_genes += 1
                unique_paths.update(paths)
            if not rxns:
                unmapped.append(gene)

        unmapped = sorted(set(unmapped))

        return {
            "model_id": model.id,
            "model_genes": len(model.genes),
            "regulatory_gene_count": len(unique_reg_genes),
            "genes_mapped_to_reactions": mapped_rxn_genes,
            "genes_mapped_to_pathways": mapped_path_genes,
            "unique_mapped_reactions": len(unique_rxns),
            "unique_mapped_pathways": len(unique_paths),
            "unmapped_gene_count": len(unmapped),
            "unmapped_genes": unmapped[:100]  # cap for display
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"iCGB21FR quality endpoint failed: {str(e)}")
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
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        folder = os.path.join(root_dir, 'data', 'reference', 'AllOrganismsFiles')
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


