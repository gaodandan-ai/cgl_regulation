from fastapi import FastAPI, HTTPException, Header, Response
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from contextlib import asynccontextmanager
import logging
import sys
import os
import json

from ai_handlers import (
    perform_summarize,
    perform_pathway_analysis,
    perform_gene_analysis,
    perform_protein_domain_analysis,
    perform_binding_site_analysis,
    perform_motif_prediction,
    call_llm_api
)

if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    BACKEND_DIR = os.path.join(sys._MEIPASS, "backend")
    PARENT_DIR = sys._MEIPASS
else:
    BACKEND_DIR = os.path.dirname(os.path.abspath(__file__))
    PARENT_DIR = os.path.dirname(BACKEND_DIR)

if BACKEND_DIR not in sys.path:
    sys.path.insert(0, BACKEND_DIR)
if PARENT_DIR not in sys.path:
    sys.path.insert(0, PARENT_DIR)

import run_server

from model_loader import get_model_status, load_model_if_needed
try:
    from thermo_pruner import get_pruning_report as _get_thermo_pruning_report
    _THERMO_PRUNER_AVAILABLE = True
except ImportError:
    _get_thermo_pruning_report = None
    _THERMO_PRUNER_AVAILABLE = False
from simulation import run_baseline_simulation, run_gene_knockout, run_gene_set_knockout, run_tf_perturbation, run_fva_analysis, run_dynamic_rfba, run_dynamic_recfba, run_ecfba_simulation, run_mfa_comparison
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
    RECFBARequest,
    RECFBAResponse,
    ECFBARequest,
    ECFBAResponse,
    MFAComparisonResponse,
    PathwayReactionsRequest
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("app")

# Mappings and cache directories will be initialized in lifespan startup event below.

ESSENTIAL_GENES = {}
PRODORIC_PWMS = {}
BRENDA_KCAT_MAPPINGS = {}
STRING_INTERACTIONS = {}
ABASY_ROLES = {}
RHEA_MAPPINGS = {}
CHEBI_MAPPINGS = {}
COG_ANNOTATIONS = {}

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

@asynccontextmanager
async def lifespan(app: FastAPI):
    global ESSENTIAL_GENES, PRODORIC_PWMS, BRENDA_KCAT_MAPPINGS, STRING_INTERACTIONS, ABASY_ROLES, RHEA_MAPPINGS, CHEBI_MAPPINGS, COG_ANNOTATIONS
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

    # Load Rhea mappings
    try:
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        rhea_path = os.path.join(root_dir, "data", "reference", "rhea_mappings.json")
        if os.path.exists(rhea_path):
            with open(rhea_path, "r", encoding="utf-8") as f:
                RHEA_MAPPINGS = json.load(f)
            logger.info(f"Loaded {len(RHEA_MAPPINGS)} Rhea mappings.")
        else:
            logger.warning(f"rhea_mappings.json not found at {rhea_path}")
    except Exception as e:
        logger.error(f"Error loading rhea_mappings.json: {e}")

    # Load ChEBI mappings
    try:
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        chebi_path = os.path.join(root_dir, "data", "reference", "chebi_mappings.json")
        if os.path.exists(chebi_path):
            with open(chebi_path, "r", encoding="utf-8") as f:
                CHEBI_MAPPINGS = json.load(f)
            logger.info(f"Loaded {len(CHEBI_MAPPINGS)} ChEBI mappings.")
        else:
            logger.warning(f"chebi_mappings.json not found at {chebi_path}")
    except Exception as e:
        logger.error(f"Error loading chebi_mappings.json: {e}")

    # Load COG annotations
    try:
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        cog_path = os.path.join(root_dir, "data", "reference", "cog_annotations.json")
        if os.path.exists(cog_path):
            with open(cog_path, "r", encoding="utf-8") as f:
                COG_ANNOTATIONS = json.load(f)
            logger.info(f"Loaded {len(COG_ANNOTATIONS)} COG annotations.")
        else:
            logger.warning(f"cog_annotations.json not found at {cog_path}")
    except Exception as e:
        logger.error(f"Error loading cog_annotations.json: {e}")

    yield

app = FastAPI(title="Cgl Regulation FBA Simulator API", version="0.5.0", lifespan=lifespan)

# Enable CORS for frontend integration across ports
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Enable Gzip compression to speed up transfer of large datasets (like 6MB metabolic mapping JSON)
app.add_middleware(GZipMiddleware, minimum_size=1000)

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
            db_links = RHEA_MAPPINGS.get(rxn.id) or RHEA_MAPPINGS.get("R_" + rxn.id)
            met_links = {}
            for met in rxn.metabolites:
                m_details = CHEBI_MAPPINGS.get(met.id) or CHEBI_MAPPINGS.get("M_" + met.id)
                if m_details:
                    met_links[met.id] = m_details
                    
            matches.append({
                "reactionId": rxn.id,
                "name": rxn.name,
                "equation": rxn.reaction,
                "lowerBound": float(rxn.lower_bound),
                "upperBound": float(rxn.upper_bound),
                "metabolites": [met.id for met in rxn.metabolites],
                "databaseLinks": db_links,
                "metaboliteLinks": met_links
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

@app.post("/api/simulation/recfba", response_model=RECFBAResponse)
def dynamic_recfba(req: RECFBARequest):
    try:
        root_dir = os.path.dirname(os.path.dirname(__file__))
        json_model_path = os.path.join(root_dir, "data", "reference", "model", "ecCGL1-main", "ecCGL1-main", "model", "iCW773_irr_enz_constraint.json")
        result = run_dynamic_recfba(
            json_model_path,
            req.tfPerturbations,
            req.proteinPoolLimit,
            req.temperature,
            req.initialGlucose,
            req.initialBiomass,
            req.timeSteps,
            brenda_kcat_mappings=BRENDA_KCAT_MAPPINGS
        )
        return result
    except Exception as e:
        logger.error(f"Dynamic recFBA simulation failed: {str(e)}")
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
    if not gene or not gene.strip():
        raise HTTPException(status_code=400, detail="Missing required query parameter: gene")
    api_key = x_ai_api_key or x_gemini_api_key or ""
    try:
        result = perform_summarize(gene, name, api_key, x_ai_provider, x_ai_model, x_ai_base_url)
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
        result = perform_pathway_analysis(pathway, api_key, x_ai_provider, x_ai_model, x_ai_base_url)
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
        result = perform_gene_analysis(query, api_key, x_ai_provider, x_ai_model, x_ai_base_url)
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
        result = perform_protein_domain_analysis(gene, api_key, x_ai_provider, x_ai_model, x_ai_base_url)
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
        result = perform_binding_site_analysis(gene, api_key, x_ai_provider, x_ai_model, x_ai_base_url)
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
                    
        result = perform_motif_prediction(tf)
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
            motif_res = perform_motif_prediction(tf)
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
        
        # Enrich result with uniprot_id, cog_id, ec_numbers
        result["uniprot_id"] = ""
        result["cog_id"] = ""
        result["ec_numbers"] = []
        
        # 1. UniProt ID lookup
        from gene_utils import GENE_TO_UNIPROT
        for tag in [cg, cgl]:
            if tag and tag.lower() in GENE_TO_UNIPROT:
                result["uniprot_id"] = GENE_TO_UNIPROT[tag.lower()]
                break
                
        # 2. COG ID lookup
        for tag in [cg, cgl]:
            if tag and tag.lower() in COG_ANNOTATIONS:
                result["cog_id"] = COG_ANNOTATIONS[tag.lower()].get("cog_id", "")
                break
                
        # 3. EC numbers lookup from metabolic model mapping
        try:
            from metabolic_mapper import load_metabolic_model_mappings
            mapping = load_metabolic_model_mappings()
            ec_set = set()
            from gene_utils import expand_gene_aliases
            for tag in [cg, cgl]:
                if tag:
                    for alias in expand_gene_aliases(tag):
                        reactions = mapping.get("gene_to_reactions", {}).get(alias, [])
                        for r in reactions:
                            ec = r.get("ec_number")
                            if ec:
                                ec_set.add(ec)
            result["ec_numbers"] = sorted(list(ec_set))
        except Exception as ex:
            logger.warning(f"Error lookup EC numbers: {ex}")
            
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
            "n_confirmed": 0,
            "n_no_data": 0,
            "total_reactions": 0,
            "lock_fraction": 0.0,
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
def metabolic_pathways(response: Response, pathway: str = "", query: str = ""):
    target = pathway or query
    try:
        result = run_server.handle_metabolic_pathways(target)
        if not target:
            # Cache default large pathways payload (6MB) for 1 hour to improve load speeds
            response.headers["Cache-Control"] = "public, max-age=3600"
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
        prompt = "Hello! Please return a single word: Success."
        response = call_llm_api(prompt, x_ai_provider, api_key, x_ai_model, x_ai_base_url)
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

@app.get("/api/go_enrichment")
def go_enrichment(tf: str = ""):
    """
    GO term enrichment analysis for a TF's target regulon.
    Fetches per-gene GO annotations from UniProt (cached), aggregates
    term frequencies, and applies hypergeometric test + BH FDR correction.
    Results are partitioned by GO namespace: BP / MF / CC.
    """
    try:
        result = run_server.handle_go_enrichment(tf)
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
def get_essential_genes(response: Response):
    """Return the database of C. glutamicum essential genes."""
    response.headers["Cache-Control"] = "public, max-age=3600"
    return ESSENTIAL_GENES

@app.get("/api/quality/brenda")
def get_brenda_mappings(response: Response):
    """Return the database of C. glutamicum BRENDA kcat mappings."""
    response.headers["Cache-Control"] = "public, max-age=3600"
    return BRENDA_KCAT_MAPPINGS

@app.get("/api/quality/abasy")
def get_abasy_roles(response: Response):
    """Return the database of C. glutamicum Abasy roles."""
    response.headers["Cache-Control"] = "public, max-age=3600"
    return ABASY_ROLES

@app.get("/api/quality/cog")
def get_cog_annotations(response: Response):
    """Return the database of C. glutamicum COG functional annotations."""
    response.headers["Cache-Control"] = "public, max-age=3600"
    return COG_ANNOTATIONS

@app.get("/api/thermo/pruning-report")
def get_thermo_pruning_report(response: Response):
    """
    Return the thermodynamic directionality pruning report for the loaded model.
    Shows how many reactions had their bounds tightened based on ΔrG' feasibility analysis.
    Triggers model loading if not already done (pruning runs at load time).
    """
    response.headers["Cache-Control"] = "public, max-age=1800"
    if not _THERMO_PRUNER_AVAILABLE or _get_thermo_pruning_report is None:
        return {
            "enabled": False,
            "message": "Thermodynamic pruning module not available."
        }
    try:
        # Trigger model loading if not yet done (pruning runs during load)
        try:
            load_model_if_needed()
        except Exception:
            pass  # Model load may fail (missing file) — still return what we have

        report = _get_thermo_pruning_report()

        # If pruning ran successfully, return the live report
        if report.get("enabled"):
            return report

        # Fallback: model not yet loaded or pruning not run — read from pre-built JSON
        thermo_path = os.path.join(os.path.dirname(BACKEND_DIR), "data", "reference", "thermo_dgr_data.json")
        if os.path.exists(thermo_path):
            with open(thermo_path, "r", encoding="utf-8") as f:
                raw = json.load(f)
            meta = raw.get("_meta", {})
            reactions = {k: v for k, v in raw.items() if not k.startswith("_") and isinstance(v, dict)}
            n_fwd = meta.get("n_forward_locked", 0)
            n_rev = meta.get("n_reverse_locked", 0)
            n_near_eq = meta.get("n_near_equilibrium", 0)
            n_no_data = meta.get("n_no_data", 0)
            total = meta.get("total_reactions", len(reactions))
            return {
                "enabled": True,
                "total_reactions": total,
                "data_coverage_pct": meta.get("coverage_pct", 0.0),
                "n_pruned": n_fwd + n_rev,
                "n_forward_locked": n_fwd,
                "n_reverse_locked": n_rev,
                "n_confirmed_forward": 0,
                "n_confirmed_reverse": 0,
                "n_thermo_consistent": 0,
                "n_near_equilibrium": n_near_eq,
                "n_no_data": n_no_data,
                "epsilon_kJ": meta.get("epsilon_kJ", 1.0),
                "conditions": meta.get("conditions", "pH 7.0, I=0.1 M, T=30°C"),
                "top_pruned": [],
                "confirmed_reactions": [],
                "all_pruned_count": n_fwd + n_rev,
                "message": "Pre-built thermo data (model not yet loaded for live pruning)",
            }

        return report

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
def get_network_centrality(response: Response, limit: int = 30, tfs_only: bool = True):
    """Return network centrality metrics for TFs in the regulatory network."""
    response.headers["Cache-Control"] = "public, max-age=1800"
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
def get_centrality_for_gene(response: Response, locus: str):
    """Return centrality metrics for a specific gene locus tag."""
    response.headers["Cache-Control"] = "public, max-age=1800"
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
def get_string_ppi(gene: str = "", min_score: int = 400, limit: int = 50):
    """Return STRING v12 protein-protein interaction partners for a gene locus tag.

    Parameters
    ----------
    gene      : CG locus tag (e.g. 'cg0001') or gene name (e.g. 'dnaA')
    min_score : minimum combined STRING score 0-1000 (default 400 = medium confidence)
    limit     : max number of partners to return (default 50)
    """
    if not gene:
        raise HTTPException(status_code=400, detail="Missing gene parameter")
    gene_lower = gene.strip().lower()
    if gene_lower == "_meta":
        raise HTTPException(status_code=400, detail="Invalid gene identifier")

    # Try direct lookup
    partners = STRING_INTERACTIONS.get(gene_lower)

    # If not found, try mapping aliases (e.g. gene name → locus)
    if not partners:
        try:
            aliases = run_server.expand_gene_aliases(gene_lower)
            for alias in aliases:
                a_lower = alias.lower()
                if a_lower in STRING_INTERACTIONS:
                    partners = STRING_INTERACTIONS[a_lower]
                    gene_lower = a_lower
                    break
        except Exception:
            pass

    if not partners:
        return {
            "gene":        gene,
            "resolved_id": gene_lower,
            "partners":    [],
            "total":       0,
            "filtered":    0,
        }

    # Apply score filter and limit
    filtered = [p for p in partners if p.get("score", 0) >= min_score]
    total    = len(filtered)
    filtered = filtered[:limit]  # already sorted by score desc

    # Attach STRING meta (gene count / edge count) from _meta block
    meta = STRING_INTERACTIONS.get("_meta", {})

    return {
        "gene":          gene,
        "resolved_id":   gene_lower,
        "partners":      filtered,
        "total":         total,
        "filtered":      total,
        "string_meta": {
            "version":       meta.get("version", "12.0"),
            "min_score":     min_score,
            "n_genes":       meta.get("n_genes", 0),
            "n_edges":       meta.get("n_edges", 0),
            "n_high_conf":   meta.get("n_high_conf", 0),
        },
    }


@app.get("/api/analysis/string_ppi/neighborhood")
def get_ppi_neighborhood(
    genes: str = "",
    min_score: int = 400,
    limit_per_gene: int = 30,
):
    """Return a merged PPI subgraph for one or more genes (comma-separated).

    Useful for rendering multi-gene PPI networks with node-expansion support.
    Each edge includes all 7 STRING channel scores.

    Parameters
    ----------
    genes           : comma-separated cg locus tags or gene names
    min_score       : minimum combined STRING score 0-1000 (default 400)
    limit_per_gene  : max partners to include per seed gene (default 30)
    """
    if not genes:
        raise HTTPException(status_code=400, detail="Missing genes parameter")

    gene_list = [g.strip().lower() for g in genes.split(",") if g.strip()]
    if not gene_list:
        raise HTTPException(status_code=400, detail="No valid gene identifiers provided")

    nodes: dict[str, dict] = {}   # id → {id, name, degree, is_seed}
    edges: dict[str, dict] = {}   # "src::tgt" → edge dict

    def resolve_locus(raw: str) -> str:
        """Return the canonical cg locus tag for a gene name/alias."""
        if raw in STRING_INTERACTIONS:
            return raw
        try:
            for alias in run_server.expand_gene_aliases(raw):
                al = alias.lower()
                if al in STRING_INTERACTIONS:
                    return al
        except Exception:
            pass
        return raw

    def gene_display_name(locus: str) -> str:
        try:
            name = run_server.normalize_gene_locus(locus)
            return name if name and name != locus else locus
        except Exception:
            return locus

    def add_node(locus: str, is_seed: bool = False):
        if locus not in nodes:
            nodes[locus] = {
                "id": locus,
                "name": gene_display_name(locus),
                "degree": 0,
                "is_seed": is_seed,
            }
        if is_seed:
            nodes[locus]["is_seed"] = True

    def add_edge(src: str, tgt: str, pdata: dict):
        # Normalise direction so A::B and B::A map to the same key
        key = f"{min(src, tgt)}::{max(src, tgt)}"
        if key not in edges:
            edges[key] = {
                "id": f"ppi-{src}-{tgt}",
                "source": src,
                "target": tgt,
                "score": pdata.get("score", 0),
                "experimental":  pdata.get("experimental", 0),
                "database":      pdata.get("database", 0),
                "coexpression":  pdata.get("coexpression", 0),
                "neighborhood":  pdata.get("neighborhood", 0),
                "cooccurrence":  pdata.get("cooccurrence", 0),
                "textmining":    pdata.get("textmining", 0),
                "fusion":        pdata.get("fusion", 0),
                "type":          pdata.get("type", ""),
            }
        else:
            # Keep higher score if duplicate
            if pdata.get("score", 0) > edges[key]["score"]:
                edges[key]["score"] = pdata["score"]
        # Track degree
        nodes[src]["degree"] = nodes[src].get("degree", 0) + 1
        nodes[tgt]["degree"] = nodes[tgt].get("degree", 0) + 1

    # Process each seed gene
    for raw_gene in gene_list:
        locus = resolve_locus(raw_gene)
        add_node(locus, is_seed=True)

        partners = STRING_INTERACTIONS.get(locus, [])
        filtered = [p for p in partners if p.get("score", 0) >= min_score]
        filtered = filtered[:limit_per_gene]

        for p in filtered:
            partner_id = p["partner"].lower()
            add_node(partner_id, is_seed=False)
            add_edge(locus, partner_id, p)

    meta = STRING_INTERACTIONS.get("_meta", {})
    return {
        "genes":   gene_list,
        "nodes":   list(nodes.values()),
        "edges":   list(edges.values()),
        "total_nodes": len(nodes),
        "total_edges": len(edges),
        "string_meta": {
            "version":     meta.get("version", "12.0"),
            "min_score":   min_score,
            "n_genes":     meta.get("n_genes", 0),
            "n_edges":     meta.get("n_edges", 0),
        },
    }


@app.get("/api/analysis/cross_network")
def get_cross_network(gene: str = "", min_ppi_score: int = 400):
    """Cross-network analysis: find regulatory targets that also physically interact with the TF.

    Returns three sets:
    - regulatory_targets : genes regulated by `gene`
    - ppi_partners       : proteins physically interacting with `gene` above min_ppi_score
    - cross_links        : the intersection — dual relationship genes
    """
    if not gene:
        raise HTTPException(status_code=400, detail="'gene' parameter is required")

    # 1. Regulatory targets
    is_tf, reg_targets = run_server.get_regulatory_targets_for_tf(gene)

    # Build locus-indexed lookup (handles gene names)
    reg_by_locus = {}
    for t in reg_targets:
        locus = (t.get("locus") or "").lower()
        if locus:
            reg_by_locus[locus] = t

    # 2. PPI partners
    # Resolve the query gene to the STRING key
    q_lower = gene.strip().lower()
    ppi_key = q_lower
    if ppi_key not in STRING_INTERACTIONS:
        for alias in run_server.expand_gene_aliases(q_lower):
            if alias in STRING_INTERACTIONS:
                ppi_key = alias
                break

    ppi_raw = STRING_INTERACTIONS.get(ppi_key, [])
    ppi_partners = [p for p in ppi_raw if p.get("score", 0) >= min_ppi_score]

    # Index ppi partners by their locus-like key
    ppi_by_locus = {}
    for p in ppi_partners:
        partner_id = (p.get("partner") or "").lower()
        # Try to map partner STRING ID to a cg locus tag
        try:
            canonical = run_server.normalize_gene_locus(partner_id)
        except Exception:
            canonical = partner_id
        key = canonical or partner_id
        ppi_by_locus[key] = dict(p, canonical=key)

    # 3. Cross-links: regulatory targets that appear in PPI partners
    cross_links = []
    for locus, reg_info in reg_by_locus.items():
        # Try both the locus and any aliases
        ppi_hit = ppi_by_locus.get(locus)
        if ppi_hit is None:
            for alias in run_server.expand_gene_aliases(locus):
                if alias in ppi_by_locus:
                    ppi_hit = ppi_by_locus[alias]
                    break
        if ppi_hit:
            cross_links.append({
                "gene":            locus,
                "name":            reg_info.get("name", locus),
                "regulation_role": reg_info.get("regulation", reg_info.get("role", "?")),
                "evidence":        reg_info.get("evidence", ""),
                "ppi_score":       ppi_hit.get("score", 0),
                "ppi_experimental":ppi_hit.get("experimental", 0),
                "ppi_database":    ppi_hit.get("database", 0),
                "ppi_coexpression":ppi_hit.get("coexpression", 0),
            })

    cross_links.sort(key=lambda x: x["ppi_score"], reverse=True)

    return {
        "query":              gene,
        "is_tf":              is_tf,
        "n_regulatory":       len(reg_targets),
        "n_ppi_partners":     len(ppi_partners),
        "n_cross_links":      len(cross_links),
        "cross_links":        cross_links,
        "regulatory_targets": reg_targets,
        "ppi_summary":        [{"gene": p.get("partner",""), "score": p.get("score",0)} for p in ppi_partners[:20]],
    }


@app.get("/api/analysis/tf_similarity")
def get_tf_similarity(min_targets: int = 3, metric: str = "jaccard", top_n: int = 40):
    """Compute pairwise TF co-regulation similarity.

    Builds a regulon (set of target genes) for each TF from regulations.csv,
    then computes pairwise Jaccard similarity: |A∩B| / |A∪B|.

    Returns the matrix rows, TF labels, and shared-target details for heatmap rendering.

    Parameters
    ----------
    min_targets : minimum regulon size to include a TF (default 3)
    metric      : 'jaccard' (default)
    top_n       : max number of TFs to include (sorted by regulon size desc, default 40)
    """
    import csv

    path = run_server.get_absolute_path("data/reference/regulations.csv")
    if not os.path.exists(path):
        raise HTTPException(status_code=404, detail="regulations.csv not found")

    # Build TF → targets map
    tf_targets: dict = {}
    tf_names: dict = {}
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            tf_locus = (row.get("TF_locusTag") or "").strip()
            tf_name  = (row.get("TF_name") or "").strip()
            tg_locus = (row.get("TG_locusTag") or "").strip()
            if not tf_locus or not tg_locus:
                continue
            tf_targets.setdefault(tf_locus, set()).add(tg_locus)
            if tf_name and tf_name != tf_locus:
                tf_names[tf_locus] = tf_name

    # Filter by min_targets and pick top_n by regulon size
    qualified = [(tf, tgts) for tf, tgts in tf_targets.items() if len(tgts) >= min_targets]
    qualified.sort(key=lambda x: len(x[1]), reverse=True)
    qualified = qualified[:top_n]

    tfs     = [q[0] for q in qualified]
    targets = [q[1] for q in qualified]

    # Pairwise Jaccard
    n = len(tfs)
    matrix = []
    shared_detail = {}

    for i in range(n):
        row_vals = []
        for j in range(n):
            if i == j:
                row_vals.append(1.0)
            else:
                inter = targets[i] & targets[j]
                union = targets[i] | targets[j]
                jac = round(len(inter) / len(union), 4) if union else 0.0
                row_vals.append(jac)
                if jac > 0 and i < j:
                    key = f"{tfs[i]}__{tfs[j]}"
                    shared_detail[key] = sorted(inter)
        matrix.append(row_vals)

    tf_labels = [{"locus": tf, "name": tf_names.get(tf, tf), "n_targets": len(targets[i])} for i, tf in enumerate(tfs)]

    return {
        "n_tfs":   n,
        "tfs":     tf_labels,
        "matrix":  matrix,
        "shared":  shared_detail,
        "metric":  metric,
    }


@app.get("/api/analysis/string_ppi/shortest_path")
def get_ppi_shortest_path(
    source: str = "",
    target: str = "",
    min_score: int = 400,
    max_hops: int = 6,
):
    """Find the shortest protein–protein interaction path between two genes.

    Uses BFS on the STRING v12 interaction graph.

    Parameters
    ----------
    source    : start gene (cg locus tag or gene name)
    target    : end gene (cg locus tag or gene name)
    min_score : minimum combined STRING score (default 400)
    max_hops  : maximum path length to search (default 6)
    """
    if not source or not target:
        raise HTTPException(status_code=400, detail="Both 'source' and 'target' parameters are required")

    def resolve(raw: str) -> str:
        r = raw.strip().lower()
        if r in STRING_INTERACTIONS:
            return r
        try:
            for alias in run_server.expand_gene_aliases(r):
                al = alias.lower()
                if al in STRING_INTERACTIONS:
                    return al
        except Exception:
            pass
        return r

    def display_name(locus: str) -> str:
        try:
            n = run_server.normalize_gene_locus(locus)
            return n if n and n != locus else locus
        except Exception:
            return locus

    src_id = resolve(source)
    tgt_id = resolve(target)

    if src_id == tgt_id:
        return {"found": False, "message": "Source and target are the same gene", "path": [], "edges": []}

    if src_id not in STRING_INTERACTIONS:
        return {"found": False, "message": f"No PPI data for source '{source}'", "path": [], "edges": []}
    if tgt_id not in STRING_INTERACTIONS:
        return {"found": False, "message": f"No PPI data for target '{target}'", "path": [], "edges": []}

    # BFS
    from collections import deque
    queue    = deque([[src_id]])
    visited  = {src_id}

    while queue:
        path = queue.popleft()
        current = path[-1]
        if len(path) > max_hops + 1:
            break
        partners = [p for p in STRING_INTERACTIONS.get(current, []) if p.get("score", 0) >= min_score]
        for p in partners:
            pid = p["partner"].lower()
            if pid in visited:
                continue
            new_path = path + [pid]
            if pid == tgt_id:
                # Build response
                nodes = [{"id": n, "name": display_name(n), "is_seed": n in (src_id, tgt_id)} for n in new_path]
                edges = []
                for i in range(len(new_path) - 1):
                    a, b = new_path[i], new_path[i + 1]
                    partners_a = STRING_INTERACTIONS.get(a, [])
                    edge_data = next((x for x in partners_a if x["partner"].lower() == b), None)
                    if edge_data:
                        edges.append({
                            "id": f"path-{a}-{b}",
                            "source": a,
                            "target": b,
                            "score": edge_data.get("score", 0),
                            "experimental":  edge_data.get("experimental", 0),
                            "database":      edge_data.get("database", 0),
                            "coexpression":  edge_data.get("coexpression", 0),
                            "textmining":    edge_data.get("textmining", 0),
                            "neighborhood":  edge_data.get("neighborhood", 0),
                            "cooccurrence":  edge_data.get("cooccurrence", 0),
                            "fusion":        edge_data.get("fusion", 0),
                        })
                return {
                    "found": True,
                    "hops": len(new_path) - 1,
                    "path_ids": new_path,
                    "nodes": nodes,
                    "edges": edges,
                    "source": src_id,
                    "target": tgt_id,
                }
            visited.add(pid)
            queue.append(new_path)

    return {"found": False, "message": f"No path found within {max_hops} hops", "path": [], "edges": []}


@app.get("/api/analysis/string_ppi/hub_ranking")
def get_ppi_hub_ranking(min_score: int = 700, limit: int = 50):
    """Return the top hub proteins ranked by PPI degree (number of partners above min_score).

    Parameters
    ----------
    min_score : minimum combined STRING score to count a partner (default 700 = high conf)
    limit     : number of top proteins to return (default 50)
    """
    def display_name(locus: str) -> str:
        try:
            n = run_server.normalize_gene_locus(locus)
            return n if n and n != locus else locus
        except Exception:
            return locus

    rows = []
    for gene, partners in STRING_INTERACTIONS.items():
        if gene == "_meta":
            continue
        filtered = [p for p in partners if p.get("score", 0) >= min_score]
        if not filtered:
            continue
        # Compute average per-channel scores
        def avg_ch(ch):
            vals = [p.get(ch, 0) for p in filtered]
            return round(sum(vals) / len(vals)) if vals else 0
        rows.append({
            "gene":          gene,
            "name":          display_name(gene),
            "degree":        len(filtered),
            "avg_score":     round(sum(p.get("score", 0) for p in filtered) / len(filtered)),
            "experimental":  avg_ch("experimental"),
            "database":      avg_ch("database"),
            "coexpression":  avg_ch("coexpression"),
            "textmining":    avg_ch("textmining"),
            "top_partners":  [p["partner"] for p in sorted(filtered, key=lambda x: x.get("score", 0), reverse=True)[:5]],
        })

    rows.sort(key=lambda x: x["degree"], reverse=True)
    return {
        "min_score":   min_score,
        "total_genes": len(rows),
        "hubs":        rows[:limit],
    }


@app.get("/api/analysis/network_ppi_edges")
def get_network_ppi_edges(genes: str = "", min_score: int = 400):
    """Return physical PPI edges existing among a set of genes (comma-separated locus tags).

    Used for cross-layer overlay mapping.
    """
    if not genes:
        return {"edges": []}
    
    gene_list = [g.strip().lower() for g in genes.split(",") if g.strip()]
    gene_set = set(gene_list)
    
    # Resolve aliases/synonyms
    resolved_set = set()
    for g in gene_set:
        resolved_set.add(g)
        try:
            for alias in run_server.expand_gene_aliases(g):
                resolved_set.add(alias.lower())
        except Exception:
            pass

    edges = []
    seen_pairs = set()
    for g in resolved_set:
        partners = STRING_INTERACTIONS.get(g, [])
        for p in partners:
            partner = p["partner"].lower()
            if partner in resolved_set:
                score = p.get("score", 0)
                if score >= min_score:
                    pair = (min(g, partner), max(g, partner))
                    if pair not in seen_pairs:
                        seen_pairs.add(pair)
                        edges.append({
                            "source": pair[0],
                            "target": pair[1],
                            "score": score,
                            "experimental": p.get("experimental", 0),
                            "database": p.get("database", 0),
                            "coexpression": p.get("coexpression", 0),
                            "neighborhood": p.get("neighborhood", 0),
                            "cooccurrence": p.get("cooccurrence", 0),
                            "textmining": p.get("textmining", 0),
                            "fusion": p.get("fusion", 0),
                        })
    return {"edges": edges}


@app.get("/api/analysis/cross_motifs")
def get_cross_motifs(motif_type: str = "co_complex", min_score: int = 400):
    """Detect and return regulatory-interaction cross-layer network motifs:

    1. 'co_complex' (Motif 1: TF_A regulates B and C, and protein B physically interacts with C)
    2. 'co_tf'      (Motif 2: TF_A physically interacts with TF_B, and both regulate C)
    3. 'feedback'   (Motif 3: TF_A regulates B, and TF_A physically interacts with B)
    """
    import csv

    # 1. Load TRN regulations
    path = run_server.get_absolute_path("data/reference/regulations.csv")
    if not os.path.exists(path):
        raise HTTPException(status_code=404, detail="regulations.csv not found")

    tf_targets = {}
    tf_names = {}
    tg_names = {}
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            tf = (row.get("TF_locusTag") or "").strip().lower()
            tf_name = (row.get("TF_name") or "").strip()
            tg = (row.get("TG_locusTag") or "").strip().lower()
            tg_name = (row.get("TG_name") or "").strip()
            if not tf or not tg:
                continue
            tf_targets.setdefault(tf, set()).add(tg)
            if tf_name and tf_name != tf:
                tf_names[tf] = tf_name
            if tg_name and tg_name != tg:
                tg_names[tg] = tg_name

    def display_name(locus: str) -> str:
        if locus in tf_names:
            return tf_names[locus]
        if locus in tg_names:
            return tg_names[locus]
        try:
            n = run_server.normalize_gene_locus(locus)
            return n if n and n != locus else locus
        except Exception:
            return locus

    # 2. Build fast PPI edge set
    ppi_edges = {} # (min(u,v), max(u,v)) -> score
    for gene, partners in STRING_INTERACTIONS.items():
        if gene == "_meta":
            continue
        g = gene.lower()
        for p in partners:
            partner = p["partner"].lower()
            score = p.get("score", 0)
            if score >= min_score:
                pair = (min(g, partner), max(g, partner))
                if pair not in ppi_edges or score > ppi_edges[pair]:
                    ppi_edges[pair] = score

    motifs = []
    
    # 3. Detect Motifs
    if motif_type == "co_complex":
        for tf, targets in tf_targets.items():
            targets_list = sorted(list(targets))
            n_targets = len(targets_list)
            for i in range(n_targets):
                b = targets_list[i]
                for j in range(i + 1, n_targets):
                    c = targets_list[j]
                    pair = (min(b, c), max(b, c))
                    if pair in ppi_edges:
                        motifs.append({
                            "tf": tf,
                            "tf_name": display_name(tf),
                            "target_b": b,
                            "target_b_name": display_name(b),
                            "target_c": c,
                            "target_c_name": display_name(c),
                            "ppi_score": ppi_edges[pair],
                        })
                        
    elif motif_type == "co_tf":
        tfs = sorted(list(tf_targets.keys()))
        n_tfs = len(tfs)
        for i in range(n_tfs):
            tf_a = tfs[i]
            for j in range(i + 1, n_tfs):
                tf_b = tfs[j]
                pair = (min(tf_a, tf_b), max(tf_a, tf_b))
                if pair in ppi_edges:
                    shared = tf_targets[tf_a].intersection(tf_targets[tf_b])
                    for c in sorted(list(shared)):
                        motifs.append({
                            "tf_a": tf_a,
                            "tf_a_name": display_name(tf_a),
                            "tf_b": tf_b,
                            "tf_b_name": display_name(tf_b),
                            "target_c": c,
                            "target_c_name": display_name(c),
                            "ppi_score": ppi_edges[pair],
                        })
                        
    elif motif_type == "feedback":
        for tf, targets in tf_targets.items():
            for b in sorted(list(targets)):
                pair = (min(tf, b), max(tf, b))
                if pair in ppi_edges:
                    motifs.append({
                        "tf": tf,
                        "tf_name": display_name(tf),
                        "target": b,
                        "target_name": display_name(b),
                        "ppi_score": ppi_edges[pair],
                    })

    # Sort results by interaction strength descending
    motifs.sort(key=lambda x: x.get("ppi_score", 0), reverse=True)
    return {
        "motif_type": motif_type,
        "min_score": min_score,
        "count": len(motifs),
        "instances": motifs[:200],
    }


@app.get("/api/quality/ppi")
def quality_ppi():
    try:
        genes = [k for k in STRING_INTERACTIONS.keys() if k != "_meta"]
        total_proteins = len(genes)
        
        edges = set()
        scores = []
        very_high = 0
        high = 0
        medium = 0
        low = 0
        
        channels = {
            "experimental": 0,
            "database": 0,
            "coexpression": 0,
            "textmining": 0,
            "neighborhood": 0,
            "cooccurrence": 0,
            "fusion": 0
        }
        
        for g in genes:
            for p in STRING_INTERACTIONS[g]:
                partner = p.get("partner")
                if not partner:
                    continue
                edge_key = tuple(sorted([g, partner.lower()]))
                if edge_key not in edges:
                    edges.add(edge_key)
                    score = p.get("score", 0)
                    scores.append(score)
                    
                    if score >= 900:
                        very_high += 1
                    elif score >= 700:
                        high += 1
                    elif score >= 400:
                        medium += 1
                    else:
                        low += 1
                        
                    for ch in channels.keys():
                        if p.get(ch, 0) > 0:
                            channels[ch] += 1
                            
        total_edges = len(edges)
        avg_score = sum(scores) / len(scores) if scores else 0
        
        return {
            "total_proteins": total_proteins,
            "total_interactions": total_edges,
            "avg_partners": round(total_edges * 2 / total_proteins, 1) if total_proteins > 0 else 0,
            "avg_score": round(avg_score, 1),
            "score_distribution": {
                "very_high": very_high,
                "high": high,
                "medium": medium,
                "low": low
            },
            "channel_support": channels
        }
    except Exception as e:
        logger.error(f"Error computing PPI quality: {e}")
        return {"error": str(e)}


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
if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    ROOT_DIR = sys._MEIPASS
else:
    ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

web_dir = os.path.join(ROOT_DIR, "web")
data_dir = os.path.join(ROOT_DIR, "data", "reference")

if os.path.exists(data_dir):
    app.mount("/data", StaticFiles(directory=data_dir), name="data")

if os.path.exists(web_dir):
    app.mount("/", StaticFiles(directory=web_dir, html=True), name="static")


