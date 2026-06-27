from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
import logging
from model_loader import get_model_status, load_model_if_needed
from simulation import run_baseline_simulation, run_gene_knockout, run_gene_set_knockout, run_tf_perturbation, run_fva_analysis
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
    FVAResponse
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

