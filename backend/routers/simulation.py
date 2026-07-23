from fastapi import APIRouter, HTTPException, Response, Header
from fastapi.responses import JSONResponse
import os
import logging

from model_loader import get_model_status, load_model_if_needed
try:
    from thermo_pruner import get_pruning_report as _get_thermo_pruning_report
    _THERMO_PRUNER_AVAILABLE = True
except ImportError:
    _get_thermo_pruning_report = None
    _THERMO_PRUNER_AVAILABLE = False

from simulation import (
    run_baseline_simulation,
    run_gene_knockout,
    run_gene_set_knockout,
    run_tf_perturbation,
    run_fva_analysis,
    run_dynamic_rfba,
    run_dynamic_recfba,
    run_ecfba_simulation,
    run_mfa_comparison
)

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

from services.reference_data import (
    BRENDA_KCAT_MAPPINGS,
    RHEA_MAPPINGS,
    CHEBI_MAPPINGS,
    check_essentiality,
    check_abasy_role
)

try:
    import run_server
except ImportError:
    import backend.run_server as run_server

router = APIRouter(tags=["Simulation & Metabolic"])
logger = logging.getLogger("app.routers.simulation")

_THERMO_DATA_CACHE = None
_METABOLIC_IMPACT_CACHE: dict = {}


def _load_thermo_gene_data():
    """Load thermo data, cached after first call."""
    global _THERMO_DATA_CACHE
    if _THERMO_DATA_CACHE is not None:
        return _THERMO_DATA_CACHE
    backend_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    root_dir = os.path.dirname(backend_dir)
    path = os.path.join(root_dir, "data", "reference", "thermo_dgr_data.json")
    if not os.path.exists(path):
        return {}
    import json
    with open(path, "r", encoding="utf-8") as f:
        _THERMO_DATA_CACHE = json.load(f)
    return _THERMO_DATA_CACHE


def classify_glutamate_reaction(rxn):
    rxn_id_lower = rxn.id.lower()
    rxn_name_lower = rxn.name.lower()
    rxn_formula = rxn.reaction
    rxn_formula_lower = rxn_formula.lower()

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


@router.get("/api/model/status", response_model=ModelStatusResponse)
def model_status():
    return get_model_status()


@router.get("/api/model/reactions/glutamate-candidates", response_model=GlutamateCandidatesResponse)
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

            essential_genes_in_rxn = []
            is_rxn_essential = False

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

    has_export_or_exchange = any(c["classification"] in ("exchange", "export") for c in candidates)
    if not has_export_or_exchange:
        warnings.append("No high-confidence glutamate export or exchange reaction was identified in the loaded model. Please select a transport or uncertain candidate for tracking manually.")

    return {"candidates": candidates, "warnings": warnings}


@router.get("/api/model/reactions/search", response_model=ReactionSearchResponse)
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

        matches_rxn = (
            query_clean in rxn_id_lower or
            query_clean in rxn_name_lower or
            query_clean in rxn_formula_lower
        )

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
            if len(matches) >= 100:
                break

    return {"query": q, "matches": matches}


@router.post("/api/model/pathway/reactions")
def get_pathway_reactions(req: PathwayReactionsRequest):
    try:
        model = load_model_if_needed()
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

    results = []
    for rxn_id in req.reactionIds:
        lookup_id = rxn_id
        if rxn_id not in model.reactions and rxn_id.startswith("R_"):
            lookup_id = rxn_id[2:]
        if lookup_id in model.reactions:
            rxn = model.reactions.get_by_id(lookup_id)
            reactants = [m.id for m in rxn.reactants]
            products = [m.id for m in rxn.products]
            results.append({
                "reactionId": rxn_id,
                "name": rxn.name,
                "equation": rxn.reaction,
                "reactants": reactants,
                "products": products
            })
    return results


@router.post("/api/simulation/baseline", response_model=BaselineSimulationResponse)
def baseline_simulation():
    try:
        model = load_model_if_needed()
        result = run_baseline_simulation(model)
        return result
    except Exception as e:
        logger.error(f"Baseline FBA failed: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/api/simulation/gene-knockout", response_model=GeneKnockoutResponse)
def gene_knockout(req: GeneKnockoutRequest):
    try:
        model = load_model_if_needed()
        result = run_gene_knockout(model, req.geneId, req.objective, req.trackReactionIds, req.method)
        return result
    except Exception as e:
        logger.error(f"Gene knockout simulation failed: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/api/simulation/gene-set-knockout", response_model=GeneSetKnockoutResponse)
def gene_set_knockout(req: GeneSetKnockoutRequest):
    try:
        model = load_model_if_needed()
        result = run_gene_set_knockout(model, req.geneIds, req.objective, req.trackReactionIds, req.method)
        return result
    except Exception as e:
        logger.error(f"Gene set knockout simulation failed: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/api/simulation/tf-perturbation", response_model=TFPerturbationResponse)
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


@router.post("/api/simulation/flux-variability", response_model=FVAResponse)
def flux_variability(req: FVARequest):
    try:
        model = load_model_if_needed()

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


@router.post("/api/simulation/rfba", response_model=RFBAResponse)
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


@router.post("/api/simulation/recfba", response_model=RECFBAResponse)
def dynamic_recfba(req: RECFBARequest):
    try:
        backend_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        root_dir = os.path.dirname(backend_dir)
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


@router.post("/api/simulation/ecfba", response_model=ECFBAResponse)
def ecfba_simulation(req: ECFBARequest):
    try:
        backend_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        root_dir = os.path.dirname(backend_dir)
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


@router.get("/api/simulation/mfa_comparison", response_model=MFAComparisonResponse)
def get_mfa_comparison():
    try:
        model = load_model_if_needed()
        result = run_mfa_comparison(model)
        return result
    except Exception as e:
        logger.error(f"MFA comparison failed: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/thermo/gene_context")
def get_thermo_gene_context(gene: str = ""):
    if not gene:
        raise HTTPException(status_code=400, detail="Missing gene parameter")

    thermo_data = _load_thermo_gene_data()
    thermo_rxns = thermo_data.get("reactions", {})

    model = load_model_if_needed()
    if model is None:
        raise HTTPException(status_code=503, detail="Model not loaded")

    gene_lower = gene.strip().lower()
    aliases = set()
    try:
        aliases = run_server.expand_gene_aliases(gene_lower)
    except Exception:
        aliases = {gene_lower, gene.strip()}

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


@router.get("/api/metabolic_impact")
def metabolic_impact(gene: str = "", query: str = ""):
    target = (gene or query).strip().lower()
    if not target:
        raise HTTPException(status_code=400, detail="Missing gene parameter")

    if target in _METABOLIC_IMPACT_CACHE:
        return JSONResponse(
            content=_METABOLIC_IMPACT_CACHE[target],
            headers={"Cache-Control": "public, max-age=3600", "X-Cache": "HIT"}
        )

    try:
        result = run_server.handle_metabolic_impact(target)
        _METABOLIC_IMPACT_CACHE[target] = result
        return JSONResponse(
            content=result,
            headers={"Cache-Control": "public, max-age=3600", "X-Cache": "MISS"}
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/metabolic_pathways")
def metabolic_pathways(response: Response, pathway: str = "", query: str = ""):
    target = pathway or query
    try:
        result = run_server.handle_metabolic_pathways(target)
        if not target:
            response.headers["Cache-Control"] = "public, max-age=3600"
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/imodulon/reactions")
def imodulon_reactions(imodulon: str = ""):
    try:
        result = run_server.handle_imodulon_reactions(imodulon)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/imodulon/simulation")
def imodulon_simulation(imodulon: str = ""):
    try:
        result = run_server.handle_imodulon_simulation(imodulon)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/engineering/simulation")
def engineering_simulation(tf: str = ""):
    try:
        result = run_server.handle_tf_simulation(tf)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/thermo/pruning-report")
def get_thermo_pruning_report(response: Response):
    response.headers["Cache-Control"] = "public, max-age=1800"
    if not _THERMO_PRUNER_AVAILABLE or _get_thermo_pruning_report is None:
        return {
            "status": "unavailable",
            "message": "Thermo directionality pruner is not available."
        }
    try:
        load_model_if_needed()
        report = _get_thermo_pruning_report()
        return report
    except Exception as e:
        logger.error(f"Failed to generate thermo pruning report: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))
