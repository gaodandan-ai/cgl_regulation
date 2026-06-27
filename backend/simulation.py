import logging
from typing import List, Dict, Any, Tuple, Optional

# Configure logging
logger = logging.getLogger("simulation")

def find_gene_in_model(model, gene_id: str) -> Any:
    """
    Search for a gene in the model by id, name, or label.
    Supports prefixes like 'G_' and case-insensitive matching.
    """
    normalized = gene_id.strip().lower()
    
    # Direct lookup
    if normalized in model.genes:
        return model.genes.get_by_id(normalized)
    
    # Try with 'G_' prefix
    g_prefixed = f"g_{normalized}"
    if g_prefixed in model.genes:
        return model.genes.get_by_id(g_prefixed)
        
    # Loose iteration matching
    for gene in model.genes:
        g_id = gene.id.lower()
        g_name = (gene.name or "").lower()
        stripped_id = g_id.replace("g_", "").replace("gene_", "")
        
        if (g_id == normalized or 
            g_name == normalized or 
            stripped_id == normalized):
            return gene
            
    return None

def apply_objective_to_model(model, objective_cfg: Any) -> Tuple[str, List[str]]:
    """
    Sets the objective function of the model based on objective_cfg.
    Returns the resolved label and any warnings.
    """
    warnings = []
    if not objective_cfg:
        return "Growth / biomass objective", warnings
        
    obj_type = getattr(objective_cfg, "objectiveType", "biomass")
    rxn_id = getattr(objective_cfg, "reactionId", None)
    
    if obj_type == "reaction" and rxn_id:
        if rxn_id in model.reactions:
            try:
                model.objective = model.reactions.get_by_id(rxn_id)
                return f"Reaction objective: {rxn_id}", warnings
            except Exception as e:
                warnings.append(f"Failed to set model objective to reaction '{rxn_id}': {str(e)}. Falling back to default biomass.")
        else:
            warnings.append(f"Selected objective reaction '{rxn_id}' not found in the model. Falling back to default biomass.")
            
    return "Growth / biomass objective", warnings

def run_fba_optimization(model) -> Tuple[float, str, List[str]]:
    """
    Optimizes the model and returns objective value, status, and warnings.
    """
    warnings = []
    try:
        solution = model.optimize()
        if solution.status == 'optimal':
            return float(solution.objective_value), solution.status, warnings
        else:
            warnings.append(f"Solver returned non-optimal status: '{solution.status}'. The model solution may be infeasible or unbounded.")
            return 0.0, solution.status, warnings
    except Exception as e:
        logger.error(f"FBA Optimization error: {str(e)}")
        warnings.append(f"FBA optimization failed: {str(e)}")
        return 0.0, "error", warnings

def run_baseline_simulation(model) -> Dict[str, Any]:
    """
    Run baseline FBA simulation.
    """
    obj_val, status, warnings = run_fba_optimization(model)
    obj_expr = str(model.objective.expression) if hasattr(model, 'objective') else "none"
    return {
        "status": status,
        "objective_value": obj_val,
        "objective_expression": obj_expr,
        "warnings": warnings
    }

def run_fba_simulation_pipeline(
    model, 
    knockout_genes: List[Any], 
    objective_cfg: Any, 
    track_reaction_ids: Optional[List[str]] = None,
    method: str = "fba"
) -> Tuple[str, Dict[str, Any], float, float, float, float, List[Dict[str, Any]], List[str]]:
    """
    Executes FBA or MOMA optimization pipeline on the model:
    - Runs baseline under customized objective
    - Runs perturbed under FBA or MOMA with knockouts
    - Tracks fluxes of selected reactions under both states.
    
    Returns:
        status, objective_response_dict, baseline_obj, perturbed_obj, change, change_percent, tracked_fluxes, warnings
    """
    warnings = []
    tracked_fluxes = []
    
    # 1. Run baseline to get objective label and baseline values
    baseline_solution = None
    with model:
        label, obj_warnings = apply_objective_to_model(model, objective_cfg)
        warnings.extend(obj_warnings)
        
        try:
            baseline_solution = model.optimize()
            b_status = baseline_solution.status
            baseline_obj = float(baseline_solution.objective_value) if b_status == 'optimal' else 0.0
            if b_status != 'optimal':
                warnings.append(f"Baseline solver returned non-optimal status: '{b_status}'.")
        except Exception as e:
            b_status = "error"
            baseline_obj = 0.0
            warnings.append(f"Baseline FBA optimization failed: {str(e)}")
            
        baseline_fluxes = {}
        if track_reaction_ids:
            for rxn_id in track_reaction_ids:
                if rxn_id in model.reactions:
                    try:
                        if b_status == 'optimal':
                            baseline_fluxes[rxn_id] = float(model.reactions.get_by_id(rxn_id).flux)
                        else:
                            baseline_fluxes[rxn_id] = 0.0
                    except Exception as e:
                        baseline_fluxes[rxn_id] = 0.0
                        warnings.append(f"Failed to get baseline flux for '{rxn_id}': {str(e)}")
                else:
                    warnings.append(f"Tracked reaction '{rxn_id}' not found in the model.")
                    
    # Construct objective details response
    obj_resp = {
        "objectiveType": objective_cfg.objectiveType if objective_cfg else "biomass",
        "reactionId": objective_cfg.reactionId if (objective_cfg and objective_cfg.objectiveType == "reaction") else None,
        "label": label
    }
    
    # 2. Run perturbed
    with model:
        apply_objective_to_model(model, objective_cfg)
        
        # Apply knockouts
        for gene in knockout_genes:
            try:
                model.genes.get_by_id(gene.id).knock_out()
            except Exception as e:
                warnings.append(f"Failed to knock out gene '{gene.id}': {str(e)}")
                
        p_status = "optimal"
        perturbed_obj = 0.0
        perturbed_fluxes = {}
        
        if method.lower() == "moma":
            # Run MOMA optimization
            from cobra.flux_analysis import moma
            sol_moma = None
            try:
                logger.info("Running MOMA (linear=True)...")
                sol_moma = moma(model, solution=baseline_solution, linear=True)
            except Exception as e:
                logger.info(f"Linear MOMA failed, trying standard MOMA: {str(e)}")
                try:
                    sol_moma = moma(model, solution=baseline_solution, linear=False)
                except Exception as ex:
                    p_status = "error"
                    warnings.append(f"MOMA optimization failed: {str(ex)}")
            
            if sol_moma and sol_moma.status == 'optimal':
                p_status = "optimal"
                
                # Retrieve actual biological objective reaction flux for perturbed growth rate
                biomass_reaction = None
                for rxn in model.reactions:
                    if rxn.objective_coefficient != 0:
                        biomass_reaction = rxn
                        break
                
                if biomass_reaction and biomass_reaction.id in sol_moma.fluxes:
                    perturbed_obj = float(sol_moma.fluxes[biomass_reaction.id])
                else:
                    perturbed_obj = float(sol_moma.objective_value)
                    
                # Store tracked fluxes
                if track_reaction_ids:
                    for rxn_id in track_reaction_ids:
                        if rxn_id in sol_moma.fluxes:
                            perturbed_fluxes[rxn_id] = float(sol_moma.fluxes[rxn_id])
            else:
                if p_status != "error":
                    p_status = sol_moma.status if sol_moma else "infeasible"
                    warnings.append(f"MOMA optimization returned non-optimal status: '{p_status}'")
                perturbed_obj = 0.0
        else:
            # Standard FBA optimization
            p_val, p_status, p_warnings = run_fba_optimization(model)
            warnings.extend(p_warnings)
            perturbed_obj = p_val
            
            if p_status == 'optimal' and track_reaction_ids:
                for rxn_id in track_reaction_ids:
                    if rxn_id in model.reactions:
                        try:
                            perturbed_fluxes[rxn_id] = float(model.reactions.get_by_id(rxn_id).flux)
                        except Exception as e:
                            warnings.append(f"Failed to get perturbed flux for '{rxn_id}': {str(e)}")
                            
        final_status = p_status if p_status == 'optimal' else b_status
        
        # Tracked fluxes
        if track_reaction_ids:
            for rxn_id in track_reaction_ids:
                if rxn_id in model.reactions:
                    b_flux = baseline_fluxes.get(rxn_id, 0.0)
                    p_flux = perturbed_fluxes.get(rxn_id, 0.0)
                    
                    flux_change = float(p_flux - b_flux)
                    flux_pct = 0.0
                    if abs(b_flux) > 1e-7:
                        flux_pct = float(flux_change / b_flux * 100)
                        
                    tracked_fluxes.append({
                        "reactionId": rxn_id,
                        "baselineFlux": b_flux,
                        "perturbedFlux": p_flux,
                        "fluxChange": flux_change,
                        "fluxChangePercent": flux_pct
                    })
                    
        # Calculate objective changes
        change = float(perturbed_obj - baseline_obj)
        change_pct = 0.0
        if abs(baseline_obj) > 1e-7:
            change_pct = float(change / baseline_obj * 100)
            
        return final_status, obj_resp, baseline_obj, perturbed_obj, change, change_pct, tracked_fluxes, warnings

def run_gene_knockout(
    model, 
    gene_id: str, 
    objective_cfg: Any = None, 
    track_reaction_ids: List[str] = None,
    method: str = "fba"
) -> Dict[str, Any]:
    """
    Simulates a single gene knockout.
    """
    gene = find_gene_in_model(model, gene_id)
    warnings = []
    
    if not gene:
        # Resolve fallback without knockout
        status, obj_resp, baseline_obj, perturbed_obj, change, change_pct, tracked_fluxes, pipeline_warnings = run_fba_simulation_pipeline(
            model, [], objective_cfg, track_reaction_ids, method
        )
        warnings.extend(pipeline_warnings)
        warnings.append(f"Gene '{gene_id}' not found in the metabolic model. Simulated as no metabolic change.")
        return {
            "status": status,
            "objective": obj_resp,
            "baselineObjective": baseline_obj,
            "perturbedObjective": baseline_obj,
            "objectiveChange": 0.0,
            "objectiveChangePercent": 0.0,
            "trackedFluxes": tracked_fluxes,
            "warnings": warnings
        }
        
    status, obj_resp, baseline_obj, perturbed_obj, change, change_pct, tracked_fluxes, pipeline_warnings = run_fba_simulation_pipeline(
        model, [gene], objective_cfg, track_reaction_ids, method
    )
    return {
        "status": status,
        "objective": obj_resp,
        "baselineObjective": baseline_obj,
        "perturbedObjective": perturbed_obj,
        "objectiveChange": change,
        "objectiveChangePercent": change_pct,
        "trackedFluxes": tracked_fluxes,
        "warnings": pipeline_warnings
    }

def run_gene_set_knockout(
    model, 
    gene_ids: List[str], 
    objective_cfg: Any = None, 
    track_reaction_ids: List[str] = None,
    method: str = "fba"
) -> Dict[str, Any]:
    """
    Simulates knockouts for multiple genes.
    """
    warnings = []
    missing_genes = []
    mapped_genes = []
    
    if not gene_ids:
        warnings.append("Empty gene list provided.")
        status, obj_resp, baseline_obj, perturbed_obj, change, change_pct, tracked_fluxes, pipeline_warnings = run_fba_simulation_pipeline(
            model, [], objective_cfg, track_reaction_ids, method
        )
        return {
            "status": status,
            "objective": obj_resp,
            "baselineObjective": baseline_obj,
            "perturbedObjective": baseline_obj,
            "objectiveChange": 0.0,
            "objectiveChangePercent": 0.0,
            "trackedFluxes": tracked_fluxes,
            "missingGenes": [],
            "warnings": warnings + pipeline_warnings
        }
        
    for gid in gene_ids:
        gene = find_gene_in_model(model, gid)
        if gene:
            mapped_genes.append(gene)
        else:
            missing_genes.append(gid)
            
    if missing_genes:
        warnings.append(f"{len(missing_genes)} genes were not found in the metabolic model: {', '.join(missing_genes)}")
        
    status, obj_resp, baseline_obj, perturbed_obj, change, change_pct, tracked_fluxes, pipeline_warnings = run_fba_simulation_pipeline(
        model, mapped_genes, objective_cfg, track_reaction_ids, method
    )
    warnings.extend(pipeline_warnings)
    
    return {
        "status": status,
        "objective": obj_resp,
        "baselineObjective": baseline_obj,
        "perturbedObjective": perturbed_obj,
        "objectiveChange": change,
        "objectiveChangePercent": change_pct,
        "trackedFluxes": tracked_fluxes,
        "missingGenes": missing_genes,
        "warnings": warnings
    }

def run_tf_perturbation(
    model, 
    tf_id: str, 
    target_gene_ids: List[str], 
    objective_cfg: Any = None, 
    track_reaction_ids: List[str] = None,
    method: str = "fba"
) -> Dict[str, Any]:
    """
    Simulates a transcription factor knockout by knocking out its downstream targets.
    """
    result = run_gene_set_knockout(model, target_gene_ids, objective_cfg, track_reaction_ids, method)
    
    return {
        "tfId": tf_id,
        "status": result["status"],
        "targetGeneCount": len(target_gene_ids),
        "mappedGeneCount": len(target_gene_ids) - len(result["missingGenes"]),
        "missingGenes": result["missingGenes"],
        "objective": result["objective"],
        "baselineObjective": result["baselineObjective"],
        "perturbedObjective": result["perturbedObjective"],
        "objectiveChange": result["objectiveChange"],
        "objectiveChangePercent": result["objectiveChangePercent"],
        "trackedFluxes": result["trackedFluxes"],
        "warnings": result["warnings"]
    }

def run_fva_analysis(
    model,
    knockout_genes: List[str],
    objective_cfg: Any,
    track_reaction_ids: Optional[List[str]] = None,
    fraction_of_optimum: float = 0.95
) -> Tuple[str, List[Dict[str, Any]], List[str]]:
    """
    Runs Flux Variability Analysis (FVA) under baseline and perturbed states.
    For each reaction in track_reaction_ids (and the active objective reaction):
    - Calculates maximum and minimum feasible fluxes under fraction_of_optimum.
    
    Returns:
        status, ranges_list, warnings
    """
    from cobra.flux_analysis import flux_variability_analysis
    warnings = []
    fva_ranges = []
    
    # 1. Resolve reactions list
    # If no track reactions specified, track the current objective reaction
    reactions_to_track = []
    if track_reaction_ids:
        for rid in track_reaction_ids:
            if rid in model.reactions:
                reactions_to_track.append(rid)
            else:
                warnings.append(f"Tracked reaction '{rid}' not found in metabolic model during FVA.")
                
    # Always include current objective reaction in FVA tracking
    with model:
        _, obj_warn = apply_objective_to_model(model, objective_cfg)
        warnings.extend(obj_warn)
        # Find active objective reaction ID(s)
        for rxn in model.reactions:
            if rxn.objective_coefficient != 0:
                if rxn.id not in reactions_to_track:
                    reactions_to_track.append(rxn.id)

    if not reactions_to_track:
        return "error", [], warnings + ["No valid reactions to track for FVA."]
        
    try:
        # 2. Run Baseline FVA
        baseline_min = {}
        baseline_max = {}
        
        with model:
            apply_objective_to_model(model, objective_cfg)
            logger.info(f"Running baseline FVA on {reactions_to_track} at fraction {fraction_of_optimum}...")
            df_baseline = flux_variability_analysis(
                model, 
                reaction_list=reactions_to_track, 
                fraction_of_optimum=fraction_of_optimum
            )
            for rid in reactions_to_track:
                if rid in df_baseline.index:
                    baseline_min[rid] = float(df_baseline.at[rid, "minimum"])
                    baseline_max[rid] = float(df_baseline.at[rid, "maximum"])
                else:
                    baseline_min[rid] = 0.0
                    baseline_max[rid] = 0.0
                    
        # 3. Run Perturbed FVA (with knockouts)
        perturbed_min = {}
        perturbed_max = {}
        
        with model:
            apply_objective_to_model(model, objective_cfg)
            
            # Apply gene knockouts
            mapped_ko_genes = []
            for g_id in knockout_genes:
                gene = find_gene_in_model(model, g_id)
                if gene:
                    mapped_ko_genes.append(gene)
                    gene.knock_out()
                    
            # Optimize to see if model is still feasible under perturbation
            sol = model.optimize()
            if sol.status != "optimal":
                warnings.append("Metabolic model became infeasible under perturbation. FVA perturbed ranges set to 0.")
                for rid in reactions_to_track:
                    perturbed_min[rid] = 0.0
                    perturbed_max[rid] = 0.0
            else:
                logger.info(f"Running perturbed FVA on {reactions_to_track} at fraction {fraction_of_optimum}...")
                df_perturbed = flux_variability_analysis(
                    model, 
                    reaction_list=reactions_to_track, 
                    fraction_of_optimum=fraction_of_optimum
                )
                for rid in reactions_to_track:
                    if rid in df_perturbed.index:
                        perturbed_min[rid] = float(df_perturbed.at[rid, "minimum"])
                        perturbed_max[rid] = float(df_perturbed.at[rid, "maximum"])
                    else:
                        perturbed_min[rid] = 0.0
                        perturbed_max[rid] = 0.0
                        
        # 4. Construct response list
        for rid in reactions_to_track:
            fva_ranges.append({
                "reactionId": rid,
                "baselineMin": baseline_min.get(rid, 0.0),
                "baselineMax": baseline_max.get(rid, 0.0),
                "perturbedMin": perturbed_min.get(rid, 0.0),
                "perturbedMax": perturbed_max.get(rid, 0.0)
            })
            
        return "optimal", fva_ranges, warnings
        
    except Exception as e:
        logger.error(f"FVA calculation error: {str(e)}")
        return "error", [], warnings + [f"FVA calculation failed: {str(e)}"]

