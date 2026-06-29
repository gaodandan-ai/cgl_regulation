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

def run_dynamic_rfba(
    model,
    tf_perturbations: Dict[str, str],  # tf_locus_tag -> "knockout" / "overexpress" / "normal"
    initial_glucose: float = 100.0,
    initial_biomass: float = 0.1,
    time_steps: int = 24,
) -> Dict[str, Any]:
    import os
    import json
    import csv
    import numpy as np

    warnings = []
    
    # 1. Resolve paths
    ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    DATA_DIR = os.path.join(ROOT_DIR, "data", "reference")
    RESULTS_PATH = os.path.join(DATA_DIR, "rna_seq_analysis_results.json")
    REGULATIONS_PATH = os.path.join(DATA_DIR, "regulations.csv")
    MAPPING_PATH = os.path.join(DATA_DIR, "gene_mapping.csv")

    # 2. Build mapping and network dictionaries
    name_to_cg = {}
    cg_to_name = {}
    cg_to_cgl = {}
    cgl_to_cg = {}

    if os.path.exists(MAPPING_PATH):
        try:
            with open(MAPPING_PATH, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    cg = row.get("cg_locus", "").strip()
                    cgl = row.get("cgl_locus", "").strip()
                    name = row.get("gene_name", "").strip()
                    if cg:
                        if name:
                            name_to_cg[name.lower()] = cg
                            name_to_cg[name.upper()] = cg
                            cg_to_name[cg] = name
                        if cgl:
                            cg_to_cgl[cg] = cgl
                            cgl_to_cg[cgl] = cg
        except Exception as e:
            warnings.append(f"Failed to load gene mapping: {str(e)}")

    tf_to_tg = {}  # TF cg_locus -> list of (tg_cgl, role)
    if os.path.exists(REGULATIONS_PATH):
        try:
            with open(REGULATIONS_PATH, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    tf_cg = row.get("TF_locusTag", "").strip()
                    tg_cg = row.get("TG_locusTag", "").strip()
                    role = row.get("Role", "").strip()
                    if tf_cg and tg_cg:
                        tg_cgl = cg_to_cgl.get(tg_cg, tg_cg)
                        tf_to_tg.setdefault(tf_cg, []).append((tg_cgl, role))
        except Exception as e:
            warnings.append(f"Failed to load regulations CSV: {str(e)}")

    # 3. Load baseline trajectories
    trajectories = {}
    if os.path.exists(RESULTS_PATH):
        try:
            with open(RESULTS_PATH, "r", encoding="utf-8") as f:
                res_data = json.load(f)
            trajectories = res_data.get("dynamic_grn", {}).get("trajectories", {})
        except Exception as e:
            warnings.append(f"Failed to load baseline trajectories: {str(e)}")

    # 4. Resolve users input perturbations (map from TF names to cg locus tags)
    resolved_perturbations = {}
    for tf_input, mode in tf_perturbations.items():
        tf_input_clean = tf_input.strip()
        tf_cg = tf_input_clean
        if tf_input_clean.lower() in name_to_cg:
            tf_cg = name_to_cg[tf_input_clean.lower()]
        elif tf_input_clean.upper() in name_to_cg:
            tf_cg = name_to_cg[tf_input_clean.upper()]
        
        resolved_perturbations[tf_cg] = mode

    # 5. Compute target gene expression multipliers based on TF perturbations
    tg_multipliers = {}  # target_cgl_locus -> float multiplier
    for tf_cg, mode in resolved_perturbations.items():
        if mode == "normal":
            continue
            
        targets = tf_to_tg.get(tf_cg, [])
        for tg_cgl, role in targets:
            mult = tg_multipliers.get(tg_cgl, 1.0)
            
            if role == "A" or "activat" in role.lower():
                if mode == "knockout":
                    mult *= 0.1
                elif mode == "overexpress":
                    mult *= 2.5
            elif role == "R" or "repress" in role.lower():
                if mode == "knockout":
                    mult *= 2.0
                elif mode == "overexpress":
                    mult *= 0.05
            
            tg_multipliers[tg_cgl] = mult

    # 6. Calculate temporal gene expression trajectories (0h to 24h)
    gene_ts = {}
    metabolic_genes = [g.id.replace("g_", "").replace("gene_", "") for g in model.genes]
    
    for g_cg in metabolic_genes:
        g_cgl = cg_to_cgl.get(g_cg, g_cg)
        traj = trajectories.get(g_cgl)
        multiplier = tg_multipliers.get(g_cgl, 1.0)
        
        if traj:
            try:
                traj_25 = np.array(traj)
                base_val = traj_25[0] if traj_25[0] > 1e-5 else 1.0
                gene_ts[g_cg] = (traj_25 / base_val) * multiplier
            except Exception:
                gene_ts[g_cg] = np.ones(25) * multiplier
        else:
            gene_ts[g_cg] = np.ones(25) * multiplier

    # 7. Find key reactions in the model
    biomass_rxn = None
    for rxn in model.reactions:
        if rxn.objective_coefficient != 0:
            biomass_rxn = rxn
            break
            
    glu_ex_rxn = None
    for rxn in model.reactions:
        if rxn.id == "EX_glu_L_e":
            glu_ex_rxn = rxn
            break
            
    if not glu_ex_rxn:
        for rxn in model.reactions:
            if rxn.id.lower().startswith("ex_glu") or ("exchange" in rxn.name.lower() and "glutamate" in rxn.name.lower()):
                glu_ex_rxn = rxn
                break

    glc_ex_rxn = None
    for rxn in model.reactions:
        if rxn.id == "EX_glc_e":
            glc_ex_rxn = rxn
            break

    # 8. Setup dynamic state and histories
    glucose = initial_glucose
    biomass = initial_biomass
    
    time_history = []
    growth_rate_history = []
    glutamate_export_history = []
    glucose_uptake_history = []
    glucose_history = [glucose]
    biomass_history = [biomass]
    
    # Store original bounds
    orig_bounds = {r.id: (r.lower_bound, r.upper_bound) for r in model.reactions}
    
    V_max_glc = 10.0  # mmol/gDW/h
    K_m_glc = 1.0  # mM
    volume = 1.0  # L
    dt = 1.0  # hour
    
    for t in range(time_steps):
        time_history.append(float(t))
        
        # Calculate glucose maximum uptake rate based on Michaelis-Menten
        if glucose <= 1e-5:
            v_uptake_max = 0.0
        else:
            S = glucose / volume
            v_uptake_max = V_max_glc * (S / (K_m_glc + S))
            
        with model:
            if glc_ex_rxn:
                glc_ex_rxn.lower_bound = - v_uptake_max
                glc_ex_rxn.upper_bound = 1000.0
                
            for rxn in model.reactions:
                if rxn.id.startswith("EX_"):
                    continue
                    
                rxn_genes = [g.id.replace("g_", "").replace("gene_", "") for g in rxn.genes]
                if not rxn_genes:
                    continue
                    
                ratios = [gene_ts[g][min(t, 24)] for g in rxn_genes if g in gene_ts]
                if ratios:
                    ratio = np.mean(ratios)
                    ratio = np.clip(ratio, 0.0, 10.0)
                    
                    lower_orig, upper_orig = orig_bounds[rxn.id]
                    if upper_orig > 0:
                        rxn.upper_bound = upper_orig * ratio
                    if lower_orig < 0:
                        rxn.lower_bound = lower_orig * ratio
            
            sol = model.optimize()
            
            if sol.status == "optimal":
                mu = float(sol.objective_value) if biomass_rxn else 0.0
                v_glu = float(sol.fluxes[glu_ex_rxn.id]) if glu_ex_rxn else 0.0
                v_glc = - float(sol.fluxes[glc_ex_rxn.id]) if (glc_ex_rxn and sol.fluxes[glc_ex_rxn.id] < 0) else 0.0
            else:
                mu = 0.0
                v_glu = 0.0
                v_glc = 0.0
                
        delta_biomass = mu * biomass * dt
        biomass_next = biomass + delta_biomass
        X_avg = 0.5 * (biomass + biomass_next)
        
        delta_glucose = - v_glc * X_avg * dt
        glucose_next = max(0.0, glucose + delta_glucose)
        
        biomass = biomass_next
        glucose = glucose_next
        
        growth_rate_history.append(mu)
        glutamate_export_history.append(v_glu)
        glucose_uptake_history.append(v_glc)
        glucose_history.append(glucose)
        biomass_history.append(biomass)

    time_history.append(float(time_steps))
    
    return {
        "status": "success",
        "time": time_history,
        "growth_rate": growth_rate_history,
        "glutamate_export": glutamate_export_history,
        "glucose_uptake": glucose_uptake_history,
        "glucose_concentration": glucose_history,
        "biomass_concentration": biomass_history,
        "warnings": warnings
    }


def run_ecfba_simulation(
    json_model_path: str,
    protein_pool_limit: float,
    enzyme_perturbations: Dict[str, float],
    target_product: str,
    temperature: float = 30.0,
    calibrate_timepoint: Optional[str] = None
) -> Dict[str, Any]:
    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from enzyme_thermal_params import get_params, compute_alpha
    import json
    import csv
    import cobra
    import math

    warnings = []
    
    if not os.path.exists(json_model_path):
        return {
            "status": "error",
            "flux": 0.0,
            "pool_limit": protein_pool_limit,
            "pool_usage": 0.0,
            "warnings": [f"Enzyme constrained JSON model not found at path: {json_model_path}"]
        }
        
    try:
        # Load and parse experimental calibration timepoint data if specified
        # NOTE: All timepoints (1h, 4h, 24h) were sampled at constant 40°C heat stress.
        # They represent the temporal adaptation response (acute → early → chronic) at 40°C.
        # Source: RNA-Seq + proteomics (limma differential analysis vs 30°C control)
        calibrated_perturbations = None
        if calibrate_timepoint in ["1h", "4h", "24h"]:
            try:
                root_dir = os.path.dirname(os.path.dirname(__file__))

                # --- Priority 1: Proteomics (direct protein abundance) ---
                prot_path = os.path.join(root_dir, "data", "raw", "protemic",
                                         f"{calibrate_timepoint}_all_results.csv")
                # --- Priority 2: RNA-Seq (transcriptomic proxy) ---
                rna_path = os.path.join(root_dir, "data", "raw", "rna_seq",
                                        f"01_{calibrate_timepoint}_differential_result.csv")

                # Gene-to-enzyme mapping for ec-FBA perturbations
                # Keys: CGL gene locus → (perturbation_name, logFC_column, padj_column)
                PROT_GENE_MAP = {
                    "Cgl2079": "gdh",    # GDH (glutamate dehydrogenase)
                    "Cgl0851": "pgi",    # PGI (phosphoglucose isomerase)
                    "Cgl2089": "pyk",    # PYK1 (pyruvate kinase)
                    "Cgl2380": "mdh",    # MDH (malate dehydrogenase)
                    "Cgl0664": "icdh",   # ICDH (isocitrate dehydrogenase)
                    "Cgl0937": "gapdh",  # GAPDH
                }
                # LysC locus not found in ecCGL1 kcat reactions; keep RNA fallback
                RNA_GENE_MAP = {
                    "Cgl0251": "lysC",   # LysC (aspartokinase)
                }

                fc_map = {}
                data_source = "default"

                # Try proteomics first
                if os.path.exists(prot_path):
                    with open(prot_path, "r", encoding="utf-8") as csvfile:
                        reader = csv.DictReader(csvfile)
                        for row in reader:
                            gene_id = row.get("cgl_id", "").strip()
                            if gene_id in PROT_GENE_MAP:
                                try:
                                    logfc = float(row.get("logFC", 0.0))
                                    padj  = float(row.get("adj.P.Val", 1.0))
                                    name  = PROT_GENE_MAP[gene_id]
                                    # Use proteomics FC; if not significant, still use
                                    # the measured value (not zero) — biological signal
                                    fc_map[name] = 2.0 ** logfc
                                except (ValueError, TypeError):
                                    pass
                    data_source = "proteomics"

                # Always supplement LysC from RNA-Seq (not in proteomics map)
                if os.path.exists(rna_path):
                    with open(rna_path, "r", encoding="utf-8") as csvfile:
                        reader = csv.DictReader(csvfile)
                        for row in reader:
                            gene_id = row.get("Geneid", "").strip()
                            if gene_id in RNA_GENE_MAP:
                                try:
                                    log2fc = float(row.get("log2FC", 0.0))
                                    name   = RNA_GENE_MAP[gene_id]
                                    if name not in fc_map:  # don't overwrite proteomics
                                        fc_map[name] = 2.0 ** log2fc
                                except (ValueError, TypeError):
                                    pass
                    if data_source == "default":
                        # Proteomics not found; fall back to RNA for GDH too
                        with open(rna_path, "r", encoding="utf-8") as csvfile:
                            reader = csv.DictReader(csvfile)
                            for row in reader:
                                gene_id = row.get("Geneid", "").strip()
                                if gene_id == "Cgl2079":
                                    try:
                                        fc_map["gdh"] = 2.0 ** float(row.get("log2FC", 0.0))
                                    except (ValueError, TypeError):
                                        pass
                        data_source = "rna_seq"

                # Apply bounds and build perturbations dict
                enzyme_perturbations = {
                    name: max(0.05, min(3.0, fc))
                    for name, fc in fc_map.items()
                }
                calibrated_perturbations = enzyme_perturbations

                # Build warning message summarizing calibration
                summary_parts = [f"{n}={fc:.2f}x" for n, fc in enzyme_perturbations.items()]
                warnings.append(
                    f"Simulation calibrated with {calibrate_timepoint} {data_source} data. "
                    f"Enzyme adjustments: {', '.join(summary_parts)}."
                )

            except Exception as e:
                warnings.append(f"Failed to read calibration data: {str(e)}. Using manual inputs.")
        # Calculate temperature-dependent HSP proteome allocation cost (Tug-of-war)
        try:
            if temperature > 30.0:
                t_diff = max(0.0, temperature - 30.0)
                # HSP pool cost (proteomics-calibrated Hill function)
                # GroEL/GroES fold-changes from proteomics: 1h(40C)=6.5x, 4h(39C)=122x, 24h(37C)=506x
                # Hill fit: k=11.0 (half-saturation at +11C above 30C = 41C), n=3
                # h_max=0.20: at saturation, HSPs occupy max 20% of total proteome
                # Reference: proteomics data (this study); Takeno et al. 2010 AEM
                hsp_fraction = 0.20 * (t_diff**3.0) / (11.0**3.0 + t_diff**3.0)
                adjusted_pool_limit = protein_pool_limit * (1.0 - hsp_fraction)
                warnings.append(
                    f"Heat stress at {temperature:.1f}°C: SigH/SigE-regulated HSPs "
                    f"(GroEL/GroES calibrated from proteomics) occupy {hsp_fraction*100:.1f}% "
                    f"of proteome. Metabolic pool: {protein_pool_limit:.3f} → {adjusted_pool_limit:.3f} g/gDW."
                )
            else:
                adjusted_pool_limit = protein_pool_limit
        except Exception as e:
            adjusted_pool_limit = protein_pool_limit
            warnings.append(f"HSP allocation calculation error: {str(e)}.")

        with open(json_model_path, "r", encoding="utf-8") as f:
            dictionary_model = json.load(f)
            
        model = cobra.io.load_json_model(json_model_path)

        # Walk up from model path to find data/gene_mapping.csv
        curr = os.path.dirname(os.path.abspath(json_model_path))
        mapping_path = None
        for _ in range(10):
            test_path = os.path.join(curr, "data", "reference", "gene_mapping.csv")
            if os.path.exists(test_path):
                mapping_path = test_path
                break
            curr = os.path.dirname(curr)
            
        name_to_cgl = {}
        cg_to_cgl = {}
        
        if mapping_path and os.path.exists(mapping_path):
            try:
                with open(mapping_path, "r", encoding="utf-8") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        cg = row.get("cg_locus", "").strip()
                        cgl = row.get("cgl_locus", "").strip()
                        name = row.get("gene_name", "").strip()
                        if cgl:
                            if name:
                                name_to_cgl[name.lower()] = cgl
                                name_to_cgl[name.upper()] = cgl
                            if cg:
                                cg_to_cgl[cg.lower()] = cgl
                                cg_to_cgl[cg.upper()] = cgl
            except Exception as e:
                warnings.append(f"Failed to load gene mapping for ecFBA: {str(e)}")
        else:
            warnings.append("gene_mapping.csv not found, name-based resolution disabled.")

        # 2. Build coefficients and constraint
        coefficients = dict()
        rxn_to_kcatmw = {}
        for rxn in model.reactions:
            for eachr in dictionary_model['reactions']:
                if rxn.id == eachr['id']:
                    kcat_mw = eachr.get('kcat_MW')
                    if kcat_mw:
                        rxn_genes = [g.id.replace("g_", "").replace("gene_", "") for g in rxn.genes]

                        # --- Per-enzyme literature-calibrated thermal parameters ---
                        # Lookup via gene locus -> enzyme_thermal_params.GENE_LOCUS_PARAMS
                        # then reaction ID -> REACTION_ID_PARAMS, then global default.
                        # All E_a/H_d/S_d values have published citations.
                        thermo_p = get_params(rxn.id, rxn_genes)
                        E_a = thermo_p["E_a"]
                        H_d = thermo_p["H_d"]
                        S_d = thermo_p["S_d"]

                        # Calculate reaction-specific alpha(T)
                        T_kelvin = temperature + 273.15
                        try:
                            rxn_alpha = compute_alpha(thermo_p, T_kelvin)
                        except Exception:
                            rxn_alpha = 1.0
                            
                        # Direction B: Thermodynamic Flux-Force Efficiency factor
                        # dG_ref and dH_rxn for near-equilibrium reactions
                        # Gene loci corrected to match ecCGL1 model (verified by inspect_ecfba_genes.py)
                        T_ref = 30.0 + 273.15
                        R = 8.314
                        dG_ref = 3000.0
                        dH_rxn = -10000.0

                        is_pgi = "Cgl0851" in rxn_genes or "PGI" in rxn.id  # corrected: Cgl0851
                        is_mdh = "Cgl2380" in rxn_genes or ("MDH" in rxn.id and "UAMDH" not in rxn.id)  # corrected: Cgl2380
                        is_icdh = "Cgl0664" in rxn_genes or "ICDHyr" in rxn.id  # corrected: Cgl0664 (was Cgl1504 wrongly mapped)

                        if is_pgi:
                            dG_ref = 1500.0
                            dH_rxn = -2500.0
                        elif is_icdh:
                            dG_ref = 1000.0
                            dH_rxn = 8400.0
                        elif is_mdh:
                            dG_ref = 800.0
                            dH_rxn = 29000.0
                            
                        try:
                            # -dG(T) = -dG_ref - dH_rxn * (1.0 - T/T_ref)
                            dg_t = dG_ref + dH_rxn * (1.0 - T_kelvin / T_ref)
                            dg_t = max(100.0, dg_t)
                            
                            eta_thermo = math.tanh(dg_t / (2.0 * R * T_kelvin))
                            eta_thermo = max(1e-4, eta_thermo)
                        except Exception:
                            eta_thermo = 0.5
                            
                        coefficients[rxn.forward_variable] = 1 / (float(kcat_mw) * rxn_alpha * eta_thermo)
                        rxn_to_kcatmw[rxn.id] = float(kcat_mw) * rxn_alpha * eta_thermo
                    break

        constraint = model.problem.Constraint(0, lb=0.0, ub=adjusted_pool_limit, name="enzyme_pool_limit")
        model.add_cons_vars(constraint)
        model.solver.update()
        constraint.set_linear_coefficients(coefficients=coefficients)

        orig_bounds = {r.id: (r.lower_bound, r.upper_bound) for r in model.reactions}

        # 3. Apply enzyme level perturbations
        for gene_locus, level in enzyme_perturbations.items():
            normalized = gene_locus.strip().lower()
            
            resolved_cgl = normalized
            if normalized in name_to_cgl:
                resolved_cgl = name_to_cgl[normalized]
            elif normalized.upper() in name_to_cgl:
                resolved_cgl = name_to_cgl[normalized.upper()]
            elif normalized in cg_to_cgl:
                resolved_cgl = cg_to_cgl[normalized]
            elif normalized.upper() in cg_to_cgl:
                resolved_cgl = cg_to_cgl[normalized.upper()]
                
            gene = None
            if resolved_cgl in model.genes:
                gene = model.genes.get_by_id(resolved_cgl)
            else:
                g_prefixed = f"g_{resolved_cgl}"
                if g_prefixed in model.genes:
                    gene = model.genes.get_by_id(g_prefixed)
            
            if gene:
                for rxn in gene.reactions:
                    if rxn.id.startswith("EX_"):
                        continue
                    lower, upper = orig_bounds[rxn.id]
                    rxn.lower_bound = lower * level
                    rxn.upper_bound = upper * level
            else:
                warnings.append(f"Gene '{gene_locus}' (resolved as '{resolved_cgl}') not found in ec-FBA model.")

        # 4. Set objective reaction
        biomass_id = "CG_biomass_cgl_ATCC13032"
        glutamate_id = "EX_glu_L_e"
        lysine_id = "EX_lys_L_e"
        
        if target_product == "growth":
            if biomass_id in model.reactions:
                model.objective = model.reactions.get_by_id(biomass_id)
            else:
                warnings.append("Growth biomass objective function not found in ecCGL1.")
        elif target_product == "glutamate":
            if glutamate_id in model.reactions:
                model.objective = model.reactions.get_by_id(glutamate_id)
            else:
                warnings.append("EX_glu_L_e exchange reaction not found in ecCGL1.")
        elif target_product == "lysine":
            if lysine_id in model.reactions:
                model.objective = model.reactions.get_by_id(lysine_id)
            else:
                warnings.append("EX_lys_L_e exchange reaction not found in ecCGL1.")

        # 5. Optimize using pFBA (parsimonious FBA)
        # pFBA simultaneously maximizes the objective AND minimizes total flux.
        # This eliminates thermodynamically infeasible high-flux loops and
        # produces flux distributions much closer to 13C-MFA observations.
        # Reference: Lewis et al. 2010, Mol. Syst. Biol. (original pFBA paper)
        try:
            from cobra.flux_analysis import pfba, flux_variability_analysis
            sol = pfba(model, fraction_of_optimum=1.0)
            sol_status = "optimal"
        except Exception as pfba_err:
            warnings.append(f"pFBA failed ({pfba_err}), falling back to standard FBA.")
            sol = model.optimize()
            sol_status = sol.status

        if sol_status == "optimal" or (hasattr(sol, 'status') and sol.status == "optimal"):
            pool_usage = 0.0
            for rxn_id, kcatmw in rxn_to_kcatmw.items():
                if rxn_id in sol.fluxes:
                    pool_usage += sol.fluxes[rxn_id] * (1 / kcatmw)

            # Extract biomass growth rate from the flux vector (not objective_value,
            # which in pFBA represents minimized total flux, not the growth rate).
            biomass_flux = 0.0
            for bio_id in ["CG_biomass_cgl_ATCC13032", "BIOMASS_Cgl_ATCC13032", "Growth"]:
                if bio_id in sol.fluxes:
                    biomass_flux = float(sol.fluxes[bio_id])
                    break
            # Fallback: use objective_value only if not pFBA
            if biomass_flux == 0.0 and not (sol_status == "optimal" and hasattr(sol, 'fluxes')):
                biomass_flux = float(sol.objective_value)

            # Extract top bottlenecks
            bottlenecks = []
            try:
                pool_dual = 0.0
                if hasattr(model, 'constraints') and "enzyme_pool_limit" in model.constraints:
                    pool_dual = float(model.constraints["enzyme_pool_limit"].dual)
            except Exception:
                pool_dual = 0.0

            for rxn_id, kcatmw in rxn_to_kcatmw.items():
                if rxn_id in sol.fluxes:
                    flux = sol.fluxes[rxn_id]
                    if abs(flux) > 1e-5:
                        usage = abs(flux) / kcatmw
                        try:
                            rxn_obj = model.reactions.get_by_id(rxn_id)
                            rxn_name = rxn_obj.name
                            genes = [g.id.replace("g_", "") for g in rxn_obj.genes]
                            gene_names_str = ", ".join(genes)
                        except Exception:
                            rxn_name = "Unknown reaction"
                            gene_names_str = ""
                        
                        shadow_price = pool_dual / kcatmw if kcatmw > 0 else 0.0
                        bottlenecks.append({
                            "reaction_id": rxn_id,
                            "reaction_name": rxn_name,
                            "genes": gene_names_str,
                            "flux": float(flux),
                            "usage": float(usage),
                            "shadow_price": float(shadow_price)
                        })
            
            # Sort bottlenecks by usage descending
            bottlenecks.sort(key=lambda x: x["usage"], reverse=True)
            top_bottlenecks = bottlenecks[:8]

            return {
                "status": "success",
                "flux": biomass_flux,
                "pool_limit": adjusted_pool_limit,
                "pool_usage": float(pool_usage),
                "warnings": warnings,
                "calibratedPerturbations": calibrated_perturbations,
                "bottlenecks": top_bottlenecks
            }
        else:
            return {
                "status": "infeasible",
                "flux": 0.0,
                "pool_limit": adjusted_pool_limit,
                "pool_usage": 0.0,
                "warnings": warnings + ["Model is infeasible under current constraints."],
                "calibratedPerturbations": calibrated_perturbations,
                "bottlenecks": []
            }
            
    except Exception as e:
        logger.error(f"ec-FBA simulation error: {str(e)}")
        return {
            "status": "error",
            "flux": 0.0,
            "pool_limit": adjusted_pool_limit,
            "pool_usage": 0.0,
            "warnings": warnings + [f"ec-FBA calculation failed: {str(e)}"],
            "calibratedPerturbations": None
        }


def run_mfa_comparison(json_model_path: str) -> Dict[str, Any]:
    """
    Compare FBA simulated fluxes against published 13C-MFA literature data
    for C. glutamicum ATCC 13032 wild-type (aerobic glucose minimal medium).
    Reference: Cheng et al. 2017, Becker & Wittmann 2011.
    """
    import os
    import json
    import cobra
    import math
    import sys

    sys.path.insert(0, os.path.dirname(__file__))
    from mfa_reference import MFA_LITERATURE_DATASET, REACTION_ID_ALIASES

    warnings_list = []

    if not os.path.exists(json_model_path):
        return {
            "status": "error",
            "items": [],
            "pearson_r": 0.0,
            "rmse": 0.0,
            "mean_deviation_pct": 0.0,
            "warnings": [f"Model file not found: {json_model_path}"]
        }

    try:
        model = cobra.io.load_json_model(json_model_path)

        # Freeze all secretion-positive exchanges but leave essential uptakes open
        # Strategy: only close exchanges that are already closed in model defaults,
        # then explicitly open key nutrients and set glucose to 10 mmol/gDW/h
        with model:
            # Allow all secretion freely (positive direction)
            for rxn in model.exchanges:
                if rxn.lower_bound >= 0:
                    pass  # already not consuming
                # Force product secretions to allow free export
                if rxn.upper_bound < 0:
                    rxn.upper_bound = 1000.0

            # Open minimal medium nutrients with generous limits
            for nut_id in ["EX_o2_e", "EX_nh4_e", "EX_pi_e", "EX_so4_e",
                            "EX_mg2_e", "EX_k_e", "EX_h2o_e", "EX_h_e", "EX_fe2_e", "EX_fe3_e",
                            "EX_mn2_e", "EX_zn2_e", "EX_cobalt2_e", "EX_cu2_e"]:
                if nut_id in model.reactions:
                    model.reactions.get_by_id(nut_id).lower_bound = -1000.0

            # Set glucose uptake rate = 10 mmol/gDW/h
            glc_set = False
            for glc_id in ["EX_glc_e", "EX_glc__D_e", "EX_glc_D_e"]:
                if glc_id in model.reactions:
                    model.reactions.get_by_id(glc_id).lower_bound = -10.0
                    glc_set = True
                    break
            if not glc_set:
                for rxn in model.exchanges:
                    if "glc" in rxn.id.lower() and "g1p" not in rxn.id.lower() and "g6p" not in rxn.id.lower():
                        rxn.lower_bound = -10.0
                        glc_set = True
                        break

            # Set growth objective (maximize biomass)
            for bid in ["CG_biomass_cgl_ATCC13032", "BIOMASS_Cgl_ATCC13032", "Growth"]:
                if bid in model.reactions:
                    model.objective = model.reactions.get_by_id(bid)
                    break

            # Use pFBA for parsimonious (biologically realistic) flux distribution
            # pFBA minimizes total flux while maintaining maximum growth,
            # giving flux distributions much closer to 13C-MFA observations.
            try:
                from cobra.flux_analysis import pfba, flux_variability_analysis
                sol = pfba(model, fraction_of_optimum=1.0)
                used_pfba = True
            except Exception:
                sol = model.optimize()
                used_pfba = False

            if sol.status != "optimal":
                return {
                    "status": "infeasible",
                    "items": [],
                    "pearson_r": 0.0,
                    "rmse": 0.0,
                    "mean_deviation_pct": 0.0,
                    "warnings": [f"Model optimization failed (status: {sol.status}). Ensure nutrients are open."]
                }

            # FVA: compute min/max flux range for MFA reactions
            # fraction_of_optimum=0.95 allows 5% suboptimality for wider ranges
            fva_results = {}
            try:
                from cobra.flux_analysis import flux_variability_analysis
                mfa_rxn_ids = []
                for entry in MFA_LITERATURE_DATASET:
                    aliases = REACTION_ID_ALIASES.get(entry["reaction_id"], [entry["reaction_id"]])
                    for alias in aliases:
                        for cand in [alias, alias+"_num1", alias+"_num2"]:
                            if cand in model.reactions:
                                mfa_rxn_ids.append(cand)
                                break

                if mfa_rxn_ids:
                    fva_df = flux_variability_analysis(
                        model, reaction_list=mfa_rxn_ids,
                        fraction_of_optimum=0.95, processes=1
                    )
                    for rxn_id, row in fva_df.iterrows():
                        fva_results[rxn_id] = {
                            "min": float(row["minimum"]),
                            "max": float(row["maximum"])
                        }
            except Exception as fva_err:
                warnings_list.append(f"FVA skipped: {fva_err}")

            # Check if actual glucose uptake flux is available for scaling
            glc_flux_actual = None
            for glc_id in ["EX_glc_e", "EX_glc__D_e"]:
                if glc_id in sol.fluxes:
                    glc_flux_actual = abs(float(sol.fluxes[glc_id]))
                    break

            items = []
            sim_fluxes = []
            mfa_fluxes = []

            for entry in MFA_LITERATURE_DATASET:
                rxn_id = entry["reaction_id"]
                aliases = REACTION_ID_ALIASES.get(rxn_id, [rxn_id])

                sim_flux = None
                matched_id = None
                for alias in aliases:
                    if alias in model.reactions:
                        raw = sol.fluxes.get(alias, None)
                        if raw is not None:
                            sim_flux = abs(float(raw))
                            matched_id = alias
                            break

                if sim_flux is None:
                    sim_flux = 0.0
                    warnings_list.append(f"Reaction '{rxn_id}' not found in model.")

                mfa_val = entry["mfa_flux"]
                deviation_pct = ((sim_flux - mfa_val) / mfa_val * 100.0) if mfa_val > 0 else 0.0

                # Look up FVA range for this reaction
                fva_min, fva_max, mfa_in_range = None, None, None
                if matched_id and matched_id in fva_results:
                    fva_min = round(abs(fva_results[matched_id]["min"]), 4)
                    fva_max = round(abs(fva_results[matched_id]["max"]), 4)
                    lo, hi = min(fva_min, fva_max), max(fva_min, fva_max)
                    mfa_in_range = bool(lo <= mfa_val <= hi)

                items.append({
                    "reaction_id": rxn_id,
                    "reaction_name": entry["reaction_name"],
                    "pathway": entry["pathway"],
                    "mfa_flux": mfa_val,
                    "mfa_std": entry["mfa_std"],
                    "sim_flux": round(sim_flux, 4),
                    "deviation_pct": round(deviation_pct, 1),
                    "matched_model_id": matched_id,
                    "reference": entry["reference"],
                    "fva_min": fva_min,
                    "fva_max": fva_max,
                    "mfa_in_range": mfa_in_range
                })

                sim_fluxes.append(sim_flux)
                mfa_fluxes.append(mfa_val)

            # Compute Pearson's R
            n = len(sim_fluxes)
            mean_s = sum(sim_fluxes) / n
            mean_m = sum(mfa_fluxes) / n
            numerator = sum((s - mean_s) * (m - mean_m) for s, m in zip(sim_fluxes, mfa_fluxes))
            denom_s = math.sqrt(sum((s - mean_s) ** 2 for s in sim_fluxes))
            denom_m = math.sqrt(sum((m - mean_m) ** 2 for m in mfa_fluxes))
            pearson_r = numerator / (denom_s * denom_m) if (denom_s * denom_m) > 0 else 0.0

            rmse = math.sqrt(sum((s - m) ** 2 for s, m in zip(sim_fluxes, mfa_fluxes)) / n)
            mean_dev = sum(abs((s - m) / m * 100.0) for s, m in zip(sim_fluxes, mfa_fluxes) if m > 0) / n

            fva_covered = sum(1 for it in items if it.get("mfa_in_range") is True)
            fva_total   = sum(1 for it in items if it.get("mfa_in_range") is not None)
            fva_coverage = round(100.0 * fva_covered / fva_total, 1) if fva_total > 0 else None

            return {
                "status": "success",
                "items": items,
                "pearson_r": round(pearson_r, 4),
                "rmse": round(rmse, 4),
                "mean_deviation_pct": round(mean_dev, 2),
                "fva_coverage_pct": fva_coverage,
                "used_pfba": used_pfba,
                "warnings": warnings_list
            }

    except Exception as e:
        logger.error(f"MFA comparison failed: {str(e)}")
        return {
            "status": "error",
            "items": [],
            "pearson_r": 0.0,
            "rmse": 0.0,
            "mean_deviation_pct": 0.0,
            "warnings": [f"MFA comparison failed: {str(e)}"]
        }
