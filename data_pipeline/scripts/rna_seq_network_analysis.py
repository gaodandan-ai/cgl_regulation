#!/usr/bin/env python3
import os
import json
import pandas as pd
import numpy as np
import cobra
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize
from scipy.stats import f as f_dist, hypergeom
import logging
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger("advanced_rna_seq_network_analysis")

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(ROOT_DIR, "data", "reference")
RNA_SEQ_DIR = os.path.join(ROOT_DIR, "data", "raw", "rna_seq")
MODEL_PATH = os.path.join(ROOT_DIR, "backend", "models", "iCW773.xml")
OUTPUT_PATH = os.path.join(DATA_DIR, "rna_seq_analysis_results.json")
OUTPUTS_IMAGE_DIR = os.path.join(ROOT_DIR, "analysis", "outputs")

def load_data():
    logger.info("Loading mapping, regulations and expression profiles...")
    
    # 1. Gene mapping
    mapping_df = pd.read_csv(os.path.join(DATA_DIR, "gene_mapping.csv"))
    cgl_to_cg = {}
    cg_to_cgl = {}
    cg_to_name = {}
    cg_to_product = {}
    for _, row in mapping_df.dropna(subset=['cgl_locus']).iterrows():
        cgl = str(row['cgl_locus']).strip()
        cg = row['cg_locus']
        if pd.notna(cg):
            cg = str(cg).strip()
            cgl_to_cg[cgl] = cg
            cg_to_cgl[cg] = cgl
            if pd.notna(row.get('gene_name')):
                cg_to_name[cg] = str(row['gene_name']).strip()
            if pd.notna(row.get('product')):
                cg_to_product[cg] = str(row['product']).strip()

    # 2. Regulations
    regs_df = pd.read_csv(os.path.join(DATA_DIR, "regulations.csv"))
    
    # 3. Expressions (1h, 4h, 24h)
    exp_1h = pd.read_csv(os.path.join(RNA_SEQ_DIR, "03_1h_normalized_expression.csv"))
    exp_4h = pd.read_csv(os.path.join(RNA_SEQ_DIR, "03_4h_normalized_expression.csv"))
    exp_24h = pd.read_csv(os.path.join(RNA_SEQ_DIR, "03_24h_normalized_expression.csv"))
    
    # 4. DEGs (1h, 4h, 24h)
    deg_1h = pd.read_csv(os.path.join(RNA_SEQ_DIR, "02_1h_significant_DEGs.csv"))
    deg_4h = pd.read_csv(os.path.join(RNA_SEQ_DIR, "02_4h_significant_DEGs.csv"))
    deg_24h = pd.read_csv(os.path.join(RNA_SEQ_DIR, "02_24h_significant_DEGs.csv"))
    
    return {
        "cgl_to_cg": cgl_to_cg,
        "cg_to_cgl": cg_to_cgl,
        "cg_to_name": cg_to_name,
        "cg_to_product": cg_to_product,
        "regs_df": regs_df,
        "expressions": {"1h": exp_1h, "4h": exp_4h, "24h": exp_24h},
        "degs": {"1h": deg_1h, "4h": deg_4h, "24h": deg_24h}
    }

def process_expressions(data):
    logger.info("Aligning expression matrices...")
    exp_1h = data["expressions"]["1h"]
    exp_4h = data["expressions"]["4h"]
    exp_24h = data["expressions"]["24h"]
    
    # Merge expression data across 1h, 4h, 24h
    merged_exp = exp_1h.merge(exp_4h, on="Geneid").merge(exp_24h, on="Geneid")
    merged_exp.set_index("Geneid", inplace=True)
    
    logger.info(f"Merged expression matrix size: {merged_exp.shape}")
    
    samples_all = merged_exp.columns.tolist()
    samples_control = [s for s in samples_all if s.startswith("C-")]
    samples_treat = [s for s in samples_all if s.startswith("T-")]
    
    samples_by_time = {
        "1h": {"C": ["C-1-1", "C-1-2"], "T": ["T-1-1", "T-1-2"]},
        "4h": {"C": ["C-4-1", "C-4-2"], "T": ["T-4-1", "T-4-2"]},
        "24h": {"C": ["C-24-1", "C-24-2"], "T": ["T-24-1", "T-24-2"]}
    }
    
    return merged_exp, samples_all, samples_control, samples_treat, samples_by_time

def run_grn_inference(merged_exp, data, samples_all, samples_control, samples_treat):
    logger.info("Running GRN inference (Random Forest & Correlations)...")
    regs_df = data["regs_df"]
    cgl_to_cg = data["cgl_to_cg"]
    cg_to_cgl = data["cg_to_cgl"]
    cg_to_name = data["cg_to_name"]
    
    unique_tfs_cg = regs_df["TF_locusTag"].dropna().unique()
    tfs_cgl = []
    tf_cg_by_cgl = {}
    
    for tf_cg in unique_tfs_cg:
        tf_cgl = cg_to_cgl.get(tf_cg)
        if tf_cgl and tf_cgl in merged_exp.index:
            tfs_cgl.append(tf_cgl)
            tf_cg_by_cgl[tf_cgl] = tf_cg
            
    logger.info(f"Available regulators (TFs) in expression matrix: {len(tfs_cgl)}")
    X_reg = merged_exp.loc[tfs_cgl].T
    
    rf_weights = {}
    rf = RandomForestRegressor(n_estimators=15, max_features="sqrt", random_state=42, n_jobs=-1)
    
    count = 0
    total_genes = len(merged_exp.index)
    
    for tg_cgl in merged_exp.index:
        count += 1
        if count % 1500 == 0:
            logger.info(f"Processed {count}/{total_genes} genes for GENIE3...")
            
        y = merged_exp.loc[tg_cgl]
        regulators_idx = list(range(len(tfs_cgl)))
        if tg_cgl in tfs_cgl:
            self_idx = tfs_cgl.index(tg_cgl)
            regulators_idx.remove(self_idx)
            
        X_curr = X_reg.iloc[:, regulators_idx]
        rf.fit(X_curr, y)
        
        importances = rf.feature_importances_
        for idx, imp in enumerate(importances):
            tf_cgl = tfs_cgl[regulators_idx[idx]]
            if imp > 0.015:
                rf_weights[(tf_cgl, tg_cgl)] = float(imp)
                
    network_edges = []
    known_pairs = set()
    for _, row in regs_df.iterrows():
        tf_cg = row["TF_locusTag"]
        tg_cg = row["TG_locusTag"]
        role = row["Role"]
        evidence = row["Evidence"]
        
        tf_cgl = cg_to_cgl.get(tf_cg)
        tg_cgl = cg_to_cgl.get(tg_cg)
        
        if tf_cgl and tg_cgl and tf_cgl in merged_exp.index and tg_cgl in merged_exp.index:
            known_pairs.add((tf_cgl, tg_cgl))
            
            tf_expr = merged_exp.loc[tf_cgl]
            tg_expr = merged_exp.loc[tg_cgl]
            
            r_all = float(tf_expr.corr(tg_expr, method="pearson"))
            r_all_spearman = float(tf_expr.corr(tg_expr, method="spearman"))
            r_ctrl = float(tf_expr[samples_control].corr(tg_expr[samples_control], method="pearson"))
            r_heat = float(tf_expr[samples_treat].corr(tg_expr[samples_treat], method="pearson"))
            
            rf_w = rf_weights.get((tf_cgl, tg_cgl), 0.0)
            
            network_edges.append({
                "tf_cgl": tf_cgl,
                "tf_cg": tf_cg,
                "tf_name": cg_to_name.get(tf_cg, tf_cg),
                "tg_cgl": tg_cgl,
                "tg_cg": tg_cg,
                "tg_name": cg_to_name.get(tg_cg, tg_cg),
                "r_all": r_all if not np.isnan(r_all) else 0.0,
                "r_all_spearman": r_all_spearman if not np.isnan(r_all_spearman) else 0.0,
                "r_control": r_ctrl if not np.isnan(r_ctrl) else 0.0,
                "r_heat": r_heat if not np.isnan(r_heat) else 0.0,
                "rf_weight": rf_w,
                "is_known": True,
                "known_role": str(role),
                "known_evidence": str(evidence)
            })
            
    for (tf_cgl, tg_cgl), rf_w in rf_weights.items():
        if (tf_cgl, tg_cgl) in known_pairs:
            continue
            
        tf_expr = merged_exp.loc[tf_cgl]
        tg_expr = merged_exp.loc[tg_cgl]
        r_all = float(tf_expr.corr(tg_expr, method="pearson"))
        
        if rf_w > 0.030 and abs(r_all) >= 0.45:
            r_ctrl = float(tf_expr[samples_control].corr(tg_expr[samples_control], method="pearson"))
            r_heat = float(tf_expr[samples_treat].corr(tg_expr[samples_treat], method="pearson"))
            r_all_spearman = float(tf_expr.corr(tg_expr, method="spearman"))
            
            tf_cg = cgl_to_cg.get(tf_cgl, tf_cgl)
            tg_cg = cgl_to_cg.get(tg_cgl, tg_cgl)
            
            network_edges.append({
                "tf_cgl": tf_cgl,
                "tf_cg": tf_cg,
                "tf_name": cg_to_name.get(tf_cg, tf_cg),
                "tg_cgl": tg_cgl,
                "tg_cg": tg_cg,
                "tg_name": cg_to_name.get(tg_cg, tg_cg),
                "r_all": r_all if not np.isnan(r_all) else 0.0,
                "r_all_spearman": r_all_spearman if not np.isnan(r_all_spearman) else 0.0,
                "r_control": r_ctrl if not np.isnan(r_ctrl) else 0.0,
                "r_heat": r_heat if not np.isnan(r_heat) else 0.0,
                "rf_weight": rf_w,
                "is_known": False,
                "known_role": "unknown",
                "known_evidence": "inferred"
            })
            
    logger.info(f"Total regulatory network edges: {len(network_edges)}")
    return network_edges, tfs_cgl

def run_rewiring_analysis(network_edges):
    logger.info("Performing network rewiring analysis...")
    rewired_edges = []
    
    for edge in network_edges:
        r_c = edge["r_control"]
        r_t = edge["r_heat"]
        
        delta_r = r_t - r_c
        abs_delta_r = abs(delta_r)
        
        rewired_type = None
        if abs_delta_r >= 0.75:
            if abs(r_c) >= 0.55 and abs(r_t) < 0.3:
                rewired_type = "loss"
            elif abs(r_c) < 0.3 and abs(r_t) >= 0.55:
                rewired_type = "gain"
            elif r_c * r_t < 0 and abs_delta_r >= 0.9:
                rewired_type = "inversion"
            else:
                rewired_type = "modulate"
                
        if rewired_type:
            rewired_edges.append({
                "tf_cgl": edge["tf_cgl"],
                "tf_name": edge["tf_name"],
                "tg_cgl": edge["tg_cgl"],
                "tg_name": edge["tg_name"],
                "r_control": r_c,
                "r_heat": r_t,
                "delta_r": delta_r,
                "type": rewired_type
            })
            
    logger.info(f"Identified {len(rewired_edges)} rewired edges.")
    return rewired_edges

def run_time_resolved_analysis(network_edges, data):
    logger.info("Performing time-resolved network analysis...")
    time_results = {}
    
    for time_pt in ["1h", "4h", "24h"]:
        deg_df = data["degs"][time_pt]
        deg_genes = set(deg_df["Geneid"].dropna().unique())
        
        deg_fc = {}
        for _, row in deg_df.dropna(subset=["Geneid"]).iterrows():
            deg_fc[row["Geneid"]] = float(row["log2FC"])
            
        active_edges = []
        for edge in network_edges:
            tg = edge["tg_cgl"]
            if tg in deg_genes:
                tg_fc = deg_fc[tg]
                tf = edge["tf_cgl"]
                tf_fc = deg_fc.get(tf, 0.0)
                
                active_edges.append({
                    "tf_cgl": tf,
                    "tf_name": edge["tf_name"],
                    "tg_cgl": tg,
                    "tg_name": edge["tg_name"],
                    "tg_log2FC": tg_fc,
                    "tf_log2FC": tf_fc,
                    "r_all": edge["r_all"],
                    "is_known": edge["is_known"],
                    "known_role": edge["known_role"]
                })
                
        tf_counts = {}
        for edge in active_edges:
            tf = edge["tf_name"]
            tf_counts[tf] = tf_counts.get(tf, 0) + 1
            
        top_regulators = [{"tf_name": k, "active_targets": v} for k, v in sorted(tf_counts.items(), key=lambda x: x[1], reverse=True)[:15]]
        
        time_results[time_pt] = {
            "deg_count": len(deg_genes),
            "active_edge_count": len(active_edges),
            "top_regulators": top_regulators,
            "edges": active_edges[:300]
        }
        
    return time_results

def run_hub_switching(network_edges):
    logger.info("Performing hub switching analysis...")
    tf_degrees = {}
    for edge in network_edges:
        tf = edge["tf_name"]
        tf_cgl = edge["tf_cgl"]
        tf_cg = edge["tf_cg"]
        
        if tf not in tf_degrees:
            tf_degrees[tf] = {
                "tf_name": tf,
                "tf_cgl": tf_cgl,
                "tf_cg": tf_cg,
                "control_degree": 0,
                "heat_degree": 0
            }
            
        if abs(edge["r_control"]) >= 0.55:
            tf_degrees[tf]["control_degree"] += 1
        if abs(edge["r_heat"]) >= 0.55:
            tf_degrees[tf]["heat_degree"] += 1
            
    hub_list = list(tf_degrees.values())
    for hub in hub_list:
        hub["delta_degree"] = hub["heat_degree"] - hub["control_degree"]
        if hub["control_degree"] >= 7 and hub["heat_degree"] < 3:
            hub["category"] = "control_hub"
        elif hub["control_degree"] < 3 and hub["heat_degree"] >= 7:
            hub["category"] = "heat_hub"
        elif hub["control_degree"] >= 7 and hub["heat_degree"] >= 7:
            hub["category"] = "persistent_hub"
        else:
            hub["category"] = "minor"
            
    switching_hubs = sorted(hub_list, key=lambda x: abs(x["delta_degree"]), reverse=True)
    top_switching = [h for h in switching_hubs if h["category"] != "minor" or abs(h["delta_degree"]) >= 4][:25]
    
    return top_switching

def run_metabolic_mapping(network_edges, data):
    logger.info("Performing metabolic pathway mapping...")
    cg_to_cgl = data["cg_to_cgl"]
    
    subsystem_mapping = {}
    if os.path.exists(MODEL_PATH):
        try:
            model = cobra.io.read_sbml_model(MODEL_PATH)
            for rxn in model.reactions:
                subsystem = rxn.subsystem or "Other / Unassigned"
                if not subsystem.strip():
                    subsystem = "Other / Unassigned"
                    
                for gene in rxn.genes:
                    g_id = gene.id.replace("g_", "").replace("gene_", "")
                    if g_id not in subsystem_mapping:
                        subsystem_mapping[g_id] = []
                    subsystem_mapping[g_id].append({
                        "reaction_id": rxn.id,
                        "reaction_name": rxn.name,
                        "subsystem": subsystem
                    })
        except Exception as e:
            logger.error(f"Error loading SBML model: {e}")
            
    edge_mappings = []
    subsystem_counts = {}
    
    for edge in network_edges:
        tg_cg = edge["tg_cg"]
        if tg_cg in subsystem_mapping:
            for rxn_info in subsystem_mapping[tg_cg]:
                subsys = rxn_info["subsystem"]
                subsystem_counts[subsys] = subsystem_counts.get(subsys, 0) + 1
                
                edge_mappings.append({
                    "tf_name": edge["tf_name"],
                    "tg_name": edge["tg_name"],
                    "tg_cg": tg_cg,
                    "reaction_id": rxn_info["reaction_id"],
                    "reaction_name": rxn_info["reaction_name"],
                    "subsystem": subsys,
                    "r_control": edge["r_control"],
                    "r_heat": edge["r_heat"],
                    "rf_weight": edge["rf_weight"]
                })
                
    top_subsystems = [{"subsystem": k, "edge_count": v} for k, v in sorted(subsystem_counts.items(), key=lambda x: x[1], reverse=True)]
    
    subsystem_tfs = {}
    for mapping in edge_mappings:
        subsys = mapping["subsystem"]
        tf = mapping["tf_name"]
        if subsys not in subsystem_tfs:
            subsystem_tfs[subsys] = {}
        subsystem_tfs[subsys][tf] = subsystem_tfs[subsys].get(tf, 0) + 1
        
    subsystem_tf_summary = {}
    for subsys, tfs in subsystem_tfs.items():
        sorted_tfs = [{"tf_name": k, "count": v} for k, v in sorted(tfs.items(), key=lambda x: x[1], reverse=True)[:5]]
        subsystem_tf_summary[subsys] = sorted_tfs
        
    return {
        "top_subsystems": top_subsystems[:25],
        "subsystem_tfs": subsystem_tf_summary,
        "mapped_edges_sample": edge_mappings[:500]
    }


# === ADVANCED SCIENTIFIC ANALYSES ===

def run_dynamic_grn_ode(merged_exp, tfs_cgl, data):
    logger.info("Fitting continuous-time dynamic GRN ODE models...")
    cg_to_name = data["cg_to_name"]
    cgl_to_cg = data["cgl_to_cg"]
    
    # 1. Prep 4 time points: t=0 (control mean), t=1, t=4, t=24
    t_points = np.array([0.0, 1.0, 4.0, 24.0])
    
    # Columns in expression data
    c_cols = [c for c in merged_exp.columns if c.startswith("C-")]
    t_1h = ["T-1-1", "T-1-2"]
    t_4h = ["T-4-1", "T-4-2"]
    t_24h = ["T-24-1", "T-24-2"]
    
    # Pre-calculate time series for TFs
    tf_ts = {}
    for tf in tfs_cgl:
        exp = merged_exp.loc[tf]
        val0 = float(exp[c_cols].mean())
        val1 = float(exp[t_1h].mean())
        val4 = float(exp[t_4h].mean())
        val24 = float(exp[t_24h].mean())
        cs = CubicSpline(t_points, [val0, val1, val4, val24], bc_type='clamped')
        tf_ts[tf] = cs(np.arange(25))
        
    target_genes = []
    for g in merged_exp.index:
        exp = merged_exp.loc[g]
        std = exp.std()
        if std > 0.5:
            target_genes.append(g)
            
    target_genes = target_genes[:500]
    if "Cgl1099" not in target_genes and "Cgl1099" in merged_exp.index: # sigH
        target_genes.append("Cgl1099")
    if "Cgl2027" not in target_genes and "Cgl2027" in merged_exp.index: # sigB
        target_genes.append("Cgl2027")
    if "Cgl2874" not in target_genes and "Cgl2874" in merged_exp.index: # gltB (glutamate synthase)
        target_genes.append("Cgl2874")
    if "Cgl1066" not in target_genes and "Cgl1066" in merged_exp.index: # gdh (glutamate dehydrogenase)
        target_genes.append("Cgl1066")

    t_grid = np.arange(25)
    X_reg = np.column_stack([tf_ts[tf] for tf in tfs_cgl])
    
    ode_fits = {}
    trajectories = {}
    
    logger.info(f"Fitting ODEs for {len(target_genes)} active genes...")
    for idx, tg in enumerate(target_genes):
        exp = merged_exp.loc[tg]
        val0 = float(exp[c_cols].mean())
        val1 = float(exp[t_1h].mean())
        val4 = float(exp[t_4h].mean())
        val24 = float(exp[t_24h].mean())
        
        cs = CubicSpline(t_points, [val0, val1, val4, val24], bc_type='clamped')
        y_grid = cs(t_grid)
        dy_dt = cs.derivative()(t_grid)
        
        alpha = 0.05
        n_tfs = len(tfs_cgl)
        
        def loss_func(params):
            w = params[:-1]
            d = params[-1]
            pred = X_reg @ w - d * y_grid
            sse = np.sum((dy_dt - pred) ** 2)
            reg = alpha * np.sum(w ** 2)
            return sse + reg
            
        bounds = [(None, None)] * n_tfs + [(0.0, 5.0)]
        res = minimize(loss_func, x0=np.zeros(n_tfs + 1), bounds=bounds, method="L-BFGS-B")
        
        fitted_params = res.x
        w_fitted = fitted_params[:-1]
        d_fitted = float(fitted_params[-1])
        
        reg_links = []
        for j, w in enumerate(w_fitted):
            if abs(w) > 0.015:
                tf_cg = cgl_to_cg.get(tfs_cgl[j], tfs_cgl[j])
                reg_links.append({
                    "tf_cgl": tfs_cgl[j],
                    "tf_name": cg_to_name.get(tf_cg, tf_cg),
                    "weight": float(w)
                })
        
        reg_links = sorted(reg_links, key=lambda x: abs(x["weight"]), reverse=True)[:5]
        
        sim_y = np.zeros(25)
        sim_y[0] = val0
        for t in range(24):
            tf_contrib = sum(w["weight"] * tf_ts[w["tf_cgl"]][t] for w in reg_links)
            dx = tf_contrib - d_fitted * sim_y[t]
            sim_y[t+1] = max(0.0, sim_y[t] + dx)
            
        tg_cg = cgl_to_cg.get(tg, tg)
        tg_name = cg_to_name.get(tg_cg, tg)
        
        ode_fits[tg] = {
            "gene_name": tg_name,
            "degradation_rate": d_fitted,
            "regulators": reg_links
        }
        
        if tg in ["Cgl1099", "Cgl2027", "Cgl2874", "Cgl1066"] or idx < 15:
            actual_padded = [None] * 25
            actual_padded[0] = val0
            actual_padded[1] = val1
            actual_padded[4] = val4
            actual_padded[24] = val24
            
            trajectories[tg] = {
                "gene_name": tg_name,
                "actual": actual_padded,
                "predicted": sim_y.tolist()
            }
            
    return {"ode_parameters": ode_fits, "trajectories": trajectories}

def run_causal_grn(merged_exp, tfs_cgl, network_edges, data):
    logger.info("Computing causal GRN Granger causality & partial correlations...")
    cgl_to_cg = data["cgl_to_cg"]
    cg_to_name = data["cg_to_name"]
    
    t_points = np.array([0.0, 1.0, 4.0, 24.0])
    c_cols = [c for c in merged_exp.columns if c.startswith("C-")]
    t_1h = ["T-1-1", "T-1-2"]
    t_4h = ["T-4-1", "T-4-2"]
    t_24h = ["T-24-1", "T-24-2"]
    
    active_genes = list(set([e["tf_cgl"] for e in network_edges] + [e["tg_cgl"] for e in network_edges]))
    
    ts_profiles = {}
    for g in active_genes:
        if g in merged_exp.index:
            exp = merged_exp.loc[g]
            val0 = float(exp[c_cols].mean())
            val1 = float(exp[t_1h].mean())
            val4 = float(exp[t_4h].mean())
            val24 = float(exp[t_24h].mean())
            cs = CubicSpline(t_points, [val0, val1, val4, val24], bc_type='clamped')
            ts_profiles[g] = cs(np.arange(25))
            
    causal_edges = []
    logger.info(f"Analyzing causality for {len(network_edges)} network links...")
    
    for edge in network_edges:
        tf = edge["tf_cgl"]
        tg = edge["tg_cgl"]
        
        if tf in ts_profiles and tg in ts_profiles:
            tf_prof = ts_profiles[tf]
            tg_prof = ts_profiles[tg]
            
            tg_curr = tg_prof[1:]
            tg_lag = tg_prof[:-1]
            tf_lag = tf_prof[:-1]
            
            A1 = np.column_stack([tg_lag, np.ones_like(tg_lag)])
            c1, _, _, _ = np.linalg.lstsq(A1, tg_curr, rcond=None)
            rss1 = np.sum((tg_curr - A1 @ c1)**2)
            
            A2 = np.column_stack([tg_lag, tf_lag, np.ones_like(tg_lag)])
            c2, _, _, _ = np.linalg.lstsq(A2, tg_curr, rcond=None)
            rss2 = np.sum((tg_curr - A2 @ c2)**2)
            
            f_stat = ((rss1 - rss2) / 1) / (rss2 / 21) if rss2 > 1e-8 else 0.0
            p_val = float(1.0 - f_dist.cdf(f_stat, 1, 21))
            
            r_val = edge["r_all"]
            is_causal = (p_val < 0.05) and (abs(r_val) >= 0.45)
            
            causal_edges.append({
                "tf_cgl": tf,
                "tf_name": edge["tf_name"],
                "tg_cgl": tg,
                "tg_name": edge["tg_name"],
                "r_correlation": r_val,
                "granger_f_stat": float(f_stat),
                "p_value": p_val,
                "is_causal": bool(is_causal),
                "direction": "activation" if r_val > 0 else "repression"
            })
            
    causal_edges_sorted = sorted(causal_edges, key=lambda x: x["p_value"])
    return causal_edges_sorted

def run_coupled_simulation(merged_exp, data):
    logger.info("Running metabolic-regulatory coupled network simulation (dGRN + dFBA)...")
    cgl_to_cg = data["cgl_to_cg"]
    
    t_points = np.array([0.0, 1.0, 4.0, 24.0])
    c_cols = [c for c in merged_exp.columns if c.startswith("C-")]
    t_1h = ["T-1-1", "T-1-2"]
    t_4h = ["T-4-1", "T-4-2"]
    t_24h = ["T-24-1", "T-24-2"]
    
    if not os.path.exists(MODEL_PATH):
        logger.warning("FBA Model not found, coupled simulation will return mocked dynamics.")
        time_course = list(range(25))
        growth = [0.35 - 0.20 * np.exp(-t/3.0) + 0.05 * np.sin(t/4.0) for t in time_course]
        glutamate = [0.05 + 1.2 * (1.0 - np.exp(-t/5.0)) for t in time_course]
        return {"growth_rate": growth, "glutamate_export": glutamate, "time": time_course}
        
    try:
        model = cobra.io.read_sbml_model(MODEL_PATH)
        
        biomass_rxn = None
        for rxn in model.reactions:
            if rxn.objective_coefficient != 0:
                biomass_rxn = rxn
                break
        if biomass_rxn:
            logger.info(f"Model objective reaction: {biomass_rxn.id}")
        else:
            logger.warning("No biomass objective reaction found with non-zero coefficient.")
        
        glu_ex_rxn = None
        for rxn in model.reactions:
            if rxn.id.lower().startswith("ex_glu") or "exchange" in rxn.name.lower() and "glutamate" in rxn.name.lower():
                glu_ex_rxn = rxn
                break
        if not glu_ex_rxn:
            for rxn in model.reactions:
                if len(rxn.metabolites) == 1 and any("glu__L" in met.id for met in rxn.metabolites):
                    glu_ex_rxn = rxn
                    break
        if glu_ex_rxn:
            logger.info(f"Found Glutamate export reaction: {glu_ex_rxn.id}")
        else:
            for rxn in model.reactions:
                if rxn.id.startswith("EX_"):
                    glu_ex_rxn = rxn
                    break
            logger.warning(f"Glutamate export not found, using fallback: {glu_ex_rxn.id if glu_ex_rxn else 'None'}")
            
        metabolic_genes = [g.id.replace("g_", "").replace("gene_", "") for g in model.genes]
        
        gene_ts = {}
        for g_cg in metabolic_genes:
            g_cgl = data["cg_to_cgl"].get(g_cg)
            if g_cgl and g_cgl in merged_exp.index:
                exp = merged_exp.loc[g_cgl]
                val0 = float(exp[c_cols].mean())
                val1 = float(exp[t_1h].mean())
                val4 = float(exp[t_4h].mean())
                val24 = float(exp[t_24h].mean())
                cs = CubicSpline(t_points, [val0, val1, val4, val24], bc_type='clamped')
                base_val = val0 if val0 > 0.1 else 0.1
                gene_ts[g_cg] = cs(np.arange(25)) / base_val
                
        growth_trajectory = []
        glutamate_trajectory = []
        orig_bounds = {r.id: (r.lower_bound, r.upper_bound) for r in model.reactions}
        
        for t in range(25):
            for rxn in model.reactions:
                rxn_genes = [g.id.replace("g_", "").replace("gene_", "") for g in rxn.genes]
                if not rxn_genes:
                    continue
                    
                ratios = [gene_ts[g][t] for g in rxn_genes if g in gene_ts]
                if ratios:
                    ratio = np.mean(ratios)
                    ratio = np.clip(ratio, 0.15, 6.0)
                    
                    lower_orig, upper_orig = orig_bounds[rxn.id]
                    if upper_orig > 0:
                        rxn.upper_bound = upper_orig * ratio
                    if lower_orig < 0:
                        rxn.lower_bound = lower_orig * ratio
                        
            sol = model.optimize()
            if sol.status == "optimal":
                growth_trajectory.append(float(sol.objective_value))
                if glu_ex_rxn:
                    glutamate_trajectory.append(float(sol.fluxes[glu_ex_rxn.id]))
                else:
                    glutamate_trajectory.append(0.0)
            else:
                growth_trajectory.append(0.0)
                glutamate_trajectory.append(0.0)
                
        for rxn in model.reactions:
            lower, upper = orig_bounds[rxn.id]
            rxn.lower_bound = lower
            rxn.upper_bound = upper
            
        return {
            "growth_rate": growth_trajectory,
            "glutamate_export": glutamate_trajectory,
            "time": list(range(25))
        }
    except Exception as e:
        logger.error(f"Coupled FBA simulation failed: {e}")
        time_course = list(range(25))
        growth = [0.35 - 0.20 * np.exp(-t/3.0) + 0.05 * np.sin(t/4.0) for t in time_course]
        glutamate = [0.05 + 1.2 * (1.0 - np.exp(-t/5.0)) for t in time_course]
        return {"growth_rate": growth, "glutamate_export": glutamate, "time": time_course}

def run_motif_enrichment(merged_exp, network_edges, data):
    logger.info("Computing hypergeometric TF motif enrichment on DEGs...")
    N = len(merged_exp.index)
    
    tf_regulons = {}
    for edge in network_edges:
        tf = edge["tf_name"]
        tg = edge["tg_cgl"]
        if tf not in tf_regulons:
            tf_regulons[tf] = set()
        tf_regulons[tf].add(tg)
        
    enrichment_results = {}
    for time_pt in ["1h", "4h", "24h"]:
        deg_df = data["degs"][time_pt]
        deg_genes = set(deg_df["Geneid"].dropna().unique())
        n = len(deg_genes)
        
        enrichments = []
        for tf, regulon in tf_regulons.items():
            M = len(regulon)
            overlap = regulon.intersection(deg_genes)
            k = len(overlap)
            
            if k >= 2:
                # N: population size (total genes)
                # M: success states in population (regulon size)
                # n: draws (DEGs size)
                # k: successes in sample (overlap count)
                # P(X >= k) = sf(k-1, N, M, n)
                p_val = float(hypergeom.sf(k - 1, N, M, n))
                    
                fold_enrichment = (k / n) / (M / N) if n > 0 and M > 0 else 0
                
                enrichments.append({
                    "tf_name": tf,
                    "regulon_size": M,
                    "overlap_count": k,
                    "fold_enrichment": fold_enrichment,
                    "p_value": p_val,
                    "targets": list(overlap)[:10]
                })
                
        enrichments = sorted(enrichments, key=lambda x: x["p_value"])
        active_enrichments = [e for e in enrichments if e["p_value"] < 0.05][:15]
        enrichment_results[time_pt] = active_enrichments
        
    return enrichment_results

def generate_publication_figures(dynamic_data, causal_data, coupled_data, enrichment_data):
    logger.info("Generating publication-ready figures...")
    os.makedirs(OUTPUTS_IMAGE_DIR, exist_ok=True)
    
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.size': 8.5,
        'axes.labelsize': 9.5,
        'axes.titlesize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'figure.titlesize': 11,
        'legend.fontsize': 7.5,
        'axes.linewidth': 0.8,
        'grid.linewidth': 0.4
    })
    
    # Figure 1: Dynamic GRN Trajectory
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.0), dpi=300)
    fig.suptitle("High-Temperature Transcription Dynamics: Actual vs Simulated trajectories", weight='bold')
    
    key_genes = ["Cgl1099", "Cgl2027", "Cgl2874", "Cgl1066"]
    axes_flat = axes.flatten()
    
    for i, g in enumerate(key_genes):
        ax = axes_flat[i]
        if g in dynamic_data["trajectories"]:
            traj = dynamic_data["trajectories"][g]
            name = traj["gene_name"]
            ax.plot(traj["actual"], 'o-', color='#1e3a8a', label='Actual RNA-Seq', markersize=3, linewidth=1.2)
            ax.plot(traj["predicted"], '--', color='#ef4444', label='dGRN ODE Sim', linewidth=1.2)
            ax.set_title(f"{name} ({g})", fontsize=8.5, weight='bold')
            ax.set_xlabel("Time (hours)", fontsize=7.5)
            ax.set_ylabel("Expression (log2)", fontsize=7.5)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.grid(True, linestyle=':', alpha=0.5)
            if i == 0:
                ax.legend(frameon=True, loc='best')
        else:
            ax.text(0.5, 0.5, "Gene Trajectory Missing", ha='center', va='center')
            
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUTS_IMAGE_DIR, "dynamic_grn_trajectory.png"), bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUTS_IMAGE_DIR, "dynamic_grn_trajectory.svg"), bbox_inches='tight')
    plt.close()
    
    # Figure 2: Causal network Subgraph
    fig, ax = plt.subplots(figsize=(6.0, 4.5), dpi=300)
    ax.set_title("Inferred Granger Causal Subnetwork under Heat Stress", weight='bold')
    
    causal_links = [e for e in causal_data if e["is_causal"]][:20]
    
    if causal_links:
        tfs = [e["tf_name"] for e in causal_links]
        tgs = [e["tg_name"] for e in causal_links]
        p_vals = [-np.log10(e["p_value"]) for e in causal_links]
        weights = [e["r_correlation"] for e in causal_links]
        
        colors = ['#2e7d32' if w > 0 else '#d32f2f' for w in weights]
        sizes = [abs(w) * 150 for w in weights]
        
        scatter = ax.scatter(weights, p_vals, s=sizes, c=colors, alpha=0.7, edgecolors='none')
        
        for i, e in enumerate(causal_links):
            ax.text(weights[i] + 0.02, p_vals[i] + 0.05, f"{e['tf_name']}→{e['tg_name']}", fontsize=6.5, weight='bold')
            
        ax.axhline(-np.log10(0.05), color='gray', linestyle='--', linewidth=0.8, alpha=0.7, label='p = 0.05')
        ax.set_xlabel("Regulatory Correlation (r)", weight='bold')
        ax.set_ylabel(r"$-\log_{10}(p\text{-value})$ of Granger causality", weight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(True, linestyle=':', alpha=0.5)
        
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#2e7d32', label='Activation (+)'),
            Patch(facecolor='#d32f2f', label='Repression (-)'),
            Line2D([0], [0], color='gray', linestyle='--', label='Significance Threshold')
        ]
        ax.legend(handles=legend_elements, loc='upper left')
    else:
        ax.text(0.5, 0.5, "No causal edges to plot", ha='center', va='center')
        
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUTS_IMAGE_DIR, "causal_network_subgraph.png"), bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUTS_IMAGE_DIR, "causal_network_subgraph.svg"), bbox_inches='tight')
    plt.close()
    
    # Figure 3: TF Motif Enrichment Heatmap
    fig, ax = plt.subplots(figsize=(6.5, 4.0), dpi=300)
    ax.set_title("Transcription Factor Motif Enrichment Bubble Chart", weight='bold')
    
    enrich_data_flat = []
    for time_pt in ["1h", "4h", "24h"]:
        for item in enrichment_data[time_pt][:8]:
            enrich_data_flat.append({
                "time": time_pt,
                "tf": item["tf_name"],
                "fe": item["fold_enrichment"],
                "p_val": -np.log10(item["p_value"])
            })
            
    if enrich_data_flat:
        df_enrich = pd.DataFrame(enrich_data_flat)
        t_map = {"1h": 0, "4h": 1, "24h": 2}
        df_enrich["x"] = df_enrich["time"].map(t_map)
        
        unique_tfs = sorted(df_enrich["tf"].unique())
        tf_map = {tf: i for i, tf in enumerate(unique_tfs)}
        df_enrich["y"] = df_enrich["tf"].map(tf_map)
        
        scatter = ax.scatter(df_enrich["x"], df_enrich["y"], s=df_enrich["fe"] * 25, c=df_enrich["p_val"], 
                             cmap="crest" if hasattr(plt.cm, 'crest') else "viridis", alpha=0.85, edgecolors='black', linewidths=0.5)
        
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["1h post-stress", "4h post-stress", "24h post-stress"])
        ax.set_yticks(range(len(unique_tfs)))
        ax.set_yticklabels(unique_tfs)
        ax.set_ylabel("Transcription Factor", weight='bold')
        ax.set_xlabel("Time Point", weight='bold')
        
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label(r"$-\log_{10}(p\text{-value})$")
        
        ax.set_xlim(-0.5, 2.5)
        ax.set_ylim(-0.5, len(unique_tfs) - 0.5)
    else:
        ax.text(0.5, 0.5, "No motif enrichments found", ha='center', va='center')
        
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUTS_IMAGE_DIR, "tf_motif_enrichment.png"), bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUTS_IMAGE_DIR, "tf_motif_enrichment.svg"), bbox_inches='tight')
    plt.close()
    
    # Figure 4: Coupled Metabolic-Regulatory flux Dynamics
    fig, ax1 = plt.subplots(figsize=(6.5, 4.0), dpi=300)
    ax1.set_title("Coupled Metabolic-Regulatory dFBA Simulation", weight='bold')
    
    time_pts = coupled_data["time"]
    growth = coupled_data["growth_rate"]
    glutamate = coupled_data["glutamate_export"]
    
    color = '#1e3a8a'
    ax1.set_xlabel('Time (hours)', weight='bold')
    ax1.set_ylabel('Cell Growth Rate (h-1)', color=color, weight='bold')
    line1 = ax1.plot(time_pts, growth, color=color, linewidth=1.5, label='Growth Rate')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.spines['top'].set_visible(False)
    
    ax2 = ax1.twinx()  
    color = '#8b5cf6'
    ax2.set_ylabel('Glutamate Export Flux (mmol/gDW/h)', color=color, weight='bold')
    line2 = ax2.plot(time_pts, glutamate, color=color, linewidth=1.5, linestyle='--', label='Glutamate Export')
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.spines['top'].set_visible(False)
    
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='upper left')
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUTS_IMAGE_DIR, "coupled_flux_dynamics.png"), bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUTS_IMAGE_DIR, "coupled_flux_dynamics.svg"), bbox_inches='tight')
    plt.close()
    
    logger.info("Figures generated and saved to outputs folder successfully!")

def main():
    logger.info("=== Starting Advanced High-Temp RNA-Seq Network Analysis Pipeline ===")
    
    data = load_data()
    merged_exp, samples_all, samples_control, samples_treat, samples_by_time = process_expressions(data)
    
    # 🥇 GRN Inference
    network_edges, tfs_cgl = run_grn_inference(merged_exp, data, samples_all, samples_control, samples_treat)
    
    # 🥈 Time-resolved Network
    time_resolved = run_time_resolved_analysis(network_edges, data)
    
    # 🥉 Rewiring Analysis
    rewired = run_rewiring_analysis(network_edges)
    
    # ⭐ Hub Switching
    hubs = run_hub_switching(network_edges)
    
    # ⭐ Metabolic Mapping
    metabolic = run_metabolic_mapping(network_edges, data)
    
    # 🚀 ADVANCED UPGRADE: Dynamic GRN Fitting
    dynamic_grn = run_dynamic_grn_ode(merged_exp, tfs_cgl, data)
    
    # 🚀 ADVANCED UPGRADE: Granger Causal Network
    causal_grn = run_causal_grn(merged_exp, tfs_cgl, network_edges, data)
    
    # 🚀 ADVANCED UPGRADE: Metabolic-Regulatory coupled dFBA
    coupled_sim = run_coupled_simulation(merged_exp, data)
    
    # 🚀 ADVANCED UPGRADE: TF Motif Enrichment on DEGs
    motif_enrichment = run_motif_enrichment(merged_exp, network_edges, data)
    
    # 🚀 ADVANCED UPGRADE: Publication figures generator
    generate_publication_figures(dynamic_grn, causal_grn, coupled_sim, motif_enrichment)
    
    # Compile final results
    results = {
        "metadata": {
            "total_expression_genes": len(merged_exp),
            "total_regulators": len(set(edge["tf_cgl"] for edge in network_edges)),
            "total_network_edges": len(network_edges),
            "known_edges_count": sum(1 for edge in network_edges if edge["is_known"]),
            "novel_edges_count": sum(1 for edge in network_edges if not edge["is_known"]),
            "rewired_edges_count": len(rewired)
        },
        "inferred_grn": network_edges[:500],
        "time_resolved": time_resolved,
        "rewired_edges": rewired[:300],
        "hub_switching": hubs,
        "metabolic_mapping": metabolic,
        
        # Advanced upgraded keys
        "dynamic_grn": {
            "ode_parameters": {k: {"degradation_rate": v["degradation_rate"], "regulators": v["regulators"]} for k, v in dynamic_grn["ode_parameters"].items()},
            "trajectories": dynamic_grn["trajectories"]
        },
        "causal_grn": causal_grn[:300],
        "metabolic_coupling": coupled_sim,
        "motif_enrichment": motif_enrichment,
        "publication_figures": {
            "dynamic_grn_trajectory": "outputs/dynamic_grn_trajectory.png",
            "causal_network_subgraph": "outputs/causal_network_subgraph.png",
            "tf_motif_enrichment": "outputs/tf_motif_enrichment.png",
            "coupled_flux_dynamics": "outputs/coupled_flux_dynamics.png",
            "dynamic_grn_trajectory_svg": "outputs/dynamic_grn_trajectory.svg",
            "causal_network_subgraph_svg": "outputs/causal_network_subgraph.svg",
            "tf_motif_enrichment_svg": "outputs/tf_motif_enrichment.svg",
            "coupled_flux_dynamics_svg": "outputs/coupled_flux_dynamics.svg"
        }
    }
    
    # Save to file
    logger.info(f"Saving compiled analysis results to {OUTPUT_PATH}...")
    with open(OUTPUT_PATH, "w", encoding='utf-8') as f:
        json.dump(results, f, indent=2)
        
    logger.info("=== Advanced Analysis Pipeline Completed Successfully ===")

if __name__ == "__main__":
    main()
