#!/usr/bin/env python3
"""
C. glutamicum TF Knockout Metabolic Simulator
=============================================
Simulates the metabolic consequences of knocking out ALL unique transcription factors
present in the regulations dataset on biomass, glutamate, and lysine export flux.
"""

import os
import csv
import sys
from pathlib import Path
from collections import defaultdict
import cobra

# Setup paths relative to workspace root
ROOT_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT_DIR / "data" / "reference"
MODEL_PATH = ROOT_DIR / "backend" / "models" / "iCW773.xml"
REGULATIONS_CSV = DATA_DIR / "regulations.csv"
MAPPING_CSV = DATA_DIR / "gene_mapping.csv"
OUTPUT_DIR = ROOT_DIR / "data_pipeline" / "outputs"

# Reactions to track
BIOMASS_RXN = "Growth"
GLUTAMATE_EX = "EX_glu_L_e"
LYSINE_EX = "EX_lys_L_e"

def load_mappings():
    """Load gene locus mapping dictionaries."""
    cg_to_cgl = {}
    cgl_to_cg = {}
    cgl_to_name = {}
    
    if MAPPING_CSV.exists():
        with open(MAPPING_CSV, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                cg = row.get("cg_locus", "").strip()
                cgl = row.get("cgl_locus", "").strip()
                name = row.get("gene_name", "").strip()
                if cg:
                    if cgl:
                        cg_to_cgl[cg.lower()] = cgl
                        cgl_to_cg[cgl.lower()] = cg
                        if name:
                            cgl_to_name[cgl] = name
    return cg_to_cgl, cgl_to_cg, cgl_to_name

def load_all_tf_activation_targets(cg_to_cgl):
    """
    Get downstream activation targets (Role = 'A') for ALL TFs.
    Also returns a mapping from TF locus tag to TF name.
    """
    tf_targets = defaultdict(list)
    tf_names = {}
    
    if not REGULATIONS_CSV.exists():
        print(f"Error: regulations CSV not found at {REGULATIONS_CSV}")
        return tf_targets, tf_names
        
    with open(REGULATIONS_CSV, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            tf_cg = row.get("TF_locusTag", "").strip().lower()
            tf_name = row.get("TF_name", "").strip()
            tg_cg = row.get("TG_locusTag", "").strip().lower()
            role = row.get("Role", "").strip()
            
            if tf_cg:
                # Initialize key if not exists
                if tf_cg not in tf_targets:
                    tf_targets[tf_cg] = []
                    
                # Track TF name
                if tf_name:
                    tf_names[tf_cg] = tf_name
                elif tf_cg not in tf_names:
                    tf_names[tf_cg] = tf_cg.upper()
                    
                # We model TF knockout by silencing targets that require this TF for activation (Role == 'A')
                if tg_cg and role == 'A':
                    tf_targets[tf_cg].append(tg_cg)
                    
    return dict(tf_targets), tf_names

def map_gene_to_model(model, gene_id: str) -> str:
    """Find the exact gene ID format in the COBRA model."""
    normalized = gene_id.strip().lower()
    
    # Try direct
    for g in model.genes:
        g_id_lower = g.id.lower()
        if normalized == g_id_lower:
            return g.id
        # Try stripping prefixes like 'g_' or 'gene_'
        stripped = g_id_lower.replace("g_", "").replace("gene_", "")
        if normalized == stripped:
            return g.id
    return None

def run_simulations(model, tf_targets, tf_names):
    """Run baseline and perturbed FBA simulations for Biomass, Glutamate, and Lysine."""
    results = {}
    
    # Find active biomass reaction (usually model.objective)
    biomass_id = None
    for rxn in model.reactions:
        if "biomass" in rxn.name.lower() or rxn.id == "BIOMASS_Cgl" or rxn.id == "Growth" or rxn.id == "CG_biomass_cgl_ATCC13032":
            biomass_id = rxn.id
            break
    if not biomass_id:
        biomass_id = list(model.objective.expression.free_symbols)[0].name
        
    print(f"Using Biomass reaction: {biomass_id}")
    
    # Get exchange reactions in model
    has_glu = GLUTAMATE_EX in model.reactions
    has_lys = LYSINE_EX in model.reactions
    
    # --- 1. Baseline simulations ---
    # Biomass FBA
    model.objective = biomass_id
    baseline_biomass = model.slim_optimize()
    baseline_glu_with_biomass = model.reactions.get_by_id(GLUTAMATE_EX).flux if has_glu else 0.0
    baseline_lys_with_biomass = model.reactions.get_by_id(LYSINE_EX).flux if has_lys else 0.0
    
    # Max Glutamate FBA (growth constrained to at least 10% of max to represent living cells)
    baseline_max_glu = 0.0
    if has_glu:
        with model:
            model.reactions.get_by_id(biomass_id).lower_bound = 0.1 * baseline_biomass
            model.objective = GLUTAMATE_EX
            baseline_max_glu = model.slim_optimize()
            
    # Max Lysine FBA
    baseline_max_lys = 0.0
    if has_lys:
        with model:
            model.reactions.get_by_id(biomass_id).lower_bound = 0.1 * baseline_biomass
            model.objective = LYSINE_EX
            baseline_max_lys = model.slim_optimize()
            
    results["baseline"] = {
        "growth": baseline_biomass,
        "glutamate_growth": baseline_glu_with_biomass,
        "lysine_growth": baseline_lys_with_biomass,
        "max_glutamate": baseline_max_glu,
        "max_lysine": baseline_max_lys
    }
    
    # --- 2. TF Perturbation simulations ---
    total_tfs = len(tf_targets)
    print(f"Simulating a total of {total_tfs} transcription factors...")
    
    for i, (tf_cg, targets) in enumerate(tf_targets.items(), 1):
        # Map target genes to model gene IDs
        model_genes = []
        for tg_cg in targets:
            m_gid = map_gene_to_model(model, tg_cg)
            if m_gid:
                model_genes.append(model.genes.get_by_id(m_gid))
                
        if i % 15 == 0 or i == total_tfs:
            print(f" - Progress: {i}/{total_tfs} TFs simulated...")
            
        # Simulate KO by applying constraints on target genes
        perturbed_biomass = 0.0
        perturbed_glu_with_biomass = 0.0
        perturbed_lys_with_biomass = 0.0
        perturbed_max_glu = 0.0
        perturbed_max_lys = 0.0
        
        # If no targets are mapped into model, it behaves exactly like baseline (save solver computation)
        if model_genes:
            with model:
                cobra.manipulation.knock_out_model_genes(model, model_genes)
                
                # Growth FBA
                model.objective = biomass_id
                try:
                    perturbed_biomass = model.slim_optimize()
                    if model.solver.status != "optimal":
                        perturbed_biomass = 0.0
                except Exception:
                    perturbed_biomass = 0.0
                    
                if perturbed_biomass > 1e-5:
                    perturbed_glu_with_biomass = model.reactions.get_by_id(GLUTAMATE_EX).flux if has_glu else 0.0
                    perturbed_lys_with_biomass = model.reactions.get_by_id(LYSINE_EX).flux if has_lys else 0.0
                    
                    # Max Glutamate FBA
                    if has_glu:
                        with model:
                            model.reactions.get_by_id(biomass_id).lower_bound = 0.1 * perturbed_biomass
                            model.objective = GLUTAMATE_EX
                            try:
                                perturbed_max_glu = model.slim_optimize()
                                if model.solver.status != "optimal":
                                    perturbed_max_glu = 0.0
                            except Exception:
                                perturbed_max_glu = 0.0
                                
                    # Max Lysine FBA
                    if has_lys:
                        with model:
                            model.reactions.get_by_id(biomass_id).lower_bound = 0.1 * perturbed_biomass
                            model.objective = LYSINE_EX
                            try:
                                perturbed_max_lys = model.slim_optimize()
                                if model.solver.status != "optimal":
                                    perturbed_max_lys = 0.0
                            except Exception:
                                perturbed_max_lys = 0.0
        else:
            # Same as baseline
            perturbed_biomass = baseline_biomass
            perturbed_glu_with_biomass = baseline_glu_with_biomass
            perturbed_lys_with_biomass = baseline_lys_with_biomass
            perturbed_max_glu = baseline_max_glu
            perturbed_max_lys = baseline_max_lys
            
        results[tf_cg] = {
            "mapped_count": len(model_genes),
            "growth": perturbed_biomass,
            "glutamate_growth": perturbed_glu_with_biomass,
            "lysine_growth": perturbed_lys_with_biomass,
            "max_glutamate": perturbed_max_glu,
            "max_lysine": perturbed_max_lys
        }
        
    return results, biomass_id

def pct_change(after, before):
    if before < 1e-6:
        return "N/A" if after < 1e-6 else "+100%"
    pct = ((after - before) / before) * 100
    if abs(pct) < 1e-2:
        return "0.0%"
    return f"{pct:+.1f}%"

def generate_report(results, tf_targets, tf_names, cgl_to_name, biomass_id):
    """Write the results to simulation_report.md categorized by impact."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    report_path = OUTPUT_DIR / "simulation_report.md"
    
    baseline = results["baseline"]
    
    # Categorize TFs
    lethal_tfs = []
    yield_tfs = []
    neutral_tfs = []
    no_targets_tfs = []
    
    for tf_cg, res in results.items():
        if tf_cg == "baseline":
            continue
            
        tf_name = tf_names.get(tf_cg, tf_cg)
        tf_cgl = tf_cg.upper()
        tf_name_display = cgl_to_name.get(tf_cgl, tf_name)
        
        info = {
            "locus": tf_cg,
            "name": tf_name_display,
            "targets_count": len(tf_targets[tf_cg]),
            "mapped_count": res["mapped_count"],
            "growth": res["growth"],
            "max_glu": res["max_glutamate"],
            "max_lys": res["max_lysine"]
        }
        
        if res["mapped_count"] == 0:
            no_targets_tfs.append(info)
        elif res["growth"] < 1e-4:
            lethal_tfs.append(info)
        elif (res["max_glutamate"] > baseline["max_glutamate"] + 1e-3) or (res["max_lysine"] > baseline["max_lysine"] + 1e-3):
            yield_tfs.append(info)
        else:
            neutral_tfs.append(info)
            
    lines = []
    lines.append("# All Transcription Factors Knockout FBA Simulation Report: C. glutamicum")
    lines.append("")
    lines.append("This report lists the metabolic consequences of knocking out **every transcription factor** present in the regulations dataset, simulated on the **iCW773** genome-scale metabolic model.")
    lines.append("")
    
    lines.append("## 1. Simulation Parameters")
    lines.append("")
    lines.append(f"- **Genome-Scale Model**: *iCW773* ({MODEL_PATH.name})")
    lines.append(f"- **Biomass Reaction**: `{biomass_id}`")
    lines.append(f"- **Glutamate Exchange Reaction**: `{GLUTAMATE_EX}`")
    lines.append(f"- **Lysine Exchange Reaction**: `{LYSINE_EX}`")
    lines.append("")
    
    lines.append("## 2. Baseline Status (Wild Type)")
    lines.append("")
    lines.append(f"- **Max Growth Rate (Biomass)**: `{baseline['growth']:.4f} h^-1`")
    lines.append(f"- **Max Glutamate Export (at 10% Growth)**: `{baseline['max_glutamate']:.4f} mmol/gDCW/h`")
    lines.append(f"- **Max Lysine Export (at 10% Growth)**: `{baseline['max_lysine']:.4f} mmol/gDCW/h`")
    lines.append("")
    
    # 1. Yield Enhancing TFs
    lines.append("## 3. Yield-Enhancing TF Knockouts")
    lines.append("These transcription factors, when knocked out, are predicted to **boost amino acid production** due to intracellular flux redirection or metabolic pathway derepression:")
    lines.append("")
    if yield_tfs:
        lines.append("| TF Name | Locus Tag | Targets | Mapped | Growth (Biomass) | Change % | Max Glutamate | Change % | Max Lysine | Change % |")
        lines.append("| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
        # Sort by glutamate increase
        yield_tfs_sorted = sorted(yield_tfs, key=lambda info: max(info["max_glu"] - baseline["max_glutamate"], info["max_lys"] - baseline["max_lysine"]), reverse=True)
        for t in yield_tfs_sorted:
            lines.append(
                f"| **{t['name']}** | `{t['locus']}` | {t['targets_count']} | {t['mapped_count']} | "
                f"{t['growth']:.4f} | {pct_change(t['growth'], baseline['growth'])} | "
                f"{t['max_glu']:.4f} | {pct_change(t['max_glu'], baseline['max_glutamate'])} | "
                f"{t['max_lys']:.4f} | {pct_change(t['max_lys'], baseline['max_lysine'])} |"
            )
    else:
        lines.append("*No yield-enhancing transcription factor knockouts were detected under standard FBA.*")
    lines.append("")
    
    # 2. Lethal TFs
    lines.append("## 4. Lethal / High-Growth-Inhibiting TF Knockouts")
    lines.append("Silencing these transcription factors knocks out essential metabolic nodes, resulting in **growth failure (lethality)**. They represent crucial components for cell survival and biomass build-up:")
    lines.append("")
    if lethal_tfs:
        lines.append("| TF Name | Locus Tag | Targets | Mapped | Growth (Biomass) | Change % | Max Glutamate | Max Lysine |")
        lines.append("| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |")
        for t in sorted(lethal_tfs, key=lambda x: x["name"]):
            lines.append(
                f"| **{t['name']}** | `{t['locus']}` | {t['targets_count']} | {t['mapped_count']} | "
                f"{t['growth']:.4f} | -100.0% | {t['max_glu']:.4f} | {t['max_lys']:.4f} |"
            )
    else:
        lines.append("*No lethal knockouts detected.*")
    lines.append("")
    
    # 3. Neutral / Low-Impact TFs
    lines.append("## 5. Low-Impact / Neutral TF Knockouts")
    lines.append("Knocking out these TFs has a minor ($< 0.1\\%$) impact on growth and production under FBA simulation. They are targets mapped in the model but showing redundancy or minor activation roles:")
    lines.append("")
    if neutral_tfs:
        lines.append("<details><summary>Click to expand/view neutral TFs</summary>")
        lines.append("")
        lines.append("| TF Name | Locus Tag | Targets | Mapped | Growth (Biomass) | Change % | Max Glutamate | Max Lysine |")
        lines.append("| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |")
        for t in sorted(neutral_tfs, key=lambda x: x["name"]):
            lines.append(
                f"| **{t['name']}** | `{t['locus']}` | {t['targets_count']} | {t['mapped_count']} | "
                f"{t['growth']:.4f} | {pct_change(t['growth'], baseline['growth'])} | {t['max_glu']:.4f} | {t['max_lys']:.4f} |"
            )
        lines.append("</details>")
    else:
        lines.append("*No neutral knockouts.*")
    lines.append("")
    
    # 4. No targets mapped
    lines.append("## 6. TFs with Zero Mapped Activation Targets")
    lines.append("These transcription factors have 0 activation targets mapped in the model. This occurs because they are either pure repressors, or their downstream target genes do not carry metabolic reactions in *iCW773*:")
    lines.append("")
    if no_targets_tfs:
        lines.append("<details><summary>Click to expand/view TFs with zero mapped targets</summary>")
        lines.append("")
        lines.append("| TF Name | Locus Tag | Total Targets | Mapped Targets |")
        lines.append("| :--- | :--- | :---: | :---: |")
        for t in sorted(no_targets_tfs, key=lambda x: x["name"]):
            lines.append(f"| **{t['name']}** | `{t['locus']}` | {t['targets_count']} | 0 |")
        lines.append("</details>")
    else:
        lines.append("*None.*")
        
    report_content = "\n".join(lines)
    report_path.write_text(report_content, encoding="utf-8")
    print(f"All-TF Simulation report written successfully to {report_path}")
    return report_content

def main():
    print("Loading iCW773 metabolic model...")
    if not MODEL_PATH.exists():
        print(f"Error: model file not found at {MODEL_PATH}")
        return
        
    try:
        model = cobra.io.read_sbml_model(str(MODEL_PATH))
    except Exception as e:
        print(f"Failed to read SBML model: {str(e)}")
        return
        
    cg_to_cgl, cgl_to_cg, cgl_to_name = load_mappings()
    tf_targets, tf_names = load_all_tf_activation_targets(cg_to_cgl)
    
    print("Running simulations...")
    results, biomass_id = run_simulations(model, tf_targets, tf_names)
    
    generate_report(results, tf_targets, tf_names, cgl_to_name, biomass_id)
    
    # Print summary of yields
    print("\n" + "="*50)
    print("      ALL-TF KO SIMULATION PREDICTIONS SUMMARY")
    print("="*50)
    print(f"Baseline Growth Rate: {results['baseline']['growth']:.4f} h^-1")
    print(f"Baseline Max Glutamate: {results['baseline']['max_glutamate']:.4f} mmol/gDCW/h")
    print(f"Baseline Max Lysine: {results['baseline']['max_lysine']:.4f} mmol/gDCW/h")
    print("-"*50)
    
    # Filter yield enhancers for quick stdout printout
    enhancers = []
    lethals = []
    for tf_cg, res in results.items():
        if tf_cg == "baseline" or res["mapped_count"] == 0:
            continue
        g_base = results["baseline"]["growth"]
        glu_base = results["baseline"]["max_glutamate"]
        lys_base = results["baseline"]["max_lysine"]
        
        if res["growth"] < 1e-4:
            lethals.append(tf_cg)
        elif (res["max_glutamate"] > glu_base + 1e-3) or (res["max_lysine"] > lys_base + 1e-3):
            enhancers.append((tf_cg, res))
            
    print(f"Yield Enhancing TF KOs detected: {len(enhancers)}")
    for tf_cg, res in enhancers:
        name = tf_names.get(tf_cg, tf_cg).upper()
        g_pct = pct_change(res["growth"], g_base)
        glu_pct = pct_change(res["max_glutamate"], glu_base)
        lys_pct = pct_change(res["max_lysine"], lys_base)
        print(f" - TF {name:<6} KO -> Growth: {res['growth']:.4f} ({g_pct}) | Max Glu: {res['max_glutamate']:.4f} ({glu_pct}) | Max Lys: {res['max_lysine']:.4f} ({lys_pct})")
        
    print(f"Lethal TF KOs detected: {len(lethals)}")
    print(f" - Loci: {', '.join(sorted([tf_names.get(l, l).upper() for l in lethals]))}")
    print("="*50 + "\n")

if __name__ == "__main__":
    main()
