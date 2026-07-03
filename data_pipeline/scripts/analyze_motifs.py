#!/usr/bin/env python3
"""
C. glutamicum Transcriptional Regulatory Network Motif Analyzer
=============================================================
Identifies all Feed-Forward Loops (FFLs) in the TRN, classifies them
into coherent/incoherent types, and isolates those controlling key
industrial metabolic genes (glutamate & lysine pathways).
"""

import os
import csv
from pathlib import Path
from collections import defaultdict

# Setup paths relative to workspace root
ROOT_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT_DIR / "data" / "reference"
REGULATIONS_CSV = DATA_DIR / "regulations.csv"
MAPPING_CSV = DATA_DIR / "gene_mapping.csv"
OUTPUT_DIR = ROOT_DIR / "data_pipeline" / "outputs"

# Target metabolic genes of interest in C. glutamicum
KEY_GENES = {
    # Lysine biosynthesis branch
    "Cgl0251": "lysC",   # Aspartokinase
    "Cgl0252": "asd",    # Aspartate-semialdehyde dehydrogenase
    "Cgl1971": "dapA",   # Dihydrodipicolinate synthase
    "Cgl1973": "dapB",   # Dihydrodipicolinate reductase
    "Cgl1106": "dapD",   # Tetrahydrodipicolinate succinylase
    "Cgl1109": "dapE",   # Succinyl-diaminopimelate desuccinylase
    "Cgl1943": "dapF",   # Diaminopimelate epimerase
    "Cgl1180": "lysA",   # Diaminopimelate decarboxylase
    # Glutamate / Nitrogen assimilation branch
    "Cgl2079": "gdh",    # Glutamate dehydrogenase
    "Cgl0184": "gltB",   # Glutamate synthase (large)
    "Cgl0185": "gltD",   # Glutamate synthase (small)
    "Cgl2214": "glnA",   # Glutamine synthetase
    # TCA Cycle / AKG supply
    "Cgl0829": "gltA",   # Citrate synthase
    "Cgl0664": "icd",    # Isocitrate dehydrogenase
    "Cgl1129": "odhA",   # 2-oxoglutarate dehydrogenase E1
    "Cgl2380": "mdh",    # Malate dehydrogenase
}

def load_mappings():
    """Load gene locus mappings (cg_locus -> cgl_locus -> gene_name)."""
    cg_to_cgl = {}
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
                        cg_to_cgl[cg] = cgl
                        if name:
                            cgl_to_name[cgl] = name
    return cg_to_cgl, cgl_to_name

def load_network(cg_to_cgl):
    """Load regulatory network into adjacency representations."""
    # Adj list: source -> set of (target, role)
    adj = defaultdict(set)
    # Role lookups: (source, target) -> role
    edge_roles = {}
    # Track nodes
    all_nodes = set()
    tfs = set()
    targets = set()
    
    if not REGULATIONS_CSV.exists():
        print(f"Error: regulations file not found at {REGULATIONS_CSV}")
        return None, None, None, None, None
        
    with open(REGULATIONS_CSV, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            tf_cg = row.get("TF_locusTag", "").strip()
            tg_cg = row.get("TG_locusTag", "").strip()
            role = row.get("Role", "").strip() # 'A' (activation), 'R' (repression), 'D' (dual)
            
            if tf_cg and tg_cg:
                # Normalize locus tags to Cgl format if mapping exists
                tf_cgl = cg_to_cgl.get(tf_cg.lower(), tf_cg)
                tg_cgl = cg_to_cgl.get(tg_cg.lower(), tg_cg)
                
                # Exclude self-loops from FFL analysis
                if tf_cgl == tg_cgl:
                    continue
                    
                adj[tf_cgl].add((tg_cgl, role))
                edge_roles[(tf_cgl, tg_cgl)] = role
                
                all_nodes.add(tf_cgl)
                all_nodes.add(tg_cgl)
                tfs.add(tf_cgl)
                targets.add(tg_cgl)
                
    return adj, edge_roles, all_nodes, tfs, targets

def classify_ffl(role_xy, role_yz, role_xz):
    """
    Classify Feed-Forward Loop (FFL) based on the signs of its 3 edges.
    Signs: A/activation = +1, R/repression = -1, D/dual = 0
    Returns FFL category (Coherent/Incoherent/Dual) and type details.
    """
    def get_sign(role):
        if role == 'A': return 1
        if role == 'R': return -1
        return 0 # dual / unknown
        
    s_xy = get_sign(role_xy)
    s_yz = get_sign(role_yz)
    s_xz = get_sign(role_xz)
    
    # If any edge is dual/unknown (0), classify as Dual/Mixed
    if s_xy == 0 or s_yz == 0 or s_xz == 0:
        return "Dual / Mixed FFL", "Contains dual-role regulator(s)"
        
    # Net sign of indirect path = s_xy * s_yz
    indirect_sign = s_xy * s_yz
    
    if indirect_sign == s_xz:
        # Coherent FFL (indirect path has same sign as direct path)
        if s_xy == 1 and s_yz == 1 and s_xz == 1:
            return "Coherent Type 1 (C1-FFL)", "All edges activating (+) -> delays activation, accelerates shut-off"
        elif s_xy == -1 and s_yz == -1 and s_xz == 1:
            return "Coherent Type 2 (C2-FFL)", "Double negative indirect path (-) * (-) = (+) -> activates target"
        elif s_xy == 1 and s_yz == -1 and s_xz == -1:
            return "Coherent Type 3 (C3-FFL)", "Activating + Repressing indirect path (+) * (-) = (-) -> represses target"
        elif s_xy == -1 and s_yz == 1 and s_xz == -1:
            return "Coherent Type 4 (C4-FFL)", "Repressing + Activating indirect path (-) * (+) = (-) -> represses target"
    else:
        # Incoherent FFL (indirect path has opposite sign to direct path)
        if s_xy == 1 and s_yz == -1 and s_xz == 1:
            return "Incoherent Type 1 (I1-FFL)", "Activating direct, repressing indirect (+) * (-) vs (+) -> generates pulses / accelerates response"
        elif s_xy == -1 and s_yz == -1 and s_xz == -1:
            return "Incoherent Type 2 (I2-FFL)", "Double repressing indirect path vs repressing direct (-) * (-) vs (-)"
        elif s_xy == 1 and s_yz == 1 and s_xz == -1:
            return "Incoherent Type 3 (I3-FFL)", "All activating indirect vs repressing direct (+) * (+) vs (-)"
        elif s_xy == -1 and s_yz == 1 and s_xz == 1:
            return "Incoherent Type 4 (I4-FFL)", "Repressing + Activating indirect vs activating direct (-) * (+) vs (+)"
            
    return "Unknown FFL", "Unclassified role combination"

def find_ffls(adj, edge_roles):
    """Find and classify all FFLs (X -> Y, Y -> Z, X -> Z)."""
    ffls = []
    
    # X is the master regulator, Y is the intermediate, Z is the target
    for x in adj:
        # Get all Ys regulated by X
        ys = [y_pair[0] for y_pair in adj[x]]
        
        for y in ys:
            # Check if Y also regulates target genes (Y is a TF)
            if y not in adj:
                continue
                
            # Get all Zs regulated by Y
            zs = [z_pair[0] for z_pair in adj[y]]
            
            for z in zs:
                # X-Y-Z must be distinct nodes
                if z == x or z == y:
                    continue
                    
                # To form an FFL, X must also directly regulate Z
                if (x, z) in edge_roles:
                    role_xy = edge_roles[(x, y)]
                    role_yz = edge_roles[(y, z)]
                    role_xz = edge_roles[(x, z)]
                    
                    category, description = classify_ffl(role_xy, role_yz, role_xz)
                    ffls.append({
                        "X": x,
                        "Y": y,
                        "Z": z,
                        "role_xy": role_xy,
                        "role_yz": role_yz,
                        "role_xz": role_xz,
                        "category": category,
                        "description": description
                    })
                    
    return ffls

def generate_report(ffls, all_nodes, tfs, targets, edge_roles, cgl_to_name):
    """Generate a Markdown report of FFLs and save it."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    report_path = OUTPUT_DIR / "motifs_report.md"
    
    # Count categories
    counts = defaultdict(int)
    for f in ffls:
        counts[f["category"]] += 1
        
    # Find FFLs regulating key metabolic genes
    metabolic_ffls = []
    for f in ffls:
        z = f["Z"]
        if z in KEY_GENES:
            metabolic_ffls.append(f)
            
    # Format node helper
    def node_fmt(locus):
        name = cgl_to_name.get(locus, locus)
        if name != locus:
            return f"**{name}** ({locus})"
        return f"`{locus}`"
        
    lines = []
    lines.append("# Transcriptional Regulatory Network Motif Report: C. glutamicum")
    lines.append("")
    lines.append(f"This report presents a topological analysis of **Feed-Forward Loops (FFLs)** within the transcriptional regulatory network of *Corynebacterium glutamicum*, compiled from literature and curated databases.")
    lines.append("")
    
    lines.append("## 1. Network Statistics")
    lines.append("")
    lines.append(f"- **Total Regulatory Nodes**: {len(all_nodes)}")
    lines.append(f"- **Transcription Factors (TFs)**: {len(tfs)}")
    lines.append(f"- **Target Genes**: {len(targets)}")
    lines.append(f"- **Total Regulatory Relationships**: {len(edge_roles)}")
    lines.append(f"- **Total Feed-Forward Loops (FFLs) Detected**: {len(ffls)}")
    lines.append("")
    
    lines.append("## 2. FFL Motif Distribution")
    lines.append("")
    lines.append("| FFL Category | Count | Description |")
    lines.append("| :--- | :---: | :--- |")
    
    sorted_categories = sorted(counts.items(), key=lambda pair: pair[1], reverse=True)
    for cat, count in sorted_categories:
        # Find description from the first match
        desc = next(f["description"] for f in ffls if f["category"] == cat)
        lines.append(f"| {cat} | {count} | {desc} |")
    lines.append("")
    
    lines.append("## 3. Key Metabolic Pathway FFLs (Glutamate, Lysine, TCA)")
    lines.append("")
    lines.append(f"Detected **{len(metabolic_ffls)}** FFLs regulating key genes in lysine biosynthesis, glutamate synthesis, and TCA cycle intermediate supplies. These loops act as metabolic tuning knobs under environmental stresses.")
    lines.append("")
    
    if metabolic_ffls:
        lines.append("| Target Gene | Role | Master TF (X) | Intermediate TF (Y) | FFL Type | Loop Structure |")
        lines.append("| :--- | :--- | :--- | :--- | :--- | :--- |")
        
        # Sort by target gene name
        metabolic_ffls_sorted = sorted(metabolic_ffls, key=lambda f: KEY_GENES[f["Z"]])
        
        for f in metabolic_ffls_sorted:
            z_locus = f["Z"]
            z_func = KEY_GENES[z_locus]
            z_name = cgl_to_name.get(z_locus, z_locus)
            
            x_name = cgl_to_name.get(f["X"], f["X"])
            y_name = cgl_to_name.get(f["Y"], f["Y"])
            
            role_str = f"X --({f['role_xy']})--> Y, Y --({f['role_yz']})--> Z, X --({f['role_xz']})--> Z"
            
            lines.append(f"| {z_name} (`{z_locus}`) | {z_func} | {x_name} | {y_name} | {f['category'].split(' (')[0]} | `{role_str}` |")
    else:
        lines.append("*No FFLs regulating the specified key genes were found in the current regulations dataset.*")
    lines.append("")
    
    lines.append("## 4. Biological Interpretations")
    lines.append("")
    lines.append("> [!NOTE]")
    lines.append("> **Coherent Type 1 (C1-FFL)** is the most common motif. It acts as a **sign-sensitive delay filter**, ensuring that a target gene (like biosynthetic enzymes) is only expressed when the master signal is sustained. Transient fluctuations will not trigger expression, saving protein synthesis resources.")
    lines.append("")
    lines.append("> [!TIP]")
    lines.append("> **Incoherent Type 1 (I1-FFL)** acts as an **accelerator and pulse generator**. When the master regulator X turns on, it quickly activates target Z. Shortly after, it activates Y (repressor), which shuts down Z. This creates a quick pulse of gene expression and speeds up cell response to rapid changes.")
    
    report_content = "\n".join(lines)
    report_path.write_text(report_content, encoding="utf-8")
    print(f"Motifs report written successfully to {report_path}")
    return report_content

def main():
    print("Starting network motif analysis...")
    cg_to_cgl, cgl_to_name = load_mappings()
    adj, edge_roles, all_nodes, tfs, targets = load_network(cg_to_cgl)
    
    if adj is None:
        return
        
    ffls = find_ffls(adj, edge_roles)
    report = generate_report(ffls, all_nodes, tfs, targets, edge_roles, cgl_to_name)
    
    # Print summary to stdout
    print("\n" + "="*40)
    print("          FFL MOTIF ANALYSIS RESULTS")
    print("="*40)
    print(f"Total regulatory nodes: {len(all_nodes)}")
    print(f"Total regulatory edges: {len(edge_roles)}")
    print(f"Total FFLs found: {len(ffls)}")
    print("-"*40)
    
    counts = defaultdict(int)
    for f in ffls:
        counts[f["category"]] += 1
    for cat, count in sorted(counts.items(), key=lambda pair: pair[1], reverse=True):
        print(f" - {cat}: {count}")
    print("="*40 + "\n")

if __name__ == "__main__":
    main()
