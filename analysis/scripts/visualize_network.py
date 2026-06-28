#!/usr/bin/env python3
import os
import sys
import argparse

# Check requirements and prompt user if missing
try:
    import pandas as pd
    import networkx as nx
    import matplotlib.pyplot as plt
except ImportError:
    print("Warning: Missing required libraries. Installing pandas, networkx, matplotlib...")
    import subprocess
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pandas", "networkx", "matplotlib"])
        import pandas as pd
        import networkx as nx
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"Error installing dependencies: {e}")
        print("Please install them manually: pip install pandas networkx matplotlib pyvis")
        sys.exit(1)

# Pyvis is optional for interactive HTML plots
pyvis_available = False
try:
    from pyvis.network import Network
    pyvis_available = True
except ImportError:
    print("Tip: Install 'pyvis' to generate interactive HTML network views (pip install pyvis).")

# Data directories - resolved relative to this script's directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(ROOT_DIR, "data")
REGULATIONS_FILE = os.path.join(DATA_DIR, "regulations.csv")
RNA_REG_FILE = os.path.join(DATA_DIR, "rna_regulation.csv")
MAPPING_FILE = os.path.join(DATA_DIR, "gene_mapping.csv")

def load_data():
    if not os.path.exists(REGULATIONS_FILE):
        print(f"Error: {REGULATIONS_FILE} not found. Looked in {REGULATIONS_FILE}.")
        sys.exit(1)
        
    print("Loading regulatory data...")
    # Load TF-TG regulations
    df_reg = pd.read_csv(REGULATIONS_FILE)
    
    # Load sRNA-mRNA regulations if file exists
    df_rna = None
    if os.path.exists(RNA_REG_FILE):
        try:
            df_rna = pd.read_csv(RNA_REG_FILE, sep='\t')
            # Filter sRNA-mRNA: keep high confidence ones (e.g. rank <= 15 or p-value < 0.05)
            # You can customize this threshold
            df_rna = df_rna[df_rna['rank'] <= 10].copy()
        except Exception as e:
            print(f"Warning: Failed to load sRNA data: {e}")
            df_rna = None

    # Load gene mapping if file exists
    df_map = None
    if os.path.exists(MAPPING_FILE):
        try:
            df_map = pd.read_csv(MAPPING_FILE)
        except Exception as e:
            print(f"Warning: Failed to load gene mapping data: {e}")
            df_map = None
            
    return df_reg, df_rna, df_map

def build_network(gene_query, df_reg, df_rna, df_map, steps=1, only_coregulated=False):
    gene_query = gene_query.strip().lower()
    
    # Create directed graph
    G = nx.DiGraph()
    
    # Helper to clean strings
    def clean_str(val):
        return str(val).strip() if pd.notna(val) else ""

    # Resolve maps
    cgl_to_cg = {}
    cg_to_cgl = {}
    name_to_cg = {}
    cg_to_product = {}
    if df_map is not None:
        for _, row in df_map.iterrows():
            cgl = clean_str(row.get('cgl_locus'))
            cg = clean_str(row.get('cg_locus'))
            name = clean_str(row.get('gene_name'))
            product = clean_str(row.get('product'))
            if cgl and cg and cgl.lower() != 'nan' and cg.lower() != 'nan':
                cgl_to_cg[cgl.lower()] = cg
                cg_to_cgl[cg.lower()] = cgl
            if name and name != '--' and name.lower() != 'nan' and cg and cg.lower() != 'nan':
                name_to_cg[name.lower()] = cg
            if cg and product:
                cg_to_product[cg.lower()] = product
            if cgl and product:
                cg_to_product[cgl.lower()] = product

    # Store mapping in Graph metadata
    G.graph['cg_to_cgl'] = cg_to_cgl
    G.graph['cg_to_product'] = cg_to_product

    # Map genes to their details (name, role)
    gene_metadata = {}
    
    # Pre-process all genes for matching (case-insensitive)
    # We want to match TF_locusTag, TF_name, TG_locusTag, TG_name, or srna, mrna
    # Find all regulations involving the queried gene
    visited_nodes = set()
    
    # Let's collect a mapping of lowercase tags/names to their standard names
    all_genes = set()
    for _, row in df_reg.iterrows():
        tf_tag, tf_name = clean_str(row['TF_locusTag']), clean_str(row['TF_name'])
        tg_tag, tg_name = clean_str(row['TG_locusTag']), clean_str(row['TG_name'])
        
        all_genes.add((tf_tag, tf_name))
        all_genes.add((tg_tag, tg_name))
        
    if df_rna is not None:
        for _, row in df_rna.iterrows():
            srna, mrna = clean_str(row['srna']), clean_str(row['mrna'])
            all_genes.add((srna, srna))
            all_genes.add((mrna, mrna))
            
    # Resolve lowercase to standard casing
    casing_map = {}
    for tag, name in all_genes:
        if tag:
            casing_map[tag.lower()] = tag
            gene_metadata[tag] = {"name": name if name else tag, "type": "gene"}
        if name:
            casing_map[name.lower()] = name
            gene_metadata[name] = {"name": name, "type": "gene"}

    # Resolve all queries (supporting list/comma-separated strings)
    if isinstance(gene_query, str):
        queries = [g.strip().lower() for g in gene_query.split(",") if g.strip()]
    else:
        queries = [g.strip().lower() for g in gene_query]

    resolved_queries = []
    for q in queries:
        resolved = None
        if q in cgl_to_cg:
            resolved = casing_map.get(cgl_to_cg[q].lower(), cgl_to_cg[q])
        elif q in name_to_cg:
            resolved = casing_map.get(name_to_cg[q].lower(), name_to_cg[q])
        elif q in casing_map:
            resolved = casing_map[q]
            
        if resolved:
            resolved_queries.append(resolved)
        else:
            print(f"Warning: Gene '{q}' not found in dataset. Checking for substring matches...")
            matches = [val for key, val in casing_map.items() if q in key]
            if not matches and df_map is not None:
                matches = [cgl_to_cg[key] for key in cgl_to_cg if q in key]
                matches = [casing_map.get(m.lower(), m) for m in matches if m.lower() in casing_map or m in casing_map]
            if matches:
                resolved_queries.append(matches[0])
                print(f"Using match: {matches[0]}")
            else:
                print(f"No matches found for '{q}'. Skipping.")

    if not resolved_queries:
        print("No valid genes found to visualize.")
        return None, None

    query_set = set(resolved_queries)
    
    # We start traversal from all resolved queries
    queue = [(g, 0) for g in resolved_queries]
    
    # Helper to get cgl-prioritized label
    def get_node_label(tag, name):
        if not tag:
            return name if name else ""
        tag_lower = tag.lower()
        if tag_lower in cg_to_cgl:
            return cg_to_cgl[tag_lower]
        if name and name != '--' and name.lower() != 'nan':
            return name
        return tag

    # Add queried nodes metadata
    for g in resolved_queries:
        g_name = gene_metadata.get(g, {}).get("name", g)
        G.add_node(g, label=get_node_label(g, g_name), type="query")

    while queue:
        curr, d = queue.pop(0)
        if curr in visited_nodes or d >= steps:
            continue
        visited_nodes.add(curr)
        
        curr_lower = curr.lower()
        
        # 1. Look for TF regulations
        # Case A: curr is TF
        tf_matches = df_reg[
            (df_reg['TF_locusTag'].str.lower() == curr_lower) | 
            (df_reg['TF_name'].str.lower() == curr_lower)
        ]
        for _, row in tf_matches.iterrows():
            tf_tag = casing_map.get(clean_str(row['TF_locusTag']).lower(), clean_str(row['TF_locusTag']))
            tg_tag = casing_map.get(clean_str(row['TG_locusTag']).lower(), clean_str(row['TG_locusTag']))
            role = clean_str(row['Role'])
            
            # Label
            tf_label = get_node_label(tf_tag, gene_metadata.get(tf_tag, {}).get("name"))
            tg_label = get_node_label(tg_tag, gene_metadata.get(tg_tag, {}).get("name"))
            
            G.add_node(tf_tag, label=tf_label, type="TF" if tf_tag not in query_set else "query")
            G.add_node(tg_tag, label=tg_label, type="gene" if tg_tag not in query_set else "query")
            G.add_edge(tf_tag, tg_tag, role=role, type="TF-TG")
            
            if tg_tag not in visited_nodes and d + 1 < steps:
                queue.append((tg_tag, d + 1))
                
        # Case B: curr is Target Gene (TG)
        tg_matches = df_reg[
            (df_reg['TG_locusTag'].str.lower() == curr_lower) | 
            (df_reg['TG_name'].str.lower() == curr_lower)
        ]
        for _, row in tg_matches.iterrows():
            tf_tag = casing_map.get(clean_str(row['TF_locusTag']).lower(), clean_str(row['TF_locusTag']))
            tg_tag = casing_map.get(clean_str(row['TG_locusTag']).lower(), clean_str(row['TG_locusTag']))
            role = clean_str(row['Role'])
            
            tf_label = get_node_label(tf_tag, gene_metadata.get(tf_tag, {}).get("name"))
            tg_label = get_node_label(tg_tag, gene_metadata.get(tg_tag, {}).get("name"))
            
            G.add_node(tf_tag, label=tf_label, type="TF" if tf_tag not in query_set else "query")
            G.add_node(tg_tag, label=tg_label, type="gene" if tg_tag not in query_set else "query")
            G.add_edge(tf_tag, tg_tag, role=role, type="TF-TG")
            
            if tf_tag not in visited_nodes and d + 1 < steps:
                queue.append((tf_tag, d + 1))
 
        # 2. Look for sRNA regulations
        if df_rna is not None:
            # Case A: curr is sRNA
            srna_matches = df_rna[df_rna['srna'].str.lower() == curr_lower]
            for _, row in srna_matches.iterrows():
                srna = casing_map.get(clean_str(row['srna']).lower(), clean_str(row['srna']))
                mrna = casing_map.get(clean_str(row['mrna']).lower(), clean_str(row['mrna']))
                energy = float(row.get('energy', 0))
                pvalue = float(row.get('copra_pvalue', 0))
                
                srna_label = get_node_label(srna, gene_metadata.get(srna, {}).get("name"))
                mrna_label = get_node_label(mrna, gene_metadata.get(mrna, {}).get("name"))
                
                G.add_node(srna, label=srna_label, type="sRNA" if srna not in query_set else "query")
                G.add_node(mrna, label=mrna_label, type="gene" if mrna not in query_set else "query")
                G.add_edge(srna, mrna, role="sRNA", type="sRNA-mRNA", energy=energy, pvalue=pvalue)
                
                if mrna not in visited_nodes and d + 1 < steps:
                    queue.append((mrna, d + 1))
                    
            # Case B: curr is mRNA (target of sRNA)
            mrna_matches = df_rna[df_rna['mrna'].str.lower() == curr_lower]
            for _, row in mrna_matches.iterrows():
                srna = casing_map.get(clean_str(row['srna']).lower(), clean_str(row['srna']))
                mrna = casing_map.get(clean_str(row['mrna']).lower(), clean_str(row['mrna']))
                energy = float(row.get('energy', 0))
                pvalue = float(row.get('copra_pvalue', 0))
                
                srna_label = get_node_label(srna, gene_metadata.get(srna, {}).get("name"))
                mrna_label = get_node_label(mrna, gene_metadata.get(mrna, {}).get("name"))
                
                G.add_node(srna, label=srna_label, type="sRNA" if srna not in query_set else "query")
                G.add_node(mrna, label=mrna_label, type="gene" if mrna not in query_set else "query")
                G.add_edge(srna, mrna, role="sRNA", type="sRNA-mRNA", energy=energy, pvalue=pvalue)
                
                if srna not in visited_nodes and d + 1 < steps:
                    queue.append((srna, d + 1))
 
    if only_coregulated:
        # Filter G to only co-regulated target genes and their connections
        target_nodes = [node for node in G.nodes() if G.nodes[node].get("type") == "gene"]
        coreg_targets = {node for node in target_nodes if G.in_degree(node) >= 2}
        
        # Nodes to keep
        nodes_to_keep = set(resolved_queries) | coreg_targets
        edges_to_keep = []
        for u, v in G.edges():
            v_type = G.nodes[v].get("type", "")
            if v_type == "gene":
                if v in coreg_targets:
                    edges_to_keep.append((u, v))
                    nodes_to_keep.add(u)
                    nodes_to_keep.add(v)
            else:
                edges_to_keep.append((u, v))
                nodes_to_keep.add(u)
                nodes_to_keep.add(v)
                
        # Rebuild filtered graph
        G_filtered = nx.DiGraph()
        G_filtered.graph.update(G.graph)
        
        for node in nodes_to_keep:
            if node in G:
                G_filtered.add_node(node, **G.nodes[node])
                
        for u, v in edges_to_keep:
            if u in G_filtered and v in G_filtered:
                G_filtered.add_edge(u, v, **G[u][v])
                
        G = G_filtered

    resolved_query = "_".join(resolved_queries) if len(resolved_queries) <= 3 else "multi_gene"
    return G, resolved_query

def draw_matplotlib(G, query_gene, output_file):
    if G.number_of_nodes() == 0:
        print("No nodes in the graph to draw.")
        return
        
    # Standard academic figure size and white facecolor
    fig = plt.figure(figsize=(10, 8), facecolor='white')
    ax = fig.add_subplot(111)
    ax.set_facecolor('white')
    
    # Layout - spring layout with customized spacing
    pos = nx.spring_layout(G, k=0.6, iterations=50)
    
    # Academic Colors mapping & Borders
    node_colors = []
    node_edge_colors = []
    for node in G.nodes():
        node_type = G.nodes[node].get("type", "gene")
        if node_type == "query":
            node_colors.append("#ffe0b2")       # Muted soft orange background
            node_edge_colors.append("#f57c00")  # Dark orange border
        elif node_type == "TF":
            node_colors.append("#e3f2fd")       # Muted soft blue background
            node_edge_colors.append("#1976d2")  # Dark blue border
        elif node_type == "sRNA":
            node_colors.append("#f3e5f5")       # Muted soft purple background
            node_edge_colors.append("#8e24aa")  # Dark purple border
        else:
            if G.in_degree(node) > 1:
                node_colors.append("#e0f2f1")       # Soft Teal/Mint background for shared target
                node_edge_colors.append("#00897b")  # Dark Teal border
            else:
                node_colors.append("#f5f5f5")       # Muted soft gray background
                node_edge_colors.append("#757575")  # Dark gray border
            
    # Node sizes dictionary mapping for better layout calculations
    node_sizes_dict = {node: (1200 if node == query_gene or G.nodes[node].get("type") == "query" else 650) for node in G.nodes()}
    node_sizes = [node_sizes_dict[node] for node in G.nodes()]
    
    # Edges styling for Academic theme
    edge_colors = []
    edge_styles = []
    for u, v in G.edges():
        role = G[u][v].get("role", "")
        edge_type = G[u][v].get("type", "")
        
        if edge_type == "sRNA-mRNA":
            edge_colors.append("#7b1fa2")  # Academic purple
            edge_styles.append("dashed")
        elif role == "A":
            edge_colors.append("#2e7d32")  # Academic green
            edge_styles.append("solid")
        elif role == "R":
            edge_colors.append("#d32f2f")  # Academic red
            edge_styles.append("solid")
        else:  # Dual or Sigma or other
            edge_colors.append("#e65100")  # Academic dark orange
            edge_styles.append("solid")

    # Draw Nodes with outlines
    nx.draw_networkx_nodes(
        G, pos, 
        node_color=node_colors, 
        edgecolors=node_edge_colors, 
        linewidths=1.5,
        node_size=node_sizes, 
        alpha=1.0,
        ax=ax
    )
    
    # Draw labels - black, crisp, sans-serif
    labels = {node: G.nodes[node].get("label", node) for node in G.nodes()}
    nx.draw_networkx_labels(
        G, pos, labels, 
        font_size=8.5, 
        font_family="sans-serif", 
        font_color="#000000",
        font_weight="normal",
        ax=ax
    )
    
    # Draw edges with correct clipping using node_size dict to avoid overlapping borders
    for idx, (u, v) in enumerate(G.edges()):
        nx.draw_networkx_edges(
            G, pos, edgelist=[(u, v)],
            edge_color=edge_colors[idx],
            style=edge_styles[idx],
            width=1.5,
            arrows=True,
            arrowsize=12,
            node_size=node_sizes,
            ax=ax
        )
        
    # Academic styled Legend
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    legend_elements = [
        Patch(facecolor='#ffe0b2', edgecolor='#f57c00', linewidth=1.5, label='Query Gene'),
        Patch(facecolor='#e3f2fd', edgecolor='#1976d2', linewidth=1.5, label='Transcription Factor (TF)'),
        Patch(facecolor='#f3e5f5', edgecolor='#8e24aa', linewidth=1.5, label='sRNA'),
        Patch(facecolor='#f5f5f5', edgecolor='#757575', linewidth=1.5, label='Target Gene'),
        Patch(facecolor='#e0f2f1', edgecolor='#00897b', linewidth=1.5, label='Shared Target (co-regulated)'),
        Line2D([0], [0], color='#2e7d32', lw=2, label='Activation (+)'),
        Line2D([0], [0], color='#d32f2f', lw=2, label='Repression (-)'),
        Line2D([0], [0], color='#e65100', lw=2, label='Dual/Sigma regulation'),
        Line2D([0], [0], color='#7b1fa2', lw=2, linestyle='--', label='sRNA regulation')
    ]
    ax.legend(
        handles=legend_elements, 
        loc='upper right', 
        fontsize=8, 
        frameon=True, 
        facecolor='white', 
        edgecolor='#e2e8f0'
    )
    
    plt.title(f"Regulatory Network (Academic Theme) centered on {query_gene}", fontsize=12, color="#0f172a", pad=15)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, facecolor='white')
    plt.close()
    print(f"Static plot saved as {output_file}")

def draw_pyvis(G, query_gene, output_file):
    if not pyvis_available:
        return
        
    # Academic white background and black fonts for Pyvis
    net = Network(height="750px", width="100%", bgcolor="#ffffff", font_color="#000000", directed=True)
    
    # Configuration options for Pyvis
    net.set_options("""
    var options = {
      "nodes": {
        "borderWidth": 2,
        "borderWidthSelected": 4,
        "size": 25,
        "font": {
          "size": 13,
          "face": "arial",
          "color": "#000000"
        }
      },
      "edges": {
        "color": {
          "inherit": false
        },
        "smooth": {
          "type": "continuous",
          "forceDirection": "none"
        }
      },
      "physics": {
        "forceAtlas2Based": {
          "gravitationalConstant": -50,
          "centralGravity": 0.01,
          "springLength": 100,
          "springConstant": 0.08
        },
        "maxVelocity": 50,
        "solver": "forceAtlas2Based",
        "timestep": 0.35,
        "stabilization": {
          "iterations": 150
        }
      }
    }
    """)
    
    # Add nodes with academic style
    for node in G.nodes():
        node_type = G.nodes[node].get("type", "gene")
        label = G.nodes[node].get("label", node)
        
        # Set node visual properties
        cg_to_cgl = G.graph.get('cg_to_cgl', {})
        cg_to_product = G.graph.get('cg_to_product', {})
        
        cgl_tag = cg_to_cgl.get(node.lower(), "N/A")
        product = cg_to_product.get(node.lower(), "")
        if not product and cgl_tag != "N/A":
            product = cg_to_product.get(cgl_tag.lower(), "")
            
        title = f"Locus Tag: {node}<br>Cgl Locus Tag: {cgl_tag}<br>Name: {label}<br>Type: {node_type}"
        if product:
            title += f"<br>Function: {product}"
        
        if node_type == "query":
            color = {"background": "#ffe0b2", "border": "#f57c00", "highlight": {"background": "#ffe0b2", "border": "#f57c00"}}
            size = 35
        elif node_type == "TF":
            color = {"background": "#e3f2fd", "border": "#1976d2", "highlight": {"background": "#e3f2fd", "border": "#1976d2"}}
            size = 28
        elif node_type == "sRNA":
            color = {"background": "#f3e5f5", "border": "#8e24aa", "highlight": {"background": "#f3e5f5", "border": "#8e24aa"}}
            size = 28
        else:
            if G.in_degree(node) > 1:
                # Shared target node
                color = {"background": "#e0f2f1", "border": "#00897b", "highlight": {"background": "#e0f2f1", "border": "#00897b"}}
                size = 24
            else:
                color = {"background": "#f5f5f5", "border": "#757575", "highlight": {"background": "#f5f5f5", "border": "#757575"}}
                size = 22
            
        net.add_node(node, label=label, title=title, color=color, size=size)
        
    # Add edges with academic color palette
    for u, v in G.edges():
        role = G[u][v].get("role", "")
        edge_type = G[u][v].get("type", "")
        
        if edge_type == "sRNA-mRNA":
            color = "#7b1fa2"  # Purple
            label = "sRNA prediction"
            dashes = True
        elif role == "A":
            color = "#2e7d32"  # Green
            label = "Activation"
            dashes = False
        elif role == "R":
            color = "#d32f2f"  # Red
            label = "Repression"
            dashes = False
        else:
            color = "#e65100"  # Dark orange
            label = "Dual/Sigma"
            dashes = False
            
        title = f"Type: {edge_type}<br>Role: {label}"
        if "energy" in G[u][v]:
            title += f"<br>Energy: {G[u][v]['energy']} kcal/mol<br>p-value: {G[u][v]['pvalue']}"
            
        net.add_edge(u, v, color=color, title=title, dashes=dashes, width=2)
        
    net.save_graph(output_file)
    print(f"Interactive interactive plot saved as {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize Corynebacterium glutamicum Regulatory Networks")
    parser.add_argument("--gene", type=str, help="Gene Locus Tag or Name (e.g., cg0350, whiB4)")
    parser.add_argument("--steps", type=int, default=1, help="Steps/degree of neighborhood connections (default: 1)")
    parser.add_argument("--output-img", type=str, help="Path to save Matplotlib PNG (default: <gene>_network.png)")
    parser.add_argument("--output-html", type=str, help="Path to save interactive HTML (default: <gene>_network.html)")
    parser.add_argument("--only-coregulated", action="store_true", help="Only show target genes co-regulated by 2 or more TFs/sRNAs")
    
    args = parser.parse_args()
    
    # Prompt for gene if not provided
    gene = args.gene
    if not gene:
        gene = input("Enter Gene Locus Tag or Name (e.g. cg0012, whiB4, sigH, cspA): ").strip()
        if not gene:
            print("No gene specified. Exiting.")
            sys.exit(0)
            
    df_reg, df_rna, df_map = load_data()
    G, resolved_query = build_network(gene, df_reg, df_rna, df_map, steps=args.steps, only_coregulated=args.only_coregulated)
    
    if G is not None and G.number_of_nodes() > 0:
        print(f"Built network with {G.number_of_nodes()} nodes and {G.number_of_edges()} interactions centered on '{resolved_query}'.")
        
        # Save plots to outputs/ directory
        SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
        ROOT_DIR = os.path.dirname(SCRIPT_DIR)
        OUTPUT_DIR = os.path.join(ROOT_DIR, "outputs")
        os.makedirs(OUTPUT_DIR, exist_ok=True)

        img_out = args.output_img if args.output_img else os.path.join(OUTPUT_DIR, f"{resolved_query}_network.png")
        draw_matplotlib(G, resolved_query, img_out)
        
        if pyvis_available:
            html_out = args.output_html if args.output_html else os.path.join(OUTPUT_DIR, f"{resolved_query}_network_interactive.html")
            draw_pyvis(G, resolved_query, html_out)
    else:
        print("Failed to build network. Please check the gene spelling.")
