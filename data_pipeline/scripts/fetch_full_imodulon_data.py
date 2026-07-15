#!/usr/bin/env python3
"""
fetch_full_imodulon_data.py
===========================
Fetches all 87 robust iModulons for C. glutamicum from iModulonDB,
maps gene IDs to cg-locus format, computes regulon overlap statistics
and KEGG pathway enrichments, and saves outputs to the reference directories.
"""

import os
import sys
import json
import csv
import urllib.request
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import run_server

OUT_DIR = ROOT / "data" / "reference" / "imodulon"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# iModulonDB APIs for C. glutamicum modulome263
BASE_URL = "https://imodulondb.org/api/datasets/c_glutamicum/modulome263"
IMODULONS_URL = f"{BASE_URL}/imodulons"
WEIGHTS_URL = f"{BASE_URL}/weights"
GENE_INFO_URL = f"{BASE_URL}/gene-info"
THRESHOLDS_URL = f"{BASE_URL}/thresholds"

def download_json(url):
    print(f"Downloading from {url}...")
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=30) as resp:
        return json.loads(resp.read().decode('utf-8'))

def sanitize_id(name):
    # Sanitize name to make a safe ID
    safe = "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in name)
    return safe.replace("-", "_")

def main():
    # 1. Load mappings from Excel
    print("Loading gene mappings...")
    df_map = pd.read_excel(str(ROOT / 'data' / 'reference' / 'gene_mapping.xlsx'))
    
    # UniProt ID mapping
    uniprot_map = df_map.set_index('Uniprot_protein_ID')['cgb_gene_ID'].dropna().to_dict()

    # Gene name mapping
    name_map = {}
    for _, row in df_map.iterrows():
        name = row.get('Gene_Name')
        cg = row.get('cgb_gene_ID')
        if pd.notna(name) and name != '--' and pd.notna(cg):
            name_map[str(name).lower()] = cg

    # Coordinate mapping (start, end, strand)
    coord_map = {}
    for _, row in df_map.iterrows():
        start = row.get('cgb_start')
        end = row.get('cgb_end')
        strand = row.get('cgb_strand')
        cg = row.get('cgb_gene_ID')
        if pd.notna(start) and start != '--' and pd.notna(end) and end != '--' and pd.notna(strand) and pd.notna(cg):
            try:
                coord_map[(int(start), int(end), strand)] = cg
            except ValueError:
                pass

    def find_overlap_cg(start, end, strand):
        best_cg = None
        best_overlap = 0
        for (s_start, s_end, s_strand), cg in coord_map.items():
            if s_strand != strand:
                continue
            overlap_start = max(start, s_start)
            overlap_end = min(end, s_end)
            if overlap_start < overlap_end:
                overlap = overlap_end - overlap_start
                len1 = end - start
                len2 = s_end - s_start
                min_len = min(len1, len2)
                fraction = overlap / min_len if min_len > 0 else 0
                if fraction > 0.8 and fraction > best_overlap:
                    best_overlap = fraction
                    best_cg = cg
        return best_cg

    # 2. Download raw datasets
    imodulons_raw = download_json(IMODULONS_URL)
    weights_raw = download_json(WEIGHTS_URL)
    gene_info_raw = download_json(GENE_INFO_URL)
    thresholds_raw = download_json(THRESHOLDS_URL)

    # Convert thresholds list to a map
    thresholds_map = {t['imodulon_name']: float(t['threshold']) for t in thresholds_raw}

    # 3. Build gene mapping dictionary (gene_id -> cg_locus)
    gene_to_cg = {}
    for g in gene_info_raw:
        gid = g['gene_id']
        cg = None
        sid = g.get('string_id')
        uniprot = g.get('uniprot_id')
        name = g.get('gene_name')
        start = g.get('start')
        end = g.get('end')
        strand = g.get('strand')
        
        if sid and '.cg' in sid:
            parts = sid.split('.cg')
            if len(parts) > 1:
                cg_num = "".join(c for c in parts[1] if c.isdigit())
                if cg_num:
                    cg = f"cg{cg_num.zfill(4)}"
        
        if not cg and uniprot in uniprot_map:
            cg = uniprot_map[uniprot]
            
        if not cg and name:
            lname = name.lower()
            if lname in name_map:
                cg = name_map[lname]
                
        if not cg and start is not None and end is not None and strand:
            try:
                key = (int(start), int(end), strand)
                if key in coord_map:
                    cg = coord_map[key]
            except ValueError:
                pass

        if not cg and start is not None and end is not None and strand:
            try:
                cg = find_overlap_cg(int(start), int(end), strand)
            except ValueError:
                pass
                
        if cg:
            gene_to_cg[gid] = cg

    # 4. Load KEGG pathway links and names
    print("Loading KEGG pathways and names...")
    run_server.load_organism_kegg_links()
    
    # 5. Load regulations to compare regulons
    print("Loading TF regulations...")
    # Map a TF name/locus to its set of target genes (cg_locus)
    tf_regulons = {}
    reg_path = ROOT / 'data' / 'reference' / 'regulations.csv'
    if reg_path.exists():
        with open(reg_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                tf_locus = row.get('TF_locusTag', '').strip().lower()
                tf_name = row.get('TF_name', '').strip().lower()
                tg_locus = row.get('TG_locusTag', '').strip().lower()
                if tf_locus and tg_locus:
                    tf_regulons.setdefault(tf_locus, set()).add(tg_locus)
                    if tf_name and tf_name != "--":
                        tf_regulons.setdefault(tf_name, set()).add(tg_locus)

    # 6. Process iModulons
    imodulon_gene_weights = {}
    imodulon_metadata = []
    
    # Build pathway metadata details
    all_pathway_genes = set()
    for genes in run_server.PATHWAY_TO_GENES.values():
        for g in genes:
            canonical_g = run_server.CGL_TO_CG.get(g, g).lower()
            all_pathway_genes.add(canonical_g)
    N = len(all_pathway_genes)

    print("Processing 87 iModulons...")
    for im in imodulons_raw:
        k = im['k']
        raw_name = im['name']
        im_id = f"iM_{k}_{sanitize_id(raw_name)}"
        category = im.get('category', 'Uncharacterized')
        var_explained = float(im.get('explained_variance') or im.get('exp_var') or 0.0)
        raw_regulator = im.get('regulator')
        if raw_regulator:
            parts = []
            for p in raw_regulator.split('/'):
                p_clean = p.strip()
                mapped = gene_to_cg.get(p_clean)
                if mapped:
                    gene_name = run_server.GENE_NAMES.get(mapped)
                    parts.append(gene_name if (gene_name and gene_name != mapped) else mapped)
                else:
                    parts.append(p_clean)
            linked_regulator = "/".join(parts)
        else:
            linked_regulator = None
        description = im.get('function_description') or ""
        
        # Check active threshold
        thresh = thresholds_map.get(raw_name, 0.3)
        
        # Collect member genes and weights
        member_genes = {}
        for row in weights_raw:
            gid = row['gene_id']
            cg = gene_to_cg.get(gid)
            if not cg:
                continue
            val = float(row.get(raw_name, 0.0))
            if abs(val) >= thresh:
                member_genes[cg] = val
                
        # Pathway enrichment
        enriched_pathways = []
        if len(member_genes) > 0 and N > 0:
            canonical_members = {g.lower() for g in member_genes.keys()}
            k_size = len(canonical_members.intersection(all_pathway_genes))
            if k_size > 0:
                for pid, pathway_genes in run_server.PATHWAY_TO_GENES.items():
                    canonical_pathway = {run_server.CGL_TO_CG.get(g, g).lower() for g in pathway_genes}
                    hits = canonical_members.intersection(canonical_pathway)
                    x = len(hits)
                    if x > 0:
                        M = len(canonical_pathway)
                        p_val = run_server.hypergeom_sf(x, N, M, k_size)
                        fold_enrichment = (x / k_size) / (M / N) if M > 0 else 0
                        pathway_name = run_server.KEGG_PATHWAY_NAMES.get(pid, pid)
                        
                        enriched_pathways.append({
                            "pathway_id": pid,
                            "pathway_name": pathway_name,
                            "p_value": p_val,
                            "fold_enrichment": fold_enrichment,
                            "hits": sorted(list(hits))
                        })
                # Sort by p-value
                enriched_pathways.sort(key=lambda x: x['p_value'])

        # Regulon overlap metrics
        regulon_overlap = None
        if linked_regulator:
            reg_key = linked_regulator.strip().lower()
            reg_targets = set()
            for k in reg_key.split('/'):
                k_clean = k.strip()
                targets = tf_regulons.get(k_clean)
                if targets:
                    reg_targets.update(targets)
            
            if reg_targets:
                member_set = {g.lower() for g in member_genes.keys()}
                intersection = member_set.intersection(reg_targets)
                overlap_size = len(intersection)
                
                precision = overlap_size / len(member_set) if len(member_set) > 0 else 0.0
                recall = overlap_size / len(reg_targets) if len(reg_targets) > 0 else 0.0
                f1 = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
                
                regulon_overlap = {
                    "regulator": linked_regulator,
                    "regulon_size": len(reg_targets),
                    "overlap_size": overlap_size,
                    "precision": precision,
                    "recall": recall,
                    "f1_score": f1
                }

        # Save to imodulon_gene_weights dict
        imodulon_gene_weights[im_id] = {
            "name": raw_name,
            "linked_regulator": linked_regulator,
            "category": category,
            "variance_explained": var_explained,
            "description": description,
            "threshold": thresh,
            "genes": member_genes,
            "regulon_overlap": regulon_overlap,
            "enriched_pathways": enriched_pathways[:10]  # top 10
        }
        
        # Save to metadata list
        imodulon_metadata.append({
            "id": im_id,
            "k": k,
            "name": raw_name,
            "linked_regulator": linked_regulator,
            "category": category,
            "variance_explained": var_explained,
            "gene_count": len(member_genes),
            "description": description,
            "threshold": thresh,
            "f1_score": regulon_overlap["f1_score"] if regulon_overlap else 0.0,
            "top_pathways": [p['pathway_name'] for p in enriched_pathways[:3]]
        })

    # Sort metadata by variance explained
    imodulon_metadata.sort(key=lambda x: -x['variance_explained'])

    # Build inverted imodulon_by_gene mapping
    imodulon_by_gene = {}
    for im_id, im in imodulon_gene_weights.items():
        for cg in im["genes"]:
            imodulon_by_gene.setdefault(cg, []).append(im_id)

    # 7. Write to files
    print(f"Writing output files to {OUT_DIR}...")
    
    with open(OUT_DIR / "imodulon_gene_weights.json", "w", encoding="utf-8") as f:
        json.dump(imodulon_gene_weights, f, indent=2, ensure_ascii=False)
        
    with open(OUT_DIR / "imodulon_metadata.json", "w", encoding="utf-8") as f:
        json.dump(imodulon_metadata, f, indent=2, ensure_ascii=False)
        
    with open(OUT_DIR / "imodulon_by_gene.json", "w", encoding="utf-8") as f:
        json.dump(imodulon_by_gene, f, indent=2, ensure_ascii=False)

    print(f"Successfully processed {len(imodulon_metadata)} iModulons!")
    print(f"Total mapped unique genes in iModulons: {len(imodulon_by_gene)}")

if __name__ == "__main__":
    main()
