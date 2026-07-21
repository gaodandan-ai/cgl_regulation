import os
import xml.etree.ElementTree as ET
import urllib.request
import json
import re

# Resolve directories
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
MODEL_PATH = os.path.join(ROOT_DIR, "backend", "models", "iCW773.xml")
RHEA_OUT = os.path.join(ROOT_DIR, "data", "reference", "rhea_mappings.json")
CHEBI_OUT = os.path.join(ROOT_DIR, "data", "reference", "chebi_mappings.json")

def clean_id(raw_id):
    cleaned = raw_id.upper()
    if cleaned.startswith("R_") or cleaned.startswith("M_"):
        cleaned = cleaned[2:]
    # Remove compartment suffixes
    if cleaned.endswith("_C") or cleaned.endswith("_E"):
        cleaned = cleaned[:-2]
    return cleaned

def parse_model_ids():
    print("Parsing model IDs...")
    tree = ET.parse(MODEL_PATH)
    root = tree.getroot()
    core_ns = "http://www.sbml.org/sbml/level3/version1/core"
    
    metabolite_ids = set()
    reaction_ids = set()
    
    metabolite_map = {} # cleaned_id -> list of model_ids
    reaction_map = {}   # cleaned_id -> list of model_ids

    # Find species
    for species in root.findall(f".//{{{core_ns}}}species"):
        mid = species.attrib.get("id", "")
        if mid:
            cleaned = clean_id(mid)
            metabolite_ids.add(cleaned)
            metabolite_map.setdefault(cleaned, []).append(mid)
            
    # Find reactions
    for reaction in root.findall(f".//{{{core_ns}}}reaction"):
        rid = reaction.attrib.get("id", "")
        if rid:
            cleaned = clean_id(rid)
            reaction_ids.add(cleaned)
            reaction_map.setdefault(cleaned, []).append(rid)
            
    print(f"Parsed {len(metabolite_ids)} metabolites and {len(reaction_ids)} reactions.")
    return metabolite_ids, metabolite_map, reaction_ids, reaction_map

def parse_links(links_str):
    links = {}
    if not links_str or links_str == "\\N":
        return links
    
    parts = re.split(r';\s*', links_str)
    for part in parts:
        if ":" not in part:
            continue
        reg_name, url = part.split(":", 1)
        reg_name = reg_name.strip()
        url = url.strip()
        
        target_id = ""
        if "identifiers.org" in url:
            target_id = url.split("/")[-1]
        else:
            target_id = url
            
        links.setdefault(reg_name, []).append(target_id)
    return links

def build_metabolite_mappings(metabolite_ids, metabolite_map):
    print("Downloading and parsing BiGG metabolites namespace...")
    url = "http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt"
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    
    chebi_mappings = {}
    
    with urllib.request.urlopen(req) as resp:
        header = resp.readline().decode('utf-8').strip().split('\t')
        bigg_idx = header.index("bigg_id")
        univ_idx = header.index("universal_bigg_id")
        name_idx = header.index("name")
        links_idx = header.index("database_links")
        
        for line_bytes in resp:
            line = line_bytes.decode('utf-8').strip('\r\n')
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) <= max(bigg_idx, univ_idx, links_idx):
                continue
                
            bigg_id = cols[bigg_idx]
            univ_id = cols[univ_idx]
            links_str = cols[links_idx]
            name = cols[name_idx] if len(cols) > name_idx else ""
            
            # Clean for matching
            c_bigg = clean_id(bigg_id)
            c_univ = clean_id(univ_id)
            
            matched_key = None
            if c_univ in metabolite_ids:
                matched_key = c_univ
            elif c_bigg in metabolite_ids:
                matched_key = c_bigg
                    
            if matched_key:
                parsed = parse_links(links_str)
                chebi_list = parsed.get("CHEBI", [])
                kegg_list = parsed.get("KEGG Compound", [])
                biocyc_list = parsed.get("BioCyc", [])
                
                for model_id in metabolite_map[matched_key]:
                    mapping = {
                        "name": name,
                        "chebi": chebi_list[0] if chebi_list else None,
                        "chebi_all": chebi_list,
                        "kegg": kegg_list[0] if kegg_list else None,
                        "biocyc": biocyc_list[0] if biocyc_list else None
                    }
                    if mapping["chebi"] or mapping["kegg"] or mapping["biocyc"]:
                        # Merge if already exists (keep first non-null values)
                        if model_id in chebi_mappings:
                            existing = chebi_mappings[model_id]
                            for k in ["chebi", "kegg", "biocyc"]:
                                if not existing[k]:
                                    existing[k] = mapping[k]
                            # Union chebi_all
                            existing["chebi_all"] = list(set(existing["chebi_all"] + mapping["chebi_all"]))
                        else:
                            chebi_mappings[model_id] = mapping
                        
    print(f"Mapped {len(chebi_mappings)} metabolites.")
    return chebi_mappings

def build_reaction_mappings(reaction_ids, reaction_map):
    print("Downloading and parsing BiGG reactions namespace...")
    url = "http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt"
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    
    rhea_mappings = {}
    
    with urllib.request.urlopen(req) as resp:
        header = resp.readline().decode('utf-8').strip().split('\t')
        bigg_idx = header.index("bigg_id")
        name_idx = header.index("name")
        links_idx = header.index("database_links")
        
        for line_bytes in resp:
            line = line_bytes.decode('utf-8').strip('\r\n')
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) <= max(bigg_idx, links_idx):
                continue
                
            bigg_id = cols[bigg_idx]
            links_str = cols[links_idx]
            name = cols[name_idx] if len(cols) > name_idx else ""
            
            c_bigg = clean_id(bigg_id)
            
            if c_bigg in reaction_ids:
                parsed = parse_links(links_str)
                rhea_raw = parsed.get("RHEA", [])
                
                rhea_ids = []
                for r in rhea_raw:
                    m = re.search(r'\d+', r)
                    if m:
                        rhea_ids.append(m.group(0))
                        
                kegg_list = parsed.get("KEGG Reaction", [])
                biocyc_list = parsed.get("BioCyc", [])
                ec_list = parsed.get("EC Number", [])
                
                for model_id in reaction_map[c_bigg]:
                    mapping = {
                        "name": name,
                        "rhea": rhea_ids[0] if rhea_ids else None,
                        "rhea_all": rhea_ids,
                        "kegg": kegg_list[0] if kegg_list else None,
                        "biocyc": biocyc_list[0] if biocyc_list else None,
                        "ec": ec_list[0] if ec_list else None
                    }
                    if mapping["rhea"] or mapping["kegg"] or mapping["biocyc"] or mapping["ec"]:
                        if model_id in rhea_mappings:
                            existing = rhea_mappings[model_id]
                            for k in ["rhea", "kegg", "biocyc", "ec"]:
                                if not existing[k]:
                                    existing[k] = mapping[k]
                            existing["rhea_all"] = list(set(existing["rhea_all"] + mapping["rhea_all"]))
                        else:
                            rhea_mappings[model_id] = mapping
                        
    print(f"Mapped {len(rhea_mappings)} reactions.")
    return rhea_mappings

def main():
    metabolite_ids, metabolite_map, reaction_ids, reaction_map = parse_model_ids()
    
    # Generate mappings
    chebi_mappings = build_metabolite_mappings(metabolite_ids, metabolite_map)
    rhea_mappings = build_reaction_mappings(reaction_ids, reaction_map)
    
    # Create directory if not exists
    os.makedirs(os.path.dirname(RHEA_OUT), exist_ok=True)
    
    # Save to file
    with open(CHEBI_OUT, "w", encoding="utf-8") as f:
        json.dump(chebi_mappings, f, ensure_ascii=False, indent=2)
    print("Saved ChEBI mappings to", CHEBI_OUT)
        
    with open(RHEA_OUT, "w", encoding="utf-8") as f:
        json.dump(rhea_mappings, f, ensure_ascii=False, indent=2)
    print("Saved Rhea mappings to", RHEA_OUT)
    
    print("Completed building mappings!")

if __name__ == "__main__":
    main()
