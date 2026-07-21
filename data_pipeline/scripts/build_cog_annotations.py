import os
import urllib.request
import urllib.parse
import json
import re
import sys
import ssl

# Resolve directories
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
COG_OUT = os.path.join(ROOT_DIR, "data", "reference", "cog_annotations.json")

# Standard COG category descriptions fallback
COG_CATEGORIES_FALLBACK = {
    "A": "RNA processing and modification",
    "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
    "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "V": "Defense mechanisms",
    "W": "Extracellular structures",
    "Y": "Nuclear structure",
    "Z": "Cytoskeleton",
    "X": "Mobilome: prophages, transposons"
}

def log(msg):
    print(msg)
    sys.stdout.flush()

def download_file(url, fallback_url=None):
    context = ssl._create_unverified_context()
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    try:
        with urllib.request.urlopen(req, context=context, timeout=30) as resp:
            return resp.read().decode('utf-8')
    except Exception as e:
        log(f"Failed to download {url}: {e}")
        if fallback_url:
            log(f"Trying fallback URL: {fallback_url}")
            fallback_req = urllib.request.Request(fallback_url, headers={'User-Agent': 'Mozilla/5.0'})
            try:
                with urllib.request.urlopen(fallback_req, context=context, timeout=30) as resp:
                    return resp.read().decode('utf-8')
            except Exception as e2:
                log(f"Fallback failed: {e2}")
        return None

def parse_cog_categories():
    log("Fetching COG categories...")
    url = "https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.fun.tab"
    fallback = "https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.fun.tab"
    content = download_file(url, fallback)
    
    cat_map = COG_CATEGORIES_FALLBACK.copy()
    if not content:
        log("Using default COG category descriptions fallback.")
        return cat_map
        
    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        if len(cols) >= 4:
            letter = cols[0].strip()
            desc = cols[3].strip()
            if len(letter) == 1 and letter.isupper():
                cat_map[letter] = desc
                
    log(f"Loaded {len(cat_map)} COG category descriptions.")
    return cat_map

def parse_cog_definitions():
    log("Fetching COG definitions...")
    url = "https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.def.tab"
    fallback = "https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab"
    content = download_file(url, fallback)
    
    def_map = {}
    if not content:
        log("Error: COG definitions could not be retrieved.")
        return def_map
        
    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
        cols = line.split('\t')
        if len(cols) >= 2:
            cog_id = cols[0].strip()
            categories = cols[1].strip()
            if cog_id.startswith("COG") and categories:
                def_map[cog_id] = categories[0]
                
    log(f"Loaded {len(def_map)} COG definitions.")
    return def_map

def fetch_uniprot_annotations():
    log("Streaming UniProt annotations for taxon 196627 via paginated search...")
    query = "taxonomy_id:196627"
    columns = "accession,gene_names,xref_eggnog"
    
    url = f"https://rest.uniprot.org/uniprotkb/search?query={urllib.parse.quote(query)}&format=tsv&fields={columns}&size=500"
    records = []
    
    context = ssl._create_unverified_context()
    while url:
        log(f"Fetching UniProt page: {url}")
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        try:
            with urllib.request.urlopen(req, context=context, timeout=20) as resp:
                content = resp.read().decode('utf-8')
                lines = content.splitlines()
                if not lines:
                    break
                    
                header = lines[0].split('\t')
                acc_idx = header.index("Entry")
                names_idx = header.index("Gene Names")
                egg_idx = header.index("eggNOG")
                
                for line in lines[1:]:
                    if not line.strip():
                        continue
                    cols = line.split('\t')
                    if len(cols) <= max(acc_idx, names_idx, egg_idx):
                        continue
                    acc = cols[acc_idx]
                    names = cols[names_idx]
                    eggnog_str = cols[egg_idx]
                    
                    cogs = []
                    if eggnog_str and eggnog_str != ";":
                        cogs = [c.strip() for c in eggnog_str.split(";") if c.strip().startswith("COG")]
                        
                    records.append({
                        "accession": acc,
                        "names": names,
                        "cogs": cogs
                    })
                    
                # Parse Link header to get next page
                link_header = resp.info().get("Link")
                url = None
                if link_header:
                    match = re.search(r'<(https://[^>]+)>;\s*rel="next"', link_header)
                    if match:
                        url = match.group(1)
        except Exception as e:
            log(f"Failed to fetch UniProt page: {e}")
            break
            
    log(f"Processed {len(records)} UniProt protein records in total.")
    return records

def main():
    cat_map = parse_cog_categories()
    def_map = parse_cog_definitions()
    uniprot_records = fetch_uniprot_annotations()
    
    if not def_map or not uniprot_records:
        log("Fatal error: Mapping databases are empty. Aborting.")
        return
        
    cog_annotations = {}
    cg_pattern = re.compile(r'\bcg\d{4}\b', re.IGNORECASE)
    cgl_pattern = re.compile(r'\bcgl\d{4}\b', re.IGNORECASE)
    
    for r in uniprot_records:
        if not r["cogs"]:
            continue
            
        primary_cog = r["cogs"][0]
        category_letter = def_map.get(primary_cog)
        if not category_letter:
            continue
            
        category_desc = cat_map.get(category_letter, "Function unknown")
        
        mapping = {
            "cog_id": primary_cog,
            "category": category_letter,
            "description": category_desc
        }
        
        names_str = r["names"]
        cg_matches = cg_pattern.findall(names_str)
        cgl_matches = cgl_pattern.findall(names_str)
        
        for match in cg_matches + cgl_matches:
            locus = match.lower()
            cog_annotations[locus] = mapping
            
    log(f"Successfully generated COG annotations for {len(cog_annotations)} genes.")
    
    os.makedirs(os.path.dirname(COG_OUT), exist_ok=True)
    with open(COG_OUT, "w", encoding="utf-8") as f:
        json.dump(cog_annotations, f, ensure_ascii=False, indent=2)
    log(f"Saved COG annotations to {COG_OUT}")

if __name__ == "__main__":
    main()
