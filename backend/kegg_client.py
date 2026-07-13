"""
backend/kegg_client.py
======================
KEGG REST API client: local cache management, pathway name loading,
and organism-level gene→pathway link loading.
"""
import os
import json
import time
import threading
import urllib.request
import urllib.parse

from gene_utils import (
    get_absolute_path,
    ROOT_DIR,
    load_gene_mappings,
    CG_TO_CGL, CGL_TO_CG, GENE_NAMES, NAME_TO_CG,
)

# ── Global caches ─────────────────────────────────────────────────────────────
KEGG_PATHWAY_NAMES:      dict  = {}   # pathway ID -> clean name
PATHWAY_NAMES_MUTEX             = threading.Lock()
GENE_TO_PATHWAYS:        dict  = {}   # gene locus -> set of pathway IDs
PATHWAY_TO_GENES:        dict  = {}   # pathway ID -> set of gene loci
ORGANISM_PATHWAYS_LOADED: bool = False
GENE_PATHWAYS_CACHE:     dict  = {}   # (cg, cgl) -> {pathways, go_terms}
PATHWAY_REGULATION_CACHE: dict = {}
KEGG_CACHE_LOADED:       bool  = False
KEGG_CACHE_HIT:          bool  = False

KEGG_CACHE_DIR  = get_absolute_path(os.path.join("data", "reference", "kegg_cache"))
KEGG_CACHE_FILE = get_absolute_path(os.path.join("data", "reference", "kegg_cache", "kegg_cgl_cgb.json"))

# ── Cache persistence ─────────────────────────────────────────────────────────

def load_kegg_cache() -> None:
    global KEGG_CACHE_LOADED, KEGG_CACHE_HIT, ORGANISM_PATHWAYS_LOADED
    if KEGG_CACHE_LOADED:
        return
    KEGG_CACHE_LOADED = True
    if not os.path.exists(KEGG_CACHE_FILE):
        return
    try:
        with open(KEGG_CACHE_FILE, "r", encoding="utf-8") as f:
            data = json.load(f)
        KEGG_PATHWAY_NAMES.update(data.get("pathway_names", {}))
        for gene, pathways in data.get("gene_to_pathways", {}).items():
            GENE_TO_PATHWAYS[gene] = set(pathways)
        for pathway, genes in data.get("pathway_to_genes", {}).items():
            PATHWAY_TO_GENES[pathway] = set(genes)
        KEGG_CACHE_HIT = True
        if GENE_TO_PATHWAYS and PATHWAY_TO_GENES:
            ORGANISM_PATHWAYS_LOADED = True
    except Exception as e:
        print("Error loading KEGG cache:", e)


def save_kegg_cache() -> None:
    try:
        os.makedirs(KEGG_CACHE_DIR, exist_ok=True)
        data = {
            "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "organisms": ["cgb", "cgl"],
            "pathway_names": KEGG_PATHWAY_NAMES,
            "gene_to_pathways": {k: sorted(v) for k, v in GENE_TO_PATHWAYS.items()},
            "pathway_to_genes": {k: sorted(v) for k, v in PATHWAY_TO_GENES.items()},
        }
        with open(KEGG_CACHE_FILE, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
    except Exception as e:
        print("Error saving KEGG cache:", e)

# ── Pathway name loading ──────────────────────────────────────────────────────

def load_kegg_pathway_names() -> None:
    global KEGG_PATHWAY_NAMES
    load_kegg_cache()
    if KEGG_PATHWAY_NAMES:
        return

    with PATHWAY_NAMES_MUTEX:
        if KEGG_PATHWAY_NAMES:
            return

        for org in ("cgb", "cgl"):
            try:
                url = f"https://rest.kegg.jp/list/pathway/{org}"
                req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
                with urllib.request.urlopen(req, timeout=10) as resp:
                    for line in resp.read().decode("utf-8").splitlines():
                        if "\t" in line:
                            pid, pname = line.split("\t", 1)
                            pname_clean = pname.split(" - Corynebacterium")[0].strip()
                            pid_clean   = pid.strip()
                            KEGG_PATHWAY_NAMES[pid_clean] = pname_clean
                            if not pid_clean.startswith("path:"):
                                KEGG_PATHWAY_NAMES[f"path:{pid_clean}"] = pname_clean
            except Exception as e:
                print(f"Error loading {org} pathway names:", e)

        if KEGG_PATHWAY_NAMES:
            save_kegg_cache()


def load_organism_kegg_links() -> None:
    global ORGANISM_PATHWAYS_LOADED
    load_kegg_cache()
    if ORGANISM_PATHWAYS_LOADED:
        return

    load_kegg_pathway_names()

    for org in ("cgb", "cgl"):
        try:
            url = f"https://rest.kegg.jp/link/pathway/{org}"
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=10) as resp:
                for line in resp.read().decode("utf-8").splitlines():
                    if "\t" in line:
                        gene_raw, path_raw = line.split("\t", 1)
                        gene = gene_raw.replace(f"{org}:", "").strip().lower()
                        path = path_raw.replace("path:", "").strip()
                        GENE_TO_PATHWAYS.setdefault(gene, set()).add(path)
                        PATHWAY_TO_GENES.setdefault(path, set()).add(gene)
        except Exception as e:
            print(f"Error loading {org} pathway links:", e)

    ORGANISM_PATHWAYS_LOADED = True
    if GENE_TO_PATHWAYS and PATHWAY_TO_GENES:
        save_kegg_cache()


def find_matching_kegg_pathways(query: str) -> list:
    load_organism_kegg_links()
    q = (query or "").strip().lower()
    if not q:
        return []

    q_digits = "".join(ch for ch in q if ch.isdigit())
    matches, seen = [], set()

    for pid, name in KEGG_PATHWAY_NAMES.items():
        clean_pid  = pid.replace("path:", "")
        pid_lower  = clean_pid.lower()
        name_lower = name.lower()
        pid_digits = "".join(ch for ch in clean_pid if ch.isdigit())
        is_match = (
            q == pid_lower
            or q in name_lower
            or (q_digits and q_digits == pid_digits)
            or (q_digits and pid_lower.endswith(q_digits))
        )
        if not is_match or clean_pid in seen:
            continue
        seen.add(clean_pid)
        matches.append({
            "id":   clean_pid,
            "name": name,
            "link": f"https://www.kegg.jp/kegg-bin/show_pathway?{clean_pid}",
        })

    if not matches and q_digits:
        for prefix in ("cgl", "cgb"):
            pid = f"{prefix}{q_digits}"
            if pid in PATHWAY_TO_GENES and pid not in seen:
                seen.add(pid)
                matches.append({
                    "id":   pid,
                    "name": KEGG_PATHWAY_NAMES.get(pid, pid),
                    "link": f"https://www.kegg.jp/kegg-bin/show_pathway?{pid}",
                })

    matches.sort(key=lambda p: (0 if p["id"].lower().endswith(q_digits) and q_digits else 1, p["name"]))
    return matches


def get_gene_pathways_and_go(cg: str, cgl: str) -> dict:
    """Fetch KEGG pathways and UniProt GO terms for a gene (cached)."""
    import urllib.request, json as _json
    cg  = (cg  or "").strip()
    cgl = (cgl or "").strip()
    cache_key = (cg.lower(), cgl.lower())
    if cache_key in GENE_PATHWAYS_CACHE:
        return GENE_PATHWAYS_CACHE[cache_key]

    load_kegg_pathway_names()
    pathways, seen_pids = [], set()

    # cgb (Bielefeld genome) paths via cg_locus
    if cg:
        try:
            url = f"https://rest.kegg.jp/link/pathway/cgb:{cg}"
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=10) as resp:
                for line in resp.read().decode("utf-8").splitlines():
                    if "\t" in line:
                        _, pid_raw = line.split("\t", 1)
                        pid_clean = pid_raw.replace("path:", "").strip()
                        pid_num   = "".join(c for c in pid_clean if c.isdigit())
                        if pid_num not in seen_pids:
                            seen_pids.add(pid_num)
                            pathways.append({
                                "id":     pid_clean,
                                "name":   KEGG_PATHWAY_NAMES.get(pid_clean, pid_clean),
                                "link":   f"https://www.kegg.jp/kegg-bin/show_pathway?{pid_clean}+cgb:{cg}",
                                "source": "KEGG",
                            })
        except Exception as e:
            print(f"Error querying cgb pathways for {cg}:", e)

    # cgl (Kyowa Hakko genome) paths via cgl_locus
    if cgl:
        cgl_norm = ("Cgl" + cgl[3:]) if len(cgl) > 3 and cgl.lower().startswith("cgl") else cgl
        try:
            url = f"https://rest.kegg.jp/link/pathway/cgl:{cgl_norm}"
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=10) as resp:
                for line in resp.read().decode("utf-8").splitlines():
                    if "\t" in line:
                        _, pid_raw = line.split("\t", 1)
                        pid_clean = pid_raw.replace("path:", "").strip()
                        pid_num   = "".join(c for c in pid_clean if c.isdigit())
                        if pid_num not in seen_pids:
                            seen_pids.add(pid_num)
                            pathways.append({
                                "id":     pid_clean,
                                "name":   KEGG_PATHWAY_NAMES.get(pid_clean, pid_clean),
                                "link":   f"https://www.kegg.jp/kegg-bin/show_pathway?{pid_clean}+cgl:{cgl_norm}",
                                "source": "KEGG",
                            })
        except Exception as e:
            print(f"Error querying cgl pathways for {cgl_norm}:", e)

    # GO terms from UniProt
    go_terms, seen_gos = [], set()
    query_tag = cg or cgl
    if query_tag:
        try:
            url = (
                f"https://rest.uniprot.org/uniprotkb/search"
                f"?query=gene:{query_tag}+AND+organism_id:196627"
                f"&fields=id,accession,go&format=json"
            )
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=10) as resp:
                data = _json.loads(resp.read().decode("utf-8"))
            for ref in (data.get("results") or [{}])[0].get("uniProtKBCrossReferences", []):
                if ref.get("database") != "GO":
                    continue
                go_id = ref.get("id")
                go_term_val = next(
                    (p.get("value") for p in ref.get("properties", []) if p.get("key") == "GoTerm"),
                    "",
                )
                if go_id and go_term_val and go_id not in seen_gos:
                    seen_gos.add(go_id)
                    go_type, go_name = "GO", go_term_val
                    if ":" in go_term_val:
                        t_code, t_name = go_term_val.split(":", 1)
                        go_type = {"P": "GO Process", "F": "GO Function", "C": "GO Component"}.get(t_code, "GO")
                        go_name = t_name.strip()
                    go_terms.append({
                        "id":   go_id,
                        "name": go_name,
                        "type": go_type,
                        "link": f"https://www.ebi.ac.uk/QuickGO/term/{go_id}",
                    })
        except Exception as e:
            print(f"Error querying UniProt GO terms for {query_tag}:", e)

    result = {"pathways": pathways, "go_terms": go_terms}
    GENE_PATHWAYS_CACHE[cache_key] = result
    return result
