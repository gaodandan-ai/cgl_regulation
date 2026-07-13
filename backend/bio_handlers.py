"""
backend/bio_handlers.py
=======================
Business-logic query handlers formerly embedded in run_server.CustomHTTPRequestHandler.
All functions are now plain callables that can be imported directly by app.py.
"""
import os
import csv
import json
import math
import urllib.parse
from collections import Counter

from gene_utils import (
    get_absolute_path,
    load_gene_mappings,
    normalize_gene_locus,
    expand_gene_aliases,
    GENE_NAMES, CG_TO_CGL, CGL_TO_CG,
)
from kegg_client import (
    load_organism_kegg_links,
    find_matching_kegg_pathways,
    KEGG_PATHWAY_NAMES, KEGG_CACHE_FILE, KEGG_CACHE_HIT,
    PATHWAY_TO_GENES, PATHWAY_REGULATION_CACHE,
)
from metabolic_mapper import load_metabolic_model_mappings

# ── Statistics helpers ────────────────────────────────────────────────────────

def hypergeom_sf(x: int, N: int, M: int, k: int) -> float:
    """P(X >= x) for the hypergeometric distribution."""
    total_comb = math.comb(N, k)
    if total_comb == 0:
        return 1.0
    total_prob = sum(math.comb(M, i) * math.comb(N - M, k - i) for i in range(x, min(k, M) + 1))
    return min(1.0, total_prob / total_comb)


def evidence_weight(evidence: str) -> float:
    text = (evidence or "").lower()
    if "experimental" in text and "predicted" in text:
        return 2.5
    if "experimental" in text:
        return 3.0
    if "predicted" in text:
        return 1.0
    return 0.5


def calculate_tf_pathway_impact(stat: dict, pathway_gene_count: int):
    edge_count  = max(1, stat["edge_count"])
    target_count = len(stat["target_genes"])
    coverage     = target_count / pathway_gene_count if pathway_gene_count else 0
    evidence_total = sum(evidence_weight(k) * v for k, v in stat["evidence"].items())
    evidence_avg   = evidence_total / edge_count
    binding_fraction    = stat["binding_site_edges"] / edge_count
    dominant_role_count = stat["roles"].most_common(1)[0][1] if stat["roles"] else 0
    direction_consistency = dominant_role_count / edge_count

    components = {
        "coverage":             round(coverage * 40, 2),
        "evidence":             round((evidence_avg / 3.0) * 25, 2),
        "binding_site":         round(binding_fraction * 20, 2),
        "direction_consistency": round(direction_consistency * 10, 2),
        "edge_support":         round(min(edge_count, 10) / 10 * 5, 2),
    }
    score      = round(sum(components.values()), 2)
    confidence = "high" if score >= 70 else "medium" if score >= 45 else "low"
    return score, components, confidence

# ── Regulon / TF helpers ──────────────────────────────────────────────────────

def get_regulatory_targets_for_tf(query: str):
    load_gene_mappings()
    q        = (query or "").strip().lower()
    resolved = normalize_gene_locus(q)
    targets, is_tf = {}, False

    path = get_absolute_path("data/reference/regulations.csv")
    if not os.path.exists(path):
        return is_tf, []
    try:
        with open(path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                tf_locus = (row.get("TF_locusTag") or "").strip()
                tf_name  = (row.get("TF_name") or "").strip()
                tf_aliases = expand_gene_aliases(tf_locus)
                if tf_name:
                    tf_aliases.add(tf_name.lower())
                if q not in tf_aliases and resolved not in tf_aliases:
                    continue
                is_tf  = True
                target = normalize_gene_locus(row.get("TG_locusTag", ""))
                if target and target not in targets:
                    targets[target] = {
                        "locus":      target,
                        "name":       row.get("TG_name", "").strip() or GENE_NAMES.get(target, target),
                        "regulation": row.get("Role", "").strip(),
                        "evidence":   row.get("Evidence", "").strip(),
                    }
    except Exception as e:
        print("Error reading regulations for metabolic impact:", e)
    return is_tf, list(targets.values())


def handle_regulon_enrichment(tf: str) -> dict:
    load_gene_mappings()
    load_organism_kegg_links()
    tf_lower    = tf.strip().lower()
    resolved_cg = CGL_TO_CG.get(tf_lower, tf)

    targets = []
    path = get_absolute_path("data/reference/regulations.csv")
    if os.path.exists(path):
        try:
            with open(path, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    tf_row  = row.get("TF_locusTag", "").strip().lower()
                    tf_name = row.get("TF_name",     "").strip().lower()
                    if tf_row == resolved_cg.lower() or (tf_name and tf_name == tf_lower):
                        tg = row.get("TG_locusTag", "").strip()
                        if tg and tg not in targets:
                            targets.append(tg)
        except Exception as e:
            print("Error reading regulations for enrichment:", e)

    if not targets:
        return {"error": f"No target genes found for transcription factor {tf}"}

    expanded = set()
    for tg in targets:
        tg_lower = tg.lower()
        expanded.add(tg_lower)
        if tg_lower in CG_TO_CGL:
            expanded.add(CG_TO_CGL[tg_lower].lower())
        if tg_lower in CGL_TO_CG:
            expanded.add(CGL_TO_CG[tg_lower].lower())

    from kegg_client import GENE_TO_PATHWAYS, PATHWAY_TO_GENES
    all_annotated = set(GENE_TO_PATHWAYS.keys())
    with_paths = expanded.intersection(all_annotated)
    canonical_reg = {CGL_TO_CG.get(g, g).lower() for g in with_paths}
    k = len(canonical_reg)

    if k == 0:
        return {"tf": tf, "regulon_size": len(targets), "annotated_regulon_size": 0, "pathways": []}

    canonical_pathway_to_genes = {
        pid: {CGL_TO_CG.get(g, g).lower() for g in genes}
        for pid, genes in PATHWAY_TO_GENES.items()
    }
    all_canonical = set()
    for genes in canonical_pathway_to_genes.values():
        all_canonical.update(genes)
    N = len(all_canonical)

    pathway_enrichments = []
    for pid, pathway_genes in canonical_pathway_to_genes.items():
        M    = len(pathway_genes)
        hits = canonical_reg.intersection(pathway_genes)
        x    = len(hits)
        if x > 0:
            pathway_enrichments.append({
                "pathway_id":     pid,
                "pathway_name":   KEGG_PATHWAY_NAMES.get(pid, pid),
                "p_value":        hypergeom_sf(x, N, M, k),
                "fold_enrichment": (x / k) / (M / N) if M > 0 else 0,
                "hits":           x,
                "total_genes":    M,
                "target_genes":   [{"locus": g, "name": GENE_NAMES.get(g, g.upper())} for g in hits],
            })
    pathway_enrichments.sort(key=lambda p: p["p_value"])
    return {
        "tf": tf,
        "regulon_size":          len(targets),
        "annotated_regulon_size": k,
        "total_annotated_genome": N,
        "pathways":              pathway_enrichments,
    }

# ── Metabolic impact / pathway handlers ──────────────────────────────────────

def handle_metabolic_impact(query: str) -> dict:
    mapping   = load_metabolic_model_mappings()
    q         = (query or "").strip()
    canonical = normalize_gene_locus(q)
    is_tf, targets = get_regulatory_targets_for_tf(q)

    seed_genes = targets if is_tf else [{"locus": canonical or q.lower(), "name": GENE_NAMES.get(canonical, q), "regulation": "", "evidence": ""}]
    mode       = "tf" if is_tf else "gene"

    affected_genes, pathway_stats, reaction_seen = [], {}, set()
    nodes, edges = {}, []

    query_id = canonical or q.lower()
    nodes[query_id] = {"id": query_id, "type": "TF" if is_tf else "gene", "label": GENE_NAMES.get(query_id, q)}

    for gene in seed_genes:
        locus = normalize_gene_locus(gene.get("locus", ""))
        if not locus:
            continue
        reactions = []
        for alias in expand_gene_aliases(locus):
            reactions.extend(mapping["gene_to_reactions"].get(alias, []))

        unique_rxns, local_seen = [], set()
        for r in reactions:
            key = f"{r['model']}:{r['id']}"
            if key not in local_seen:
                local_seen.add(key)
                unique_rxns.append(r)

        if is_tf:
            edges.append({
                "source":     query_id,
                "target":     locus,
                "type":       "regulates",
                "regulation": {"A": "activation", "R": "repression"}.get(gene.get("regulation"), "unknown"),
                "confidence": evidence_weight(gene.get("evidence", "")) / 3.0,
            })
        nodes[locus] = {"id": locus, "type": "gene", "label": GENE_NAMES.get(locus, gene.get("name") or locus)}

        affected_genes.append({**gene, "locus": locus, "name": GENE_NAMES.get(locus, gene.get("name") or locus),
                                "mapped_reaction_count": len(unique_rxns), "reactions": unique_rxns[:12]})

        for r in unique_rxns:
            rnode_id    = f"reaction:{r['model']}:{r['id']}"
            pway_label  = r.get("pathway_name") or "Unassigned pathway"
            pway_id     = r.get("pathway_id")   or pway_label
            pway_nid    = f"pathway:{r['model']}:{pway_id}"
            nodes[rnode_id] = {"id": rnode_id, "type": "reaction", "label": r.get("label") or r["id"],
                                "equation": r.get("equation",""), "model": r.get("model",""),
                                "ec_number": r.get("ec_number",""), "kcat": r.get("kcat"),
                                "molecular_weight": r.get("molecular_weight"), "kcat_MW": r.get("kcat_MW"),
                                "uniprot_ids": r.get("uniprot_ids",[]), "reaction_variant": r.get("reaction_variant",""),
                                "enzyme_constraint": r.get("enzyme_constraint",{})}
            nodes[pway_nid]  = {"id": pway_nid, "type": "pathway", "label": pway_label, "model": r.get("model","")}
            edges.append({"source": locus, "target": rnode_id, "type": "associated_with_reaction", "gpr_rule": r.get("gpr_rule","")})
            edges.append({"source": rnode_id, "target": pway_nid, "type": "belongs_to_pathway"})
            pkey = f"{r.get('model')}::{pway_label}"
            stat = pathway_stats.setdefault(pkey, {"id": pway_id, "name": pway_label, "model": r.get("model",""),
                                                    "gene_count": 0, "reaction_count": 0, "genes": set(), "reactions": set()})
            stat["genes"].add(locus)
            stat["reactions"].add(r["id"])
            reaction_seen.add(f"{r['model']}:{r['id']}")

    pathways = sorted([
        {"id": s["id"], "name": s["name"], "model": s["model"],
         "gene_count": len(s["genes"]), "reaction_count": len(s["reactions"]),
         "genes": sorted(s["genes"])[:20], "reactions": sorted(s["reactions"])[:20]}
        for s in pathway_stats.values()
    ], key=lambda p: (-p["gene_count"], -p["reaction_count"], p["name"]))

    mapped = [g for g in affected_genes if g["mapped_reaction_count"] > 0]
    return {
        "query": q, "mode": mode, "is_tf": is_tf,
        "model_mapping": {"loaded": mapping["loaded"], "models": mapping["models"],
                          "files": mapping["files"], "warnings": mapping["warnings"]},
        "summary": {"target_gene_count": len(seed_genes), "mapped_gene_count": len(mapped),
                    "reaction_count": len(reaction_seen), "pathway_count": len(pathways)},
        "affected_genes": affected_genes,
        "pathways":       pathways,
        "graph":          {"nodes": list(nodes.values()), "edges": edges},
    }


def handle_metabolic_pathways(query: str = "") -> dict:
    mapping    = load_metabolic_model_mappings()
    query_text = (query or "").strip().lower()
    pathway_index: dict = {}

    for alias, reactions in mapping.get("gene_to_reactions", {}).items():
        locus      = normalize_gene_locus(alias) or (alias or "").strip().lower()
        if not locus:
            continue
        gene_label = GENE_NAMES.get(locus, locus)
        for r in reactions or []:
            model      = r.get("model", "")
            pway_id    = r.get("pathway_id")   or r.get("pathway_name") or "Unassigned pathway"
            pway_name  = r.get("pathway_name") or pway_id or "Unassigned pathway"
            pway_key   = f"{model}::{pway_id or pway_name}"
            entry      = pathway_index.setdefault(pway_key, {"id": pway_id, "name": pway_name, "model": model, "genes": {}, "reactions": {}})
            rid        = r.get("id") or r.get("label") or "reaction"
            rxn_data   = {"reactionId": rid, "reactionName": r.get("label") or rid, "model": model,
                          "ecNumber": r.get("ec_number",""), "kcat": r.get("kcat"),
                          "molecularWeight": r.get("molecular_weight"), "kcatMW": r.get("kcat_MW"),
                          "uniprotIds": r.get("uniprot_ids",[]), "reactionVariant": r.get("reaction_variant",""),
                          "enzymeConstraint": r.get("enzyme_constraint",{})}
            entry["reactions"][rid] = rxn_data
            gene_entry = entry["genes"].setdefault(locus, {"geneId": locus, "geneLabel": gene_label, "reactions": {}})
            gene_entry["reactions"][rid] = {k: v for k, v in rxn_data.items() if k != "enzymeConstraint"}

    pathways = []
    for entry in pathway_index.values():
        genes = sorted([{"geneId": g["geneId"], "geneLabel": g["geneLabel"],
                         "reactions": sorted(g["reactions"].values(), key=lambda r: r["reactionId"])}
                        for g in entry["genes"].values()], key=lambda g: g["geneId"])
        reactions = sorted(entry["reactions"].values(), key=lambda r: r["reactionId"])
        item = {"pathwayId": entry["id"], "pathwayName": entry["name"], "model": entry["model"],
                "totalGenes": len(genes), "totalReactions": len(reactions), "genes": genes, "reactions": reactions}
        if query_text and query_text not in f"{item['pathwayId']} {item['pathwayName']}".lower():
            continue
        pathways.append(item)

    pathways.sort(key=lambda p: (-p["totalGenes"], -p["totalReactions"], p["pathwayName"]))
    return {
        "model_mapping": {"loaded": mapping.get("loaded", False), "models": mapping.get("models", []),
                          "files": mapping.get("files", []), "warnings": mapping.get("warnings", [])},
        "pathways": pathways,
    }


def handle_imodulon_reactions(imodulon_id: str) -> dict:
    path = get_absolute_path(os.path.join("data", "reference", "imodulon", "imodulon_gene_weights.json"))
    if not os.path.exists(path):
        return {"error": "iModulon weights file not found."}
    with open(path, "r", encoding="utf-8") as f:
        imod_data = json.load(f)
    if imodulon_id not in imod_data:
        return {"error": f"iModulon {imodulon_id} not found."}

    genes   = imod_data[imodulon_id].get("genes", {})
    mapping = load_metabolic_model_mappings()
    reactions, seen = [], set()
    for gene_locus in genes:
        for alias in expand_gene_aliases(gene_locus):
            for r in mapping.get("gene_to_reactions", {}).get(alias, []):
                key = (r["model"], r["id"])
                if key not in seen:
                    seen.add(key)
                    reactions.append({"reactionId": r["id"], "name": r["label"], "model": r["model"],
                                      "equation": r.get("equation",""), "gpr_rule": r.get("gpr_rule",""),
                                      "pathway_id": r.get("pathway_id",""), "pathway_name": r.get("pathway_name",""),
                                      "ec_number": r.get("ec_number",""), "kcat": r.get("kcat"),
                                      "molecular_weight": r.get("molecular_weight"), "kcat_MW": r.get("kcat_MW")})
    return {"imodulon_id": imodulon_id, "reactions": reactions}


def handle_imodulon_simulation(imodulon_id: str) -> dict:
    path = get_absolute_path(os.path.join("data", "reference", "imodulon", "imodulon_gene_weights.json"))
    if not os.path.exists(path):
        return {"error": "iModulon weights file not found."}
    with open(path, "r", encoding="utf-8") as f:
        imod_data = json.load(f)
    if imodulon_id not in imod_data:
        return {"error": f"iModulon {imodulon_id} not found."}

    genes = list(imod_data[imodulon_id].get("genes", {}).keys())
    if not genes:
        return {"error": "No genes in this iModulon."}

    fba_res, ecfba_res, warnings = {}, {}, []
    json_model_path = get_absolute_path(os.path.join(
        "data", "reference", "model", "ecCGL1-main", "ecCGL1-main", "model",
        "iCW773_irr_enz_constraint.json"))
    enzyme_perturbations = {g: 0.0 for g in genes}

    try:
        from backend.model_loader import load_model_if_needed
        from backend.simulation   import run_gene_set_knockout
        model   = load_model_if_needed()
        fba_res = run_gene_set_knockout(model, genes, track_reaction_ids=["EX_lys_L_e", "EX_glu_L_e"])
    except Exception as e:
        warnings.append(f"Standard FBA failed: {e}")

    try:
        from backend.simulation import run_ecfba_simulation
        ecfba_res = {
            obj: (run_ecfba_simulation(json_model_path, 0.2, enzyme_perturbations, obj, 30.0).get("flux") or 0.0)
            for obj in ("growth", "lysine", "glutamate")
        }
        ecfba_res["status"] = "success"
    except Exception as e:
        warnings.append(f"ecFBA failed: {e}")

    return {"imodulon_id": imodulon_id, "fba": fba_res, "ecfba": ecfba_res, "warnings": warnings}


def handle_tf_simulation(tf_id: str) -> dict:
    load_gene_mappings()
    tf_lower    = tf_id.strip().lower()
    resolved_cg = CGL_TO_CG.get(tf_lower, tf_id)

    targets = []
    path = get_absolute_path("data/reference/regulations.csv")
    if os.path.exists(path):
        try:
            with open(path, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    tf_row  = row.get("TF_locusTag", "").strip().lower()
                    tf_name = row.get("TF_name",     "").strip().lower()
                    if tf_row == resolved_cg.lower() or (tf_name and tf_name == tf_lower):
                        tg = row.get("TG_locusTag", "").strip()
                        if tg and tg not in targets:
                            targets.append(tg)
        except Exception as e:
            print("Error reading regulations for TF simulation:", e)

    if not targets:
        return {"error": f"No target genes found for transcription factor {tf_id}"}

    fba_res, ecfba_res, warnings = {}, {}, []
    json_model_path = get_absolute_path(os.path.join(
        "data", "reference", "model", "ecCGL1-main", "ecCGL1-main", "model",
        "iCW773_irr_enz_constraint.json"))
    enzyme_perturbations = {g: 0.0 for g in targets}

    try:
        from backend.model_loader import load_model_if_needed
        from backend.simulation   import run_gene_set_knockout
        model   = load_model_if_needed()
        fba_res = run_gene_set_knockout(model, targets, track_reaction_ids=["EX_lys_L_e", "EX_glu_L_e"])
    except Exception as e:
        warnings.append(f"Standard FBA failed: {e}")

    try:
        from backend.simulation import run_ecfba_simulation
        ecfba_res = {
            obj: (run_ecfba_simulation(json_model_path, 0.2, enzyme_perturbations, obj, 30.0).get("flux") or 0.0)
            for obj in ("growth", "lysine", "glutamate")
        }
        ecfba_res["status"] = "success"
    except Exception as e:
        warnings.append(f"ecFBA failed: {e}")

    return {"tf_id": tf_id, "targets_count": len(targets), "fba": fba_res, "ecfba": ecfba_res, "warnings": warnings}


def handle_pathway_regulation(query: str) -> dict:
    cache_key = (query or "").strip().lower()
    if cache_key in PATHWAY_REGULATION_CACHE:
        return PATHWAY_REGULATION_CACHE[cache_key]

    load_gene_mappings()
    load_organism_kegg_links()

    matches  = find_matching_kegg_pathways(query)
    selected = matches[:4]
    pathway_genes, pathway_ids = set(), set()
    for pw in selected:
        pid = pw["id"]
        pathway_ids.update({pid, f"path:{pid}"})
        pathway_genes.update(PATHWAY_TO_GENES.get(pid, set()))
        pathway_genes.update(PATHWAY_TO_GENES.get(f"path:{pid}", set()))

    canonical_pw_genes = {normalize_gene_locus(g) for g in pathway_genes if normalize_gene_locus(g)}

    tf_stats, regulated, edge_examples = {}, set(), []
    path = get_absolute_path("data/reference/regulations.csv")
    if os.path.exists(path) and canonical_pw_genes:
        try:
            with open(path, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    target = normalize_gene_locus(row.get("TG_locusTag"))
                    if target not in canonical_pw_genes:
                        continue
                    tf        = (row.get("TF_locusTag") or "").strip()
                    if not tf:
                        continue
                    tf_key    = tf.lower()
                    tf_name   = (row.get("TF_name") or tf).strip() or tf
                    role      = (row.get("Role")     or "").strip() or "unknown"
                    evidence  = (row.get("Evidence") or "").strip() or "unknown"
                    source    = (row.get("Source")   or "").strip() or "local"
                    tg_name   = (row.get("TG_name")  or row.get("TG_locusTag") or target).strip()
                    binding   = (row.get("Binding_site") or "").strip()

                    if tf_key not in tf_stats:
                        tf_stats[tf_key] = {"tf_locus": tf, "tf_name": tf_name, "edge_count": 0,
                                            "target_genes": set(), "roles": Counter(), "evidence": Counter(),
                                            "sources": Counter(), "binding_site_edges": 0, "examples": []}
                    stat = tf_stats[tf_key]
                    stat["edge_count"] += 1
                    stat["target_genes"].add(target)
                    stat["roles"][role]      += 1
                    stat["evidence"][evidence] += 1
                    stat["sources"][source]    += 1
                    if binding:
                        stat["binding_site_edges"] += 1
                    if len(stat["examples"]) < 5:
                        stat["examples"].append({"target_locus": target, "target_name": tg_name,
                                                  "role": role, "evidence": evidence, "has_binding_site": bool(binding)})
                    regulated.add(target)
                    if len(edge_examples) < 12:
                        edge_examples.append({"tf": tf, "tf_name": tf_name, "target_locus": target,
                                               "target_name": tg_name, "role": role, "evidence": evidence,
                                               "source": source, "has_binding_site": bool(binding)})
        except Exception as e:
            print("Error projecting pathway genes onto regulatory network:", e)

    regulators = sorted([
        {
            "tf_locus": s["tf_locus"], "tf_name": s["tf_name"],
            "impact_score":    (sc := calculate_tf_pathway_impact(s, len(canonical_pw_genes)))[0],
            "score_components": sc[1], "confidence": sc[2],
            "edge_count":    s["edge_count"],
            "target_count":  len(s["target_genes"]),
            "coverage":      round(len(s["target_genes"]) / len(canonical_pw_genes), 4) if canonical_pw_genes else 0,
            "target_genes":  sorted(s["target_genes"]),
            "roles":         dict(s["roles"].most_common()),
            "evidence":      dict(s["evidence"].most_common()),
            "sources":       dict(s["sources"].most_common()),
            "binding_site_edges": s["binding_site_edges"],
            "examples":      s["examples"],
        }
        for s in tf_stats.values()
    ], key=lambda r: (r["impact_score"], r["target_count"], r["edge_count"]), reverse=True)

    result = {
        "query":           query,
        "matched_pathways": selected,
        "all_matches_count": len(matches),
        "pathway_ids":     sorted(pathway_ids),
        "pathway_gene_count": len(canonical_pw_genes),
        "regulated_gene_count": len(regulated),
        "regulator_count": len(regulators),
        "coverage":        round(len(regulated) / len(canonical_pw_genes), 4) if canonical_pw_genes else 0,
        "pathway_genes":   [{"locus": g, "name": GENE_NAMES.get(g, g.upper()),
                              "cgl_locus": CG_TO_CGL.get(g, "")} for g in sorted(canonical_pw_genes)],
        "regulators":      regulators[:25],
        "edge_examples":   edge_examples,
        "external_resources": {
            "kegg":         [p["link"] for p in selected],
            "biocyc_search": f"https://biocyc.org/gene-search.shtml?orgid=CORYNE&query={urllib.parse.quote(query or '')}",
            "note":         "BioCyc and genome-scale model overlays can be added when local reaction/SBML files are supplied.",
        },
        "cache": {"enabled": True, "path": KEGG_CACHE_FILE, "loaded_from_disk": KEGG_CACHE_HIT},
    }
    PATHWAY_REGULATION_CACHE[cache_key] = result
    return result
