from fastapi import APIRouter, HTTPException, Response
import os
import json
import logging
from collections import deque

try:
    from db_manager import get_db_manager
    _DB_MANAGER_AVAILABLE = True
except ImportError:
    try:
        from backend.db_manager import get_db_manager
        _DB_MANAGER_AVAILABLE = True
    except ImportError:
        get_db_manager = None
        _DB_MANAGER_AVAILABLE = False

try:
    from graph_engine import get_graph_engine
except ImportError:
    try:
        from backend.graph_engine import get_graph_engine
    except ImportError:
        get_graph_engine = None

try:
    from schemas import (
        GraphCascadeResponse,
        GraphMotifResponse,
        InterventionTargetsResponse,
        InterventionTargetSchema,
        HTTPError
    )
except ImportError:
    from backend.schemas import (
        GraphCascadeResponse,
        GraphMotifResponse,
        InterventionTargetsResponse,
        InterventionTargetSchema,
        HTTPError
    )
from services.reference_data import STRING_INTERACTIONS

try:
    import run_server
except ImportError:
    import backend.run_server as run_server

router = APIRouter(tags=["Network & Regulatory"])
logger = logging.getLogger("app.routers.network")

_CENTRALITY_DATA = None


def _load_centrality():
    global _CENTRALITY_DATA
    if _CENTRALITY_DATA is not None:
        return _CENTRALITY_DATA
    backend_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    root_dir = os.path.dirname(backend_dir)
    path = os.path.join(root_dir, "data", "reference", "network_centrality.json")
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        _CENTRALITY_DATA = json.load(f)
    return _CENTRALITY_DATA


@router.get("/api/analysis/rna-seq")
def get_rna_seq_analysis():
    backend_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    root_dir = os.path.dirname(backend_dir)
    file_path = os.path.join(root_dir, "data", "reference", "rna_seq_analysis_results.json")
    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="Heat stress RNA-Seq analysis data is not publicly available yet. It will be released upon publication.")
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        return data
    except Exception as e:
        logger.error(f"Failed to read RNA-Seq analysis file: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Failed to read analysis results: {str(e)}")


@router.get("/api/analysis/dynamic-grn")
def get_dynamic_grn():
    data = get_rna_seq_analysis()
    return data.get("dynamic_grn", {})


@router.get("/api/analysis/causal-grn")
def get_causal_grn():
    data = get_rna_seq_analysis()
    return data.get("causal_grn", [])


@router.get("/api/analysis/metabolic-coupling")
def get_metabolic_coupling():
    data = get_rna_seq_analysis()
    return data.get("metabolic_coupling", {})


@router.get("/api/analysis/tf-motif-enrichment")
def get_tf_motif_enrichment():
    data = get_rna_seq_analysis()
    return data.get("motif_enrichment", {})


@router.get("/api/kegg_pathways")
def kegg_pathways(gene_name: str = "", accession: str = ""):
    try:
        result = run_server.handle_kegg_pathways(gene_name, accession)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/pathway_regulation")
def pathway_regulation(pathway_id: str = ""):
    try:
        result = run_server.handle_pathway_regulation(pathway_id)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/regulon_enrichment")
def regulon_enrichment(tf: str = ""):
    try:
        result = run_server.handle_regulon_enrichment(tf)
        return result
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/go_enrichment")
def go_enrichment(tf: str = ""):
    try:
        result = run_server.handle_go_enrichment(tf)
        return result
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/network/centrality")
def get_network_centrality(response: Response, limit: int = 30, tfs_only: bool = True):
    response.headers["Cache-Control"] = "public, max-age=1800"
    data = _load_centrality()
    if data is None:
        raise HTTPException(status_code=503, detail="Centrality data not available. Run scripts/network_centrality.py first.")
    nodes = data.get("nodes", {})
    result = [v for v in nodes.values() if (not tfs_only or v.get("is_tf"))]
    result.sort(key=lambda x: x.get("importance", 0), reverse=True)
    return {
        "_meta": data.get("_meta", {}),
        "top_tfs": result[:limit],
        "total_tfs": sum(1 for v in nodes.values() if v.get("is_tf")),
        "total_nodes": len(nodes),
    }


@router.get("/api/network/centrality/{locus}")
def get_centrality_for_gene(response: Response, locus: str):
    response.headers["Cache-Control"] = "public, max-age=1800"
    data = _load_centrality()
    if data is None:
        raise HTTPException(status_code=503, detail="Centrality data not available.")
    nodes = data.get("nodes", {})
    locus_lower = locus.strip().lower()
    entry = nodes.get(locus_lower) or nodes.get(locus)
    if entry is None:
        for k, v in nodes.items():
            if k.lower() == locus_lower:
                entry = v
                break
    if entry is None:
        raise HTTPException(status_code=404, detail="Gene not found in centrality data.")
    return entry


@router.get("/api/analysis/string_ppi")
def get_string_ppi(gene: str = "", min_score: int = 400, limit: int = 50):
    if not gene:
        raise HTTPException(status_code=400, detail="Missing gene parameter")
    gene_lower = gene.strip().lower()
    if gene_lower == "_meta":
        raise HTTPException(status_code=400, detail="Invalid gene identifier")

    partners = STRING_INTERACTIONS.get(gene_lower)
    if not partners:
        try:
            aliases = run_server.expand_gene_aliases(gene_lower)
            for alias in aliases:
                a_lower = alias.lower()
                if a_lower in STRING_INTERACTIONS:
                    partners = STRING_INTERACTIONS[a_lower]
                    gene_lower = a_lower
                    break
        except Exception:
            pass

    if not partners:
        return {
            "gene": gene,
            "resolved_id": gene_lower,
            "partners": [],
            "total": 0,
            "filtered": 0,
        }

    filtered = [p for p in partners if p.get("score", 0) >= min_score]
    total = len(filtered)
    filtered = filtered[:limit]

    meta = STRING_INTERACTIONS.get("_meta", {})
    return {
        "gene": gene,
        "resolved_id": gene_lower,
        "partners": filtered,
        "total": total,
        "filtered": total,
        "string_meta": {
            "version": meta.get("version", "12.0"),
            "min_score": min_score,
            "n_genes": meta.get("n_genes", 0),
            "n_edges": meta.get("n_edges", 0),
            "n_high_conf": meta.get("n_high_conf", 0),
        },
    }


@router.get("/api/analysis/string_ppi/neighborhood")
def get_ppi_neighborhood(
    genes: str = "",
    min_score: int = 400,
    limit_per_gene: int = 30,
):
    if not genes:
        raise HTTPException(status_code=400, detail="Missing genes parameter")

    gene_list = [g.strip().lower() for g in genes.split(",") if g.strip()]
    if not gene_list:
        raise HTTPException(status_code=400, detail="No valid gene identifiers provided")

    nodes: dict[str, dict] = {}
    edges: dict[str, dict] = {}

    def resolve_locus(raw: str) -> str:
        if raw in STRING_INTERACTIONS:
            return raw
        try:
            for alias in run_server.expand_gene_aliases(raw):
                al = alias.lower()
                if al in STRING_INTERACTIONS:
                    return al
        except Exception:
            pass
        return raw

    def gene_display_name(locus: str) -> str:
        try:
            name = run_server.normalize_gene_locus(locus)
            return name if name and name != locus else locus
        except Exception:
            return locus

    def add_node(locus: str, is_seed: bool = False):
        if locus not in nodes:
            nodes[locus] = {
                "id": locus,
                "name": gene_display_name(locus),
                "degree": 0,
                "is_seed": is_seed,
            }
        if is_seed:
            nodes[locus]["is_seed"] = True

    def add_edge(src: str, tgt: str, pdata: dict):
        key = f"{min(src, tgt)}::{max(src, tgt)}"
        if key not in edges:
            edges[key] = {
                "id": f"ppi-{src}-{tgt}",
                "source": src,
                "target": tgt,
                "score": pdata.get("score", 0),
                "experimental": pdata.get("experimental", 0),
                "database": pdata.get("database", 0),
                "coexpression": pdata.get("coexpression", 0),
                "neighborhood": pdata.get("neighborhood", 0),
                "cooccurrence": pdata.get("cooccurrence", 0),
                "textmining": pdata.get("textmining", 0),
                "fusion": pdata.get("fusion", 0),
                "type": pdata.get("type", ""),
            }
        else:
            if pdata.get("score", 0) > edges[key]["score"]:
                edges[key]["score"] = pdata["score"]
        nodes[src]["degree"] = nodes[src].get("degree", 0) + 1
        nodes[tgt]["degree"] = nodes[tgt].get("degree", 0) + 1

    for raw_gene in gene_list:
        locus = resolve_locus(raw_gene)
        add_node(locus, is_seed=True)

        partners = STRING_INTERACTIONS.get(locus, [])
        filtered = [p for p in partners if p.get("score", 0) >= min_score]
        filtered = filtered[:limit_per_gene]

        for p in filtered:
            partner_id = p["partner"].lower()
            add_node(partner_id, is_seed=False)
            add_edge(locus, partner_id, p)

    meta = STRING_INTERACTIONS.get("_meta", {})
    return {
        "genes": gene_list,
        "nodes": list(nodes.values()),
        "edges": list(edges.values()),
        "total_nodes": len(nodes),
        "total_edges": len(edges),
        "string_meta": {
            "version": meta.get("version", "12.0"),
            "min_score": min_score,
            "n_genes": meta.get("n_genes", 0),
            "n_edges": meta.get("n_edges", 0),
        },
    }


@router.get("/api/analysis/cross_network")
def get_cross_network(gene: str = "", min_ppi_score: int = 400):
    if not gene:
        raise HTTPException(status_code=400, detail="'gene' parameter is required")

    is_tf, reg_targets = run_server.get_regulatory_targets_for_tf(gene)

    reg_by_locus = {}
    for t in reg_targets:
        locus = (t.get("locus") or "").lower()
        if locus:
            reg_by_locus[locus] = t

    q_lower = gene.strip().lower()
    ppi_key = q_lower
    if ppi_key not in STRING_INTERACTIONS:
        for alias in run_server.expand_gene_aliases(q_lower):
            if alias in STRING_INTERACTIONS:
                ppi_key = alias
                break

    ppi_raw = STRING_INTERACTIONS.get(ppi_key, [])
    ppi_partners = [p for p in ppi_raw if p.get("score", 0) >= min_ppi_score]

    ppi_by_locus = {}
    for p in ppi_partners:
        partner_id = (p.get("partner") or "").lower()
        try:
            canonical = run_server.normalize_gene_locus(partner_id)
        except Exception:
            canonical = partner_id
        key = canonical or partner_id
        ppi_by_locus[key] = dict(p, canonical=key)

    cross_links = []
    for locus, reg_info in reg_by_locus.items():
        ppi_hit = ppi_by_locus.get(locus)
        if ppi_hit is None:
            for alias in run_server.expand_gene_aliases(locus):
                if alias in ppi_by_locus:
                    ppi_hit = ppi_by_locus[alias]
                    break
        if ppi_hit:
            cross_links.append({
                "gene": locus,
                "name": reg_info.get("name", locus),
                "regulation_role": reg_info.get("regulation", reg_info.get("role", "?")),
                "evidence": reg_info.get("evidence", ""),
                "ppi_score": ppi_hit.get("score", 0),
                "ppi_experimental": ppi_hit.get("experimental", 0),
                "ppi_database": ppi_hit.get("database", 0),
                "ppi_coexpression": ppi_hit.get("coexpression", 0),
            })

    cross_links.sort(key=lambda x: x["ppi_score"], reverse=True)

    return {
        "query": gene,
        "is_tf": is_tf,
        "n_regulatory": len(reg_targets),
        "n_ppi_partners": len(ppi_partners),
        "n_cross_links": len(cross_links),
        "cross_links": cross_links,
        "regulatory_targets": reg_targets,
        "ppi_summary": [{"gene": p.get("partner", ""), "score": p.get("score", 0)} for p in ppi_partners[:20]],
    }


@router.get("/api/analysis/tf_similarity")
def get_tf_similarity(min_targets: int = 3, metric: str = "jaccard", top_n: int = 40):
    import csv

    path = run_server.get_absolute_path("data/reference/regulations.csv")
    if not os.path.exists(path):
        raise HTTPException(status_code=404, detail="regulations.csv not found")

    tf_targets: dict = {}
    tf_names: dict = {}
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            tf_locus = (row.get("TF_locusTag") or "").strip()
            tf_name = (row.get("TF_name") or "").strip()
            tg_locus = (row.get("TG_locusTag") or "").strip()
            if not tf_locus or not tg_locus:
                continue
            tf_targets.setdefault(tf_locus, set()).add(tg_locus)
            if tf_name and tf_name != tf_locus:
                tf_names[tf_locus] = tf_name

    qualified = [(tf, tgts) for tf, tgts in tf_targets.items() if len(tgts) >= min_targets]
    qualified.sort(key=lambda x: len(x[1]), reverse=True)
    qualified = qualified[:top_n]

    tfs = [q[0] for q in qualified]
    targets = [q[1] for q in qualified]

    n = len(tfs)
    matrix = []
    shared_detail = {}

    for i in range(n):
        row_vals = []
        for j in range(n):
            if i == j:
                row_vals.append(1.0)
            else:
                inter = targets[i] & targets[j]
                union = targets[i] | targets[j]
                jac = round(len(inter) / len(union), 4) if union else 0.0
                row_vals.append(jac)
                if jac > 0 and i < j:
                    key = f"{tfs[i]}__{tfs[j]}"
                    shared_detail[key] = sorted(inter)
        matrix.append(row_vals)

    tf_labels = [{"locus": tf, "name": tf_names.get(tf, tf), "n_targets": len(targets[i])} for i, tf in enumerate(tfs)]

    return {
        "n_tfs": n,
        "tfs": tf_labels,
        "matrix": matrix,
        "shared": shared_detail,
        "metric": metric,
    }


@router.get("/api/analysis/string_ppi/shortest_path")
def get_ppi_shortest_path(
    source: str = "",
    target: str = "",
    min_score: int = 400,
    max_hops: int = 6,
):
    if not source or not target:
        raise HTTPException(status_code=400, detail="Both 'source' and 'target' parameters are required")

    def resolve(raw: str) -> str:
        r = raw.strip().lower()
        if r in STRING_INTERACTIONS:
            return r
        try:
            for alias in run_server.expand_gene_aliases(r):
                al = alias.lower()
                if al in STRING_INTERACTIONS:
                    return al
        except Exception:
            pass
        return r

    def display_name(locus: str) -> str:
        try:
            n = run_server.normalize_gene_locus(locus)
            return n if n and n != locus else locus
        except Exception:
            return locus

    src_id = resolve(source)
    tgt_id = resolve(target)

    if src_id == tgt_id:
        return {"found": False, "message": "Source and target are the same gene", "path": [], "edges": []}

    if src_id not in STRING_INTERACTIONS:
        return {"found": False, "message": f"No PPI data for source '{source}'", "path": [], "edges": []}
    if tgt_id not in STRING_INTERACTIONS:
        return {"found": False, "message": f"No PPI data for target '{target}'", "path": [], "edges": []}

    queue = deque([[src_id]])
    visited = {src_id}

    while queue:
        path = queue.popleft()
        current = path[-1]
        if len(path) > max_hops + 1:
            break
        partners = [p for p in STRING_INTERACTIONS.get(current, []) if p.get("score", 0) >= min_score]
        for p in partners:
            pid = p["partner"].lower()
            if pid in visited:
                continue
            new_path = path + [pid]
            if pid == tgt_id:
                nodes = [{"id": n, "name": display_name(n), "is_seed": n in (src_id, tgt_id)} for n in new_path]
                edges = []
                for i in range(len(new_path) - 1):
                    a, b = new_path[i], new_path[i + 1]
                    partners_a = STRING_INTERACTIONS.get(a, [])
                    edge_data = next((x for x in partners_a if x["partner"].lower() == b), None)
                    if edge_data:
                        edges.append({
                            "id": f"path-{a}-{b}",
                            "source": a,
                            "target": b,
                            "score": edge_data.get("score", 0),
                            "experimental": edge_data.get("experimental", 0),
                            "database": edge_data.get("database", 0),
                            "coexpression": edge_data.get("coexpression", 0),
                            "textmining": edge_data.get("textmining", 0),
                            "neighborhood": edge_data.get("neighborhood", 0),
                            "cooccurrence": edge_data.get("cooccurrence", 0),
                            "fusion": edge_data.get("fusion", 0),
                        })
                return {
                    "found": True,
                    "hops": len(new_path) - 1,
                    "path_ids": new_path,
                    "nodes": nodes,
                    "edges": edges,
                    "source": src_id,
                    "target": tgt_id,
                }
            visited.add(pid)
            queue.append(new_path)

    return {"found": False, "message": f"No path found within {max_hops} hops", "path": [], "edges": []}


@router.get("/api/analysis/string_ppi/hub_ranking")
def get_ppi_hub_ranking(min_score: int = 700, limit: int = 50):
    def display_name(locus: str) -> str:
        try:
            n = run_server.normalize_gene_locus(locus)
            return n if n and n != locus else locus
        except Exception:
            return locus

    rows = []
    for gene, partners in STRING_INTERACTIONS.items():
        if gene == "_meta":
            continue
        filtered = [p for p in partners if p.get("score", 0) >= min_score]
        if not filtered:
            continue
        def avg_ch(ch):
            vals = [p.get(ch, 0) for p in filtered]
            return round(sum(vals) / len(vals)) if vals else 0
        rows.append({
            "gene": gene,
            "name": display_name(gene),
            "degree": len(filtered),
            "avg_score": round(sum(p.get("score", 0) for p in filtered) / len(filtered)),
            "experimental": avg_ch("experimental"),
            "database": avg_ch("database"),
            "coexpression": avg_ch("coexpression"),
            "textmining": avg_ch("textmining"),
            "top_partners": [p["partner"] for p in sorted(filtered, key=lambda x: x.get("score", 0), reverse=True)[:5]],
        })

    rows.sort(key=lambda x: x["degree"], reverse=True)
    return {
        "min_score": min_score,
        "total_genes": len(rows),
        "hubs": rows[:limit],
    }


@router.get("/api/analysis/network_ppi_edges")
def get_network_ppi_edges(genes: str = "", min_score: int = 400):
    if not genes:
        return {"edges": []}

    gene_list = [g.strip().lower() for g in genes.split(",") if g.strip()]
    gene_set = set(gene_list)

    resolved_set = set()
    for g in gene_set:
        resolved_set.add(g)
        try:
            for alias in run_server.expand_gene_aliases(g):
                resolved_set.add(alias.lower())
        except Exception:
            pass

    edges = []
    seen_pairs = set()
    for g in resolved_set:
        partners = STRING_INTERACTIONS.get(g, [])
        for p in partners:
            partner = p["partner"].lower()
            if partner in resolved_set:
                score = p.get("score", 0)
                if score >= min_score:
                    pair = (min(g, partner), max(g, partner))
                    if pair not in seen_pairs:
                        seen_pairs.add(pair)
                        edges.append({
                            "source": pair[0],
                            "target": pair[1],
                            "score": score,
                            "experimental": p.get("experimental", 0),
                            "database": p.get("database", 0),
                            "coexpression": p.get("coexpression", 0),
                            "neighborhood": p.get("neighborhood", 0),
                            "cooccurrence": p.get("cooccurrence", 0),
                            "textmining": p.get("textmining", 0),
                            "fusion": p.get("fusion", 0),
                        })
    return {"edges": edges}


@router.get("/api/analysis/cross_motifs")
def get_cross_motifs(motif_type: str = "co_complex", min_score: int = 400):
    import csv

    path = run_server.get_absolute_path("data/reference/regulations.csv")
    if not os.path.exists(path):
        raise HTTPException(status_code=404, detail="regulations.csv not found")

    tf_targets = {}
    tf_names = {}
    tg_names = {}
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            tf = (row.get("TF_locusTag") or "").strip().lower()
            tf_name = (row.get("TF_name") or "").strip()
            tg = (row.get("TG_locusTag") or "").strip().lower()
            tg_name = (row.get("TG_name") or "").strip()
            if not tf or not tg:
                continue
            tf_targets.setdefault(tf, set()).add(tg)
            if tf_name and tf_name != tf:
                tf_names[tf] = tf_name
            if tg_name and tg_name != tg:
                tg_names[tg] = tg_name

    def display_name(locus: str) -> str:
        if locus in tf_names:
            return tf_names[locus]
        if locus in tg_names:
            return tg_names[locus]
        try:
            n = run_server.normalize_gene_locus(locus)
            return n if n and n != locus else locus
        except Exception:
            return locus

    ppi_edges = {}
    for gene, partners in STRING_INTERACTIONS.items():
        if gene == "_meta":
            continue
        g = gene.lower()
        for p in partners:
            partner = p["partner"].lower()
            score = p.get("score", 0)
            if score >= min_score:
                pair = (min(g, partner), max(g, partner))
                if pair not in ppi_edges or score > ppi_edges[pair]:
                    ppi_edges[pair] = score

    motifs = []

    if motif_type == "co_complex":
        for tf, targets in tf_targets.items():
            targets_list = sorted(list(targets))
            n_targets = len(targets_list)
            for i in range(n_targets):
                b = targets_list[i]
                for j in range(i + 1, n_targets):
                    c = targets_list[j]
                    pair = (min(b, c), max(b, c))
                    if pair in ppi_edges:
                        motifs.append({
                            "tf": tf,
                            "tf_name": display_name(tf),
                            "target_b": b,
                            "target_b_name": display_name(b),
                            "target_c": c,
                            "target_c_name": display_name(c),
                            "ppi_score": ppi_edges[pair],
                        })

    elif motif_type == "co_tf":
        tfs = sorted(list(tf_targets.keys()))
        n_tfs = len(tfs)
        for i in range(n_tfs):
            tf_a = tfs[i]
            for j in range(i + 1, n_tfs):
                tf_b = tfs[j]
                pair = (min(tf_a, tf_b), max(tf_a, tf_b))
                if pair in ppi_edges:
                    shared = tf_targets[tf_a].intersection(tf_targets[tf_b])
                    for c in sorted(list(shared)):
                        motifs.append({
                            "tf_a": tf_a,
                            "tf_a_name": display_name(tf_a),
                            "tf_b": tf_b,
                            "tf_b_name": display_name(tf_b),
                            "target_c": c,
                            "target_c_name": display_name(c),
                            "ppi_score": ppi_edges[pair],
                        })

    elif motif_type == "feedback":
        for tf, targets in tf_targets.items():
            for b in sorted(list(targets)):
                pair = (min(tf, b), max(tf, b))
                if pair in ppi_edges:
                    motifs.append({
                        "tf": tf,
                        "tf_name": display_name(tf),
                        "target": b,
                        "target_name": display_name(b),
                        "ppi_score": ppi_edges[pair],
                    })

    motifs.sort(key=lambda x: x.get("ppi_score", 0), reverse=True)
    return {
        "motif_type": motif_type,
        "min_score": min_score,
        "count": len(motifs),
        "instances": motifs[:200],
    }


@router.get("/api/graph/cascade", response_model=GraphCascadeResponse)
def get_graph_cascade(source: str = "", target: str = "", max_depth: int = 3):
    if not source or not target:
        raise HTTPException(status_code=400, detail="source and target parameters are required")
    ge = get_graph_engine() if get_graph_engine else None
    if not ge:
        raise HTTPException(status_code=503, detail="Graph engine unavailable")
    paths = ge.find_cascade_paths(source=source, target=target, max_depth=max_depth)
    return GraphCascadeResponse(source=source, target=target, max_depth=max_depth, paths=paths)


@router.get("/api/graph/motifs", response_model=GraphMotifResponse)
def get_graph_motifs(type: str = "ffl", limit: int = 50):
    ge = get_graph_engine() if get_graph_engine else None
    if not ge:
        raise HTTPException(status_code=503, detail="Graph engine unavailable")
    items = ge.detect_motifs(motif_type=type, limit=limit)
    return GraphMotifResponse(motif_type=type, count=len(items), items=items)


@router.get("/api/tf/effectors/{tf_id}")
def get_tf_effectors_api(tf_id: str):
    db = get_db_manager()
    res = db.get_tf_effector_info(tf_id)
    if not res:
        raise HTTPException(status_code=404, detail=f"TF Effector info not found for {tf_id}")
    return res


@router.get("/api/network/extended")
def get_extended_network_api(locus: str = "", mode: str = "all", edge_type: str = None):
    db = get_db_manager()
    edges = db.get_extended_edges(locus, mode=mode, edge_type=edge_type)
    return {"locus": locus, "mode": mode, "edge_type": edge_type, "count": len(edges), "edges": edges}


@router.get("/api/network/allosteric-feedback")
def get_allosteric_feedback_api(query: str = None):
    db = get_db_manager()
    loops = db.get_allosteric_feedback_loops(tf_or_metabolite=query)
    return {"query": query, "count": len(loops), "loops": loops}


@router.get("/api/network/srna-competition")
def get_srna_competition_api(srna_id: str = None):
    db = get_db_manager()
    targets = db.get_srna_target_competition(srna_id=srna_id)
    return {"srna_id": srna_id, "count": len(targets), "targets": targets}


@router.get("/api/imodulon/gene/{gene_id}")
def get_imodulon_gene_api(gene_id: str):
    db = get_db_manager()
    imodulons = db.get_imodulons_for_gene(gene_id)
    return {"gene_id": gene_id, "count": len(imodulons), "imodulons": imodulons}


@router.get("/api/network/rf-scores")
def get_rf_scores_api(locus: str = "", min_confidence: float = 0.3):
    db = get_db_manager()
    scores = db.get_rf_edge_scores(locus, min_confidence=min_confidence)
    return {"locus": locus, "min_confidence": min_confidence, "count": len(scores), "scores": scores}


@router.get("/api/tf/hierarchy-rankings")
def get_tf_hierarchy_rankings_api():
    db = get_db_manager()
    rankings = db.get_tf_hierarchy_rankings()
    return {"count": len(rankings), "rankings": rankings}


@router.get("/api/network/rewired")
def get_network_rewired_api(locus: str = None):
    db = get_db_manager()
    edges = db.get_rewired_edges(locus=locus)
    return {"locus": locus, "count": len(edges), "edges": edges}


@router.get("/api/tfbs/collectf")
def get_collectf_tfbs_api(locus: str = None):
    db = get_db_manager()
    sites = db.get_collectf_tfbs(locus=locus)
    return {"locus": locus, "count": len(sites), "sites": sites}


@router.get("/api/pathway/gene/{gene_id}")
def get_pathways_for_gene_api(gene_id: str):
    db = get_db_manager()
    pathways = db.get_pathways_for_gene(gene_id)
    return {"gene_id": gene_id, "count": len(pathways), "pathways": pathways}


@router.get("/api/pathway/info/{pathway_id}")
def get_genes_in_pathway_api(pathway_id: str):
    db = get_db_manager()
    genes = db.get_genes_in_pathway(pathway_id)
    return {"pathway_id": pathway_id, "count": len(genes), "genes": genes}


@router.get("/api/imodulon/condition")
def get_condition_specific_regulons_api(condition: str = None):
    db = get_db_manager()
    activities = db.get_condition_specific_regulons(condition_name=condition)
    return {"condition": condition, "count": len(activities), "activities": activities}


@router.get("/api/imodulon/overlap")
def get_imodulon_regulon_overlap_api(imodulon_id: str = None):
    db = get_db_manager()
    overlaps = db.get_imodulon_regulon_overlap(imodulon_id=imodulon_id)
    return {"imodulon_id": imodulon_id, "count": len(overlaps), "overlaps": overlaps}


@router.get("/api/condition-regulation/runs")
def get_condition_regulation_runs_api():
    db = get_db_manager()
    runs = db.get_condition_regulation_runs()
    return {"count": len(runs), "runs": runs}


@router.get("/api/condition-regulation/conditions")
def get_condition_regulation_conditions_api(run_id: str = "iron_regulon_v1"):
    db = get_db_manager()
    conditions = db.get_condition_regulation_conditions(run_id=run_id)
    return {"run_id": run_id, "count": len(conditions), "conditions": conditions}


@router.get("/api/condition-regulation/summary")
def get_condition_regulation_summary_api(
    comparison_id: str = None, tf: str = None,
    run_id: str = "iron_regulon_v1",
):
    db = get_db_manager()
    summaries = db.get_condition_regulation_summary(
        comparison_id=comparison_id, tf_name=tf, run_id=run_id,
    )
    return {
        "run_id": run_id,
        "comparison_id": comparison_id,
        "tf": tf,
        "count": len(summaries),
        "summaries": summaries,
    }


@router.get("/api/condition-regulation/edges")
def get_condition_regulation_edges_api(
    comparison_id: str,
    tf: str = None,
    state: str = None,
    min_score: float = 0.0,
    limit: int = 100,
    offset: int = 0,
    run_id: str = "iron_regulon_v1",
):
    if not 0.0 <= min_score <= 1.0:
        raise HTTPException(status_code=400, detail="min_score must be between 0 and 1")
    db = get_db_manager()
    result = db.get_condition_regulation_edges(
        comparison_id=comparison_id,
        tf_name=tf,
        support_state=state,
        min_score=min_score,
        limit=limit,
        offset=offset,
        run_id=run_id,
    )
    return {
        "run_id": run_id,
        "comparison_id": comparison_id,
        "tf": tf,
        "state": state,
        "min_score": min_score,
        "total": result["total"],
        "count": len(result["edges"]),
        "edges": result["edges"],
    }


@router.get(
    "/api/intervention-targets",
    response_model=InterventionTargetsResponse,
    responses={400: {"model": HTTPError}}
)
def get_intervention_targets_api(
    q: str = None,
    strategy: str = None,
    min_modules: int = 1,
    max_risk: float = 1.0,
    grade: str = None,
    include_known_essential: bool = True,
    limit: int = 100,
    offset: int = 0,
):
    if not 0.0 <= max_risk <= 1.0:
        raise HTTPException(status_code=400, detail="max_risk must be between 0 and 1")
    db = get_db_manager()
    result = db.get_intervention_targets(
        query=q, strategy=strategy, min_modules=min_modules, max_risk=max_risk,
        evidence_grade=grade, include_known_essential=include_known_essential,
        limit=limit, offset=offset,
    )
    return {"total": result["total"], "limit": limit, "targets": result["targets"]}


@router.get(
    "/api/intervention-targets/{locus}",
    response_model=InterventionTargetSchema,
    responses={404: {"model": HTTPError}}
)
def get_intervention_target_detail_api(locus: str):
    db = get_db_manager()
    target = db.get_intervention_target_detail(locus)
    if not target:
        raise HTTPException(status_code=404, detail="Target priority not found")
    return target
