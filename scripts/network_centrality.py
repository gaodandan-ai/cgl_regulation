"""
scripts/network_centrality.py
================================
Computes network centrality metrics for the C. glutamicum regulatory network:
  - Out-degree (regulon size)
  - In-degree (how many TFs regulate this gene)
  - Betweenness centrality (how often on shortest paths)
  - PageRank (authority score)
  - Closeness centrality
  - Hub score (HITS algorithm)
  - Activation ratio (fraction of activating edges)
  - Evidence quality score

Output: data/reference/network_centrality.json
"""

import os, json, math, logging
import pandas as pd
import networkx as nx
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("network_centrality")

ROOT_DIR  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REG_CSV   = os.path.join(ROOT_DIR, "data", "reference", "regulations.csv")
TCS_CSV   = os.path.join(ROOT_DIR, "data", "reference", "tcs_regulations.csv")
OUT_PATH  = os.path.join(ROOT_DIR, "data", "reference", "network_centrality.json")


def evidence_score(ev_str: str) -> float:
    """Map evidence type to quality weight [0,1]."""
    ev = str(ev_str).lower()
    if "experimental + predicted" in ev:
        return 1.0
    if "experimental" in ev:
        return 0.85
    if "predicted" in ev:
        return 0.35
    return 0.2


def main():
    log.info("Loading regulatory network...")
    df = pd.read_csv(REG_CSV)

    # Optional: merge TCS regulations
    if os.path.exists(TCS_CSV):
        try:
            tcs_df = pd.read_csv(TCS_CSV)
            if set(["TF_locusTag", "TG_locusTag"]).issubset(tcs_df.columns):
                # Add missing columns with defaults
                for col in df.columns:
                    if col not in tcs_df.columns:
                        tcs_df[col] = None
                df = pd.concat([df, tcs_df[df.columns]], ignore_index=True)
                log.info(f"Merged TCS regulations. Total edges: {len(df)}")
        except Exception as e:
            log.warning(f"Could not merge TCS data: {e}")

    # --- Build directed graph ---
    G = nx.DiGraph()

    tf_meta    = {}   # TF metadata
    gene_meta  = {}   # all gene metadata
    tf_act_cnt = defaultdict(int)
    tf_rep_cnt = defaultdict(int)
    tf_ev_sum  = defaultdict(float)
    tf_ev_cnt  = defaultdict(int)

    for _, row in df.iterrows():
        tf   = str(row["TF_locusTag"]).strip()
        tg   = str(row["TG_locusTag"]).strip()
        role = str(row.get("Role", "")).strip().upper()
        ev   = row.get("Evidence", "predicted")
        ev_w = evidence_score(ev)

        # Skip self-loops for cleaner betweenness
        if tf == tg:
            continue

        G.add_edge(tf, tg, role=role, evidence=ev, weight=ev_w)

        if role == "A":
            tf_act_cnt[tf] += 1
        elif role == "R":
            tf_rep_cnt[tf] += 1

        tf_ev_sum[tf] += ev_w
        tf_ev_cnt[tf]  += 1

        # Collect metadata
        if tf not in tf_meta:
            tf_meta[tf] = {
                "locus": tf,
                "alt_locus": str(row.get("TF_altLocusTag", "")) if pd.notna(row.get("TF_altLocusTag")) else "",
                "name": str(row.get("TF_name", tf)).strip() if pd.notna(row.get("TF_name")) else tf,
                "is_sigma": bool(row.get("Is_sigma_factor", False)),
            }

    log.info(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # --- Compute centrality metrics ---
    log.info("Computing out-degree (regulon size)...")
    out_degree = dict(G.out_degree())

    log.info("Computing in-degree...")
    in_degree  = dict(G.in_degree())

    log.info("Computing betweenness centrality (this may take ~10s)...")
    betweenness = nx.betweenness_centrality(G, normalized=True, weight=None)

    log.info("Computing PageRank...")
    try:
        pagerank = nx.pagerank(G, alpha=0.85, weight="weight", max_iter=200)
    except nx.PowerIterationFailedConvergence:
        pagerank = nx.pagerank(G, alpha=0.85, weight=None, max_iter=500)

    log.info("Computing HITS (hub/authority scores)...")
    try:
        hub_scores, auth_scores = nx.hits(G, max_iter=200, normalized=True)
    except nx.PowerIterationFailedConvergence:
        hub_scores  = {n: 0.0 for n in G.nodes()}
        auth_scores = {n: 0.0 for n in G.nodes()}

    log.info("Computing closeness centrality...")
    # Use undirected version for closeness (connectivity measure)
    G_undir = G.to_undirected()
    closeness = nx.closeness_centrality(G_undir)

    # --- Compute composite TF importance score ---
    # Formula: 0.30*norm(out_degree) + 0.25*norm(betweenness) + 0.20*norm(pagerank)
    #          + 0.10*ev_quality + 0.15*norm(hub)
    log.info("Computing composite importance scores...")

    def safe_norm(d: dict) -> dict:
        """Normalize dict values to [0, 1]."""
        vals = list(d.values())
        vmin, vmax = min(vals), max(vals)
        if vmax == vmin:
            return {k: 0.5 for k in d}
        return {k: (v - vmin) / (vmax - vmin) for k, v in d.items()}

    norm_out    = safe_norm(out_degree)
    norm_btw    = safe_norm(betweenness)
    norm_pr     = safe_norm(pagerank)
    norm_hub    = safe_norm(hub_scores)
    norm_close  = safe_norm(closeness)

    all_nodes   = list(G.nodes())
    tf_nodes    = list(tf_meta.keys())

    composite   = {}
    for n in all_nodes:
        ev_q = (tf_ev_sum.get(n, 0) / max(tf_ev_cnt.get(n, 1), 1))  # [0,1]
        composite[n] = (
            0.30 * norm_out.get(n, 0)
          + 0.25 * norm_btw.get(n, 0)
          + 0.20 * norm_pr.get(n, 0)
          + 0.10 * ev_q
          + 0.15 * norm_hub.get(n, 0)
        )

    # --- Assemble output ---
    result_nodes = {}
    for n in all_nodes:
        act  = tf_act_cnt.get(n, 0)
        rep  = tf_rep_cnt.get(n, 0)
        total = act + rep
        result_nodes[n] = {
            "locus":         n,
            "name":          tf_meta.get(n, {}).get("name", n),
            "is_tf":         n in tf_meta,
            "is_sigma":      tf_meta.get(n, {}).get("is_sigma", False),
            "out_degree":    int(out_degree.get(n, 0)),
            "in_degree":     int(in_degree.get(n, 0)),
            "betweenness":   round(betweenness.get(n, 0.0), 6),
            "pagerank":      round(pagerank.get(n, 0.0), 6),
            "hub_score":     round(hub_scores.get(n, 0.0), 6),
            "auth_score":    round(auth_scores.get(n, 0.0), 6),
            "closeness":     round(closeness.get(n, 0.0), 6),
            "importance":    round(composite.get(n, 0.0), 4),
            "n_activations": int(act),
            "n_repressions": int(rep),
            "activation_ratio": round(act / total, 3) if total > 0 else 0.0,
            "ev_quality":    round(tf_ev_sum.get(n, 0) / max(tf_ev_cnt.get(n, 1), 1), 3),
        }

    # Rank TFs by composite importance
    tf_ranked = sorted(
        [v for v in result_nodes.values() if v["is_tf"]],
        key=lambda x: x["importance"],
        reverse=True
    )
    for rank, tf in enumerate(tf_ranked, 1):
        result_nodes[tf["locus"]]["rank"] = rank

    output = {
        "_meta": {
            "n_nodes":       G.number_of_nodes(),
            "n_edges":       G.number_of_edges(),
            "n_tfs":         len(tf_nodes),
            "n_target_genes": len(set(all_nodes) - set(tf_nodes)),
            "top_tf":        tf_ranked[0]["locus"] if tf_ranked else None,
        },
        "nodes": result_nodes,
        "top_tfs_by_importance": [
            {"locus": tf["locus"], "name": tf["name"], "importance": tf["importance"],
             "out_degree": tf["out_degree"], "betweenness": tf["betweenness"],
             "pagerank": tf["pagerank"], "activation_ratio": tf["activation_ratio"]}
            for tf in tf_ranked[:20]
        ],
    }

    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)
    log.info(f"Saved: {OUT_PATH}")

    # Print top 20 TFs
    print("\n" + "="*80)
    print(f"{'Rank':>4}  {'TF':>10}  {'Name':<12}  {'OutDeg':>6}  {'Betw':>8}  {'PageR':>8}  {'Score':>7}")
    print("-"*80)
    for i, tf in enumerate(tf_ranked[:20], 1):
        print(f"{i:>4}  {tf['locus']:>10}  {tf['name']:<12}  {tf['out_degree']:>6}  "
              f"{tf['betweenness']:.4f}  {tf['pagerank']:.4f}  {tf['importance']:.4f}")
    print("="*80)


if __name__ == "__main__":
    main()
