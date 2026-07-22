#!/usr/bin/env python3
"""
graph_engine.py
===============
High-performance NetworkX graph analysis engine over cgl_regulation.db.
Provides fast multi-hop regulatory cascade tracing, motif detection (FFL, Feedback),
and sub-graph centrality calculations.
"""

import os
import sys
import logging
import networkx as nx

logger = logging.getLogger("cgl_graph")

_backend_dir = os.path.dirname(os.path.abspath(__file__))
_root_dir = os.path.dirname(_backend_dir)
if _backend_dir not in sys.path:
    sys.path.insert(0, _backend_dir)
if _root_dir not in sys.path:
    sys.path.insert(0, _root_dir)

from db_manager import get_db_manager

class RegulatoryNetworkGraph:
    """
    In-memory NetworkX Graph wrapper for Corynebacterium glutamicum regulatory network.
    """
    _instance = None

    def __init__(self):
        self.graph = nx.DiGraph()
        self.is_loaded = False
        self.load_graph()

    @classmethod
    def get_instance(cls):
        if cls._instance is None:
            cls._instance = RegulatoryNetworkGraph()
        return cls._instance

    def load_graph(self):
        db = get_db_manager()
        conn = db.get_connection()
        if not conn:
            logger.warning("Database connection unavailable for Graph Engine.")
            return

        cursor = conn.cursor()
        cursor.execute("SELECT TF_locusTag, TF_name, TG_locusTag, TG_name, Role, evidence_score FROM regulations")
        rows = cursor.fetchall()

        self.graph.clear()
        for r in rows:
            tf_id = (r["TF_locusTag"] or r["TF_name"] or "").strip()
            tg_id = (r["TG_locusTag"] or r["TG_name"] or "").strip()
            if tf_id and tg_id:
                self.graph.add_node(tf_id.lower(), name=r["TF_name"] or tf_id, is_tf=True)
                self.graph.add_node(tg_id.lower(), name=r["TG_name"] or tg_id)
                self.graph.add_edge(
                    tf_id.lower(),
                    tg_id.lower(),
                    role=r["Role"] or "A",
                    score=float(r["evidence_score"] or 1.0)
                )

        self.is_loaded = True
        logger.info(f"Loaded Regulatory Graph: {self.graph.number_of_nodes()} nodes, {self.graph.number_of_edges()} edges.")

    def find_cascade_paths(self, source: str, target: str, max_depth: int = 3):
        """
        Find all directed paths from source to target up to max_depth.
        """
        if not self.is_loaded:
            self.load_graph()

        src = source.lower().strip()
        tgt = target.lower().strip()

        if src not in self.graph or tgt not in self.graph:
            return []

        try:
            paths = list(nx.all_simple_paths(self.graph, source=src, target=tgt, cutoff=max_depth))
            formatted_paths = []
            for path in paths[:20]: # Return top 20 paths
                path_details = []
                for i in range(len(path) - 1):
                    u, v = path[i], path[i+1]
                    edge_data = self.graph.get_edge_data(u, v, {})
                    path_details.append({
                        "from": u,
                        "to": v,
                        "role": edge_data.get("role", "A"),
                        "score": edge_data.get("score", 1.0)
                    })
                formatted_paths.append({
                    "nodes": path,
                    "length": len(path) - 1,
                    "edges": path_details
                })
            return formatted_paths
        except Exception as e:
            logger.error(f"Error finding cascade paths: {e}")
            return []

    def detect_motifs(self, motif_type: str = "ffl", limit: int = 50):
        """
        Detect structural regulatory motifs:
          - 'ffl': Feed-Forward Loop (A -> B, A -> C, B -> C)
          - 'feedback': Mutual / Feedback Loop (A -> B, B -> A)
        """
        if not self.is_loaded:
            self.load_graph()

        results = []
        g = self.graph

        if motif_type == "ffl":
            for a in g.nodes():
                succ_a = set(g.successors(a))
                for b in succ_a:
                    if a == b:
                        continue
                    succ_b = set(g.successors(b))
                    # Common targets regulated by both A and B
                    common_targets = succ_a.intersection(succ_b)
                    for c in common_targets:
                        if c != a and c != b:
                            results.append({
                                "motif": "Feed-Forward Loop (FFL)",
                                "tf_a": a,
                                "tf_b": b,
                                "target_c": c
                            })
                            if len(results) >= limit:
                                return results

        elif motif_type == "feedback":
            seen_pairs = set()
            for u, v in g.edges():
                if g.has_edge(v, u) and u != v:
                    pair = tuple(sorted([u, v]))
                    if pair not in seen_pairs:
                        seen_pairs.add(pair)
                        results.append({
                            "motif": "Feedback Loop",
                            "gene_a": u,
                            "gene_b": v
                        })
                        if len(results) >= limit:
                            return results

        return results

    def compute_subgraph_centrality(self, gene_list: list):
        """
        Compute centrality metrics for a given subset of genes.
        """
        if not self.is_loaded:
            self.load_graph()

        valid_genes = [g.lower().strip() for g in gene_list if g.lower().strip() in self.graph]
        if not valid_genes:
            return {}

        subgraph = self.graph.subgraph(valid_genes).copy()
        if subgraph.number_of_nodes() == 0:
            return {}

        deg = dict(subgraph.degree())
        in_deg = dict(subgraph.in_degree())
        out_deg = dict(subgraph.out_degree())

        try:
            between = nx.betweenness_centrality(subgraph)
        except Exception:
            between = {n: 0.0 for n in subgraph.nodes()}

        try:
            pr = nx.pagerank(subgraph, max_iter=200)
        except Exception:
            pr = {n: 0.0 for n in subgraph.nodes()}

        metrics = {}
        for n in valid_genes:
            metrics[n] = {
                "degree": deg.get(n, 0),
                "in_degree": in_deg.get(n, 0),
                "out_degree": out_deg.get(n, 0),
                "betweenness": round(between.get(n, 0.0), 4),
                "pagerank": round(pr.get(n, 0.0), 4)
            }
        return metrics

# Export helper functions
def get_graph_engine():
    return RegulatoryNetworkGraph.get_instance()
