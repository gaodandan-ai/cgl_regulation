import os
import sys
import logging
from collections import deque
from typing import List, Dict, Any, Optional

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
    import run_server
except ImportError:
    import backend.run_server as run_server

logger = logging.getLogger("app.services.graph_rag")


class GraphRAGService:
    def __init__(self):
        self._graph_cache = None

    def get_causal_paths(self, source: str, target: str, max_depth: int = 3) -> List[Dict[str, Any]]:
        """
        Find multi-hop causal regulatory & metabolic paths between source and target entities.
        """
        if not source or not target:
            return []

        src_clean = source.strip().lower()
        tgt_clean = target.strip().lower()

        # Expand aliases to canonical locus tags
        try:
            src_locus = run_server.normalize_gene_locus(src_clean) or src_clean
            tgt_locus = run_server.normalize_gene_locus(tgt_clean) or tgt_clean
        except Exception:
            src_locus = src_clean
            tgt_locus = tgt_clean

        paths = []

        if _DB_MANAGER_AVAILABLE and get_db_manager():
            db = get_db_manager()

            # 1. Try DB extended network path search
            try:
                edges_src = db.get_extended_edges(src_locus, mode="all")
                # Direct edge check
                for e in edges_src:
                    t_id = (e.get("target") or "").lower()
                    t_alias = (e.get("target_locus") or "").lower()
                    if tgt_locus in (t_id, t_alias) or tgt_clean in (t_id, t_alias):
                        paths.append({
                            "hops": 1,
                            "nodes": [src_locus, tgt_locus],
                            "mechanism": e.get("edge_type", "regulatory"),
                            "details": [f"{src_locus} --[{e.get('interaction_type', 'regulates')}]--> {tgt_locus}"]
                        })
            except Exception as ex:
                logger.warning(f"GraphRAG direct path query failed: {ex}")

        # 2. Try graph engine if available
        try:
            from graph_engine import get_graph_engine
            ge = get_graph_engine()
            if ge:
                ge_paths = ge.find_cascade_paths(source=src_locus, target=tgt_locus, max_depth=max_depth)
                for p in ge_paths:
                    nodes = p.get("nodes", [])
                    if len(nodes) > 1:
                        paths.append({
                            "hops": len(nodes) - 1,
                            "nodes": nodes,
                            "mechanism": "cascade_regulatory",
                            "details": [f"{nodes[i]} --> {nodes[i+1]}" for i in range(len(nodes)-1)]
                        })
        except Exception as ex:
            logger.warning(f"GraphRAG graph_engine cascade failed: {ex}")

        # Deduplicate paths
        seen = set()
        unique_paths = []
        for p in paths:
            key = tuple(p["nodes"])
            if key not in seen:
                seen.add(key)
                unique_paths.append(p)

        return unique_paths

    def query_graph_rag_reasoning(
        self,
        source: str = "",
        target: str = "",
        query: str = "",
        max_depth: int = 3
    ) -> Dict[str, Any]:
        """
        Build a GraphRAG knowledge graph causal reasoning report combining graph paths with literature evidence.
        """
        from rag_service import RAGService
        rag_svc = RAGService()

        # If source/target not explicitly provided, try to extract from query text
        if not source and query:
            # Simple locus tag extraction regex (e.g. cg0001, cg0350)
            import re
            loci = re.findall(r'\bcg\d{4}\b', query.lower())
            if len(loci) >= 2:
                source, target = loci[0], loci[1]
            elif len(loci) == 1:
                source = loci[0]

        causal_paths = []
        if source and target:
            causal_paths = self.get_causal_paths(source=source, target=target, max_depth=max_depth)

        # Retrieve relevant literature evidence via RAG
        lit_query = f"{source} {target} {query}".strip()
        lit_evidence = rag_svc.query_similarity(query=lit_query, provider="google", api_key="", model_name="", base_url="", top_n=3)

        # Format GraphRAG prompt context snippet
        context_lines = []
        context_lines.append(f"### GraphRAG Causal Network Evidence for '{source}' -> '{target}'")
        if causal_paths:
            context_lines.append(f"Found {len(causal_paths)} causal graph path(s):")
            for idx, p in enumerate(causal_paths, 1):
                nodes_str = " -> ".join(p["nodes"])
                context_lines.append(f"  Path {idx} ({p['hops']} hops): {nodes_str} [{p['mechanism']}]")
        else:
            context_lines.append("No direct multi-hop graph path detected in SQLite regulatory topology.")

        if lit_evidence:
            context_lines.append("\n### Supporting Literature Extracts:")
            for item in lit_evidence:
                context_lines.append(f"- **Source [{item.get('file')}]** (Relevance: {item.get('score', 0):.2f}): {item.get('text')}")

        formatted_context = "\n".join(context_lines)

        return {
            "source": source,
            "target": target,
            "query": query,
            "causal_paths_count": len(causal_paths),
            "causal_paths": causal_paths,
            "literature_evidence_count": len(lit_evidence),
            "literature_evidence": lit_evidence,
            "graph_rag_context": formatted_context
        }


_graph_rag_service_instance = None

def get_graph_rag_service() -> GraphRAGService:
    global _graph_rag_service_instance
    if _graph_rag_service_instance is None:
        _graph_rag_service_instance = GraphRAGService()
    return _graph_rag_service_instance
