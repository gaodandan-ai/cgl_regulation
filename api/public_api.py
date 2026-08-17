"""Lightweight public endpoints that do not bundle the desktop FBA stack."""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import JSONResponse

from backend.db_manager import get_db_manager
from backend.services.provenance import build_provenance


LOGGER = logging.getLogger("cgl.public_api")
ROOT = Path(__file__).resolve().parents[1]


def _version() -> str:
    try:
        return str(json.loads((ROOT / "web" / "version.json").read_text(encoding="utf-8"))["version"])
    except (OSError, KeyError, TypeError, ValueError):
        return "0.0.0"


def _db():
    manager = get_db_manager()
    if manager.get_connection() is None:
        raise HTTPException(status_code=503, detail="Public data service is unavailable")
    return manager


def configure_public_api(app: FastAPI) -> None:
    origins = [
        item.strip()
        for item in os.environ.get(
            "CGL_ALLOWED_ORIGINS", "https://cgl-regulation.vercel.app"
        ).split(",")
        if item.strip()
    ]
    app.add_middleware(
        CORSMiddleware,
        allow_origins=origins,
        allow_credentials=False,
        allow_methods=["GET", "OPTIONS"],
        allow_headers=["Accept", "Content-Type", "If-None-Match"],
    )
    app.add_middleware(GZipMiddleware, minimum_size=1000)

    @app.middleware("http")
    async def security_headers(request: Request, call_next):
        response = await call_next(request)
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["X-Frame-Options"] = "DENY"
        response.headers["Referrer-Policy"] = "strict-origin-when-cross-origin"
        response.headers["Permissions-Policy"] = "camera=(), microphone=(), geolocation=()"
        response.headers["Content-Security-Policy"] = (
            "default-src 'self'; object-src 'none'; base-uri 'self'; frame-ancestors 'none'"
        )
        return response

    @app.exception_handler(Exception)
    async def unhandled_exception(request: Request, exc: Exception):
        LOGGER.exception("Public API error on %s", request.url.path)
        return JSONResponse(status_code=500, content={"detail": "Internal server error"})

    @app.get("/api/health")
    def health():
        manager = get_db_manager()
        db_path = getattr(manager, "_db_path", "")
        available = bool(db_path and os.path.isfile(db_path))
        return {
            "status": "ok" if available else "degraded",
            "app": "cgl-regulation",
            "version": _version(),
            "database": "available" if available else "unavailable",
            "deployment": "public",
            "capabilities": {
                "regulatory_data": available,
                "condition_regulation": available,
                "target_priorities": available,
                "metabolic_simulation": False,
                "desktop_application": True,
            },
        }

    @app.get("/api/provenance")
    def provenance():
        return build_provenance(
            manager=get_db_manager(), root=ROOT, version=_version(), deployment="public"
        )

    @app.get("/api/gene/coordinates/{gene_id}")
    def gene_coordinates(gene_id: str):
        result = _db().get_gene_coordinates(gene_id)
        if not result:
            raise HTTPException(status_code=404, detail="Gene coordinates not found")
        return result

    @app.get("/api/tf/effectors/{tf_id}")
    def tf_effectors(tf_id: str):
        result = _db().get_tf_effector_info(tf_id)
        if not result:
            raise HTTPException(status_code=404, detail="TF effector information not found")
        return result

    @app.get("/api/network/extended")
    def extended_network(locus: str = "", mode: str = "all", edge_type: str = None):
        edges = _db().get_extended_edges(locus, mode=mode, edge_type=edge_type)
        return {"locus": locus, "mode": mode, "edge_type": edge_type, "count": len(edges), "edges": edges}

    @app.get("/api/gene/profile/{gene_id}")
    def gene_profile(gene_id: str):
        result = _db().get_full_gene_profile(gene_id)
        if not result:
            raise HTTPException(status_code=404, detail="Gene profile not found")
        return result

    @app.get("/api/gene/neighborhood/{gene_id}")
    def gene_neighborhood(gene_id: str, window_bp: int = 20000):
        genes = _db().get_genomic_neighborhood(gene_id, window_bp=window_bp)
        return {"center_gene": gene_id, "window_bp": window_bp, "count": len(genes), "genes": genes}

    @app.get("/api/network/allosteric-feedback")
    def allosteric_feedback(query: str = None):
        rows = _db().get_allosteric_feedback_loops(query)
        return {"query": query, "count": len(rows), "loops": rows}

    @app.get("/api/network/srna-competition")
    def srna_competition(srna_id: str = None):
        rows = _db().get_srna_target_competition(srna_id)
        return {"srna_id": srna_id, "count": len(rows), "targets": rows}

    @app.get("/api/imodulon/gene/{gene_id}")
    def imodulon_gene(gene_id: str):
        rows = _db().get_imodulons_for_gene(gene_id)
        return {"gene_id": gene_id, "count": len(rows), "imodulons": rows}

    @app.get("/api/network/rf-scores")
    def rf_scores(locus: str = "", min_confidence: float = 0.3):
        rows = _db().get_rf_edge_scores(locus, min_confidence)
        return {"locus": locus, "min_confidence": min_confidence, "count": len(rows), "scores": rows}

    @app.get("/api/tf/hierarchy-rankings")
    def hierarchy_rankings():
        rows = _db().get_tf_hierarchy_rankings()
        return {"count": len(rows), "rankings": rows}

    @app.get("/api/network/rewired")
    def rewired_network(locus: str = None):
        rows = _db().get_rewired_edges(locus)
        return {"locus": locus, "count": len(rows), "edges": rows}

    @app.get("/api/tfbs/collectf")
    def collectf_tfbs(locus: str = None):
        rows = _db().get_collectf_tfbs(locus)
        return {"locus": locus, "count": len(rows), "sites": rows}

    @app.get("/api/pathway/gene/{gene_id}")
    def gene_pathways(gene_id: str):
        rows = _db().get_pathways_for_gene(gene_id)
        return {"gene_id": gene_id, "count": len(rows), "pathways": rows}

    @app.get("/api/pathway/info/{pathway_id}")
    def pathway_genes(pathway_id: str):
        rows = _db().get_genes_in_pathway(pathway_id)
        return {"pathway_id": pathway_id, "count": len(rows), "genes": rows}

    @app.get("/api/ncrna/list")
    def ncrna_list(rna_type: str = None):
        rows = _db().get_ncrnas(rna_type)
        return {"rna_type": rna_type, "count": len(rows), "ncrnas": rows}

    @app.get("/api/ncrna/targets")
    def ncrna_targets(locus: str = None):
        rows = _db().get_srna_targets(locus)
        return {"locus": locus, "count": len(rows), "targets": rows}

    @app.get("/api/imodulon/condition")
    def imodulon_condition(condition: str = None):
        rows = _db().get_condition_specific_regulons(condition)
        return {"condition": condition, "count": len(rows), "activities": rows}

    @app.get("/api/imodulon/overlap")
    def imodulon_overlap(imodulon_id: str = None):
        rows = _db().get_imodulon_regulon_overlap(imodulon_id)
        return {"imodulon_id": imodulon_id, "count": len(rows), "overlaps": rows}

    @app.get("/api/condition-regulation/runs")
    def condition_runs():
        rows = _db().get_condition_regulation_runs()
        return {"count": len(rows), "runs": rows}

    @app.get("/api/condition-regulation/conditions")
    def condition_conditions(run_id: str = "iron_regulon_v1"):
        rows = _db().get_condition_regulation_conditions(run_id)
        return {"run_id": run_id, "count": len(rows), "conditions": rows}

    @app.get("/api/condition-regulation/summary")
    def condition_summary(comparison_id: str = None, tf: str = None, run_id: str = "iron_regulon_v1"):
        rows = _db().get_condition_regulation_summary(comparison_id, tf, run_id)
        return {"run_id": run_id, "comparison_id": comparison_id, "tf": tf, "count": len(rows), "summaries": rows}

    @app.get("/api/condition-regulation/edges")
    def condition_edges(
        comparison_id: str, tf: str = None, state: str = None,
        min_score: float = 0.0, limit: int = 100, offset: int = 0,
        run_id: str = "iron_regulon_v1",
    ):
        if not 0 <= min_score <= 1:
            raise HTTPException(status_code=400, detail="min_score must be between 0 and 1")
        result = _db().get_condition_regulation_edges(
            comparison_id, tf, state, min_score, limit, offset, run_id
        )
        return {
            "run_id": run_id, "comparison_id": comparison_id, "tf": tf,
            "state": state, "min_score": min_score, "total": result["total"],
            "count": len(result["edges"]), "edges": result["edges"],
        }

    @app.get("/api/intervention-targets")
    def intervention_targets(
        q: str = None, strategy: str = None, min_modules: int = 1,
        max_risk: float = 1.0, grade: str = None,
        include_known_essential: bool = True, limit: int = 100, offset: int = 0,
    ):
        if not 0 <= max_risk <= 1:
            raise HTTPException(status_code=400, detail="max_risk must be between 0 and 1")
        result = _db().get_intervention_targets(
            q, strategy, min_modules, max_risk, grade,
            include_known_essential, limit, offset,
        )
        return {"total": result["total"], "count": len(result["targets"]), "targets": result["targets"]}

    @app.get("/api/intervention-targets/{locus}")
    def intervention_target(locus: str):
        result = _db().get_intervention_target_detail(locus)
        if not result:
            raise HTTPException(status_code=404, detail="Target priority not found")
        return result

    @app.get("/api/{path:path}")
    def desktop_only(path: str):
        raise HTTPException(
            status_code=503,
            detail="This analysis requires the desktop application",
        )
