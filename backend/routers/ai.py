from fastapi import APIRouter, HTTPException, Header
from pydantic import BaseModel
import os
import sys
import logging
import traceback

backend_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if backend_dir not in sys.path:
    sys.path.insert(0, backend_dir)

try:
    from ai_handlers import (
        perform_summarize,
        perform_pathway_analysis,
        perform_gene_analysis,
        perform_protein_domain_analysis,
        perform_binding_site_analysis,
        perform_motif_prediction,
        handle_ai_engineering_command,
        call_llm_api
    )
except ImportError:
    from backend.ai_handlers import (
        perform_summarize,
        perform_pathway_analysis,
        perform_gene_analysis,
        perform_protein_domain_analysis,
        perform_binding_site_analysis,
        perform_motif_prediction,
        handle_ai_engineering_command,
        call_llm_api
    )

try:
    from services.reference_data import PRODORIC_PWMS
except ImportError:
    from backend.services.reference_data import PRODORIC_PWMS

try:
    from sequence_tools import compute_binding_affinity
    _SEQ_TOOLS_AVAILABLE = True
except ImportError:
    try:
        from backend.sequence_tools import compute_binding_affinity
        _SEQ_TOOLS_AVAILABLE = True
    except ImportError:
        _SEQ_TOOLS_AVAILABLE = False

router = APIRouter(tags=["AI Copilot"])
logger = logging.getLogger("app.routers.ai")


class AIEngineeringCommandRequest(BaseModel):
    command: str
    gene: str = ""
    provider: str = "google"
    api_key: str = ""
    model_name: str = ""
    base_url: str = ""


@router.post("/api/ai/engineering_command")
def ai_engineering_command_endpoint(req: AIEngineeringCommandRequest):
    return handle_ai_engineering_command(
        command=req.command,
        gene=req.gene,
        provider=req.provider,
        api_key=req.api_key,
        model_name=req.model_name,
        base_url=req.base_url
    )


@router.get("/api/summarize")
def summarize(
    gene: str = "",
    name: str = "",
    x_ai_api_key: str = Header(None, alias="X-AI-API-Key"),
    x_gemini_api_key: str = Header(None, alias="X-Gemini-API-Key"),
    x_ai_provider: str = Header("google", alias="X-AI-Provider"),
    x_ai_model: str = Header("", alias="X-AI-Model"),
    x_ai_base_url: str = Header("", alias="X-AI-Base-URL"),
):
    if not gene or not gene.strip():
        raise HTTPException(status_code=400, detail="Missing required query parameter: gene")
    api_key = x_ai_api_key or x_gemini_api_key or ""
    try:
        result = perform_summarize(gene, name, api_key, x_ai_provider, x_ai_model, x_ai_base_url)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/pathway")
def pathway(
    pathway: str = "",
    x_ai_api_key: str = Header(None, alias="X-AI-API-Key"),
    x_gemini_api_key: str = Header(None, alias="X-Gemini-API-Key"),
    x_ai_provider: str = Header("google", alias="X-AI-Provider"),
    x_ai_model: str = Header("", alias="X-AI-Model"),
    x_ai_base_url: str = Header("", alias="X-AI-Base-URL"),
):
    api_key = x_ai_api_key or x_gemini_api_key or ""
    try:
        result = perform_pathway_analysis(pathway, api_key, x_ai_provider, x_ai_model, x_ai_base_url)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/gene_assistant")
def gene_assistant(
    query: str = "",
    x_ai_api_key: str = Header(None, alias="X-AI-API-Key"),
    x_gemini_api_key: str = Header(None, alias="X-Gemini-API-Key"),
    x_ai_provider: str = Header("google", alias="X-AI-Provider"),
    x_ai_model: str = Header("", alias="X-AI-Model"),
    x_ai_base_url: str = Header("", alias="X-AI-Base-URL"),
):
    api_key = x_ai_api_key or x_gemini_api_key or ""
    try:
        result = perform_gene_analysis(query, api_key, x_ai_provider, x_ai_model, x_ai_base_url)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/protein_domain")
def protein_domain(
    gene: str = "",
    x_ai_api_key: str = Header(None, alias="X-AI-API-Key"),
    x_gemini_api_key: str = Header(None, alias="X-Gemini-API-Key"),
    x_ai_provider: str = Header("google", alias="X-AI-Provider"),
    x_ai_model: str = Header("", alias="X-AI-Model"),
    x_ai_base_url: str = Header("", alias="X-AI-Base-URL"),
):
    api_key = x_ai_api_key or x_gemini_api_key or ""
    try:
        result = perform_protein_domain_analysis(gene, api_key, x_ai_provider, x_ai_model, x_ai_base_url)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/binding_site")
def binding_site(
    gene: str = "",
    x_ai_api_key: str = Header(None, alias="X-AI-API-Key"),
    x_gemini_api_key: str = Header(None, alias="X-Gemini-API-Key"),
    x_ai_provider: str = Header("google", alias="X-AI-Provider"),
    x_ai_model: str = Header("", alias="X-AI-Model"),
    x_ai_base_url: str = Header("", alias="X-AI-Base-URL"),
):
    api_key = x_ai_api_key or x_gemini_api_key or ""
    try:
        result = perform_binding_site_analysis(gene, api_key, x_ai_provider, x_ai_model, x_ai_base_url)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/predict_motif")
def predict_motif(tf: str = ""):
    try:
        tf_lower = tf.strip().lower()

        if tf_lower in PRODORIC_PWMS:
            pwm_data = PRODORIC_PWMS[tf_lower]
            return {
                "tf": tf_lower,
                "tf_name": pwm_data.get("tf_name", tf),
                "consensus": pwm_data.get("consensus", ""),
                "pwm": pwm_data.get("pwm"),
                "nsites": pwm_data.get("targets_count", 0),
                "source": "PRODORIC (Local DB)"
            }
        else:
            for k, v in PRODORIC_PWMS.items():
                if v.get("tf_name", "").lower() == tf_lower:
                    return {
                        "tf": k,
                        "tf_name": v.get("tf_name", tf),
                        "consensus": v.get("consensus", ""),
                        "pwm": v.get("pwm"),
                        "nsites": v.get("targets_count", 0),
                        "source": "PRODORIC (Local DB)"
                    }

        result = perform_motif_prediction(tf)
        return result
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/predict_binding_affinity")
def predict_binding_affinity(tf: str = "", sequence: str = "", temperature: float = 30.0):
    if not tf or not sequence:
        raise HTTPException(status_code=400, detail="Missing tf or sequence parameter")
    try:
        tf_lower = tf.strip().lower()
        pwm = None
        tf_name = tf
        consensus = ""
        targets_count = 0

        if tf_lower in PRODORIC_PWMS:
            pwm_data = PRODORIC_PWMS[tf_lower]
            pwm = pwm_data.get("pwm")
            tf_name = pwm_data.get("tf_name", tf)
            consensus = pwm_data.get("consensus", "")
            targets_count = pwm_data.get("targets_count", 0)
        else:
            for k, v in PRODORIC_PWMS.items():
                if v.get("tf_name", "").lower() == tf_lower:
                    pwm = v.get("pwm")
                    tf_name = v.get("tf_name", tf)
                    consensus = v.get("consensus", "")
                    targets_count = v.get("targets_count", 0)
                    break

        if pwm and _SEQ_TOOLS_AVAILABLE:
            res = compute_binding_affinity(sequence, pwm, temperature)
            res["tf"] = tf
            res["tf_name"] = tf_name
            res["consensus"] = consensus
            res["targets_count"] = targets_count
            res["source"] = "PRODORIC (Local PWM)"
            return res
        else:
            raise HTTPException(status_code=404, detail=f"PWM for TF '{tf}' not found")
    except HTTPException:
        raise
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/test_ai")
def test_ai(
    x_ai_api_key: str = Header(None, alias="X-AI-API-Key"),
    x_gemini_api_key: str = Header(None, alias="X-Gemini-API-Key"),
    x_ai_provider: str = Header("google", alias="X-AI-Provider"),
    x_ai_model: str = Header("", alias="X-AI-Model"),
    x_ai_base_url: str = Header("", alias="X-AI-Base-URL"),
):
    api_key = x_ai_api_key or x_gemini_api_key or ""
    try:
        prompt = "Hello! Please return a single word: Success."
        response = call_llm_api(prompt, x_ai_provider, api_key, x_ai_model, x_ai_base_url)
        return {"status": "success", "message": f"连接成功！AI 响应：{response}"}
    except Exception as e:
        return {"status": "error", "message": str(e)}


class FetchPubMedRequest(BaseModel):
    query: str
    max_results: int = 5
    auto_ingest: bool = True


@router.post("/api/ai/literature/fetch_pubmed")
def fetch_pubmed_endpoint(req: FetchPubMedRequest):
    try:
        from rag_service import RAGService
        rag_svc = RAGService()
        articles = rag_svc.fetch_pubmed_abstracts(req.query, max_results=req.max_results)
        ingested = 0
        if req.auto_ingest and articles:
            ingested = rag_svc.ingest_fetched_articles(articles, tag=req.query)
        return {
            "query": req.query,
            "count": len(articles),
            "ingested_count": ingested,
            "articles": articles
        }
    except Exception as e:
        logger.error(f"PubMed fetch failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/ai/graph_rag/reasoning")
def graph_rag_reasoning_endpoint(
    source: str = "",
    target: str = "",
    query: str = "",
    max_depth: int = 3
):
    try:
        from services.graph_rag_service import get_graph_rag_service
        graph_rag_svc = get_graph_rag_service()
        result = graph_rag_svc.query_graph_rag_reasoning(
            source=source,
            target=target,
            query=query,
            max_depth=max_depth
        )
        return result
    except Exception as e:
        logger.error(f"GraphRAG reasoning failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))
