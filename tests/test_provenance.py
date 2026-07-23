from pathlib import Path

from fastapi.testclient import TestClient

from backend.app import app
from backend.services.provenance import build_provenance


ROOT = Path(__file__).resolve().parents[1]


def test_local_provenance_is_auditable():
    with TestClient(app) as client:
        response = client.get("/api/provenance")
    assert response.status_code == 200
    payload = response.json()
    assert payload["app"] == "cgl-regulation"
    assert payload["database"]["schema_version"] >= 1
    if payload["database"]["dataset_count"]:
        assert all(item["sha256"] for item in payload["database"]["datasets"])
    else:
        assert payload["status"] == "partial"
        assert payload["release_manifest"]["row_counts"]
        assert payload["warnings"]
    assert payload["evidence_policy"]["predicted"]
    assert payload["method_versions"]["evidence_scoring"] == "heuristic-confidence-v1.0.0"
    assert payload["method_versions"]["network_normalization"] == "normalized-network-v1.0.0"
    assert payload["method_versions"]["query_navigation"] == "query-navigation-v1.0.0"
    assert payload["method_versions"]["gene_identifier_index"] == "gene-identifier-index-v1.0.0"
    assert payload["method_versions"]["network_render_session"] == "network-render-session-v1.0.0"
    assert payload["method_versions"]["network_interaction_binder"] == "network-interaction-binder-v1.0.0"
    assert payload["method_versions"]["network_styles"] == "network-styles-v2.0.0"
    assert payload["method_versions"]["network_graph"] == "network-graph-v1.0.0"
    assert payload["method_versions"]["network_ppi_loader"] == "network-ppi-loader-v1.0.0"
    assert payload["interpretation_limits"]


class _UnavailableManager:
    def get_connection(self):
        return None


def test_provenance_reports_degraded_metadata_without_leaking_paths():
    payload = build_provenance(
        manager=_UnavailableManager(), root=ROOT, version="test", deployment="public"
    )
    assert payload["status"] == "partial"
    assert payload["database"]["status"] == "unavailable"
    assert "F:\\" not in str(payload)
