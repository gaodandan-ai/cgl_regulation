from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INDEX = (ROOT / "web" / "index.html").read_text(encoding="utf-8")
APP = (ROOT / "web" / "app.js").read_text(encoding="utf-8")


def test_api_client_loads_before_application():
    api_position = INDEX.index('src="lib/apiClient.js')
    loader_position = INDEX.index('src="lib/dataLoader.js')
    export_position = INDEX.index('src="lib/exportUtils.js')
    scoring_position = INDEX.index('src="lib/evidenceScoring.js')
    normalizer_position = INDEX.index('src="lib/networkNormalizer.js')
    navigation_position = INDEX.index('src="lib/queryNavigation.js')
    identifier_position = INDEX.index('src="lib/geneIdentifierIndex.js')
    render_session_position = INDEX.index('src="lib/networkRenderSession.js')
    interaction_position = INDEX.index('src="lib/networkInteractionBinder.js')
    styles_position = INDEX.index('src="lib/networkStyles.js')
    graph_position = INDEX.index('src="lib/networkGraph.js')
    ppi_loader_position = INDEX.index('src="lib/networkPpiLoader.js')
    app_position = INDEX.index('src="app.js')
    assert api_position < loader_position < export_position < scoring_position < normalizer_position < navigation_position < identifier_position < render_session_position < interaction_position < styles_position < graph_position < ppi_loader_position < app_position


def test_migrated_calls_use_shared_api_client():
    expected_calls = [
        "CglApiClient.postJson('/api/ai/engineering_command'",
        "CglApiClient.getJson(TCS_SYSTEMS_URL)",
        "CglApiClient.getJson(SIGMA_ANNOTATIONS_URL)",
        "CglApiClient.getJson('/api/network/centrality?limit=50&tfs_only=true')",
        "CglApiClient.getJson(`/api/thermo/gene_context",
    ]
    for call in expected_calls:
        assert call in APP


def test_ai_response_is_escaped_before_html_rendering():
    assert "${escapeHtml(data.result || '')}" in APP
    assert "${data.result}" not in APP


def test_startup_assets_use_shared_parallel_loader():
    assert APP.count("CglDataLoader.loadAssets") >= 3
    assert "Optional data asset unavailable" in APP
    assert "const fetches =" not in APP
    assert "fetch(EDGE_CONFIDENCE_SCORES_URL)" not in APP
    assert "fetch(IMODULON_WEIGHTS_URL)" not in APP


def test_network_and_quality_exports_use_shared_escaping():
    assert "CglExportUtils.toCsv(rows" in APP
    assert "CglExportUtils.download" in APP


def test_evidence_scoring_rules_are_imported_not_redeclared():
    assert "} = window.CglEvidenceScoring;" in APP
    assert "function combineConfidenceScores(" not in APP
    assert "function confidenceFromEvidence(" not in APP


def test_network_normalization_is_injected_and_not_redeclared():
    assert "CglNetworkNormalizer.createNormalizer" in APP
    assert "networkNormalizer.normalizeNetwork(regulations, rnaRegulations)" in APP
    assert "function normalizeTfEdge(" not in APP
    assert "function normalizeSrnaEdge(" not in APP


def test_query_navigation_uses_shared_state_machine():
    assert "CglQueryNavigation.createHistory()" in APP
    assert "queryNavigationHistory.suspend()" in APP
    assert "await renderNetwork(target)" in APP
    assert "let queryHistory =" not in APP
    assert "function normalizeQueryList(" not in APP


def test_gene_identifier_search_uses_shared_auditable_index():
    assert "CglGeneIdentifierIndex.createIndex" in APP
    assert "geneIdentifierIndex?.resolve(q)" in APP
    assert "geneIdentifierIndex.searchSuggestions(query, 15)" in APP
    assert "escapeHtml(item.display)" in APP
    assert "// 1. Process gene mappings first" not in APP


def test_network_rendering_rejects_stale_async_results():
    assert "CglNetworkRenderSession.createSession()" in APP
    assert "networkRenderSession.begin(locusTag)" in APP
    assert "networkRenderSession.isActive(renderTransaction.id)" in APP
    assert "renderTransaction.signal" in APP
    assert "const rendered = await renderNetwork(resolvedLoci)" in APP


def test_cytoscape_interactions_use_instance_scoped_binder():
    assert "window.CglNetworkInteractionBinder" in APP
    assert "networkInteractionBinder.bindLevelOfDetail(renderedCy)" in APP
    assert "networkInteractionBinder.markSharedTargets(renderedCy)" in APP
    assert "networkInteractionBinder.bindInteractions(renderedCy" in APP
    assert "let lastTapNode = null" not in APP


def test_rnaseq_network_styles_are_imported_with_threshold_dependency():
    assert "window.CglNetworkStyles" in APP
    assert "...networkStyles.createBaseNodeStyles()" in APP
    assert "...networkStyles.createBaseEdgeStyles()" in APP
    assert "...networkStyles.createRegulationEdgeStyles()" in APP
    assert "...networkStyles.createInteractionStateStyles()" in APP
    assert "...networkStyles.createRnaSeqStyles({" in APP
    assert "thresholdValue: (id, fallback) => document.getElementById(id)?.value ?? fallback" in APP
    assert "selector: 'node.rnaseq-node'" not in APP


def test_network_instance_and_ppi_edges_use_shared_graph_module():
    assert "window.CglNetworkGraph" in APP
    assert "networkGraph.createGraph({" in APP
    assert "cytoscapeImpl: cytoscape" in APP
    assert "networkGraph.addPpiEdges(renderedCy, ppiEdges)" in APP
    assert "hideEdgesOnViewport: elements.length" not in APP


def test_ppi_requests_use_parallel_abortable_loader():
    assert "window.CglNetworkPpiLoader" in APP
    assert "networkPpiLoader.loadQueryInteractions({" in APP
    assert "networkPpiLoader.loadVisibleEdges({" in APP
    assert "signal: renderTransaction.signal" in APP
    assert "for (const locus of queryList)" not in APP


def test_header_exposes_stable_windows_app_download():
    download_url = (
        "https://github.com/gaodandan-ai/cgl_regulation/"
        "releases/latest/download/cgl_regulation.exe"
    )
    assert download_url in INDEX
    assert "Download App" in INDEX
    assert INDEX.index("header-github-text") < INDEX.index("header-download-text")
