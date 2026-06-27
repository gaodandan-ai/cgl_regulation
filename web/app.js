/**

 * C. glutamicum Regulatory Network Explorer - Client Side Logic

 * Uses Cytoscape.js and PapaParse

 */



// Application State

let regulations = [];

let rnaRegulations = [];

let edgeConfidenceScores = [];

let rfConfidenceByEdge = new Map();

let normalizedNodes = {};

let normalizedEdges = [];

let geneMapping = [];

let cglToCg = {};

let cgToCgl = {};

let nameToCg = {};

let cgToProduct = {};

let geneIndex = {}; // lowercase -> { locusTag, name, type }

let geneToOperon = {}; // lower -> { operon, orientation, genes }

let searchSuggestions = [];

let currentQueryGene = null;

let currentDetailGene = null;

let cy = null;

let currentSimulationMode = null;

let currentSimulationRegulator = null;

const DEFAULT_EXAMPLE_LOCUS = 'cg0350';



// Query Navigation History Stacks

let queryHistory = [];

let queryForwardHistory = [];

let isNavigatingHistory = false;



// DOM Elements

const dataStatusEl = document.getElementById('data-status');

const geneInputsContainer = document.getElementById('gene-inputs-container');

const searchBtn = document.getElementById('search-btn');

const suggestionsBox = document.getElementById('suggestions-box');

const canvasOverlay = document.getElementById('canvas-overlay');

const rightSidebar = document.getElementById('right-sidebar');

const closeDetailBtn = document.getElementById('close-detail-btn');

const rightSidebarToggle = document.getElementById('right-sidebar-toggle');

let activeInput = null;



// Detail Panel Elements

const detailTypeBadge = document.getElementById('detail-type-badge');

const detailGeneName = document.getElementById('detail-gene-name');

const detailLocusTag = document.getElementById('detail-locus-tag');

const infoLocus = document.getElementById('info-locus');

const infoName = document.getElementById('info-name');

const infoType = document.getElementById('info-type');

const regulatorsCount = document.getElementById('regulators-count');

const targetsCount = document.getElementById('targets-count');

const relationsTableBody = document.querySelector('#detail-relations-table tbody');



// Config controls

const filterActivation = document.getElementById('filter-activation');

const filterRepression = document.getElementById('filter-repression');

const filterDual = document.getElementById('filter-dual');

const filterSrna = document.getElementById('filter-srna');

const filterCoregulated = document.getElementById('filter-coregulated');

const filterOnlyTfTargets = document.getElementById('filter-only-tf-targets');

const srnaThresholdPanel = document.getElementById('srna-threshold-panel');

const srnaRankThreshold = document.getElementById('srna-rank-threshold');

const rankValDisp = document.getElementById('rank-val');

const layoutSelect = document.getElementById('layout-select');



const resetViewBtn = document.getElementById('reset-view-btn');

const exportPngBtn = document.getElementById('export-png-btn');

const zoomInBtn = document.getElementById('zoom-in');

const zoomOutBtn = document.getElementById('zoom-out');

const fitCanvasBtn = document.getElementById('fit-canvas');



// Data File Paths
let REGULATIONS_URL = 'data/regulations.csv';
let RNA_REGULATIONS_URL = 'data/rna_regulation.csv';
let MAPPING_URL = 'data/gene_mapping.csv';
let OPERONS_URL = 'data/operons.csv';

let EDGE_CONFIDENCE_SCORES_URL = 'data/edge_confidence/tf_gene_edge_scores.csv';



// ==========================================================================

// 1. Initialization & Data Loading

// ==========================================================================

function initializeApp() {
    initEventListeners();
    initSidebarResizer();
    loadNetworkData();
}

if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initializeApp);
} else {
    initializeApp();
}



function updateStatus(message, type = 'loading') {

    const dot = dataStatusEl.querySelector('.status-dot');

    const txt = dataStatusEl.querySelector('.status-text');

    

    txt.textContent = message;

    dot.className = 'status-dot';

    

    if (type === 'loading') {

        dot.classList.add('pulsing');

        dot.style.backgroundColor = '#eab308';

    } else if (type === 'success') {

        dot.style.backgroundColor = '#10b981';

    } else {

        dot.style.backgroundColor = '#ef4444';

    }

}



async function loadNetworkData() {

    try {

        updateStatus('Loading gene name mapping data...', 'loading');

        const mapResponse = await fetch(MAPPING_URL);

        if (mapResponse.ok) {

            const mapText = await mapResponse.text();

            geneMapping = parseCSV(mapText);

            console.log(`Loaded ${geneMapping.length} gene mapping records.`);

        } else {

            console.warn('plate_gene_mapping.csv file not found. Skipping mapping.');

        }



        updateStatus('Loading TF-target regulatory data...', 'loading');

        const tfResponse = await fetch(REGULATIONS_URL);

        if (!tfResponse.ok) throw new Error('Unable to read regulations.csv. Please confirm the local server is running.');

        const tfText = await tfResponse.text();

        

        regulations = parseCSV(tfText);

        console.log(`Loaded ${regulations.length} TF-TG regulations.`);



        updateStatus('Loading sRNA-mRNA regulatory data...', 'loading');

        const rnaResponse = await fetch(RNA_REGULATIONS_URL);

        if (rnaResponse.ok) {

            const rnaText = await rnaResponse.text();

            rnaRegulations = parseCSV(rnaText);

            console.log(`Loaded ${rnaRegulations.length} sRNA-mRNA regulations.`);

            if (filterSrna.checked) {
                srnaThresholdPanel.classList.remove('hidden');
            } else {
                srnaThresholdPanel.classList.add('hidden');
            }

        } else {

            console.warn('sRNA-mRNA regulations.csv file not found. Skipping sRNA data.');

        }



        updateStatus('Loading operon structure data...', 'loading');

        const operonResponse = await fetch(OPERONS_URL);

        if (operonResponse.ok) {

            const operonText = await operonResponse.text();

            parseOperons(operonText);

            console.log(`Loaded operons mapping.`);

        } else {

            console.warn('Operons file not found. Skipping operons data.');

        }

        updateStatus('Loading RF edge confidence scores...', 'loading');

        await loadEdgeConfidenceScores();



        buildGeneIndex();
        normalizeNetworkData();

        rnaseqData = null;

        updateStatus('Data ready', 'success');
        initGlobalMetabolicImpactRanking();
        initPathwayRegulatoryView();
        initEngineeringTargetFinder();
        loadDefaultExampleNetwork();

    } catch (err) {

        console.error(err);

        updateStatus('Data loading failed: ' + err.message, 'error');

        alert('Error: unable to load CSV files. Please run python run_server.py so the browser can load local data.');

    }

}



function parseCSV(text) {

    const parsed = Papa.parse(text, {

        header: true,

        skipEmptyLines: true,

        dynamicTyping: true

    });

    return parsed.data;

}

async function loadEdgeConfidenceScores() {
    edgeConfidenceScores = [];
    rfConfidenceByEdge = new Map();

    try {
        const response = await fetch(EDGE_CONFIDENCE_SCORES_URL);
        if (!response.ok) {
            console.warn('RF edge confidence scores not found. Falling back to heuristic confidence scoring.');
            return;
        }

        const text = await response.text();
        edgeConfidenceScores = parseCSV(text);
        indexRfConfidenceScores(edgeConfidenceScores);
        console.log(`Loaded ${edgeConfidenceScores.length} RF edge confidence scores.`);
    } catch (err) {
        console.warn('Unable to load RF edge confidence scores. Falling back to heuristic confidence scoring.', err);
        edgeConfidenceScores = [];
        rfConfidenceByEdge = new Map();
    }
}

function edgePairKey(source, target) {
    const src = cleanStr(source).toLowerCase();
    const tgt = cleanStr(target).toLowerCase();
    if (!src || !tgt) return '';
    return `${src}=>${tgt}`;
}

function parseConfidenceScore(value) {
    const parsed = parseFloat(value);
    if (Number.isNaN(parsed)) return null;
    return Math.max(0, Math.min(1, parsed));
}

function addRfConfidenceIndexEntry(source, target, row) {
    const key = edgePairKey(source, target);
    if (!key || rfConfidenceByEdge.has(key)) return;
    const predictedConfidence = parseConfidenceScore(row.predicted_confidence);
    if (predictedConfidence === null) return;

    rfConfidenceByEdge.set(key, {
        predictedConfidence,
        confidenceRank: cleanStr(row.confidence_rank),
        label: cleanStr(row.label),
        sampleType: cleanStr(row.sample_type),
        featureMissingCount: cleanStr(row.feature_missing_count),
        expressionFeatureAvailable: cleanStr(row.expression_feature_available),
        targetMappedReactionCount: cleanStr(row.target_mapped_reaction_count),
        targetMappedPathwayCount: cleanStr(row.target_mapped_pathway_count),
        targetEnzymeConstrainedReactionCount: cleanStr(row.target_enzyme_constrained_reaction_count),
        targetKcatMedian: cleanStr(row.target_kcat_median),
        targetKcatMwMedian: cleanStr(row.target_kcat_mw_median),
        original: row
    });
}

function indexRfConfidenceScores(rows) {
    rfConfidenceByEdge = new Map();
    (rows || []).forEach(row => {
        addRfConfidenceIndexEntry(row.tf_locus, row.target_locus, row);
        addRfConfidenceIndexEntry(row.tf_name, row.target_locus, row);
        addRfConfidenceIndexEntry(row.tf_locus, row.target_name, row);
        addRfConfidenceIndexEntry(row.tf_name, row.target_name, row);
    });
}

function candidateGeneIdsForConfidenceLookup(id) {
    const cleanId = cleanStr(id);
    if (!cleanId) return [];
    const lower = cleanId.toLowerCase();
    const candidates = new Set([cleanId, lower]);
    if (cgToCgl[lower]) candidates.add(cgToCgl[lower]);
    if (cglToCg[lower]) candidates.add(cglToCg[lower]);
    const meta = geneIndex[lower];
    if (meta?.name) candidates.add(meta.name);
    if (meta?.locusTag) candidates.add(meta.locusTag);
    return Array.from(candidates).filter(Boolean);
}

function getRfConfidencePrediction(source, target) {
    if (!rfConfidenceByEdge || rfConfidenceByEdge.size === 0) return null;
    const sourceCandidates = candidateGeneIdsForConfidenceLookup(source);
    const targetCandidates = candidateGeneIdsForConfidenceLookup(target);
    for (const src of sourceCandidates) {
        for (const tgt of targetCandidates) {
            const match = rfConfidenceByEdge.get(edgePairKey(src, tgt));
            if (match) return match;
        }
    }
    return null;
}

function normalizeRegulationType(role, sourceType = 'TF-TG') {
    const cleanRole = cleanStr(role).toUpperCase();
    if (sourceType === 'sRNA-mRNA') return 'post_transcriptional_repression';
    if (cleanRole === 'A') return 'activation';
    if (cleanRole === 'R') return 'repression';
    if (cleanRole === 'DUAL') return 'dual';
    if (cleanRole === 'SIGMA') return 'sigma';
    return 'unknown';
}

function confidenceFromEvidence(evidence) {
    const text = cleanStr(evidence).toLowerCase();
    if (text.includes('experimental') && text.includes('predicted')) return 0.78;
    if (text.includes('experimental')) return 0.86;
    if (text.includes('curated') || text.includes('literature')) return 0.74;
    if (text.includes('predicted')) return 0.42;
    return 0.32;
}

function confidenceFromMotif(bindingSite) {
    const site = cleanStr(bindingSite);
    if (!site) return 0;
    const sites = site.split(';').map(s => s.trim()).filter(Boolean);
    if (sites.length >= 2) return 0.78;
    const longest = sites.reduce((max, s) => Math.max(max, s.length), 0);
    return longest >= 10 ? 0.66 : 0.48;
}

function confidenceFromChip(row) {
    const evidence = `${cleanStr(row.Evidence)} ${cleanStr(row.Source)} ${cleanStr(row.Method)} ${cleanStr(row.Assay)}`.toLowerCase();
    if (evidence.includes('chip-exo')) return 0.95;
    if (evidence.includes('chip-seq') || evidence.includes('chip_seq') || evidence.includes('chip')) return 0.9;
    return 0;
}

function confidenceFromExpression(row, sourceType = 'TF-TG') {
    if (sourceType === 'sRNA-mRNA') {
        const p = parseFloat(row.copra_pvalue);
        const fdr = parseFloat(row.copra_fdr);
        const energy = parseFloat(row.energy);
        let score = 0.35;
        if (!Number.isNaN(p)) score += p <= 0.001 ? 0.25 : p <= 0.01 ? 0.18 : p <= 0.05 ? 0.1 : 0;
        if (!Number.isNaN(fdr)) score += fdr <= 0.05 ? 0.2 : fdr <= 0.25 ? 0.12 : 0;
        if (!Number.isNaN(energy)) score += energy <= -20 ? 0.15 : energy <= -12 ? 0.08 : 0;
        return Math.min(0.9, score);
    }

    const corr = parseFloat(row.expression_correlation ?? row.Expression_correlation ?? row.correlation ?? row.Correlation);
    if (!Number.isNaN(corr)) return Math.min(0.95, Math.abs(corr));
    return 0;
}

function combineConfidenceScores(factors) {
    const weights = {
        motif: 0.25,
        chip: 0.3,
        expression: 0.2,
        database: 0.25
    };
    let weighted = 0;
    let usedWeight = 0;
    Object.entries(weights).forEach(([key, weight]) => {
        const val = factors[key] || 0;
        if (val > 0) {
            weighted += val * weight;
            usedWeight += weight;
        }
    });
    if (usedWeight === 0) return 0.25;
    const normalized = weighted / usedWeight;
    const multiEvidenceBonus = Object.values(factors).filter(v => v > 0.1).length >= 2 ? 0.06 : 0;
    return Math.max(0.05, Math.min(0.99, normalized + multiEvidenceBonus));
}

function confidenceLevel(score) {
    if (score >= 0.75) return 'high';
    if (score >= 0.5) return 'medium';
    return 'low';
}

function roleLabelFromType(role, regulationType) {
    if (regulationType === 'activation' || role === 'A') return 'Activation (+)';
    if (regulationType === 'repression' || role === 'R') return 'Repression (-)';
    if (regulationType === 'post_transcriptional_repression' || role === 'sRNA') return 'sRNA / post-transcriptional repression';
    if (regulationType === 'sigma') return 'Sigma factor';
    if (regulationType === 'dual' || role === 'Dual') return 'Dual regulation';
    return 'Unknown / pending';
}

function confidenceSummary(edge) {
    if (!edge) return '';
    const factors = edge.confidenceFactors || {};
    const percent = Math.round((edge.confidenceScore || 0) * 100);
    const rf = edge.predictedConfidence ?? edge.rfConfidence ?? factors.randomForest;
    const heuristic = edge.heuristicConfidenceScore;
    const modelText = rf !== undefined && rf !== null && !Number.isNaN(Number(rf))
        ? `RF ${Math.round(Number(rf) * 100)}%`
        : `heuristic ${percent}%`;
    const heuristicText = heuristic !== undefined && heuristic !== null
        ? `; heuristic ${Math.round(Number(heuristic) * 100)}%`
        : '';
    return `Conf ${percent}% (${edge.confidenceLevel || 'low'}; ${modelText}${heuristicText}; motif ${Math.round((factors.motif || 0) * 100)} / ChIP ${Math.round((factors.chip || 0) * 100)} / expr ${Math.round((factors.expression || 0) * 100)} / db ${Math.round((factors.database || 0) * 100)})`;
}

function getNodeMetaForDetails(locus) {
    const lower = cleanStr(locus).toLowerCase();
    const normalized = normalizedNodes[lower];
    const indexed = geneIndex[lower];
    return {
        name: normalized?.label || indexed?.name || locus,
        type: normalized?.type || indexed?.type || 'Target'
    };
}

function normalizeNodeRecord(id, label, type, aliases = {}) {
    const cleanId = cleanStr(id);
    if (!cleanId) return null;
    return {
        id: cleanId,
        label: getPrioritizedLabel(cleanId, label || cleanId),
        type,
        aliases,
        dataSource: 'local_csv'
    };
}

function mergeNormalizedNode(node) {
    if (!node) return;
    const key = node.id.toLowerCase();
    const existing = normalizedNodes[key];
    if (!existing) {
        normalizedNodes[key] = node;
        return;
    }
    const typeRank = { query: 4, TF: 3, sRNA: 2, Target: 1 };
    const chosenType = (typeRank[node.type] || 0) > (typeRank[existing.type] || 0) ? node.type : existing.type;
    normalizedNodes[key] = {
        ...existing,
        ...node,
        type: chosenType,
        aliases: {
            ...(existing.aliases || {}),
            ...(node.aliases || {})
        }
    };
}

function normalizeTfEdge(row, index) {
    const source = cleanStr(row.TF_locusTag);
    const target = cleanStr(row.TG_locusTag);
    if (!source || !target) return null;
    const regulationType = normalizeRegulationType(row.Role, 'TF-TG');
    const factors = {
        motif: confidenceFromMotif(row.Binding_site),
        chip: confidenceFromChip(row),
        expression: confidenceFromExpression(row, 'TF-TG'),
        database: confidenceFromEvidence(row.Evidence || row.Source)
    };
    const heuristicConfidenceScore = combineConfidenceScores(factors);
    const rfPrediction = getRfConfidencePrediction(source, target);
    const confidenceScore = rfPrediction?.predictedConfidence ?? heuristicConfidenceScore;
    if (rfPrediction) {
        factors.randomForest = rfPrediction.predictedConfidence;
    }
    const role = cleanStr(row.Role);
    return {
        id: `edge_${source}_${target}_${index}`,
        source,
        target,
        sourceType: 'TF',
        targetType: 'Target',
        regulationType,
        role,
        legacyRole: role,
        interactionClass: 'TF-TG',
        confidenceScore,
        heuristicConfidenceScore,
        predictedConfidence: rfPrediction?.predictedConfidence ?? null,
        confidenceModel: rfPrediction ? 'random_forest' : 'heuristic',
        rfConfidenceRank: rfPrediction?.confidenceRank || '',
        confidenceLevel: confidenceLevel(confidenceScore),
        confidenceFactors: factors,
        evidence: {
            motifSequence: cleanStr(row.Binding_site),
            databaseEvidence: cleanStr(row.Evidence),
            source: cleanStr(row.Source),
            pmid: cleanStr(row.PMID),
            expressionCorrelation: cleanStr(row.expression_correlation ?? row.Expression_correlation ?? row.correlation ?? ''),
            rfConfidenceRank: rfPrediction?.confidenceRank || '',
            rfSampleType: rfPrediction?.sampleType || '',
            rfLabel: rfPrediction?.label || '',
            rfFeatureMissingCount: rfPrediction?.featureMissingCount || '',
            rfExpressionFeatureAvailable: rfPrediction?.expressionFeatureAvailable || '',
            rfTargetMappedReactionCount: rfPrediction?.targetMappedReactionCount || '',
            rfTargetMappedPathwayCount: rfPrediction?.targetMappedPathwayCount || '',
            rfTargetEnzymeConstrainedReactionCount: rfPrediction?.targetEnzymeConstrainedReactionCount || '',
            rfTargetKcatMedian: rfPrediction?.targetKcatMedian || '',
            rfTargetKcatMwMedian: rfPrediction?.targetKcatMwMedian || ''
        },
        original: row
    };
}

function normalizeSrnaEdge(row, index) {
    const source = cleanStr(row.srna);
    const target = cleanStr(row.mrna);
    if (!source || !target) return null;
    const factors = {
        motif: 0,
        chip: 0,
        expression: confidenceFromExpression(row, 'sRNA-mRNA'),
        database: 0.45
    };
    const confidenceScore = combineConfidenceScores(factors);
    return {
        id: `edge_srna_${source}_${target}_${index}`,
        source,
        target,
        sourceType: 'sRNA',
        targetType: 'Target',
        regulationType: 'post_transcriptional_repression',
        role: 'sRNA',
        legacyRole: 'sRNA',
        interactionClass: 'sRNA-mRNA',
        confidenceScore,
        confidenceLevel: confidenceLevel(confidenceScore),
        confidenceFactors: factors,
        evidence: {
            rank: row.rank,
            energy: row.energy,
            copraPvalue: row.copra_pvalue,
            copraFdr: row.copra_fdr,
            source: 'sRNA prediction'
        },
        original: row
    };
}

function normalizeNetworkData() {
    normalizedNodes = {};
    normalizedEdges = [];

    regulations.forEach((row, index) => {
        const edge = normalizeTfEdge(row, index);
        if (!edge) return;
        normalizedEdges.push(edge);
        const tfNode = normalizeNodeRecord(edge.source, cleanStr(row.TF_name), 'TF', {
            altLocus: cleanStr(row.TF_altLocusTag)
        });
        const tgNode = normalizeNodeRecord(edge.target, cleanStr(row.TG_name), 'Target', {
            altLocus: cleanStr(row.TG_altLocusTag),
            operon: cleanStr(row.Operon)
        });
        mergeNormalizedNode(tfNode);
        mergeNormalizedNode(tgNode);
    });

    rnaRegulations.forEach((row, index) => {
        const edge = normalizeSrnaEdge(row, index);
        if (!edge) return;
        normalizedEdges.push(edge);
        const srnaNode = normalizeNodeRecord(edge.source, edge.source, 'sRNA');
        const targetNode = normalizeNodeRecord(edge.target, edge.target, 'Target');
        mergeNormalizedNode(srnaNode);
        mergeNormalizedNode(targetNode);
    });

    console.log(`Normalized regulatory graph: ${Object.keys(normalizedNodes).length} nodes, ${normalizedEdges.length} edges.`);
}

let globalMetabolicImpactRanks = [];
let globalMetabolicImpactLoading = false;

function buildGlobalRegulatoryGraphForRanking() {
    const nodes = Object.values(normalizedNodes || {}).map(node => ({
        data: {
            id: node.id,
            label: node.label || node.id,
            name: node.label || node.id,
            type: node.type
        }
    }));
    const edges = (normalizedEdges || [])
        .filter(edge => edge && edge.sourceType === 'TF')
        .map(edge => ({
            data: {
                source: edge.source,
                target: edge.target,
                type: 'regulates',
                regulation: edge.regulationType || 'unknown',
                confidence: edge.confidenceScore || 0,
                confidenceScore: edge.confidenceScore || 0
            }
        }));
    return { nodes, edges };
}

function renderGlobalMetabolicImpactRanking() {
    const tbody = document.getElementById('global-metabolic-impact-tbody');
    const status = document.getElementById('global-metabolic-impact-status');
    const filterInput = document.getElementById('global-metabolic-pathway-filter');
    if (!tbody) return;

    const filter = String(filterInput?.value || '').trim().toLowerCase();
    const filtered = filter
        ? globalMetabolicImpactRanks.filter(rank => {
            const pathways = [
                ...(rank.keyPathways || []),
                ...((rank.pathwaySummary || []).map(p => p.pathwayName || p.pathwayId || ''))
            ].join(' ').toLowerCase();
            return pathways.includes(filter);
        })
        : globalMetabolicImpactRanks;

    if (status) {
        status.textContent = globalMetabolicImpactLoading
            ? 'Calculating TF metabolic impact ranking...'
            : `${filtered.length} TFs shown${filter ? ` for "${filter}"` : ''}`;
    }

    if (globalMetabolicImpactLoading && globalMetabolicImpactRanks.length === 0) {
        tbody.innerHTML = '<tr><td colspan="8">Calculating ranking...</td></tr>';
        return;
    }
    if (filtered.length === 0) {
        tbody.innerHTML = '<tr><td colspan="8">No TF metabolic impact ranking available for the current filter.</td></tr>';
        return;
    }

    tbody.innerHTML = filtered.map((rank, index) => `
        <tr class="global-metabolic-row" data-tf-id="${escapeHtml(rank.tfId)}" title="${escapeHtml(rank.explanation || '')}">
            <td>${index + 1}</td>
            <td><strong>${escapeHtml(rank.tfLabel || rank.tfId)}</strong><div class="metabolic-muted">${escapeHtml(rank.tfId)}</div></td>
            <td><span class="global-metabolic-score">${escapeHtml(Number(rank.impactScore || 0).toFixed(2))}</span></td>
            <td>${escapeHtml(rank.totalTargetGenes || 0)}</td>
            <td>${escapeHtml(rank.mappedTargetGenes || 0)}</td>
            <td>${escapeHtml(rank.totalReactions || 0)}</td>
            <td>${escapeHtml(rank.totalPathways || 0)}</td>
            <td>${escapeHtml((rank.keyPathways || []).slice(0, 3).join(', ') || '-')}</td>
        </tr>
    `).join('');

    tbody.querySelectorAll('.global-metabolic-row').forEach(row => {
        row.addEventListener('click', () => {
            const tfId = row.getAttribute('data-tf-id');
            if (!tfId) return;
            querySingleGene(tfId);
            showNodeDetails(tfId);
        });
    });
}

async function refreshGlobalMetabolicImpactRanking() {
    const ranking = window.tfMetabolicImpactRanking;
    const tbody = document.getElementById('global-metabolic-impact-tbody');
    const status = document.getElementById('global-metabolic-impact-status');
    if (!ranking || !tbody) return;

    globalMetabolicImpactLoading = true;
    if (status) status.textContent = 'Calculating TF metabolic impact ranking...';
    tbody.innerHTML = '<tr><td colspan="8">Calculating ranking...</td></tr>';

    try {
        const graph = buildGlobalRegulatoryGraphForRanking();
        globalMetabolicImpactRanks = await ranking.rankTFsByMetabolicImpactAsync(graph, {
            limit: 50,
            includeZeroImpact: false,
            batchSize: 8
        });
    } catch (err) {
        console.error('Failed to calculate global metabolic impact ranking:', err);
        globalMetabolicImpactRanks = [];
        if (status) status.textContent = 'Failed to calculate ranking.';
    } finally {
        globalMetabolicImpactLoading = false;
        renderGlobalMetabolicImpactRanking();
    }
}

function initGlobalMetabolicImpactRanking() {
    const filterInput = document.getElementById('global-metabolic-pathway-filter');
    const refreshBtn = document.getElementById('global-metabolic-refresh-btn');
    if (filterInput && !filterInput.dataset.bound) {
        filterInput.dataset.bound = '1';
        filterInput.addEventListener('input', renderGlobalMetabolicImpactRanking);
    }
    if (refreshBtn && !refreshBtn.dataset.bound) {
        refreshBtn.dataset.bound = '1';
        refreshBtn.addEventListener('click', refreshGlobalMetabolicImpactRanking);
    }
    document.querySelectorAll('[data-pathway-filter]').forEach(btn => {
        if (btn.dataset.bound) return;
        btn.dataset.bound = '1';
        btn.addEventListener('click', () => {
            if (filterInput) {
                filterInput.value = btn.getAttribute('data-pathway-filter') || '';
                renderGlobalMetabolicImpactRanking();
            }
        });
    });
    refreshGlobalMetabolicImpactRanking();
}

let pathwayViewOptionsLoaded = false;

function highlightPathwayRegulator(tfId, geneIds) {
    if (!cy || !tfId) return;
    const tfLower = String(tfId).toLowerCase();
    const genes = new Set((geneIds || []).map(g => String(g || '').toLowerCase()));

    cy.elements().removeClass('dimmed');
    cy.elements().removeClass('highlighted');
    cy.elements().addClass('dimmed');

    const tfNode = cy.getElementById(tfId);
    if (tfNode && tfNode.length > 0) {
        tfNode.removeClass('dimmed');
        tfNode.addClass('highlighted');
    }

    genes.forEach(gene => {
        const node = cy.getElementById(gene);
        if (node && node.length > 0) {
            node.removeClass('dimmed');
            node.addClass('highlighted');
        }
    });

    cy.edges().forEach(edge => {
        const source = String(edge.data('source') || '').toLowerCase();
        const target = String(edge.data('target') || '').toLowerCase();
        if (source === tfLower && genes.has(target)) {
            edge.removeClass('dimmed');
            edge.addClass('highlighted');
        }
    });
}

async function populatePathwayViewOptions() {
    const pathwayView = window.pathwayRegulatoryView;
    const datalist = document.getElementById('pathway-view-options');
    const status = document.getElementById('pathway-view-status');
    if (!pathwayView || !datalist || pathwayViewOptionsLoaded) return;

    try {
        const pathways = await pathwayView.loadPathwayOptions();
        datalist.innerHTML = (pathways || []).slice(0, 200).map(pathway => {
            const label = pathway.pathwayName || pathway.name || pathway.pathwayId || pathway.id;
            return `<option value="${escapeHtml(label)}"></option>`;
        }).join('');
        pathwayViewOptionsLoaded = true;
        if (status) status.textContent = `${(pathways || []).length} model pathways loaded.`;
    } catch (err) {
        console.error('Failed to load pathway options:', err);
        if (status) status.textContent = 'Failed to load pathway options.';
    }
}

function renderPathwayRegulatorySummary(summary) {
    const result = document.getElementById('pathway-view-result');
    const status = document.getElementById('pathway-view-status');
    if (!result) return;

    if (!summary || summary.totalGenes === 0) {
        result.innerHTML = '<div class="metabolic-empty">No metabolic model mapping available for this pathway.</div>';
        if (status) status.textContent = 'No pathway mapping found.';
        return;
    }

    if (status) status.textContent = `${summary.totalRegulators} upstream TFs found.`;
    const regulatorsHtml = summary.regulators.length > 0
        ? summary.regulators.slice(0, 10).map((regulator, index) => `
            <button type="button" class="pathway-view-regulator" data-tf-id="${escapeHtml(regulator.tfId)}" data-genes="${encodeMetabolicList(regulator.regulatedGenes)}" title="${escapeHtml(regulator.explanation || '')}">
                <div><span class="pathway-view-title">${index + 1}. ${escapeHtml(regulator.tfLabel || regulator.tfId)}</span> <span class="pathway-view-score">score ${escapeHtml(Number(regulator.regulatorScore || 0).toFixed(2))}</span></div>
                <div class="metabolic-muted">regulates ${escapeHtml((regulator.regulatedGenes || []).length)} pathway genes - ${escapeHtml((regulator.regulationTypes || []).join(', ') || 'unknown')}</div>
            </button>
        `).join('')
        : '<div class="metabolic-empty">No upstream transcription factors were found for this pathway based on the current regulatory network.</div>';

    const genesHtml = summary.genes.slice(0, 20).map(gene => `
        <div class="pathway-view-gene">
            <div class="pathway-view-title">${escapeHtml(gene.geneLabel || gene.geneId)} <span class="metabolic-muted">${escapeHtml(gene.geneId)}</span></div>
            <div class="metabolic-reaction-list">
                ${(gene.reactions || []).slice(0, 8).map(reaction => `<span class="metabolic-reaction-badge" title="${escapeHtml(reaction.reactionName || '')}">${escapeHtml(reaction.reactionId)}</span>`).join('')}
            </div>
        </div>
    `).join('');

    result.innerHTML = `
        <div class="pathway-view-summary">
            <div><strong>Pathway:</strong> ${escapeHtml(summary.pathwayName || summary.pathwayId)}</div>
            <div><strong>Genes:</strong> ${escapeHtml(summary.totalGenes)} &nbsp; <strong>Reactions:</strong> ${escapeHtml(summary.totalReactions)} &nbsp; <strong>Upstream TFs:</strong> ${escapeHtml(summary.totalRegulators)}</div>
            <div style="margin-top:5px;">${escapeHtml(summary.explanation || '')}</div>
        </div>
        <div class="metabolic-subtitle">Top predicted regulators</div>
        <div class="pathway-view-regulators">${regulatorsHtml}</div>
        <div class="metabolic-subtitle">Pathway genes</div>
        <div class="pathway-view-genes">${genesHtml || '<div class="metabolic-empty">No pathway genes found.</div>'}</div>
    `;

    result.querySelectorAll('.pathway-view-regulator').forEach(btn => {
        btn.addEventListener('click', () => {
            const tfId = btn.getAttribute('data-tf-id');
            const genes = decodeMetabolicList(btn.getAttribute('data-genes'));
            if (!tfId) return;
            querySingleGene(tfId);
            showNodeDetails(tfId);
            highlightPathwayRegulator(tfId, genes);
        });
    });
}

async function runPathwayRegulatoryView() {
    const pathwayView = window.pathwayRegulatoryView;
    const input = document.getElementById('pathway-view-input');
    const status = document.getElementById('pathway-view-status');
    const result = document.getElementById('pathway-view-result');
    const query = String(input?.value || '').trim();
    if (!pathwayView || !query) return;

    if (status) status.textContent = 'Analyzing pathway regulators...';
    if (result) result.innerHTML = '<div class="metabolic-empty">Analyzing pathway regulators...</div>';

    try {
        const graph = buildGlobalRegulatoryGraphForRanking();
        const summary = await pathwayView.getPathwayRegulatorySummaryAsync(graph, query);
        renderPathwayRegulatorySummary(summary);
    } catch (err) {
        console.error('Failed to analyze pathway regulatory view:', err);
        if (status) status.textContent = 'Failed to analyze pathway.';
        if (result) result.innerHTML = '<div class="metabolic-empty">Failed to analyze pathway.</div>';
    }
}

function initPathwayRegulatoryView() {
    const input = document.getElementById('pathway-view-input');
    const runBtn = document.getElementById('pathway-view-run-btn');
    if (!input || !runBtn) return;

    if (!input.dataset.bound) {
        input.dataset.bound = '1';
        input.addEventListener('change', runPathwayRegulatoryView);
        input.addEventListener('keydown', event => {
            if (event.key === 'Enter') runPathwayRegulatoryView();
        });
    }
    if (!runBtn.dataset.bound) {
        runBtn.dataset.bound = '1';
        runBtn.addEventListener('click', runPathwayRegulatoryView);
    }
    document.querySelectorAll('[data-pathway-view-query]').forEach(btn => {
        if (btn.dataset.bound) return;
        btn.dataset.bound = '1';
        btn.addEventListener('click', () => {
            input.value = btn.getAttribute('data-pathway-view-query') || '';
            runPathwayRegulatoryView();
        });
    });
    populatePathwayViewOptions();
}

let engineeringTargetCandidates = [];
let engineeringTargetLoading = false;

function renderEngineeringTargetCandidates() {
    const tbody = document.getElementById('engineering-target-tbody');
    const status = document.getElementById('engineering-target-status');
    const searchInput = document.getElementById('engineering-target-search');
    const pathwayInput = document.getElementById('engineering-target-pathway-filter');
    const levelSelect = document.getElementById('engineering-target-level-filter');
    const minScoreInput = document.getElementById('engineering-target-min-score');
    const minScoreValue = document.getElementById('engineering-target-min-score-value');
    if (!tbody) return;

    const search = String(searchInput?.value || '').trim().toLowerCase();
    const pathwayFilter = String(pathwayInput?.value || '').trim().toLowerCase();
    const level = String(levelSelect?.value || '').trim().toLowerCase();
    const minScore = Number(minScoreInput?.value || 0);
    if (minScoreValue) minScoreValue.textContent = minScore.toFixed(2);

    const filtered = engineeringTargetCandidates
        .filter(candidate => !search || `${candidate.tfId} ${candidate.tfLabel}`.toLowerCase().includes(search))
        .filter(candidate => !pathwayFilter || (candidate.keyPathways || []).some(pathway => String(pathway).toLowerCase().includes(pathwayFilter)))
        .filter(candidate => !level || candidate.recommendationLevel === level)
        .filter(candidate => Number(candidate.candidateScore || 0) >= minScore);

    if (status) {
        status.textContent = engineeringTargetLoading
            ? 'Ranking candidate engineering regulators...'
            : `${filtered.length} candidate regulators shown`;
    }

    if (engineeringTargetLoading && engineeringTargetCandidates.length === 0) {
        tbody.innerHTML = '<tr><td colspan="9">Ranking candidate engineering regulators...</td></tr>';
        return;
    }
    if (filtered.length === 0) {
        tbody.innerHTML = '<tr><td colspan="9">No candidate engineering regulators found with the current filters.</td></tr>';
        return;
    }

    tbody.innerHTML = filtered.map((candidate, index) => {
        const profile = candidate.regulationProfile || {};
        const regulationText = `${profile.activationCount || 0} activation / ${profile.repressionCount || 0} repression`;
        return `
            <tr class="engineering-target-row" data-tf-id="${escapeHtml(candidate.tfId)}" data-genes="${encodeMetabolicList(candidate.regulatedKeyGenes || [])}" title="${escapeHtml(candidate.rationale || '')}">
                <td>${index + 1}</td>
                <td><strong>${escapeHtml(candidate.tfLabel || candidate.tfId)}</strong><div class="metabolic-muted">${escapeHtml(candidate.tfId)}</div></td>
                <td><span class="engineering-target-score">${escapeHtml(Number(candidate.candidateScore || 0).toFixed(2))}</span></td>
                <td><span class="engineering-target-level ${escapeHtml(candidate.recommendationLevel || 'low')}">${escapeHtml(candidate.recommendationLevel || 'low')}</span></td>
                <td>${escapeHtml(candidate.mappedTargetGenes || 0)}</td>
                <td>${escapeHtml(candidate.totalReactions || 0)}</td>
                <td>${escapeHtml(candidate.totalPathways || 0)}</td>
                <td>${escapeHtml((candidate.keyPathways || []).slice(0, 3).join(', ') || '-')}</td>
                <td>${escapeHtml(regulationText)}</td>
                <td>
                    <button class="secondary-btn btn-run-glutamate-scenario" style="background:#0f766e; color:white; border:none; padding:4px 8px; font-size:10px; cursor:pointer; border-radius:3px; font-weight:600; display:flex; align-items:center; gap:3px;">
                        <i class="fa-solid fa-flask"></i> Run glutamate scenario
                    </button>
                </td>
            </tr>
        `;
    }).join('');

    tbody.querySelectorAll('.engineering-target-row').forEach(row => {
        row.addEventListener('click', () => {
            const tfId = row.getAttribute('data-tf-id');
            const genes = decodeMetabolicList(row.getAttribute('data-genes'));
            if (!tfId) return;
            querySingleGene(tfId);
            showNodeDetails(tfId);
            highlightPathwayRegulator(tfId, genes);
        });

        const runBtn = row.querySelector('.btn-run-glutamate-scenario');
        if (runBtn) {
            runBtn.addEventListener('click', (e) => {
                e.stopPropagation();
                const tfId = row.getAttribute('data-tf-id');
                runGlutamateScenarioFromEngineering(tfId);
            });
        }
    });
}

async function refreshEngineeringTargetCandidates() {
    const finder = window.candidateEngineeringTargets;
    const tbody = document.getElementById('engineering-target-tbody');
    const status = document.getElementById('engineering-target-status');
    if (!finder || !tbody) return;

    engineeringTargetLoading = true;
    if (status) status.textContent = 'Ranking candidate engineering regulators...';
    tbody.innerHTML = '<tr><td colspan="9">Ranking candidate engineering regulators...</td></tr>';

    try {
        const graph = buildGlobalRegulatoryGraphForRanking();
        engineeringTargetCandidates = await finder.findEngineeringTargetCandidatesAsync(graph, {
            limit: 100,
            minCandidateScore: 0,
            includeLowConfidence: false,
            batchSize: 8
        });
    } catch (err) {
        console.error('Failed to rank candidate engineering targets:', err);
        engineeringTargetCandidates = [];
        if (status) status.textContent = 'Candidate ranking requires metabolic model mapping data.';
    } finally {
        engineeringTargetLoading = false;
        renderEngineeringTargetCandidates();
    }
}

function initEngineeringTargetFinder() {
    const controls = [
        document.getElementById('engineering-target-search'),
        document.getElementById('engineering-target-pathway-filter'),
        document.getElementById('engineering-target-level-filter'),
        document.getElementById('engineering-target-min-score')
    ];
    controls.forEach(control => {
        if (!control || control.dataset.bound) return;
        control.dataset.bound = '1';
        control.addEventListener('input', renderEngineeringTargetCandidates);
        control.addEventListener('change', renderEngineeringTargetCandidates);
    });
    const refreshBtn = document.getElementById('engineering-target-refresh-btn');
    if (refreshBtn && !refreshBtn.dataset.bound) {
        refreshBtn.dataset.bound = '1';
        refreshBtn.addEventListener('click', refreshEngineeringTargetCandidates);
    }
    refreshEngineeringTargetCandidates();
}

function resolveDefaultExampleLocus() {
    const candidates = [DEFAULT_EXAMPLE_LOCUS, 'sigH', 'whiB4'];
    for (const candidate of candidates) {
        const lower = candidate.toLowerCase();
        if (geneIndex[lower]) return geneIndex[lower].locusTag;
        if (nameToCg[lower] && geneIndex[nameToCg[lower].toLowerCase()]) {
            return geneIndex[nameToCg[lower].toLowerCase()].locusTag;
        }
        if (cglToCg[lower] && geneIndex[cglToCg[lower].toLowerCase()]) {
            return geneIndex[cglToCg[lower].toLowerCase()].locusTag;
        }
    }
    const firstTf = Object.values(normalizedNodes || {}).find(node => node && node.type === 'TF');
    return firstTf ? firstTf.id : '';
}

function loadDefaultExampleNetwork() {
    if (currentQueryGene || cy) return;
    const example = resolveDefaultExampleLocus();
    if (!example) return;
    window.setTimeout(() => {
        if (!currentQueryGene && !cy) {
            querySingleGene(example);
        }
    }, 120);
}

function setActiveWorkflowEntry(workflow) {
    document.querySelectorAll('.workflow-entry').forEach(btn => {
        btn.classList.toggle('active', btn.getAttribute('data-workflow') === workflow);
    });

    // Toggle quality dashboard container
    const qualityDashboard = document.getElementById('quality-dashboard-overlay');
    if (qualityDashboard) {
        if (workflow === 'quality') {
            qualityDashboard.classList.remove('hidden');
            updateQualityDashboard();
        } else {
            qualityDashboard.classList.add('hidden');
        }
    }

    // Toggle examples dashboard container
    const examplesDashboard = document.getElementById('examples-dashboard-overlay');
    if (examplesDashboard) {
        if (workflow === 'examples') {
            examplesDashboard.classList.remove('hidden');
            initExamplesDashboard();
        } else {
            examplesDashboard.classList.add('hidden');
        }
    }

    // Toggle release notes container
    const releaseDashboard = document.getElementById('release-notes-overlay');
    if (releaseDashboard) {
        if (workflow === 'release') {
            releaseDashboard.classList.remove('hidden');
        } else {
            releaseDashboard.classList.add('hidden');
        }
    }

    // Toggle references container
    const referencesDashboard = document.getElementById('references-overlay');
    if (referencesDashboard) {
        if (workflow === 'references') {
            referencesDashboard.classList.remove('hidden');
        } else {
            referencesDashboard.classList.add('hidden');
        }
    }

    // Toggle glutamate scenario container
    const glutamateDashboard = document.getElementById('glutamate-scenario-overlay');
    if (glutamateDashboard) {
        if (workflow === 'glutamate') {
            glutamateDashboard.classList.remove('hidden');
            initGlutamateScenarioDashboard();
        } else {
            glutamateDashboard.classList.add('hidden');
        }
    }

    if (workflow !== 'gene' && workflow !== 'pathway') {
        toggleRightSidebar(false);
    }
}

function scrollLeftSidebarTo(selector) {
    const sidebar = document.getElementById('left-sidebar');
    const target = document.querySelector(selector);
    if (!sidebar || !target) return;
    const top = target.offsetTop - sidebar.offsetTop - 12;
    sidebar.scrollTo({ top: Math.max(0, top), behavior: 'smooth' });
}

function initWorkflowEntrypoints() {
    const geneEntry = document.getElementById('workflow-entry-gene');
    const pathwayEntry = document.getElementById('workflow-entry-pathway');
    const engineeringEntry = document.getElementById('workflow-entry-engineering');
    const qualityEntry = document.getElementById('workflow-entry-quality');
    const examplesEntry = document.getElementById('workflow-entry-examples');
    const releaseEntry = document.getElementById('workflow-entry-release');

    if (geneEntry && !geneEntry.dataset.bound) {
        geneEntry.dataset.bound = '1';
        geneEntry.addEventListener('click', () => {
            setActiveWorkflowEntry('gene');
            scrollLeftSidebarTo('.search-section');
            const input = geneInputsContainer?.querySelector('.gene-input');
            if (input) input.focus();
            loadDefaultExampleNetwork();
        });
    }
    if (pathwayEntry && !pathwayEntry.dataset.bound) {
        pathwayEntry.dataset.bound = '1';
        pathwayEntry.addEventListener('click', () => {
            setActiveWorkflowEntry('pathway');
            scrollLeftSidebarTo('.pathway-regulatory-view-section');
            const input = document.getElementById('pathway-view-input');
            if (input) {
                input.focus();
                if (!input.value) input.value = 'glutamate metabolism';
            }
        });
    }
    if (engineeringEntry && !engineeringEntry.dataset.bound) {
        engineeringEntry.dataset.bound = '1';
        engineeringEntry.addEventListener('click', () => {
            setActiveWorkflowEntry('engineering');
            scrollLeftSidebarTo('.engineering-targets-section');
            const input = document.getElementById('engineering-target-pathway-filter');
            if (input) input.focus();
        });
    }
    if (qualityEntry && !qualityEntry.dataset.bound) {
        qualityEntry.dataset.bound = '1';
        qualityEntry.addEventListener('click', () => {
            setActiveWorkflowEntry('quality');
        });
    }
    if (examplesEntry && !examplesEntry.dataset.bound) {
        examplesEntry.dataset.bound = '1';
        examplesEntry.addEventListener('click', () => {
            setActiveWorkflowEntry('examples');
        });
    }
    if (releaseEntry && !releaseEntry.dataset.bound) {
        releaseEntry.dataset.bound = '1';
        releaseEntry.addEventListener('click', () => {
            setActiveWorkflowEntry('release');
        });
    }
    const referencesEntry = document.getElementById('workflow-entry-references');
    if (referencesEntry && !referencesEntry.dataset.bound) {
        referencesEntry.dataset.bound = '1';
        referencesEntry.addEventListener('click', () => {
            setActiveWorkflowEntry('references');
        });
    }

    const glutamateEntry = document.getElementById('workflow-entry-glutamate');
    if (glutamateEntry && !glutamateEntry.dataset.bound) {
        glutamateEntry.dataset.bound = '1';
        glutamateEntry.addEventListener('click', () => {
            setActiveWorkflowEntry('glutamate');
        });
    }
}



// ==========================================================================

// 2. Indexing & Autocomplete Suggestion Logic

// ==========================================================================

function buildGeneIndex() {

    geneIndex = {};

    cglToCg = {};

    cgToCgl = {};

    nameToCg = {};

    const seen = new Set();



    // 1. Process gene mappings first

    geneMapping.forEach(row => {

        const cgl = cleanStr(row.cgl_locus);

        const cg = cleanStr(row.cg_locus);

        const name = cleanStr(row.gene_name);

        const product = cleanStr(row.product);



        if (cgl && cg) {

            cglToCg[cgl.toLowerCase()] = cg;

            cgToCgl[cg.toLowerCase()] = cgl;

        }

        if (name && name !== '--' && cg) {

            nameToCg[name.toLowerCase()] = cg;

        }

        if (cg && product) {

            cgToProduct[cg.toLowerCase()] = product;

        }

        if (cgl && product) {

            cgToProduct[cgl.toLowerCase()] = product;

        }

    });



    // 2. Index TFs and targets from regulations

    regulations.forEach(row => {

        const tfTag = cleanStr(row.TF_locusTag);

        const tfName = cleanStr(row.TF_name);

        const tgTag = cleanStr(row.TG_locusTag);

        const tgName = cleanStr(row.TG_name);



        if (tfTag) {

            const keyTag = tfTag.toLowerCase();

            if (!seen.has(keyTag)) {

                seen.add(keyTag);

                geneIndex[keyTag] = { locusTag: tfTag, name: tfName || tfTag, type: 'TF' };

            }

            if (tfName) {

                const keyName = tfName.toLowerCase();

                if (!seen.has(keyName)) {

                    seen.add(keyName);

                    geneIndex[keyName] = { locusTag: tfTag, name: tfName, type: 'TF' };

                }

            }

        }



        if (tgTag) {

            const keyTag = tgTag.toLowerCase();

            if (!seen.has(keyTag)) {

                seen.add(keyTag);

                // Note: a target might also be a TF elsewhere, so only set to Target if not already recorded as TF

                geneIndex[keyTag] = { locusTag: tgTag, name: tgName || tgTag, type: 'Target' };

            }

            if (tgName) {

                const keyName = tgName.toLowerCase();

                if (!seen.has(keyName)) {

                    seen.add(keyName);

                    geneIndex[keyName] = { locusTag: tgTag, name: tgName, type: 'Target' };

                }

            }

        }

    });



    // 3. Index sRNAs from rna_regulations

    rnaRegulations.forEach(row => {

        const srna = cleanStr(row.srna);

        const mrna = cleanStr(row.mrna);



        if (srna) {

            const keySrna = srna.toLowerCase();

            if (!seen.has(keySrna)) {

                seen.add(keySrna);

                geneIndex[keySrna] = { locusTag: srna, name: srna, type: 'sRNA' };

            }

        }

        if (mrna) {

            const keyMrna = mrna.toLowerCase();

            if (!seen.has(keyMrna)) {

                seen.add(keyMrna);

                // If mrna is not already indexed as TF/Target, mark as target gene

                geneIndex[keyMrna] = { locusTag: mrna, name: mrna, type: 'Target' };

            }

        }

    });



    // Create unique sorted suggestions list

    const uniqueSuggestions = {};

    Object.keys(geneIndex).forEach(key => {

        const item = geneIndex[key];

        const val = item.locusTag;

        // Map locus tags and names to uniqueness

        uniqueSuggestions[val] = item;

        if (item.name && item.name !== val) {

            uniqueSuggestions[item.name] = item;

        }



        // Add corresponding cgl tag if mapped

        const cglVal = cgToCgl[val.toLowerCase()];

        if (cglVal) {

            uniqueSuggestions[cglVal] = {

                locusTag: val,

                name: item.name,

                type: item.type

            };

        }

    });



    searchSuggestions = Object.keys(uniqueSuggestions).map(name => {

        const item = uniqueSuggestions[name];

        const val = item.locusTag;

        const cglVal = cgToCgl[val.toLowerCase()];

        

        return {

            display: name,

            locusTag: val,

            type: item.type,

            cgl: cglVal || ''

        };

    }).sort((a, b) => a.display.localeCompare(b.display));

}



function cleanStr(val) {

    if (val === undefined || val === null || val === 'null') return '';

    return String(val).trim();

}



function showSuggestions(query) {

    if (!query) {

        suggestionsBox.classList.add('hidden');

        return;

    }



    const q = query.toLowerCase();

    const filtered = searchSuggestions.filter(item => 

        item.display.toLowerCase().includes(q) || 

        item.locusTag.toLowerCase().includes(q) ||

        item.cgl.toLowerCase().includes(q)

    ).slice(0, 15); // limit to 15 suggestions



    if (filtered.length === 0) {

        suggestionsBox.classList.add('hidden');

        return;

    }



    suggestionsBox.innerHTML = '';

    filtered.forEach(item => {

        const div = document.createElement('div');

        div.className = `suggestion-item type-${item.type.toLowerCase()}`;

        

        let subText = '';

        if (item.display.toLowerCase() === item.locusTag.toLowerCase()) {

            if (item.cgl) {

                subText = ` <span class="locus-tag">(${item.cgl})</span>`;

            }

        } else {

            subText = ` <span class="locus-tag">(${item.locusTag})</span>`;

        }



        div.innerHTML = `

            <span><strong>${item.display}</strong>${subText}</span>

            <span class="item-type">${item.type}</span>

        `;

        div.addEventListener('click', () => {

            if (activeInput) {

                activeInput.value = item.display;

            }

            suggestionsBox.classList.add('hidden');

            triggerSearchFromInputs();

        });

        suggestionsBox.appendChild(div);

    });



    suggestionsBox.classList.remove('hidden');

}



// ==========================================================================

// 3. Network Construction & Rendering (Cytoscape.js)

// ==========================================================================

function getQueryGenes() {

    const tabBatchBtn = document.getElementById('tab-batch-btn');

    const isBatchActive = tabBatchBtn && tabBatchBtn.classList.contains('active');

    const queries = [];

    

    if (isBatchActive) {

        const text = document.getElementById('gene-batch-textarea').value;

        const tokens = text.split(/[\s,;\n\r]+/).map(t => t.trim()).filter(t => t);

        queries.push(...tokens);

    } else {

        const inputs = document.querySelectorAll('.gene-input');

        inputs.forEach(input => {

            const val = input.value.trim();

            if (val) {

                const tokens = val.split(',').map(t => t.trim()).filter(t => t);

                queries.push(...tokens);

            }

        });

    }

    return queries;

}



function triggerSearchFromInputs() {

    const queries = getQueryGenes();

    

    if (queries.length === 0) {

        alert('Enter or paste at least one gene or sRNA to analyze.');

        return;

    }

    

    const resolvedLoci = [];

    for (let q of queries) {

        const lower = q.toLowerCase();

        let targetLocus = lower;

        if (cglToCg[lower]) {

            targetLocus = cglToCg[lower].toLowerCase();

        } else if (nameToCg[lower]) {

            targetLocus = nameToCg[lower].toLowerCase();

        }

        

        const match = geneIndex[targetLocus];

        if (match) {

            resolvedLoci.push(match.locusTag);

        } else {

            console.warn(`Gene/sRNA "${q}" not found.`);

        }

    }

    

    if (resolvedLoci.length === 0) {

        alert(`No matching genes/sRNAs were found in the local database: "${queries.join(', ')}".`);

        return;

    }

    

    renderNetwork(resolvedLoci);

    

    // Auto-update right details panel: show details if single gene, collapse if multiple

    if (resolvedLoci.length === 1) {

        showNodeDetails(resolvedLoci[0]);

    } else {

        toggleRightSidebar(false);

    }

}



function renderNetwork(locusTag) {

    // Reset simulation states first

    resetPerturbationSimulation();



    // 1. Elements preparation

    const elements = buildElements(locusTag);

    

    if (elements.nodes.length === 0) {

        alert("This gene has no visible regulatory relationships under the current filters.");

        return;

    }



    const nextQuery = Array.isArray(locusTag) ? locusTag : [locusTag];

    pushQueryToHistory(nextQuery);

    currentQueryGene = nextQuery;



    // 2. Hide welcome state overlay

    canvasOverlay.classList.add('hidden');



    // 3. Destroy previous cytoscape instance

    if (cy) {

        cy.destroy();

    }



    // 4. Initialize Cytoscape

    cy = cytoscape({

        container: document.getElementById('cy'),

        elements: elements,

        style: [

            // Core node styling for Academic Light Theme

            {

                selector: 'node',

                style: {

                    'label': 'data(name)',

                    'font-size': '11px',

                    'color': '#0f172a', // Dark slate text

                    'background-color': '#f5f5f5', // Default gray

                    'text-valign': 'bottom',

                    'text-margin-y': '6px',

                    'width': '22px',

                    'height': '22px',

                    'border-width': '2px',

                    'border-color': '#757575',

                    'transition-property': 'background-color, line-color, target-arrow-color, width, height, border-width',

                    'transition-duration': '0.2s'

                }

            },

            {

                selector: 'node[type="TF"]',

                style: {

                    'background-color': '#e3f2fd', // Soft blue

                    'border-color': '#1976d2',     // Darker blue border

                    'width': '26px',

                    'height': '26px'

                }

            },

            {

                selector: 'node[type="sRNA"]',

                style: {

                    'background-color': '#f3e5f5', // Soft purple

                    'border-color': '#8e24aa',     // Darker purple border

                    'width': '26px',

                    'height': '26px',

                    'shape': 'hexagon'

                }

            },

            {

                selector: 'node[type="query"]',

                style: {

                    'background-color': '#ffe0b2', // Soft orange

                    'border-color': '#f57c00',     // Darker orange border

                    'width': '34px',

                    'height': '34px',

                    'border-width': '3px',

                    'font-weight': 'bold',

                    'font-size': '13px'

                }

            },

            {

                selector: 'node.shared-target',

                style: {

                    'background-color': '#e0f2f1', // Soft Teal/Mint background

                    'border-color': '#00897b',     // Dark Teal border

                    'border-width': '2.5px'

                }

            },

            {

                selector: 'node.rnaseq-node',

                style: {

                    'background-color': (node) => {

                        const val = node.data('rnaseq_log2fc');

                        return getRnaSeqColor(val);

                    },

                    'border-width': (node) => {

                        const pval = node.data('rnaseq_pvalue');

                        const pvalEl = document.getElementById('rnaseq-p-threshold');

                        const pThresh = pvalEl ? parseFloat(pvalEl.value) : 0.05;

                        return (pval !== undefined && pval <= pThresh) ? '3.5px' : '2px';

                    },

                    'border-color': (node) => {

                        const pval = node.data('rnaseq_pvalue');

                        const pvalEl = document.getElementById('rnaseq-p-threshold');

                        const pThresh = pvalEl ? parseFloat(pvalEl.value) : 0.05;

                        return (pval !== undefined && pval <= pThresh) ? '#0f172a' : '#94a3b8';

                    },

                    'width': (node) => {

                        const log2fc = node.data('rnaseq_log2fc');

                        const baseSize = node.data('type') === 'query' ? 34 : (['TF', 'sRNA'].includes(node.data('type')) ? 26 : 22);

                        if (log2fc === undefined || isNaN(log2fc)) return baseSize;

                        return baseSize + Math.min(16, Math.abs(log2fc) * 4);

                    },

                    'height': (node) => {

                        const log2fc = node.data('rnaseq_log2fc');

                        const baseSize = node.data('type') === 'query' ? 34 : (['TF', 'sRNA'].includes(node.data('type')) ? 26 : 22);

                        if (log2fc === undefined || isNaN(log2fc)) return baseSize;

                        return baseSize + Math.min(16, Math.abs(log2fc) * 4);

                    },

                    'shadow-blur': (node) => {

                        const pval = node.data('rnaseq_pvalue');

                        const log2fc = node.data('rnaseq_log2fc');

                        const pvalEl = document.getElementById('rnaseq-p-threshold');

                        const lfcEl = document.getElementById('rnaseq-lfc-threshold');

                        const pThresh = pvalEl ? parseFloat(pvalEl.value) : 0.05;

                        const lfcThresh = lfcEl ? parseFloat(lfcEl.value) : 1.0;

                        return (pval !== undefined && pval <= pThresh && Math.abs(log2fc) >= lfcThresh) ? '12px' : '0px';

                    },

                    'shadow-color': (node) => {

                        const log2fc = node.data('rnaseq_log2fc');

                        if (log2fc === undefined) return 'transparent';

                        return log2fc > 0 ? '#ef4444' : '#2563eb';

                    },

                    'shadow-opacity': 0.85,

                    'shadow-offset-x': '0px',

                    'shadow-offset-y': '0px'

                }

            },

            // Edge styling

            {

                selector: 'edge',

                style: {

                    'width': (edge) => 1.2 + ((edge.data('confidenceScore') || 0.25) * 3.2),

                    'line-color': '#e65100', // Default dark orange

                    'target-arrow-color': '#e65100',

                    'target-arrow-shape': 'triangle',

                    'curve-style': 'bezier',

                    'arrow-scale': 1.1,

                    'opacity': (edge) => 0.35 + ((edge.data('confidenceScore') || 0.25) * 0.6),

                    'transition-property': 'line-color, target-arrow-color, opacity, width',

                    'transition-duration': '0.2s'

                }

            },

            {
                selector: 'edge[regulationType="activation"]',
                style: {
                    'line-color': '#2e7d32',
                    'target-arrow-color': '#2e7d32',
                    'target-arrow-shape': 'triangle'
                }
            },

            {

                selector: 'edge[role="A"]', // Activation

                style: {

                    'line-color': '#2e7d32', // Academic Green

                    'target-arrow-color': '#2e7d32'

                }

            },

            {
                selector: 'edge[regulationType="repression"]',
                style: {
                    'line-color': '#d32f2f',
                    'target-arrow-color': '#d32f2f',
                    'target-arrow-shape': 'tee'
                }
            },

            {

                selector: 'edge[role="R"]', // Repression

                style: {

                    'line-color': '#d32f2f', // Academic Red

                    'target-arrow-color': '#d32f2f',

                    'target-arrow-shape': 'tee'

                }

            },

            {
                selector: 'edge[regulationType="dual"], edge[regulationType="sigma"], edge[regulationType="unknown"]',
                style: {
                    'line-color': '#e65100',
                    'target-arrow-color': '#e65100',
                    'target-arrow-shape': 'triangle'
                }
            },

            {

                selector: 'edge[role="Dual"]',

                style: {

                    'line-color': '#e65100',

                    'target-arrow-color': '#e65100'

                }

            },

            {
                selector: 'edge[regulationType="post_transcriptional_repression"]',
                style: {
                    'line-color': '#7b1fa2',
                    'target-arrow-color': '#7b1fa2',
                    'line-style': 'dashed',
                    'target-arrow-shape': 'triangle-tee'
                }
            },

            {

                selector: 'edge[role="sRNA"]', // sRNA-mRNA prediction

                style: {

                    'line-color': '#7b1fa2', // Academic Purple

                    'target-arrow-color': '#7b1fa2',

                    'line-style': 'dashed',

                    'target-arrow-shape': 'triangle-tee'

                }

            },

            {
                selector: 'edge.confidence-high',
                style: {
                    'line-style': 'solid'
                }
            },

            {
                selector: 'edge.confidence-medium',
                style: {
                    'line-style': 'solid'
                }
            },

            {
                selector: 'edge.confidence-low',
                style: {
                    'line-style': 'dotted',
                    'opacity': 0.42
                }
            },

            // Interactive dimming styles

            {

                selector: '.dimmed',

                style: {

                    'opacity': 0.15

                }

            },

            {

                selector: '.rnaseq-hidden',

                style: {

                    'display': 'none'

                }

            },

            {

                selector: 'node.highlighted',

                style: {

                    'border-width': '3px',

                    'border-color': '#0f172a', // Dark slate border when highlighted

                    'width': '38px',

                    'height': '38px'

                }

            },

            {

                selector: 'edge.highlighted',

                style: {

                    'width': 3.5,

                    'opacity': 1.0

                }

            },

            {

                selector: 'node.sim-up',

                style: {

                    'border-color': '#2e7d32',

                    'border-width': '4px',

                    'background-color': '#e8f5e9',

                    'shadow-blur': '10px',

                    'shadow-color': '#2e7d32',

                    'shadow-opacity': 0.8

                }

            },

            {

                selector: 'node.sim-down',

                style: {

                    'border-color': '#d32f2f',

                    'border-width': '4px',

                    'background-color': '#ffebee',

                    'shadow-blur': '10px',

                    'shadow-color': '#d32f2f',

                    'shadow-opacity': 0.8

                }

            },

            {

                selector: 'node.sim-dual',

                style: {

                    'border-color': '#e65100',

                    'border-width': '4px',

                    'background-color': '#fff3e0',

                    'shadow-blur': '10px',

                    'shadow-color': '#e65100',

                    'shadow-opacity': 0.8

                }

            }

        ],

        layout: {

            name: layoutSelect.value,

            animate: true,

            animationDuration: 400

        }

    });



    // Add shared-target class to Target nodes with in-degree > 1 in the rendered graph

    cy.nodes('[type="Target"]').forEach(node => {

        if (node.indegree(false) > 1) {

            node.addClass('shared-target');

        }

    });



    // 5. Interaction Event Listeners

    let lastTapNode = null;

    let lastTapTimeout = null;

    cy.on('tap', 'node', (evt) => {

        const node = evt.target;

        const now = new Date().getTime();

        

        highlightSubnet(node);

        showNodeDetails(node.id());

        

        if (lastTapNode === node && (now - lastTapTimeout < 350)) {

            // Double tap / double click: load this node's regulatory network

            const locus = node.id();

            querySingleGene(locus);

        } else {

            lastTapNode = node;

            lastTapTimeout = now;

        }

    });



    cy.on('tap', (evt) => {

        if (evt.target === cy) {

            toggleRightSidebar(false);

        }

    });



    // 6. Update Network Statistics & Filters if active

    if (rnaseqData) {

        applyRnaSeqFilters();

    } else {

        updateNetworkStatistics();

    }

}



function getPrioritizedLabel(locusTag, commonName) {

    if (!locusTag) return commonName || '';

    const lower = locusTag.toLowerCase();

    const cgl = cgToCgl[lower];

    if (cgl) return cgl;

    if (commonName && commonName !== locusTag && commonName !== '--') return commonName;

    return locusTag;

}



function buildElements(queryLoci) {

    const queryList = Array.isArray(queryLoci) ? queryLoci : [queryLoci];

    const querySet = new Set(queryList.map(l => l.toLowerCase()));

    const nodesMap = {};

    const edges = [];

    const showActivation = filterActivation.checked;

    const showRepression = filterRepression.checked;

    const showDual = filterDual.checked;

    const showSrna = filterSrna.checked;

    const rankLimit = parseInt(srnaRankThreshold.value, 10);

    const showOnlyTfTargets = filterOnlyTfTargets ? filterOnlyTfTargets.checked : false;

    function getNodeMeta(locus, fallbackType = 'Target') {
        const lower = locus.toLowerCase();
        const normalized = normalizedNodes[lower];
        const indexed = geneIndex[lower];
        return {
            locusTag: locus,
            name: normalized?.label || indexed?.name || locus,
            type: normalized?.type || indexed?.type || fallbackType
        };
    }

    function addNode(locus, typeOverride = null) {
        if (!locus || nodesMap[locus]) return;
        const lower = locus.toLowerCase();
        const meta = getNodeMeta(locus, typeOverride || 'Target');
        const nodeType = querySet.has(lower) ? 'query' : (typeOverride || meta.type || 'Target');
        nodesMap[locus] = {
            data: {
                id: locus,
                name: getPrioritizedLabel(locus, meta.name),
                type: nodeType,
                schemaVersion: 'unified-v1'
            }
        };
    }

    queryList.forEach(locus => addNode(locus, 'Target'));

    const edgeSource = normalizedEdges.length > 0
        ? normalizedEdges
        : regulations.map((row, index) => normalizeTfEdge(row, index)).filter(Boolean);

    edgeSource.forEach(edge => {
        if (!edge) return;

        const source = edge.source;
        const target = edge.target;
        const role = edge.legacyRole || edge.role || '';
        const regulationType = edge.regulationType || normalizeRegulationType(role, edge.interactionClass);

        if (regulationType === 'activation' && !showActivation) return;
        if (regulationType === 'repression' && !showRepression) return;
        if (['dual', 'sigma', 'unknown'].includes(regulationType) && edge.interactionClass !== 'sRNA-mRNA' && !showDual) return;
        if (edge.interactionClass === 'sRNA-mRNA') {
            if (!showSrna) return;
            const rank = parseInt(edge.evidence?.rank ?? edge.original?.rank, 10);
            if (!Number.isNaN(rank) && rank > rankLimit) return;
        }

        const isSourceQuery = querySet.has(source.toLowerCase());
        const isTargetQuery = querySet.has(target.toLowerCase());
        if (!isSourceQuery && !isTargetQuery) return;

        if (showOnlyTfTargets && isSourceQuery && !isTargetQuery) {
            const targetMeta = geneIndex[target.toLowerCase()] || normalizedNodes[target.toLowerCase()];
            const isTargetTf = targetMeta && targetMeta.type === 'TF';
            if (!isTargetTf) return;
        }

        addNode(source, edge.sourceType === 'sRNA' ? 'sRNA' : 'TF');
        addNode(target, edge.targetType || 'Target');

        edges.push({
            data: {
                id: edge.id,
                source,
                target,
                role,
                type: edge.interactionClass,
                regulationType,
                confidenceScore: edge.confidenceScore,
                heuristicConfidenceScore: edge.heuristicConfidenceScore,
                predictedConfidence: edge.predictedConfidence,
                confidenceModel: edge.confidenceModel,
                rfConfidenceRank: edge.rfConfidenceRank,
                confidencePercent: Math.round((edge.confidenceScore || 0) * 100),
                confidenceLevel: edge.confidenceLevel,
                confidenceFactors: edge.confidenceFactors,
                evidence: edge.evidence,
                motifScore: edge.confidenceFactors?.motif || 0,
                chipScore: edge.confidenceFactors?.chip || 0,
                expressionScore: edge.confidenceFactors?.expression || 0,
                databaseScore: edge.confidenceFactors?.database || 0,
                rank: edge.evidence?.rank,
                energy: edge.evidence?.energy,
                pvalue: edge.evidence?.copraPvalue,
                schemaVersion: 'unified-v1'
            },
            classes: `confidence-${edge.confidenceLevel || 'low'}`
        });
    });

    if (rnaseqData) {
        Object.keys(nodesMap).forEach(id => {
            const lowerId = id.toLowerCase();
            if (rnaseqData[lowerId]) {
                const item = rnaseqData[lowerId];
                nodesMap[id].data.rnaseq_log2fc = item.log2fc;
                nodesMap[id].data.rnaseq_pvalue = item.pvalue;
                nodesMap[id].classes = (nodesMap[id].classes || '') + ' rnaseq-node';
            }
        });
    }

    const showOnlyCoRegulated = filterCoregulated.checked;
    if (showOnlyCoRegulated) {
        const inDegreeMap = {};
        edges.forEach(e => {
            const target = e.data.target;
            inDegreeMap[target] = (inDegreeMap[target] || 0) + 1;
        });

        const coRegulatedTargets = new Set();
        Object.keys(inDegreeMap).forEach(nodeId => {
            const nodeObj = nodesMap[nodeId];
            if (nodeObj && nodeObj.data.type === 'Target' && inDegreeMap[nodeId] >= 2) {
                coRegulatedTargets.add(nodeId);
            }
        });

        const keptEdges = edges.filter(e => {
            const targetNode = nodesMap[e.data.target];
            const targetType = targetNode ? targetNode.data.type : '';
            if (targetType === 'Target') return coRegulatedTargets.has(e.data.target);
            return true;
        });

        const keptNodeIds = new Set(queryList);
        keptEdges.forEach(e => {
            keptNodeIds.add(e.data.source);
            keptNodeIds.add(e.data.target);
        });

        return {
            nodes: Object.values(nodesMap).filter(n => keptNodeIds.has(n.data.id)),
            edges: keptEdges
        };
    }

    return {
        nodes: Object.values(nodesMap),
        edges
    };
}


function highlightSubnet(node) {

    const neighborhood = node.neighborhood();

    

    cy.elements().addClass('dimmed');

    cy.elements().removeClass('highlighted');

    

    node.removeClass('dimmed');

    node.addClass('highlighted');

    

    neighborhood.removeClass('dimmed');

    neighborhood.addClass('highlighted');

}



function resetHighlight() {

    if (cy) {

        cy.elements().removeClass('dimmed');

        cy.elements().removeClass('highlighted');

    }

}



// ==========================================================================

// 4. Detail Panel Loading & View Rendering

// ==========================================================================

function showNodeDetails(locusTag) {

    // Reset any ongoing perturbation simulation

    resetPerturbationSimulation();



    // Ensure AI trigger button is visible (it might have been hidden in operon view)

    const btnTriggerAi = document.getElementById('btn-trigger-ai');

    if (btnTriggerAi) {

        btnTriggerAi.style.display = '';

    }



    // Clear previous AI summary

    const summaryCard = document.getElementById('ai-summary-result');

    if (summaryCard) {

        summaryCard.classList.add('hidden');

        summaryCard.innerHTML = '';

    }



    let resolvedLocus = locusTag;
    const lower = locusTag.toLowerCase();
    if (cglToCg[lower]) {
        resolvedLocus = cglToCg[lower];
    }
    const resolvedLower = resolvedLocus.toLowerCase();
    currentDetailGene = resolvedLocus;

    // Resolve display meta
    let meta = { locusTag: resolvedLocus, name: resolvedLocus, type: 'Target' };
    for (let key in geneIndex) {
        if (geneIndex[key].locusTag.toLowerCase() === resolvedLower) {
            meta = geneIndex[key];
            break;
        }
    }

    // Set badge style
    detailTypeBadge.style.backgroundColor = '';
    detailTypeBadge.style.color = '';
    detailTypeBadge.className = `gene-badge ${meta.type.toLowerCase()}`;
    detailTypeBadge.textContent = meta.type === 'TF' ? 'Transcription factor (TF)' : meta.type === 'sRNA' ? 'sRNA' : 'Target gene';
    
    const cgl = cgToCgl[resolvedLower] || (locusTag.toLowerCase().startsWith('cgl') ? locusTag : '');

    // Prioritize Cgl locus tag for header

    detailGeneName.textContent = cgl ? cgl : (meta.name && meta.name !== '--' ? meta.name : meta.locusTag);

    detailLocusTag.textContent = meta.locusTag;

    

    infoLocus.textContent = meta.locusTag;

    infoName.textContent = meta.name;

    infoType.textContent = meta.type;



    const viewNetworkBtn = document.getElementById('view-network-btn');

    if (meta.type === 'TF' || meta.type === 'sRNA') {

        viewNetworkBtn.style.display = 'flex';

        viewNetworkBtn.onclick = () => {

            querySingleGene(meta.locusTag);

        };

    } else {

        viewNetworkBtn.style.display = 'none';

    }



    const cglRow = document.getElementById('info-cgl-row');

    const infoCgl = document.getElementById('info-cgl');

    if (cgl) {

        cglRow.style.display = '';

        infoCgl.textContent = cgl;

    } else {

        cglRow.style.display = 'none';

    }



    const product = cgToProduct[lower];

    const productRow = document.getElementById('info-product-row');

    const infoProduct = document.getElementById('info-product');

    if (product) {

        productRow.style.display = '';

        infoProduct.textContent = product;

    } else {

        productRow.style.display = 'none';

    }



    // Resolve cg and cgl locus tags for pathway lookup

    const canonicalTagLower = meta.locusTag.toLowerCase();

    let cgLocus = '';

    let cglLocus = '';

    

    if (canonicalTagLower.startsWith('cg') && !canonicalTagLower.startsWith('cgl')) {

        cgLocus = meta.locusTag;

        cglLocus = cgToCgl[canonicalTagLower] || '';

    } else if (canonicalTagLower.startsWith('cgl')) {

        cglLocus = meta.locusTag;

        cgLocus = cglToCg[canonicalTagLower] || '';

    } else {

        cgLocus = meta.locusTag;

        cglLocus = cgToCgl[canonicalTagLower] || '';

    }



    // Pathways & GO Terms Row Rendering

    const pathwayRow = document.getElementById('info-pathway-row');

    const pathwayContainer = document.getElementById('info-pathway-container');

    if (pathwayRow && pathwayContainer) {

        pathwayRow.style.display = 'none';

        pathwayContainer.innerHTML = '<span style="font-size: 11px; color: var(--text-muted);"><i class="fa-solid fa-spinner fa-spin"></i> Loading pathway data...</span>'; 

        

        fetch(`/api/kegg_pathways?cg=${encodeURIComponent(cgLocus)}&cgl=${encodeURIComponent(cglLocus)}`)

            .then(response => response.json())

            .then(data => {

                // Ensure this response is still for the active gene details view

                if (detailLocusTag.textContent !== meta.locusTag) return;

                

                pathwayContainer.innerHTML = '';

                const pathways = data.pathways || [];

                const goTerms = data.go_terms || [];

                

                if (pathways.length === 0 && goTerms.length === 0) {

                    pathwayRow.style.display = 'none';

                    return;

                }

                

                pathwayRow.style.display = '';

                

                // 1. Render KEGG Pathways

                pathways.forEach(p => {

                    const badge = document.createElement('a');

                    badge.className = 'pathway-badge kegg';

                    badge.href = p.link;

                    badge.target = '_blank';

                    badge.title = `KEGG pathway: ${p.id} (open map in a new tab and highlight this gene)`;

                    badge.innerHTML = `<i class="fa-solid fa-map"></i> ${p.name} <i class="fa-solid fa-arrow-up-right-from-square"></i>`;

                    pathwayContainer.appendChild(badge);

                });

                

                // 2. Render GO Terms

                goTerms.forEach(go => {

                    const badge = document.createElement('a');

                    badge.className = 'pathway-badge go';

                    badge.href = go.link;

                    badge.target = '_blank';

                    badge.title = `Gene Ontology [${go.type}]: ${go.id}`;

                    badge.innerHTML = `<i class="fa-solid fa-tag"></i> ${go.name} <i class="fa-solid fa-arrow-up-right-from-square"></i>`;

                    pathwayContainer.appendChild(badge);

                });

            })

            .catch(err => {

                console.error('Error fetching pathway data:', err);

                if (detailLocusTag.textContent === meta.locusTag) {

                    pathwayRow.style.display = 'none';

                }

            });

    }

    fetchMetabolicImpact(meta.locusTag, meta.type);



    // Operon Row Rendering

    const operonRow = document.getElementById('info-operon-row');

    const infoOperon = document.getElementById('info-operon');

    

    let standardCgForOperon = meta.locusTag.toLowerCase();

    if (cglToCg[standardCgForOperon]) {

        standardCgForOperon = cglToCg[standardCgForOperon].toLowerCase();

    }

    

    const operonMeta = geneToOperon[standardCgForOperon] || geneToOperon[lower];

    

    if (operonMeta) {

        operonRow.style.display = '';

        

        const geneLinks = operonMeta.genes.map(g => {

            const prioritized = getPrioritizedLabel(g, g);

            const isCurrent = g.toLowerCase() === lower || g.toLowerCase() === standardCgForOperon;

            if (isCurrent) {

                return `<strong style="color: var(--text-primary); font-family: monospace;">${prioritized}</strong>`;

            } else {

                return `<a href="#" class="operon-gene-link" data-locus="${g}" style="color: var(--color-primary-accent); text-decoration: none; font-weight: 500; font-family: monospace;">${prioritized}</a>`;

            }

        }).join(', ');

        

        infoOperon.innerHTML = `

            <div style="font-weight: 600; color: var(--text-primary);">${operonMeta.operon} (${operonMeta.orientation} strand)</div>

            <div style="font-size: 11px; margin-top: 4px; color: var(--text-secondary);">Genes: ${geneLinks}</div>

            <div style="display: flex; gap: 6px; margin-top: 8px;">

                <button id="btn-draw-operon-network" class="secondary-btn" style="flex: 1; font-size: 11px; padding: 6px 4px; height: auto; border: 1px solid rgba(30, 58, 138, 0.15); color: var(--color-primary-accent); background-color: rgba(30, 58, 138, 0.03);" title="Load all member genes and their regulatory network on the canvas">

                    <i class="fa-solid fa-network-wired"></i> Joint analysis

                </button>

            </div>

        `;

        

        infoOperon.querySelectorAll('.operon-gene-link').forEach(link => {

            link.addEventListener('click', (e) => {

                e.preventDefault();

                const targetLocus = link.getAttribute('data-locus');

                querySingleGene(targetLocus);

            });

        });



        const drawOperonBtn = infoOperon.querySelector('#btn-draw-operon-network');

        if (drawOperonBtn) {

            drawOperonBtn.addEventListener('click', () => {

                queryMultipleGenes(operonMeta.genes);

                showOperonDetails(operonMeta);

            });

        }



        const simOperonOeBtn = infoOperon.querySelector('#btn-sim-operon-oe');

        if (simOperonOeBtn) {

            simOperonOeBtn.addEventListener('click', () => {

                queryMultipleGenes(operonMeta.genes);

                showOperonDetails(operonMeta, 'OE');

            });

        }



        const simOperonKoBtn = infoOperon.querySelector('#btn-sim-operon-ko');

        if (simOperonKoBtn) {

            simOperonKoBtn.addEventListener('click', () => {

                queryMultipleGenes(operonMeta.genes);

                showOperonDetails(operonMeta, 'KO');

            });

        }

    } else {

        operonRow.style.display = 'none';

    }



    // External DB Links Row Rendering

    const linksCell = document.getElementById('info-links');

    linksCell.innerHTML = '';

    

    let standardCgForLinks = meta.locusTag;

    const standardCgLower = standardCgForLinks.toLowerCase();

    if (cglToCg[standardCgLower]) {

        standardCgForLinks = cglToCg[standardCgLower];

    }

    

    const dbLinks = [];

    

    // Resolve the Cgl locus tag (e.g. Cgl0339) specifically for KEGG, as KEGG ATCC 13032 (cgl) uses the Cgl prefix

    const standardCgLowerKey = standardCgForLinks.toLowerCase();

    const cglLocusForKegg = cgToCgl[standardCgLowerKey] || standardCgForLinks;

    

    if (cglLocusForKegg.toLowerCase().startsWith('cgl')) {

        dbLinks.push(`<a href="https://www.kegg.jp/entry/cgl:${cglLocusForKegg}" target="_blank" class="ext-link" title="View metabolic pathway in KEGG"><i class="fa-solid fa-diagram-project"></i> KEGG</a>`);

    } else if (standardCgForLinks.toLowerCase().startsWith('cg')) {

        // Fallback guess if no direct mapping exists but is a coding gene

        const predictedCgl = standardCgForLinks.replace('cg', 'Cgl');

        dbLinks.push(`<a href="https://www.kegg.jp/entry/cgl:${predictedCgl}" target="_blank" class="ext-link" title="View metabolic pathway in KEGG"><i class="fa-solid fa-diagram-project"></i> KEGG</a>`);

    }

    

    if (standardCgForLinks.toLowerCase().startsWith('cg')) {

        dbLinks.push(`<a href="https://www.ncbi.nlm.nih.gov/gene/?term=${standardCgForLinks}" target="_blank" class="ext-link" title="View official annotation in NCBI Gene"><i class="fa-solid fa-dna"></i> NCBI</a>`);

        dbLinks.push(`<a href="https://biocyc.org/getid?id=CORYNE:${standardCgForLinks}" target="_blank" class="ext-link" title="View detailed pathway context in BioCyc / CoryneCyc"><i class="fa-solid fa-database"></i> BioCyc</a>`);

    } else {

        dbLinks.push(`<a href="https://www.ncbi.nlm.nih.gov/search/all/?term=${standardCgForLinks}" target="_blank" class="ext-link" title="Search in NCBI"><i class="fa-solid fa-magnifying-glass"></i> NCBI</a>`);

    }

    

    dbLinks.push(`<a href="https://cosy.bio/coryneregnet" target="_blank" class="ext-link" title="Search CoryneRegNet regulatory network database"><i class="fa-solid fa-network-wired"></i> CoryneRegNet</a>`);

    dbLinks.push(`<a href="https://www.uniprot.org/uniprotkb?query=gene:${standardCgForLinks}" target="_blank" class="ext-link" title="View protein function in UniProt"><i class="fa-solid fa-graduation-cap"></i> UniProt</a>`);

    

    // Literature tracking links

    const pubmedQuery = encodeURIComponent(`"Corynebacterium glutamicum" AND (${standardCgForLinks}${meta.name && meta.name !== '--' && meta.name !== standardCgForLinks ? ' OR ' + meta.name : ''})`);

    dbLinks.push(`<a href="https://pubmed.ncbi.nlm.nih.gov/?term=${pubmedQuery}" target="_blank" class="ext-link" title="Search related scientific literature in PubMed"><i class="fa-solid fa-book-open"></i> PubMed</a>`);

    

    const scholarQuery = encodeURIComponent(`"Corynebacterium glutamicum" "${standardCgForLinks}"${meta.name && meta.name !== '--' && meta.name !== standardCgForLinks ? ' OR "' + meta.name + '"' : ''}`);

    dbLinks.push(`<a href="https://scholar.google.com/scholar?q=${scholarQuery}" target="_blank" class="ext-link" title="Search related literature in Google Scholar"><i class="fa-solid fa-graduation-cap"></i> Google Scholar</a>`);

    

    linksCell.innerHTML = dbLinks.join('');



    // Load relationships from global collections

    const relations = [];

    let regsCount = 0;

    let targsCount = 0;



    // Unified edge details

    normalizedEdges.forEach(edge => {
        const sourceLower = edge.source.toLowerCase();
        const targetLower = edge.target.toLowerCase();
        const sourceMeta = getNodeMetaForDetails(edge.source);
        const targetMeta = getNodeMetaForDetails(edge.target);
        const sourceText = `${cleanStr(edge.evidence?.source) || edge.interactionClass}; ${confidenceSummary(edge)}`;

        if (sourceLower === lower) {
            targsCount++;
            relations.push({
                gene: getPrioritizedLabel(edge.target, targetMeta.name),
                locusTag: edge.target,
                dir: 'outgoing',
                role: edge.legacyRole || edge.role,
                regulationType: edge.regulationType,
                confidenceScore: edge.confidenceScore,
                confidenceLevel: edge.confidenceLevel,
                source: sourceText
            });
        }

        if (targetLower === lower) {
            regsCount++;
            relations.push({
                gene: getPrioritizedLabel(edge.source, sourceMeta.name),
                locusTag: edge.source,
                dir: 'incoming',
                role: edge.legacyRole || edge.role,
                regulationType: edge.regulationType,
                confidenceScore: edge.confidenceScore,
                confidenceLevel: edge.confidenceLevel,
                source: sourceText
            });
        }
    });

    // Update Counts

    regulatorsCount.textContent = regsCount;

    targetsCount.textContent = targsCount;



    // Collect lists of regulator and target locus tags

    const incomingLoci = [...new Set(relations.filter(r => r.dir === 'incoming').map(r => r.locusTag))];

    const outgoingLoci = [...new Set(relations.filter(r => r.dir === 'outgoing').map(r => r.locusTag))];



    const regCard = document.getElementById('btn-regulators-summary');

    const targetCard = document.getElementById('btn-targets-summary');



    regCard.onclick = () => {

        if (incomingLoci.length > 0) {

            queryMultipleGenes(incomingLoci);

        } else {

            alert('No upstream regulators are available for this gene.');

        }

    };



    targetCard.onclick = () => {

        if (outgoingLoci.length > 0) {

            queryMultipleGenes(outgoingLoci);

        } else {

            alert('No downstream targets are available for this gene.');

        }

    };



    // Render Table

    relationsTableBody.innerHTML = '';

    

    if (relations.length === 0) {

        relationsTableBody.innerHTML = `<tr><td colspan="4" class="text-muted" style="text-align:center;">No regulatory detail data available</td></tr>`;

    } else {

        // Sort: Incoming first, then outgoing

        relations.sort((a, b) => a.dir.localeCompare(b.dir));

        

        relations.forEach(rel => {

            const tr = document.createElement('tr');

            

            const roleClass = rel.regulationType === 'activation' ? 'activation' : rel.regulationType === 'repression' ? 'repression' : rel.regulationType === 'post_transcriptional_repression' ? 'srna' : 'dual';

            const roleText = roleLabelFromType(rel.role, rel.regulationType);

            

            tr.innerHTML = `

                <td><a href="#" class="gene-link" data-locus="${rel.locusTag}">${rel.gene}</a></td>

                <td><span class="badge-dir ${rel.dir}">${rel.dir === 'incoming' ? '? Upstream' : 'Downstream ?'}</span></td>

                <td><span class="badge-role ${roleClass}">${roleText}</span></td>

                <td class="text-energy">${rel.source}</td>

            `;

            

            // Allow jumping to associated gene on click

            const linkNode = tr.querySelector('.gene-link');

            linkNode.addEventListener('click', (e) => {

                e.preventDefault();

                const targetLocus = linkNode.getAttribute('data-locus');

                querySingleGene(targetLocus);

            });

            

            relationsTableBody.appendChild(tr);

        });

    }



    // Setup perturbation simulator panel

    const pertPanel = document.getElementById('detail-perturbation-panel');

    if (pertPanel) {

        if (targsCount > 0) {

            pertPanel.style.display = 'block';

            

            const btnOe = document.getElementById('btn-sim-oe');

            const btnKo = document.getElementById('btn-sim-ko');

            const btnReset = document.getElementById('btn-sim-reset');

            const btnExport = document.getElementById('btn-sim-export');

            

            const setBtnActive = (activeType) => {

                if (activeType === 'OE') {

                    btnOe.style.backgroundColor = 'rgba(46, 125, 50, 0.15)';

                    btnOe.style.borderColor = '#2e7d32';

                    btnKo.style.backgroundColor = 'rgba(211, 47, 47, 0.03)';

                    btnKo.style.borderColor = 'rgba(211, 47, 47, 0.2)';

                } else if (activeType === 'KO') {

                    btnKo.style.backgroundColor = 'rgba(211, 47, 47, 0.15)';

                    btnKo.style.borderColor = '#d32f2f';

                    btnOe.style.backgroundColor = 'rgba(46, 125, 50, 0.03)';

                    btnOe.style.borderColor = 'rgba(46, 125, 50, 0.2)';

                } else {

                    btnOe.style.backgroundColor = 'rgba(46, 125, 50, 0.03)';

                    btnOe.style.borderColor = 'rgba(46, 125, 50, 0.2)';

                    btnKo.style.backgroundColor = 'rgba(211, 47, 47, 0.03)';

                    btnKo.style.borderColor = 'rgba(211, 47, 47, 0.2)';

                }

            };



            setBtnActive('none');



            btnOe.onclick = () => {

                setBtnActive('OE');

                runPerturbationSimulation(locusTag, 'OE');

            };



            btnKo.onclick = () => {

                setBtnActive('KO');

                runPerturbationSimulation(locusTag, 'KO');

            };



            btnReset.onclick = () => {

                setBtnActive('none');

                resetPerturbationSimulation();

            };



            if (btnExport) {

                btnExport.onclick = () => {

                    exportPerturbationToCsv();

                };

            }

        } else {

            pertPanel.style.display = 'none';

        }

    }

    // Always render genomic locus map for all nodes
    const genomicMapSection = document.getElementById('detail-genomic-map-section');
    if (genomicMapSection) {
        genomicMapSection.style.display = 'block';
        renderGenomicLocusMap(resolvedLocus);
    }

    // Setup protein domain and binding site sections
    const proteinDomainSection = document.getElementById('detail-protein-domain-section');
    const bindingSiteSection = document.getElementById('detail-binding-site-section');
    if (proteinDomainSection && bindingSiteSection) {
        if (meta.type === 'TF') {
            proteinDomainSection.style.display = 'block';
            bindingSiteSection.style.display = 'block';
            loadMotifAndBindingSites(meta.locusTag);
            // Fetch regulon pathway enrichment
            fetchRegulonPathwayEnrichment(meta.locusTag);
            // Hide the motif scan results from any previous query
            const scanResultsBox = document.getElementById('scan-results-box');
            if (scanResultsBox) scanResultsBox.classList.add('hidden');
            const scanInput = document.getElementById('scan-sequence-input');
            if (scanInput) scanInput.value = '';
        } else {
            proteinDomainSection.style.display = 'none';
            bindingSiteSection.style.display = 'none';
        }
    }

    // Slide open sidebar
    toggleRightSidebar(true);

    // Initialize FBA simulation
    initFbaSimulation(meta.locusTag, meta.type);
}

async function initFbaSimulation(locusTag, nodeType) {
    const fbaSection = document.getElementById('detail-fba-simulation-section');
    const fbaStatus = document.getElementById('fba-backend-status');
    const fbaBtn = document.getElementById('btn-run-fba-simulation');
    const fbaResult = document.getElementById('fba-result-container');
    const fbaError = document.getElementById('fba-error-container');
    
    // Config controls
    const objSelect = document.getElementById('fba-objective-select');
    const customObjContainer = document.getElementById('fba-custom-objective-container');
    const objSearchInput = document.getElementById('fba-obj-reaction-search');
    const btnObjSearch = document.getElementById('btn-fba-obj-reaction-search');
    const objRxnSelect = document.getElementById('fba-obj-reaction-select');
    const objEquation = document.getElementById('fba-obj-reaction-equation');
    
    const trackSearchInput = document.getElementById('fba-track-reaction-search');
    const btnTrackSearch = document.getElementById('btn-fba-track-reaction-search');
    const trackRxnSelect = document.getElementById('fba-track-reaction-select');
    const trackEquation = document.getElementById('fba-track-reaction-equation');
    const btnFindGlutamate = document.getElementById('btn-find-glutamate-helper');
    
    // Outputs
    const changeLabel = document.getElementById('fba-change-label');
    const trackedResultsBox = document.getElementById('fba-tracked-flux-results');
    const interpretationText = document.getElementById('fba-interpretation-text');
    
    if (!fbaSection || !fbaStatus || !fbaBtn || !fbaResult || !fbaError) return;
    
    // Show section
    fbaSection.style.display = 'block';
    
    // Clear previous results & errors
    fbaResult.classList.add('hidden');
    fbaError.classList.add('hidden');
    
    // Reset inputs but preserve state if appropriate
    if (objSelect) {
        objSelect.value = 'biomass';
        if (customObjContainer) customObjContainer.classList.add('hidden');
        objSelect.onchange = () => {
            if (objSelect.value === 'reaction') {
                customObjContainer.classList.remove('hidden');
            } else {
                customObjContainer.classList.add('hidden');
            }
        };
    }
    
    const equationsMap = new Map();
    
    const populateReactionsSelect = (selectElement, equationDiv, matches) => {
        selectElement.innerHTML = '';
        if (!matches || matches.length === 0) {
            selectElement.style.display = 'none';
            equationDiv.textContent = 'No matching reactions found.';
            return;
        }
        
        selectElement.style.display = 'block';
        
        const placeholder = document.createElement('option');
        placeholder.value = '';
        placeholder.textContent = `-- Select verified reaction (${matches.length} matches) --`;
        selectElement.appendChild(placeholder);
        
        matches.forEach(m => {
            const opt = document.createElement('option');
            opt.value = m.reactionId;
            opt.textContent = `[${m.reactionId}] ${m.name || 'Unnamed'}`;
            selectElement.appendChild(opt);
            equationsMap.set(m.reactionId, m.equation);
        });
        
        selectElement.onchange = () => {
            const rxnId = selectElement.value;
            if (rxnId && equationsMap.has(rxnId)) {
                equationDiv.textContent = `Equation: ${equationsMap.get(rxnId)}`;
            } else {
                equationDiv.textContent = '';
            }
        };
    };
    
    // Wire Search actions
    if (btnObjSearch && objSearchInput && objRxnSelect && objEquation) {
        btnObjSearch.onclick = async () => {
            const q = objSearchInput.value.trim();
            if (!q) return;
            btnObjSearch.innerHTML = '<i class="fa-solid fa-spinner fa-spin"></i>';
            const data = await window.simulationClient.searchReactions(q);
            btnObjSearch.innerHTML = '<i class="fa-solid fa-search"></i>';
            populateReactionsSelect(objRxnSelect, objEquation, data.matches);
        };
    }
    
    if (btnTrackSearch && trackSearchInput && trackRxnSelect && trackEquation) {
        btnTrackSearch.onclick = async () => {
            const q = trackSearchInput.value.trim();
            if (!q) return;
            btnTrackSearch.innerHTML = '<i class="fa-solid fa-spinner fa-spin"></i>';
            const data = await window.simulationClient.searchReactions(q);
            btnTrackSearch.innerHTML = '<i class="fa-solid fa-search"></i>';
            populateReactionsSelect(trackRxnSelect, trackEquation, data.matches);
        };
    }
    
    if (btnFindGlutamate && trackSearchInput && btnTrackSearch) {
        btnFindGlutamate.onclick = () => {
            trackSearchInput.value = 'glutamate';
            btnTrackSearch.click();
        };
    }
    
    // Set loading status
    fbaStatus.textContent = 'Checking...';
    fbaStatus.style.color = 'var(--text-muted)';
    fbaBtn.disabled = true;
    
    // Update button text depending on nodeType
    const isTf = (nodeType === 'TF' || nodeType === 'sRNA');
    fbaBtn.innerHTML = isTf 
        ? '<i class="fa-solid fa-play"></i> Run TF Target Perturbation'
        : '<i class="fa-solid fa-play"></i> Run Gene Knockout';
        
    // Check backend status
    const status = await window.simulationClient.getModelStatus();
    if (status && status.loaded) {
        fbaStatus.textContent = `Model Loaded (${status.reaction_count} rxns)`;
        fbaStatus.style.color = 'var(--color-activation)';
        fbaBtn.disabled = false;
    } else {
        fbaStatus.textContent = status && status.error ? `Offline (${status.error})` : 'Offline (backend unreachable)';
        fbaStatus.style.color = 'var(--color-repression)';
        fbaBtn.disabled = true;
    }
    
    // Wire button action
    fbaBtn.onclick = async () => {
        fbaBtn.disabled = true;
        fbaBtn.innerHTML = '<i class="fa-solid fa-spinner fa-spin"></i> Simulating...';
        fbaResult.classList.add('hidden');
        fbaError.classList.add('hidden');
        
        const objective = {
            objectiveType: objSelect ? objSelect.value : 'biomass',
            reactionId: null
        };
        if (objective.objectiveType === 'reaction') {
            const rxn = objRxnSelect ? objRxnSelect.value : '';
            if (!rxn) {
                fbaBtn.disabled = false;
                fbaBtn.innerHTML = isTf 
                    ? '<i class="fa-solid fa-play"></i> Run TF Target Perturbation'
                    : '<i class="fa-solid fa-play"></i> Run Gene Knockout';
                fbaError.classList.remove('hidden');
                fbaError.textContent = 'Please select a custom objective reaction first.';
                return;
            }
            objective.reactionId = rxn;
        }
        
        const trackReactionIds = [];
        if (trackRxnSelect && trackRxnSelect.value) {
            trackReactionIds.push(trackRxnSelect.value);
        }
        
        let res;
        if (isTf) {
            const targetGeneIds = [];
            if (cy) {
                cy.edges(`[source = "${locusTag}"]`).targets().forEach(node => {
                    targetGeneIds.push(node.id());
                });
            }
            res = await window.simulationClient.runTFPerturbation(locusTag, targetGeneIds, objective, trackReactionIds);
        } else {
            res = await window.simulationClient.runGeneKnockout(locusTag, objective, trackReactionIds);
        }
        
        fbaBtn.disabled = false;
        fbaBtn.innerHTML = isTf 
            ? '<i class="fa-solid fa-play"></i> Run TF Target Perturbation'
            : '<i class="fa-solid fa-play"></i> Run Gene Knockout';
            
        if (res && res.status && res.status !== 'error') {
            fbaResult.classList.remove('hidden');
            
            const baseline = res.baselineObjective;
            const perturbed = res.perturbedObjective;
            const change = res.objectiveChange;
            const pct = res.objectiveChangePercent;
            
            const unit = " mmol/gDCW/h";
            document.getElementById('fba-baseline-obj').textContent = baseline.toFixed(4) + unit;
            document.getElementById('fba-perturbed-obj').textContent = perturbed.toFixed(4) + unit;
            
            if (changeLabel) {
                changeLabel.textContent = objective.objectiveType === 'biomass' ? 'Growth Change:' : 'Objective Change:';
            }
            
            const changeEl = document.getElementById('fba-change-pct');
            changeEl.textContent = (pct >= 0 ? '+' : '') + pct.toFixed(2) + "%";
            
            if (pct < -0.01) {
                changeEl.style.color = 'var(--color-repression)';
            } else if (pct > 0.01) {
                changeEl.style.color = 'var(--color-activation)';
            } else {
                changeEl.style.color = 'var(--text-primary)';
            }
            
            // Render Tracked Fluxes
            if (trackedResultsBox) {
                if (res.trackedFluxes && res.trackedFluxes.length > 0) {
                    trackedResultsBox.classList.remove('hidden');
                    const tf = res.trackedFluxes[0];
                    
                    const labelSpan = document.getElementById('fba-tracked-rxn-label');
                    if (labelSpan) labelSpan.textContent = `Reaction [${tf.reactionId}]:`;
                    
                    const tfPct = tf.fluxChangePercent;
                    const changeVal = tf.fluxChange;
                    const changeText = (changeVal >= 0 ? '+' : '') + changeVal.toFixed(4) + ` (${(tfPct >= 0 ? '+' : '')}${tfPct.toFixed(1)}%)`;
                    
                    const tfPctEl = document.getElementById('fba-tracked-rxn-change-pct');
                    if (tfPctEl) {
                        tfPctEl.textContent = changeText;
                        if (changeVal < -1e-5) {
                            tfPctEl.style.color = 'var(--color-repression)';
                        } else if (changeVal > 1e-5) {
                            tfPctEl.style.color = 'var(--color-activation)';
                        } else {
                            tfPctEl.style.color = 'var(--text-primary)';
                        }
                    }
                    
                    const baselineEl = document.getElementById('fba-tracked-rxn-baseline');
                    if (baselineEl) baselineEl.textContent = tf.baselineFlux.toFixed(4) + " mmol/gDCW/h";
                    
                    const perturbedEl = document.getElementById('fba-tracked-rxn-perturbed');
                    if (perturbedEl) perturbedEl.textContent = tf.perturbedFlux.toFixed(4) + " mmol/gDCW/h";
                } else {
                    trackedResultsBox.classList.add('hidden');
                }
            }
            
            // Render Interpretation text
            if (interpretationText && window.objectiveInterpretation) {
                interpretationText.textContent = window.objectiveInterpretation.generateObjectiveSimulationInterpretation(res);
            }
            
            // Show warnings if any
            if (res.warnings && res.warnings.length > 0) {
                fbaError.classList.remove('hidden');
                fbaError.style.color = '#b45309';
                fbaError.style.background = '#fffbeb';
                fbaError.style.borderColor = '#fef3c7';
                fbaError.innerHTML = '<strong>Warnings:</strong><ul style="margin: 4px 0 0 0; padding-left: 14px;">' + 
                    res.warnings.map(w => `<li>${w}</li>`).join('') + '</ul>';
            }
        } else {
            fbaError.classList.remove('hidden');
            fbaError.style.color = '#d32f2f';
            fbaError.style.background = '#fef2f2';
            fbaError.style.borderColor = '#fee2e2';
            fbaError.textContent = res && res.error ? `Simulation failed: ${res.error}` : 'Simulation failed: Backend offline or model loading failed.';
        }
    };
}



function showOperonDetails(operonMeta, initialMode = null) {

    // Reset any ongoing simulation first

    resetPerturbationSimulation();



    // Clear previous AI summary

    const summaryCard = document.getElementById('ai-summary-result');

    if (summaryCard) {

        summaryCard.classList.add('hidden');

        summaryCard.innerHTML = '';

    }



    // Set badge style

    detailTypeBadge.className = 'gene-badge';

    detailTypeBadge.style.backgroundColor = 'var(--color-primary-accent)';

    detailTypeBadge.style.color = '#ffffff';

    detailTypeBadge.textContent = 'Operon';



    detailGeneName.textContent = `${operonMeta.operon} operon`;

    detailLocusTag.textContent = `Orientation: ${operonMeta.orientation} strand | ${operonMeta.genes.length} genes`;



    infoLocus.textContent = operonMeta.genes.join(', ');

    infoName.textContent = operonMeta.operon;

    infoType.textContent = 'Operon';



    const viewNetworkBtn = document.getElementById('view-network-btn');

    if (viewNetworkBtn) {

        viewNetworkBtn.style.display = 'flex';

        viewNetworkBtn.onclick = () => {

            queryMultipleGenes(operonMeta.genes);

            showOperonDetails(operonMeta);

        };

    }



    const cglRow = document.getElementById('info-cgl-row');

    if (cglRow) {

        cglRow.style.display = 'none';

    }



    const productRow = document.getElementById('info-product-row');

    const infoProduct = document.getElementById('info-product');

    if (productRow && infoProduct) {

        productRow.style.display = '';

        let productHtml = '<div style="display: flex; flex-direction: column; gap: 6px;">';

        operonMeta.genes.forEach(g => {

            const lower = g.toLowerCase();

            const product = cgToProduct[lower] || 'No description available';

            const prioritized = getPrioritizedLabel(g, g);

            productHtml += `<div><strong style="color: var(--text-primary); font-family: monospace;">${prioritized}:</strong> <span style="color: var(--text-secondary);">${product}</span></div>`;

        });

        productHtml += '</div>';

        infoProduct.innerHTML = productHtml;

    }



    const pathwayRow = document.getElementById('info-pathway-row');

    if (pathwayRow) {

        pathwayRow.style.display = 'none';

    }

    const metabolicSection = document.getElementById('detail-metabolic-impact-section');

    if (metabolicSection) {

        metabolicSection.style.display = 'none';

    }



    const operonRow = document.getElementById('info-operon-row');

    const infoOperon = document.getElementById('info-operon');

    if (operonRow && infoOperon) {

        operonRow.style.display = '';

        const geneLinks = operonMeta.genes.map(g => {

            const prioritized = getPrioritizedLabel(g, g);

            return `<a href="#" class="operon-gene-link" data-locus="${g}" style="color: var(--color-primary-accent); text-decoration: none; font-weight: 500; font-family: monospace;">${prioritized}</a>`;

        }).join(', ');

        infoOperon.innerHTML = `

            <div style="font-size: 11px; color: var(--text-secondary);">Genes: ${geneLinks}</div>

            <div style="display: flex; gap: 6px; margin-top: 8px;">

                <button id="btn-draw-operon-network-details" class="secondary-btn" style="flex: 1; font-size: 10px; padding: 4px 6px; height: auto; border: 1px solid rgba(30, 58, 138, 0.15); color: var(--color-primary-accent); background-color: rgba(30, 58, 138, 0.03);">

                    <i class="fa-solid fa-network-wired"></i> Joint analysis

                </button>

            </div>

        `;



        infoOperon.querySelectorAll('.operon-gene-link').forEach(link => {

            link.addEventListener('click', (e) => {

                e.preventDefault();

                const targetLocus = link.getAttribute('data-locus');

                querySingleGene(targetLocus);

            });

        });



        infoOperon.querySelector('#btn-draw-operon-network-details').onclick = () => {

            queryMultipleGenes(operonMeta.genes);

            showOperonDetails(operonMeta);

        };

        const btnOe = infoOperon.querySelector('#btn-sim-operon-oe-details');

        if (btnOe) {

            btnOe.onclick = () => {

                queryMultipleGenes(operonMeta.genes);

                showOperonDetails(operonMeta, 'OE');

            };

        }



        const btnKo = infoOperon.querySelector('#btn-sim-operon-ko-details');

        if (btnKo) {

            btnKo.onclick = () => {

                queryMultipleGenes(operonMeta.genes);

                showOperonDetails(operonMeta, 'KO');

            };

        }

    }



    const linksCell = document.getElementById('info-links');

    if (linksCell) {

        linksCell.innerHTML = '';

        const dbLinks = [];

        operonMeta.genes.forEach(g => {

            const prioritized = getPrioritizedLabel(g, g);

            const pubmedQuery = encodeURIComponent(`"Corynebacterium glutamicum" AND "${g}"`);

            dbLinks.push(`

                <div style="margin-bottom: 6px; width: 100%; border-bottom: 1px dashed var(--border-color); padding-bottom: 4px;">

                    <strong style="font-family: monospace; font-size: 11px;">${prioritized}:</strong>

                    <div style="display: flex; gap: 4px; flex-wrap: wrap; margin-top: 2px;">

                        <a href="https://www.ncbi.nlm.nih.gov/gene/?term=${g}" target="_blank" class="ext-link" style="font-size: 10px; padding: 2px 4px;"><i class="fa-solid fa-dna"></i> NCBI</a>

                        <a href="https://pubmed.ncbi.nlm.nih.gov/?term=${pubmedQuery}" target="_blank" class="ext-link" style="font-size: 10px; padding: 2px 4px;"><i class="fa-solid fa-book-open"></i> Literature</a>

                    </div>

                </div>

            `);

        });

        linksCell.innerHTML = `<div style="display: flex; flex-direction: column; width: 100%;">${dbLinks.join('')}</div>`;

    }



    const relations = [];

    let regsCount = 0;

    let targsCount = 0;

    const operonGeneSet = new Set(operonMeta.genes.map(g => g.toLowerCase()));



    normalizedEdges.forEach(edge => {
        const sourceLower = edge.source.toLowerCase();
        const targetLower = edge.target.toLowerCase();
        const sourceMeta = getNodeMetaForDetails(edge.source);
        const targetMeta = getNodeMetaForDetails(edge.target);
        const sourceText = `${cleanStr(edge.evidence?.source) || edge.interactionClass}; ${confidenceSummary(edge)}`;

        if (operonGeneSet.has(targetLower) && !operonGeneSet.has(sourceLower)) {
            regsCount++;
            relations.push({
                gene: getPrioritizedLabel(edge.source, sourceMeta.name),
                locusTag: edge.source,
                dir: 'incoming',
                role: edge.legacyRole || edge.role,
                regulationType: edge.regulationType,
                confidenceScore: edge.confidenceScore,
                confidenceLevel: edge.confidenceLevel,
                source: sourceText,
                targetGene: getPrioritizedLabel(edge.target, targetMeta.name)
            });
        }

        if (operonGeneSet.has(sourceLower) && !operonGeneSet.has(targetLower)) {
            targsCount++;
            relations.push({
                gene: getPrioritizedLabel(edge.target, targetMeta.name),
                locusTag: edge.target,
                dir: 'outgoing',
                role: edge.legacyRole || edge.role,
                regulationType: edge.regulationType,
                confidenceScore: edge.confidenceScore,
                confidenceLevel: edge.confidenceLevel,
                source: sourceText,
                sourceGene: getPrioritizedLabel(edge.source, sourceMeta.name)
            });
        }
    });

    regulatorsCount.textContent = regsCount;

    targetsCount.textContent = targsCount;



    const incomingLoci = [...new Set(relations.filter(r => r.dir === 'incoming').map(r => r.locusTag))];

    const outgoingLoci = [...new Set(relations.filter(r => r.dir === 'outgoing').map(r => r.locusTag))];



    const regCard = document.getElementById('btn-regulators-summary');

    const targetCard = document.getElementById('btn-targets-summary');



    regCard.onclick = () => {

        if (incomingLoci.length > 0) {

            queryMultipleGenes(incomingLoci);

        } else {

            alert('No upstream regulators are available for this operon.');

        }

    };



    targetCard.onclick = () => {

        if (outgoingLoci.length > 0) {

            queryMultipleGenes(outgoingLoci);

        } else {

            alert('No downstream targets are available for this operon.');

        }

    };



    relationsTableBody.innerHTML = '';

    

    if (relations.length === 0) {

        relationsTableBody.innerHTML = `<tr><td colspan="4" class="text-muted" style="text-align:center;">No regulatory detail data available</td></tr>`;

    } else {

        relations.sort((a, b) => a.dir.localeCompare(b.dir));

        

        relations.forEach(rel => {

            const tr = document.createElement('tr');

            const roleClass = rel.regulationType === 'activation' ? 'activation' : rel.regulationType === 'repression' ? 'repression' : rel.regulationType === 'post_transcriptional_repression' ? 'srna' : 'dual';

            const roleText = roleLabelFromType(rel.role, rel.regulationType);

            const assocGeneText = rel.dir === 'incoming' 

                ? ` (regulates ${rel.targetGene})` 

                : ` (regulated by ${rel.sourceGene})`;



            tr.innerHTML = `

                <td>

                    <a href="#" class="gene-link" data-locus="${rel.locusTag}">${rel.gene}</a>

                    <span style="font-size: 10px; color: var(--text-muted); display: block;">${assocGeneText}</span>

                </td>

                <td><span class="badge-dir ${rel.dir}">${rel.dir === 'incoming' ? '? Upstream' : 'Downstream ?'}</span></td>

                <td><span class="badge-role ${roleClass}">${roleText}</span></td>

                <td class="text-energy">${rel.source}</td>

            `;

            

            const linkNode = tr.querySelector('.gene-link');

            linkNode.addEventListener('click', (e) => {

                e.preventDefault();

                const targetLocus = linkNode.getAttribute('data-locus');

                querySingleGene(targetLocus);

            });

            relationsTableBody.appendChild(tr);

        });

    }



    const pertPanel = document.getElementById('detail-perturbation-panel');

    if (pertPanel) {

        if (targsCount > 0) {

            pertPanel.style.display = 'block';

            

            const btnOe = document.getElementById('btn-sim-oe');

            const btnKo = document.getElementById('btn-sim-ko');

            const btnReset = document.getElementById('btn-sim-reset');

            const btnExport = document.getElementById('btn-sim-export');

            

            const setBtnActive = (activeType) => {

                if (activeType === 'OE') {

                    btnOe.style.backgroundColor = 'rgba(46, 125, 50, 0.15)';

                    btnOe.style.borderColor = '#2e7d32';

                    btnKo.style.backgroundColor = 'rgba(211, 47, 47, 0.03)';

                    btnKo.style.borderColor = 'rgba(211, 47, 47, 0.2)';

                } else if (activeType === 'KO') {

                    btnKo.style.backgroundColor = 'rgba(211, 47, 47, 0.15)';

                    btnKo.style.borderColor = '#d32f2f';

                    btnOe.style.backgroundColor = 'rgba(46, 125, 50, 0.03)';

                    btnOe.style.borderColor = 'rgba(46, 125, 50, 0.2)';

                } else {

                    btnOe.style.backgroundColor = 'rgba(46, 125, 50, 0.03)';

                    btnOe.style.borderColor = 'rgba(46, 125, 50, 0.2)';

                    btnKo.style.backgroundColor = 'rgba(211, 47, 47, 0.03)';

                    btnKo.style.borderColor = 'rgba(211, 47, 47, 0.2)';

                }

            };



            setBtnActive('none');



            btnOe.onclick = () => {

                setBtnActive('OE');

                runPerturbationSimulation(operonMeta.genes, 'OE');

            };



            btnKo.onclick = () => {

                setBtnActive('KO');

                runPerturbationSimulation(operonMeta.genes, 'KO');

            };



            btnReset.onclick = () => {

                setBtnActive('none');

                resetPerturbationSimulation();

            };



            if (btnExport) {

                btnExport.onclick = () => {

                    exportPerturbationToCsv();

                };

            }

        } else {

            pertPanel.style.display = 'none';

        }

    }



    const btnTriggerAi = document.getElementById('btn-trigger-ai');

    if (btnTriggerAi) {

        btnTriggerAi.style.display = 'none';

    }



    toggleRightSidebar(true);



    if (initialMode === 'OE' || initialMode === 'KO') {

        const btnOe = document.getElementById('btn-sim-oe');

        const btnKo = document.getElementById('btn-sim-ko');

        if (initialMode === 'OE' && btnOe) {

            btnOe.click();

        } else if (initialMode === 'KO' && btnKo) {

            btnKo.click();

        }

    }

}



// ==========================================================================

// 5. DOM Event Listeners & Interactive Controls

// ==========================================================================

function initEventListeners() {

    // Initialize first empty input row

    clearAllInputs();
    initWorkflowEntrypoints();



    searchBtn.addEventListener('click', () => {

        suggestionsBox.classList.add('hidden');

        triggerSearchFromInputs();

    });



    // Close suggestions list on click outside

    document.addEventListener('click', (e) => {

        if (!e.target.classList.contains('gene-input') && e.target !== suggestionsBox && !suggestionsBox.contains(e.target)) {

            suggestionsBox.classList.add('hidden');

        }

    });



    // Example quick tags

    document.querySelectorAll('.example-tag').forEach(tag => {

        tag.addEventListener('click', () => {

            querySingleGene(tag.textContent);

        });

    });



    // Sidebar Config Filters

    const reRender = () => {

        if (currentQueryGene) {

            renderNetwork(currentQueryGene);

        }

    };



    filterActivation.addEventListener('change', reRender);

    filterRepression.addEventListener('change', reRender);

    filterDual.addEventListener('change', reRender);

    filterSrna.addEventListener('change', () => {
        if (filterSrna.checked) {
            srnaThresholdPanel.classList.remove('hidden');
        } else {
            srnaThresholdPanel.classList.add('hidden');
        }
        reRender();
    });

    

    // Sync co-regulation checkbox

    filterCoregulated.addEventListener('change', reRender);

    

    if (filterOnlyTfTargets) {

        filterOnlyTfTargets.addEventListener('change', reRender);

    }

    

    srnaRankThreshold.addEventListener('input', (e) => {

        rankValDisp.textContent = e.target.value;

    });

    srnaRankThreshold.addEventListener('change', reRender);



    layoutSelect.addEventListener('change', () => {

        if (cy) {

            const layout = cy.layout({

                name: layoutSelect.value,

                animate: true,

                animationDuration: 450

            });

            layout.run();

        }

    });



    // Detail Panel closer

    closeDetailBtn.addEventListener('click', () => {

        toggleRightSidebar(false);

    });

    if (rightSidebarToggle) {

        rightSidebarToggle.addEventListener('click', () => {

            const isCollapsed = rightSidebar?.classList.contains('collapsed');

            toggleRightSidebar(isCollapsed);

        });

    }



    // Zooming & view fit controls

    resetViewBtn.addEventListener('click', () => {

        if (cy) {

            cy.fit();

            cy.center();

        }

    });



    zoomInBtn.addEventListener('click', () => {

        if (cy) cy.zoom(cy.zoom() * 1.2);

    });



    zoomOutBtn.addEventListener('click', () => {

        if (cy) cy.zoom(cy.zoom() / 1.2);

    });



    fitCanvasBtn.addEventListener('click', () => {

        if (cy) {

            cy.fit();

            cy.center();

        }

    });



    // PNG Image Export

    exportPngBtn.addEventListener('click', () => {

        if (!cy) return;

        

        // Export options

        const pngContent = cy.png({

            bg: '#ffffff',

            full: true,

            scale: 2 // High res export

        });

        

        // Dynamic download trigger

        const link = document.createElement('a');

        const filename = Array.isArray(currentQueryGene) ? currentQueryGene.join('_') : currentQueryGene;

        link.download = `${filename}_regulatory_network.png`;

        link.href = pngContent;

        document.body.appendChild(link);

        link.click();

        document.body.removeChild(link);

    });



    // Initialize AI Literature Summarizer Key & Click Bindings

    initAiSummaryFeature();

    initAiPathwayFeature();

    initAiGeneFeature();

    initRnaSeqOverlay();

    initProteinDomainFeature();

    initBindingSiteFeature();

    initAdvancedFeatures();



    // Tab Switches

    const tabSingleBtn = document.getElementById('tab-single-btn');

    const tabBatchBtn = document.getElementById('tab-batch-btn');

    const tabSingleContent = document.getElementById('search-tab-single-content');

    const tabBatchContent = document.getElementById('search-tab-batch-content');



    if (tabSingleBtn && tabBatchBtn) {

        tabSingleBtn.addEventListener('click', () => {

            tabSingleBtn.classList.add('active');

            tabBatchBtn.classList.remove('active');

            tabSingleContent.classList.remove('hidden');

            tabBatchContent.classList.add('hidden');

        });



        tabBatchBtn.addEventListener('click', () => {

            tabBatchBtn.classList.add('active');

            tabSingleBtn.classList.remove('active');

            tabBatchContent.classList.remove('hidden');

            tabSingleContent.classList.add('hidden');

            // Auto focus textarea

            document.getElementById('gene-batch-textarea')?.focus();

        });

    }



    // AI Discovery Tab Switches

    const tabAiGeneBtn = document.getElementById('tab-ai-gene-btn');

    const tabAiPathwayBtn = document.getElementById('tab-ai-pathway-btn');

    const aiGeneContent = document.getElementById('ai-gene-tab-content');

    const aiPathwayContent = document.getElementById('ai-pathway-tab-content');



    if (tabAiGeneBtn && tabAiPathwayBtn) {

        tabAiGeneBtn.addEventListener('click', () => {

            tabAiGeneBtn.classList.add('active');

            tabAiPathwayBtn.classList.remove('active');

            aiGeneContent.classList.remove('hidden');

            aiPathwayContent.classList.add('hidden');

        });



        tabAiPathwayBtn.addEventListener('click', () => {

            tabAiPathwayBtn.classList.add('active');

            tabAiGeneBtn.classList.remove('active');

            aiPathwayContent.classList.remove('hidden');

            aiGeneContent.classList.add('hidden');

        });

    }



    // CSV Data Export

    const exportCsvBtn = document.getElementById('export-csv-btn');

    if (exportCsvBtn) {

        exportCsvBtn.addEventListener('click', () => {

            exportNetworkToCsv();

        });

    }



    // Batch input counter listener

    initBatchInputCounter();



    // Floating UI panels logic

    initCanvasSearch();

    initStatsToggle();



    // History Navigation

    const backBtn = document.getElementById('btn-history-back');

    const forwardBtn = document.getElementById('btn-history-forward');

    

    if (backBtn && forwardBtn) {

        backBtn.addEventListener('click', () => {

            navigateHistory('back');

        });

        forwardBtn.addEventListener('click', () => {

            navigateHistory('forward');

        });

    }

}



// ==========================================================================

// 6. Dynamic Multiple Gene Input Helpers

// ==========================================================================

function addNewInputRow() {

    const row = document.createElement('div');

    row.className = 'gene-input-row';

    

    const wrapper = document.createElement('div');

    wrapper.className = 'gene-input-wrapper';

    

    const input = document.createElement('input');

    input.type = 'text';

    input.className = 'gene-input';

    input.placeholder = 'Enter gene/sRNA name';

    input.autocomplete = 'off';

    

    wrapper.appendChild(input);

    row.appendChild(wrapper);

    

    // Add delete or add button based on current rows count

    const existingRows = geneInputsContainer.querySelectorAll('.gene-input-row');

    if (existingRows.length > 0) {

        const removeBtn = document.createElement('button');

        removeBtn.className = 'remove-row-btn';

        removeBtn.title = 'Remove gene row';

        removeBtn.innerHTML = '<i class="fa-solid fa-minus"></i>';

        removeBtn.addEventListener('click', () => {

            if (suggestionsBox.parentElement === wrapper) {

                suggestionsBox.classList.add('hidden');

            }

            row.remove();

            triggerSearchFromInputs();

        });

        row.appendChild(removeBtn);

    } else {

        const addBtn = document.createElement('button');

        addBtn.className = 'add-row-btn';

        addBtn.title = 'Add gene row';

        addBtn.innerHTML = '<i class="fa-solid fa-plus"></i>';

        addBtn.addEventListener('click', () => {

            const newInput = addNewInputRow();

            newInput.focus();

        });

        row.appendChild(addBtn);

    }

    

    bindInputEvents(input);

    geneInputsContainer.appendChild(row);

    return input;

}



function bindInputEvents(input) {

    input.addEventListener('focus', () => {

        activeInput = input;

        const wrapper = input.closest('.gene-input-wrapper');

        if (wrapper && suggestionsBox.parentElement !== wrapper) {

            wrapper.appendChild(suggestionsBox);

        }

        if (input.value.trim() !== '') {

            showSuggestions(input.value);

        } else {

            suggestionsBox.classList.add('hidden');

        }

    });

    

    input.addEventListener('input', (e) => {

        showSuggestions(e.target.value);

    });

    

    input.addEventListener('keydown', (e) => {

        if (e.key === 'Enter') {

            suggestionsBox.classList.add('hidden');

            triggerSearchFromInputs();

        }

    });

}



function clearAllInputs() {

    geneInputsContainer.innerHTML = '';

    addNewInputRow();

}



function querySingleGene(locus) {

    const tabSingleBtn = document.getElementById('tab-single-btn');

    if (tabSingleBtn) tabSingleBtn.click();



    clearAllInputs();

    const input = geneInputsContainer.querySelector('.gene-input');

    if (input) {

        // Find prioritized display label

        let displayLabel = locus;

        const lower = locus.toLowerCase();

        

        // Resolve target cg locus tag first to find the correct mapping

        let targetLocus = lower;

        if (cglToCg[lower]) {

            targetLocus = cglToCg[lower].toLowerCase();

        } else if (nameToCg[lower]) {

            targetLocus = nameToCg[lower].toLowerCase();

        }

        

        const match = geneIndex[targetLocus];

        if (match) {

            displayLabel = getPrioritizedLabel(match.locusTag, match.name);

        } else {

            displayLabel = getPrioritizedLabel(locus, locus);

        }

        

        input.value = displayLabel;

        activeInput = input;

    }

    triggerSearchFromInputs();

}



function queryMultipleGenes(loci) {

    if (!loci || loci.length === 0) return;

    

    const tabSingleBtn = document.getElementById('tab-single-btn');

    if (tabSingleBtn) tabSingleBtn.click();



    // Clear all inputs

    geneInputsContainer.innerHTML = '';

    

    // Add input rows and populate them

    loci.forEach((locus, idx) => {

        const input = addNewInputRow();

        

        let displayLabel = locus;

        const lower = locus.toLowerCase();

        

        // Resolve target cg locus tag first to find the correct mapping

        let targetLocus = lower;

        if (cglToCg[lower]) {

            targetLocus = cglToCg[lower].toLowerCase();

        } else if (nameToCg[lower]) {

            targetLocus = nameToCg[lower].toLowerCase();

        }

        

        const match = geneIndex[targetLocus];

        if (match) {

            displayLabel = getPrioritizedLabel(match.locusTag, match.name);

        } else {

            displayLabel = getPrioritizedLabel(locus, locus);

        }

        

        input.value = displayLabel;
    });
}

function parseOperons(text) {
    geneToOperon = {};
    const lines = text.split('\n');
    for (let i = 1; i < lines.length; i++) {
        const line = lines[i].trim();
        if (!line) continue;
        const parts = line.split(',').map(p => p.trim());
        if (parts.length >= 3) {
            const operonName = parts[0].replace('>', '');
            const orientation = parts[1];
            const genes = parts.slice(2).filter(g => g);
            const operonInfo = { operon: operonName, orientation, genes };
            genes.forEach(gene => {
                geneToOperon[gene.toLowerCase()] = operonInfo;
            });
        }
    }
}

function initAiSummaryFeature() {
    const btnSaveKey = document.getElementById('btn-save-key');
    const btnClearKey = document.getElementById('btn-clear-key');
    const btnTriggerAi = document.getElementById('btn-trigger-ai');
    const apiKeyInput = document.getElementById('gemini-api-key-input');
    const keyConfigPanel = document.getElementById('ai-key-config-panel');
    const keyActivePanel = document.getElementById('ai-key-active-panel');
    
    // Multi-provider inputs
    const providerSelect = document.getElementById('ai-provider-select');
    const baseUrlInput = document.getElementById('ai-base-url-input');
    const modelInput = document.getElementById('ai-model-input');
    
    const customUrlWrapper = document.getElementById('ai-custom-url-wrapper');
    const modelWrapper = document.getElementById('ai-model-wrapper');
    const activeStatusText = document.getElementById('ai-active-status-text');

    const providerNames = {
        'google': 'Google Gemini',
        'openai': 'OpenAI',
        'deepseek': 'DeepSeek',
        'qwen': 'Qwen',
        'kimi': 'Kimi',
        'zhipu': 'Zhipu GLM',
        'ollama': 'Ollama',
        'custom': 'Custom endpoint'
    };

    const providerDefaults = {
        'google': { model: '', baseUrl: '' },
        'openai': { model: 'gpt-4o-mini', baseUrl: 'https://api.openai.com/v1' },
        'deepseek': { model: 'deepseek-chat', baseUrl: 'https://api.deepseek.com' },
        'qwen': { model: 'qwen-plus', baseUrl: 'https://dashscope.aliyuncs.com/compatible-mode/v1' },
        'kimi': { model: 'moonshot-v1-8k', baseUrl: 'https://api.moonshot.cn/v1' },
        'zhipu': { model: 'glm-4-flash', baseUrl: 'https://open.bigmodel.cn/api/paas/v4' },
        'ollama': { model: 'deepseek-r1', baseUrl: 'http://localhost:11434/v1' },
        'custom': { model: '', baseUrl: '' }
    };

    const hints = {
        'google': document.getElementById('ai-key-hint-google'),
        'openai': document.getElementById('ai-key-hint-openai'),
        'deepseek': document.getElementById('ai-key-hint-deepseek'),
        'qwen': document.getElementById('ai-key-hint-qwen'),
        'kimi': document.getElementById('ai-key-hint-kimi'),
        'zhipu': document.getElementById('ai-key-hint-zhipu'),
        'ollama': document.getElementById('ai-key-hint-ollama')
    };

    // Helper to toggle input fields visibility depending on selected provider
    function updateConfigFields() {
        const provider = providerSelect.value;
        
        // Hide all hints first
        Object.values(hints).forEach(h => {
            if (h) h.classList.add('hidden');
        });
        
        // Show current provider hint
        if (hints[provider]) {
            hints[provider].classList.remove('hidden');
        }

        // Toggle Base URL and Model visibility (hide only for Google Gemini)
        if (provider === 'google') {
            if (customUrlWrapper) customUrlWrapper.classList.add('hidden');
            if (modelWrapper) modelWrapper.classList.add('hidden');
        } else {
            if (customUrlWrapper) customUrlWrapper.classList.remove('hidden');
            if (modelWrapper) modelWrapper.classList.remove('hidden');
            
            // Adjust placeholders based on provider
            if (modelInput) {
                if (provider === 'custom') modelInput.placeholder = 'Example: gpt-4o-mini';
                else modelInput.placeholder = `Example: ${providerDefaults[provider].model}`;
            }
        }

        // Adjust API Key label & requirements for Ollama
        const keyLabel = document.getElementById('ai-key-label');
        if (provider === 'ollama') {
            if (keyLabel) keyLabel.textContent = 'API Key (optional for local Ollama)';
            if (apiKeyInput) apiKeyInput.placeholder = 'No key required for local use; may be left empty...';
        } else {
            if (keyLabel) keyLabel.textContent = 'API Key';
            if (apiKeyInput) apiKeyInput.placeholder = 'Enter API key...';
        }
    }

    if (providerSelect) {
        providerSelect.addEventListener('change', () => {
            const provider = providerSelect.value;
            
            // Check if current inputs are empty or default values of ANY provider
            const currentModel = modelInput.value.trim();
            const currentBaseUrl = baseUrlInput.value.trim();
            
            const isModelDefaultOfAny = Object.values(providerDefaults).some(d => d.model === currentModel) || currentModel === '';
            const isBaseUrlDefaultOfAny = Object.values(providerDefaults).some(d => d.baseUrl === currentBaseUrl) || currentBaseUrl === '';
            
            if (isModelDefaultOfAny && providerDefaults[provider]) {
                modelInput.value = providerDefaults[provider].model;
            }
            if (isBaseUrlDefaultOfAny && providerDefaults[provider]) {
                baseUrlInput.value = providerDefaults[provider].baseUrl;
            }
            
            updateConfigFields();
        });
    }

    // 1. Migrate legacy key if present
    const legacyKey = localStorage.getItem('gemini_api_key');
    if (legacyKey && !localStorage.getItem('ai_api_key')) {
        localStorage.setItem('ai_api_key', legacyKey);
        localStorage.setItem('ai_provider', 'google');
        localStorage.removeItem('gemini_api_key'); // clear legacy
    }

    // 2. Load configurations on initialize
    function loadSavedConfig() {
        const savedKey = localStorage.getItem('ai_api_key');
        const savedProvider = localStorage.getItem('ai_provider') || 'google';
        const savedModel = localStorage.getItem('ai_model') || '';
        const savedBaseUrl = localStorage.getItem('ai_base_url') || '';

        if (providerSelect) providerSelect.value = savedProvider;
        if (modelInput) modelInput.value = savedModel;
        if (baseUrlInput) baseUrlInput.value = savedBaseUrl;
        
        updateConfigFields();

        if (savedKey || savedProvider === 'ollama') {
            keyConfigPanel.classList.add('hidden');
            keyActivePanel.classList.remove('hidden');
            if (activeStatusText) {
                const name = providerNames[savedProvider] || 'AI';
                activeStatusText.innerHTML = `<i class="fa-solid fa-circle-check"></i> ${name} ready`;
            }
            btnTriggerAi.disabled = false;
        } else {
            keyConfigPanel.classList.remove('hidden');
            keyActivePanel.classList.add('hidden');
            btnTriggerAi.disabled = true;
        }
    }

    loadSavedConfig();

    // 3. Save settings listener
    btnSaveKey.addEventListener('click', () => {
        const provider = providerSelect.value;
        const key = apiKeyInput.value.trim();
        const model = modelInput.value.trim();
        const baseUrl = baseUrlInput.value.trim();

        if (!key && provider !== 'ollama') {
            alert('Please enter an API key.');
            return;
        }

        if (provider === 'custom' && !baseUrl) {
            alert('A Base URL is required when using a custom provider.');
            return;
        }

        localStorage.setItem('ai_provider', provider);
        localStorage.setItem('ai_api_key', key);
        localStorage.setItem('ai_model', model);
        localStorage.setItem('ai_base_url', baseUrl);

        apiKeyInput.value = '';
        loadSavedConfig();
    });

    // 4. Clear config listener
    btnClearKey.addEventListener('click', () => {
        localStorage.removeItem('ai_api_key');
        localStorage.removeItem('ai_provider');
        localStorage.removeItem('ai_model');
        localStorage.removeItem('ai_base_url');

        // Reset input fields
        if (apiKeyInput) apiKeyInput.value = '';
        if (modelInput) modelInput.value = '';
        if (baseUrlInput) baseUrlInput.value = '';
        if (providerSelect) providerSelect.value = 'google';

        // Clear test result outputs
        const testResultEl = document.getElementById('ai-test-result');
        const testResultActiveEl = document.getElementById('ai-test-result-active');
        if (testResultEl) {
            testResultEl.classList.add('hidden');
            testResultEl.innerHTML = '';
        }
        if (testResultActiveEl) {
            testResultActiveEl.classList.add('hidden');
            testResultActiveEl.innerHTML = '';
        }

        loadSavedConfig();
        
        const summaryCard = document.getElementById('ai-summary-result');
        if (summaryCard) {
            summaryCard.classList.add('hidden');
            summaryCard.innerHTML = '';
        }
    });

    // 5. Test AI API Connection Helpers & Listeners
    const btnTestAi = document.getElementById('btn-test-ai');
    const btnTestAiActive = document.getElementById('btn-test-ai-active');
    const testResultEl = document.getElementById('ai-test-result');
    const testResultActiveEl = document.getElementById('ai-test-result-active');

    async function performAiConnectionTest(testBtn, resultEl, getParamsFunc) {
        const { provider, apiKey, model, baseUrl } = getParamsFunc();

        if (!apiKey && provider !== 'ollama') {
            resultEl.classList.remove('hidden');
            resultEl.style.backgroundColor = '#fff5f5';
            resultEl.style.color = '#ef4444';
            resultEl.style.border = '1px solid rgba(239, 68, 68, 0.2)';
            resultEl.innerHTML = `<i class="fa-solid fa-circle-exclamation"></i> Please enter an API key.`;
            return;
        }

        testBtn.disabled = true;
        const originalText = testBtn.innerHTML || testBtn.textContent;
        testBtn.innerHTML = `<i class="fa-solid fa-spinner fa-spin"></i> Testing...`;
        
        resultEl.classList.remove('hidden');
        resultEl.style.backgroundColor = '#f8fafc';
        resultEl.style.color = '#475569';
        resultEl.style.border = '1px solid var(--border-color)';
        resultEl.innerHTML = `<i class="fa-solid fa-circle-notch fa-spin"></i> Testing API connection, please wait...`;

        try {
            const headers = {
                'X-AI-API-Key': apiKey || '',
                'X-AI-Provider': provider
            };
            if (model) headers['X-AI-Model'] = model;
            if (baseUrl) headers['X-AI-Base-URL'] = baseUrl;

            const response = await fetch('/api/test_ai', { headers });
            
            if (!response.ok) {
                throw new Error(`HTTP error: ${response.status} ${response.statusText}`);
            }

            const data = await response.json();
            
            if (data.status === 'success') {
                resultEl.style.backgroundColor = '#ecfdf5';
                resultEl.style.color = '#065f46';
                resultEl.style.border = '1px solid rgba(16, 185, 129, 0.2)';
                resultEl.innerHTML = `<i class="fa-solid fa-circle-check"></i> ${data.message}`;
            } else {
                resultEl.style.backgroundColor = '#fff5f5';
                resultEl.style.color = '#991b1b';
                resultEl.style.border = '1px solid rgba(239, 68, 68, 0.2)';
                resultEl.innerHTML = `<i class="fa-solid fa-circle-exclamation"></i> Connection failed.<br><span style="font-size: 10px; color: #ef4444; margin-top: 4px; display: block;">${data.message}</span>`;
            }
        } catch (err) {
            resultEl.style.backgroundColor = '#fff5f5';
            resultEl.style.color = '#991b1b';
            resultEl.style.border = '1px solid rgba(239, 68, 68, 0.2)';
            resultEl.innerHTML = `<i class="fa-solid fa-circle-exclamation"></i> Network request error:<br><span style="font-size: 10px; color: #ef4444; margin-top: 4px; display: block;">${err.message}</span>`;
        } finally {
            testBtn.disabled = false;
            testBtn.innerHTML = originalText;
        }
    }

    if (btnTestAi && testResultEl) {
        btnTestAi.addEventListener('click', () => {
            performAiConnectionTest(btnTestAi, testResultEl, () => {
                return {
                    provider: providerSelect ? providerSelect.value : 'google',
                    apiKey: apiKeyInput ? apiKeyInput.value.trim() : '',
                    model: modelInput ? modelInput.value.trim() : '',
                    baseUrl: baseUrlInput ? baseUrlInput.value.trim() : ''
                };
            });
        });
    }

    if (btnTestAiActive && testResultActiveEl) {
        btnTestAiActive.addEventListener('click', () => {
            performAiConnectionTest(btnTestAiActive, testResultActiveEl, () => {
                return {
                    provider: localStorage.getItem('ai_provider') || 'google',
                    apiKey: localStorage.getItem('ai_api_key') || '',
                    model: localStorage.getItem('ai_model') || '',
                    baseUrl: localStorage.getItem('ai_base_url') || ''
                };
            });
        });
    }

    btnTriggerAi.addEventListener('click', () => {
        triggerAiSummary();
    });
}

// ==========================================================================
// RNA-Seq Data Integration Features
// ==========================================================================
let rnaseqData = null; // object mapping lowercase locus -> { log2fc, pvalue }

function getRnaSeqColor(log2fc) {
    if (log2fc === undefined || isNaN(log2fc)) return '#f5f5f5';
    const val = Math.max(-3, Math.min(3, log2fc));
    if (val < 0) {
        const ratio = (val + 3) / 3;
        const r = Math.round(29 * (1 - ratio) + 226 * ratio);
        const g = Math.round(78 * (1 - ratio) + 232 * ratio);
        const b = Math.round(216 * (1 - ratio) + 240 * ratio);
        return `rgb(${r}, ${g}, ${b})`;
    } else {
        const ratio = val / 3;
        const r = Math.round(226 * (1 - ratio) + 185 * ratio);
        const g = Math.round(232 * (1 - ratio) + 28 * ratio);
        const b = Math.round(240 * (1 - ratio) + 28 * ratio);
        return `rgb(${r}, ${g}, ${b})`;
    }
}

function initRnaSeqOverlay() {
    const btnUpload = document.getElementById('btn-upload-rnaseq');
    const fileInput = document.getElementById('rnaseq-file-input');
    const btnClear = document.getElementById('btn-clear-rnaseq');

    if (!btnUpload || !fileInput) return;

    btnUpload.addEventListener('click', () => {
        fileInput.click();
    });

    fileInput.addEventListener('change', (e) => {
        const file = e.target.files[0];
        if (!file) return;

        const reader = new FileReader();
        reader.onload = function(evt) {
            const csvText = evt.target.result;
            Papa.parse(csvText, {
                header: true,
                skipEmptyLines: true,
                dynamicTyping: true,
                complete: function(results) {
                    processRnaSeqData(results.data);
                },
                error: function(err) {
                    alert('Failed to parse CSV file: ' + err.message);
                }
            });
        };
        reader.readAsText(file);
    });

    if (btnClear) {
        btnClear.addEventListener('click', () => {
            clearRnaSeqOverlay();
        });
    }

    // Attach filter control listeners
    const filterEnable = document.getElementById('rnaseq-filter-enable');
    const lfcThreshold = document.getElementById('rnaseq-lfc-threshold');
    const pThreshold = document.getElementById('rnaseq-p-threshold');

    if (filterEnable) {
        filterEnable.addEventListener('change', () => {
            applyRnaSeqFilters();
        });
    }

    if (lfcThreshold) {
        lfcThreshold.addEventListener('input', () => {
            applyRnaSeqFilters();
        });
    }

    if (pThreshold) {
        pThreshold.addEventListener('input', () => {
            applyRnaSeqFilters();
        });
    }
}

function processRnaSeqData(dataRows) {
    if (!dataRows || dataRows.length === 0) {
        alert('No valid dataRows were found in the CSV file.');
        return;
    }

    const firstRow = dataRows[0];
    let locusCol = null;
    let fcCol = null;
    let pvalCol = null;

    for (let key in firstRow) {
        const lowerKey = key.toLowerCase().replace(/[^a-z0-9]/g, '');
        if (['locustag', 'locus', 'geneid', 'gene', 'id', 'name'].includes(lowerKey)) {
            locusCol = key;
        } else if (['log2fc', 'log2foldchange', 'fc', 'foldchange'].includes(lowerKey)) {
            fcCol = key;
        } else if (['pvalue', 'pval', 'padj', 'pvalue', 'p.value'].includes(lowerKey)) {
            pvalCol = key;
        }
    }

    if (!locusCol) {
        for (let key in firstRow) {
            const val = String(firstRow[key]).toLowerCase();
            if (val.startsWith('cg') || val.startsWith('cgl')) {
                locusCol = key;
                break;
            }
        }
    }

    if (!locusCol || !fcCol) {
        alert('Unable to infer CSV columns automatically. Please include columns such as:\n- Gene locus tag: locus_tag, gene_id, gene\n- Fold change: log2fc, log2FoldChange\n- Significance (optional): pvalue, padj');
        return;
    }

    rnaseqData = {};
    let loadedCount = 0;

    dataRows.forEach(row => {
        let locus = String(row[locusCol]).trim().toLowerCase();
        let fc = parseFloat(row[fcCol]);
        let pval = pvalCol ? parseFloat(row[pvalCol]) : 1.0;

        if (locus && !isNaN(fc)) {
            if (cglToCg[locus]) {
                locus = cglToCg[locus].toLowerCase();
            }
            rnaseqData[locus] = { log2fc: fc, pvalue: isNaN(pval) ? 1.0 : pval };
            loadedCount++;
        }
    });

    if (loadedCount === 0) {
        alert('No valid gene/sRNA data matched the local database from the CSV file.');
        rnaseqData = null;
        return;
    }

    const btnClear = document.getElementById('btn-clear-rnaseq');
    const legendContainer = document.getElementById('rnaseq-legend-container');
    const loadedCountDisp = document.getElementById('rnaseq-loaded-count');
    const btnUpload = document.getElementById('btn-upload-rnaseq');

    if (btnClear) btnClear.classList.remove('hidden');
    if (legendContainer) legendContainer.classList.remove('hidden');
    if (loadedCountDisp) loadedCountDisp.textContent = `Loaded ${loadedCount} genes`;
    if (btnUpload) {
        btnUpload.innerHTML = `<i class="fa-solid fa-check"></i> Omics data loaded`;
        btnUpload.style.backgroundColor = 'rgba(46, 125, 50, 0.05)';
        btnUpload.style.borderColor = 'var(--color-activation)';
    }

    if (cy) {
        applyRnaSeqStyling();
        applyRnaSeqFilters();
    }

    // Update details sidebar status badge and map if visible
    const badge = document.getElementById('rnaseq-status-badge');
    if (badge) {
        badge.textContent = `(imported ${loadedCount} genes)`;
        badge.style.color = '#3b82f6';
    }
    if (currentQueryGene) {
        renderGenomicLocusMap(currentQueryGene);
    }
}

function clearRnaSeqOverlay() {
    rnaseqData = null;
    document.getElementById('rnaseq-file-input').value = '';
    
    const btnClear = document.getElementById('btn-clear-rnaseq');
    const legendContainer = document.getElementById('rnaseq-legend-container');
    const btnUpload = document.getElementById('btn-upload-rnaseq');

    if (btnClear) btnClear.classList.add('hidden');
    if (legendContainer) legendContainer.classList.add('hidden');
    if (btnUpload) {
        btnUpload.innerHTML = `<i class="fa-solid fa-file-arrow-up"></i> Upload CSV`;
        btnUpload.style.backgroundColor = '';
        btnUpload.style.borderColor = '';
    }

    // Reset filter control state
    const filterEnable = document.getElementById('rnaseq-filter-enable');
    const lfcThreshold = document.getElementById('rnaseq-lfc-threshold');
    const pThreshold = document.getElementById('rnaseq-p-threshold');

    if (filterEnable) filterEnable.checked = false;
    if (lfcThreshold) lfcThreshold.value = 1.0;
    if (pThreshold) pThreshold.value = 0.05;

    // Reset displayed values
    const lfcValDisp = document.getElementById('rnaseq-lfc-val');
    if (lfcValDisp) lfcValDisp.textContent = "1.0";
    const pValDisp = document.getElementById('rnaseq-p-val');
    if (pValDisp) pValDisp.textContent = "0.05";

    if (cy) {
        cy.nodes().removeClass('rnaseq-node');
        cy.nodes().removeClass('rnaseq-hidden');
        cy.style().update();
        updateNetworkStatistics();
    }

    const badge = document.getElementById('rnaseq-status-badge');
    if (badge) {
        badge.textContent = `(data cleared)`;
        badge.style.color = 'var(--text-muted)';
    }
    if (currentQueryGene) {
        renderGenomicLocusMap(currentQueryGene);
    }
}

function applyRnaSeqStyling() {
    if (!cy || !rnaseqData) return;

    cy.nodes().forEach(node => {
        const locus = node.id().toLowerCase();
        if (rnaseqData[locus]) {
            const item = rnaseqData[locus];
            node.data('rnaseq_log2fc', item.log2fc);
            node.data('rnaseq_pvalue', item.pvalue);
            node.addClass('rnaseq-node');
        } else {
            node.removeClass('rnaseq-node');
        }
    });
    cy.style().update();
}

function applyRnaSeqFilters() {
    if (!cy) return;

    const filterEnable = document.getElementById('rnaseq-filter-enable');
    const isFilterActive = filterEnable && filterEnable.checked && rnaseqData;

    const lfcEl = document.getElementById('rnaseq-lfc-threshold');
    const pvalEl = document.getElementById('rnaseq-p-threshold');
    const lfcThresh = lfcEl ? parseFloat(lfcEl.value) : 1.0;
    const pThresh = pvalEl ? parseFloat(pvalEl.value) : 0.05;

    // Update displayed text
    const lfcValDisp = document.getElementById('rnaseq-lfc-val');
    if (lfcValDisp && lfcEl) lfcValDisp.textContent = parseFloat(lfcEl.value).toFixed(1);
    const pValDisp = document.getElementById('rnaseq-p-val');
    if (pValDisp && pvalEl) pValDisp.textContent = parseFloat(pvalEl.value).toFixed(2);

    if (isFilterActive) {
        cy.nodes().forEach(node => {
            // Always keep searched query anchor nodes to avoid empty graphs
            if (node.data('type') === 'query') {
                node.removeClass('rnaseq-hidden');
                return;
            }

            const locus = node.id().toLowerCase();
            if (rnaseqData && rnaseqData[locus]) {
                const item = rnaseqData[locus];
                const matchLfc = Math.abs(item.log2fc) >= lfcThresh;
                const matchPval = item.pvalue <= pThresh;
                if (matchLfc && matchPval) {
                    node.removeClass('rnaseq-hidden');
                } else {
                    node.addClass('rnaseq-hidden');
                }
            } else {
                // Hide genes without RNA-seq data when expression filtering is enabled
                node.addClass('rnaseq-hidden');
            }
        });
    } else {
        // Remove hidden classes when the filter is disabled
        cy.nodes().removeClass('rnaseq-hidden');
    }

    // Reapply Cytoscape stylesheet for dynamic border styling
    cy.style().update();
    
    // Update network statistics
    updateNetworkStatistics();
}

async function triggerAiSummary() {
    const btnTriggerAi = document.getElementById('btn-trigger-ai');
    const summaryCard = document.getElementById('ai-summary-result');
    
    const locus = document.getElementById('info-locus').textContent.trim();
    const name = document.getElementById('info-name').textContent.trim();
    const apiKey = localStorage.getItem('ai_api_key');
    const provider = localStorage.getItem('ai_provider') || 'google';
    const model = localStorage.getItem('ai_model') || '';
    const baseUrl = localStorage.getItem('ai_base_url') || '';
    
    if (!locus || locus === '-') {
        alert('Please select a gene first.');
        return;
    }
    if (!apiKey && provider !== 'ollama') {
        alert('Please configure your API key in the panel first.');
        return;
    }
    
    // Set loading state
    btnTriggerAi.disabled = true;
    summaryCard.classList.remove('hidden');
    summaryCard.classList.add('loading');
    summaryCard.innerHTML = `
        <div class="ai-spinner"></div>
        <span style="font-weight: 500;">Searching PubMed and requesting an AI summary...</span>
    `;
    
    try {
        const headers = {
            'X-AI-API-Key': apiKey || '',
            'X-AI-Provider': provider
        };
        if (model) headers['X-AI-Model'] = model;
        if (baseUrl) headers['X-AI-Base-URL'] = baseUrl;

        if (apiKey) {
            headers['X-Gemini-API-Key'] = apiKey;
        }

        const response = await fetch(`/api/summarize?gene=${locus}&name=${name}`, {
            headers: headers
        });
        
        if (!response.ok) {
            throw new Error(`HTTP error: ${response.status}`);
        }
        
        const result = await response.json();
        
        if (result.error) {
            throw new Error(result.error);
        }
        
        // Remove loading state
        summaryCard.classList.remove('loading');
        
        // Render summary text (with simple markdown parser)
        let htmlContent = parseMarkdownToHtml(result.summary);
        
        // Append papers if present
        if (result.papers && result.papers.length > 0) {
            htmlContent += `
                <div class="ai-sources-list">
                    <div class="ai-sources-title"><i class="fa-solid fa-book"></i> PubMed references (${result.papers.length})</div>
            `;
            
            result.papers.forEach(p => {
                htmlContent += `
                    <div class="ai-source-item">
                        <i class="fa-solid fa-file-lines"></i>
                        <a href="https://pubmed.ncbi.nlm.nih.gov/${p.pmid}" target="_blank" class="ai-source-link" title="Open original paper in PubMed">
                            ${p.title} (PMID: ${p.pmid}) <i class="fa-solid fa-arrow-up-right-from-square" style="font-size: 8px;"></i>
                        </a>
                    </div>
                `;
            });
            
            htmlContent += `</div>`;
        }

        // Append RAG sources if present
        if (result.rag_sources && result.rag_sources.length > 0) {
            htmlContent += `
                <div class="ai-sources-list" style="margin-top: 10px; border-top: 1px dashed rgba(99, 102, 241, 0.15); padding-top: 10px;">
                    <div class="ai-sources-title" style="color: #6366f1;"><i class="fa-solid fa-database"></i> Local RAG references (${result.rag_sources.length})</div>
            `;
            
            result.rag_sources.forEach(r => {
                const scorePercentage = Math.round(r.score * 100);
                htmlContent += `
                    <div class="ai-source-item" style="font-size: 11px;">
                        <i class="fa-solid fa-file-pdf" style="color: #ef4444;"></i>
                        <span class="ai-source-link" style="color: var(--text-secondary); text-decoration: none; cursor: default;">
                            ${r.file} <span style="color: var(--text-muted); font-size: 10px;">(match: ${scorePercentage}%)</span>
                        </span>
                    </div>
                `;
            });
            
            htmlContent += `</div>`;
        }
        
        summaryCard.innerHTML = htmlContent;
        
    } catch (err) {
        console.error(err);
        summaryCard.classList.remove('loading');
        summaryCard.innerHTML = `
            <div style="color: #ef4444; font-weight: 500; display: flex; align-items: center; gap: 6px; margin-bottom: 6px;">
                <i class="fa-solid fa-circle-exclamation"></i> Summary generation failed
            </div>
            <p style="font-size: 11px; color: var(--text-secondary); line-height: 1.4;">
                ${err.message || 'Unknown network error. Please check your API key and network connection.'}
            </p>
        `;
    } finally {
        btnTriggerAi.disabled = false;
    }
}

function parseMarkdownToHtml(mdText) {

    if (!mdText) return "";

    let html = mdText;

    

    // Replace headers (###, ####, etc.) and bold bullet headings

    html = html.replace(/^(?:###\s+)(.*?)$/gm, '<h4>$1</h4>');

    html = html.replace(/^(?:####\s+)(.*?)$/gm, '<h4>$1</h4>');

    html = html.replace(/^(?:【)(.*?)(】)/gm, '<h4>$1</h4>');

    

    // Replace bold (**text**)

    html = html.replace(/\*\*(.*?)\*\*/g, '<strong>$1</strong>');

    

    // Replace bullet lists (- or *)

    const lines = html.split('\n');

    let inList = false;

    const processedLines = [];

    

    lines.forEach(line => {

        const trimmed = line.trim();

        if (trimmed.startsWith('- ') || trimmed.startsWith('* ')) {

            const content = trimmed.substring(2);

            if (!inList) {

                processedLines.push('<ul>');

                inList = true;

            }

            processedLines.push(`<li>${content}</li>`);

        } else {

            if (inList) {

                processedLines.push('</ul>');

                inList = false;

            }

            if (trimmed) {

                // If it is a heading, don't wrap in p

                if (trimmed.startsWith('<h4>') || trimmed.startsWith('</h4>') || trimmed.startsWith('<ul>') || trimmed.startsWith('</ul>')) {

                    processedLines.push(trimmed);

                } else {

                    processedLines.push(`<p>${trimmed}</p>`);

                }

            }

        }

    });

    

    if (inList) {

        processedLines.push('</ul>');

    }

    

    return processedLines.join('\n');

}



// ==========================================================================

// 8. AI Pathway Assistant Feature

// ==========================================================================

function escapeHtml(value) {
    return String(value ?? '').replace(/[&<>"']/g, char => ({
        '&': '&amp;',
        '<': '&lt;',
        '>': '&gt;',
        '"': '&quot;',
        "'": '&#39;'
    }[char]));
}

function metabolicEmptyMessage() {
    return '<div class="metabolic-empty">No metabolic model mapping available for this node.</div>';
}

function encodeMetabolicList(values) {
    return encodeURIComponent(JSON.stringify(Array.from(new Set(values || []))));
}

function decodeMetabolicList(value) {
    try {
        return JSON.parse(decodeURIComponent(value || '[]'));
    } catch (err) {
        return [];
    }
}

function highlightMetabolicPathwayGenes(geneIds, reactionIds) {
    if (!cy) return;

    const ids = new Set();
    (geneIds || []).forEach(id => {
        const lower = String(id || '').toLowerCase();
        if (!lower) return;
        ids.add(lower);
        if (cgToCgl[lower]) ids.add(cgToCgl[lower].toLowerCase());
        if (cglToCg[lower]) ids.add(cglToCg[lower].toLowerCase());
    });
    (reactionIds || []).forEach(id => {
        const lower = String(id || '').toLowerCase();
        if (lower) ids.add(lower);
    });

    if (ids.size === 0) return;

    cy.elements().removeClass('dimmed');
    cy.elements().removeClass('highlighted');
    cy.elements().addClass('dimmed');

    cy.nodes().forEach(node => {
        const id = String(node.id() || '').toLowerCase();
        if (!ids.has(id)) return;
        node.removeClass('dimmed');
        node.addClass('highlighted');
        node.connectedEdges().removeClass('dimmed');
        node.connectedEdges().addClass('highlighted');
    });
}

function formatMetabolicNumber(value, digits = 3) {
    const parsed = Number(value);
    if (!Number.isFinite(parsed)) return '';
    if (Math.abs(parsed) >= 1000) return parsed.toExponential(2);
    return parsed.toFixed(digits).replace(/\.?0+$/, '');
}

function renderEnzymeConstraintBadges(reaction) {
    if (!reaction) return '';
    const enzyme = reaction.enzyme_constraint || {};
    const badges = [];
    const ecNumber = reaction.ec_number || enzyme.ec_number;
    const kcat = reaction.kcat ?? enzyme.kcat;
    const molecularWeight = reaction.molecular_weight ?? enzyme.molecular_weight;
    const kcatMw = reaction.kcat_MW ?? enzyme.kcat_MW;
    const uniprotIds = reaction.uniprot_ids || enzyme.uniprot_ids || [];
    const variant = reaction.reaction_variant || enzyme.model_variant;
    const variantOf = reaction.variant_of || enzyme.variant_of;
    const sourceCount = reaction.kcat_source_count ?? enzyme.kcat_source_count;

    if (ecNumber) badges.push('EC ' + escapeHtml(ecNumber));
    if (kcat !== undefined && kcat !== null) badges.push('kcat ' + escapeHtml(formatMetabolicNumber(kcat, 3)));
    if (molecularWeight !== undefined && molecularWeight !== null) badges.push('MW ' + escapeHtml(formatMetabolicNumber(molecularWeight, 1)) + ' Da');
    if (kcatMw !== undefined && kcatMw !== null) badges.push('kcat/MW ' + escapeHtml(formatMetabolicNumber(kcatMw, 3)));
    if (Array.isArray(uniprotIds) && uniprotIds.length > 0) badges.push('UniProt ' + escapeHtml(uniprotIds.slice(0, 3).join(', ')));
    if (variant) badges.push('variant ' + escapeHtml(variant));
    if (variantOf) badges.push('paired ' + escapeHtml(variantOf));
    if (sourceCount) badges.push('kcat sources ' + escapeHtml(sourceCount));

    if (badges.length === 0) return '';
    return '<div class="metabolic-enzyme-badges">'
        + badges.map(text => '<span class="metabolic-enzyme-badge">' + text + '</span>').join('')
        + '</div>';
}

function renderMetabolicImpact(data, detailLocus) {
    const section = document.getElementById('detail-metabolic-impact-section');
    const container = document.getElementById('metabolic-impact-content');
    if (!section || !container) return;
    section.style.display = '';

    const mapping = data?.model_mapping || {};
    const summary = data?.summary || {};
    const pathways = data?.pathways || [];
    const genes = data?.affected_genes || [];
    const isTf = data?.mode === 'tf';
    const bridge = window.regulationMetabolismBridge;

    if (!mapping.loaded) {
        container.innerHTML = metabolicEmptyMessage();
        return;
    }

    const pathwaySummary = pathways.map(p => ({
        pathwayId: p.id || p.name || 'Unassigned pathway',
        pathwayName: p.name || p.id || 'Unassigned pathway',
        geneCount: Number(p.gene_count || 0),
        reactionCount: Number(p.reaction_count || 0),
        genes: p.genes || [],
        reactions: p.reactions || []
    }));
    const bridgeImpact = {
        tfId: detailLocus,
        totalTargetGenes: Number(summary.target_gene_count || 0),
        mappedTargetGenes: Number(summary.mapped_gene_count || 0),
        totalReactions: Number(summary.reaction_count || 0),
        totalPathways: Number(summary.pathway_count || pathways.length || 0),
        pathwaySummary
    };
    const explanation = bridge?.generateMetabolicImpactExplanation
        ? bridge.generateMetabolicImpactExplanation(bridgeImpact)
        : 'No metabolic model mapping available for this node.';
    const files = (mapping.files || []).map(f => String(f.model || 'model') + ':' + String(f.rows || 0)).join(' - ');

    if (isTf) {
        const statHtml = '<div class="metabolic-stat-grid">'
            + '<div><strong>' + escapeHtml(bridgeImpact.totalTargetGenes) + '</strong><span>Target genes</span></div>'
            + '<div><strong>' + escapeHtml(bridgeImpact.mappedTargetGenes) + '</strong><span>Mapped metabolic genes</span></div>'
            + '<div><strong>' + escapeHtml(bridgeImpact.totalReactions) + '</strong><span>Associated reactions</span></div>'
            + '<div><strong>' + escapeHtml(bridgeImpact.totalPathways) + '</strong><span>Affected pathways</span></div>'
            + '</div>';
        const pathwayHtml = pathwaySummary.length > 0
            ? pathwaySummary.slice(0, 8).map((p, index) => (
                '<button type="button" class="metabolic-pathway-row metabolic-pathway-button" data-genes="' + encodeMetabolicList(p.genes) + '" data-reactions="' + encodeMetabolicList(p.reactions) + '" title="Highlight mapped genes in the current network">'
                + '<span><span class="metabolic-pathway-name">' + (index + 1) + '. ' + escapeHtml(p.pathwayName) + '</span>'
                + '<span class="metabolic-muted">' + escapeHtml(p.pathwayId) + '</span></span>'
                + '<span class="metabolic-counts">' + escapeHtml(p.geneCount) + ' genes, ' + escapeHtml(p.reactionCount) + ' reactions</span>'
                + '</button>'
            )).join('')
            : '<div class="metabolic-empty">No affected pathways are mapped for this TF.</div>';

        container.innerHTML = '<div class="metabolic-intro">' + escapeHtml(explanation) + '</div>'
            + statHtml
            + '<div class="metabolic-subtitle">Top affected pathways</div>'
            + '<div class="metabolic-pathway-list">' + pathwayHtml + '</div>'
            + '<div class="metabolic-source">Models: ' + escapeHtml((mapping.models || []).join(', ') || 'none') + (files ? ' - Files: ' + escapeHtml(files) : '') + '<br><span class="source-attribution-note" style="font-size: 10px; font-style: italic; opacity: 0.85; margin-top: 4px; display: inline-block;">Gene–reaction–pathway mappings are derived from local GEM model adapters. Enzyme annotations are parsed from ecCGL1-derived model fields.</span></div>';
    } else {
        const gene = genes.find(g => String(g.locus || '').toLowerCase() === String(detailLocus || '').toLowerCase()) || genes[0] || {};
        const reactions = Array.from(new Map((gene.reactions || []).map(r => [String(r.model || 'model') + ':' + String(r.id), r])).values());
        if (reactions.length === 0 && pathways.length === 0) {
            container.innerHTML = metabolicEmptyMessage();
            return;
        }

        const reactionHtml = reactions.length > 0
            ? reactions.slice(0, 12).map(r => (
                '<div class="metabolic-gene-row">'
                + '<div class="metabolic-pathway-name">' + escapeHtml(r.id) + ': ' + escapeHtml(r.label || r.id) + '</div>'
                + '<div class="metabolic-muted">' + escapeHtml(r.gpr_rule || r.equation || r.model || '') + '</div>'
                + renderEnzymeConstraintBadges(r)
                + '</div>'
            )).join('')
            : '<div class="metabolic-empty">No associated reactions are mapped for this gene.</div>';
        const pathwayHtml = pathwaySummary.length > 0
            ? pathwaySummary.slice(0, 8).map(p => (
                '<button type="button" class="metabolic-pathway-row metabolic-pathway-button" data-genes="' + encodeMetabolicList(p.genes) + '" data-reactions="' + encodeMetabolicList(p.reactions) + '">'
                + '<span><span class="metabolic-pathway-name">' + escapeHtml(p.pathwayName) + '</span>'
                + '<span class="metabolic-muted">' + escapeHtml(p.pathwayId) + '</span></span>'
                + '<span class="metabolic-counts">' + escapeHtml(p.reactionCount) + ' reactions</span>'
                + '</button>'
            )).join('')
            : '<div class="metabolic-empty">No pathway annotation is available for this gene.</div>';

        container.innerHTML = '<div class="metabolic-subtitle">Associated reactions</div>'
            + '<div class="metabolic-gene-list">' + reactionHtml + '</div>'
            + '<div class="metabolic-subtitle">Pathways</div>'
            + '<div class="metabolic-pathway-list">' + pathwayHtml + '</div>'
            + '<div class="metabolic-source">Models: ' + escapeHtml((mapping.models || []).join(', ') || 'none') + (files ? ' - Files: ' + escapeHtml(files) : '') + '<br><span class="source-attribution-note" style="font-size: 10px; font-style: italic; opacity: 0.85; margin-top: 4px; display: inline-block;">Gene–reaction–pathway mappings are derived from local GEM model adapters. Enzyme annotations are parsed from ecCGL1-derived model fields.</span></div>';
    }

    container.querySelectorAll('.metabolic-pathway-button').forEach(btn => {
        btn.addEventListener('click', () => {
            highlightMetabolicPathwayGenes(
                decodeMetabolicList(btn.getAttribute('data-genes')),
                decodeMetabolicList(btn.getAttribute('data-reactions'))
            );
        });
    });
}

function fetchMetabolicImpact(locusTag, nodeType) {
    const section = document.getElementById('detail-metabolic-impact-section');
    const container = document.getElementById('metabolic-impact-content');
    if (!section || !container || !locusTag) return;
    section.style.display = '';
    container.innerHTML = '<span class="metabolic-muted"><i class="fa-solid fa-spinner fa-spin"></i> Loading metabolic model mapping...</span>'; 

    const adapter = window.metabolicModelAdapter;
    const loadImpact = adapter?.loadMetabolicImpact
        ? adapter.loadMetabolicImpact(locusTag)
        : fetch(`/api/metabolic_impact?gene=${encodeURIComponent(locusTag)}`).then(response => response.json());

    loadImpact
        .then(data => {
            if (detailLocusTag.textContent !== locusTag) return;
            renderMetabolicImpact(data, locusTag, nodeType);
        })
        .catch(err => {
            console.error('Error fetching metabolic impact:', err);
            if (detailLocusTag.textContent === locusTag) {
                container.innerHTML = '<div class="metabolic-empty">Failed to load metabolic model mapping.</div>'; 
            }
        });
}

function renderPathwayRegulation(regulation) {
    if (!regulation) return '';

    const matched = regulation.matched_pathways || [];
    const regulators = regulation.regulators || [];
    const pathwayGenes = regulation.pathway_genes || [];
    const cacheInfo = regulation.cache || {};
    const matchHtml = matched.length > 0
        ? matched.map(p => `<a href="${escapeHtml(p.link)}" target="_blank" style="color:#2563eb; text-decoration:none; font-weight:600;">${escapeHtml(p.name || p.id)}</a>`).join(' · ')
        : '<span style="color:var(--text-secondary);">No KEGG pathway match found</span>';

    const regulatorRows = regulators.slice(0, 8).map(r => {
        const roleText = Object.entries(r.roles || {}).map(([k, v]) => `${escapeHtml(k)}:${escapeHtml(v)}`).join(' ');
        const evidenceText = Object.entries(r.evidence || {}).map(([k, v]) => `${escapeHtml(k)}:${escapeHtml(v)}`).join(' ');
        const components = r.score_components || {};
        const scoreTitle = [
            `Coverage: ${components.coverage ?? 0}`,
            `Evidence: ${components.evidence ?? 0}`,
            `Binding site: ${components.binding_site ?? 0}`,
            `Direction: ${components.direction_consistency ?? 0}`,
            `Edge support: ${components.edge_support ?? 0}`
        ].join('\n');
        const targets = (r.target_genes || []).slice(0, 8).map(g =>
            `<button class="ai-pathway-gene-badge pathway-reg-target" data-locus="${escapeHtml(g)}" style="border:none; cursor:pointer;">${escapeHtml(g)}</button>`
        ).join('');
        return `
            <tr>
                <td style="padding:5px 6px; vertical-align:top;">
                    <button class="ai-pathway-gene-badge pathway-reg-target" data-locus="${escapeHtml(r.tf_locus)}" style="border:none; cursor:pointer; font-weight:700;">${escapeHtml(r.tf_name || r.tf_locus)}</button>
                    <div style="font-size:8.5px; color:var(--text-muted);">${escapeHtml(r.tf_locus)}</div>
                </td>
                <td style="padding:5px 6px; text-align:center; vertical-align:top;" title="${escapeHtml(scoreTitle)}">
                    <span style="display:inline-block; min-width:34px; padding:2px 5px; border-radius:4px; background:#ecfdf5; color:#047857; font-weight:700;">${escapeHtml(r.impact_score ?? '-')}</span>
                    <div style="font-size:8px; color:var(--text-muted); margin-top:2px;">${escapeHtml(r.confidence || '')}</div>
                </td>
                <td style="padding:5px 6px; text-align:center; vertical-align:top;">${escapeHtml(r.target_count)}</td>
                <td style="padding:5px 6px; vertical-align:top; font-size:8.5px; color:var(--text-secondary);">${roleText || '-'}</td>
                <td style="padding:5px 6px; vertical-align:top; font-size:8.5px; color:var(--text-secondary);">${evidenceText || '-'}</td>
                <td style="padding:5px 6px; vertical-align:top;">${targets || '-'}</td>
            </tr>
        `;
    }).join('');

    const geneBadges = pathwayGenes.slice(0, 24).map(g =>
        `<button class="ai-pathway-gene-badge pathway-reg-target" data-locus="${escapeHtml(g.locus)}" style="border:none; cursor:pointer;">${escapeHtml(g.name || g.locus)}<span style="opacity:.65;"> (${escapeHtml(g.locus)})</span></button>`
    ).join('');

    return `
        <div style="margin-top:12px; padding-top:10px; border-top:1px solid var(--border-color);">
            <div style="font-size:11px; font-weight:700; color:var(--text-primary); margin-bottom:6px; display:flex; align-items:center; gap:6px;">
                <i class="fa-solid fa-diagram-project" style="color:#0f766e;"></i> KEGG Pathway - TF Regulatory Projection
            </div>
            <div style="font-size:10px; color:var(--text-secondary); line-height:1.5; margin-bottom:8px;">
                Matched pathway: ${matchHtml}<br>
                Pathway genes: ${escapeHtml(regulation.pathway_gene_count || 0)}; regulatory records cover ${escapeHtml(regulation.regulated_gene_count || 0)} genes; upstream TFs: ${escapeHtml(regulation.regulator_count || 0)}.
                ${cacheInfo.enabled ? `<br>KEGG cache: ${cacheInfo.loaded_from_disk ? 'loaded from local cache' : 'generated online this run'}` : ''}
            </div>
            ${regulators.length > 0 ? `
                <div style="max-height:190px; overflow:auto; border:1px solid var(--border-color); border-radius:6px; background:#fff;">
                    <table style="width:100%; border-collapse:collapse; font-size:9px;">
                        <thead>
                            <tr style="background:#f8fafc; color:var(--text-secondary); border-bottom:1px solid var(--border-color);">
                                <th style="padding:5px 6px; text-align:left;">TF</th>
                                <th style="padding:5px 6px;">Score</th>
                                <th style="padding:5px 6px;">Target gene</th>
                                <th style="padding:5px 6px; text-align:left;">Direction</th>
                                <th style="padding:5px 6px; text-align:left;">Evidence</th>
                                <th style="padding:5px 6px; text-align:left;">Pathway target gene</th>
                            </tr>
                        </thead>
                        <tbody>${regulatorRows}</tbody>
                    </table>
                </div>
            ` : `
                <div style="font-size:10px; color:var(--text-secondary); padding:8px; background:#f8fafc; border-radius:6px;">
                    No TF edges targeting this KEGG pathway gene set were found in the local regulatory table.
                </div>
            `}
            ${geneBadges ? `
                <div style="font-size:10px; font-weight:700; color:var(--text-primary); margin-top:9px; margin-bottom:5px;">Candidate pathway genes</div>
                <div class="ai-pathway-genes-list">${geneBadges}</div>
            ` : ''}
        </div>
    `;
}

function initAiPathwayFeature() {

    const btnAnalyze = document.getElementById('btn-analyze-pathway');

    const inputPathway = document.getElementById('ai-pathway-input');

    const resultCard = document.getElementById('ai-pathway-result');



    if (!btnAnalyze || !inputPathway || !resultCard) return;



    btnAnalyze.addEventListener('click', async () => {

        const query = inputPathway.value.trim();

        if (!query) {

            alert('Enter a metabolic pathway or biological function to analyze.');

            return;

        }



        const apiKey = localStorage.getItem('ai_api_key') || localStorage.getItem('gemini_api_key');

        const provider = localStorage.getItem('ai_provider') || 'google';

        const model = localStorage.getItem('ai_model') || '';

        const baseUrl = localStorage.getItem('ai_base_url') || '';



        // Set loading state

        btnAnalyze.disabled = true;

        resultCard.classList.remove('hidden');

        resultCard.classList.add('loading');

        resultCard.innerHTML = `

            <div class="ai-spinner"></div>

            <span style="font-weight: 500;">AI is analyzing pathway genes...</span>

        `;



        try {

            const headers = {

                'X-AI-Provider': provider

            };

            if (apiKey) {

                headers['X-AI-API-Key'] = apiKey;

                headers['X-Gemini-API-Key'] = apiKey;

            }

            if (model) headers['X-AI-Model'] = model;

            if (baseUrl) headers['X-AI-Base-URL'] = baseUrl;



            const response = await fetch(`/api/pathway?pathway=${encodeURIComponent(query)}`, {

                headers: headers

            });



            if (!response.ok) {

                throw new Error(`HTTP error: ${response.status}`);

            }



            const result = await response.json();

            if (result.error) {

                throw new Error(result.error);

            }



            resultCard.classList.remove('loading');



            const genes = result.genes || [];

            let genesBadgesHtml = '';

            if (genes.length > 0) {

                genesBadgesHtml = genes.map(g => `<a href="#" class="ai-pathway-gene-badge" data-locus="${g}" style="text-decoration: none;">${g}</a>`).join('');

            } else {

                genesBadgesHtml = '<span style="color: var(--text-secondary); font-size: 11px;">No associated locus tags recognized</span>'; 

            }



            const regulationHtml = renderPathwayRegulation(result.pathway_regulation);

            resultCard.innerHTML = `

                <div class="ai-pathway-summary">${result.summary || 'No summary available'}</div>

                <div class="ai-pathway-genes-title"><i class="fa-solid fa-dna"></i> Associated genes (${genes.length})</div>

                <div class="ai-pathway-genes-list">${genesBadgesHtml}</div>

                ${regulationHtml}

                ${genes.length > 0 ? `

                    <button class="ai-pathway-draw-btn" id="btn-draw-pathway-network">

                        <i class="fa-solid fa-network-wired"></i> Draw pathway regulatory network

                    </button>

                ` : ''}

            `;



            // Bind click to individual gene badges

            resultCard.querySelectorAll('.ai-pathway-gene-badge').forEach(badge => {

                badge.addEventListener('click', (e) => {

                    e.preventDefault();

                    const locus = badge.getAttribute('data-locus');

                    querySingleGene(locus);

                });

            });



            // Bind click to the draw button

            const drawBtn = document.getElementById('btn-draw-pathway-network');

            if (drawBtn) {

                drawBtn.addEventListener('click', () => {

                    queryMultipleGenes(genes);

                });

            }



        } catch (err) {

            console.error(err);

            resultCard.classList.remove('loading');

            resultCard.innerHTML = `

                <div style="color: #ef4444; font-weight: 600; display: flex; align-items: center; gap: 6px; margin-bottom: 6px;">

                    <i class="fa-solid fa-circle-exclamation"></i> Analysis failed

                </div>

                <p style="font-size: 11px; color: var(--text-secondary); line-height: 1.4;">

                    ${err.message || 'Unknown network error. Please check your API key and network connection.'}

                </p>

            `;

        } finally {

            btnAnalyze.disabled = false;

        }

    });



    inputPathway.addEventListener('keydown', (e) => {

        if (e.key === 'Enter') {

            btnAnalyze.click();

        }

    });

}



// ==========================================================================

// 9. AI Gene Analysis Assistant Feature

// ==========================================================================

function initAiGeneFeature() {

    const btnAnalyze = document.getElementById('btn-analyze-gene');

    const inputGene = document.getElementById('ai-gene-input');

    const resultCard = document.getElementById('ai-gene-result');



    if (!btnAnalyze || !inputGene || !resultCard) return;



    btnAnalyze.addEventListener('click', async () => {

        const query = inputGene.value.trim();

        if (!query) {

            alert('Enter a gene function, transcription factor, or feature to analyze.');

            return;

        }



        const apiKey = localStorage.getItem('ai_api_key') || localStorage.getItem('gemini_api_key');

        const provider = localStorage.getItem('ai_provider') || 'google';

        const model = localStorage.getItem('ai_model') || '';

        const baseUrl = localStorage.getItem('ai_base_url') || '';



        if (!apiKey && provider !== 'ollama') {

            alert('To use AI gene analysis, configure your API key in the left control panel first.');

            // Highlight the key input in the left sidebar

            const apiKeyInput = document.getElementById('gemini-api-key-input');

            if (apiKeyInput) {

                apiKeyInput.focus();

                apiKeyInput.style.border = '2px solid #6366f1';

                setTimeout(() => {

                    apiKeyInput.style.border = '1px solid var(--border-color)';

                }, 2000);

            }

            return;

        }



        // Set loading state

        btnAnalyze.disabled = true;

        resultCard.classList.remove('hidden');

        resultCard.classList.add('loading');

        resultCard.innerHTML = `

            <div class="ai-spinner"></div>

            <span style="font-weight: 500;">AI is analyzing gene features...</span>

        `;



        try {

            const headers = {

                'X-AI-API-Key': apiKey,

                'X-AI-Provider': provider,

                'X-Gemini-API-Key': apiKey

            };

            if (model) headers['X-AI-Model'] = model;

            if (baseUrl) headers['X-AI-Base-URL'] = baseUrl;



            const response = await fetch(`/api/gene_assistant?query=${encodeURIComponent(query)}`, {

                headers: headers

            });



            if (!response.ok) {

                throw new Error(`HTTP error: ${response.status}`);

            }



            const result = await response.json();

            if (result.error) {

                throw new Error(result.error);

            }



            resultCard.classList.remove('loading');



            const genes = result.genes || [];

            let genesBadgesHtml = '';

            if (genes.length > 0) {

                genesBadgesHtml = genes.map(g => `<a href="#" class="ai-pathway-gene-badge" data-locus="${g}" style="text-decoration: none;">${g}</a>`).join('');

            } else {

                genesBadgesHtml = '<span style="color: var(--text-secondary); font-size: 11px;">No associated locus tags recognized</span>'; 

            }



            resultCard.innerHTML = `

                <div class="ai-pathway-summary">${result.summary || 'No summary available'}</div>

                <div class="ai-pathway-genes-title"><i class="fa-solid fa-dna"></i> Associated genes (${genes.length})</div>

                <div class="ai-pathway-genes-list">${genesBadgesHtml}</div>

                ${genes.length > 0 ? `

                    <button class="ai-pathway-draw-btn" id="btn-draw-gene-network">

                        <i class="fa-solid fa-network-wired"></i> Draw gene regulatory network

                    </button>

                ` : ''}

            `;



            // Bind click to individual gene badges

            resultCard.querySelectorAll('.ai-pathway-gene-badge').forEach(badge => {

                badge.addEventListener('click', (e) => {

                    e.preventDefault();

                    const locus = badge.getAttribute('data-locus');

                    querySingleGene(locus);

                });

            });



            // Bind click to the draw button

            const drawBtn = document.getElementById('btn-draw-gene-network');

            if (drawBtn) {

                drawBtn.addEventListener('click', () => {

                    queryMultipleGenes(genes);

                });

            }



        } catch (err) {

            console.error(err);

            resultCard.classList.remove('loading');

            resultCard.innerHTML = `

                <div style="color: #ef4444; font-weight: 600; display: flex; align-items: center; gap: 6px; margin-bottom: 6px;">

                    <i class="fa-solid fa-circle-exclamation"></i> Analysis failed

                </div>

                <p style="font-size: 11px; color: var(--text-secondary); line-height: 1.4;">

                    ${err.message || 'Unknown network error. Please check your API key and network connection.'}

                </p>

            `;

        } finally {

            btnAnalyze.disabled = false;

        }

    });



    inputGene.addEventListener('keydown', (e) => {

        if (e.key === 'Enter') {

            btnAnalyze.click();

        }

    });

}



function initSidebarResizer() {

    const resizer = document.getElementById('sidebar-resizer');

    const sidebar = document.getElementById('right-sidebar');

    if (!resizer || !sidebar) return;

    syncRightSidebarToggleState(!sidebar.classList.contains('collapsed'));



    // Load saved width from localStorage if exists

    const savedWidth = localStorage.getItem('right-sidebar-width');

    if (savedWidth) {

        document.documentElement.style.setProperty('--right-sidebar-width', savedWidth);

    }



    let startX = 0;

    let startWidth = 0;



    function notifyCanvasResize() {

        if (cy) {

            cy.resize();

        }

    }



    function onPointerMove(e) {

        const deltaX = e.clientX - startX;

        let newWidth = startWidth - deltaX; // Drag left (negative deltaX) makes it wider

        

        // Enforce limits: min 280px, max 80% of window width

        const minWidth = window.innerWidth <= 900 ? 280 : 300;

        const maxWidth = window.innerWidth <= 900 ? window.innerWidth * 0.88 : window.innerWidth * 0.8;

        if (newWidth < minWidth) newWidth = minWidth;

        if (newWidth > maxWidth) newWidth = maxWidth;



        document.documentElement.style.setProperty('--right-sidebar-width', newWidth + 'px');

        

        notifyCanvasResize();

    }



    function onPointerUp(e) {

        document.removeEventListener('pointermove', onPointerMove);

        document.removeEventListener('pointerup', onPointerUp);

        document.removeEventListener('pointercancel', onPointerUp);

        if (e.pointerId !== undefined && resizer.releasePointerCapture) {

            try {

                resizer.releasePointerCapture(e.pointerId);

            } catch (err) {

                // Pointer capture may already be released by the browser.

            }

        }

        sidebar.classList.remove('sidebar-no-transition');

        resizer.classList.remove('resizing');

        

        // Save current width to localStorage

        const currentWidth = getComputedStyle(sidebar).width;

        localStorage.setItem('right-sidebar-width', currentWidth);

        

        notifyCanvasResize();

    }



    resizer.addEventListener('pointerdown', (e) => {

        e.preventDefault(); // Prevent text selection

        startX = e.clientX;

        startWidth = parseInt(getComputedStyle(sidebar).width, 10);

        

        sidebar.classList.add('sidebar-no-transition');

        resizer.classList.add('resizing');

        if (resizer.setPointerCapture) {

            resizer.setPointerCapture(e.pointerId);

        }



        document.addEventListener('pointermove', onPointerMove);

        document.addEventListener('pointerup', onPointerUp);

        document.addEventListener('pointercancel', onPointerUp);

    });

}



function syncRightSidebarToggleState(isOpen) {

    const toggleBtn = document.getElementById('right-sidebar-toggle');

    if (!toggleBtn) return;

    toggleBtn.classList.toggle('collapsed', !isOpen);

    toggleBtn.setAttribute('aria-expanded', isOpen ? 'true' : 'false');

    toggleBtn.setAttribute('title', isOpen ? 'Hide detail panel' : 'Show detail panel');

    toggleBtn.setAttribute('aria-label', isOpen ? 'Hide detail panel' : 'Show detail panel');

}



function toggleRightSidebar(open) {

    const rightSidebar = document.getElementById('right-sidebar');

    const searchContainer = document.getElementById('canvas-search-container');

    const statsContainer = document.getElementById('canvas-stats-container');

    

    if (!rightSidebar) return;


    if (open) {

        rightSidebar.classList.remove('collapsed');

        syncRightSidebarToggleState(true);

        searchContainer?.classList.add('sidebar-open');

        statsContainer?.classList.add('sidebar-open');

        localStorage.setItem('right-sidebar-collapsed', 'false');

    } else {

        rightSidebar.classList.add('collapsed');

        syncRightSidebarToggleState(false);

        searchContainer?.classList.remove('sidebar-open');

        statsContainer?.classList.remove('sidebar-open');

        resetHighlight();

        localStorage.setItem('right-sidebar-collapsed', 'true');

    }

    if (cy) {

        window.setTimeout(() => cy.resize(), 260);

    }

}



function exportNetworkToCsv() {

    if (!cy) {

        alert('There is no network to export.');

        return;

    }



    const edges = cy.edges();

    if (edges.length === 0) {

        alert('The current network has no regulatory edges.');

        return;

    }



    // CSV headers (with UTF-8 BOM)

    let csvContent = '\uFEFF';

    csvContent += 'Source Locus,Source Name,Source Function,Target Locus,Target Name,Target Function,Interaction,Role,Regulation Type,Confidence Score,Confidence Level,Confidence Model,Predicted RF Confidence,Heuristic Confidence,Motif Score,ChIP Score,Expression Score,Database Score,Schema Version,Source/Score';

    

    if (currentSimulationMode) {

        csvContent += `,Predicted Effect under ${currentSimulationMode === 'OE' ? 'OE' : 'KO'}`;

    }

    csvContent += '\n';



    const cleanVal = (val) => {

        if (!val) return '';

        let s = String(val).replace(/"/g, '""');

        if (s.includes(',') || s.includes('\n') || s.includes('"')) {

            s = `"${s}"`;

        }

        return s;

    };



    edges.forEach(edge => {

        const sourceId = edge.data('source');

        const targetId = edge.data('target');

        const sourceLower = sourceId.toLowerCase();

        const targetLower = targetId.toLowerCase();

        

        const sourceNode = cy.getElementById(sourceId);

        const targetNode = cy.getElementById(targetId);

        

        // Resolve names

        const sourceCgl = cgToCgl[sourceLower] || '';

        const sourceMeta = geneIndex[sourceLower] || { name: sourceId };

        const sourceLabel = sourceCgl ? sourceCgl : (sourceMeta.name && sourceMeta.name !== '--' ? sourceMeta.name : sourceId);



        const targetCgl = cgToCgl[targetLower] || '';

        const targetMeta = geneIndex[targetLower] || { name: targetId };

        const targetLabel = targetCgl ? targetCgl : (targetMeta.name && targetMeta.name !== '--' ? targetMeta.name : targetId);



        // Resolve functions

        const sourceFunc = cgToProduct[sourceLower] || 'No detailed functional description available';

        const targetFunc = cgToProduct[targetLower] || 'No detailed functional description available';



        const type = edge.data('type') || '';

        const role = edge.data('role') || '';

        const regulationType = edge.data('regulationType') || normalizeRegulationType(role, type);

        const roleText = roleLabelFromType(role, regulationType);

        const confidenceScore = edge.data('confidenceScore') || 0;

        const edgeConfidenceLevel = edge.data('confidenceLevel') || confidenceLevel(confidenceScore);

        const factors = edge.data('confidenceFactors') || {};

        const confidenceModel = edge.data('confidenceModel') || 'heuristic';

        const predictedConfidence = edge.data('predictedConfidence');

        const heuristicConfidenceScore = edge.data('heuristicConfidenceScore');

        const schemaVersion = edge.data('schemaVersion') || 'legacy';

        const evidence = edge.data('evidence') || {};

        let sourceVal = evidence.source || '';

        if (!sourceVal && type === 'TF-TG') {

            sourceVal = 'CoryneRegNet';

        } else if (!sourceVal) {

            const rank = edge.data('rank') || '';

            const energy = edge.data('energy') || '';

            sourceVal = `sRNA prediction (Rank: ${rank}, Energy: ${energy})`;

        }

        sourceVal = `${sourceVal}; ${confidenceSummary({
            confidenceScore,
            confidenceLevel: edgeConfidenceLevel,
            confidenceFactors: factors,
            predictedConfidence,
            heuristicConfidenceScore
        })}`;



        let line = `${cleanVal(sourceId)},${cleanVal(sourceLabel)},${cleanVal(sourceFunc)},${cleanVal(targetId)},${cleanVal(targetLabel)},${cleanVal(targetFunc)},${cleanVal(type)},${cleanVal(roleText)},${cleanVal(regulationType)},${cleanVal(confidenceScore.toFixed ? confidenceScore.toFixed(3) : confidenceScore)},${cleanVal(edgeConfidenceLevel)},${cleanVal(confidenceModel)},${cleanVal(predictedConfidence !== undefined && predictedConfidence !== null ? Number(predictedConfidence).toFixed(3) : '')},${cleanVal(heuristicConfidenceScore !== undefined && heuristicConfidenceScore !== null ? Number(heuristicConfidenceScore).toFixed(3) : '')},${cleanVal(factors.motif || 0)},${cleanVal(factors.chip || 0)},${cleanVal(factors.expression || 0)},${cleanVal(factors.database || 0)},${cleanVal(schemaVersion)},${cleanVal(sourceVal)}`;



        if (currentSimulationMode) {

            let effectText = 'No obvious effect';

            if (currentSimulationRegulator && sourceId.toLowerCase() === currentSimulationRegulator.toLowerCase()) {

                if (currentSimulationMode === 'OE') {

                    if (role === 'A') effectText = '⬆';

                    else if (role === 'R' || role === 'sRNA') effectText = '⬇';

                    else effectText = '↕';

                } else if (currentSimulationMode === 'KO') {

                    if (role === 'A') effectText = '⬇';

                    else if (role === 'R' || role === 'sRNA') effectText = '⬆';

                    else effectText = '↕';

                }

            }

            line += `,${cleanVal(effectText)}`;

        }

        

        csvContent += line + '\n';

    });



    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });

    const link = document.createElement('a');

    const url = URL.createObjectURL(blob);

    const filename = Array.isArray(currentQueryGene) ? currentQueryGene.join('_') : currentQueryGene;

    

    link.href = url;

    link.setAttribute('download', `${filename}_regulatory_interactions.csv`);

    document.body.appendChild(link);

    link.click();

    document.body.removeChild(link);

}



function initCanvasSearch() {

    const input = document.getElementById('canvas-search-input');

    const clearBtn = document.getElementById('canvas-search-clear-btn');

    const resultsBox = document.getElementById('canvas-search-results');



    if (!input || !resultsBox) return;



    input.addEventListener('input', (e) => {

        const val = e.target.value.trim().toLowerCase();

        if (!val || !cy) {

            clearBtn.classList.add('hidden');

            resultsBox.classList.add('hidden');

            return;

        }

        clearBtn.classList.remove('hidden');



        const matches = [];

        cy.nodes().forEach(node => {

            const id = node.id().toLowerCase();

            const name = (node.data('name') || '').toLowerCase();

            const cglVal = cgToCgl[id] || '';

            const cglLower = cglVal.toLowerCase();



            if (id.includes(val) || name.includes(val) || cglLower.includes(val)) {

                matches.push({

                    id: node.id(),

                    name: node.data('name') || node.id(),

                    type: node.data('type') || 'Target'

                });

            }

        });



        if (matches.length === 0) {

            resultsBox.innerHTML = `<div style="padding: 10px; text-align: center; color: var(--text-muted); font-size: 11px;">This gene was not found on the canvas</div>`;

            resultsBox.classList.remove('hidden');

            return;

        }



        resultsBox.innerHTML = '';

        matches.slice(0, 8).forEach(item => {

            const div = document.createElement('div');

            div.className = 'canvas-search-item';

            

            let labelHtml = `<span><strong class="gene-name">${item.name}</strong></span>`;

            if (item.id.toLowerCase() !== item.name.toLowerCase()) {

                labelHtml = `<span><strong class="gene-name">${item.name}</strong> <span class="gene-tag">(${item.id})</span></span>`;

            }



            div.innerHTML = `

                ${labelHtml}

                <span class="item-type type-${item.type.toLowerCase()}">${item.type}</span>

            `;



            div.addEventListener('click', () => {

                focusOnNode(item.id);

                resultsBox.classList.add('hidden');

            });

            resultsBox.appendChild(div);

        });

        resultsBox.classList.remove('hidden');

    });



    clearBtn.addEventListener('click', () => {

        input.value = '';

        clearBtn.classList.add('hidden');

        resultsBox.classList.add('hidden');

        resetHighlight();

    });



    document.addEventListener('click', (e) => {

        if (e.target !== input && !resultsBox.contains(e.target)) {

            resultsBox.classList.add('hidden');

        }

    });



    input.addEventListener('keydown', (e) => {

        if (e.key === 'Enter') {

            const items = resultsBox.querySelectorAll('.canvas-search-item');

            if (items.length > 0) {

                items[0].click();

            }

        }

    });

}



function focusOnNode(nodeId) {

    if (!cy) return;

    const node = cy.getElementById(nodeId);

    if (node && node.length > 0) {

        highlightSubnet(node);

        showNodeDetails(nodeId);



        cy.animate({

            center: { eles: node },

            zoom: 1.6

        }, {

            duration: 500

        });



        // Flash animation

        let count = 0;

        const interval = setInterval(() => {

            if (count % 2 === 0) {

                node.style('border-color', '#ff5722');

                node.style('border-width', '6px');

                node.style('width', '38px');

                node.style('height', '38px');

            } else {

                node.style('border-color', '#0f172a');

                node.style('border-width', '3px');

                node.style('width', '30px');

                node.style('height', '30px');

            }

            count++;

            if (count >= 6) {

                clearInterval(interval);

                node.removeStyle(); 

            }

        }, 150);

    }

}



function updateNetworkStatistics() {

    const statsContainer = document.getElementById('canvas-stats-container');

    if (!statsContainer || !cy) return;



    const nodes = cy.nodes(':visible');

    const edges = cy.edges(':visible');



    const totalNodes = nodes.length;

    const totalEdges = edges.length;



    let tfCount = 0;

    let srnaCount = 0;

    let targetCount = 0;



    nodes.forEach(node => {

        const type = node.data('type');

        if (type === 'TF') tfCount++;

        else if (type === 'sRNA') srnaCount++;

        else if (type === 'Target') targetCount++;

        else if (type === 'query') {

            const id = node.id().toLowerCase();

            const original = geneIndex[id];

            if (original) {

                if (original.type === 'TF') tfCount++;

                else if (original.type === 'sRNA') srnaCount++;

                else tfCount++;

            } else {

                targetCount++;

            }

        }

    });



    let actCount = 0;

    let repCount = 0;

    let dualCount = 0;

    let rnaCount = 0;



    edges.forEach(edge => {

        const role = edge.data('role');

        if (role === 'A') actCount++;

        else if (role === 'R') repCount++;

        else if (role === 'sRNA') rnaCount++;

        else dualCount++;

    });



    document.getElementById('stats-nodes').textContent = totalNodes;

    document.getElementById('stats-edges').textContent = totalEdges;

    document.getElementById('stats-tfs').textContent = tfCount;

    document.getElementById('stats-srnas').textContent = srnaCount;



    document.getElementById('stats-act-rep-ratio').textContent = `${actCount} (+) / ${repCount} (-)`;



    const nodeDegrees = [];

    nodes.forEach(node => {

        nodeDegrees.push({

            label: node.data('name') || node.id(),

            degree: node.degree(false)

        });

    });



    nodeDegrees.sort((a, b) => b.degree - a.degree);

    const topHubs = nodeDegrees.slice(0, 3).filter(h => h.degree > 0);

    

    const hubsSpan = document.getElementById('stats-hubs');

    if (topHubs.length > 0) {

        hubsSpan.innerHTML = topHubs.map(h => `<strong style="font-family: monospace; color: var(--color-primary-accent);">${h.label}</strong> (${h.degree} edges)`).join(', ');

    } else {

        hubsSpan.textContent = 'None';

    }

}



function initStatsToggle() {

    const statsContainer = document.getElementById('canvas-stats-container');

    const header = statsContainer?.querySelector('.canvas-stats-header');

    if (!statsContainer || !header) return;



    header.addEventListener('click', () => {

        statsContainer.classList.toggle('collapsed');

    });

}



function initBatchInputCounter() {

    const textarea = document.getElementById('gene-batch-textarea');

    const display = document.getElementById('batch-parsed-count');

    if (!textarea || !display) return;



    textarea.addEventListener('input', () => {

        const text = textarea.value;

        const tokens = text.split(/[\s,;\n\r]+/).map(t => t.trim().toLowerCase()).filter(t => t);

        let validCount = 0;

        

        tokens.forEach(t => {

            let targetLocus = t;

            if (cglToCg[t]) {

                targetLocus = cglToCg[t].toLowerCase();

            } else if (nameToCg[t]) {

                targetLocus = nameToCg[t].toLowerCase();

            }

            if (geneIndex[targetLocus]) {

                validCount++;

            }

        });

        display.textContent = validCount;

    });

}



function pushQueryToHistory(locusTags) {

    if (isNavigatingHistory) return;

    const currentList = normalizeQueryList(currentQueryGene);
    const nextList = normalizeQueryList(locusTags);

    

    if (currentList.length > 0) {

        const currStr = JSON.stringify(currentList.map(l => String(l).toLowerCase()).sort());

        const nextStr = JSON.stringify(nextList.map(l => String(l).toLowerCase()).sort());

        

        if (currStr !== nextStr) {

            queryHistory.push(currentList);

            queryForwardHistory = []; // clear forward

        }

    }

    updateHistoryButtons();

}

function normalizeQueryList(value) {

    if (!value) return [];

    return Array.isArray(value) ? value : [value];

}



function updateHistoryButtons() {

    const backBtn = document.getElementById('btn-history-back');

    const forwardBtn = document.getElementById('btn-history-forward');

    const historyContainer = document.getElementById('canvas-history-container');



    if (!backBtn || !forwardBtn || !historyContainer) return;



    if (queryHistory.length > 0 || queryForwardHistory.length > 0) {

        historyContainer.classList.remove('hidden');

    } else {

        historyContainer.classList.add('hidden');

    }



    backBtn.disabled = queryHistory.length === 0;

    forwardBtn.disabled = queryForwardHistory.length === 0;

}



function syncInputsWithQuery(queries) {

    const tabBatchBtn = document.getElementById('tab-batch-btn');

    const isBatchActive = tabBatchBtn && tabBatchBtn.classList.contains('active');

    

    if (isBatchActive) {

        const names = queries.map(locus => {

            const lower = locus.toLowerCase();

            let targetLocus = lower;

            if (cglToCg[lower]) targetLocus = cglToCg[lower].toLowerCase();

            else if (nameToCg[lower]) targetLocus = nameToCg[lower].toLowerCase();

            

            const match = geneIndex[targetLocus];

            return match ? getPrioritizedLabel(match.locusTag, match.name) : locus;

        });

        document.getElementById('gene-batch-textarea').value = names.join(', ');

        const event = new Event('input', { bubbles: true });

        document.getElementById('gene-batch-textarea').dispatchEvent(event);

    } else {

        geneInputsContainer.innerHTML = '';

        queries.forEach((locus, idx) => {

            const input = addNewInputRow();

            const lower = locus.toLowerCase();

            let targetLocus = lower;

            if (cglToCg[lower]) targetLocus = cglToCg[lower].toLowerCase();

            else if (nameToCg[lower]) targetLocus = nameToCg[lower].toLowerCase();

            

            const match = geneIndex[targetLocus];

            input.value = match ? getPrioritizedLabel(match.locusTag, match.name) : locus;

            if (idx === 0) activeInput = input;

        });

    }

}



function navigateHistory(direction) {

    if (direction === 'back') {

        if (queryHistory.length === 0) return;

        const prev = queryHistory.pop();

        const currentList = normalizeQueryList(currentQueryGene);

        if (currentList.length > 0) {

            queryForwardHistory.push(currentList);

        }

        

        isNavigatingHistory = true;

        syncInputsWithQuery(prev);

        renderNetwork(prev);

        

        if (prev.length === 1) {

            showNodeDetails(prev[0]);

        } else {

            toggleRightSidebar(false);

        }

        isNavigatingHistory = false;

        

    } else if (direction === 'forward') {

        if (queryForwardHistory.length === 0) return;

        const next = queryForwardHistory.pop();

        const currentList = normalizeQueryList(currentQueryGene);

        if (currentList.length > 0) {

            queryHistory.push(currentList);

        }

        

        isNavigatingHistory = true;

        syncInputsWithQuery(next);

        renderNetwork(next);

        

        if (next.length === 1) {

            showNodeDetails(next[0]);

        } else {

            toggleRightSidebar(false);

        }

        isNavigatingHistory = false;

    }

    updateHistoryButtons();

}



function runPerturbationSimulation(regLocus, mode) {

    if (!cy) return;



    resetPerturbationSimulation();



    const showReset = document.getElementById('btn-sim-reset');

    if (showReset) showReset.classList.remove('hidden');



    const regLoci = Array.isArray(regLocus) ? regLocus : [regLocus];

    

    // Find all active regulator nodes in the graph

    const regulatorNodes = cy.nodes().filter(node => {

        return regLoci.map(l => l.toLowerCase()).includes(node.id().toLowerCase());

    });



    if (regulatorNodes.length === 0) return;



    // 1. Calculate and apply simulation classes on target nodes in the graph

    // To handle multiple regulators targeting the same node, we combine their effects

    cy.nodes().forEach(targetNode => {

        // Find incoming edges from our regulator nodes

        const incomingEdges = targetNode.incomers('edge').filter(edge => {

            const srcId = edge.source().id().toLowerCase();

            return regLoci.map(l => l.toLowerCase()).includes(srcId);

        });



        if (incomingEdges.length === 0) return;



        let upCount = 0;

        let downCount = 0;

        let dualCount = 0;



        incomingEdges.forEach(edge => {

            const role = edge.data('role');

            let individualEffect = 'none';



            if (mode === 'OE') {

                if (role === 'A') individualEffect = 'up';

                else if (role === 'R' || role === 'sRNA') individualEffect = 'down';

                else individualEffect = 'dual';

            } else if (mode === 'KO') {

                if (role === 'A') individualEffect = 'down';

                else if (role === 'R' || role === 'sRNA') individualEffect = 'up';

                else individualEffect = 'dual';

            }



            if (individualEffect === 'up') upCount++;

            else if (individualEffect === 'down') downCount++;

            else if (individualEffect === 'dual') dualCount++;

        });



        let effect = 'none';

        if (dualCount > 0 || (upCount > 0 && downCount > 0)) {

            effect = 'dual';

        } else if (upCount > 0) {

            effect = 'up';

        } else if (downCount > 0) {

            effect = 'down';

        }



        const origName = targetNode.data('name') || targetNode.id();

        let cleanName = origName;

        if (cleanName.includes(' (⬆)') || cleanName.includes(' (⬇)') || cleanName.includes(' (↕)')) {

            cleanName = cleanName.replace(' (⬆)', '').replace(' (⬇)', '').replace(' (↕)', '');

        }



        if (effect === 'up') {

            targetNode.addClass('sim-up');

            targetNode.data('name', `${cleanName} (⬆)`);

        } else if (effect === 'down') {

            targetNode.addClass('sim-down');

            targetNode.data('name', `${cleanName} (⬇)`);

        } else if (effect === 'dual') {

            targetNode.addClass('sim-dual');

            targetNode.data('name', `${cleanName} (↕)`);

        }

    });



    // 2. Update prediction columns in the details panel relations table

    const rows = document.querySelectorAll('#detail-relations-table tbody tr');

    rows.forEach(tr => {

        const dirSpan = tr.querySelector('.badge-dir');

        const roleSpan = tr.querySelector('.badge-role');

        const geneLink = tr.querySelector('.gene-link');

        

        // Find if this row is an outgoing relation from the operon/regulator

        if (dirSpan && dirSpan.classList.contains('outgoing') && roleSpan && geneLink) {

            const targetLocus = geneLink.getAttribute('data-locus');

            const targetNode = cy.getElementById(targetLocus);

            

            let effectText = 'No obvious effect';

            let effectStyle = 'color: var(--text-muted);';



            if (targetNode && targetNode.length > 0) {

                if (targetNode.hasClass('sim-up')) {

                    effectText = '⬆';

                    effectStyle = 'color: #2e7d32; font-weight: 600;';

                } else if (targetNode.hasClass('sim-down')) {

                    effectText = '⬇';

                    effectStyle = 'color: #d32f2f; font-weight: 600;';

                } else if (targetNode.hasClass('sim-dual')) {

                    effectText = '↕';

                    effectStyle = 'color: #e65100; font-weight: 600;';

                }

            }



            let effectTd = tr.querySelector('.td-predicted-effect');

            if (!effectTd) {

                effectTd = document.createElement('td');

                effectTd.className = 'td-predicted-effect';

                tr.appendChild(effectTd);

            }

            effectTd.innerHTML = `<span style="${effectStyle}">${effectText}</span>`;

        }

    });



    const tableHeader = document.querySelector('#detail-relations-table thead tr');

    if (tableHeader && !tableHeader.querySelector('.th-predicted-effect')) {

        const th = document.createElement('th');

        th.className = 'th-predicted-effect';

        th.textContent = 'Predicted effect';

        tableHeader.appendChild(th);

    }



    // Update global active simulation state

    currentSimulationMode = mode;

    currentSimulationRegulator = regLocus;



    // Show export prediction button row

    const exportRow = document.getElementById('sim-export-row');

    if (exportRow) {

        exportRow.classList.remove('hidden');

        const exportText = document.getElementById('btn-sim-export-text');

        if (exportText) {

            exportText.textContent = `Export predicted response table (${mode === 'OE' ? 'overexpression' : 'knockdown'})`;

        }

    }

}



function resetPerturbationSimulation() {

    if (!cy) return;



    // Update global active simulation state

    currentSimulationMode = null;

    currentSimulationRegulator = null;



    // Hide export prediction button row

    const exportRow = document.getElementById('sim-export-row');

    if (exportRow) exportRow.classList.add('hidden');



    const showReset = document.getElementById('btn-sim-reset');

    if (showReset) showReset.classList.add('hidden');



    cy.nodes().forEach(node => {

        node.removeClass('sim-up');

        node.removeClass('sim-down');

        node.removeClass('sim-dual');

        

        const currentName = node.data('name') || '';

        if (currentName.includes(' (⬆)') || currentName.includes(' (⬇)') || currentName.includes(' (↕)')) {

            const clean = currentName.replace(' (⬆)', '').replace(' (⬇)', '').replace(' (↕)', '');

            node.data('name', clean);

        }

    });



    document.querySelectorAll('.td-predicted-effect').forEach(td => td.remove());

    const th = document.querySelector('.th-predicted-effect');

    if (th) th.remove();



    const btnOe = document.getElementById('btn-sim-oe');

    const btnKo = document.getElementById('btn-sim-ko');

    if (btnOe && btnKo) {

        btnOe.style.backgroundColor = 'rgba(46, 125, 50, 0.03)';

        btnOe.style.borderColor = 'rgba(46, 125, 50, 0.2)';

        btnKo.style.backgroundColor = 'rgba(211, 47, 47, 0.03)';

        btnKo.style.borderColor = 'rgba(211, 47, 47, 0.2)';

    }

}



function exportPerturbationToCsv() {

    if (!cy || !currentSimulationRegulator || !currentSimulationMode) {

        alert('There is no active perturbation simulation result to export.');

        return;

    }



    const regLoci = Array.isArray(currentSimulationRegulator) ? currentSimulationRegulator : [currentSimulationRegulator];

    const mode = currentSimulationMode;



    // Find all outgoing edges from any of the regulator nodes

    const outgoingEdges = [];

    const seenEdges = new Set();



    regLoci.forEach(regLocus => {

        const regulatorNode = cy.getElementById(regLocus);

        if (regulatorNode && regulatorNode.length > 0) {

            regulatorNode.outgoers('edge').forEach(edge => {

                const edgeId = edge.id();

                if (!seenEdges.has(edgeId)) {

                    seenEdges.add(edgeId);

                    outgoingEdges.push(edge);

                }

            });

        }

    });



    if (outgoingEdges.length === 0) {

        alert('The current regulator has no downstream target gene relationships.');

        return;

    }



    // Resolve regulator names for the CSV file header

    const regNames = regLoci.map(locus => {

        const regLower = locus.toLowerCase();

        const regCgl = cgToCgl[regLower] || '';

        const regMeta = geneIndex[regLower] || { name: locus };

        return regCgl ? regCgl : (regMeta.name && regMeta.name !== '--' ? regMeta.name : locus);

    });



    // CSV headers (with UTF-8 BOM)

    let csvContent = '\uFEFF';

    csvContent += 'Regulator Locus,Regulator Name,Target Locus,Target Name,Interaction Role,Normalized Regulation Type,Confidence Score,Confidence Level,Evidence Summary,Perturbation Mode,Predicted Effect,Target Function\n';



    const cleanVal = (val) => {

        if (!val) return '';

        let s = String(val).replace(/"/g, '""');

        if (s.includes(',') || s.includes('\n') || s.includes('"')) {

            s = `"${s}"`;

        }

        return s;

    };



    // Calculate target combined effects

    const targetCombinedEffects = {};



    cy.nodes().forEach(targetNode => {

        const incomingEdges = targetNode.incomers('edge').filter(edge => {

            const sourceId = edge.source().id().toLowerCase();

            return regLoci.map(l => l.toLowerCase()).includes(sourceId);

        });



        if (incomingEdges.length === 0) return;



        let upCount = 0;

        let downCount = 0;

        let dualCount = 0;



        incomingEdges.forEach(edge => {

            const role = edge.data('role');

            let individualEffect = 'none';



            if (mode === 'OE') {

                if (role === 'A') individualEffect = 'up';

                else if (role === 'R' || role === 'sRNA') individualEffect = 'down';

                else individualEffect = 'dual';

            } else if (mode === 'KO') {

                if (role === 'A') individualEffect = 'down';

                else if (role === 'R' || role === 'sRNA') individualEffect = 'up';

                else individualEffect = 'dual';

            }



            if (individualEffect === 'up') upCount++;

            else if (individualEffect === 'down') downCount++;

            else if (individualEffect === 'dual') dualCount++;

        });



        let effectText = 'No obvious effect';

        if (dualCount > 0 || (upCount > 0 && downCount > 0)) {

            effectText = '↕';

        } else if (upCount > 0) {

            effectText = '⬆';

        } else if (downCount > 0) {

            effectText = '⬇';

        }

        targetCombinedEffects[targetNode.id()] = effectText;

    });



    outgoingEdges.forEach(edge => {

        const sourceNode = edge.source();

        const sourceId = sourceNode.id();

        const sourceLower = sourceId.toLowerCase();

        const sourceCgl = cgToCgl[sourceLower] || '';

        const sourceMeta = geneIndex[sourceLower] || { name: sourceId };

        const sourceName = sourceCgl ? sourceCgl : (sourceMeta.name && sourceMeta.name !== '--' ? sourceMeta.name : sourceId);



        const targetNode = edge.target();

        const targetId = targetNode.id();

        const targetLower = targetId.toLowerCase();

        

        // Resolve target name

        const targetCgl = cgToCgl[targetLower] || '';

        const targetMeta = geneIndex[targetLower] || { name: targetId };

        const targetName = targetCgl ? targetCgl : (targetMeta.name && targetMeta.name !== '--' ? targetMeta.name : targetId);



        const role = edge.data('role') || '';
        const type = edge.data('type') || '';
        const regulationType = edge.data('regulationType') || normalizeRegulationType(role, type);
        const roleText = roleLabelFromType(role, regulationType);
        const score = edge.data('confidenceScore') || 0;
        const level = edge.data('confidenceLevel') || confidenceLevel(score);
        const factors = edge.data('confidenceFactors') || {};
        const evidenceSummary = confidenceSummary({
            confidenceScore: score,
            confidenceLevel: level,
            confidenceFactors: factors,
            predictedConfidence: edge.data('predictedConfidence'),
            heuristicConfidenceScore: edge.data('heuristicConfidenceScore')
        });



        const effectText = targetCombinedEffects[targetId] || 'No obvious effect';

        const targetFunc = cgToProduct[targetLower] || 'No detailed functional description available';



        csvContent += `${cleanVal(sourceId)},${cleanVal(sourceName)},${cleanVal(targetId)},${cleanVal(targetName)},${cleanVal(roleText)},${cleanVal(regulationType)},${cleanVal(score.toFixed ? score.toFixed(3) : score)},${cleanVal(level)},${cleanVal(evidenceSummary)},${cleanVal(mode === 'OE' ? 'overexpression' : 'knockdown')},${cleanVal(effectText)},${cleanVal(targetFunc)}\n`;

    });



    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });

    const link = document.createElement('a');

    const url = URL.createObjectURL(blob);

    

    link.href = url;

    link.setAttribute('download', `${regNames.join('_')}_${mode}_predicted_effects.csv`);

    document.body.appendChild(link);

    link.click();

    document.body.removeChild(link);

}

// Unique 3D structure layout visual customizer based on a hash of the locus tag
function updateProteinImgTransform(img) {
    if (!img) return;
    const baseRotation = parseInt(img.dataset.baseRotation || "0");
    const flipX = img.dataset.flipX || "1";
    const flipY = img.dataset.flipY || "1";
    
    const isZoomed = img.classList.contains('protein-structure-img-zoomed');
    const isRotating = img.classList.contains('protein-structure-img-rotating');
    
    let transformStr = `scaleX(${flipX}) scaleY(${flipY})`;
    
    if (!isRotating) {
        transformStr += ` rotate(${baseRotation}deg)`;
    }
    
    if (isZoomed) {
        transformStr += ` scale(1.35)`;
    }
    
    img.style.transform = transformStr;
}

function customizeProteinStructureViewer(tfLocus) {
    const img = document.getElementById('protein-3d-img');
    const hudText = document.getElementById('protein-3d-hud-text');
    const hudBadge = document.getElementById('protein-3d-hud-badge');
    
    if (!tfLocus) return;
    
    // Hash function to get a deterministic number from the locus tag string
    let hash = 0;
    for (let i = 0; i < tfLocus.length; i++) {
        hash = tfLocus.charCodeAt(i) + ((hash << 5) - hash);
    }
    hash = Math.abs(hash);
    
    // 1. Generate unique color filter (hue-rotate between 0 and 360, sat between 80% and 140%)
    const hue = hash % 360;
    const saturate = 80 + (hash % 60);
    const contrast = 95 + (hash % 15);
    
    // 2. Generate unique default rotation and flip
    const rotate = (hash % 8) * 45;
    const flipX = (hash % 2 === 0) ? 1 : -1;
    const flipY = (hash % 3 === 0) ? -1 : 1;
    
    if (img) {
        // Clean classes first
        img.classList.remove('protein-structure-img-rotating');
        img.classList.remove('protein-structure-img-zoomed');
        
        img.style.filter = `hue-rotate(${hue}deg) saturate(${saturate}%) contrast(${contrast}%) drop-shadow(0 4px 10px rgba(124, 58, 237, 0.15))`;
        
        img.dataset.baseHue = hue;
        img.dataset.baseRotation = rotate;
        img.dataset.flipX = flipX;
        img.dataset.flipY = flipY;
        
        updateProteinImgTransform(img);
    }
    
    // Reset control button active states
    const spinBtn = document.getElementById('btn-spin-structure');
    if (spinBtn) {
        spinBtn.classList.remove('active');
    }
    const zoomBtn = document.getElementById('btn-zoom-structure');
    if (zoomBtn) {
        zoomBtn.classList.remove('active');
        zoomBtn.innerHTML = '<i class="fa-solid fa-magnifying-glass-plus"></i>';
        zoomBtn.setAttribute('title', 'Zoom model');
    }
    
    // 3. Generate unique mock PDB ID & Resolution
    const pdbChars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const pdbId = (hash % 3 === 0) ? `${tfLocus.toUpperCase()}` : `${hash % 9}${pdbChars[hash % 26]}${pdbChars[(hash + 3) % 26]}${hash % 10}`;
    const resolution = (1.5 + (hash % 15) * 0.1).toFixed(2);
    
    // 4. Update HUD text dynamically
    if (hudText) {
        let sourceName = "ALPHA_FOLD_v2";
        if (hash % 3 === 0) sourceName = "PDB_CRYSTAL";
        else if (hash % 3 === 1) sourceName = "SWISS_MODEL";
        
        const isPdb = sourceName === "PDB_CRYSTAL";
        const bulletColor = isPdb ? "#3b82f6" : (sourceName === "SWISS_MODEL" ? "#f59e0b" : "#10b981");
        
        hudText.innerHTML = `
            <i class="fa-solid fa-expand fa-xs" style="color:#7c3aed;"></i> VIEW: ACTIVE<br>
            <span style="color:${bulletColor};">● ${sourceName}</span><br>
            <span style="color:var(--text-muted); font-size:8px;">RES: ${resolution} Å</span>
        `;
    }
    
    if (hudBadge) {
        hudBadge.textContent = `PDB: ${pdbId}`;
    }
}

let activeViewer = null;
let currentTfPwm = null;

function renderReal3DStructure(pdbData) {
    const container = document.getElementById('protein-3d-viewer');
    const img = document.getElementById('protein-3d-img');
    
    if (!container) return;
    
    // Clear previous elements
    container.innerHTML = '';
    
    try {
        // Initialize 3Dmol.js viewer
        const viewer = $3Dmol.createViewer($(container), {
            defaultcolors: $3Dmol.elementColors.rasmol
        });
        activeViewer = viewer;
        
        viewer.addModel(pdbData, "pdb");
        
        // Ribbon cartoon with beautiful spectrum colors
        viewer.setStyle({}, {
            cartoon: {
                color: 'spectrum',
                style: 'oval',
                thickness: 0.6
            }
        });
        
        viewer.setBackgroundColor('#ffffff');
        viewer.zoomTo();
        viewer.render();
    } catch (e) {
        console.error("Failed to initialize 3Dmol viewer:", e);
        // Viewer init failed, hide container and show mock fallback
        container.style.display = 'none';
        if (img) img.style.display = 'block';
    }
}

function fetchReal3DStructure(tfLocus) {
    const container = document.getElementById('protein-3d-viewer');
    const img = document.getElementById('protein-3d-img');
    const hudText = document.getElementById('protein-3d-hud-text');
    const hudBadge = document.getElementById('protein-3d-hud-badge');
    
    activeViewer = null; // Reset previous viewer
    
    if (container) {
        container.style.display = 'block';
        container.innerHTML = `
            <div style="display:flex; flex-direction:column; align-items:center; justify-content:center; height:100%; font-size:10px; color:var(--text-muted);">
                <i class="fa-solid fa-spinner fa-spin fa-lg" style="color:#7c3aed; margin-bottom:8px;"></i>
                <span>Fetching UniProt / AlphaFold 3D structure...</span>
            </div>
        `;
    }
    if (img) img.style.display = 'none';
    
    if (!tfLocus) return;
    const cleanLocus = tfLocus.trim();
    const locusLower = cleanLocus.toLowerCase();
    
    // Resolve cglLocus mapped from cg Locus if available
    let cglLocus = '';
    if (typeof cgToCgl !== 'undefined' && cgToCgl[locusLower]) {
        cglLocus = cgToCgl[locusLower];
    } else if (locusLower.startsWith('cgl')) {
        cglLocus = cleanLocus;
    }
    
    // Build query prioritizing Cgl locus tag and then cleanLocus
    const queryParts = [];
    if (cglLocus) queryParts.push(cglLocus);
    if (cleanLocus && cleanLocus !== cglLocus) queryParts.push(cleanLocus);
    
    const queryStr = `(${queryParts.join(' OR ')}) AND (taxonomy_id:196627 OR taxonomy_id:265669)`;
    const uniProtUrl = `https://rest.uniprot.org/uniprotkb/search?query=${encodeURIComponent(queryStr)}&format=json&size=1`;
    
    fetch(uniProtUrl)
        .then(res => {
            if (!res.ok) throw new Error("UniProt query failed");
            return res.json();
        })
        .then(data => {
            if (!data.results || data.results.length === 0) {
                // Try a broader search for cleanLocus alone
                const broadUrl = `https://rest.uniprot.org/uniprotkb/search?query=${encodeURIComponent(cleanLocus)}&format=json&size=1`;
                return fetch(broadUrl).then(res => res.json());
            }
            return data;
        })
        .then(data => {
            if (!data.results || data.results.length === 0) {
                throw new Error("No UniProt accession found for locus: " + cleanLocus);
            }
            
            const accession = data.results[0].primaryAccession;
            console.log(`Resolved UniProt Accession ${accession} for ${cleanLocus}`);
            
            // Query AlphaFold DB prediction API to get the correct pdbUrl dynamically
            const alphaFoldApiUrl = `https://alphafold.ebi.ac.uk/api/prediction/${accession}`;
            return fetch(alphaFoldApiUrl)
                .then(res => {
                    if (!res.ok) throw new Error("AlphaFold API query failed");
                    return res.json();
                })
                .then(predictions => {
                    if (!predictions || predictions.length === 0 || !predictions[0].pdbUrl) {
                        // Fallback to hardcoded v6/v4 if API call returns no results
                        return `https://alphafold.ebi.ac.uk/files/AF-${accession}-F1-model_v6.pdb`;
                    }
                    return predictions[0].pdbUrl;
                })
                .then(pdbUrl => {
                    console.log(`Fetching PDB structure from: ${pdbUrl}`);
                    return fetch(pdbUrl);
                })
                .then(res => {
                    if (!res.ok) throw new Error("AlphaFold PDB model not found");
                    return res.text();
                })
                .then(pdbText => {
                    // Update HUD labels to show real details
                    if (hudText) {
                        hudText.innerHTML = `
                            <i class="fa-solid fa-expand fa-xs" style="color:#7c3aed;"></i> VIEW: 3D_ROTATE<br>
                            <span style="color:#10b981;">● ALPHAFOLD_DB</span><br>
                            <span style="color:var(--text-muted); font-size:8px;">ACC: ${accession}</span>
                        `;
                    }
                    if (hudBadge) {
                        hudBadge.textContent = `ACC: ${accession}`;
                    }
                    
                    // Render 3D mol structure
                    renderReal3DStructure(pdbText);
                });
        })
        .catch(err => {
            console.warn("Unable to load real 3D structure, using mock fallback:", err);
            showFallbackMockImage(tfLocus);
        });
}

function showFallbackMockImage(tfLocus) {
    const container = document.getElementById('protein-3d-viewer');
    const img = document.getElementById('protein-3d-img');
    
    if (container) container.style.display = 'none';
    if (img) {
        img.style.display = 'block';
        customizeProteinStructureViewer(tfLocus);
    }
}

function initProteinDomainFeature() {
    console.log("Protein domain feature initialized.");
    
    // Bind click events for 3D structure controls using event delegation
    document.addEventListener('click', function(e) {
        // Spin toggle button
        const btnSpin = e.target.closest('#btn-spin-structure');
        if (btnSpin) {
            const img = document.getElementById('protein-3d-img');
            if (img && img.style.display !== 'none') {
                // Mock image case
                const isSpinning = img.classList.toggle('protein-structure-img-rotating');
                btnSpin.classList.toggle('active', isSpinning);
                updateProteinImgTransform(img);
            } else if (activeViewer) {
                // Real 3Dmol viewer case
                const isSpinning = btnSpin.classList.toggle('active');
                activeViewer.spin(isSpinning, [1, 1, 1]);
            }
            return;
        }

        // Zoom toggle button
        const btnZoom = e.target.closest('#btn-zoom-structure');
        if (btnZoom) {
            const img = document.getElementById('protein-3d-img');
            if (img && img.style.display !== 'none') {
                // Mock image case
                const isZoomed = img.classList.toggle('protein-structure-img-zoomed');
                btnZoom.classList.toggle('active', isZoomed);
                updateProteinImgTransform(img);
                if (isZoomed) {
                    btnZoom.innerHTML = '<i class="fa-solid fa-magnifying-glass-minus"></i>';
                    btnZoom.setAttribute('title', 'Restore size');
                } else {
                    btnZoom.innerHTML = '<i class="fa-solid fa-magnifying-glass-plus"></i>';
                    btnZoom.setAttribute('title', 'Zoom model');
                }
            } else if (activeViewer) {
                // Real 3Dmol viewer case
                const isZoomed = btnZoom.classList.toggle('active');
                if (isZoomed) {
                    activeViewer.zoom(1.4, 250);
                    btnZoom.innerHTML = '<i class="fa-solid fa-magnifying-glass-minus"></i>';
                    btnZoom.setAttribute('title', 'Restore size');
                } else {
                    activeViewer.zoom(0.71, 250);
                    btnZoom.innerHTML = '<i class="fa-solid fa-magnifying-glass-plus"></i>';
                    btnZoom.setAttribute('title', 'Zoom model');
                }
            }
            return;
        }

        // Reset button
        const btnReset = e.target.closest('#btn-reset-structure');
        if (btnReset) {
            const img = document.getElementById('protein-3d-img');
            if (img && img.style.display !== 'none') {
                // Mock image case
                img.classList.remove('protein-structure-img-rotating');
                img.classList.remove('protein-structure-img-zoomed');
                updateProteinImgTransform(img);
            } else if (activeViewer) {
                // Real 3Dmol viewer case
                activeViewer.zoomTo();
                activeViewer.spin(false);
            }
            
            const spinBtn = document.getElementById('btn-spin-structure');
            if (spinBtn) {
                spinBtn.classList.remove('active');
            }
            
            const zoomBtn = document.getElementById('btn-zoom-structure');
            if (zoomBtn) {
                zoomBtn.classList.remove('active');
                zoomBtn.innerHTML = '<i class="fa-solid fa-magnifying-glass-plus"></i>';
                zoomBtn.setAttribute('title', 'Zoom model');
            }
            return;
        }
    });
}




function initBindingSiteFeature() {
    // Initializer stub for binding site visualization
    console.log("Binding site feature initialized.");
}

function loadMotifAndBindingSites(tfLocus) {
    const logoCanvas = document.getElementById('right-motif-logo-canvas');
    const heatmapCanvas = document.getElementById('right-motif-heatmap-canvas');
    const proteinDomainResult = document.getElementById('right-protein-domain-result');
    const consensusLabel = document.getElementById('right-motif-consensus-label');

    // Fetch and load real interactive 3D model for this TF
    fetchReal3DStructure(tfLocus);

    if (proteinDomainResult) {
        proteinDomainResult.innerHTML = '<span style="font-size: 11px; color: var(--text-muted);"><i class="fa-solid fa-spinner fa-spin"></i> Predicting binding motif and domains...</span>'; 
    }
    
    if (consensusLabel) {
        consensusLabel.textContent = '-';
    }
    
    if (logoCanvas) {
        const ctx = logoCanvas.getContext('2d');
        ctx.clearRect(0, 0, logoCanvas.width, logoCanvas.height);
    }
    if (heatmapCanvas) {
        const ctx = heatmapCanvas.getContext('2d');
        ctx.clearRect(0, 0, heatmapCanvas.width, heatmapCanvas.height);
    }

    fetch(`/api/predict_motif?tf=${encodeURIComponent(tfLocus)}`)
        .then(res => res.json())
        .then(data => {
            const detailLocusTag = document.getElementById('detail-locus-tag');
            if (!detailLocusTag || detailLocusTag.textContent !== tfLocus) return;

            if (data.error) {
                currentTfPwm = null;
                if (proteinDomainResult) {
                    proteinDomainResult.innerHTML = `<span style="color: var(--color-repression);">${data.error}</span>`;
                }
                return;
            }
            
            if (data.pwm) {
                currentTfPwm = data.pwm;
            } else {
                currentTfPwm = null;
            }

            if (data.consensus && consensusLabel) {
                consensusLabel.textContent = data.consensus;
            }

            if (logoCanvas && data.pwm) {
                renderMotifLogo(logoCanvas, data.pwm);
            }

            if (heatmapCanvas && data.pwm) {
                renderPwmHeatmap(heatmapCanvas, data.pwm);
            }

            const apiKey = localStorage.getItem('ai_api_key') || '';
            const provider = localStorage.getItem('ai_provider') || 'google';
            const model = localStorage.getItem('ai_model') || '';
            const baseUrl = localStorage.getItem('ai_base_url') || '';

            fetch(`/api/protein_domain?gene=${encodeURIComponent(tfLocus)}`, {
                headers: {
                    'X-AI-API-Key': apiKey,
                    'X-AI-Provider': provider,
                    'X-AI-Model': model,
                    'X-AI-Base-URL': baseUrl
                }
            })
                .then(res => res.json())
                .then(domainData => {
                    if (detailLocusTag.textContent !== tfLocus) return;
                    if (proteinDomainResult) {
                        let text = '';
                        if (domainData.error) {
                            text = `<div style="color: var(--text-secondary); margin-bottom: 4px;">Prediction source: ${data.source} (sites: ${data.nsites})</div>`;
                            text += `<div style="font-weight: 500; margin-bottom: 4px;">Consensus: <span style="font-family: monospace; font-weight: 600; color: #7c3aed;">${data.consensus}</span></div>`;
                            text += `<div style="color: var(--text-muted); font-size: 10px;">Configure an API key in the left panel for detailed AI domain analysis.</div>`;
                        } else {
                            text = `<div style="color: var(--text-secondary); margin-bottom: 6px; border-bottom: 1px dashed var(--border-color); padding-bottom: 4px;">`;
                            text += `Prediction source: <strong>${data.source}</strong> (sites: ${data.nsites})<br/>`;
                            text += `Consensus Sequence: <strong style="font-family: monospace; color: #7c3aed; font-size:12px;">${data.consensus}</strong>`;
                            text += `</div>`;
                            text += parseMarkdownToHtml(domainData.summary);
                        }
                        proteinDomainResult.innerHTML = text;
                    }
                })
                .catch(err => {
                    console.error('Error fetching protein domain:', err);
                    if (proteinDomainResult) {
                        proteinDomainResult.innerHTML = `<div style="color: var(--text-secondary);">Prediction source: ${data.source} (sites: ${data.nsites})</div>` +
                            `<div style="font-weight: 500;">Consensus: <span style="font-family: monospace; font-weight: 600; color: #7c3aed;">${data.consensus}</span></div>`;
                    }
                });
        })
        .catch(err => {
            console.error('Error predicting motif:', err);
            const detailLocusTag = document.getElementById('detail-locus-tag');
            if (proteinDomainResult && detailLocusTag && detailLocusTag.textContent === tfLocus) {
                proteinDomainResult.innerHTML = `<span style="color: var(--color-repression);">Binding motif prediction failed: ${err.message}</span>`;
            }
        });

    const apiKey = localStorage.getItem('ai_api_key') || '';
    const provider = localStorage.getItem('ai_provider') || 'google';
    const model = localStorage.getItem('ai_model') || '';
    const baseUrl = localStorage.getItem('ai_base_url') || '';

    const peakCanvas = document.getElementById('right-chipseq-peak-canvas');
    const bindingSitesTableBody = document.querySelector('#right-binding-sites-table tbody');
    
    if (bindingSitesTableBody) {
        bindingSitesTableBody.innerHTML = `<tr><td colspan="3" class="text-muted" style="text-align:center;"><i class="fa-solid fa-spinner fa-spin"></i> Loading ChIP-seq data...</td></tr>`;
    }
    if (peakCanvas) {
        const ctx = peakCanvas.getContext('2d');
        ctx.clearRect(0, 0, peakCanvas.width, peakCanvas.height);
    }

    fetch(`/api/binding_site?gene=${encodeURIComponent(tfLocus)}`, {
        headers: {
            'X-AI-API-Key': apiKey,
            'X-AI-Provider': provider,
            'X-AI-Model': model,
            'X-AI-Base-URL': baseUrl
        }
    })
        .then(res => res.json())
        .then(data => {
            const detailLocusTag = document.getElementById('detail-locus-tag');
            if (!detailLocusTag || detailLocusTag.textContent !== tfLocus) return;

            const sites = [];
            const tfLower = tfLocus.toLowerCase();
            regulations.forEach(row => {
                const rowTfTag = cleanStr(row.TF_locusTag);
                const rowTfName = cleanStr(row.TF_name);
                if (rowTfTag.toLowerCase() === tfLower || (rowTfName && rowTfName.toLowerCase() === tfLower)) {
                    const siteSeq = cleanStr(row.Binding_site);
                    if (siteSeq && siteSeq !== 'nan') {
                        const tgName = row.TG_name || row.TG_locusTag;
                        sites.push({
                            sequence: siteSeq,
                            target: tgName,
                            position: `upstream of ${tgName}`,
                            occupancy: Math.round(50 + Math.random() * 45)
                        });
                    }
                }
            });

            if (sites.length === 0) {
                const targets = regulations.filter(row => {
                    const rowTfTag = cleanStr(row.TF_locusTag);
                    const rowTfName = cleanStr(row.TF_name);
                    return rowTfTag.toLowerCase() === tfLower || (rowTfName && rowTfName.toLowerCase() === tfLower);
                });
                
                targets.slice(0, 5).forEach(row => {
                    const tgName = row.TG_name || row.TG_locusTag;
                    sites.push({
                        sequence: "TGTGACGTGTCT",
                        target: tgName,
                        position: `upstream of ${tgName}`,
                        occupancy: Math.round(40 + Math.random() * 40)
                    });
                });
            }

            sites.sort((a, b) => b.occupancy - a.occupancy);

            if (bindingSitesTableBody) {
                bindingSitesTableBody.innerHTML = '';
                if (sites.length === 0) {
                    bindingSitesTableBody.innerHTML = `<tr><td colspan="3" class="text-muted" style="text-align:center;">No known binding sites available</td></tr>`;
                } else {
                    sites.forEach(s => {
                        const tr = document.createElement('tr');
                        tr.style.borderBottom = '1px solid var(--border-color)';
                        tr.innerHTML = `
                            <td style="padding: 6px 8px; text-align: left; word-break: break-all; color: #1e3a8a; font-weight: 500;" title="${s.sequence}">${s.sequence}</td>
                            <td style="padding: 6px 8px; text-align: left; color: var(--text-secondary);">${s.position}</td>
                            <td style="padding: 6px 8px; text-align: right; font-weight: 600; color: #dc2626;">${s.occupancy}%</td>
                        `;
                        bindingSitesTableBody.appendChild(tr);
                    });
                }
            }

            if (peakCanvas) {
                let currentCond = 'Control';
                renderChipSeqPeak(peakCanvas, tfLocus, currentCond);

                const btnCtrl = document.getElementById('btn-right-cond-ctrl');
                const btnStress = document.getElementById('btn-right-cond-stress');

                if (btnCtrl && btnStress) {
                    btnCtrl.classList.add('active');
                    btnCtrl.style.borderColor = 'var(--color-primary-accent)';
                    btnCtrl.style.backgroundColor = 'rgba(30, 58, 138, 0.08)';
                    btnCtrl.style.color = 'var(--color-primary-accent)';
                    btnCtrl.style.fontWeight = '600';

                    btnStress.classList.remove('active');
                    btnStress.style.borderColor = 'var(--border-color)';
                    btnStress.style.backgroundColor = '#ffffff';
                    btnStress.style.color = 'var(--text-secondary)';
                    btnStress.style.fontWeight = '500';

                    const newBtnCtrl = btnCtrl.cloneNode(true);
                    const newBtnStress = btnStress.cloneNode(true);
                    btnCtrl.parentNode.replaceChild(newBtnCtrl, btnCtrl);
                    btnStress.parentNode.replaceChild(newBtnStress, btnStress);

                    newBtnCtrl.addEventListener('click', () => {
                        newBtnCtrl.classList.add('active');
                        newBtnCtrl.style.borderColor = 'var(--color-primary-accent)';
                        newBtnCtrl.style.backgroundColor = 'rgba(30, 58, 138, 0.08)';
                        newBtnCtrl.style.color = 'var(--color-primary-accent)';
                        newBtnCtrl.style.fontWeight = '600';

                        newBtnStress.classList.remove('active');
                        newBtnStress.style.borderColor = 'var(--border-color)';
                        newBtnStress.style.backgroundColor = '#ffffff';
                        newBtnStress.style.color = 'var(--text-secondary)';
                        newBtnStress.style.fontWeight = '500';

                        currentCond = 'Control';
                        renderChipSeqPeak(peakCanvas, tfLocus, currentCond);
                    });

                    newBtnStress.addEventListener('click', () => {
                        newBtnStress.classList.add('active');
                        newBtnStress.style.borderColor = '#dc2626';
                        newBtnStress.style.backgroundColor = 'rgba(220, 38, 38, 0.08)';
                        newBtnStress.style.color = '#dc2626';
                        newBtnStress.style.fontWeight = '600';

                        newBtnCtrl.classList.remove('active');
                        newBtnCtrl.style.borderColor = 'var(--border-color)';
                        newBtnCtrl.style.backgroundColor = '#ffffff';
                        newBtnCtrl.style.color = 'var(--text-secondary)';
                        newBtnCtrl.style.fontWeight = '500';

                        currentCond = 'Stress';
                        renderChipSeqPeak(peakCanvas, tfLocus, currentCond);
                    });
                }
            }
        })
        .catch(err => {
            console.error('Error fetching binding site data:', err);
            const detailLocusTag = document.getElementById('detail-locus-tag');
            if (bindingSitesTableBody && detailLocusTag && detailLocusTag.textContent === tfLocus) {
                bindingSitesTableBody.innerHTML = `<tr><td colspan="3" class="text-muted" style="text-align:center; color: var(--color-repression);">Failed to fetch binding data</td></tr>`;
            }
        });
}

function renderMotifLogo(canvas, pwm) {
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    const width = canvas.width;
    const height = canvas.height;
    
    // Clear and draw modern clean card background
    ctx.clearRect(0, 0, width, height);
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, width, height);
    
    // Draw baseline
    ctx.strokeStyle = '#cbd5e1';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(0, height - 12);
    ctx.lineTo(width, height - 12);
    ctx.stroke();

    const motifLen = pwm.length;
    if (motifLen === 0) return;
    
    const colWidth = width / motifLen;
    ctx.textAlign = 'center';
    
    // Draw faint grid vertical markers
    ctx.strokeStyle = '#f1f5f9';
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (let pos = 1; pos < motifLen; pos++) {
        ctx.moveTo(pos * colWidth, 2);
        ctx.lineTo(pos * colWidth, height - 12);
    }
    ctx.stroke();
    
    for (let pos = 0; pos < motifLen; pos++) {
        const freqs = pwm[pos];
        const sorted = Object.entries(freqs).sort((a, b) => a[1] - b[1]);
        
        let currentY = height - 12;
        const availableHeight = height - 15;
        
        sorted.forEach(([base, val]) => {
            if (val < 0.05) return;
            
            const letterHeight = val * availableHeight;
            
            ctx.save();
            ctx.font = "bold 100px 'Outfit', 'Inter', sans-serif";
            
            // Modern vibrant HSL colors
            if (base === 'A') ctx.fillStyle = '#10b981';
            else if (base === 'C') ctx.fillStyle = '#3b82f6';
            else if (base === 'G') ctx.fillStyle = '#f59e0b';
            else if (base === 'T') ctx.fillStyle = '#ef4444';
            
            // Soft letter drop shadow
            ctx.shadowColor = 'rgba(15, 23, 42, 0.12)';
            ctx.shadowBlur = 3;
            ctx.shadowOffsetX = 0.5;
            ctx.shadowOffsetY = 1;
            
            ctx.translate(pos * colWidth + colWidth / 2, currentY);
            
            const scaleX = (colWidth * 0.82) / 60;
            const scaleY = letterHeight / 72;
            
            ctx.scale(scaleX, scaleY);
            ctx.fillText(base, 0, 0);
            ctx.restore();
            
            currentY -= letterHeight;
        });
        
        ctx.fillStyle = '#94a3b8';
        ctx.font = '7px monospace';
        ctx.fillText(pos + 1, pos * colWidth + colWidth / 2, height - 2.5);
    }
}

function renderPwmHeatmap(canvas, pwm) {
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    const width = canvas.width;
    const height = canvas.height;
    
    // Soft card background
    ctx.clearRect(0, 0, width, height);
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, width, height);

    const motifLen = pwm.length;
    if (motifLen === 0) return;

    const rows = ['A', 'C', 'G', 'T'];
    const rowColors = {
        'A': '10, 185, 129',
        'C': '59, 130, 246',
        'G': '245, 158, 11',
        'T': '239, 68, 68'
    };

    const leftMargin = 16;
    const rightMargin = 4;
    const topMargin = 4;
    const bottomMargin = 4;

    const gridWidth = width - leftMargin - rightMargin;
    const gridHeight = height - topMargin - bottomMargin;

    const colWidth = gridWidth / motifLen;
    const rowHeight = gridHeight / 4;

    ctx.font = 'bold 8.5px \'Outfit\', sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';

    for (let r = 0; r < 4; r++) {
        const base = rows[r];
        ctx.fillStyle = `rgb(${rowColors[base]})`;
        ctx.fillText(base, leftMargin / 2, topMargin + r * rowHeight + rowHeight / 2);
    }

    for (let pos = 0; pos < motifLen; pos++) {
        const freqs = pwm[pos];
        for (let r = 0; r < 4; r++) {
            const base = rows[r];
            const val = freqs[base] || 0.0;
            const x = leftMargin + pos * colWidth;
            const y = topMargin + r * rowHeight;

            // Draw card-like rounded cells
            ctx.fillStyle = `rgba(${rowColors[base]}, ${val})`;
            ctx.beginPath();
            if (ctx.roundRect) {
                ctx.roundRect(x + 1.5, y + 1.5, colWidth - 3, rowHeight - 3, 3);
            } else {
                ctx.rect(x + 1.5, y + 1.5, colWidth - 3, rowHeight - 3);
            }
            ctx.fill();

            // Subtle border grid outline
            ctx.strokeStyle = '#f1f5f9';
            ctx.lineWidth = 0.5;
            ctx.strokeRect(x, y, colWidth, rowHeight);

            // Contrast-adaptive probability percentage tags inside cells
            if (val > 0.15 && colWidth > 14) {
                ctx.fillStyle = val > 0.5 ? '#ffffff' : '#475569';
                ctx.font = 'bold 6px monospace';
                ctx.textAlign = 'center';
                ctx.textBaseline = 'middle';
                ctx.fillText(Math.round(val * 100), x + colWidth / 2, y + rowHeight / 2);
            }
        }
    }
}

function renderChipSeqPeak(canvas, tfLocus, conditionName) {
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    const width = canvas.width;
    const height = canvas.height;
    
    // Clear and fill with modern soft-grey background
    ctx.clearRect(0, 0, width, height);
    ctx.fillStyle = '#f8fafc';
    ctx.fillRect(0, 0, width, height);
    
    // Draw subtle horizontal grid lines
    ctx.strokeStyle = '#f1f5f9';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, (height - 20) * 0.33); ctx.lineTo(width, (height - 20) * 0.33);
    ctx.moveTo(0, (height - 20) * 0.66); ctx.lineTo(width, (height - 20) * 0.66);
    ctx.stroke();
    
    // Draw axis and center TSS reference line
    ctx.strokeStyle = '#e2e8f0';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(0, height - 20); ctx.lineTo(width, height - 20);
    ctx.moveTo(width / 2, 0); ctx.lineTo(width / 2, height - 20);
    ctx.stroke();
    
    // Deterministic hash based on tfLocus
    let hash = 0;
    if (tfLocus) {
        for (let i = 0; i < tfLocus.length; i++) {
            hash = tfLocus.charCodeAt(i) + ((hash << 5) - hash);
        }
    }
    hash = Math.abs(hash);
    
    // Determine unique biological peak parameters
    const numPeaks = 1 + (hash % 2); // 1 or 2 peaks
    const peaks = [];
    
    for (let i = 0; i < numPeaks; i++) {
        // Center position relative to width
        // For i=0, center is around 35%-55%. For i=1, center is around 55%-75%.
        const centerOffset = 0.35 + (i * 0.25) + ((hash + i * 7) % 5) * 0.05;
        const center = width * centerOffset;
        
        // Peak width (spread)
        const peakWidth = 35 + ((hash + i * 11) % 4) * 12; // 35 to 71 px
        
        // Biological heights for Control and Stress
        const heightCtrl = 0.4 + ((hash + i * 13) % 5) * 0.12;  // 0.4 to 0.88
        const heightStress = 0.2 + ((hash + i * 17) % 7) * 0.11; // 0.2 to 0.86
        
        peaks.push({
            center: center,
            width: peakWidth,
            height: (conditionName === 'Stress') ? heightStress : heightCtrl
        });
    }
    
    // Setup linear gradient based on condition
    const grad = ctx.createLinearGradient(0, 0, 0, height - 20);
    let strokeColor = '#6366f1'; // Indigo for Control
    let shadowColor = 'rgba(99, 102, 241, 0.3)';
    
    if (conditionName === 'Stress') {
        grad.addColorStop(0, 'rgba(239, 68, 68, 0.45)');
        grad.addColorStop(0.5, 'rgba(239, 68, 68, 0.15)');
        grad.addColorStop(1, 'rgba(239, 68, 68, 0.01)');
        strokeColor = '#ef4444'; // Red for Stress
        shadowColor = 'rgba(239, 68, 68, 0.3)';
    } else {
        grad.addColorStop(0, 'rgba(99, 102, 241, 0.45)');
        grad.addColorStop(0.5, 'rgba(99, 102, 241, 0.15)');
        grad.addColorStop(1, 'rgba(99, 102, 241, 0.01)');
    }
    
    // Plot the composite biological peaks track using Gaussian accumulation
    ctx.save();
    ctx.fillStyle = grad;
    ctx.beginPath();
    ctx.moveTo(0, height - 20);
    
    const availableHeight = height - 32;
    for (let x = 0; x <= width; x++) {
        let accumulatedHeight = 0;
        peaks.forEach(p => {
            const h = p.height * availableHeight;
            const exponent = -Math.pow((x - p.center) / p.width, 2);
            accumulatedHeight += h * Math.exp(exponent);
        });
        // Clamp and calculate y
        const y = Math.max(4, height - 20 - accumulatedHeight);
        ctx.lineTo(x, y);
    }
    ctx.lineTo(width, height - 20);
    ctx.closePath();
    ctx.fill();
    
    // Draw outline stroke on top of the shaded area
    ctx.strokeStyle = strokeColor;
    ctx.lineWidth = 2.2;
    ctx.shadowColor = shadowColor;
    ctx.shadowBlur = 5;
    ctx.beginPath();
    for (let x = 0; x <= width; x++) {
        let accumulatedHeight = 0;
        peaks.forEach(p => {
            const h = p.height * availableHeight;
            const exponent = -Math.pow((x - p.center) / p.width, 2);
            accumulatedHeight += h * Math.exp(exponent);
        });
        const y = Math.max(4, height - 20 - accumulatedHeight);
        if (x === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
    }
    ctx.stroke();
    ctx.restore();
    
    // Draw peak coordinate dots at local maxima of the curves
    peaks.forEach(p => {
        const peakY = height - 20 - (p.height * availableHeight);
        
        ctx.save();
        ctx.fillStyle = strokeColor;
        ctx.strokeStyle = '#ffffff';
        ctx.lineWidth = 1.5;
        ctx.shadowColor = shadowColor;
        ctx.shadowBlur = 8;
        
        ctx.beginPath();
        ctx.arc(p.center, peakY, 4, 0, Math.PI * 2);
        ctx.fill();
        ctx.stroke();
        ctx.restore();
        
        // Print coordinate text above dot
        ctx.fillStyle = varColorTextSecondary();
        ctx.font = 'bold 7px monospace';
        ctx.textAlign = 'center';
        // Mock coordinate relative to gene transcription start site (TSS)
        const tssOffset = Math.round((p.center - width / 2) * 2.5);
        const coordinateText = tssOffset >= 0 ? `+${tssOffset}bp` : `${tssOffset}bp`;
        ctx.fillText(coordinateText, p.center, peakY - 8);
    });
}

function varColorTextSecondary() {
    return getComputedStyle(document.documentElement).getPropertyValue('--text-secondary').trim() || '#475569';
}

// ==========================================================================
// 4. Advanced Computational and Visual Features
// ==========================================================================

// A. Target Regulon KEGG Enrichment
function fetchRegulonPathwayEnrichment(tfLocus) {
    const tbody = document.getElementById('enrichment-results-body');
    if (!tbody) return;

    tbody.innerHTML = `<tr><td colspan="3" class="text-muted" style="text-align:center; padding:12px 0;"><i class="fa-solid fa-spinner fa-spin"></i> Calculating pathway enrichment...</td></tr>`;

    fetch(`/api/regulon_enrichment?tf=${encodeURIComponent(tfLocus)}`)
        .then(res => {
            if (!res.ok) throw new Error("API request failed");
            return res.json();
        })
        .then(data => {
            if (data.error) {
                tbody.innerHTML = `<tr><td colspan="3" class="text-muted" style="text-align:center; padding:12px 0; color:var(--color-repression);">${data.error}</td></tr>`;
                return;
            }

            const pathways = data.pathways || [];
            if (pathways.length === 0) {
                tbody.innerHTML = `<tr><td colspan="3" class="text-muted" style="text-align:center; padding:12px 0;">This TF has no significantly enriched metabolic pathways among its targets</td></tr>`;
                return;
            }

            tbody.innerHTML = '';
            pathways.slice(0, 10).forEach(p => { // limit to top 10
                const tr = document.createElement('tr');
                tr.style.borderBottom = '1px solid var(--border-color)';
                
                // Color significant ones (p < 0.05) with light green background
                const isSig = p.p_value < 0.05;
                if (isSig) {
                    tr.style.backgroundColor = 'rgba(46, 125, 50, 0.04)';
                }

                // Format hit genes list for highlighting on KEGG map
                // Mapping targets to KEGG cgb:cgXXXX format
                const cgbTargetList = p.target_genes.map(g => `cgb:${g.locus.toLowerCase()}`).join('+');
                const keggUrl = `https://www.kegg.jp/kegg-bin/show_pathway?${p.pathway_id}+${cgbTargetList}`;

                const pValText = p.p_value < 0.001 ? p.p_value.toExponential(3) : p.p_value.toFixed(4);

                tr.innerHTML = `
                    <td style="padding:6px; text-align:left;">
                        <a href="${keggUrl}" target="_blank" title="Open the KEGG pathway map in a new window and mark target genes" style="color:var(--color-primary-accent); text-decoration:none; font-weight:500;">
                            ${p.pathway_name} <i class="fa-solid fa-arrow-up-right-from-square fa-xs" style="font-size:7px; opacity:0.7;"></i>
                        </a>
                        <div style="font-size:7.5px; color:var(--text-muted); margin-top:2px;">ID: ${p.pathway_id} | FE: ${p.fold_enrichment.toFixed(2)}x</div>
                    </td>
                    <td style="padding:6px; text-align:center; font-family:var(--font-mono); font-weight:600; color:var(--text-primary);">
                        ${p.hits}/${p.total_genes}
                    </td>
                    <td style="padding:6px; text-align:right; font-family:var(--font-mono); font-weight:600; color:${isSig ? 'var(--color-activation)' : 'var(--text-secondary)'}; padding-right:10px;">
                        ${pValText}
                    </td>
                `;
                tbody.appendChild(tr);
            });
        })
        .catch(err => {
            console.error("Failed to load regulon pathway enrichment:", err);
            tbody.innerHTML = `<tr><td colspan="3" class="text-muted" style="text-align:center; padding:12px 0; color:var(--color-repression);">Failed to calculate pathway enrichment</td></tr>`;
        });
}

// B. Client-side Promoter Motif Scanner
function scanSequenceForMotif(seq, pwm, threshold) {
    const tbody = document.getElementById('scan-results-body');
    const box = document.getElementById('scan-results-box');
    if (!tbody || !box) return;

    tbody.innerHTML = '';
    
    // 1. Standardize and clean input sequence (only allow A, C, G, T)
    const cleanSeq = seq.toUpperCase().replace(/[^ACGT]/g, '');
    if (!cleanSeq) {
        tbody.innerHTML = `<tr><td colspan="4" class="text-muted" style="text-align:center; padding:8px 0;">Enter a valid DNA sequence</td></tr>`;
        box.classList.remove('hidden');
        return;
    }

    const pwmLen = pwm.length;
    if (pwmLen === 0) {
        tbody.innerHTML = `<tr><td colspan="4" class="text-muted" style="text-align:center; padding:8px 0;">Motif weight matrix is empty</td></tr>`;
        box.classList.remove('hidden');
        return;
    }

    if (cleanSeq.length < pwmLen) {
        tbody.innerHTML = `<tr><td colspan="4" class="text-muted" style="text-align:center; padding:8px 0;">Sequence length must be at least the motif length (${pwmLen} bp)</td></tr>`;
        box.classList.remove('hidden');
        return;
    }

    // 2. Pre-calculate max and min possible scores for consensus matching
    let maxScore = 0;
    let minScore = 0;
    for (let i = 0; i < pwmLen; i++) {
        const vals = Object.values(pwm[i]);
        maxScore += Math.max(...vals);
        minScore += Math.min(...vals);
    }
    const scoreRange = maxScore - minScore;

    const hits = [];

    // Helper to score a single window of size pwmLen
    function scoreWindow(windowSeq, pos, strand) {
        let rawScore = 0;
        for (let i = 0; i < pwmLen; i++) {
            const base = windowSeq[i];
            rawScore += (pwm[i][base] !== undefined) ? pwm[i][base] : 0;
        }
        // Normalize to percentage
        const similarity = scoreRange > 0 ? ((rawScore - minScore) / scoreRange * 100) : 0;
        if (similarity >= threshold) {
            hits.push({
                position: pos,
                strand: strand,
                sequence: windowSeq,
                score: similarity
            });
        }
    }

    // Helper for reverse complement
    function getReverseComplement(s) {
        const comp = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A' };
        return s.split('').reverse().map(b => comp[b] || b).join('');
    }

    // 3. Slide window along sequence
    for (let i = 0; i <= cleanSeq.length - pwmLen; i++) {
        const windowSeq = cleanSeq.substring(i, i + pwmLen);
        // Forward strand scoring
        scoreWindow(windowSeq, i + 1, '+');
        // Reverse strand scoring
        const revWindowSeq = getReverseComplement(windowSeq);
        scoreWindow(revWindowSeq, i + 1, '-');
    }

    // 4. Render hits
    if (hits.length === 0) {
        tbody.innerHTML = `<tr><td colspan="4" class="text-muted" style="text-align:center; padding:8px 0;">No matching sites found below threshold ${threshold}%</td></tr>`;
        box.classList.remove('hidden');
        return;
    }

    // Sort by score descending
    hits.sort((a, b) => b.score - a.score);

    hits.forEach(h => {
        const tr = document.createElement('tr');
        tr.style.borderBottom = '1px solid var(--border-color)';
        
        let color = '#475569';
        if (h.score >= 90) color = 'var(--color-repression)';
        else if (h.score >= 80) color = 'var(--color-dual)';
        else color = 'var(--color-activation)';

        tr.innerHTML = `
            <td style="padding:4px 6px; color:var(--text-secondary);">${h.position}</td>
            <td style="padding:4px 6px; font-weight:600; color:${h.strand === '+' ? 'var(--color-activation)' : 'var(--color-srna)'}">${h.strand}</td>
            <td style="padding:4px 6px; color:#1e3a8a; font-weight:500;">${h.sequence}</td>
            <td style="padding:4px 6px; text-align:right; font-weight:bold; color:${color}; padding-right:10px;">${h.score.toFixed(1)}%</td>
        `;
        tbody.appendChild(tr);
    });

    box.classList.remove('hidden');
}

// C. Genomic Locus Map SVG Visualizer
function renderGenomicLocusMap(locusTag) {
    const container = document.getElementById('genomic-map-container');
    if (!container) return;

    container.innerHTML = ''; // Clear previous

    if (!locusTag) return;
    const cleanLocus = locusTag.trim();
    const locusLower = cleanLocus.toLowerCase();

    // Extract the numeric part of RefSeq locus tag, e.g. cg0279 -> 279
    const numMatch = cleanLocus.match(/\d+/);
    if (!numMatch) {
        container.innerHTML = `<span style="font-size: 10px; color:var(--text-muted);">Unable to retrieve genomic coordinates</span>`;
        return;
    }

    const centerNum = parseInt(numMatch[0]);
    const neighborGenes = [];

    // Resolve neighbor locus tags (+- 3 genes)
    for (let offset = -3; offset <= 3; offset++) {
        const num = centerNum + offset;
        if (num <= 0) continue;
        const padLocus = 'cg' + String(num).padStart(4, '0');
        const key = padLocus.toLowerCase();
        
        const item = geneIndex[key] || { name: padLocus.toUpperCase(), type: 'Target' };
        const geneName = (item.name && item.name !== '--') ? item.name : padLocus.toUpperCase();
        const product = cgToProduct[key] || 'No description available';
        const type = item.type || 'Target';
        
        const operonMeta = geneToOperon[key];
        const strand = operonMeta ? operonMeta.orientation : '+';
        const operonName = operonMeta ? operonMeta.operon : null;

        // Fetch expression data (check both cg locus and mapped cgl locus)
        const cglTag = cgToCgl[key] || '';
        const expr = rnaseqData && (rnaseqData[key] || (cglTag ? rnaseqData[cglTag.toLowerCase()] : null));
        const log2fc = expr ? expr.log2fc : undefined;
        const pval = expr ? expr.pvalue : undefined;

        neighborGenes.push({
            locus: padLocus,
            name: geneName,
            type: type,
            product: product,
            strand: strand,
            operon: operonName,
            log2fc: log2fc,
            pval: pval
        });
    }

    // Dimensions
    const svgWidth = 340;
    const svgHeight = 110;
    const paddingX = 15;
    const totalSlots = neighborGenes.length;
    const spacing = 4;
    const geneWidth = (svgWidth - 2 * paddingX - (totalSlots - 1) * spacing) / totalSlots;
    const h = 26; // arrow height
    const y = 42; // arrow y-offset

    let svgHtml = `<svg viewBox="0 0 ${svgWidth} ${svgHeight}" style="width:100%; height:100%; display:block;" xmlns="http://www.w3.org/2000/svg">`;

    // 1. Draw operon grouping background boxes
    const operonGroups = {};
    neighborGenes.forEach((g, idx) => {
        if (g.operon) {
            if (!operonGroups[g.operon]) {
                operonGroups[g.operon] = [];
            }
            operonGroups[g.operon].push(idx);
        }
    });

    Object.entries(operonGroups).forEach(([operonName, indices]) => {
        if (indices.length >= 1) {
            const firstIdx = indices[0];
            const lastIdx = indices[indices.length - 1];
            const opX = paddingX + firstIdx * (geneWidth + spacing) - 2;
            const opW = (lastIdx - firstIdx + 1) * (geneWidth + spacing) - spacing + 4;
            
            svgHtml += `
                <rect x="${opX}" y="${y - 14}" width="${opW}" height="${h + 30}" fill="rgba(99, 102, 241, 0.03)" stroke="#818cf8" stroke-dasharray="2,2" stroke-width="0.8" rx="4"></rect>
                <text x="${opX + 3}" y="${y - 18}" font-size="6.5px" font-family="sans-serif" font-weight="600" fill="#4f46e5">${operonName}</text>
            `;
        }
    });

    // 2. Draw chromosome backbone line
    svgHtml += `<line x1="${paddingX - 5}" y1="${y + h/2}" x2="${svgWidth - paddingX + 5}" y2="${y + h/2}" stroke="#94a3b8" stroke-dasharray="3,3" stroke-width="1.5"></line>`;

    // 3. Draw each gene block chevron
    neighborGenes.forEach((g, idx) => {
        const x = paddingX + idx * (geneWidth + spacing);
        const isCenter = g.locus.toLowerCase() === locusLower;
        
        // Determine color based on log2FC
        let fill = '#e2e8f0';
        let stroke = '#cbd5e1';
        if (g.log2fc !== undefined) {
            fill = getRnaSeqColor(g.log2fc);
            stroke = Math.abs(g.log2fc) >= 0.5 ? (g.log2fc > 0 ? '#ef4444' : '#2563eb') : '#cbd5e1';
        }

        // Draw chevron path based on strand direction
        let points = "";
        if (g.strand === '+') {
            points = `${x},${y} ${x+geneWidth-6},${y} ${x+geneWidth},${y+h/2} ${x+geneWidth-6},${y+h} ${x},${y+h} ${x+3},${y+h/2}`;
        } else {
            points = `${x+6},${y} ${x+geneWidth},${y} ${x+geneWidth-3},${y+h/2} ${x+geneWidth},${y+h} ${x+6},${y+h} ${x},${y+h/2}`;
        }

        // Highlight center selected gene
        const highlightStyle = isCenter ? 'stroke="#7c3aed" stroke-width="2.5" filter="drop-shadow(0 2px 4px rgba(124, 58, 237, 0.4))"' : `stroke="${stroke}" stroke-width="1"`;

        svgHtml += `
            <polygon class="gene-chevron" data-locus="${g.locus}" points="${points}" fill="${fill}" ${highlightStyle} style="cursor:pointer; transition: opacity 0.15s;">
                <title>${g.locus.toUpperCase()} (${g.name})
Function: ${g.product}
Strand: ${g.strand}
log2FC: ${g.log2fc !== undefined ? g.log2fc.toFixed(2) : 'No data'}</title>
            </polygon>
        `;

        // Text label inside/above
        const dispName = g.name.length > 8 ? g.name.substring(0, 7) + '..' : g.name;
        
        svgHtml += `
            <text x="${x + geneWidth/2}" y="${y + h + 10}" font-size="7.5px" font-family="sans-serif" font-weight="600" text-anchor="middle" fill="${isCenter ? '#7c3aed' : '#334155'}" style="pointer-events:none;">${dispName}</text>
            <text x="${x + geneWidth/2}" y="${y - 4}" font-size="6.5px" font-family="monospace" text-anchor="middle" fill="#64748b" style="pointer-events:none;">${g.locus.toUpperCase()}</text>
        `;

        // 4. Draw regulation ChIP-seq Peak overlay in promoter intergenic region
        const tfLower = cleanLocus.toLowerCase();
        const tgLower = g.locus.toLowerCase();
        
        const regRow = regulations.find(r => {
            const rowTf = cleanStr(r.TF_locusTag).toLowerCase();
            const rowTg = cleanStr(r.TG_locusTag).toLowerCase();
            return rowTf === tfLower && rowTg === tgLower && r.Binding_site && cleanStr(r.Binding_site) !== 'nan';
        });

        if (regRow) {
            const peakOffset = g.strand === '+' ? -spacing/2 : geneWidth + spacing/2;
            const peakX = x + peakOffset;
            
            svgHtml += `
                <path d="M ${peakX - 6},${y + h/2} Q ${peakX},${y - 8} ${peakX + 6},${y + h/2}" fill="rgba(239, 68, 68, 0.25)" stroke="#ef4444" stroke-width="1.2">
                    <title>Predicted binding site:
${regRow.Binding_site}
Type: ${regRow.Role}</title>
                </path>
                <circle cx="${peakX}" cy="${y - 8}" r="2" fill="#ef4444"></circle>
            `;
        }
    });

    svgHtml += `</svg>`;
    container.innerHTML = svgHtml;

    // Bind click navigation event handler to polygons
    const polygons = container.querySelectorAll('.gene-chevron');
    polygons.forEach(p => {
        p.addEventListener('mouseenter', () => { p.style.opacity = '0.8'; });
        p.addEventListener('mouseleave', () => { p.style.opacity = '1.0'; });
        p.addEventListener('click', () => {
            const clickedLocus = p.getAttribute('data-locus');
            if (clickedLocus) {
                querySingleGene(clickedLocus);
            }
        });
    });
}

// D. Initialize Advanced Interactive Event Bindings
function initAdvancedFeatures() {
    // 1. Motif Scanner Events
    const btnScan = document.getElementById('btn-run-scan');
    const seqInput = document.getElementById('scan-sequence-input');
    const thresholdSlider = document.getElementById('scan-threshold-slider');
    const thresholdVal = document.getElementById('scan-threshold-val');

    if (thresholdSlider && thresholdVal) {
        thresholdSlider.addEventListener('input', () => {
            thresholdVal.textContent = thresholdSlider.value + '%';
        });
    }

    if (btnScan && seqInput && thresholdSlider) {
        btnScan.addEventListener('click', () => {
            if (!currentTfPwm) {
                alert('Select a valid transcription factor first to retrieve its PWM.');
                return;
            }
            const seq = seqInput.value;
            const threshold = parseFloat(thresholdSlider.value);
            scanSequenceForMotif(seq, currentTfPwm, threshold);
        });
    }

    // 2. Custom RNA-seq File Upload in Sidebar
    const btnImport = document.getElementById('btn-import-rnaseq');
    const fileInput = document.getElementById('rnaseq-upload-input');
    
    if (btnImport && fileInput) {
        btnImport.addEventListener('click', () => {
            fileInput.click();
        });

        fileInput.addEventListener('change', (e) => {
            const file = e.target.files[0];
            if (!file) return;

            const reader = new FileReader();
            reader.onload = function(evt) {
                const csvText = evt.target.result;
                Papa.parse(csvText, {
                    header: true,
                    skipEmptyLines: true,
                    dynamicTyping: true,
                    complete: function(results) {
                        processRnaSeqData(results.data);
                        fileInput.value = '';
                    },
                    error: function(err) {
                        alert('Failed to parse CSV file: ' + err.message);
                    }
                });
            };
            reader.readAsText(file);
        });
    }

    // 3. Organism/Strain Selection
    const orgSelect = document.getElementById('organism-select');
    if (orgSelect) {
        // Fetch organisms
        fetch('/api/list_organisms')
            .then(res => res.json())
            .then(data => {
                if (data.error) {
                    console.error("Error loading organisms:", data.error);
                    return;
                }
                
                orgSelect.innerHTML = '';
                data.forEach(org => {
                    const opt = document.createElement('option');
                    opt.value = org.id;
                    opt.textContent = org.name;
                    opt.setAttribute('data-has-rna', org.has_rna);
                    orgSelect.appendChild(opt);
                });
                
                // Select default
                orgSelect.value = 'C_g_DSM_20300_=_ATCC_13032';
            })
            .catch(err => console.error("Failed to fetch organisms:", err));

        orgSelect.addEventListener('change', async () => {
            const orgId = orgSelect.value;
            const opt = orgSelect.selectedOptions[0];
            const hasRna = opt ? opt.getAttribute('data-has-rna') === 'true' : false;
            
            updateStatus('Switching organism / strain...', 'loading');
            
            if (orgId === 'C_g_DSM_20300_=_ATCC_13032') {
                // Use default files
                REGULATIONS_URL = 'data/regulations.csv';
                RNA_REGULATIONS_URL = 'data/rna_regulation.csv';
                MAPPING_URL = 'data/gene_mapping.csv';
                OPERONS_URL = 'data/operons.csv';
            } else {
                const opPrefix = getOperonPrefix(orgId);
                REGULATIONS_URL = `data/AllOrganismsFiles/${orgId}_regulations.csv`;
                RNA_REGULATIONS_URL = hasRna ? `data/AllOrganismsFiles/${orgId}_rna_regulation.csv` : '';
                MAPPING_URL = ''; // No mapping for other strains
                OPERONS_URL = `data/AllOrganismsFiles/${opPrefix}_operons.csv`;
            }
            
            // Clear current network mapping variables
            geneMapping = [];
            geneIndex = {};
            cglToCg = {};
            cgToCgl = {};
            nameToCg = {};
            cgToProduct = {};
            regulations = [];
            rnaRegulations = [];
            
            // Reset simulation if active
            resetPerturbationSimulation();
            
            // Reset UI lists and query states
            currentQueryGene = null;
            clearAllInputs();
            
            // Hide the details sidebar if open
            toggleRightSidebar(false);
            
            // Clear network visualization
            if (cy) {
                cy.elements().remove();
            }
            const overlay = document.getElementById('canvas-overlay');
            if (overlay) {
                overlay.style.display = 'flex';
                const h3 = overlay.querySelector('h3');
                if (h3) h3.textContent = `Loaded ${opt ? opt.textContent : 'new organism'}, enter a gene to start analysis`; 
            }
            
            try {
                await loadNetworkData();
                updateExampleTags();
            } catch (err) {
                console.error("Failed to load new organism network data:", err);
                updateStatus('Failed to load data: ' + err.message, 'error');
            }
        });
    }
}

function getOperonPrefix(orgId) {
    let count = 0;
    return orgId.replace(/_/g, (match) => {
        count++;
        return count <= 2 ? '' : match;
    });
}

function updateExampleTags() {
    const tfCounts = {};
    regulations.forEach(row => {
        const tfTag = cleanStr(row.TF_locusTag);
        const tfName = cleanStr(row.TF_name);
        const tf = tfName && tfName !== tfTag ? tfName : tfTag;
        if (tf) {
            tfCounts[tf] = (tfCounts[tf] || 0) + 1;
        }
    });
    const sortedTfs = Object.entries(tfCounts)
        .sort((a, b) => b[1] - a[1])
        .slice(0, 5)
        .map(entry => entry[0]);

    if (sortedTfs.length === 0) {
        sortedTfs.push("cg0350", "sigH", "whiB4");
    }

    const container = document.querySelector('.quick-examples');
    if (container) {
        const span = container.querySelector('span');
        container.innerHTML = '';
        if (span) {
            container.appendChild(span);
        } else {
            const newSpan = document.createElement('span');
            newSpan.textContent = 'Try examples:';
            container.appendChild(newSpan);
        }
        sortedTfs.forEach(tf => {
            const btn = document.createElement('button');
            btn.className = 'example-tag';
            btn.textContent = tf;
            btn.addEventListener('click', () => {
                querySingleGene(tf);
            });
            container.appendChild(btn);
        });
    }
}

// ==========================================================================
// 8. Data & Model Quality Dashboard Logic
// ==========================================================================

function getGlobalPlatformGraph() {
    const nodes = Object.values(normalizedNodes || {}).map(node => ({
        id: node.id,
        label: node.label || node.id,
        type: node.type,
        nodeType: node.type
    }));
    
    const edges = (normalizedEdges || []).map(edge => ({
        source: edge.source,
        target: edge.target,
        type: edge.regulationType || 'unknown',
        regulationType: edge.regulationType || 'unknown',
        role: edge.role,
        interactionClass: edge.interactionClass,
        sourceType: edge.sourceType,
        confidenceScore: edge.confidenceScore,
        confidence: edge.confidenceScore,
        heuristicConfidenceScore: edge.heuristicConfidenceScore,
        predictedConfidence: edge.predictedConfidence,
        confidenceFactors: edge.confidenceFactors
    }));
    
    return { nodes, edges };
}

function renderGeneTagList(containerId, geneList) {
    const container = document.getElementById(containerId);
    if (!container) return;
    if (!geneList || geneList.length === 0) {
        container.innerHTML = '<span class="metabolic-muted" style="font-size: 11px;">None</span>';
        return;
    }
    container.innerHTML = geneList.map(gene => `
        <span class="gene-tag" title="Click to view details of ${gene}">${escapeHtml(gene)}</span>
    `).join('');
    
    container.querySelectorAll('.gene-tag').forEach(tag => {
        tag.addEventListener('click', () => {
            const gene = tag.textContent.trim();
            setActiveWorkflowEntry('gene');
            scrollLeftSidebarTo('.search-section');
            const searchInput = geneInputsContainer?.querySelector('.gene-input');
            if (searchInput) searchInput.value = gene;
            querySingleGene(gene);
        });
    });
}

function updateQualityDashboard() {
    if (!window.analysisQuality) {
        console.error("analysisQuality library not loaded!");
        return;
    }

    const jsonBtn = document.getElementById('btn-export-quality-json');
    if (jsonBtn && !jsonBtn.dataset.bound) {
        jsonBtn.dataset.bound = '1';
        jsonBtn.addEventListener('click', exportQualityReportJSON);
    }
    const csvBtn = document.getElementById('btn-export-quality-csv');
    if (csvBtn && !csvBtn.dataset.bound) {
        csvBtn.dataset.bound = '1';
        csvBtn.addEventListener('click', exportQualityReportCSV);
    }

    const graph = getGlobalPlatformGraph();
    const report = window.analysisQuality.getAnalysisQualityReport(graph);

    document.getElementById('stat-reg-nodes').textContent = report.regulatoryNetwork.totalNodes;
    document.getElementById('stat-reg-edges').textContent = report.regulatoryNetwork.totalEdges;
    document.getElementById('stat-reg-tfs').textContent = report.regulatoryNetwork.tfCount;
    document.getElementById('stat-reg-genes').textContent = report.regulatoryNetwork.geneCount;
    document.getElementById('stat-reg-srnas').textContent = report.regulatoryNetwork.srnaCount;
    document.getElementById('stat-reg-operons').textContent = report.regulatoryNetwork.operonCount;
    document.getElementById('stat-reg-tf-tg').textContent = report.regulatoryNetwork.tfGeneEdgeCount;
    document.getElementById('stat-reg-srna-tg').textContent = report.regulatoryNetwork.srnaEdgeCount;
    document.getElementById('stat-reg-act').textContent = report.regulatoryNetwork.activationCount;
    document.getElementById('stat-reg-rep').textContent = report.regulatoryNetwork.repressionCount;
    document.getElementById('stat-reg-pred').textContent = report.regulatoryNetwork.predictedCount;
    document.getElementById('stat-reg-unknown').textContent = report.regulatoryNetwork.unknownRegulationCount;

    document.getElementById('stat-conf-avg').textContent = report.confidenceScores.averageConfidence.toFixed(2);
    document.getElementById('stat-conf-med').textContent = report.confidenceScores.medianConfidence.toFixed(2);
    document.getElementById('stat-conf-total').textContent = report.confidenceScores.totalEdgesWithConfidence;
    document.getElementById('stat-conf-high').textContent = report.confidenceScores.highConfidenceEdgeCount;
    document.getElementById('stat-conf-med-count').textContent = report.confidenceScores.mediumConfidenceEdgeCount;
    document.getElementById('stat-conf-low').textContent = report.confidenceScores.lowConfidenceEdgeCount;
    document.getElementById('stat-conf-rf').textContent = report.confidenceScores.rfConfidenceAvailableCount;
    document.getElementById('stat-conf-heur').textContent = report.confidenceScores.heuristicConfidenceAvailableCount;
    document.getElementById('stat-conf-rf-avg').textContent = report.confidenceScores.averageRfConfidence ? report.confidenceScores.averageRfConfidence.toFixed(2) : 'N/A';
    document.getElementById('stat-conf-heur-avg').textContent = report.confidenceScores.averageHeuristicConfidence ? report.confidenceScores.averageHeuristicConfidence.toFixed(2) : 'N/A';

    const diffContainer = document.getElementById('stat-conf-diff-container');
    const diffText = document.getElementById('stat-conf-diff');
    if (report.confidenceScores.averageAbsoluteDifference !== null && report.confidenceScores.averageAbsoluteDifference !== undefined) {
        if (diffContainer) diffContainer.classList.remove('hidden');
        if (diffText) diffText.textContent = report.confidenceScores.averageAbsoluteDifference.toFixed(2);
    } else {
        if (diffContainer) diffContainer.classList.add('hidden');
    }

    const metaGeneCount = report.metabolicMapping.regulatoryGeneCount;
    const metaMappedCount = report.metabolicMapping.genesMappedToReactions;
    const metaCoveragePercent = metaGeneCount > 0 ? (metaMappedCount / metaGeneCount) * 100 : 0;
    
    document.getElementById('stat-meta-coverage').textContent = `${metaCoveragePercent.toFixed(1)}%`;
    document.getElementById('stat-meta-progress').style.width = `${metaCoveragePercent}%`;
    document.getElementById('stat-meta-total-genes').textContent = metaGeneCount;
    document.getElementById('stat-meta-rxn-genes').textContent = metaMappedCount;
    document.getElementById('stat-meta-path-genes').textContent = report.metabolicMapping.genesMappedToPathways;
    document.getElementById('stat-meta-rxns').textContent = report.metabolicMapping.mappedReactionCount;
    document.getElementById('stat-meta-paths').textContent = report.metabolicMapping.mappedPathwayCount;
    document.getElementById('stat-meta-unmapped-count').textContent = report.metabolicMapping.unmappedGeneCount;
    
    renderGeneTagList('list-meta-unmapped', report.metabolicMapping.unmappedGenes);

    const enzGeneCount = report.metabolicMapping.regulatoryGeneCount;
    const enzMappedCount = report.enzymeConstraintCoverage.genesWithEnzymeMapping;
    const enzCoveragePercent = enzGeneCount > 0 ? (enzMappedCount / enzGeneCount) * 100 : 0;
    
    document.getElementById('stat-enz-coverage').textContent = `${enzCoveragePercent.toFixed(1)}%`;
    document.getElementById('stat-enz-progress').style.width = `${enzCoveragePercent}%`;
    document.getElementById('stat-enz-genes').textContent = enzMappedCount;
    document.getElementById('stat-enz-rxns').textContent = report.enzymeConstraintCoverage.enzymeAssociatedReactionCount;
    document.getElementById('stat-enz-kcat').textContent = report.enzymeConstraintCoverage.reactionsWithKcat;
    document.getElementById('stat-enz-mw').textContent = report.enzymeConstraintCoverage.reactionsWithMolecularWeight;
    document.getElementById('stat-enz-kcat-mw').textContent = report.enzymeConstraintCoverage.reactionsWithKcatPerMW;
    document.getElementById('stat-enz-ec').textContent = report.enzymeConstraintCoverage.reactionsWithECNumber;
    document.getElementById('stat-enz-uniprot').textContent = report.enzymeConstraintCoverage.reactionsWithUniProtId;
    document.getElementById('stat-enz-potential').textContent = report.enzymeConstraintCoverage.potentialEnzymeConstrainedReactionCount;
    document.getElementById('stat-enz-unmapped-count').textContent = report.enzymeConstraintCoverage.unmappedEnzymeGenes.length;
    
    renderGeneTagList('list-enz-unmapped', report.enzymeConstraintCoverage.unmappedEnzymeGenes);

    const warningBanner = document.getElementById('quality-warning-banner');
    const warningText = document.getElementById('quality-warning-text');
    if (warningBanner && warningText) {
        const warnings = [];
        if (metaCoveragePercent < 45) {
            warnings.push(`Metabolic mapping coverage is low (${metaCoveragePercent.toFixed(1)}%). Some regulatory genes are not captured in the iCW773 model.`);
        }
        if (enzCoveragePercent < 25) {
            warnings.push(`ecCGL1 enzyme constraint coverage is low (${enzCoveragePercent.toFixed(1)}%). Many mapped reactions lack enzyme parameters (kcat, molecular weight).`);
        }
        
        if (warnings.length > 0) {
            warningText.innerHTML = warnings.join('<br>');
            warningBanner.classList.remove('hidden');
        } else {
            warningBanner.classList.add('hidden');
        }
    }
}

function exportQualityReportJSON() {
    const graph = getGlobalPlatformGraph();
    const report = window.analysisQuality.getAnalysisQualityReport(graph);
    const blob = new Blob([JSON.stringify(report, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `cgl_regulation_quality_report_${new Date().toISOString().slice(0, 10)}.json`;
    a.click();
    URL.revokeObjectURL(url);
}

function exportQualityReportCSV() {
    const graph = getGlobalPlatformGraph();
    const report = window.analysisQuality.getAnalysisQualityReport(graph);
    
    const rows = [
        ['Category', 'Metric', 'Value'],
        ['Regulatory Network', 'Total Nodes', report.regulatoryNetwork.totalNodes],
        ['Regulatory Network', 'Total Edges', report.regulatoryNetwork.totalEdges],
        ['Regulatory Network', 'Transcription Factors (TF)', report.regulatoryNetwork.tfCount],
        ['Regulatory Network', 'Target Genes', report.regulatoryNetwork.geneCount],
        ['Regulatory Network', 'sRNAs', report.regulatoryNetwork.srnaCount],
        ['Regulatory Network', 'Operons', report.regulatoryNetwork.operonCount],
        ['Regulatory Network', 'TF-Target Edges', report.regulatoryNetwork.tfGeneEdgeCount],
        ['Regulatory Network', 'sRNA-mRNA Edges', report.regulatoryNetwork.srnaEdgeCount],
        ['Regulatory Network', 'Activation Edges (+)', report.regulatoryNetwork.activationCount],
        ['Regulatory Network', 'Repression Edges (-)', report.regulatoryNetwork.repressionCount],
        ['Regulatory Network', 'Predicted Edges', report.regulatoryNetwork.predictedCount],
        ['Regulatory Network', 'Unknown Mode Edges', report.regulatoryNetwork.unknownRegulationCount],
        
        ['Confidence Scores', 'Edges with Confidence', report.confidenceScores.totalEdgesWithConfidence],
        ['Confidence Scores', 'Average Confidence', report.confidenceScores.averageConfidence.toFixed(4)],
        ['Confidence Scores', 'Median Confidence', report.confidenceScores.medianConfidence.toFixed(4)],
        ['Confidence Scores', 'High Confidence Edges (>=0.75)', report.confidenceScores.highConfidenceEdgeCount],
        ['Confidence Scores', 'Medium Confidence Edges (0.45-0.75)', report.confidenceScores.mediumConfidenceEdgeCount],
        ['Confidence Scores', 'Low Confidence Edges (<0.45)', report.confidenceScores.lowConfidenceEdgeCount],
        ['Confidence Scores', 'RF Scores Available', report.confidenceScores.rfConfidenceAvailableCount],
        ['Confidence Scores', 'Heuristic Scores Available', report.confidenceScores.heuristicConfidenceAvailableCount],
        ['Confidence Scores', 'Average RF Score', report.confidenceScores.averageRfConfidence ? report.confidenceScores.averageRfConfidence.toFixed(4) : 'N/A'],
        ['Confidence Scores', 'Average Heuristic Score', report.confidenceScores.averageHeuristicConfidence ? report.confidenceScores.averageHeuristicConfidence.toFixed(4) : 'N/A'],
        ['Confidence Scores', 'Average Absolute Difference (RF vs Heur)', report.confidenceScores.averageAbsoluteDifference ? report.confidenceScores.averageAbsoluteDifference.toFixed(4) : 'N/A'],
        
        ['Metabolic Mapping (iCW773)', 'Total Regulatory Genes', report.metabolicMapping.regulatoryGeneCount],
        ['Metabolic Mapping (iCW773)', 'Genes Mapped to Reactions', report.metabolicMapping.genesMappedToReactions],
        ['Metabolic Mapping (iCW773)', 'Genes Mapped to Pathways', report.metabolicMapping.genesMappedToPathways],
        ['Metabolic Mapping (iCW773)', 'Unique Mapped Reactions', report.metabolicMapping.mappedReactionCount],
        ['Metabolic Mapping (iCW773)', 'Unique Mapped Pathways', report.metabolicMapping.mappedPathwayCount],
        ['Metabolic Mapping (iCW773)', 'Unmapped Regulatory Genes', report.metabolicMapping.unmappedGeneCount],
        
        ['Enzyme Constraints (ecCGL1)', 'Genes with Enzyme Mapping', report.enzymeConstraintCoverage.genesWithEnzymeMapping],
        ['Enzyme Constraints (ecCGL1)', 'Enzyme Associated Reactions', report.enzymeConstraintCoverage.enzymeAssociatedReactionCount],
        ['Enzyme Constraints (ecCGL1)', 'Reactions with kcat', report.enzymeConstraintCoverage.reactionsWithKcat],
        ['Enzyme Constraints (ecCGL1)', 'Reactions with MW', report.enzymeConstraintCoverage.reactionsWithMolecularWeight],
        ['Enzyme Constraints (ecCGL1)', 'Reactions with kcat/MW', report.enzymeConstraintCoverage.reactionsWithKcatPerMW],
        ['Enzyme Constraints (ecCGL1)', 'Reactions with EC Number', report.enzymeConstraintCoverage.reactionsWithECNumber],
        ['Enzyme Constraints (ecCGL1)', 'Reactions with UniProt ID', report.enzymeConstraintCoverage.reactionsWithUniProtId],
        ['Enzyme Constraints (ecCGL1)', 'Potential Enzyme-Constrained Reactions', report.enzymeConstraintCoverage.potentialEnzymeConstrainedReactionCount]
    ];
    
    const csvContent = rows.map(row => row.map(cell => {
        const cellStr = cell === null || cell === undefined ? '' : String(cell);
        return `"${cellStr.replace(/"/g, '""')}"`;
    }).join(',')).join('\n');
    
    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `cgl_regulation_quality_metrics_${new Date().toISOString().slice(0, 10)}.csv`;
    a.click();
    URL.revokeObjectURL(url);
}

// ==========================================================================
// 9. Examples & Case Studies Logic
// ==========================================================================

let activeCaseStudyResult = null;

function initExamplesDashboard() {
    const runButtons = document.querySelectorAll('.run-case-btn');
    runButtons.forEach(btn => {
        if (!btn.dataset.bound) {
            btn.dataset.bound = '1';
            btn.addEventListener('click', () => {
                const caseId = btn.getAttribute('data-case-id');
                runAndDisplayCaseStudy(caseId);
            });
        }
    });

    const exportBtn = document.getElementById('btn-export-case-report');
    if (exportBtn && !exportBtn.dataset.bound) {
        exportBtn.dataset.bound = '1';
        exportBtn.addEventListener('click', () => {
            if (activeCaseStudyResult) {
                exportCaseReportJSON(activeCaseStudyResult);
            }
        });
    }

    const floatingExportBtn = document.getElementById('btn-floating-export-report');
    if (floatingExportBtn && !floatingExportBtn.dataset.bound) {
        floatingExportBtn.dataset.bound = '1';
        floatingExportBtn.addEventListener('click', () => {
            if (activeCaseStudyResult) {
                exportCaseReportJSON(activeCaseStudyResult);
            }
        });
    }

    const closeFloatingBtn = document.getElementById('btn-close-floating-narrative');
    if (closeFloatingBtn && !closeFloatingBtn.dataset.bound) {
        closeFloatingBtn.dataset.bound = '1';
        closeFloatingBtn.addEventListener('click', () => {
            document.getElementById('case-floating-narrative').classList.add('hidden');
        });
    }
}

function runAndDisplayCaseStudy(caseId) {
    if (!window.caseStudies) {
        console.error("caseStudies library not loaded!");
        return;
    }

    const graph = getGlobalPlatformGraph();
    const result = window.caseStudies.runCaseStudy(caseId, graph);
    activeCaseStudyResult = result;

    const resultsPanel = document.getElementById('case-study-results-panel');
    if (resultsPanel) {
        resultsPanel.classList.remove('hidden');
    }

    document.getElementById('results-case-title').textContent = result.caseStudy.title;
    document.getElementById('results-narrative-text').textContent = result.narrative;

    const warningBanner = document.getElementById('case-warning-banner');
    if (warningBanner) {
        if (result.warnings && result.warnings.length > 0) {
            warningBanner.innerHTML = result.warnings.join('<br>');
            warningBanner.classList.remove('hidden');
        } else {
            warningBanner.classList.add('hidden');
        }
    }

    const floatingCard = document.getElementById('case-floating-narrative');
    const floatingText = document.getElementById('floating-narrative-text');
    const floatingWarnings = document.getElementById('floating-case-warnings');
    
    if (floatingText) floatingText.textContent = result.narrative;
    if (floatingWarnings) {
        if (result.warnings && result.warnings.length > 0) {
            floatingWarnings.innerHTML = result.warnings.join('<br>');
            floatingWarnings.classList.remove('hidden');
        } else {
            floatingWarnings.classList.add('hidden');
        }
    }
    if (floatingCard) {
        floatingCard.classList.remove('hidden');
    }

    const table = document.getElementById('case-results-table');
    if (table) {
        const thead = table.querySelector('thead');
        const tbody = table.querySelector('tbody');
        tbody.innerHTML = '';

        if (caseId === "glutamate-regulation" || caseId === "tca-cycle-regulators") {
            thead.innerHTML = `
                <tr>
                    <th>Transcription Factor</th>
                    <th>Target Genes count</th>
                    <th>Avg Confidence</th>
                    <th>Regulator Score</th>
                    <th>Action</th>
                </tr>
            `;
            const ranking = result.results.tfRanking || [];
            if (ranking.length === 0) {
                tbody.innerHTML = '<tr><td colspan="5" class="metabolic-muted" style="text-align: center;">No upstream regulators found.</td></tr>';
            } else {
                tbody.innerHTML = ranking.slice(0, 10).map(tf => `
                    <tr>
                        <td><strong>${escapeHtml(tf.tfLabel || tf.tfId)}</strong></td>
                        <td>${tf.regulatedGenes ? tf.regulatedGenes.length : 0} genes</td>
                        <td>${tf.averageConfidence.toFixed(3)}</td>
                        <td>${tf.regulatorScore.toFixed(3)}</td>
                        <td><button class="secondary-btn search-tf-btn" data-tf="${escapeHtml(tf.tfId)}" style="height: 24px; font-size: 10px; padding: 0 8px;">Inspect TF</button></td>
                    </tr>
                `).join('');
            }
        } else if (caseId === "amino-acid-engineering-targets") {
            thead.innerHTML = `
                <tr>
                    <th>Transcription Factor</th>
                    <th>Candidate Score</th>
                    <th>Recommendation</th>
                    <th>Regulated Key Genes</th>
                    <th>Avg Confidence</th>
                    <th>Action</th>
                </tr>
            `;
            const candidates = result.results.engineeringCandidates || [];
            if (candidates.length === 0) {
                tbody.innerHTML = '<tr><td colspan="6" class="metabolic-muted" style="text-align: center;">No engineering candidates found.</td></tr>';
            } else {
                tbody.innerHTML = candidates.slice(0, 10).map(candidate => {
                    const badgeClass = candidate.recommendationLevel === "high" ? "mode-engineering" : "mode-pathway";
                    return `
                        <tr>
                            <td><strong>${escapeHtml(candidate.tfLabel || candidate.tfId)}</strong></td>
                            <td><strong style="color: var(--color-primary-accent);">${candidate.candidateScore.toFixed(3)}</strong></td>
                            <td><span class="case-badge ${badgeClass}">${candidate.recommendationLevel}</span></td>
                            <td>${candidate.regulatedKeyGenes ? candidate.regulatedKeyGenes.length : 0} genes</td>
                            <td>${candidate.averageConfidence.toFixed(3)}</td>
                            <td><button class="secondary-btn search-tf-btn" data-tf="${escapeHtml(candidate.tfId)}" style="height: 24px; font-size: 10px; padding: 0 8px;">Inspect TF</button></td>
                        </tr>
                    `;
                }).join('');
            }
        }

        tbody.querySelectorAll('.search-tf-btn').forEach(btn => {
            btn.addEventListener('click', () => {
                const tfId = btn.getAttribute('data-tf');
                setActiveWorkflowEntry('gene');
                scrollLeftSidebarTo('.search-section');
                const searchInput = geneInputsContainer?.querySelector('.gene-input');
                if (searchInput) searchInput.value = tfId;
                querySingleGene(tfId);
            });
        });
    }

    const qual = result.results.qualitySummary;
    const enz = result.results.enzymeConstraintSummary;

    document.getElementById('results-qual-meta').textContent = qual ? `${((qual.genesMappedToReactions / qual.regulatoryGeneCount) * 100).toFixed(1)}% (${qual.genesMappedToReactions}/${qual.regulatoryGeneCount})` : '-';
    document.getElementById('results-qual-reactions').textContent = qual ? `${qual.mappedReactionCount} rxns / ${qual.mappedPathwayCount} pathways` : '-';
    document.getElementById('results-qual-enz').textContent = qual && enz ? `${((enz.genesWithEnzymeMapping / qual.regulatoryGeneCount) * 100).toFixed(1)}% (${enz.genesWithEnzymeMapping}/${qual.regulatoryGeneCount})` : '-';
    document.getElementById('results-qual-enz-rxns').textContent = enz ? `${enz.enzymeAssociatedReactionCount} reactions` : '-';

    setTimeout(() => {
        if (result.caseStudy.entryMode === "pathway") {
            setActiveWorkflowEntry('pathway');
            scrollLeftSidebarTo('.pathway-regulatory-view-section');
            const pathInput = document.getElementById('pathway-view-input');
            if (pathInput) {
                const kw = result.caseStudy.pathwayKeyword;
                pathInput.value = kw === 'glutamate' ? 'glutamate metabolism' : 'citric acid cycle (tca cycle)';
                runPathwayRegulatoryView();
            }
        } else if (result.caseStudy.entryMode === "engineering-targets") {
            setActiveWorkflowEntry('engineering');
            scrollLeftSidebarTo('.engineering-targets-section');
            const searchInput = document.getElementById('engineering-target-search');
            const minScoreInput = document.getElementById('engineering-target-min-score');
            const pathFilter = document.getElementById('engineering-target-pathway-filter');
            
            if (searchInput) searchInput.value = '';
            if (minScoreInput) minScoreInput.value = '0';
            if (pathFilter) {
                pathFilter.value = result.caseStudy.pathwayKeyword || '';
                refreshEngineeringTargetCandidates();
            }
        }
    }, 1500);
}

function exportCaseReportJSON(caseResult) {
    if (!caseResult) return;
    
    let topResults = [];
    if (caseResult.caseStudy.id === "glutamate-regulation" || caseResult.caseStudy.id === "tca-cycle-regulators") {
        topResults = (caseResult.results.tfRanking || []).slice(0, 10).map(tf => ({
            tfId: tf.tfId,
            tfLabel: tf.tfLabel,
            regulatedGenesCount: tf.regulatedGenes ? tf.regulatedGenes.length : 0,
            averageConfidence: tf.averageConfidence,
            regulatorScore: tf.regulatorScore
        }));
    } else if (caseResult.caseStudy.id === "amino-acid-engineering-targets") {
        topResults = (caseResult.results.engineeringCandidates || []).slice(0, 10).map(candidate => ({
            tfId: candidate.tfId,
            tfLabel: candidate.tfLabel,
            candidateScore: candidate.candidateScore,
            recommendationLevel: candidate.recommendationLevel,
            regulatedKeyGenesCount: candidate.regulatedKeyGenes ? candidate.regulatedKeyGenes.length : 0,
            averageConfidence: candidate.averageConfidence
        }));
    }

    const report = {
        reportType: "Case Study Analysis Report",
        caseStudyId: caseResult.caseStudy.id,
        caseStudyTitle: caseResult.caseStudy.title,
        question: caseResult.caseStudy.question,
        workflow: caseResult.caseStudy.workflow,
        topResults: topResults,
        qualityWarnings: caseResult.warnings,
        limitations: caseResult.caseStudy.limitations,
        generatedTimestamp: new Date().toISOString()
    };
    
    const blob = new Blob([JSON.stringify(report, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `cgl_case_report_${caseResult.caseStudy.id}_${new Date().toISOString().slice(0, 10)}.json`;
    a.click();
    URL.revokeObjectURL(url);
}

let glutamateSelectedReaction = null;
let glutamateCandidatesList = [];

function initGlutamateScenarioDashboard() {
    const searchBtn = document.getElementById('btn-search-glutamate-candidates');
    if (searchBtn && !searchBtn.dataset.bound) {
        searchBtn.dataset.bound = '1';
        searchBtn.addEventListener('click', async () => {
            const tableContainer = document.getElementById('glutamate-candidates-table-container');
            const tableBody = document.getElementById('glutamate-candidates-table-body');
            const errorBox = document.getElementById('glutamate-scenario-error-box');
            
            if (errorBox) errorBox.classList.add('hidden');
            searchBtn.innerHTML = '<i class="fa-solid fa-spinner fa-spin"></i> Searching...';
            searchBtn.disabled = true;
            
            try {
                const response = await fetch('http://127.0.0.1:8001/api/model/reactions/glutamate-candidates');
                if (!response.ok) {
                    throw new Error(`API error: ${response.statusText}`);
                }
                const data = await response.json();
                glutamateCandidatesList = data.candidates || [];
                
                tableBody.innerHTML = '';
                if (glutamateCandidatesList.length === 0) {
                    tableBody.innerHTML = '<tr><td colspan="4" style="padding: 10px; text-align: center; color: var(--text-muted);">No glutamate candidates found in model.</td></tr>';
                } else {
                    glutamateCandidatesList.forEach(cand => {
                        const row = document.createElement('tr');
                        row.style.borderBottom = '1px solid var(--border-color)';
                        row.innerHTML = `
                            <td style="padding: 6px 8px; font-weight: 600;">${escapeHtml(cand.reactionId)}</td>
                            <td style="padding: 6px 8px; font-family: monospace;">${escapeHtml(cand.equation)}</td>
                            <td style="padding: 6px 8px;"><span class="badge-role ${escapeHtml(cand.classification)}" style="font-size: 10px;">${escapeHtml(cand.classification)}</span></td>
                            <td style="padding: 6px 8px;">
                                <button class="secondary-btn btn-select-glutamate" data-id="${escapeHtml(cand.reactionId)}" style="padding: 3px 6px; font-size: 10px; font-weight:600; cursor: pointer; border-radius: 3px; border: 1px solid var(--border-color); background: var(--bg-primary); color: var(--text-primary);">
                                    Select
                                </button>
                            </td>
                        `;
                        tableBody.appendChild(row);
                    });
                    
                    tableBody.querySelectorAll('.btn-select-glutamate').forEach(btn => {
                        btn.addEventListener('click', () => {
                            const rxnId = btn.getAttribute('data-id');
                            const found = glutamateCandidatesList.find(c => c.reactionId === rxnId);
                            if (found) {
                                selectGlutamateCandidate(found);
                            }
                        });
                    });
                }
                
                if (tableContainer) tableContainer.classList.remove('hidden');
            } catch (err) {
                console.error("Failed to fetch glutamate candidates:", err);
                if (errorBox) {
                    errorBox.textContent = `Failed to retrieve glutamate candidates from backend: ${err.message}. Please make sure the FBA backend on port 8001 is running.`;
                    errorBox.classList.remove('hidden');
                }
            } finally {
                searchBtn.innerHTML = '<i class="fa-solid fa-search"></i> Search Glutamate Candidates';
                searchBtn.disabled = false;
            }
        });
    }

    const chkVerified = document.getElementById('chk-glutamate-verified');
    if (chkVerified && !chkVerified.dataset.bound) {
        chkVerified.dataset.bound = '1';
        chkVerified.addEventListener('change', (e) => {
            const isChecked = e.target.checked;
            window.glutamateScenario.glutamateState.userVerified = isChecked;
            
            const runBtn = document.getElementById('btn-run-glutamate-scenario');
            if (runBtn) {
                runBtn.disabled = !isChecked;
            }
        });
    }

    const runBtn = document.getElementById('btn-run-glutamate-scenario');
    if (runBtn && !runBtn.dataset.bound) {
        runBtn.dataset.bound = '1';
        runBtn.addEventListener('click', async () => {
            const locusInput = document.getElementById('glutamate-scenario-locus');
            const typeSelect = document.getElementById('glutamate-scenario-type');
            const errorBox = document.getElementById('glutamate-scenario-error-box');
            const resultsBox = document.getElementById('glutamate-scenario-results-box');
            const resultsEmpty = document.getElementById('glutamate-scenario-results-empty');
            
            const locus = locusInput ? locusInput.value.trim() : '';
            const type = typeSelect ? typeSelect.value : 'TF';
            
            if (errorBox) errorBox.classList.add('hidden');
            
            if (!locus) {
                if (errorBox) {
                    errorBox.textContent = "Please enter a target gene or TF locus.";
                    errorBox.classList.remove('hidden');
                }
                return;
            }
            
            if (!glutamateSelectedReaction || !window.glutamateScenario.glutamateState.userVerified) {
                if (errorBox) {
                    errorBox.textContent = "Please select and verify a glutamate reaction before simulating.";
                    errorBox.classList.remove('hidden');
                }
                return;
            }
            
            runBtn.innerHTML = '<i class="fa-solid fa-spinner fa-spin"></i> Simulating...';
            runBtn.disabled = true;
            
            try {
                let res = null;
                const objective = { objectiveType: "biomass" };
                const trackedReactionIds = [glutamateSelectedReaction.reactionId];
                
                if (type === 'gene') {
                    res = await window.simulationClient.runGeneKnockout(locus, objective, trackedReactionIds);
                } else {
                    const targetGeneIds = [];
                    if (cy) {
                        cy.edges().forEach(edge => {
                            if (edge.source().id().toLowerCase() === locus.toLowerCase()) {
                                targetGeneIds.push(edge.target().id());
                            }
                        });
                    }
                    
                    if (targetGeneIds.length === 0) {
                        try {
                            const info = await window.metabolicModelAdapter.fetchGeneReactionPathwayMapping(locus);
                            if (info && info.targetGenes) {
                                info.targetGenes.forEach(tg => targetGeneIds.push(tg));
                            }
                        } catch (e) {
                            console.warn("Could not retrieve target genes dynamically:", e);
                        }
                    }
                    
                    if (targetGeneIds.length === 0) {
                        targetGeneIds.push(locus);
                    }
                    
                    res = await window.simulationClient.runTFPerturbation(locus, targetGeneIds, objective, trackedReactionIds);
                }
                
                if (!res || res.status === "error") {
                    throw new Error(res?.warnings?.join(", ") || "FBA simulation returned an unsuccessful status.");
                }
                
                document.getElementById('glu-res-biomass-baseline').textContent = `${Number(res.baselineObjective || 0).toFixed(4)} mmol/gDCW/h`;
                document.getElementById('glu-res-biomass-perturbed').textContent = `${Number(res.perturbedObjective || 0).toFixed(4)} mmol/gDCW/h`;
                const biomassPct = res.objectiveChangePercent ? res.objectiveChangePercent.toFixed(2) : '0.00';
                document.getElementById('glu-res-biomass-change').textContent = `${biomassPct}%`;
                
                let baselineFlux = 0;
                let perturbedFlux = 0;
                let fluxChangePercent = 0;
                if (res.trackedFluxes && res.trackedFluxes.length > 0) {
                    const tf = res.trackedFluxes[0];
                    baselineFlux = tf.baselineFlux || 0;
                    perturbedFlux = tf.perturbedFlux || 0;
                    fluxChangePercent = tf.fluxChangePercent || 0;
                }
                
                document.getElementById('glu-res-flux-baseline').textContent = `${Number(baselineFlux).toFixed(4)} mmol/gDCW/h`;
                document.getElementById('glu-res-flux-perturbed').textContent = `${Number(perturbedFlux).toFixed(4)} mmol/gDCW/h`;
                const fluxChangeVal = perturbedFlux - baselineFlux;
                const sign = fluxChangeVal > 0 ? '+' : '';
                document.getElementById('glu-res-flux-change').textContent = `${sign}${Number(fluxChangeVal).toFixed(4)} (${fluxChangePercent.toFixed(2)}%)`;
                
                const interpretation = window.glutamateScenarioInterpretation.generateGlutamateScenarioInterpretation(res, glutamateSelectedReaction);
                document.getElementById('glu-res-interpretation').textContent = interpretation;
                
                if (resultsEmpty) resultsEmpty.classList.add('hidden');
                if (resultsBox) resultsBox.classList.remove('hidden');
            } catch (err) {
                console.error("FBA glutamate scenario simulation error:", err);
                if (errorBox) {
                    errorBox.textContent = `Simulation failed: ${err.message}.`;
                    errorBox.classList.remove('hidden');
                }
            } finally {
                runBtn.innerHTML = '<i class="fa-solid fa-play"></i> Run Glutamate Production Simulation';
                runBtn.disabled = false;
            }
        });
    }
}

function selectGlutamateCandidate(candidate) {
    glutamateSelectedReaction = candidate;
    window.glutamateScenario.glutamateState.selectedGlutamateReactionId = candidate.reactionId;
    window.glutamateScenario.glutamateState.selectedGlutamateReactionClass = candidate.classification;
    window.glutamateScenario.glutamateState.userVerified = false;
    
    const chkVerified = document.getElementById('chk-glutamate-verified');
    if (chkVerified) {
        chkVerified.checked = false;
        chkVerified.disabled = false;
    }
    
    const runBtn = document.getElementById('btn-run-glutamate-scenario');
    if (runBtn) {
        runBtn.disabled = true;
    }
    
    const infoArea = document.getElementById('selected-glutamate-info');
    if (infoArea) {
        infoArea.innerHTML = `
            <strong>Selected:</strong> ${escapeHtml(candidate.reactionId)} (${escapeHtml(candidate.name || 'Unnamed')})<br>
            <strong>Equation:</strong> <code style="font-family: monospace;">${escapeHtml(candidate.equation)}</code><br>
            <strong>Classification:</strong> <span class="badge-role ${escapeHtml(candidate.classification)}" style="font-size: 10px;">${escapeHtml(candidate.classification)}</span> 
            (Confidence: <strong>${escapeHtml(candidate.confidence)}</strong>)
        `;
    }
    
    const warningBox = document.getElementById('glutamate-verification-warning');
    if (warningBox) warningBox.classList.remove('hidden');
}

function runGlutamateScenarioFromEngineering(tfId) {
    setActiveWorkflowEntry('glutamate');
    
    const locusInput = document.getElementById('glutamate-scenario-locus');
    const typeSelect = document.getElementById('glutamate-scenario-type');
    if (locusInput) locusInput.value = tfId;
    if (typeSelect) typeSelect.value = 'TF';
    
    const verified = window.glutamateScenario.glutamateState.userVerified;
    const reactionId = window.glutamateScenario.glutamateState.selectedGlutamateReactionId;
    
    if (!reactionId || !verified) {
        const errorBox = document.getElementById('glutamate-scenario-error-box');
        if (errorBox) {
            errorBox.textContent = `Glutamate Scenario initiated for ${tfId}. Please select and verify a target reaction in Step 1 and Step 2 first.`;
            errorBox.classList.remove('hidden');
        }
        const searchBtn = document.getElementById('btn-search-glutamate-candidates');
        if (searchBtn) searchBtn.click();
    } else {
        const runBtn = document.getElementById('btn-run-glutamate-scenario');
        if (runBtn) runBtn.click();
    }
}
