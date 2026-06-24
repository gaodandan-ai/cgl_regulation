/**
 * C. glutamicum Regulatory Network Explorer - Client Side Logic
 * Uses Cytoscape.js and PapaParse
 */

// Application State
let regulations = [];
let rnaRegulations = [];
let geneMapping = [];
let cglToCg = {};
let cgToCgl = {};
let nameToCg = {};
let cgToProduct = {};
let geneIndex = {}; // lowercase -> { locusTag, name, type }
let geneToOperon = {}; // lower -> { operon, orientation, genes }
let searchSuggestions = [];
let currentQueryGene = null;
let cy = null;
let currentSimulationMode = null;
let currentSimulationRegulator = null;

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
const REGULATIONS_URL = 'data/regulations.csv';
const RNA_REGULATIONS_URL = 'data/rna_regulation.csv';
const MAPPING_URL = 'data/gene_mapping.csv';
const OPERONS_URL = 'data/operons.csv';

// ==========================================================================
// 1. Initialization & Data Loading
// ==========================================================================
document.addEventListener('DOMContentLoaded', () => {
    initEventListeners();
    initSidebarResizer();
    loadNetworkData();
});

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
        updateStatus('正在加载基因命名映射数据...', 'loading');
        const mapResponse = await fetch(MAPPING_URL);
        if (mapResponse.ok) {
            const mapText = await mapResponse.text();
            geneMapping = parseCSV(mapText);
            console.log(`Loaded ${geneMapping.length} gene mapping records.`);
        } else {
            console.warn('plate_gene_mapping.csv file not found. Skipping mapping.');
        }

        updateStatus('正在加载 TF-TG 调控数据...', 'loading');
        const tfResponse = await fetch(REGULATIONS_URL);
        if (!tfResponse.ok) throw new Error('无法读取 regulations.csv，请确认本地服务已启动。');
        const tfText = await tfResponse.text();
        
        regulations = parseCSV(tfText);
        console.log(`Loaded ${regulations.length} TF-TG regulations.`);

        updateStatus('正在加载 sRNA-mRNA 调控数据...', 'loading');
        const rnaResponse = await fetch(RNA_REGULATIONS_URL);
        if (rnaResponse.ok) {
            const rnaText = await rnaResponse.text();
            rnaRegulations = parseCSV(rnaText);
            console.log(`Loaded ${rnaRegulations.length} sRNA-mRNA regulations.`);
            srnaThresholdPanel.classList.remove('hidden');
        } else {
            console.warn('sRNA-mRNA regulations.csv file not found. Skipping sRNA data.');
        }

        updateStatus('正在加载操纵子结构数据...', 'loading');
        const operonResponse = await fetch(OPERONS_URL);
        if (operonResponse.ok) {
            const operonText = await operonResponse.text();
            parseOperons(operonText);
            console.log(`Loaded operons mapping.`);
        } else {
            console.warn('Operons file not found. Skipping operons data.');
        }

        buildGeneIndex();
        updateStatus('数据已就绪', 'success');
    } catch (err) {
        console.error(err);
        updateStatus('数据加载失败: ' + err.message, 'error');
        alert('错误：无法加载 CSV 文件。请确保运行了 python run_server.py 以便浏览器加载数据。');
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
        alert('请输入或粘贴至少一个基因或sRNA进行分析。');
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
        alert(`未在本地数据库中匹配到输入的基因/sRNA："${queries.join(', ')}"。`);
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
        alert("基于当前的过滤条件，该基因没有任何可见的调控关系。");
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
            // Edge styling
            {
                selector: 'edge',
                style: {
                    'width': 2,
                    'line-color': '#e65100', // Default dark orange
                    'target-arrow-color': '#e65100',
                    'target-arrow-shape': 'triangle',
                    'curve-style': 'bezier',
                    'arrow-scale': 1.1,
                    'transition-property': 'line-color, target-arrow-color, opacity, width',
                    'transition-duration': '0.2s'
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
                selector: 'edge[role="R"]', // Repression
                style: {
                    'line-color': '#d32f2f', // Academic Red
                    'target-arrow-color': '#d32f2f',
                    'target-arrow-shape': 'tee'
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
                selector: 'edge[role="sRNA"]', // sRNA-mRNA prediction
                style: {
                    'line-color': '#7b1fa2', // Academic Purple
                    'target-arrow-color': '#7b1fa2',
                    'line-style': 'dashed',
                    'target-arrow-shape': 'triangle-tee'
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

    // 6. Update Network Statistics
    updateNetworkStatistics();
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

    // Local filter state
    const showActivation = filterActivation.checked;
    const showRepression = filterRepression.checked;
    const showDual = filterDual.checked;
    const showSrna = filterSrna.checked;
    const rankLimit = parseInt(srnaRankThreshold.value, 10);
    const showOnlyTfTargets = filterOnlyTfTargets ? filterOnlyTfTargets.checked : false;

    // Add all query nodes first
    queryList.forEach(locus => {
        const lower = locus.toLowerCase();
        let meta = { locusTag: locus, name: locus, type: 'Target' };
        for (let key in geneIndex) {
            if (geneIndex[key].locusTag.toLowerCase() === lower) {
                meta = geneIndex[key];
                break;
            }
        }
        nodesMap[locus] = {
            data: {
                id: locus,
                name: getPrioritizedLabel(locus, meta.name),
                type: 'query'
            }
        };
    });

    // 1. Process TF regulations
    regulations.forEach(row => {
        const tfTag = cleanStr(row.TF_locusTag);
        const tfName = cleanStr(row.TF_name);
        const tgTag = cleanStr(row.TG_locusTag);
        const tgName = cleanStr(row.TG_name);
        const role = cleanStr(row.Role); // A, R, Dual, etc.

        // Edge matching role filter
        if (role === 'A' && !showActivation) return;
        if (role === 'R' && !showRepression) return;
        if ((role === 'Dual' || role === 'Sigma' || role === '') && !showDual) return;

        // Is ANY query involved?
        const isTfQuery = querySet.has(tfTag.toLowerCase());
        const isTgQuery = querySet.has(tgTag.toLowerCase());

        if (isTfQuery || isTgQuery) {
            // TF targets filter
            if (showOnlyTfTargets && isTfQuery && !isTgQuery) {
                const targetMeta = geneIndex[tgTag.toLowerCase()];
                const isTargetTf = targetMeta && targetMeta.type === 'TF';
                if (!isTargetTf) return;
            }
            
            // Add nodes
            if (!nodesMap[tfTag]) {
                nodesMap[tfTag] = {
                    data: {
                        id: tfTag,
                        name: getPrioritizedLabel(tfTag, tfName),
                        type: querySet.has(tfTag.toLowerCase()) ? 'query' : 'TF'
                    }
                };
            }
            if (!nodesMap[tgTag]) {
                nodesMap[tgTag] = {
                    data: {
                        id: tgTag,
                        name: getPrioritizedLabel(tgTag, tgName),
                        type: querySet.has(tgTag.toLowerCase()) ? 'query' : 'Target'
                    }
                };
            }
            
            // Add edge
            edges.push({
                data: {
                    id: `edge_${tfTag}_${tgTag}`,
                    source: tfTag,
                    target: tgTag,
                    role: role,
                    type: 'TF-TG'
                }
            });
        }
    });

    // 2. Process sRNA regulations
    if (showSrna && rnaRegulations.length > 0) {
        rnaRegulations.forEach(row => {
            const srna = cleanStr(row.srna);
            const mrna = cleanStr(row.mrna);
            const rank = parseInt(row.rank, 10);
            
            if (rank > rankLimit) return; // Slider filter

            const isSrnaQuery = querySet.has(srna.toLowerCase());
            const isMrnaQuery = querySet.has(mrna.toLowerCase());

            if (isSrnaQuery || isMrnaQuery) {
                // TF targets filter for sRNA targets
                if (showOnlyTfTargets && isSrnaQuery && !isMrnaQuery) {
                    const targetMeta = geneIndex[mrna.toLowerCase()];
                    const isTargetTf = targetMeta && targetMeta.type === 'TF';
                    if (!isTargetTf) return;
                }
                
                if (!nodesMap[srna]) {
                    nodesMap[srna] = {
                        data: {
                            id: srna,
                            name: getPrioritizedLabel(srna, srna),
                            type: querySet.has(srna.toLowerCase()) ? 'query' : 'sRNA'
                        }
                    };
                }
                if (!nodesMap[mrna]) {
                    nodesMap[mrna] = {
                        data: {
                            id: mrna,
                            name: getPrioritizedLabel(mrna, mrna),
                            type: querySet.has(mrna.toLowerCase()) ? 'query' : 'Target'
                        }
                    };
                }

                edges.push({
                    data: {
                        id: `edge_srna_${srna}_${mrna}`,
                        source: srna,
                        target: mrna,
                        role: 'sRNA',
                        type: 'sRNA-mRNA',
                        rank: rank,
                        energy: row.energy,
                        pvalue: row.copra_pvalue
                    }
                });
            }
        });
    }

    // 3. Filter for co-regulated target genes if checked
    const showOnlyCoRegulated = filterCoregulated.checked;
    if (showOnlyCoRegulated) {
        // Count incoming edges for Target nodes in the current edges set
        const inDegreeMap = {};
        edges.forEach(e => {
            const target = e.data.target;
            inDegreeMap[target] = (inDegreeMap[target] || 0) + 1;
        });

        // Determine which Target nodes have inDegree >= 2
        const coRegulatedTargets = new Set();
        Object.keys(inDegreeMap).forEach(nodeId => {
            const nodeObj = nodesMap[nodeId];
            if (nodeObj && (nodeObj.data.type === 'Target') && inDegreeMap[nodeId] >= 2) {
                coRegulatedTargets.add(nodeId);
            }
        });

        // Filter edges: only keep edges pointing to co-regulated target genes
        const keptEdges = edges.filter(e => {
            const targetNode = nodesMap[e.data.target];
            const targetType = targetNode ? targetNode.data.type : '';
            if (targetType === 'Target') {
                return coRegulatedTargets.has(e.data.target);
            }
            return true; // Keep edges pointing to queries/TFs/sRNAs
        });

        // Nodes must be kept if they are in the query list or involved in kept edges
        const keptNodeIds = new Set(queryList);
        keptEdges.forEach(e => {
            keptNodeIds.add(e.data.source);
            keptNodeIds.add(e.data.target);
        });

        const keptNodes = Object.values(nodesMap).filter(n => keptNodeIds.has(n.data.id));

        return {
            nodes: keptNodes,
            edges: keptEdges
        };
    }

    return {
        nodes: Object.values(nodesMap),
        edges: edges
    };
}

// Highlight subnetwork (1st degree neighbors)
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

    // Clear previous AI summary
    const summaryCard = document.getElementById('ai-summary-result');
    if (summaryCard) {
        summaryCard.classList.add('hidden');
        summaryCard.innerHTML = '';
    }

    const lower = locusTag.toLowerCase();
    
    // Resolve display meta
    let meta = { locusTag: locusTag, name: locusTag, type: 'Target' };
    for (let key in geneIndex) {
        if (geneIndex[key].locusTag.toLowerCase() === lower) {
            meta = geneIndex[key];
            break;
        }
    }

    // Set badge style
    detailTypeBadge.className = `gene-badge ${meta.type.toLowerCase()}`;
    detailTypeBadge.textContent = meta.type === 'TF' ? '转录因子 (TF)' : meta.type === 'sRNA' ? '小RNA (sRNA)' : '靶基因 (Target)';
    
    const cgl = cgToCgl[lower];
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
        pathwayContainer.innerHTML = '<span style="font-size: 11px; color: var(--text-muted);"><i class="fa-solid fa-spinner fa-spin"></i> 正在加载通路数据...</span>';
        
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
                    badge.title = `KEGG通路: ${p.id} (点击在新标签页中查看地图并高亮基因)`;
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
            <div style="font-weight: 600; color: var(--text-primary);">${operonMeta.operon} (${operonMeta.orientation} 向)</div>
            <div style="font-size: 11px; margin-top: 4px; color: var(--text-secondary);">包含基因: ${geneLinks}</div>
            <button id="btn-draw-operon-network" class="secondary-btn" style="margin-top: 8px; font-size: 11px; padding: 6px 10px; height: auto; width: 100%; border: 1px solid rgba(30, 58, 138, 0.15); color: var(--color-primary-accent); background-color: rgba(30, 58, 138, 0.03);">
                <i class="fa-solid fa-network-wired"></i> 联合分析操纵子成员 (${operonMeta.genes.length} 基因)
            </button>
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
        dbLinks.push(`<a href="https://www.kegg.jp/entry/cgl:${cglLocusForKegg}" target="_blank" class="ext-link" title="在 KEGG 中查看代谢通路"><i class="fa-solid fa-diagram-project"></i> KEGG</a>`);
    } else if (standardCgForLinks.toLowerCase().startsWith('cg')) {
        // Fallback guess if no direct mapping exists but is a coding gene
        const predictedCgl = standardCgForLinks.replace('cg', 'Cgl');
        dbLinks.push(`<a href="https://www.kegg.jp/entry/cgl:${predictedCgl}" target="_blank" class="ext-link" title="在 KEGG 中查看代谢通路"><i class="fa-solid fa-diagram-project"></i> KEGG</a>`);
    }
    
    if (standardCgForLinks.toLowerCase().startsWith('cg')) {
        dbLinks.push(`<a href="https://www.ncbi.nlm.nih.gov/gene/?term=${standardCgForLinks}" target="_blank" class="ext-link" title="在 NCBI Gene 查看官方注释"><i class="fa-solid fa-dna"></i> NCBI</a>`);
        dbLinks.push(`<a href="https://biocyc.org/getid?id=CORYNE:${standardCgForLinks}" target="_blank" class="ext-link" title="在 BioCyc / CoryneCyc 谷棒专属数据库中查看详细通路"><i class="fa-solid fa-database"></i> BioCyc</a>`);
    } else {
        dbLinks.push(`<a href="https://www.ncbi.nlm.nih.gov/search/all/?term=${standardCgForLinks}" target="_blank" class="ext-link" title="在 NCBI 中检索"><i class="fa-solid fa-magnifying-glass"></i> NCBI</a>`);
    }
    
    dbLinks.push(`<a href="https://cosy.bio/coryneregnet" target="_blank" class="ext-link" title="在 CoryneRegNet 谷棒转录调控网络数据库中检索"><i class="fa-solid fa-network-wired"></i> CoryneRegNet</a>`);
    dbLinks.push(`<a href="https://www.uniprot.org/uniprotkb?query=gene:${standardCgForLinks}" target="_blank" class="ext-link" title="在 UniProt 中查看蛋白功能"><i class="fa-solid fa-graduation-cap"></i> UniProt</a>`);
    
    // 文献追踪链接
    const pubmedQuery = encodeURIComponent(`"Corynebacterium glutamicum" AND (${standardCgForLinks}${meta.name && meta.name !== '--' && meta.name !== standardCgForLinks ? ' OR ' + meta.name : ''})`);
    dbLinks.push(`<a href="https://pubmed.ncbi.nlm.nih.gov/?term=${pubmedQuery}" target="_blank" class="ext-link" title="在 PubMed 检索该基因相关的科研文献"><i class="fa-solid fa-book-open"></i> PubMed 文献</a>`);
    
    const scholarQuery = encodeURIComponent(`"Corynebacterium glutamicum" "${standardCgForLinks}"${meta.name && meta.name !== '--' && meta.name !== standardCgForLinks ? ' OR "' + meta.name + '"' : ''}`);
    dbLinks.push(`<a href="https://scholar.google.com/scholar?q=${scholarQuery}" target="_blank" class="ext-link" title="在 Google 学术中检索该基因文献"><i class="fa-solid fa-graduation-cap"></i> 谷歌学术</a>`);
    
    linksCell.innerHTML = dbLinks.join('');

    // Load relationships from global collections
    const relations = [];
    let regsCount = 0;
    let targsCount = 0;

    // TF-TG details
    regulations.forEach(row => {
        const tfTag = cleanStr(row.TF_locusTag);
        const tfName = cleanStr(row.TF_name);
        const tgTag = cleanStr(row.TG_locusTag);
        const tgName = cleanStr(row.TG_name);
        const role = cleanStr(row.Role);

        if (tfTag.toLowerCase() === lower) {
            targsCount++;
            relations.push({
                gene: getPrioritizedLabel(tgTag, tgName),
                locusTag: tgTag,
                dir: 'outgoing',
                role: role,
                source: 'CorynebNet'
            });
        }
        if (tgTag.toLowerCase() === lower) {
            regsCount++;
            relations.push({
                gene: getPrioritizedLabel(tfTag, tfName),
                locusTag: tfTag,
                dir: 'incoming',
                role: role,
                source: 'CorynebNet'
            });
        }
    });

    // sRNA details
    rnaRegulations.forEach(row => {
        const srna = cleanStr(row.srna);
        const mrna = cleanStr(row.mrna);
        const rank = parseInt(row.rank, 10);
        const energy = row.energy;

        if (srna.toLowerCase() === lower) {
            targsCount++;
            relations.push({
                gene: getPrioritizedLabel(mrna, mrna),
                locusTag: mrna,
                dir: 'outgoing',
                role: 'sRNA',
                source: `Rank: ${rank} (E: ${energy})`
            });
        }
        if (mrna.toLowerCase() === lower) {
            regsCount++;
            relations.push({
                gene: getPrioritizedLabel(srna, srna),
                locusTag: srna,
                dir: 'incoming',
                role: 'sRNA',
                source: `Rank: ${rank} (E: ${energy})`
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
            alert('该基因暂无上游调控因子。');
        }
    };

    targetCard.onclick = () => {
        if (outgoingLoci.length > 0) {
            queryMultipleGenes(outgoingLoci);
        } else {
            alert('该基因暂无下游靶标。');
        }
    };

    // Render Table
    relationsTableBody.innerHTML = '';
    
    if (relations.length === 0) {
        relationsTableBody.innerHTML = `<tr><td colspan="4" class="text-muted" style="text-align:center;">暂无调控明细数据</td></tr>`;
    } else {
        // Sort: Incoming first, then outgoing
        relations.sort((a, b) => a.dir.localeCompare(b.dir));
        
        relations.forEach(rel => {
            const tr = document.createElement('tr');
            
            const roleClass = rel.role === 'A' ? 'activation' : rel.role === 'R' ? 'repression' : rel.role === 'sRNA' ? 'srna' : 'dual';
            const roleText = rel.role === 'A' ? '激活 (+)' : rel.role === 'R' ? '抑制 (-)' : rel.role === 'sRNA' ? 'sRNA预测' : '双重/Sigma';
            
            tr.innerHTML = `
                <td><a href="#" class="gene-link" data-locus="${rel.locusTag}">${rel.gene}</a></td>
                <td><span class="badge-dir ${rel.dir}">${rel.dir === 'incoming' ? '← 上游' : '下游 →'}</span></td>
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

    // Slide open sidebar
    toggleRightSidebar(true);
}

// ==========================================================================
// 5. DOM Event Listeners & Interactive Controls
// ==========================================================================
function initEventListeners() {
    // Initialize first empty input row
    clearAllInputs();

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
    filterSrna.addEventListener('change', reRender);
    
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
    input.placeholder = '输入基因/sRNA名称';
    input.autocomplete = 'off';
    
    wrapper.appendChild(input);
    row.appendChild(wrapper);
    
    // Add delete or add button based on current rows count
    const existingRows = geneInputsContainer.querySelectorAll('.gene-input-row');
    if (existingRows.length > 0) {
        const removeBtn = document.createElement('button');
        removeBtn.className = 'remove-row-btn';
        removeBtn.title = '删除基因栏';
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
        addBtn.title = '添加基因栏';
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
    showNodeDetails(locus);
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
        if (idx === 0) {
            activeInput = input;
        }
    });
    
    triggerSearchFromInputs();
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

// ==========================================================================
// 7. AI Literature Summary & Function Assistant
// ==========================================================================
function initAiSummaryFeature() {
    const btnSaveKey = document.getElementById('btn-save-key');
    const btnClearKey = document.getElementById('btn-clear-key');
    const btnTriggerAi = document.getElementById('btn-trigger-ai');
    const apiKeyInput = document.getElementById('gemini-api-key-input');
    const keyConfigPanel = document.getElementById('ai-key-config-panel');
    const keyActivePanel = document.getElementById('ai-key-active-panel');

    // 1. Load key from localStorage on initialize
    const savedKey = localStorage.getItem('gemini_api_key');
    if (savedKey) {
        keyConfigPanel.classList.add('hidden');
        keyActivePanel.classList.remove('hidden');
        btnTriggerAi.disabled = false;
    } else {
        keyConfigPanel.classList.remove('hidden');
        keyActivePanel.classList.add('hidden');
        btnTriggerAi.disabled = true;
    }

    // 2. Save Key listener
    btnSaveKey.addEventListener('click', () => {
        const key = apiKeyInput.value.trim();
        if (!key) {
            alert('请输入有效的 Gemini API Key。');
            return;
        }
        localStorage.setItem('gemini_api_key', key);
        apiKeyInput.value = '';
        
        keyConfigPanel.classList.add('hidden');
        keyActivePanel.classList.remove('hidden');
        btnTriggerAi.disabled = false;
    });

    // 3. Clear Key listener
    btnClearKey.addEventListener('click', () => {
        localStorage.removeItem('gemini_api_key');
        keyConfigPanel.classList.remove('hidden');
        keyActivePanel.classList.add('hidden');
        btnTriggerAi.disabled = true;
        
        const summaryCard = document.getElementById('ai-summary-result');
        summaryCard.classList.add('hidden');
        summaryCard.innerHTML = '';
    });

    // 4. AI Trigger listener
    btnTriggerAi.addEventListener('click', () => {
        triggerAiSummary();
    });
}

async function triggerAiSummary() {
    const btnTriggerAi = document.getElementById('btn-trigger-ai');
    const summaryCard = document.getElementById('ai-summary-result');
    
    const locus = document.getElementById('info-locus').textContent.trim();
    const name = document.getElementById('info-name').textContent.trim();
    const apiKey = localStorage.getItem('gemini_api_key');
    
    if (!locus || locus === '-') {
        alert('请先选择一个基因进行分析。');
        return;
    }
    if (!apiKey) {
        alert('请先在上方配置您的 Gemini API Key。');
        return;
    }
    
    // Set loading state
    btnTriggerAi.disabled = true;
    summaryCard.classList.remove('hidden');
    summaryCard.classList.add('loading');
    summaryCard.innerHTML = `
        <div class="ai-spinner"></div>
        <span style="font-weight: 500;">正在检索 PubMed 文献并请求 AI 总结中...</span>
    `;
    
    try {
        const response = await fetch(`/api/summarize?gene=${locus}&name=${name}`, {
            headers: {
                'X-Gemini-API-Key': apiKey
            }
        });
        
        if (!response.ok) {
            throw new Error(`HTTP 错误: ${response.status}`);
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
                    <div class="ai-sources-title"><i class="fa-solid fa-book"></i> 参考 PubMed 文献 (${result.papers.length} 篇)</div>
            `;
            
            result.papers.forEach(p => {
                htmlContent += `
                    <div class="ai-source-item">
                        <i class="fa-solid fa-file-lines"></i>
                        <a href="https://pubmed.ncbi.nlm.nih.gov/${p.pmid}" target="_blank" class="ai-source-link" title="点击在 PubMed 查看原始文献">
                            ${p.title} (PMID: ${p.pmid}) <i class="fa-solid fa-arrow-up-right-from-square" style="font-size: 8px;"></i>
                        </a>
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
                <i class="fa-solid fa-circle-exclamation"></i> 总结生成失败
            </div>
            <p style="font-size: 11px; color: var(--text-secondary); line-height: 1.4;">
                ${err.message || '未知网络错误，请检查您的 API Key 是否正确或网络连接状态。'}
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
function initAiPathwayFeature() {
    const btnAnalyze = document.getElementById('btn-analyze-pathway');
    const inputPathway = document.getElementById('ai-pathway-input');
    const resultCard = document.getElementById('ai-pathway-result');

    if (!btnAnalyze || !inputPathway || !resultCard) return;

    btnAnalyze.addEventListener('click', async () => {
        const query = inputPathway.value.trim();
        if (!query) {
            alert('请输入要分析的代谢通路或生理功能名称。');
            return;
        }

        const apiKey = localStorage.getItem('gemini_api_key');
        if (!apiKey) {
            alert('要使用 AI 通路分析，请先在右侧详情面板配置您的 Gemini API Key！');
            // Open the detail sidebar if it is collapsed, and highlight the key input
            toggleRightSidebar(true);
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
            <span style="font-weight: 500;">AI 正在分析通路基因...</span>
        `;

        try {
            const response = await fetch(`/api/pathway?pathway=${encodeURIComponent(query)}`, {
                headers: {
                    'X-Gemini-API-Key': apiKey
                }
            });

            if (!response.ok) {
                throw new Error(`HTTP 错误: ${response.status}`);
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
                genesBadgesHtml = '<span style="color: var(--text-secondary); font-size: 11px;">未识别到关联 locus tags</span>';
            }

            resultCard.innerHTML = `
                <div class="ai-pathway-summary">${result.summary || '无总结信息'}</div>
                <div class="ai-pathway-genes-title"><i class="fa-solid fa-dna"></i> 关联基因 (${genes.length})</div>
                <div class="ai-pathway-genes-list">${genesBadgesHtml}</div>
                ${genes.length > 0 ? `
                    <button class="ai-pathway-draw-btn" id="btn-draw-pathway-network">
                        <i class="fa-solid fa-network-wired"></i> 一键绘制该通路调控网络
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
                    <i class="fa-solid fa-circle-exclamation"></i> 分析失败
                </div>
                <p style="font-size: 11px; color: var(--text-secondary); line-height: 1.4;">
                    ${err.message || '未知网络错误，请检查您的 API Key 是否正确或网络连接状态。'}
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
            alert('请输入要分析的基因功能描述、转录因子或特征。');
            return;
        }

        const apiKey = localStorage.getItem('gemini_api_key');
        if (!apiKey) {
            alert('要使用 AI 基因分析，请先在右侧详情面板配置您的 Gemini API Key！');
            // Open the detail sidebar if it is collapsed, and highlight the key input
            toggleRightSidebar(true);
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
            <span style="font-weight: 500;">AI 正在分析基因特征...</span>
        `;

        try {
            const response = await fetch(`/api/gene_assistant?query=${encodeURIComponent(query)}`, {
                headers: {
                    'X-Gemini-API-Key': apiKey
                }
            });

            if (!response.ok) {
                throw new Error(`HTTP 错误: ${response.status}`);
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
                genesBadgesHtml = '<span style="color: var(--text-secondary); font-size: 11px;">未识别到关联 locus tags</span>';
            }

            resultCard.innerHTML = `
                <div class="ai-pathway-summary">${result.summary || '无总结信息'}</div>
                <div class="ai-pathway-genes-title"><i class="fa-solid fa-dna"></i> 关联基因 (${genes.length})</div>
                <div class="ai-pathway-genes-list">${genesBadgesHtml}</div>
                ${genes.length > 0 ? `
                    <button class="ai-pathway-draw-btn" id="btn-draw-gene-network">
                        <i class="fa-solid fa-network-wired"></i> 一键绘制该基因调控网络
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
                    <i class="fa-solid fa-circle-exclamation"></i> 分析失败
                </div>
                <p style="font-size: 11px; color: var(--text-secondary); line-height: 1.4;">
                    ${err.message || '未知网络错误，请检查您的 API Key 是否正确或网络连接状态。'}
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

// ==========================================================================
// 10. Floating UI & Helper Actions (CSV Export, Search Focus, Stats)
// ==========================================================================
function initSidebarResizer() {
    const resizer = document.getElementById('sidebar-resizer');
    const sidebar = document.getElementById('right-sidebar');
    if (!resizer || !sidebar) return;

    // Load saved width from localStorage if exists
    const savedWidth = localStorage.getItem('right-sidebar-width');
    if (savedWidth) {
        document.documentElement.style.setProperty('--right-sidebar-width', savedWidth);
    }

    let startX = 0;
    let startWidth = 0;

    function onMouseMove(e) {
        const deltaX = e.clientX - startX;
        let newWidth = startWidth - deltaX; // Drag left (negative deltaX) makes it wider
        
        // Enforce limits: min 280px, max 80% of window width
        const minWidth = 280;
        const maxWidth = window.innerWidth * 0.8;
        if (newWidth < minWidth) newWidth = minWidth;
        if (newWidth > maxWidth) newWidth = maxWidth;

        document.documentElement.style.setProperty('--right-sidebar-width', newWidth + 'px');
        
        // Notify Cytoscape of layout resize
        if (cy) {
            cy.resize();
        }
    }

    function onMouseUp() {
        document.removeEventListener('mousemove', onMouseMove);
        document.removeEventListener('mouseup', onMouseUp);
        sidebar.classList.remove('sidebar-no-transition');
        resizer.classList.remove('resizing');
        
        // Save current width to localStorage
        const currentWidth = getComputedStyle(sidebar).width;
        localStorage.setItem('right-sidebar-width', currentWidth);
        
        if (cy) {
            cy.resize();
        }
    }

    resizer.addEventListener('mousedown', (e) => {
        e.preventDefault(); // Prevent text selection
        startX = e.clientX;
        startWidth = parseInt(getComputedStyle(sidebar).width, 10);
        
        sidebar.classList.add('sidebar-no-transition');
        resizer.classList.add('resizing');

        document.addEventListener('mousemove', onMouseMove);
        document.addEventListener('mouseup', onMouseUp);
    });
}

function toggleRightSidebar(open) {
    const rightSidebar = document.getElementById('right-sidebar');
    const searchContainer = document.getElementById('canvas-search-container');
    const statsContainer = document.getElementById('canvas-stats-container');
    
    if (open) {
        rightSidebar.classList.remove('collapsed');
        searchContainer?.classList.add('sidebar-open');
        statsContainer?.classList.add('sidebar-open');
    } else {
        rightSidebar.classList.add('collapsed');
        searchContainer?.classList.remove('sidebar-open');
        statsContainer?.classList.remove('sidebar-open');
        resetHighlight();
    }
}

function exportNetworkToCsv() {
    if (!cy) {
        alert('当前没有可导出的网络。');
        return;
    }

    const edges = cy.edges();
    if (edges.length === 0) {
        alert('当前网络中无任何调控边关系。');
        return;
    }

    // CSV headers (with UTF-8 BOM)
    let csvContent = '\uFEFF';
    csvContent += '源节点Locus Tag(Source Locus),源节点名称(Source Name),源节点功能(Source Function),目标节点Locus Tag(Target Locus),目标节点名称(Target Name),目标节点功能(Target Function),调控类型(Interaction),调控作用(Role),数据来源/分数(Source)';
    
    if (currentSimulationMode) {
        csvContent += `,预测转录效应(Predicted Effect under ${currentSimulationMode === 'OE' ? 'OE' : 'KO'})`;
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
        const sourceFunc = cgToProduct[sourceLower] || '暂无详细功能描述';
        const targetFunc = cgToProduct[targetLower] || '暂无详细功能描述';

        const type = edge.data('type') || '';
        const role = edge.data('role') || '';
        let roleText = '未知';
        if (role === 'A') roleText = '激活 (+)';
        else if (role === 'R') roleText = '抑制 (-)';
        else if (role === 'sRNA') roleText = 'sRNA预测';
        else roleText = '双重/Sigma';
        
        let sourceVal = '';
        if (type === 'TF-TG') {
            sourceVal = 'CorynebNet';
        } else {
            const rank = edge.data('rank') || '';
            const energy = edge.data('energy') || '';
            sourceVal = `sRNA预测 (Rank: ${rank}, Energy: ${energy})`;
        }

        let line = `${cleanVal(sourceId)},${cleanVal(sourceLabel)},${cleanVal(sourceFunc)},${cleanVal(targetId)},${cleanVal(targetLabel)},${cleanVal(targetFunc)},${cleanVal(type)},${cleanVal(roleText)},${cleanVal(sourceVal)}`;

        if (currentSimulationMode) {
            let effectText = '无明显效应';
            if (currentSimulationRegulator && sourceId.toLowerCase() === currentSimulationRegulator.toLowerCase()) {
                if (currentSimulationMode === 'OE') {
                    if (role === 'A') effectText = '表达增强 ⬆';
                    else if (role === 'R' || role === 'sRNA') effectText = '表达减弱 ⬇';
                    else effectText = '复杂/双重 ↕';
                } else if (currentSimulationMode === 'KO') {
                    if (role === 'A') effectText = '表达减弱 ⬇';
                    else if (role === 'R' || role === 'sRNA') effectText = '表达增强 ⬆';
                    else effectText = '复杂/双重 ↕';
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
            resultsBox.innerHTML = `<div style="padding: 10px; text-align: center; color: var(--text-muted); font-size: 11px;">画布中未找到该基因</div>`;
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

    const nodes = cy.nodes();
    const edges = cy.edges();

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
        hubsSpan.innerHTML = topHubs.map(h => `<strong style="font-family: monospace; color: var(--color-primary-accent);">${h.label}</strong> (${h.degree}条)`).join(', ');
    } else {
        hubsSpan.textContent = '无';
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
    
    if (currentQueryGene) {
        const currStr = JSON.stringify(currentQueryGene.map(l => l.toLowerCase()).sort());
        const nextStr = JSON.stringify((Array.isArray(locusTags) ? locusTags : [locusTags]).map(l => l.toLowerCase()).sort());
        
        if (currStr !== nextStr) {
            queryHistory.push(currentQueryGene);
            queryForwardHistory = []; // clear forward
        }
    }
    updateHistoryButtons();
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
        if (currentQueryGene) {
            queryForwardHistory.push(currentQueryGene);
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
        if (currentQueryGene) {
            queryHistory.push(currentQueryGene);
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

    const regulatorNode = cy.getElementById(regLocus);
    if (!regulatorNode || regulatorNode.length === 0) return;

    const outgoingEdges = regulatorNode.outgoers('edge');
    
    outgoingEdges.forEach(edge => {
        const targetNode = edge.target();
        const role = edge.data('role');

        let effect = 'none'; // 'up', 'down', 'dual'
        if (mode === 'OE') {
            if (role === 'A') effect = 'up';
            else if (role === 'R' || role === 'sRNA') effect = 'down';
            else effect = 'dual';
        } else if (mode === 'KO') {
            if (role === 'A') effect = 'down';
            else if (role === 'R' || role === 'sRNA') effect = 'up';
            else effect = 'dual';
        }

        const origName = targetNode.data('name') || targetNode.id();
        
        if (effect === 'up') {
            targetNode.addClass('sim-up');
            targetNode.data('name', `${origName} (⬆ 增强)`);
        } else if (effect === 'down') {
            targetNode.addClass('sim-down');
            targetNode.data('name', `${origName} (⬇ 减弱)`);
        } else if (effect === 'dual') {
            targetNode.addClass('sim-dual');
            targetNode.data('name', `${origName} (↕ 双重)`);
        }
    });

    const rows = document.querySelectorAll('#detail-relations-table tbody tr');
    rows.forEach(tr => {
        const dirSpan = tr.querySelector('.badge-dir');
        const roleSpan = tr.querySelector('.badge-role');
        if (dirSpan && dirSpan.classList.contains('outgoing') && roleSpan) {
            const roleClass = roleSpan.className;
            let role = 'Dual';
            if (roleClass.includes('activation')) role = 'A';
            else if (roleClass.includes('repression')) role = 'R';
            else if (roleClass.includes('srna')) role = 'sRNA';

            let effectText = '';
            let effectStyle = '';
            if (mode === 'OE') {
                if (role === 'A') { effectText = '表达增强 ⬆'; effectStyle = 'color: #2e7d32; font-weight: 600;'; }
                else if (role === 'R' || role === 'sRNA') { effectText = '表达减弱 ⬇'; effectStyle = 'color: #d32f2f; font-weight: 600;'; }
                else { effectText = '复杂/双重 ↕'; effectStyle = 'color: #e65100; font-weight: 600;'; }
            } else if (mode === 'KO') {
                if (role === 'A') { effectText = '表达减弱 ⬇'; effectStyle = 'color: #d32f2f; font-weight: 600;'; }
                else if (role === 'R' || role === 'sRNA') { effectText = '表达增强 ⬆'; effectStyle = 'color: #2e7d32; font-weight: 600;'; }
                else { effectText = '复杂/双重 ↕'; effectStyle = 'color: #e65100; font-weight: 600;'; }
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
        th.textContent = '预测效应';
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
            exportText.textContent = `导出预测效应表格 (${mode === 'OE' ? '过表达' : '敲除'})`;
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
        if (currentName.includes(' (⬆ 增强)') || currentName.includes(' (⬇ 减弱)') || currentName.includes(' (↕ 双重)')) {
            const clean = currentName.replace(' (⬆ 增强)', '').replace(' (⬇ 减弱)', '').replace(' (↕ 双重)', '');
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
        alert('当前无活跃的扰动模拟结果可供导出。');
        return;
    }

    const regLocus = currentSimulationRegulator;
    const mode = currentSimulationMode;

    const regulatorNode = cy.getElementById(regLocus);
    if (!regulatorNode || regulatorNode.length === 0) {
        alert('无法在画布中找到当前的调控因子节点。');
        return;
    }

    const outgoingEdges = regulatorNode.outgoers('edge');
    if (outgoingEdges.length === 0) {
        alert('当前调控因子没有下游靶基因关系。');
        return;
    }

    // Resolve regulator metadata
    const regLower = regLocus.toLowerCase();
    const regCgl = cgToCgl[regLower] || '';
    const regMeta = geneIndex[regLower] || { name: regLocus };
    const regName = regCgl ? regCgl : (regMeta.name && regMeta.name !== '--' ? regMeta.name : regLocus);

    // CSV headers (with UTF-8 BOM)
    let csvContent = '\uFEFF';
    csvContent += '调控基因Locus Tag(Regulator Locus),调控基因名称(Regulator Name),靶基因Locus Tag(Target Locus),靶基因名称(Target Name),调控关系(Interaction Role),扰动模式(Perturbation Mode),预测表达效应(Predicted Effect),靶基因功能描述(Target Function)\n';

    const cleanVal = (val) => {
        if (!val) return '';
        let s = String(val).replace(/"/g, '""');
        if (s.includes(',') || s.includes('\n') || s.includes('"')) {
            s = `"${s}"`;
        }
        return s;
    };

    outgoingEdges.forEach(edge => {
        const targetNode = edge.target();
        const targetId = targetNode.id();
        const targetLower = targetId.toLowerCase();
        
        // Resolve target name
        const targetCgl = cgToCgl[targetLower] || '';
        const targetMeta = geneIndex[targetLower] || { name: targetId };
        const targetName = targetCgl ? targetCgl : (targetMeta.name && targetMeta.name !== '--' ? targetMeta.name : targetId);

        const role = edge.data('role') || '';
        let roleText = '未知';
        if (role === 'A') roleText = '激活 (+)';
        else if (role === 'R') roleText = '抑制 (-)';
        else if (role === 'sRNA') roleText = 'sRNA预测';
        else roleText = '双重/Sigma';

        let effectText = '无明显效应';
        if (mode === 'OE') {
            if (role === 'A') effectText = '表达增强 ⬆';
            else if (role === 'R' || role === 'sRNA') effectText = '表达减弱 ⬇';
            else effectText = '复杂/双重 ↕';
        } else if (mode === 'KO') {
            if (role === 'A') effectText = '表达减弱 ⬇';
            else if (role === 'R' || role === 'sRNA') effectText = '表达增强 ⬆';
            else effectText = '复杂/双重 ↕';
        }

        // Get target function
        const targetFunc = cgToProduct[targetLower] || '暂无详细功能描述';

        csvContent += `${cleanVal(regLocus)},${cleanVal(regName)},${cleanVal(targetId)},${cleanVal(targetName)},${cleanVal(roleText)},${cleanVal(mode === 'OE' ? '过表达 (OE)' : '敲除 (KO)')},${cleanVal(effectText)},${cleanVal(targetFunc)}\n`;
    });

    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement('a');
    const url = URL.createObjectURL(blob);
    
    link.href = url;
    link.setAttribute('download', `${regLocus}_${mode}_predicted_effects.csv`);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}


