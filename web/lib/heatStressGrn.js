/**
 * C. glutamicum Regulatory Network Explorer
 * Heat Stress RNA-Seq & GRN Analysis Controller
 */

(function () {
    const state = {
        data: null,
        activeSubtab: 'grn',
        cyInstance: null,
        corrThreshold: 0.60,
        rfThreshold: 0.020,
        currentTimePoint: '1h',
        loading: false,
        dynamicChartInstance: null,
        metabChartInstance: null
    };

    function init() {
        console.log("Initializing Heat Stress GRN Analysis panel...");

        // Subtabs event bindings
        const subtabs = {
            'btn-subtab-grn': 'grn',
            'btn-subtab-dynamic': 'dynamic',
            'btn-subtab-causal': 'causal',
            'btn-subtab-enrich': 'enrich',
            'btn-subtab-time': 'time',
            'btn-subtab-rewire': 'rewire',
            'btn-subtab-metab': 'metab',
            'btn-subtab-ecfba': 'ecfba',
            'btn-subtab-mfa': 'mfa',
            'btn-subtab-figures': 'figures'
        };

        // Select selectors bindings
        const dynSelect = document.getElementById('rnaseq-dynamic-gene-select');
        if (dynSelect) {
            dynSelect.addEventListener('change', () => {
                renderDynamicGRN();
            });
        }

        const enrichSelect = document.getElementById('rnaseq-enrich-time-select');
        if (enrichSelect) {
            enrichSelect.addEventListener('change', () => {
                renderMotifEnrichment();
            });
        }

        Object.entries(subtabs).forEach(([btnId, tabName]) => {
            const btn = document.getElementById(btnId);
            if (btn) {
                btn.addEventListener('click', () => {
                    setActiveSubtab(tabName);
                });
            }
        });

        // Slider bindings
        const corrSlider = document.getElementById('rnaseq-corr-threshold');
        const corrVal = document.getElementById('rnaseq-corr-val');
        if (corrSlider && corrVal) {
            corrSlider.addEventListener('input', () => {
                state.corrThreshold = parseFloat(corrSlider.value);
                corrVal.textContent = state.corrThreshold.toFixed(2);
                updateCurrentView();
            });
        }

        const rfSlider = document.getElementById('rnaseq-rf-threshold');
        const rfVal = document.getElementById('rnaseq-rf-val');
        if (rfSlider && rfVal) {
            rfSlider.addEventListener('input', () => {
                state.rfThreshold = parseFloat(rfSlider.value);
                rfVal.textContent = state.rfThreshold.toFixed(3);
                updateCurrentView();
            });
        }

        // Time selector
        const timeSelect = document.getElementById('rnaseq-time-select');
        if (timeSelect) {
            timeSelect.addEventListener('change', () => {
                state.currentTimePoint = timeSelect.value;
                updateCurrentView();
            });
        }

        // Canvas actions
        const btnFitFloat = document.getElementById('btn-rnaseq-fit-float');
        if (btnFitFloat) {
            btnFitFloat.addEventListener('click', () => {
                if (state.cyInstance) {
                    state.cyInstance.fit();
                }
            });
        }

        const btnLayoutFloat = document.getElementById('btn-rnaseq-layout-float');
        if (btnLayoutFloat) {
            btnLayoutFloat.addEventListener('click', () => {
                runLayout();
            });
        }

        // Additional Canvas controls
        const btnZoomIn = document.getElementById('btn-rnaseq-zoom-in');
        if (btnZoomIn) {
            btnZoomIn.addEventListener('click', () => {
                if (state.cyInstance) {
                    state.cyInstance.zoom(state.cyInstance.zoom() * 1.2);
                }
            });
        }

        const btnZoomOut = document.getElementById('btn-rnaseq-zoom-out');
        if (btnZoomOut) {
            btnZoomOut.addEventListener('click', () => {
                if (state.cyInstance) {
                    state.cyInstance.zoom(state.cyInstance.zoom() / 1.2);
                }
            });
        }

        const btnExport = document.getElementById('btn-rnaseq-export');
        if (btnExport) {
            btnExport.addEventListener('click', () => {
                if (state.cyInstance) {
                    const pngContent = state.cyInstance.png({ bg: '#ffffff', full: true });
                    const link = document.createElement('a');
                    link.href = pngContent;
                    link.download = `heat_stress_subnetwork_${state.activeSubtab}.png`;
                    document.body.appendChild(link);
                    link.click();
                    document.body.removeChild(link);
                }
            });
        }

        // Inspector Close
        const btnInspectClose = document.getElementById('btn-rnaseq-inspect-close');
        if (btnInspectClose) {
            btnInspectClose.addEventListener('click', () => {
                const panel = document.getElementById('rnaseq-inspection-panel');
                if (panel) panel.classList.add('hidden');
                if (state.cyInstance) {
                    state.cyInstance.elements().removeClass('dimmed');
                    state.cyInstance.elements().removeClass('highlighted');
                }
            });
        }

        const btnRfbaRun = document.getElementById('rfba-run-btn');
        if (btnRfbaRun) {
            btnRfbaRun.addEventListener('click', async () => {
                const tfSelect = document.getElementById('rfba-tf-select');
                const modeSelect = document.getElementById('rfba-mode-select');
                const glucoseInput = document.getElementById('rfba-glucose-input');
                const biomassInput = document.getElementById('rfba-biomass-input');
                
                if (!tfSelect || !modeSelect || !glucoseInput || !biomassInput) return;
                
                const tf = tfSelect.value;
                const mode = modeSelect.value;
                const glucose = parseFloat(glucoseInput.value) || 100.0;
                const biomass = parseFloat(biomassInput.value) || 0.1;
                
                btnRfbaRun.disabled = true;
                btnRfbaRun.innerHTML = '<i class="fa-solid fa-spinner fa-spin"></i> Simulating...';
                
                try {
                    const tfPerturbations = {};
                    tfPerturbations[tf] = mode;
                    
                    const res = await window.simulationClient.runDynamicRFBA(tfPerturbations, glucose, biomass);
                    if (res.status === "error" || res.error) {
                        alert(`Simulation error: ${res.error || res.status}`);
                    } else {
                        plotRFBASimulation(res, tf, mode);
                    }
                } catch (err) {
                    console.error("rFBA run error:", err);
                    alert(`Failed to run rFBA simulation: ${err.message}`);
                } finally {
                    btnRfbaRun.disabled = false;
                    btnRfbaRun.innerHTML = '<i class="fa-solid fa-play"></i> Run Dynamic rFBA';
                }
            });
        }

        const poolSlider = document.getElementById('ecfba-pool-slider');
        const gdhSlider = document.getElementById('ecfba-gdh-slider');
        const lyscSlider = document.getElementById('ecfba-lysc-slider');
        const tempSlider = document.getElementById('ecfba-temp-slider');
        const productSelect = document.getElementById('ecfba-product-select');

        const updateSlidersText = () => {
            if (poolSlider) setElText('ecfba-pool-val', parseFloat(poolSlider.value).toFixed(3));
            if (gdhSlider) setElText('ecfba-gdh-val', parseFloat(gdhSlider.value).toFixed(1) + 'x');
            if (lyscSlider) setElText('ecfba-lysc-val', parseFloat(lyscSlider.value).toFixed(1) + 'x');
            if (tempSlider) setElText('ecfba-temp-val', parseFloat(tempSlider.value).toFixed(1) + '°C');
        };

        if (poolSlider) {
            poolSlider.addEventListener('input', updateSlidersText);
            poolSlider.addEventListener('change', () => renderECFBASimulation());
        }
        if (gdhSlider) {
            gdhSlider.addEventListener('input', updateSlidersText);
            gdhSlider.addEventListener('change', () => renderECFBASimulation());
        }
        if (lyscSlider) {
            lyscSlider.addEventListener('input', updateSlidersText);
            lyscSlider.addEventListener('change', () => renderECFBASimulation());
        }
        if (tempSlider) {
            tempSlider.addEventListener('input', updateSlidersText);
            tempSlider.addEventListener('change', () => renderECFBASimulation());
        }
        if (productSelect) {
            productSelect.addEventListener('change', () => renderECFBASimulation());
        }
        const calSelect = document.getElementById('ecfba-calibration-select');
        if (calSelect) {
            calSelect.addEventListener('change', () => {
                const disabled = calSelect.value !== 'none';
                if (gdhSlider) gdhSlider.disabled = disabled;
                if (lyscSlider) lyscSlider.disabled = disabled;
                renderECFBASimulation();
            });
        }
    }

    async function activate() {
        if (!state.data) {
            await fetchAnalysisData();
        } else {
            updateStatsBar();
            updateCurrentView();
        }
    }

    async function fetchAnalysisData() {
        if (state.loading) return;
        state.loading = true;
        
        const findingsText = document.getElementById('rnaseq-findings-text');
        if (findingsText) findingsText.innerHTML = '<i class="fa-solid fa-spinner fa-spin"></i> Fetching precomputed transcriptomic & network metrics...';

        try {
            // Fetch from API
            const response = await fetch('/api/analysis/rna-seq');
            if (!response.ok) {
                // Differentiate between "not published yet" (403/404) and actual errors
                if (response.status === 403 || response.status === 404 || response.status === 500) {
                    const body = await response.json().catch(() => ({}));
                    throw Object.assign(new Error('not_published'), { detail: body.detail || null, status: response.status });
                }
                throw new Error(`Server error: ${response.status} ${response.statusText}`);
            }
            const data = await response.json();
            state.data = data;
            
            console.log("Successfully loaded Heat Stress GRN data!");
            
            updateStatsBar();
            updateFindingsSummary();
            updateCurrentView();
        } catch (err) {
            console.warn("Heat Stress GRN not available:", err.message);
            if (findingsText) {
                if (err.message === 'not_published') {
                    findingsText.innerHTML = `
                        <div style="text-align:center; padding: 40px 20px; color: var(--text-secondary);">
                            <i class="fa-solid fa-hourglass-half" style="font-size:36px; color:#a5b4fc; margin-bottom:16px; display:block;"></i>
                            <strong style="font-size:14px; color:var(--text-primary); display:block; margin-bottom:8px;">Coming Soon</strong>
                            <p style="font-size:12px; line-height:1.6; max-width:380px; margin:0 auto;">
                                The heat stress transcriptomics &amp; GRN analysis dataset will be publicly released upon publication of the associated research paper. Stay tuned.
                            </p>
                        </div>`;
                } else {
                    findingsText.innerHTML = `<span style="color: #ef4444;"><i class="fa-solid fa-circle-exclamation"></i> Unable to load data: ${err.message}. Check that the backend is running.</span>`;
                }
            }
        } finally {
            state.loading = false;
        }
    }

    function updateStatsBar() {
        if (!state.data) return;
        const meta = state.data.metadata;
        
        // Update stats elements in index.html
        setElText('rnaseq-stat-genes', meta.total_expression_genes ? meta.total_expression_genes.toLocaleString() : '-');
        setElText('rnaseq-stat-tfs', meta.total_regulators ? meta.total_regulators.toLocaleString() : '-');
        setElText('rnaseq-stat-edges', meta.total_network_edges ? meta.total_network_edges.toLocaleString() : '-');
        
        const rewiredCount = state.data.rewired_edges ? state.data.rewired_edges.length : 0;
        setElText('rnaseq-stat-rewired', rewiredCount.toLocaleString());
    }

    function updateFindingsSummary() {
        const textEl = document.getElementById('rnaseq-findings-text');
        if (!textEl || !state.data) return;
        
        const hubs = state.data.hub_switching || [];
        const topStressHubs = hubs
            .filter(h => h.delta_degree > 0)
            .slice(0, 3)
            .map(h => h.tf_name)
            .join(', ');
            
        const topSubsys = state.data.metabolic_mapping?.top_subsystems?.[0]?.subsystem || "amino acid biosynthesis";
        const rewiredCount = state.data.rewired_edges ? state.data.rewired_edges.length : 0;

        textEl.innerHTML = `
            Heat stress induces massive regulatory rewiring of <strong>${rewiredCount} edges</strong>. 
            Primary heat shock hubs include <strong>${topStressHubs || 'sigH, sigB'}</strong>, which heavily redirect transcription 
            towards <strong>${topSubsys}</strong> and stress response subsystems.
        `;
    }

    function setActiveSubtab(tabName) {
        state.activeSubtab = tabName;
        
        // Toggle tab active class
        const subtabBtns = ['grn', 'dynamic', 'causal', 'enrich', 'time', 'rewire', 'metab', 'ecfba', 'figures'];
        subtabBtns.forEach(t => {
            const btn = document.getElementById(`btn-subtab-${t}`);
            const content = document.getElementById(`subtab-content-${t}`);
            
            if (btn) btn.classList.toggle('active', t === tabName);
            if (content) content.classList.toggle('hidden', t !== tabName);
        });

        updateCurrentView();
    }

    function updateCurrentView() {
        if (!state.data) return;

        switch (state.activeSubtab) {
            case 'grn':
                renderOverviewGRN();
                break;
            case 'dynamic':
                renderDynamicGRN();
                break;
            case 'causal':
                renderCausalGRN();
                break;
            case 'enrich':
                renderMotifEnrichment();
                break;
            case 'time':
                renderTimeResolved();
                break;
            case 'rewire':
                renderRewiring();
                break;
            case 'metab':
                renderMetabolicMapping();
                break;
            case 'ecfba':
                renderECFBASimulation();
                break;
            case 'mfa':
                renderMFAComparison();
                break;
            case 'figures':
                // Figures subtab uses static HTML files linked to download.
                break;
        }
    }

    // --- Tab Rendering Functions ---

    function renderOverviewGRN() {
        const edges = state.data.inferred_grn || [];
        const tbody = document.getElementById('rnaseq-grn-table-body');
        
        // Filter edges based on sliders
        const filteredEdges = edges.filter(edge => 
            Math.abs(edge.r_all) >= state.corrThreshold && 
            edge.rf_weight >= state.rfThreshold
        );

        tbody.innerHTML = '';
        if (filteredEdges.length === 0) {
            tbody.innerHTML = '<tr><td colspan="4" style="padding: 12px; text-align: center; color: var(--text-muted);">No edges match current thresholds.</td></tr>';
        } else {
            filteredEdges.slice(0, 100).forEach(edge => {
                const tr = document.createElement('tr');
                tr.style.borderBottom = '1px solid var(--border-color)';
                tr.style.cursor = 'pointer';
                tr.innerHTML = `
                    <td style="padding: 6px 8px; font-weight:600; color:var(--color-primary-accent);">${escapeHtml(edge.tf_name)}</td>
                    <td style="padding: 6px 8px;">${escapeHtml(edge.tg_name)}</td>
                    <td style="padding: 6px 8px; font-family: monospace;">${edge.r_all.toFixed(3)}</td>
                    <td style="padding: 6px 8px; font-family: monospace;">${edge.rf_weight.toFixed(4)}</td>
                `;
                tr.addEventListener('click', () => {
                    highlightNodeInNetwork(edge.tf_cgl);
                });
                tbody.appendChild(tr);
            });
        }

        buildAndDrawNetwork(filteredEdges.slice(0, 80));
    }

    function renderTimeResolved() {
        const timeData = state.data.time_resolved?.[state.currentTimePoint];
        const tbody = document.getElementById('rnaseq-time-table-body');
        
        setElText('rnaseq-time-deg-count', timeData?.deg_count || '0');
        setElText('rnaseq-time-edge-count', timeData?.active_edge_count || '0');

        tbody.innerHTML = '';
        if (!timeData || !timeData.top_regulators || timeData.top_regulators.length === 0) {
            tbody.innerHTML = '<tr><td colspan="2" style="padding: 12px; text-align: center; color: var(--text-muted);">No active regulators at this time point.</td></tr>';
        } else {
            timeData.top_regulators.forEach(reg => {
                const tr = document.createElement('tr');
                tr.style.borderBottom = '1px solid var(--border-color)';
                tr.innerHTML = `
                    <td style="padding: 6px 8px; font-weight:600;">${escapeHtml(reg.tf_name)}</td>
                    <td style="padding: 6px 8px; text-align: right; font-weight:bold; color:var(--color-activation);">${reg.active_targets} targets</td>
                `;
                tbody.appendChild(tr);
            });
        }

        const edges = timeData?.edges || [];
        buildAndDrawNetwork(edges.slice(0, 80), 'time');
    }

    function renderRewiring() {
        const edges = state.data.rewired_edges || [];
        const tbody = document.getElementById('rnaseq-rewire-table-body');
        
        // Count rewiring categories
        let gained = 0, lost = 0, inverted = 0, modulated = 0;
        edges.forEach(e => {
            if (e.type === 'gain') gained++;
            else if (e.type === 'loss') lost++;
            else if (e.type === 'inversion') inverted++;
            else modulated++;
        });

        setElText('rnaseq-rewire-gained', gained);
        setElText('rnaseq-rewire-lost', lost);
        setElText('rnaseq-rewire-inverted', inverted);
        setElText('rnaseq-rewire-modulated', modulated);

        tbody.innerHTML = '';
        if (edges.length === 0) {
            tbody.innerHTML = '<tr><td colspan="4" style="padding: 12px; text-align: center; color: var(--text-muted);">No rewired edges found.</td></tr>';
        } else {
            edges.slice(0, 100).forEach(edge => {
                const badgeColor = edge.type === 'gain' ? '#10b981' : (edge.type === 'loss' ? '#3b82f6' : '#a855f7');
                const tr = document.createElement('tr');
                tr.style.borderBottom = '1px solid var(--border-color)';
                tr.innerHTML = `
                    <td style="padding: 6px 8px; font-weight:600;">${escapeHtml(edge.tf_name)}</td>
                    <td style="padding: 6px 8px;">${escapeHtml(edge.tg_name)}</td>
                    <td style="padding: 6px 8px; font-family: monospace;">${edge.r_control.toFixed(2)} &rarr; ${edge.r_heat.toFixed(2)}</td>
                    <td style="padding: 6px 8px;"><span style="background: ${badgeColor}; color:white; padding:2px 6px; border-radius:4px; font-size:9.5px; font-weight:600;">${edge.type}</span></td>
                `;
                tbody.appendChild(tr);
            });
        }

        buildAndDrawNetwork(edges.slice(0, 60), 'rewire');
    }

    function renderHubSwitching() {
        const hubs = state.data.hub_switching || [];
        const tbody = document.getElementById('rnaseq-hub-table-body');

        tbody.innerHTML = '';
        if (hubs.length === 0) {
            tbody.innerHTML = '<tr><td colspan="4" style="padding: 12px; text-align: center; color: var(--text-muted);">No hub switches identified.</td></tr>';
        } else {
            hubs.forEach(hub => {
                const sign = hub.delta_degree > 0 ? '+' : '';
                const color = hub.delta_degree > 0 ? 'var(--color-activation)' : 'var(--color-repression)';
                const tr = document.createElement('tr');
                tr.style.borderBottom = '1px solid var(--border-color)';
                tr.innerHTML = `
                    <td style="padding: 6px 8px; font-weight:700; color:var(--color-primary-accent);">${escapeHtml(hub.tf_name)}</td>
                    <td style="padding: 6px 8px; text-align: right; font-family: monospace;">${hub.control_degree}</td>
                    <td style="padding: 6px 8px; text-align: right; font-family: monospace;">${hub.heat_degree}</td>
                    <td style="padding: 6px 8px; text-align: right; font-family: monospace; font-weight:bold; color: ${color};">${sign}${hub.delta_degree}</td>
                `;
                tbody.appendChild(tr);
            });
        }

        // Draw hub subnetwork
        const hubNames = new Set(hubs.slice(0, 5).map(h => h.tf_name));
        const edges = (state.data.inferred_grn || []).filter(e => hubNames.has(e.tf_name));
        buildAndDrawNetwork(edges.slice(0, 80));
    }

    function renderMetabolicMapping() {
        const metab = state.data.metabolic_mapping;
        const tbody = document.getElementById('rnaseq-metab-table-body');
        
        tbody.innerHTML = '';
        if (!metab || !metab.top_subsystems || metab.top_subsystems.length === 0) {
            tbody.innerHTML = '<tr><td colspan="2" style="padding: 12px; text-align: center; color: var(--text-muted);">No metabolic mappings available.</td></tr>';
        } else {
            metab.top_subsystems.forEach(sub => {
                const tr = document.createElement('tr');
                tr.style.borderBottom = '1px solid var(--border-color)';
                tr.innerHTML = `
                    <td style="padding: 6px 8px; font-weight:600;">${escapeHtml(sub.subsystem)}</td>
                    <td style="padding: 6px 8px; text-align: right; font-weight:bold; color: #8b5cf6;">${sub.edge_count} links</td>
                `;
                tbody.appendChild(tr);
            });
        }

        // Render coupled FBA growth and glutamate export flux over time
        const coupling = state.data.metabolic_coupling;
        const canvas = document.getElementById('rnaseq-metab-chart');
        if (coupling && canvas && window.Chart) {
            if (state.metabChartInstance) {
                state.metabChartInstance.destroy();
            }
            const ctx = canvas.getContext('2d');
            state.metabChartInstance = new Chart(ctx, {
                type: 'line',
                data: {
                    labels: coupling.time.map(t => `${t}h`),
                    datasets: [
                        {
                            label: 'Growth Rate (h-1)',
                            data: coupling.growth_rate,
                            borderColor: '#1e3a8a',
                            backgroundColor: 'rgba(30, 58, 138, 0.05)',
                            borderWidth: 1.8,
                            yAxisID: 'y_growth',
                            tension: 0.15
                        },
                        {
                            label: 'Glu Export (mmol/gDW/h)',
                            data: coupling.glutamate_export,
                            borderColor: '#8b5cf6',
                            borderDash: [3, 3],
                            backgroundColor: 'transparent',
                            borderWidth: 1.8,
                            yAxisID: 'y_flux',
                            tension: 0.15
                        }
                    ]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {
                        legend: { labels: { font: { size: 9 } } }
                    },
                    scales: {
                        x: { ticks: { font: { size: 8 } }, grid: { display: false } },
                        y_growth: {
                            type: 'linear',
                            position: 'left',
                            ticks: { font: { size: 8 }, color: '#1e3a8a' },
                            title: { display: true, text: 'Growth', font: { size: 9, weight: 'bold' } }
                        },
                        y_flux: {
                            type: 'linear',
                            position: 'right',
                            ticks: { font: { size: 8 }, color: '#8b5cf6' },
                            title: { display: true, text: 'Glutamate Export', font: { size: 9, weight: 'bold' } },
                            grid: { drawOnChartArea: false }
                        }
                    }
                }
            });
        }

        // Draw metabolic mapped submap
        const edges = metab?.mapped_edges_sample || [];
        buildAndDrawNetwork(edges.slice(0, 60), 'metab');
    }

    function renderDynamicGRN() {
        const dynamicData = state.data.dynamic_grn;
        if (!dynamicData) return;

        const select = document.getElementById('rnaseq-dynamic-gene-select');
        if (!select) return;
        const tgLocus = select.value;

        const traj = dynamicData.trajectories[tgLocus];
        const params = dynamicData.ode_parameters[tgLocus];

        const canvas = document.getElementById('rnaseq-dynamic-chart');
        if (canvas && window.Chart) {
            if (state.dynamicChartInstance) {
                state.dynamicChartInstance.destroy();
            }
            const ctx = canvas.getContext('2d');
            state.dynamicChartInstance = new Chart(ctx, {
                type: 'line',
                data: {
                    labels: Array.from({length: 25}, (_, i) => `${i}h`),
                    datasets: [
                        {
                            label: 'Actual RNA-Seq (Discrete Points)',
                            data: traj ? traj.actual : [],
                            borderColor: '#1e3a8a',
                            backgroundColor: '#1e3a8a',
                            borderWidth: 0,
                            pointRadius: 6,
                            pointHoverRadius: 8,
                            showLine: false
                        },
                        {
                            label: 'dGRN ODE Predicted',
                            data: traj ? traj.predicted : [],
                            borderColor: '#ef4444',
                            borderDash: [4, 4],
                            backgroundColor: 'transparent',
                            borderWidth: 1.8,
                            pointRadius: 0,
                            tension: 0.15
                        }
                    ]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {
                        legend: { labels: { font: { size: 9 } } }
                    },
                    scales: {
                        x: { ticks: { font: { size: 8 } }, grid: { display: false } },
                        y: {
                            ticks: { font: { size: 8 } },
                            title: { display: true, text: 'Expr (log2)', font: { size: 9, weight: 'bold' } }
                        }
                    }
                }
            });
        }

        const degEl = document.getElementById('rnaseq-dynamic-deg');
        const regListEl = document.getElementById('rnaseq-dynamic-reg-list');
        if (degEl && regListEl) {
            if (params) {
                degEl.textContent = params.degradation_rate.toFixed(3);
                regListEl.innerHTML = '';
                if (!params.regulators || params.regulators.length === 0) {
                    regListEl.innerHTML = '<div style="color:var(--text-muted); padding:2px;">No regulators fitted.</div>';
                } else {
                    params.regulators.forEach(reg => {
                        const div = document.createElement('div');
                        div.style.display = 'flex';
                        div.style.justify = 'space-between';
                        div.style.padding = '2px 0';
                        const sign = reg.weight > 0 ? '+' : '';
                        const color = reg.weight > 0 ? 'var(--color-activation)' : 'var(--color-repression)';
                        div.innerHTML = `
                            <span style="font-weight:600; color:var(--color-primary-accent); cursor:pointer;" onclick="window.heatStressGrn.highlightNode('${reg.tf_cgl}')">${escapeHtml(reg.tf_name)}</span>
                            <span style="font-family:monospace; font-weight:bold; color:${color};">${sign}${reg.weight.toFixed(3)}</span>
                        `;
                        regListEl.appendChild(div);
                    });
                }
            } else {
                degEl.textContent = '-';
                regListEl.innerHTML = '<div style="color:var(--text-muted); padding:2px;">No parameters available.</div>';
            }
        }

        if (params && params.regulators && params.regulators.length > 0) {
            const dynamicEdges = params.regulators.map(reg => ({
                tf_cgl: reg.tf_cgl,
                tf_name: reg.tf_name,
                tg_cgl: tgLocus,
                tg_name: traj ? traj.gene_name : tgLocus,
                r_all: reg.weight,
                is_known: true
            }));
            buildAndDrawNetwork(dynamicEdges);
        }
    }

    function renderCausalGRN() {
        const causalEdges = state.data.causal_grn || [];
        const tbody = document.getElementById('rnaseq-causal-table-body');
        if (!tbody) return;
        tbody.innerHTML = '';
        if (causalEdges.length === 0) {
            tbody.innerHTML = '<tr><td colspan="4" style="padding: 12px; text-align: center; color: var(--text-muted);">No causal edges identified.</td></tr>';
        } else {
            causalEdges.slice(0, 100).forEach(edge => {
                const isAct = edge.direction === 'activation';
                const color = isAct ? 'var(--color-activation)' : 'var(--color-repression)';
                const sign = isAct ? 'Activation (+)' : 'Repression (-)';
                const tr = document.createElement('tr');
                tr.style.borderBottom = '1px solid var(--border-color)';
                tr.style.cursor = 'pointer';
                tr.innerHTML = `
                    <td style="padding: 6px 8px; font-weight:600; color:var(--color-primary-accent);">${escapeHtml(edge.tf_name)}</td>
                    <td style="padding: 6px 8px;">${escapeHtml(edge.tg_name)}</td>
                    <td style="padding: 6px 8px; font-weight:bold; color:${color};">${sign}</td>
                    <td style="padding: 6px 8px; font-family: monospace; text-align: right;">${edge.p_value.toExponential(2)}</td>
                `;
                tr.addEventListener('click', () => {
                    highlightNodeInNetwork(edge.tf_cgl);
                });
                tbody.appendChild(tr);
            });
        }

        const topCausal = causalEdges.filter(e => e.is_causal).slice(0, 50);
        buildAndDrawNetwork(topCausal.map(e => ({
            tf_cgl: e.tf_cgl,
            tf_name: e.tf_name,
            tg_cgl: e.tg_cgl,
            tg_name: e.tg_name,
            r_all: e.r_correlation,
            is_known: true
        })));
    }

    function renderMotifEnrichment() {
        const enrichTime = document.getElementById('rnaseq-enrich-time-select').value;
        const enrichList = state.data.motif_enrichment?.[enrichTime] || [];
        const tbody = document.getElementById('rnaseq-enrich-table-body');
        if (!tbody) return;
        tbody.innerHTML = '';
        if (enrichList.length === 0) {
            tbody.innerHTML = '<tr><td colspan="3" style="padding: 12px; text-align: center; color: var(--text-muted);">No enriched TFs at this time.</td></tr>';
            const box = document.getElementById('rnaseq-enrich-logo-box');
            if (box) box.classList.add('hidden');
        } else {
            enrichList.forEach((item, idx) => {
                const tr = document.createElement('tr');
                tr.style.borderBottom = '1px solid var(--border-color)';
                tr.style.cursor = 'pointer';
                if (idx === 0) tr.classList.add('selected-row');
                tr.innerHTML = `
                    <td style="padding: 6px 8px; font-weight:700; color:var(--color-primary-accent);">${escapeHtml(item.tf_name)}</td>
                    <td style="padding: 6px 8px; text-align: right; font-family: monospace;">${item.fold_enrichment.toFixed(1)}x</td>
                    <td style="padding: 6px 8px; text-align: right; font-family: monospace;">${item.p_value.toExponential(2)}</td>
                `;
                tr.addEventListener('click', () => {
                    $(tbody).find('tr').removeClass('selected-row');
                    tr.classList.add('selected-row');
                    showTFLogo(item.tf_name);
                });
                tbody.appendChild(tr);
            });
            showTFLogo(enrichList[0].tf_name);
        }

        const enrichTFs = new Set(enrichList.slice(0, 5).map(e => e.tf_name));
        const edges = (state.data.inferred_grn || []).filter(e => enrichTFs.has(e.tf_name));
        buildAndDrawNetwork(edges.slice(0, 80));
    }

    async function showTFLogo(tfName) {
        const box = document.getElementById('rnaseq-enrich-logo-box');
        const label = document.getElementById('rnaseq-enrich-logo-tf');
        const canvas = document.getElementById('rnaseq-enrich-logo-canvas');
        if (!box || !label || !canvas) return;

        label.textContent = tfName;
        box.classList.remove('hidden');

        try {
            const response = await fetch(`/api/predict_motif?tf=${encodeURIComponent(tfName)}`);
            if (response.ok) {
                const data = await response.json();
                if (data && data.sequences && data.sequences.length > 0) {
                    if (window.renderMotifLogo) {
                        window.renderMotifLogo(canvas, data.sequences);
                    } else {
                        drawFallbackLogo(canvas, data.consensus || "ACGTACGT");
                    }
                } else {
                    drawFallbackLogo(canvas, data.consensus || "ACGTACGT");
                }
            } else {
                drawFallbackLogo(canvas, "ACGTACGT");
            }
        } catch (err) {
            drawFallbackLogo(canvas, "ACGTACGT");
        }
    }

    function drawFallbackLogo(canvas, consensus) {
        const ctx = canvas.getContext('2d');
        if (!ctx) return;
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        ctx.fillStyle = '#f8fafc';
        ctx.fillRect(0, 0, canvas.width, canvas.height);
        ctx.fillStyle = '#475569';
        ctx.font = 'bold 12px sans-serif';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText(`Consensus: ${consensus}`, canvas.width / 2, canvas.height / 2);
    }

    // --- Cytoscape Network Drawing ---

    function buildAndDrawNetwork(edges, mode = 'normal') {
        const cyContainer = document.getElementById('rna-seq-cy');
        if (!cyContainer || !window.cytoscape) return;

        // Compile unique nodes
        const nodesMap = new Map();
        const cyEdges = [];

        edges.forEach((edge, idx) => {
            const sourceId = edge.tf_cgl || edge.tf_name;
            const sourceName = edge.tf_name || sourceId;
            const targetId = edge.tg_cgl || edge.tg_cg || edge.tg_name;
            const targetName = edge.tg_name || targetId;

            // Add source node
            if (!nodesMap.has(sourceId)) {
                nodesMap.set(sourceId, {
                    id: sourceId,
                    label: sourceName,
                    type: 'TF'
                });
            }

            // Add target node
            if (!nodesMap.has(targetId)) {
                nodesMap.set(targetId, {
                    id: targetId,
                    label: targetName,
                    type: 'Target'
                });
            }

            // Add edge
            let relType = 'activation';
            let role = 'A';
            if (edge.r_all < 0 || edge.r_control < 0 || edge.r_heat < 0 || edge.r_correlation < 0) {
                relType = 'repression';
                role = 'R';
            }

            let weight = 0.5;
            if (edge.rf_weight !== undefined) {
                weight = edge.rf_weight * 10;
            }

            cyEdges.push({
                data: {
                    id: `rnaseq_edge_${idx}`,
                    source: sourceId,
                    target: targetId,
                    r_val: edge.r_all ?? edge.r_control ?? 0,
                    rf_weight: edge.rf_weight ?? 0,
                    regulationType: relType,
                    role: role,
                    is_known: edge.is_known ?? true,
                    confidenceScore: edge.rf_weight ?? Math.abs(edge.r_all ?? 0.5)
                }
            });
        });

        const cyNodes = Array.from(nodesMap.values()).map(n => ({
            data: {
                id: n.id,
                label: n.label,
                type: n.type
            }
        }));

        // Initialize Cytoscape
        if (state.cyInstance) {
            state.cyInstance.destroy();
        }

        state.cyInstance = window.cytoscape({
            container: cyContainer,
            elements: {
                nodes: cyNodes,
                edges: cyEdges
            },
            style: [
                {
                    selector: 'node',
                    style: {
                        'label': 'data(label)',
                        'text-valign': 'bottom',
                        'text-margin-y': '6px',
                        'color': '#0f172a',
                        'font-size': '11px',
                        'font-family': 'var(--font-primary)',
                        'background-color': '#f5f5f5',
                        'border-color': '#757575',
                        'border-width': '2px',
                        'width': '22px',
                        'height': '22px',
                        'transition-property': 'background-color, line-color, target-arrow-color, width, height, border-width',
                        'transition-duration': '0.2s'
                    }
                },
                {
                    selector: 'node[type="TF"]',
                    style: {
                        'background-color': '#e3f2fd',
                        'border-color': '#1976d2',
                        'width': '26px',
                        'height': '26px',
                        'shape': 'ellipse'
                    }
                },
                {
                    selector: 'node[type="sRNA"]',
                    style: {
                        'background-color': '#f3e5f5',
                        'border-color': '#8e24aa',
                        'width': '26px',
                        'height': '26px',
                        'shape': 'hexagon'
                    }
                },
                {
                    selector: 'node[type="Target"]',
                    style: {
                        'background-color': '#f5f5f5',
                        'border-color': '#757575',
                        'width': '22px',
                        'height': '22px',
                        'shape': 'ellipse'
                    }
                },
                {
                    selector: 'edge',
                    style: {
                        'width': (edge) => 1.2 + ((edge.data('confidenceScore') || 0.25) * 3.2),
                        'line-color': '#e65100',
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
                    selector: 'edge[regulationType="repression"]',
                    style: {
                        'line-color': '#d32f2f',
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
                    selector: 'edge[role="A"]',
                    style: {
                        'line-color': '#2e7d32',
                        'target-arrow-color': '#2e7d32',
                        'target-arrow-shape': 'triangle'
                    }
                },
                {
                    selector: 'edge[role="R"]',
                    style: {
                        'line-color': '#d32f2f',
                        'target-arrow-color': '#d32f2f',
                        'target-arrow-shape': 'tee'
                    }
                },
                {
                    selector: 'edge[is_known=false]',
                    style: {
                        'line-style': 'dashed',
                        'line-color': '#e65100',
                        'target-arrow-color': '#e65100'
                    }
                },
                {
                    selector: 'node.dimmed',
                    style: {
                        'opacity': 0.15,
                        'text-opacity': 0.15
                    }
                },
                {
                    selector: 'edge.dimmed',
                    style: {
                        'opacity': 0.05
                    }
                },
                {
                    selector: 'node.highlighted',
                    style: {
                        'border-width': '3.5px',
                        'border-color': 'var(--color-primary-accent)',
                        'font-weight': 'bold',
                        'font-size': '11px'
                    }
                }
            ],
            layout: {
                name: 'cose',
                animate: true,
                animationDuration: 450,
                nodeRepulsion: 45000,
                idealEdgeLength: 80
            }
        });

        // --- Interactive Events & Visual Upgrades ---

        // 1. Mouse Over Node (Highlight neighbors & show tooltip)
        state.cyInstance.on('mouseover', 'node', function(e) {
            const node = e.target;
            const neighborhood = node.neighborhood().add(node);
            
            state.cyInstance.elements().addClass('dimmed');
            neighborhood.removeClass('dimmed');
            neighborhood.addClass('highlighted');
            
            // Show Tooltip
            const tooltip = document.getElementById('rnaseq-canvas-tooltip');
            if (tooltip) {
                const poPos = e.renderedPosition;
                tooltip.style.left = (poPos.x + 12) + 'px';
                tooltip.style.top = (poPos.y - 12) + 'px';
                tooltip.innerHTML = `
                    <strong>${escapeHtml(node.data('label'))}</strong> (${node.id()})<br>
                    <span style="color:var(--text-muted);">Type:</span> <span style="font-weight:600; text-transform:uppercase;">${node.data('type')}</span><br>
                    <span style="color:var(--text-muted);">Degree in View:</span> <strong>${node.degree()}</strong>
                `;
                tooltip.classList.remove('hidden');
            }
        });

        state.cyInstance.on('mouseout', 'node', function(e) {
            state.cyInstance.elements().removeClass('dimmed');
            state.cyInstance.elements().removeClass('highlighted');
            const tooltip = document.getElementById('rnaseq-canvas-tooltip');
            if (tooltip) tooltip.classList.add('hidden');
        });

        // 2. Mouse Over Edge (Show interaction details tooltip)
        state.cyInstance.on('mouseover', 'edge', function(e) {
            const edge = e.target;
            const tooltip = document.getElementById('rnaseq-canvas-tooltip');
            if (tooltip) {
                const poPos = e.renderedPosition;
                tooltip.style.left = (poPos.x + 12) + 'px';
                tooltip.style.top = (poPos.y - 12) + 'px';
                
                const typeStr = edge.data('regulationType') === 'activation' ? '<span style="color:#2e7d32; font-weight:600;">Activation (+)</span>' : '<span style="color:#d32f2f; font-weight:600;">Repression (-)</span>';
                const rVal = edge.data('r_val');
                const rfVal = edge.data('rf_weight');
                
                tooltip.innerHTML = `
                    <strong>Regulation Link:</strong><br>
                    ${escapeHtml(edge.source().data('label'))} &rarr; ${escapeHtml(edge.target().data('label'))}<br>
                    Mode: ${typeStr}<br>
                    Correlation r: <span style="font-family:monospace; font-weight:600;">${rVal.toFixed(3)}</span><br>
                    RF Weight: <span style="font-family:monospace; font-weight:600;">${rfVal.toFixed(4)}</span>
                `;
                tooltip.classList.remove('hidden');
            }
        });

        state.cyInstance.on('mouseout', 'edge', function(e) {
            const tooltip = document.getElementById('rnaseq-canvas-tooltip');
            if (tooltip) tooltip.classList.add('hidden');
        });

        // 3. Node Tap (Select and populate inspector card)
        state.cyInstance.on('tap', 'node', function(evt){
            const node = evt.target;
            const locus = node.id();
            const name = node.data('label');
            const type = node.data('type');
            
            // Populate Inspector Panel
            const panel = document.getElementById('rnaseq-inspection-panel');
            const nameEl = document.getElementById('rnaseq-inspect-name');
            const locusEl = document.getElementById('rnaseq-inspect-locus');
            const typeEl = document.getElementById('rnaseq-inspect-type');
            const connEl = document.getElementById('rnaseq-inspect-connections');
            
            if (panel && nameEl && locusEl && typeEl && connEl) {
                nameEl.textContent = name || locus;
                locusEl.textContent = locus;
                typeEl.textContent = type;
                
                // Find all connected edges in current view
                const connectedEdges = node.connectedEdges();
                connEl.innerHTML = '';
                
                if (connectedEdges.length === 0) {
                    connEl.innerHTML = '<div style="color:var(--text-muted); padding: 4px; text-align:center;">No connections in view.</div>';
                } else {
                    connectedEdges.forEach(edge => {
                        const otherNode = edge.connectedNodes().not(node);
                        const otherName = otherNode.data('label') || otherNode.id();
                        const isSource = edge.source().id() === locus;
                        const direction = isSource ? `Targets ${otherName}` : `Regulated by ${otherName}`;
                        const typeClass = edge.data('regulationType') === 'activation' ? 'rnaseq-inspect-link-activation' : 'rnaseq-inspect-link-repression';
                        
                        const div = document.createElement('div');
                        div.className = `rnaseq-inspect-link-item ${typeClass}`;
                        div.innerHTML = `
                            <span>${direction}</span>
                            <span style="font-family:monospace; font-weight:600;">r: ${edge.data('r_val').toFixed(2)}</span>
                        `;
                        div.style.cursor = 'pointer';
                        div.addEventListener('click', (e) => {
                            e.stopPropagation();
                            highlightNodeInNetwork(otherNode.id());
                        });
                        connEl.appendChild(div);
                    });
                }
                panel.classList.remove('hidden');
            }
            
            // Trigger main view lookup if possible
            if (window.querySingleGene) {
                console.log(`Inspecting node: ${locus}`);
                const searchInput = document.querySelector('.gene-input');
                if (searchInput) searchInput.value = name || locus;
                window.querySingleGene(locus);
                if (window.toggleRightSidebar) {
                    window.toggleRightSidebar(true);
                }
            }
        });
    }

    function runLayout() {
        if (state.cyInstance) {
            state.cyInstance.layout({
                name: 'cose',
                animate: true,
                animationDuration: 400,
                nodeRepulsion: 45000,
                idealEdgeLength: 80
            }).run();
        }
    }

    function highlightNodeInNetwork(nodeId) {
        if (!state.cyInstance) return;
        const node = state.cyInstance.getElementById(nodeId);
        if (node.length > 0) {
            state.cyInstance.animate({
                fit: {
                    eles: node,
                    padding: 80
                }
            }, {
                duration: 350
            });
            node.flashClass('highlighted', 1000);
        }
    }

    async function renderECFBASimulation() {
        const poolSlider = document.getElementById('ecfba-pool-slider');
        const gdhSlider = document.getElementById('ecfba-gdh-slider');
        const lyscSlider = document.getElementById('ecfba-lysc-slider');
        const tempSlider = document.getElementById('ecfba-temp-slider');
        const calSelect = document.getElementById('ecfba-calibration-select');
        const productSelect = document.getElementById('ecfba-product-select');
        
        if (!poolSlider || !gdhSlider || !lyscSlider || !tempSlider || !productSelect) return;
        
        const poolLimit = parseFloat(poolSlider.value) || 0.129;
        const gdhLevel = parseFloat(gdhSlider.value) !== undefined ? parseFloat(gdhSlider.value) : 1.0;
        const lyscLevel = parseFloat(lyscSlider.value) !== undefined ? parseFloat(lyscSlider.value) : 1.0;
        const temperature = parseFloat(tempSlider.value) || 30.0;
        const product = productSelect.value || "growth";
        const calibrateTimepoint = calSelect && calSelect.value !== 'none' ? calSelect.value : null;
        
        const fluxDisplay = document.getElementById('ecfba-flux-display');
        if (fluxDisplay) {
            fluxDisplay.innerHTML = '<i class="fa-solid fa-spinner fa-spin" style="font-size: 12px;"></i>';
        }
        
        try {
            const perturbations = {
                "gdh": gdhLevel,
                "lysC": lyscLevel
            };
            
            const res = await window.simulationClient.runECFBA(poolLimit, perturbations, product, temperature, calibrateTimepoint);
            
            if (res.status === "success") {
                const suffix = product === 'growth' ? ' h-1' : ' mmol/gDW/h';
                if (fluxDisplay) {
                    fluxDisplay.style.color = '#10b981';
                    fluxDisplay.textContent = res.flux.toFixed(4) + suffix;
                }
                
                if (res.calibratedPerturbations) {
                    if (gdhSlider && res.calibratedPerturbations.gdh !== undefined) {
                        gdhSlider.value = res.calibratedPerturbations.gdh;
                        setElText('ecfba-gdh-val', res.calibratedPerturbations.gdh.toFixed(2) + 'x');
                    }
                    if (lyscSlider && res.calibratedPerturbations.lysC !== undefined) {
                        lyscSlider.value = res.calibratedPerturbations.lysC;
                        setElText('ecfba-lysc-val', res.calibratedPerturbations.lysC.toFixed(2) + 'x');
                    }
                }
                
                const usagePercent = Math.min(100, Math.max(0, Math.round((res.poolUsage / res.poolLimit) * 100)));
                setElText('ecfba-usage-percent', usagePercent + '%');
                setElText('ecfba-usage-abs', res.poolUsage.toFixed(4));
                setElText('ecfba-limit-abs', res.poolLimit.toFixed(3));
                
                const progressBar = document.getElementById('ecfba-usage-progress');
                if (progressBar) {
                    progressBar.style.width = usagePercent + '%';
                }
            } else {
                if (fluxDisplay) {
                    fluxDisplay.style.color = '#ef4444';
                    fluxDisplay.textContent = 'Infeasible';
                }
                setElText('ecfba-usage-percent', '0%');
                setElText('ecfba-usage-abs', '0.0000');
                setElText('ecfba-limit-abs', poolLimit.toFixed(3));
                
                const progressBar = document.getElementById('ecfba-usage-progress');
                if (progressBar) progressBar.style.width = '0%';
            }
            
            const warnBox = document.getElementById('ecfba-warning-box');
            if (warnBox) {
                if (res.warnings && res.warnings.length > 0) {
                    warnBox.style.display = 'block';
                    warnBox.innerHTML = '<i class="fa-solid fa-circle-exclamation"></i> ' + res.warnings.join('<br>');
                } else {
                    warnBox.style.display = 'none';
                }
            }
        } catch (err) {
            console.error("ec-FBA run error:", err);
            if (fluxDisplay) {
                fluxDisplay.style.color = '#ef4444';
                fluxDisplay.textContent = 'Error';
            }
            const warnBox = document.getElementById('ecfba-warning-box');
            if (warnBox) {
                warnBox.style.display = 'block';
                warnBox.innerHTML = '<i class="fa-solid fa-triangle-exclamation"></i> Calculation failed: ' + err.message;
            }
        }
    }

    function plotRFBASimulation(res, tf, mode) {
        const canvas = document.getElementById('rnaseq-metab-chart');
        if (!canvas || !window.Chart) return;
        
        if (state.metabChartInstance) {
            state.metabChartInstance.destroy();
        }
        
        const modeLabel = mode === 'normal' ? 'normal' : `${tf} ${mode}`;
        const ctx = canvas.getContext('2d');
        state.metabChartInstance = new Chart(ctx, {
            type: 'line',
            data: {
                labels: res.time.map(t => `${t}h`),
                datasets: [
                    {
                        label: `Growth (${modeLabel})`,
                        data: res.growth_rate,
                        borderColor: '#10b981',
                        backgroundColor: 'rgba(16, 185, 129, 0.05)',
                        borderWidth: 1.8,
                        yAxisID: 'y_growth',
                        tension: 0.15
                    },
                    {
                        label: `Glu Export (${modeLabel})`,
                        data: res.glutamate_export,
                        borderColor: '#8b5cf6',
                        borderDash: [3, 3],
                        backgroundColor: 'transparent',
                        borderWidth: 1.8,
                        yAxisID: 'y_flux',
                        tension: 0.15
                    },
                    {
                        label: 'Glucose (mM)',
                        data: res.glucose_concentration,
                        borderColor: '#f59e0b',
                        backgroundColor: 'transparent',
                        borderWidth: 1.2,
                        yAxisID: 'y_glucose',
                        tension: 0.15
                    }
                ]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                plugins: {
                    legend: { labels: { font: { size: 8 } } }
                },
                scales: {
                    x: { ticks: { font: { size: 8 } }, grid: { display: false } },
                    y_growth: {
                        type: 'linear',
                        position: 'left',
                        ticks: { font: { size: 8 }, color: '#10b981' },
                        title: { display: true, text: 'Growth', font: { size: 9, weight: 'bold' } }
                    },
                    y_flux: {
                        type: 'linear',
                        position: 'right',
                        ticks: { font: { size: 8 }, color: '#8b5cf6' },
                        title: { display: true, text: 'Glutamate Export', font: { size: 9, weight: 'bold' } },
                        grid: { drawOnChartArea: false }
                    },
                    y_glucose: {
                        type: 'linear',
                        position: 'right',
                        ticks: { font: { size: 8 }, color: '#f59e0b' },
                        title: { display: true, text: 'Glucose', font: { size: 9, weight: 'bold' } },
                        grid: { drawOnChartArea: false }
                    }
                }
            }
        });
    }

    // --- 13C-MFA Literature Comparison ---

    const PATHWAY_COLORS = {
        'Glycolysis': '#3b82f6',
        'PPP': '#10b981',
        'TCA Cycle': '#f59e0b',
        'Anaplerosis': '#8b5cf6',
        'Amino Acid Biosynthesis': '#ef4444'
    };

    let mfaScatterChart = null;

    async function renderMFAComparison() {
        const runBtn = document.getElementById('btn-run-mfa-comparison');
        if (runBtn) {
            runBtn.addEventListener('click', executeMFAComparison);
        }
    }

    async function executeMFAComparison() {
        const runBtn = document.getElementById('btn-run-mfa-comparison');
        if (runBtn) {
            runBtn.disabled = true;
            runBtn.innerHTML = '<i class="fa-solid fa-spinner fa-spin"></i> Running...';
        }

        const tbody = document.getElementById('mfa-comparison-tbody');
        const warnBox = document.getElementById('mfa-warning-box');
        if (tbody) tbody.innerHTML = '<tr><td colspan="5" style="padding:10px; text-align:center;"><i class="fa-solid fa-spinner fa-spin"></i> Comparing...</td></tr>';
        if (warnBox) { warnBox.style.display = 'none'; warnBox.textContent = ''; }

        try {
            const res = await window.simulationClient.runMFAComparison();

            if (res.status === 'success' || res.items.length > 0) {
                // Update stats
                setElText('mfa-pearson-val', res.pearson_r.toFixed(3));
                setElText('mfa-rmse-val', res.rmse.toFixed(2));
                setElText('mfa-meandev-val', res.mean_deviation_pct.toFixed(1) + '%');

                // Color pearson value
                const pearsonEl = document.getElementById('mfa-pearson-val');
                if (pearsonEl) {
                    pearsonEl.style.color = res.pearson_r >= 0.90 ? '#10b981' : res.pearson_r >= 0.70 ? '#f59e0b' : '#ef4444';
                }

                // Render comparison table
                if (tbody) {
                    tbody.innerHTML = '';
                    res.items.forEach(item => {
                        const devAbs = Math.abs(item.deviation_pct);
                        const devColor = devAbs <= 15 ? '#10b981' : devAbs <= 30 ? '#f59e0b' : '#ef4444';
                        const pathColor = PATHWAY_COLORS[item.pathway] || '#6b7280';
                        const row = document.createElement('tr');
                        row.style.borderBottom = '1px solid var(--border-color)';
                        row.innerHTML = `
                            <td style="padding: 3px 5px; font-weight:600; font-size:9.5px;">
                                <span style="display:inline-block; width:6px; height:6px; border-radius:50%; background:${pathColor}; margin-right:4px;"></span>
                                ${escapeHtml(item.reaction_id)}
                            </td>
                            <td style="padding: 3px 5px; font-size:9px; color:var(--text-muted);">${escapeHtml(item.pathway)}</td>
                            <td style="padding: 3px 5px; font-family:monospace; font-size:9.5px; text-align:right;">${item.mfa_flux.toFixed(2)} <span style="color:var(--text-muted)">±${item.mfa_std}</span></td>
                            <td style="padding: 3px 5px; font-family:monospace; font-size:9.5px; text-align:right; font-weight:600;">${item.sim_flux.toFixed(4)}</td>
                            <td style="padding: 3px 5px; text-align:right; font-size:9.5px; font-weight:700; color:${devColor};">${item.deviation_pct > 0 ? '+' : ''}${item.deviation_pct}%</td>
                        `;
                        tbody.appendChild(row);
                    });
                }

                // Render scatter plot
                const canvas = document.getElementById('mfa-scatter-canvas');
                if (canvas && typeof Chart !== 'undefined') {
                    if (mfaScatterChart) mfaScatterChart.destroy();
                    const ctx = canvas.getContext('2d');

                    const datasets = {};
                    res.items.forEach(item => {
                        const pathway = item.pathway;
                        if (!datasets[pathway]) {
                            datasets[pathway] = { label: pathway, data: [], backgroundColor: PATHWAY_COLORS[pathway] || '#6b7280', pointRadius: 5 };
                        }
                        datasets[pathway].data.push({ x: item.mfa_flux, y: item.sim_flux, label: item.reaction_id });
                    });

                    const maxVal = Math.max(...res.items.map(i => Math.max(i.mfa_flux, i.sim_flux))) * 1.1;

                    mfaScatterChart = new Chart(ctx, {
                        type: 'scatter',
                        data: { datasets: Object.values(datasets) },
                        options: {
                            responsive: true, maintainAspectRatio: false,
                            plugins: {
                                legend: { display: true, labels: { font: { size: 9 }, boxWidth: 10 } },
                                tooltip: { callbacks: { label: ctx => `${ctx.raw.label}: (MFA ${ctx.raw.x.toFixed(2)}, FBA ${ctx.raw.y.toFixed(2)})` } },
                                annotation: {
                                    annotations: {
                                        line1: { type: 'line', xMin: 0, xMax: maxVal, yMin: 0, yMax: maxVal, borderColor: '#9ca3af', borderWidth: 1, borderDash: [4, 4], label: { display: true, content: 'y = x', font: { size: 8 } } }
                                    }
                                }
                            },
                            scales: {
                                x: { min: 0, max: maxVal, title: { display: true, text: '¹³C-MFA (mmol/gDW/h)', font: { size: 9 } }, ticks: { font: { size: 9 } } },
                                y: { min: 0, max: maxVal, title: { display: true, text: 'FBA Simulated (mmol/gDW/h)', font: { size: 9 } }, ticks: { font: { size: 9 } } }
                            }
                        }
                    });
                }

                // Show warnings
                if (warnBox && res.warnings && res.warnings.length > 0) {
                    warnBox.innerHTML = '<strong>⚠ Mapping warnings:</strong><br>' + res.warnings.map(w => escapeHtml(w)).join('<br>');
                    warnBox.style.display = 'block';
                }

            } else {
                if (tbody) tbody.innerHTML = `<tr><td colspan="5" style="padding:12px; text-align:center; color:#ef4444;">${escapeHtml(res.warnings.join('; ') || 'Comparison failed.')}</td></tr>`;
            }
        } catch (err) {
            if (tbody) tbody.innerHTML = `<tr><td colspan="5" style="padding:12px; text-align:center; color:#ef4444;">Error: ${escapeHtml(err.message)}</td></tr>`;
        } finally {
            if (runBtn) {
                runBtn.disabled = false;
                runBtn.innerHTML = '<i class="fa-solid fa-rotate-right"></i> Re-run FBA vs. MFA Comparison';
            }
        }
    }

    // Helper functions
    function setElText(id, text) {
        const el = document.getElementById(id);
        if (el) el.textContent = text;
    }

    function escapeHtml(str) {
        if (!str) return '';
        return str
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;')
            .replace(/'/g, '&#039;');
    }

    // Expose window library
    window.heatStressGrn = {
        init: init,
        activate: activate,
        highlightNode: highlightNodeInNetwork
    };

    // Auto-run init on load
    $(document).ready(() => {
        init();
    });
})();
