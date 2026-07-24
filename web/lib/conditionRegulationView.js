/** Interactive multi-module condition-specific regulation explorer. */
(function () {
    const state = {
        initialized: false, runs: [], conditions: [],
        runId: 'iron_regulon_v1', comparisonId: '', tf: 'HrrA',
        supportState: '', minScore: 0.70, network: null,
    };
    const tfByRun = {
        iron_regulon_v1: ['HrrA', 'DtxR'],
        oxygen_regulon_v1: ['ArnR', 'GlxR', 'HrrA'],
        carbon_regulon_v1: ['GlxR', 'RamA', 'RamB', 'SugR'],
        nitrogen_regulon_v1: ['AmtR', 'ArgR', 'LtbR'],
        stress_regulon_v1: ['SigH', 'OxyR', 'MtrA'],
    };
    const runLabels = {
        iron_regulon_v1: 'Iron / heme — HrrA & DtxR',
        oxygen_regulon_v1: 'Oxygen limitation — ArnR, GlxR & HrrA',
        carbon_regulon_v1: 'Carbon sources — GlxR, RamA, RamB & SugR',
        nitrogen_regulon_v1: 'Nitrogen & amino acids — AmtR, ArgR & LtbR',
        stress_regulon_v1: 'Oxidative & envelope stress — SigH, OxyR & MtrA',
    };
    const esc = value => String(value ?? '')
        .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;').replace(/'/g, '&#039;');
    const fmt = (value, digits = 3) => value == null ? '—' : Number(value).toFixed(digits);
    const stateLabel = {
        condition_supported: 'Supported', direction_conflict: 'Direction conflict',
        weak_context_support: 'Weak context', insufficient_dynamic_data: 'Insufficient data',
    };
    const stateClass = {
        condition_supported: 'cr-supported', direction_conflict: 'cr-conflict',
        weak_context_support: 'cr-weak', insufficient_dynamic_data: 'cr-insufficient',
    };

    function shell(container) {
        container.innerHTML = `
            <style>
                #condition-regulation-root{color:#0f172a;font-family:var(--font-primary,Inter,sans-serif)}
                .cr-header{display:flex;align-items:flex-start;justify-content:space-between;gap:16px;margin-bottom:14px}
                .cr-header h2{font-size:20px;margin:0}.cr-header p{font-size:11.5px;color:#64748b;margin:5px 0 0;max-width:900px;line-height:1.5}
                .cr-panel{background:#fff;border:1px solid #dbe3ee;border-radius:10px;box-shadow:0 2px 8px rgba(15,23,42,.04)}
                .cr-filters{display:grid;grid-template-columns:minmax(220px,1.2fr) minmax(280px,2fr) 135px 180px 115px auto;gap:10px;padding:12px;margin-bottom:12px;align-items:end}
                .cr-field label{display:block;font-size:10px;font-weight:700;color:#475569;margin-bottom:4px;text-transform:uppercase;letter-spacing:.03em}
                .cr-field select,.cr-field input{width:100%;height:34px;border:1px solid #cbd5e1;border-radius:6px;background:#fff;padding:0 8px;font-size:11px;color:#0f172a}
                .cr-refresh{height:34px;padding:0 14px;border:0;border-radius:6px;background:#0f766e;color:#fff;font-size:11px;font-weight:700;cursor:pointer}
                .cr-kpis{display:grid;grid-template-columns:repeat(5,minmax(120px,1fr));gap:10px;margin-bottom:12px}.cr-kpi{padding:10px 12px}.cr-kpi span{display:block;font-size:9.5px;color:#64748b;text-transform:uppercase}.cr-kpi strong{display:block;font-size:18px;margin-top:4px}.cr-kpi small{display:block;font-size:9.5px;color:#64748b;margin-top:3px}
                .cr-main{display:grid;grid-template-columns:minmax(420px,1.15fr) minmax(440px,1fr);gap:12px;min-height:430px}.cr-section-title{padding:10px 12px;border-bottom:1px solid #e2e8f0;font-size:11px;font-weight:800;display:flex;justify-content:space-between}.cr-network{height:390px}.cr-table-wrap{max-height:390px;overflow:auto}.cr-table{width:100%;border-collapse:collapse;font-size:10.5px}.cr-table th{position:sticky;top:0;background:#f8fafc;z-index:1;text-align:left;color:#475569;padding:7px;border-bottom:1px solid #cbd5e1}.cr-table td{padding:7px;border-bottom:1px solid #eef2f7}.cr-gene{border:0;background:none;color:#0369a1;font-weight:700;cursor:pointer;padding:0}.cr-badge{display:inline-block;padding:2px 6px;border-radius:10px;font-size:9px;font-weight:700}.cr-supported{background:#dcfce7;color:#166534}.cr-conflict{background:#fee2e2;color:#991b1b}.cr-weak{background:#fef3c7;color:#92400e}.cr-insufficient{background:#e2e8f0;color:#475569}.cr-score{font-family:ui-monospace,monospace;font-weight:800}.cr-note{margin-top:12px;padding:10px 12px;background:#eff6ff;border:1px solid #bfdbfe;border-radius:8px;color:#1e3a8a;font-size:10.5px;line-height:1.5}.cr-loading,.cr-error{padding:34px;text-align:center;color:#64748b}.cr-error{color:#b91c1c}
                @media(max-width:1100px){.cr-filters{grid-template-columns:1fr 1fr}.cr-kpis{grid-template-columns:repeat(2,1fr)}.cr-main{grid-template-columns:1fr}.cr-network{height:330px}}
            </style>
            <div class="cr-header"><div><h2><i class="fa-solid fa-flask-vial" style="color:#0f766e"></i> Condition-specific regulation</h2>
                <p>Compare inferred TF activity, target response, evidence completeness and direction consistency across environmental contrasts. Scores are contextual support values, not validation probabilities.</p></div>
                <button class="secondary-btn" id="cr-back"><i class="fa-solid fa-arrow-left"></i> Gene Explorer</button></div>
            <div class="cr-panel cr-filters">
                <div class="cr-field"><label for="cr-run">Analysis module</label><select id="cr-run"></select></div>
                <div class="cr-field"><label for="cr-condition">Condition contrast</label><select id="cr-condition"></select></div>
                <div class="cr-field"><label for="cr-tf">Transcription factor</label><select id="cr-tf"></select></div>
                <div class="cr-field"><label for="cr-state">Evidence state</label><select id="cr-state"><option value="">All states</option><option value="condition_supported">Condition supported</option><option value="direction_conflict">Direction conflict</option><option value="weak_context_support">Weak context</option><option value="insufficient_dynamic_data">Insufficient dynamic data</option></select></div>
                <div class="cr-field"><label for="cr-score">Minimum score</label><input id="cr-score" type="number" min="0" max="1" step="0.05" value="0.70"></div>
                <button class="cr-refresh" id="cr-refresh"><i class="fa-solid fa-rotate"></i> Apply</button>
            </div>
            <div id="cr-results"><div class="cr-panel cr-loading"><i class="fa-solid fa-spinner fa-spin"></i> Loading condition scores…</div></div>
            <div class="cr-note" id="cr-note"></div>`;
    }

    function renderRunControls() {
        document.getElementById('cr-run').innerHTML = state.runs.map(run =>
            `<option value="${esc(run.run_id)}" ${run.run_id === state.runId ? 'selected' : ''}>${esc(runLabels[run.run_id] || run.scope)}</option>`
        ).join('');
        const names = tfByRun[state.runId] || [];
        if (!names.includes(state.tf)) state.tf = names[0] || '';
        document.getElementById('cr-tf').innerHTML = names.map(name =>
            `<option ${name === state.tf ? 'selected' : ''}>${esc(name)}</option>`
        ).join('') + `<option value="">All TFs</option>`;
        const notes = {
            oxygen_regulon_v1: '<strong>Interpretation:</strong> ArnR, GlxR and HrrA activities are signed target projections and are not independent validation.',
            carbon_regulon_v1: '<strong>Interpretation:</strong> GlxR, RamA, RamB and SugR activities are signed target projections. The weakly annotated RamA iModulon is retained in the database but is not treated as independent TF activity.',
            nitrogen_regulon_v1: '<strong>Interpretation:</strong> AmtR and ArgR use regulon-oriented ICA activities; LtbR uses a signed target projection. LysG is withheld from quantitative scoring because current clear-direction coverage is insufficient.',
            stress_regulon_v1: '<strong>Interpretation:</strong> SigH, OxyR and MtrA activities are signed target projections. MprA has no mapped entry in the current ATCC 13032 layer; SigM/SigE are withheld because their imported edges lack regulatory direction.',
            iron_regulon_v1: '<strong>Interpretation:</strong> DtxR activity is ICA-derived and oriented with ΔdtxR contrasts; HrrA is target-projection-derived.',
        };
        document.getElementById('cr-note').innerHTML = `${notes[state.runId] || ''} Direction conflicts are capped at 0.55; missing dynamic evidence is capped at 0.65.`;
    }

    async function loadRuns() {
        const response = await fetch('/api/condition-regulation/runs');
        if (!response.ok) throw new Error(`Module request failed (${response.status})`);
        state.runs = (await response.json()).runs || [];
        if (!state.runs.some(run => run.run_id === state.runId)) state.runId = state.runs[0]?.run_id || '';
        renderRunControls();
    }

    async function loadConditions(reset = false) {
        const response = await fetch(`/api/condition-regulation/conditions?${new URLSearchParams({run_id: state.runId})}`);
        if (!response.ok) throw new Error(`Condition request failed (${response.status})`);
        state.conditions = (await response.json()).conditions || [];
        if (reset || !state.conditions.some(item => item.comparison_id === state.comparisonId)) {
            const preferred = state.runId === 'oxygen_regulon_v1'
                ? state.conditions.find(item => /anaerob 30 min/i.test(item.condition_label))
                : state.runId === 'carbon_regulon_v1'
                    ? state.conditions.find(item => /lactate/i.test(item.condition_label) && /glucose/i.test(item.condition_label) && !/ramA/i.test(item.condition_label))
                    : state.runId === 'nitrogen_regulon_v1'
                        ? state.conditions.find(item => /glutamine.*vs.*w\/o \(NH4\)2SO4 and urea/i.test(item.condition_label))
                        : state.runId === 'stress_regulon_v1'
                            ? state.conditions.find(item => /hydroperoxide/i.test(item.condition_label))
                            : state.conditions.find(item => /hrra/i.test(item.condition_label) && /FeSO4/i.test(item.condition_label));
            state.comparisonId = preferred?.comparison_id || state.conditions[0]?.comparison_id || '';
        }
        document.getElementById('cr-condition').innerHTML = state.conditions.map(item =>
            `<option value="${esc(item.comparison_id)}" ${item.comparison_id === state.comparisonId ? 'selected' : ''}>${item.has_fdr_signal ? '● ' : ''}${esc(item.condition_label)}</option>`
        ).join('');
    }

    async function loadResults() {
        const results = document.getElementById('cr-results');
        if (!state.comparisonId) { results.innerHTML = '<div class="cr-panel cr-error">No scored condition contrasts are available.</div>'; return; }
        results.innerHTML = '<div class="cr-panel cr-loading"><i class="fa-solid fa-spinner fa-spin"></i> Loading activity and edge evidence…</div>';
        const common = {run_id: state.runId, comparison_id: state.comparisonId};
        const edgeParams = new URLSearchParams({...common, min_score: String(state.minScore), limit: '100'});
        const summaryParams = new URLSearchParams(common);
        if (state.tf) { edgeParams.set('tf', state.tf); summaryParams.set('tf', state.tf); }
        if (state.supportState) edgeParams.set('state', state.supportState);
        try {
            const [summaryResponse, edgeResponse] = await Promise.all([
                fetch(`/api/condition-regulation/summary?${summaryParams}`),
                fetch(`/api/condition-regulation/edges?${edgeParams}`),
            ]);
            if (!summaryResponse.ok || !edgeResponse.ok) throw new Error('The condition API returned an error.');
            const summaries = (await summaryResponse.json()).summaries || [];
            const edgePayload = await edgeResponse.json();
            renderResults(summaries, edgePayload.edges || [], edgePayload.total || 0);
        } catch (error) { results.innerHTML = `<div class="cr-panel cr-error">${esc(error.message)}</div>`; }
    }

    function renderResults(summaries, edges, total) {
        const results = document.getElementById('cr-results');
        const primary = summaries[0] || {};
        const activityText = summaries.length ? summaries.map(s => `${esc(s.tf_name)} ${fmt(s.activity_score)}`).join(' · ') : '—';
        const qValues = summaries.filter(s => s.empirical_q_value != null).map(s => Number(s.empirical_q_value));
        const bestQ = qValues.length ? Math.min(...qValues) : null;
        const supported = edges.filter(edge => edge.support_state === 'condition_supported').length;
        const conflicts = edges.filter(edge => edge.support_state === 'direction_conflict').length;
        results.innerHTML = `<div class="cr-kpis">
            <div class="cr-panel cr-kpi"><span>TF activity</span><strong style="font-size:14px">${activityText}</strong><small>${esc(primary.activity_method || 'No activity estimate')}</small></div>
            <div class="cr-panel cr-kpi"><span>Best empirical q</span><strong>${fmt(bestQ, 4)}</strong><small>BH-FDR within this module</small></div>
            <div class="cr-panel cr-kpi"><span>Matching edges</span><strong>${total}</strong><small>Showing top ${edges.length}</small></div>
            <div class="cr-panel cr-kpi"><span>Supported shown</span><strong style="color:#15803d">${supported}</strong><small>Score and dynamic response pass</small></div>
            <div class="cr-panel cr-kpi"><span>Conflicts shown</span><strong style="color:#b91c1c">${conflicts}</strong><small>Direction-inconsistent evidence</small></div></div>
            <div class="cr-main"><section class="cr-panel"><div class="cr-section-title"><span>Condition-specific network</span><span>${esc(primary.condition_label || '')}</span></div><div id="cr-network" class="cr-network"></div></section>
            <section class="cr-panel"><div class="cr-section-title"><span>Ranked regulatory edges</span><span>score ≥ ${state.minScore.toFixed(2)}</span></div><div class="cr-table-wrap"><table class="cr-table"><thead><tr><th>TF → target</th><th>Role</th><th>Expression</th><th>Replicate</th><th>Score</th><th>State</th></tr></thead><tbody>${edges.map(edgeRow).join('') || '<tr><td colspan="6" style="text-align:center;padding:24px">No edges match these filters.</td></tr>'}</tbody></table></div></section></div>`;
        results.querySelectorAll('.cr-gene').forEach(button => button.addEventListener('click', () => {
            window.setActiveWorkflowEntry?.('gene'); window.querySingleGene?.(button.dataset.locus);
        }));
        drawNetwork(edges);
    }

    function edgeRow(edge) {
        const css = stateClass[edge.support_state] || 'cr-insufficient';
        return `<tr><td><strong>${esc(edge.tf_name)}</strong> → <button class="cr-gene" data-locus="${esc(edge.target_locus)}">${esc(edge.target_name || edge.target_locus)}</button><br><small style="color:#94a3b8">${esc(edge.target_locus)}</small></td><td>${esc(edge.regulation_role)}</td><td>${fmt(edge.target_expression_mean)}</td><td>${fmt(edge.replicate_consistency, 2)}</td><td class="cr-score">${fmt(edge.condition_score)}</td><td><span class="cr-badge ${css}">${esc(stateLabel[edge.support_state] || edge.support_state)}</span></td></tr>`;
    }

    function drawNetwork(edges) {
        const container = document.getElementById('cr-network');
        if (!container) return;
        if (state.network) { state.network.destroy(); state.network = null; }
        if (!window.cytoscape || !edges.length) { container.innerHTML = '<div class="cr-loading">No network edges match the active filters.</div>'; return; }
        const elements = [], seen = new Set();
        edges.slice(0, 70).forEach(edge => {
            const tfId = `tf:${edge.tf_locus}`, targetId = `gene:${edge.target_locus}`;
            if (!seen.has(tfId)) { elements.push({data:{id:tfId,label:edge.tf_name,type:'tf'}}); seen.add(tfId); }
            if (!seen.has(targetId)) { elements.push({data:{id:targetId,label:edge.target_name || edge.target_locus,type:'gene',locus:edge.target_locus,state:edge.support_state}}); seen.add(targetId); }
            elements.push({data:{id:`${tfId}>${targetId}`,source:tfId,target:targetId,role:edge.regulation_role,state:edge.support_state,score:edge.condition_score}});
        });
        state.network = window.cytoscape({container, elements, minZoom: 0.15, maxZoom: 2.0, textureOnViewport: true, hideEdgesOnViewport: elements.length > 100, pixelRatio: 'auto', style:[
            {selector:'node',style:{'label':'data(label)','font-size':8,'text-valign':'center','text-halign':'center','width':28,'height':28,'background-color':'#94a3b8','color':'#334155'}},
            {selector:'node[type="tf"]',style:{'shape':'round-rectangle','width':58,'height':34,'background-color':'#0f766e','color':'#fff','font-weight':700}},
            {selector:'node[state="condition_supported"]',style:{'background-color':'#22c55e'}},{selector:'node[state="direction_conflict"]',style:{'background-color':'#ef4444'}},{selector:'node[state="weak_context_support"]',style:{'background-color':'#f59e0b'}},
            {selector:'edge',style:{'width':'mapData(score,0,1,1,4)','curve-style':'bezier','target-arrow-shape':'triangle','line-color':'#8b5cf6','target-arrow-color':'#8b5cf6','opacity':.72}},
            {selector:'edge[role="A"]',style:{'line-color':'#16a34a','target-arrow-color':'#16a34a'}},{selector:'edge[role="R"]',style:{'line-color':'#dc2626','target-arrow-color':'#dc2626'}},{selector:'edge[state="direction_conflict"]',style:{'line-style':'dashed','opacity':.45}},
        ], layout:{name:'cose',animate:false,fit:true,padding:elements.length <= 3 ? 120 : 24,nodeRepulsion:7000,idealEdgeLength:75}});
        state.network.on('tap', 'node[type="gene"]', event => { window.setActiveWorkflowEntry?.('gene'); window.querySingleGene?.(event.target.data('locus')); });
    }

    function bind() {
        document.getElementById('cr-back').addEventListener('click', () => window.setActiveWorkflowEntry?.('gene'));
        document.getElementById('cr-run').addEventListener('change', async event => {
            state.runId = event.target.value; state.comparisonId = ''; renderRunControls();
            await loadConditions(true); await loadResults();
        });
        document.getElementById('cr-condition').addEventListener('change', event => { state.comparisonId = event.target.value; loadResults(); });
        document.getElementById('cr-tf').addEventListener('change', event => { state.tf = event.target.value; loadResults(); });
        document.getElementById('cr-state').addEventListener('change', event => { state.supportState = event.target.value; loadResults(); });
        document.getElementById('cr-refresh').addEventListener('click', () => {
            state.comparisonId = document.getElementById('cr-condition').value; state.tf = document.getElementById('cr-tf').value;
            state.supportState = document.getElementById('cr-state').value;
            state.minScore = Math.max(0, Math.min(1, Number(document.getElementById('cr-score').value) || 0)); loadResults();
        });
    }

    window.ConditionRegulationView = {async init(containerId = 'condition-regulation-root') {
        const container = document.getElementById(containerId); if (!container) return;
        if (!state.initialized) {
            shell(container); bind();
            try { await loadRuns(); await loadConditions(); state.initialized = true; }
            catch (error) { document.getElementById('cr-results').innerHTML = `<div class="cr-panel cr-error">${esc(error.message)}</div>`; return; }
        }
        await loadResults();
    }};
})();
