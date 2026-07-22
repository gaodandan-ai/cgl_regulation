/** Cross-module engineering target priority explorer. */
(function () {
    const state = {initialized:false, targets:[], total:0, selected:null};
    const esc = value => String(value ?? '').replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;').replace(/'/g,'&#039;');
    const fmt = (value, digits=3) => value == null ? '—' : Number(value).toFixed(digits);
    const strategyLabels = {
        dynamic_tuning_only:'Dynamic tuning only', careful_titration:'Careful titration',
        multi_stress_control_node:'Multi-stress control node',
        metabolic_intervention_candidate:'Metabolic intervention',
        context_specific_candidate:'Context-specific candidate',
    };
    const moduleLabels = {
        iron_regulon_v1:'Iron', oxygen_regulon_v1:'Oxygen', carbon_regulon_v1:'Carbon',
        nitrogen_regulon_v1:'Nitrogen', stress_regulon_v1:'Stress',
    };

    function shell(container) {
        container.innerHTML = `<style>
            #intervention-priority-root{color:#0f172a;font-family:var(--font-primary,Inter,sans-serif)}
            .ip-header{display:flex;justify-content:space-between;align-items:flex-start;gap:16px;margin-bottom:14px}.ip-header h2{font-size:20px;margin:0}.ip-header p{font-size:11.5px;color:#64748b;max-width:920px;line-height:1.5;margin:5px 0 0}
            .ip-panel{background:#fff;border:1px solid #dbe3ee;border-radius:10px;box-shadow:0 2px 8px rgba(15,23,42,.04)}
            .ip-filters{display:grid;grid-template-columns:minmax(220px,1.6fr) minmax(190px,1fr) 90px 100px 110px auto;gap:9px;padding:12px;margin-bottom:12px;align-items:end}.ip-field label{display:block;font-size:9.5px;font-weight:800;color:#475569;margin-bottom:4px;text-transform:uppercase}.ip-field input,.ip-field select{width:100%;height:34px;border:1px solid #cbd5e1;border-radius:6px;padding:0 8px;font-size:11px;background:#fff}.ip-check{display:flex;align-items:center;gap:6px;height:34px;font-size:10px;color:#475569}.ip-check input{width:auto}.ip-apply{height:34px;border:0;border-radius:6px;background:#4338ca;color:#fff;font-size:11px;font-weight:800;padding:0 15px;cursor:pointer}
            .ip-kpis{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin-bottom:12px}.ip-kpi{padding:10px 12px}.ip-kpi span{display:block;font-size:9px;text-transform:uppercase;color:#64748b}.ip-kpi strong{display:block;font-size:19px;margin-top:4px}.ip-layout{display:grid;grid-template-columns:minmax(720px,1.8fr) minmax(310px,.8fr);gap:12px}.ip-table-wrap{max-height:620px;overflow:auto}.ip-table{width:100%;border-collapse:collapse;font-size:10.5px}.ip-table th{position:sticky;top:0;background:#f8fafc;z-index:1;text-align:left;padding:7px;color:#475569;border-bottom:1px solid #cbd5e1}.ip-table td{padding:7px;border-bottom:1px solid #eef2f7;vertical-align:middle}.ip-table tr{cursor:pointer}.ip-table tr:hover{background:#f8fafc}.ip-gene{border:0;background:none;color:#0369a1;font-weight:800;padding:0;cursor:pointer}.ip-grade{display:inline-flex;width:20px;height:20px;border-radius:50%;align-items:center;justify-content:center;font-weight:900}.ip-grade-A{background:#dcfce7;color:#166534}.ip-grade-B{background:#dbeafe;color:#1d4ed8}.ip-grade-C{background:#fef3c7;color:#92400e}.ip-grade-D{background:#e2e8f0;color:#475569}.ip-bar{width:72px;height:6px;background:#e2e8f0;border-radius:4px;overflow:hidden}.ip-bar i{display:block;height:100%;background:#4f46e5}.ip-badge{display:inline-block;padding:2px 6px;border-radius:10px;background:#ede9fe;color:#5b21b6;font-size:9px;font-weight:700}.ip-essential{background:#fee2e2;color:#991b1b}.ip-detail{padding:13px;min-height:300px}.ip-detail h3{font-size:16px;margin:0 0 3px}.ip-detail .ip-sub{font-size:10px;color:#64748b;margin-bottom:12px}.ip-score-grid{display:grid;grid-template-columns:1fr 1fr;gap:7px;margin:10px 0}.ip-score-box{background:#f8fafc;border-radius:6px;padding:8px}.ip-score-box span{font-size:8.5px;color:#64748b;text-transform:uppercase}.ip-score-box strong{display:block;font-size:15px;margin-top:2px}.ip-module{border-top:1px solid #e2e8f0;padding:8px 0;font-size:10px}.ip-module strong{display:flex;justify-content:space-between}.ip-note{margin-top:12px;background:#fff7ed;border:1px solid #fed7aa;border-radius:8px;padding:10px;color:#9a3412;font-size:10px;line-height:1.5}.ip-loading,.ip-error{padding:34px;text-align:center;color:#64748b}.ip-error{color:#b91c1c}
            @media(max-width:1150px){.ip-filters{grid-template-columns:1fr 1fr 1fr}.ip-layout{grid-template-columns:1fr}.ip-kpis{grid-template-columns:1fr 1fr}}
        </style>
        <div class="ip-header"><div><h2><i class="fa-solid fa-bullseye" style="color:#4338ca"></i> Cross-module target priorities</h2><p>Rank targets by regulatory evidence, systems impact, engineering tractability and risk across iron, oxygen, carbon, nitrogen and stress modules. Rankings guide follow-up; they do not prescribe a perturbation without a product objective.</p></div><button class="secondary-btn" id="ip-back"><i class="fa-solid fa-arrow-left"></i> Gene Explorer</button></div>
        <div class="ip-panel ip-filters">
            <div class="ip-field"><label>Gene or product</label><input id="ip-query" placeholder="katA, cg0310, catalase…"></div>
            <div class="ip-field"><label>Strategy</label><select id="ip-strategy"><option value="">All strategies</option>${Object.entries(strategyLabels).map(([v,l])=>`<option value="${v}">${l}</option>`).join('')}</select></div>
            <div class="ip-field"><label>Min modules</label><input id="ip-modules" type="number" min="1" max="5" value="2"></div>
            <div class="ip-field"><label>Max risk</label><input id="ip-risk" type="number" min="0" max="1" step="0.05" value="1"></div>
            <div class="ip-field"><label>Evidence grade</label><select id="ip-grade"><option value="">All grades</option><option>A</option><option>B</option><option>C</option><option>D</option></select></div>
            <div><label class="ip-check"><input id="ip-essential" type="checkbox" checked> Include known essential</label><button class="ip-apply" id="ip-apply">Apply</button></div>
        </div>
        <div id="ip-results"><div class="ip-panel ip-loading"><i class="fa-solid fa-spinner fa-spin"></i> Ranking targets…</div></div>
        <div class="ip-note"><strong>Safety:</strong> “not in curated essential set” does not prove non-essentiality. Deletion, overexpression and dynamic-control decisions require product-specific flux analysis and experimental validation.</div>`;
    }

    function params() {
        const out = new URLSearchParams({
            min_modules:document.getElementById('ip-modules').value || '1',
            max_risk:document.getElementById('ip-risk').value || '1', limit:'150',
            include_known_essential:String(document.getElementById('ip-essential').checked),
        });
        const q=document.getElementById('ip-query').value.trim(), strategy=document.getElementById('ip-strategy').value, grade=document.getElementById('ip-grade').value;
        if(q)out.set('q',q);if(strategy)out.set('strategy',strategy);if(grade)out.set('grade',grade);return out;
    }

    async function load() {
        const results=document.getElementById('ip-results');results.innerHTML='<div class="ip-panel ip-loading"><i class="fa-solid fa-spinner fa-spin"></i> Ranking targets…</div>';
        try{const response=await fetch(`/api/intervention-targets?${params()}`);if(!response.ok)throw new Error(`Target API failed (${response.status})`);const payload=await response.json();state.targets=payload.targets||[];state.total=payload.total||0;render();}
        catch(error){results.innerHTML=`<div class="ip-panel ip-error">${esc(error.message)}</div>`;}
    }

    function render() {
        const aCount=state.targets.filter(x=>x.evidence_grade==='A').length, multi=state.targets.filter(x=>x.module_count>=3).length, essential=state.targets.filter(x=>x.essentiality_status==='known_essential').length;
        document.getElementById('ip-results').innerHTML=`<div class="ip-kpis"><div class="ip-panel ip-kpi"><span>Matching targets</span><strong>${state.total}</strong></div><div class="ip-panel ip-kpi"><span>A-grade shown</span><strong>${aCount}</strong></div><div class="ip-panel ip-kpi"><span>≥3 modules shown</span><strong>${multi}</strong></div><div class="ip-panel ip-kpi"><span>Known essential shown</span><strong>${essential}</strong></div></div>
        <div class="ip-layout"><section class="ip-panel ip-table-wrap"><table class="ip-table"><thead><tr><th>Target</th><th>Grade</th><th>Modules</th><th>Evidence</th><th>Impact</th><th>Risk</th><th>Priority</th><th>Strategy</th></tr></thead><tbody>${state.targets.map(row).join('')||'<tr><td colspan="8">No targets match.</td></tr>'}</tbody></table></section><aside class="ip-panel ip-detail" id="ip-detail">Select a target to inspect its score decomposition and module evidence.</aside></div>`;
        document.querySelectorAll('tr[data-locus]').forEach(tr=>tr.addEventListener('click',()=>detail(tr.dataset.locus)));
        document.querySelectorAll('.ip-gene').forEach(button=>button.addEventListener('click',event=>{event.stopPropagation();window.setActiveWorkflowEntry?.('gene');window.querySingleGene?.(button.dataset.locus);}));
        if(state.targets.length)detail(state.selected&&state.targets.some(x=>x.target_locus===state.selected)?state.selected:state.targets[0].target_locus);
    }

    function row(x){return `<tr data-locus="${esc(x.target_locus)}"><td><button class="ip-gene" data-locus="${esc(x.target_locus)}">${esc(x.target_name||x.target_locus)}</button><br><small>${esc(x.target_locus)}</small></td><td><span class="ip-grade ip-grade-${x.evidence_grade}">${x.evidence_grade}</span></td><td>${x.module_count}<br><small>${x.modules.map(m=>moduleLabels[m]||m).join(', ')}</small></td><td>${scoreBar(x.evidence_score)}</td><td>${scoreBar(x.systems_impact_score)}</td><td>${fmt(x.risk_score)}</td><td><strong>${fmt(x.priority_score)}</strong></td><td><span class="ip-badge ${x.essentiality_status==='known_essential'?'ip-essential':''}">${esc(strategyLabels[x.strategy_class]||x.strategy_class)}</span></td></tr>`;}
    function scoreBar(v){return `<strong>${fmt(v)}</strong><div class="ip-bar"><i style="width:${Math.max(0,Math.min(100,Number(v)*100))}%"></i></div>`;}

    async function detail(locus){state.selected=locus;const panel=document.getElementById('ip-detail');if(!panel)return;panel.innerHTML='<div class="ip-loading">Loading detail…</div>';try{const response=await fetch(`/api/intervention-targets/${encodeURIComponent(locus)}`);if(!response.ok)throw new Error('Detail unavailable');const x=await response.json();panel.innerHTML=`<h3>${esc(x.target_name||x.target_locus)}</h3><div class="ip-sub">${esc(x.target_locus)} · ${esc(x.product||'No product annotation')}</div><span class="ip-badge ${x.essentiality_status==='known_essential'?'ip-essential':''}">${esc(strategyLabels[x.strategy_class]||x.strategy_class)}</span><div class="ip-score-grid"><div class="ip-score-box"><span>Priority</span><strong>${fmt(x.priority_score)}</strong></div><div class="ip-score-box"><span>Risk</span><strong>${fmt(x.risk_score)}</strong></div><div class="ip-score-box"><span>Evidence</span><strong>${fmt(x.evidence_score)}</strong></div><div class="ip-score-box"><span>Impact</span><strong>${fmt(x.systems_impact_score)}</strong></div><div class="ip-score-box"><span>Tractability</span><strong>${fmt(x.engineering_tractability_score)}</strong></div><div class="ip-score-box"><span>Pathways</span><strong>${x.pathway_count}</strong></div></div><h4 style="font-size:11px;margin:12px 0 3px">Module evidence</h4>${(x.module_evidence||[]).map(m=>`<div class="ip-module"><strong><span>${esc(moduleLabels[m.module_run_id]||m.module_run_id)}</span><span>${fmt(m.mean_score)}</span></strong><div>${m.condition_count} conditions · ${m.supported_count} supported · ${m.significant_context_count} FDR contexts</div></div>`).join('')}<div class="ip-sub" style="margin-top:10px">Proteomics: ${x.proteomics_detected?'detected by conservative mapping':'not conservatively mapped'} · Essentiality: ${esc(x.essentiality_status)}</div>`;}catch(error){panel.innerHTML=`<div class="ip-error">${esc(error.message)}</div>`;}}

    function bind(){document.getElementById('ip-back').addEventListener('click',()=>window.setActiveWorkflowEntry?.('gene'));document.getElementById('ip-apply').addEventListener('click',load);document.getElementById('ip-query').addEventListener('keydown',e=>{if(e.key==='Enter')load();});}
    window.InterventionPriorityView={async init(containerId='intervention-priority-root'){const container=document.getElementById(containerId);if(!container)return;if(!state.initialized){shell(container);bind();state.initialized=true;}await load();}};
})();
