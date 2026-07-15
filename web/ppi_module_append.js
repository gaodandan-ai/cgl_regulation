
/* ============================================================
   PPI Explorer Module v2 — Network / Shortest Path / Hub Ranking
   ============================================================ */

(function () {
    'use strict';

    // ── State ─────────────────────────────────────────────────────────────────
    let _ppiCy = null;
    let _ppiInitialized = false;
    let _ppiCurrentGene = null;
    let _ppiEdgeTooltip = null;
    let _ppiNodeData = {};
    let _ppiMode = 'network';   // 'network' | 'path' | 'hub'
    let _ppiHubData = [];       // cached hub rows
    let _ppiActiveChannels = new Set(['experimental','database','coexpression','textmining','neighborhood','cooccurrence','fusion']);

    const CHANNEL_COLORS = {
        experimental: '#16a34a',
        database:     '#2563eb',
        coexpression: '#7c3aed',
        textmining:   '#d97706',
        neighborhood: '#0891b2',
        cooccurrence: '#db2777',
        fusion:       '#ea580c',
    };

    function scoreClass(s) {
        return s >= 700 ? 'high' : s >= 400 ? 'medium' : 'low';
    }
    function dominantChannel(edge) {
        const chs = ['experimental','database','coexpression','textmining','neighborhood','cooccurrence','fusion'];
        let best = 'textmining', bv = 0;
        chs.forEach(function(ch){ const v = edge[ch]||0; if(v>bv){bv=v;best=ch;} });
        return best;
    }
    function edgeVisible(ed, minScore) {
        if ((ed.score||0) < minScore) return false;
        return _ppiActiveChannels.has(dominantChannel(ed));
    }

    // ── Mode switching ────────────────────────────────────────────────────────
    function switchPpiMode(mode) {
        _ppiMode = mode;
        document.querySelectorAll('.ppi-mode-tab').forEach(function(btn){
            btn.classList.toggle('active', btn.dataset.mode === mode);
        });
        document.querySelectorAll('.ppi-toolbar-panel').forEach(function(p){ p.classList.add('hidden'); });
        var tb = document.getElementById('ppi-toolbar-' + mode);
        if (tb) tb.classList.remove('hidden');

        var cyContainer   = document.getElementById('ppi-cy-container');
        var hubContainer  = document.getElementById('ppi-hub-table-container');
        var detailPanel   = document.getElementById('ppi-detail-panel');

        if (mode === 'hub') {
            if (cyContainer)  cyContainer.style.display  = 'none';
            if (hubContainer) { hubContainer.classList.remove('hidden'); hubContainer.style.display = 'block'; }
            if (detailPanel)  detailPanel.style.display  = 'none';
            loadHubRanking();
        } else {
            if (cyContainer)  cyContainer.style.display  = '';
            if (hubContainer) { hubContainer.classList.add('hidden'); hubContainer.style.display = 'none'; }
            if (detailPanel)  detailPanel.style.display  = '';
        }
    }

    // ── Cytoscape style ───────────────────────────────────────────────────────
    function ppiStyle() {
        return [
            { selector: 'node', style: {
                'label': 'data(name)',
                'text-valign': 'bottom', 'text-halign': 'center',
                'font-size': '10px', 'font-family': 'Inter,system-ui,sans-serif',
                'color': '#334155', 'text-margin-y': 4,
                'shape': 'ellipse',
                'width':  function(e){ return Math.max(22, Math.min(46, 22+(e.data('degree')||0)*2)); },
                'height': function(e){ return Math.max(22, Math.min(46, 22+(e.data('degree')||0)*2)); },
                'background-color': '#6366f1',
                'border-width': 2, 'border-color': '#fff',
                'text-wrap': 'wrap', 'text-max-width': '80px',
            }},
            { selector: 'node.ppi-seed', style: {
                'background-color': '#4f46e5', 'border-color': '#4f46e5', 'border-width': 3,
                'width':  function(e){ return Math.max(32, Math.min(52, 32+(e.data('degree')||0)*1.5)); },
                'height': function(e){ return Math.max(32, Math.min(52, 32+(e.data('degree')||0)*1.5)); },
                'font-weight': '700', 'font-size': '11px', 'color': '#1e1b4b', 'z-index': 10,
            }},
            { selector: 'node.ppi-partner',  style: { 'background-color': '#a5b4fc', 'border-color': '#fff' }},
            { selector: 'node.ppi-expanded', style: { 'background-color': '#10b981', 'border-color': '#059669' }},
            { selector: 'node:selected',     style: { 'border-color': '#f59e0b', 'border-width': 3, 'background-color': '#fde68a' }},
            { selector: 'edge', style: {
                'width': function(e){ return Math.max(1.2, Math.min(5, (e.data('score')||400)/220)); },
                'line-color': function(e){ return CHANNEL_COLORS[e.data('channel')]||'#94a3b8'; },
                'opacity': 0.75, 'curve-style': 'bezier',
                'target-arrow-shape': 'none', 'line-style': 'solid',
            }},
            { selector: 'edge.ppi-edge-conf-low',    style: { 'line-style': 'dotted', 'opacity': 0.4 }},
            { selector: 'edge.ppi-edge-conf-medium', style: { 'opacity': 0.6 }},
            { selector: 'edge.ppi-edge-conf-high',   style: { 'opacity': 0.85 }},
            { selector: '.ppi-dim',         style: { 'opacity': 0.08 }},
            { selector: '.ppi-highlighted', style: { 'opacity': 1 }},
            // Shortest path highlight
            { selector: '.ppi-path-node', style: { 'background-color': '#fde68a', 'border-color': '#f59e0b', 'border-width': 3 }},
            { selector: '.ppi-path-edge', style: { 'line-color': '#f59e0b', 'width': 4, 'opacity': 1 }},
        ];
    }

    // ── Build & render elements ───────────────────────────────────────────────
    function buildPpiElements(data, minScore) {
        _ppiNodeData = {};
        var nodes = data.nodes.map(function(n) {
            _ppiNodeData[n.id] = n;
            return { data: { id: n.id, name: n.name||n.id, label: n.name||n.id, degree: n.degree||0, is_seed: n.is_seed||false }, classes: n.is_seed ? 'ppi-seed' : 'ppi-partner' };
        });
        var edges = data.edges.filter(function(e){ return edgeVisible(e, minScore); }).map(function(e) {
            var ch = dominantChannel(e);
            return { data: Object.assign({}, e, {channel: ch}), classes: 'ppi-edge-conf-'+scoreClass(e.score||0)+' ppi-ch-'+ch };
        });
        return { nodes: nodes, edges: edges };
    }

    function renderPpiGraph(elements) {
        var container = document.getElementById('ppi-cy-container');
        if (!container) return;
        if (_ppiCy) { _ppiCy.destroy(); _ppiCy = null; }
        if (!elements.nodes.length) return;
        _ppiCy = window.cytoscape({
            container: container,
            elements: elements,
            style: ppiStyle(),
            layout: { name: 'cose', animate: true, animationDuration: 600,
                nodeRepulsion: function(){ return 8000; }, idealEdgeLength: function(){ return 80; },
                gravity: 0.3, randomize: true, fit: true, padding: 30 },
            minZoom: 0.1, maxZoom: 5,
        });
        bindPpiEvents();
    }

    // ── Tooltip ───────────────────────────────────────────────────────────────
    function showEdgeTooltip(evt, ed) {
        if (!_ppiEdgeTooltip) { _ppiEdgeTooltip = document.createElement('div'); _ppiEdgeTooltip.className = 'ppi-edge-tooltip'; document.body.appendChild(_ppiEdgeTooltip); }
        var chs = [{k:'experimental',l:'Exp',c:'#4ade80'},{k:'database',l:'DB',c:'#60a5fa'},{k:'coexpression',l:'CoExp',c:'#c084fc'},{k:'textmining',l:'Text',c:'#fbbf24'},{k:'neighborhood',l:'Nbhd',c:'#38bdf8'},{k:'cooccurrence',l:'CoOcc',c:'#f472b6'},{k:'fusion',l:'Fuse',c:'#fb923c'}];
        var bars = chs.map(function(c){ var v=ed[c.k]||0; if(!v)return''; return '<div style="display:flex;align-items:center;gap:4px;margin-top:2px;"><span style="min-width:34px;font-size:10px;color:#94a3b8;">'+c.l+'</span><div style="flex:1;height:5px;border-radius:3px;background:#1e293b;"><div style="width:'+Math.round(v/10)+'%;height:100%;background:'+c.c+';border-radius:3px;"></div></div><span style="font-size:10px;min-width:22px;text-align:right;">'+v+'</span></div>'; }).join('');
        _ppiEdgeTooltip.innerHTML = '<strong>Score: '+(ed.score||0)+'</strong>'+bars;
        _ppiEdgeTooltip.style.display = 'block';
        moveEdgeTooltip(evt);
    }
    function moveEdgeTooltip(evt) { if(!_ppiEdgeTooltip)return; var oe=evt.originalEvent||evt; _ppiEdgeTooltip.style.left=(oe.clientX+14)+'px'; _ppiEdgeTooltip.style.top=(oe.clientY+14)+'px'; }
    function hideEdgeTooltip() { if(_ppiEdgeTooltip)_ppiEdgeTooltip.style.display='none'; }

    // ── Right detail panel ────────────────────────────────────────────────────
    function updatePpiDetailPanel(nodeId) {
        var header = document.getElementById('ppi-detail-header');
        var nameEl = document.getElementById('ppi-detail-gene-name');
        var locusEl= document.getElementById('ppi-detail-gene-locus');
        var content= document.getElementById('ppi-detail-content');
        if (!content) return;
        var meta = _ppiNodeData[nodeId]||{};
        var name = (meta.name && meta.name!==nodeId) ? meta.name : nodeId;
        if (header)  header.style.display='';
        if (nameEl)  nameEl.textContent=name;
        if (locusEl) locusEl.textContent=nodeId;
        content.innerHTML='';
        var wrap=document.createElement('div'); wrap.id='ppi-detail-string-card'; content.appendChild(wrap);
        if (typeof renderStringPpiCard==='function') renderStringPpiCard(nodeId,wrap);
        else wrap.innerHTML='<div style="color:#94a3b8;font-size:12px;text-align:center;padding:20px;">Loading…</div>';
        if (_ppiCy) {
            _ppiCy.elements().addClass('ppi-dim').removeClass('ppi-highlighted');
            var n=_ppiCy.getElementById(nodeId); n.removeClass('ppi-dim').addClass('ppi-highlighted'); n.neighborhood().removeClass('ppi-dim').addClass('ppi-highlighted');
        }
    }

    // ── Cytoscape events ──────────────────────────────────────────────────────
    function bindPpiEvents() {
        if (!_ppiCy) return;
        _ppiCy.on('tap','node',function(evt){ updatePpiDetailPanel(evt.target.id()); });
        _ppiCy.on('dblclick dbltap','node',function(evt){
            var slider=document.getElementById('ppi-score-slider');
            expandPpiNode(evt.target.id(), slider?parseInt(slider.value):400);
        });
        _ppiCy.on('mouseover','edge',function(evt){ showEdgeTooltip(evt,evt.target.data()); });
        _ppiCy.on('mousemove','edge',function(evt){ moveEdgeTooltip(evt); });
        _ppiCy.on('mouseout','edge',function(){ hideEdgeTooltip(); });
        _ppiCy.on('tap',function(evt){
            if(evt.target!==_ppiCy)return;
            _ppiCy.elements().removeClass('ppi-dim ppi-highlighted');
            var h=document.getElementById('ppi-detail-header'), c=document.getElementById('ppi-detail-content');
            if(h)h.style.display='none';
            if(c)c.innerHTML='<div style="text-align:center;color:#94a3b8;padding:40px 0;font-size:12px;"><i class="fa-solid fa-hand-pointer" style="font-size:24px;margin-bottom:8px;display:block;color:#c7d2fe;"></i>Click a node to see its interactions</div>';
        });
    }

    // ── Node expansion ────────────────────────────────────────────────────────
    async function expandPpiNode(nodeId, minScore) {
        if (!_ppiCy) return;
        var existingIds = new Set(_ppiCy.nodes().map(function(n){return n.id();}));
        showPpiLoading(true);
        try {
            var resp = await fetch('/api/analysis/string_ppi/neighborhood?genes='+encodeURIComponent(nodeId)+'&min_score='+minScore+'&limit_per_gene=20');
            if (!resp.ok) return;
            var data = await resp.json();
            data.nodes.forEach(function(n){
                if(!existingIds.has(n.id)){
                    _ppiNodeData[n.id]=n;
                    _ppiCy.add({group:'nodes',data:{id:n.id,name:n.name||n.id,label:n.name||n.id,degree:n.degree||0,is_seed:false},classes:'ppi-partner ppi-expanded'});
                    existingIds.add(n.id);
                }
            });
            data.edges.forEach(function(e){
                if(existingIds.has(e.source)&&existingIds.has(e.target)&&!_ppiCy.getElementById(e.id).length&&edgeVisible(e,minScore)){
                    var ch=dominantChannel(e);
                    _ppiCy.add({group:'edges',data:Object.assign({},e,{channel:ch}),classes:'ppi-edge-conf-'+scoreClass(e.score||0)+' ppi-ch-'+ch});
                }
            });
            _ppiCy.layout({name:'cose',animate:true,animationDuration:400,fit:false,nodeRepulsion:function(){return 6000;}}).run();
            updatePpiStats();
        } catch(err){ console.warn('[PPI]expand:',err); }
        finally { showPpiLoading(false); }
    }

    // ── Stats ─────────────────────────────────────────────────────────────────
    function updatePpiStats() {
        var el=document.getElementById('ppi-stats-label');
        if(!el||!_ppiCy)return;
        el.textContent=_ppiCy.nodes().length+' nodes · '+_ppiCy.edges(':visible').length+' edges';
    }
    function showPpiLoading(show) {
        var el=document.getElementById('ppi-loading'), emp=document.getElementById('ppi-empty-state');
        if(el)el.classList.toggle('hidden',!show);
        if(show&&emp)emp.style.display='none';
    }

    // ── Network fetch & render ────────────────────────────────────────────────
    async function fetchAndRenderPpi(gene) {
        if (!gene) return;
        _ppiCurrentGene = gene;
        showPpiLoading(true);
        var slider=document.getElementById('ppi-score-slider'), score=slider?parseInt(slider.value):400;
        try {
            var resp=await fetch('/api/analysis/string_ppi/neighborhood?genes='+encodeURIComponent(gene)+'&min_score='+score+'&limit_per_gene=30');
            if(!resp.ok)throw new Error('HTTP '+resp.status);
            var data=await resp.json();
            if(!data.nodes||!data.nodes.length){
                showPpiLoading(false);
                var emp=document.getElementById('ppi-empty-state'); if(emp)emp.style.display='flex';
                return;
            }
            var elements=buildPpiElements(data,score);
            renderPpiGraph(elements);
            var sl=document.getElementById('ppi-stats-label');
            if(sl)sl.textContent=data.total_nodes+' nodes · '+data.total_edges+' edges  ·  STRING v'+((data.string_meta&&data.string_meta.version)||'12');
        } catch(err){ console.error('[PPI]fetch:',err); }
        finally { showPpiLoading(false); }
    }

    // ── Channel filter ────────────────────────────────────────────────────────
    function applyChannelFilter() {
        if(!_ppiCy)return;
        var sl=document.getElementById('ppi-score-slider'), ms=sl?parseInt(sl.value):400;
        _ppiCy.edges().forEach(function(e){ e.style('display',edgeVisible(e.data(),ms)?'element':'none'); });
        updatePpiStats();
    }

    // ── ① Shortest Path ───────────────────────────────────────────────────────
    async function findShortestPath() {
        var src=(document.getElementById('ppi-path-source')||{}).value||'';
        var tgt=(document.getElementById('ppi-path-target')||{}).value||'';
        var scoreEl=document.getElementById('ppi-path-score');
        var label=document.getElementById('ppi-path-result-label');
        if(!src.trim()||!tgt.trim()){if(label)label.textContent='Please enter both genes.';return;}
        var ms=scoreEl?parseInt(scoreEl.value):400;
        showPpiLoading(true);
        if(label)label.textContent='Searching…';
        try {
            var resp=await fetch('/api/analysis/string_ppi/shortest_path?source='+encodeURIComponent(src.trim())+'&target='+encodeURIComponent(tgt.trim())+'&min_score='+ms+'&max_hops=6');
            var data=await resp.json();
            if(!data.found){
                showPpiLoading(false);
                if(label)label.innerHTML='<span style="color:#ef4444;">'+data.message+'</span>';
                return;
            }
            // Build elements (path nodes + their full neighborhoods for context)
            _ppiNodeData={};
            var nodes=data.nodes.map(function(n){_ppiNodeData[n.id]=n;return{data:{id:n.id,name:n.name||n.id,label:n.name||n.id,degree:1,is_seed:n.is_seed||false},classes:n.is_seed?'ppi-seed ppi-path-node':'ppi-partner ppi-path-node'};});
            var edges=data.edges.map(function(e){var ch=dominantChannel(e);return{data:Object.assign({},e,{channel:ch}),classes:'ppi-path-edge ppi-ch-'+ch};});
            renderPpiGraph({nodes:nodes,edges:edges});
            if(label)label.innerHTML='<span style="color:#16a34a;font-weight:700;">Found! Path length: '+data.hops+' hop'+(data.hops>1?'s':'')+'</span>';
            var sl=document.getElementById('ppi-stats-label');
            if(sl)sl.textContent=data.hops+' hops · '+(data.nodes.length)+' nodes';
        } catch(err){ if(label)label.textContent='Error: '+err.message; }
        finally { showPpiLoading(false); }
    }

    // ── ② Hub Ranking ─────────────────────────────────────────────────────────
    async function loadHubRanking() {
        var sel=document.getElementById('ppi-hub-score-select');
        var ms=sel?parseInt(sel.value):700;
        var tbody=document.getElementById('ppi-hub-tbody');
        var meta=document.getElementById('ppi-hub-meta-label');
        if(!tbody)return;
        tbody.innerHTML='<tr><td colspan="11" style="text-align:center;padding:30px;color:#94a3b8;"><i class="fa-solid fa-spinner fa-spin"></i> Loading…</td></tr>';
        try {
            var resp=await fetch('/api/analysis/string_ppi/hub_ranking?min_score='+ms+'&limit=50');
            var data=await resp.json();
            _ppiHubData=data.hubs||[];
            if(meta)meta.textContent=data.total_genes+' genes with ≥1 partner at score '+ms;
            renderHubTable(_ppiHubData);
        } catch(err){ if(tbody)tbody.innerHTML='<tr><td colspan="11" style="color:#ef4444;padding:20px;">Error loading hub data</td></tr>'; }
    }

    function renderHubTable(rows) {
        var tbody=document.getElementById('ppi-hub-tbody');
        if(!tbody)return;
        var maxDegree=rows.length?rows[0].degree:1;
        var filterVal=(document.getElementById('ppi-hub-search')||{}).value||'';
        var fl=filterVal.toLowerCase();
        var filtered=fl?rows.filter(function(r){return r.name.toLowerCase().includes(fl)||r.gene.includes(fl);}):rows;
        if(!filtered.length){tbody.innerHTML='<tr><td colspan="11" style="text-align:center;padding:30px;color:#94a3b8;">No matching proteins</td></tr>';return;}
        tbody.innerHTML=filtered.map(function(r,i){
            var pct=Math.round((r.degree/maxDegree)*100);
            var partners=(r.top_partners||[]).slice(0,5).map(function(p){return '<span style="font-size:10px;padding:1px 6px;border-radius:10px;background:#eef2ff;color:#4338ca;margin:1px;display:inline-block;">'+p+'</span>';}).join('');
            return '<tr>'+
                '<td style="font-weight:700;color:var(--text-secondary);">'+(i+1)+'</td>'+
                '<td style="font-weight:700;color:#1e293b;">'+r.name+'</td>'+
                '<td style="color:#94a3b8;font-size:11px;">'+r.gene+'</td>'+
                '<td><div class="ppi-hub-degree-bar"><span style="font-weight:700;color:#6366f1;min-width:28px;">'+r.degree+'</span><div style="flex:1;max-width:100px;"><div class="ppi-hub-degree-fill" style="width:'+pct+'%;"></div></div></div></td>'+
                '<td style="text-align:center;color:var(--text-secondary);">'+r.avg_score+'</td>'+
                '<td style="text-align:center;color:#16a34a;font-weight:600;">'+r.experimental+'</td>'+
                '<td style="text-align:center;color:#2563eb;font-weight:600;">'+r.database+'</td>'+
                '<td style="text-align:center;color:#7c3aed;font-weight:600;">'+r.coexpression+'</td>'+
                '<td style="text-align:center;color:#d97706;font-weight:600;">'+r.textmining+'</td>'+
                '<td>'+partners+'</td>'+
                '<td style="text-align:center;"><button onclick="window._ppiViewGene(\''+r.gene+'\')" style="padding:3px 10px;border-radius:6px;border:1px solid #c7d2fe;background:#eef2ff;color:#4338ca;font-size:10.5px;cursor:pointer;font-weight:600;">Network</button></td>'+
            '</tr>';
        }).join('');
    }

    // Global accessor for hub table "Network" buttons
    window._ppiViewGene = function(gene) {
        switchPpiMode('network');
        var inp=document.getElementById('ppi-search-input'); if(inp)inp.value=gene;
        fetchAndRenderPpi(gene);
    };

    // ── ③ Multi-gene subnetwork (handled via fetchAndRenderPpi with comma input) ──

    // ── One-time init ─────────────────────────────────────────────────────────
    window.initPpiExplorer = function () {
        if (_ppiInitialized) return;
        _ppiInitialized = true;

        // Mode tabs
        document.querySelectorAll('.ppi-mode-tab').forEach(function(btn){
            btn.addEventListener('click', function(){ switchPpiMode(btn.dataset.mode); });
        });

        // Network tab
        var si=document.getElementById('ppi-search-input'), sb=document.getElementById('ppi-search-btn');
        function doSearch(){ var g=si?si.value.trim():''; if(g)fetchAndRenderPpi(g); }
        if(sb)sb.addEventListener('click',doSearch);
        if(si)si.addEventListener('keydown',function(e){if(e.key==='Enter')doSearch();});

        document.querySelectorAll('.ppi-example-btn').forEach(function(btn){
            btn.addEventListener('click',function(){
                var g=btn.dataset.gene; if(si)si.value=g;
                switchPpiMode('network'); fetchAndRenderPpi(g);
            });
        });

        var scoreSlider=document.getElementById('ppi-score-slider'), scoreVal=document.getElementById('ppi-score-val');
        if(scoreSlider){
            scoreSlider.addEventListener('input',function(){if(scoreVal)scoreVal.textContent=scoreSlider.value;applyChannelFilter();});
            scoreSlider.addEventListener('change',function(){if(_ppiCurrentGene)fetchAndRenderPpi(_ppiCurrentGene);});
        }

        document.querySelectorAll('.ppi-ch-pill').forEach(function(pill){
            pill.addEventListener('click',function(){
                var ch=pill.dataset.channel;
                if(_ppiActiveChannels.has(ch)){_ppiActiveChannels.delete(ch);pill.classList.remove('active');}
                else{_ppiActiveChannels.add(ch);pill.classList.add('active');}
                applyChannelFilter();
            });
        });

        var lb=document.getElementById('ppi-layout-btn');
        if(lb)lb.addEventListener('click',function(){if(_ppiCy)_ppiCy.layout({name:'cose',animate:true,animationDuration:500,fit:true,nodeRepulsion:function(){return 8000;}}).run();});
        var fb=document.getElementById('ppi-fit-btn');
        if(fb)fb.addEventListener('click',function(){if(_ppiCy)_ppiCy.fit(undefined,30);});

        // Shortest Path tab
        var pathBtn=document.getElementById('ppi-path-btn');
        if(pathBtn)pathBtn.addEventListener('click',findShortestPath);
        var ps=document.getElementById('ppi-path-source'), pt=document.getElementById('ppi-path-target');
        function tryPath(e){if(e.key==='Enter')findShortestPath();}
        if(ps)ps.addEventListener('keydown',tryPath);
        if(pt)pt.addEventListener('keydown',tryPath);
        var pathScore=document.getElementById('ppi-path-score'), pathScoreVal=document.getElementById('ppi-path-score-val');
        if(pathScore)pathScore.addEventListener('input',function(){if(pathScoreVal)pathScoreVal.textContent=pathScore.value;});

        // Hub Ranking tab
        var hubLoad=document.getElementById('ppi-hub-load-btn');
        if(hubLoad)hubLoad.addEventListener('click',loadHubRanking);
        var hubSel=document.getElementById('ppi-hub-score-select');
        if(hubSel)hubSel.addEventListener('change',loadHubRanking);
        var hubSearch=document.getElementById('ppi-hub-search');
        if(hubSearch)hubSearch.addEventListener('input',function(){renderHubTable(_ppiHubData);});
    };

}());
