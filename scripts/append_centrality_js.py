"""scripts/append_centrality_js.py — appends centrality functions to app.js"""

JS_CODE = r"""
// ── Network Centrality Analysis (pre-computed) ─────────────────────────────────
let _centralityLoaded = false;

async function loadPrecomputedCentrality() {
    const tbody = document.getElementById('topo-centrality-tbody');
    if (!tbody) return;
    tbody.innerHTML = `<tr><td colspan="9" style="text-align:center;padding:20px;"><i class="fa-solid fa-spinner fa-spin"></i> Loading centrality data…</td></tr>`;

    try {
        const resp = await fetch('/api/network/centrality?limit=50&tfs_only=true');
        if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
        const data = await resp.json();
        _centralityLoaded = true;

        const tfs = data.top_tfs || [];
        const maxScore = tfs.length > 0 ? tfs[0].importance : 1;

        tbody.innerHTML = tfs.map((tf, i) => {
            const pct = maxScore > 0 ? Math.round(tf.importance / maxScore * 100) : 0;
            const actPct = Math.round((tf.activation_ratio || 0) * 100);
            const actColor = actPct > 60 ? '#16a34a' : actPct < 40 ? '#dc2626' : '#d97706';
            const goldColors = ['#f59e0b','#94a3b8','#cd7f32'];
            const rankBadge = i < 3
                ? `<span style="background:${goldColors[i]};color:#fff;border-radius:4px;padding:1px 5px;font-size:10px;font-weight:700;">#${i+1}</span>`
                : `<span style="color:var(--text-muted);">${i+1}</span>`;
            const sigmaTag = tf.is_sigma
                ? `<span style="font-size:9px;background:rgba(139,92,246,0.15);color:#8b5cf6;border-radius:3px;padding:0 4px;margin-left:4px;">σ</span>`
                : '';
            const displayName = tf.name && tf.name !== tf.locus ? tf.name : tf.locus;
            const bar = `<div style="background:#e2e8f0;border-radius:3px;height:6px;width:80px;overflow:hidden;display:inline-block;">
                <div style="background:linear-gradient(90deg,#6366f1,#8b5cf6);height:100%;width:${pct}%;transition:width 0.4s;"></div></div>`;

            return `<tr style="border-bottom:1px solid var(--border-color);cursor:pointer;"
                onmouseover="this.style.background='rgba(99,102,241,0.04)'"
                onmouseout="this.style.background=''">
                <td style="padding:7px 8px;">${rankBadge}</td>
                <td style="padding:7px 8px;">
                    <a href="#" class="topo-gene-link" data-locus="${escapeHtml(tf.locus)}"
                       style="font-weight:600;color:var(--color-primary-accent);text-decoration:none;"
                    >${escapeHtml(displayName)}${sigmaTag}</a>
                    <div style="font-size:10px;color:var(--text-muted);">${escapeHtml(tf.locus)}</div>
                </td>
                <td style="padding:7px 8px;text-align:center;font-weight:700;">${tf.out_degree}</td>
                <td style="padding:7px 8px;text-align:center;font-family:monospace;font-size:11px;"
                    title="Betweenness: ${tf.betweenness}">${(tf.betweenness * 1000).toFixed(2)}</td>
                <td style="padding:7px 8px;text-align:center;font-family:monospace;font-size:11px;"
                    title="PageRank: ${tf.pagerank}">${(tf.pagerank * 1000).toFixed(2)}</td>
                <td style="padding:7px 8px;text-align:center;font-family:monospace;font-size:11px;"
                    title="Hub score: ${tf.hub_score}">${(tf.hub_score * 1000).toFixed(2)}</td>
                <td style="padding:7px 8px;text-align:center;font-weight:600;color:${actColor};">${actPct}%</td>
                <td style="padding:7px 8px;text-align:center;font-weight:700;color:var(--color-primary-accent);"
                    title="Composite importance score">${(tf.importance * 100).toFixed(1)}</td>
                <td style="padding:7px 8px;">${bar}</td>
            </tr>`;
        }).join('');

        // Bind gene links
        tbody.querySelectorAll('.topo-gene-link').forEach(a => {
            a.addEventListener('click', e => {
                e.preventDefault();
                setActiveWorkflowEntry('gene');
                setTimeout(() => queryGene(a.getAttribute('data-locus')), 100);
            });
        });

        const meta = data._meta || {};
        showToast(`Centrality loaded: ${meta.n_tfs || tfs.length} TFs over ${meta.n_edges || '?'} edges`, 'success');
    } catch (err) {
        console.error('Centrality load error:', err);
        const tbody2 = document.getElementById('topo-centrality-tbody');
        if (tbody2) tbody2.innerHTML = `<tr><td colspan="9" style="text-align:center;padding:20px;color:#dc2626;">Error: ${escapeHtml(err.message)}</td></tr>`;
    }
}

// Auto-load when centrality tab is first clicked
document.addEventListener('click', e => {
    const btn = e.target.closest('[data-topo-tab="centrality"]');
    if (btn && !_centralityLoaded) {
        setTimeout(loadPrecomputedCentrality, 150);
    }
});

window.loadPrecomputedCentrality = loadPrecomputedCentrality;
"""

with open("web/app.js", "r", encoding="utf-8") as f:
    content = f.read()

# Append after the last window.closeEngineeringSimulationModal line
MARKER = "window.closeEngineeringSimulationModal = closeEngineeringSimulationModal;\n"
if MARKER in content:
    content = content.replace(MARKER, MARKER + JS_CODE, 1)
    with open("web/app.js", "w", encoding="utf-8") as f:
        f.write(content)
    print("SUCCESS: centrality JS appended to app.js")
else:
    print("ERROR: marker not found")
