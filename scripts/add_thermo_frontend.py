"""scripts/add_thermo_frontend.py — injects thermodynamic context UI into key frontend areas"""

THERMO_JS = r'''

// ── Thermodynamic Context UI ──────────────────────────────────────────────────
// Cache for gene thermo context (avoid repeated fetches)
const _thermoContextCache = new Map();

async function fetchThermoContext(locus) {
    if (!locus) return null;
    if (_thermoContextCache.has(locus)) return _thermoContextCache.get(locus);
    try {
        const r = await fetch(`/api/thermo/gene_context?gene=${encodeURIComponent(locus)}`);
        if (!r.ok) return null;
        const data = await r.json();
        _thermoContextCache.set(locus, data);
        return data;
    } catch (e) {
        console.warn('Thermo context fetch failed:', e);
        return null;
    }
}

function renderThermoContextCard(ctx, containerEl) {
    if (!ctx || !containerEl) return;

    const level = ctx.thermo_support_level || 'none';
    const n_locked = ctx.n_locked || 0;
    const n_total = ctx.total_reactions || 0;
    const confidence = ctx.ko_thermo_confidence || 0;
    const annotated = ctx.thermo_annotated || [];

    // Level styling
    const levelConfig = {
        'strong':   { color: '#16a34a', bg: 'rgba(22,163,74,0.08)',  icon: '🔒', label: 'Strong' },
        'moderate': { color: '#d97706', bg: 'rgba(217,119,6,0.08)',   icon: '⚠️', label: 'Moderate' },
        'weak':     { color: '#9ca3af', bg: 'rgba(156,163,175,0.08)', icon: '〰️', label: 'Weak' },
        'none':     { color: '#9ca3af', bg: 'rgba(156,163,175,0.06)', icon: '❓', label: 'No data' },
    };
    const lc = levelConfig[level] || levelConfig['none'];

    // Build reaction rows — show locked ones first, max 8
    const rowsHtml = annotated.slice(0, 8).map(r => {
        const dir = r.direction_locked;
        let badge = '';
        if (dir === 'forward') {
            badge = `<span style="background:#dcfce7;color:#16a34a;border-radius:4px;padding:1px 5px;font-size:9px;font-weight:700;">→ FWD LOCK</span>`;
        } else if (dir === 'reverse') {
            badge = `<span style="background:#fee2e2;color:#dc2626;border-radius:4px;padding:1px 5px;font-size:9px;font-weight:700;">← REV LOCK</span>`;
        } else if (r.has_thermo_data) {
            badge = `<span style="background:#fef9c3;color:#92400e;border-radius:4px;padding:1px 5px;font-size:9px;">⇌ near-eq</span>`;
        } else {
            badge = `<span style="background:#f3f4f6;color:#9ca3af;border-radius:4px;padding:1px 5px;font-size:9px;">no data</span>`;
        }
        const dgr = r.dgr_prime_0 != null
            ? `<span style="font-family:monospace;font-size:10px;color:var(--text-secondary);">ΔG'°=${r.dgr_prime_0.toFixed(1)}</span>`
            : '';
        return `<div style="display:flex;align-items:center;gap:6px;padding:3px 0;border-bottom:1px solid var(--border-color);">
            <span style="font-weight:600;font-size:11px;min-width:80px;">${escapeHtml(r.reaction_id)}</span>
            ${badge}
            ${dgr}
        </div>`;
    }).join('');

    // Confidence bar
    const confPct = Math.round(confidence * 100);
    const confColor = confidence > 0.6 ? '#16a34a' : confidence > 0.3 ? '#d97706' : '#9ca3af';

    containerEl.innerHTML = `
    <div style="margin-top:10px;border:1px solid var(--border-color);border-radius:10px;overflow:hidden;">
        <div style="display:flex;align-items:center;justify-content:space-between;padding:8px 12px;background:${lc.bg};border-bottom:1px solid var(--border-color);">
            <div style="display:flex;align-items:center;gap:8px;">
                <span style="font-size:13px;">${lc.icon}</span>
                <span style="font-weight:600;font-size:12px;">Thermodynamic Support</span>
                <span style="font-size:11px;font-weight:700;color:${lc.color};">${lc.label}</span>
            </div>
            <div style="font-size:11px;color:var(--text-secondary);">${n_locked}/${n_total} reactions locked</div>
        </div>
        <div style="padding:10px 12px;">
            <div style="margin-bottom:8px;">
                <div style="display:flex;justify-content:space-between;font-size:11px;color:var(--text-secondary);margin-bottom:3px;">
                    <span>Prediction confidence</span>
                    <span style="color:${confColor};font-weight:600;">${confPct}%</span>
                </div>
                <div style="background:#e2e8f0;border-radius:4px;height:5px;overflow:hidden;">
                    <div style="background:${confColor};height:100%;width:${confPct}%;transition:width 0.5s;"></div>
                </div>
            </div>
            ${rowsHtml || '<div style="font-size:11px;color:var(--text-muted);">No reactions in thermodynamic database.</div>'}
            ${annotated.length > 8 ? `<div style="font-size:10px;color:var(--text-muted);margin-top:4px;">+${annotated.length - 8} more reactions…</div>` : ''}
        </div>
    </div>`;
}

// Inject into the existing fetchMetabolicImpact flow
const _origFetchMetabolicImpact = typeof fetchMetabolicImpact === 'function' ? fetchMetabolicImpact : null;

async function fetchMetabolicImpactWithThermo(locusTag, nodeType) {
    if (_origFetchMetabolicImpact) _origFetchMetabolicImpact(locusTag, nodeType);

    // Fetch thermo context in parallel
    const ctx = await fetchThermoContext(locusTag);
    if (!ctx) return;

    // Find the container after the metabolic impact section renders
    setTimeout(() => {
        const container = document.getElementById('metabolic-impact-content');
        if (!container) return;
        // Add thermo card at the bottom of the metabolic impact panel
        let thermoContainer = document.getElementById('thermo-context-card');
        if (!thermoContainer) {
            thermoContainer = document.createElement('div');
            thermoContainer.id = 'thermo-context-card';
            container.parentElement.appendChild(thermoContainer);
        }
        renderThermoContextCard(ctx, thermoContainer);
    }, 500);
}

window.fetchThermoContext = fetchThermoContext;
window.renderThermoContextCard = renderThermoContextCard;

// Override fetchMetabolicImpact in the side panel fetching code
// by hooking into the detail panel update sequence
const _origDetailPanel = window._detailPanelThermoHooked;
if (!_origDetailPanel) {
    window._detailPanelThermoHooked = true;

    // Listen for gene detail load events using MutationObserver on the metabolic section
    const metabolicSection = document.getElementById('detail-metabolic-impact-section');
    if (metabolicSection) {
        const obs = new MutationObserver(() => {
            const locus = (document.getElementById('detail-locus-tag') || {}).textContent || '';
            if (!locus) return;
            let thermoContainer = document.getElementById('thermo-context-card');
            if (!thermoContainer) {
                thermoContainer = document.createElement('div');
                thermoContainer.id = 'thermo-context-card';
                thermoContainer.style.padding = '0 0 8px 0';
                metabolicSection.appendChild(thermoContainer);
            }
            thermoContainer.innerHTML = `<div style="padding:8px 0;color:var(--text-muted);font-size:11px;"><i class="fa-solid fa-spinner fa-spin"></i> Loading thermodynamic context…</div>`;
            fetchThermoContext(locus).then(ctx => {
                if (ctx) renderThermoContextCard(ctx, thermoContainer);
                else thermoContainer.innerHTML = '';
            });
        });
        obs.observe(metabolicSection, { childList: true, subtree: false, attributes: true, attributeFilter: ['style'] });
    }
}

'''

# Inject the thermo UI JS at the end of app.js
MARKER = "window.loadPrecomputedCentrality = loadPrecomputedCentrality;\n"

with open("web/app.js", "r", encoding="utf-8") as f:
    content = f.read()

if MARKER in content:
    content = content.replace(MARKER, MARKER + THERMO_JS, 1)
    with open("web/app.js", "w", encoding="utf-8") as f:
        f.write(content)
    print("SUCCESS: Thermodynamic context UI injected into app.js")
else:
    print("ERROR: marker not found")
    idx = content.find("loadPrecomputedCentrality")
    print(f"  Found at: {idx}")
