"""scripts/add_thermo_pathway_color.py — adds thermo coloring to Pathway View reaction nodes"""

# 1. Add thermo coloring after pathway cytoscape is built (line 1486)
MARKER_PATHWAY = "                if (loading) loading.classList.add('hidden');\n\n                // Update stats"

PATHWAY_THERMO_CODE = """                if (loading) loading.classList.add('hidden');

                // ── Thermodynamic Coloring of Reaction Nodes ──────────────────
                (async () => {
                    try {
                        const thermoResp = await fetch('/api/thermo/pruning-report');
                        if (!thermoResp.ok) return;
                        const thermoReport = await thermoResp.json();
                        const pruned = thermoReport.pruned_reactions || [];
                        const confirmed = thermoReport.confirmed_reactions || [];

                        // Build quick lookup maps
                        const fwdLocked = new Set(pruned.filter(r => r.direction === 'forward').map(r => r.reaction_id));
                        const revLocked = new Set(pruned.filter(r => r.direction === 'reverse').map(r => r.reaction_id));
                        const confirmedFwd = new Set(confirmed.filter(r => r.direction === 'forward').map(r => r.reaction_id));
                        const confirmedRev = new Set(confirmed.filter(r => r.direction === 'reverse').map(r => r.reaction_id));

                        // Apply colors to reaction nodes
                        if (pathwayKeggCy) {
                            pathwayKeggCy.nodes('[type="reaction"]').forEach(node => {
                                const rxnId = node.data('label');
                                if (fwdLocked.has(rxnId)) {
                                    node.style({ 'background-color': '#dcfce7', 'border-color': '#16a34a', 'border-width': '2.5px' });
                                    node.data('thermo', 'fwd-locked');
                                } else if (revLocked.has(rxnId)) {
                                    node.style({ 'background-color': '#fee2e2', 'border-color': '#dc2626', 'border-width': '2.5px' });
                                    node.data('thermo', 'rev-locked');
                                } else if (confirmedFwd.has(rxnId) || confirmedRev.has(rxnId)) {
                                    node.style({ 'background-color': '#fef9c3', 'border-color': '#d97706', 'border-width': '1.5px' });
                                    node.data('thermo', 'near-eq');
                                }
                            });
                        }
                    } catch (e) {
                        console.warn('Pathway thermo coloring failed:', e);
                    }
                })();

                // Update stats"""

# 2. Add thermo info to the reaction detail popup (line 1518-1528)
OLD_REACTION_DETAIL = """                            detailContent.innerHTML = `
                                <div style="margin-bottom:10px;">
                                    <div style="font-weight:700;color:var(--color-activation);font-size:12px;margin-bottom:4px;">${escapeHtml(data.label)}</div>
                                    <div style="color:var(--text-secondary);font-size:10px;">${escapeHtml(data.name || 'No name')}</div>
                                </div>
                                <div style="background:rgba(46,125,50,0.05);border:1px solid rgba(46,125,50,0.15);border-radius:6px;padding:8px;margin-bottom:8px;">
                                    <div style="font-size:9px;font-weight:700;color:var(--color-activation);margin-bottom:4px;text-transform:uppercase;">Equation</div>
                                    <div style="font-family:monospace;font-size:9.5px;color:#1b5e20;word-break:break-all;line-height:1.5;">${escapeHtml(data.equation || '—')}</div>
                                </div>
                                <div style="font-size:10px;margin-bottom:4px;"><span style="color:#4f46e5;font-weight:600;">↳ Substrates:</span> <span style="color:#312e81;">${escapeHtml(incoming)}</span></div>
                                <div style="font-size:10px;"><span style="color:#ea580c;font-weight:600;">↳ Products:</span> <span style="color:#7c2d12;">${escapeHtml(outgoing)}</span></div>`;"""

NEW_REACTION_DETAIL = """                            // Fetch thermo info for this reaction
                            const thermoInfo = await (async () => {
                                try {
                                    const r = await fetch('/api/thermo/pruning-report');
                                    if (!r.ok) return null;
                                    const rep = await r.json();
                                    const allRxns = [...(rep.pruned_reactions || []), ...(rep.confirmed_reactions || [])];
                                    return allRxns.find(rx => rx.reaction_id === data.label) || null;
                                } catch { return null; }
                            })();
                            const thermoTag = thermoInfo ? (() => {
                                const dir = thermoInfo.direction;
                                const dg = thermoInfo.dgr_prime_0 != null ? `ΔG'°=${thermoInfo.dgr_prime_0.toFixed(1)} kJ/mol` : '';
                                if (thermoInfo.status === 'newly_locked') {
                                    return dir === 'forward'
                                        ? `<div style="margin-top:8px;padding:6px 8px;background:#dcfce7;border:1px solid #16a34a;border-radius:6px;font-size:10px;">🔒 <strong>Forward-locked</strong> by thermodynamics${dg ? ' · ' + dg : ''}<br><span style="color:var(--text-muted)">ΔG'∈[${thermoInfo.dgr_prime_min?.toFixed(1)}, ${thermoInfo.dgr_prime_max?.toFixed(1)}] kJ/mol</span></div>`
                                        : `<div style="margin-top:8px;padding:6px 8px;background:#fee2e2;border:1px solid #dc2626;border-radius:6px;font-size:10px;">🔒 <strong>Reverse-locked</strong> by thermodynamics${dg ? ' · ' + dg : ''}</div>`;
                                }
                                return `<div style="margin-top:8px;padding:6px 8px;background:#fef9c3;border:1px solid #d97706;border-radius:6px;font-size:10px;">⚖️ Near-equilibrium · ${dg}</div>`;
                            })() : '';
                            detailContent.innerHTML = `
                                <div style="margin-bottom:10px;">
                                    <div style="font-weight:700;color:var(--color-activation);font-size:12px;margin-bottom:4px;">${escapeHtml(data.label)}</div>
                                    <div style="color:var(--text-secondary);font-size:10px;">${escapeHtml(data.name || 'No name')}</div>
                                </div>
                                <div style="background:rgba(46,125,50,0.05);border:1px solid rgba(46,125,50,0.15);border-radius:6px;padding:8px;margin-bottom:8px;">
                                    <div style="font-size:9px;font-weight:700;color:var(--color-activation);margin-bottom:4px;text-transform:uppercase;">Equation</div>
                                    <div style="font-family:monospace;font-size:9.5px;color:#1b5e20;word-break:break-all;line-height:1.5;">${escapeHtml(data.equation || '—')}</div>
                                </div>
                                <div style="font-size:10px;margin-bottom:4px;"><span style="color:#4f46e5;font-weight:600;">↳ Substrates:</span> <span style="color:#312e81;">${escapeHtml(incoming)}</span></div>
                                <div style="font-size:10px;"><span style="color:#ea580c;font-weight:600;">↳ Products:</span> <span style="color:#7c2d12;">${escapeHtml(outgoing)}</span></div>
                                ${thermoTag}`;"""

# Also need to make the tap handler async since we're using await inside it
OLD_TAP = "                pathwayKeggCy.on('tap', 'node', function(evt) {"
NEW_TAP = "                pathwayKeggCy.on('tap', 'node', async function(evt) {"

with open("web/app.js", "r", encoding="utf-8") as f:
    content = f.read()

changed = 0

if MARKER_PATHWAY in content:
    content = content.replace(MARKER_PATHWAY, PATHWAY_THERMO_CODE, 1)
    changed += 1
    print("OK: Pathway thermo coloring code injected")
else:
    print("WARN: pathway coloring marker not found")

if OLD_TAP in content:
    content = content.replace(OLD_TAP, NEW_TAP, 1)
    changed += 1
    print("OK: tap handler made async")
else:
    print("WARN: tap handler not found")

if OLD_REACTION_DETAIL in content:
    content = content.replace(OLD_REACTION_DETAIL, NEW_REACTION_DETAIL, 1)
    changed += 1
    print("OK: Reaction detail thermo info injected")
else:
    print("WARN: reaction detail marker not found; checking what we have...")
    idx = content.find("detailContent.innerHTML = `")
    print(f"  detailContent.innerHTML at char {idx}")

if changed > 0:
    with open("web/app.js", "w", encoding="utf-8") as f:
        f.write(content)
    print(f"SUCCESS: {changed} changes applied to app.js")
else:
    print("ERROR: No changes applied")
