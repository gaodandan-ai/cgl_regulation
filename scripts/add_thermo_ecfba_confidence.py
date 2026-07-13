"""scripts/add_thermo_ecfba_confidence.py — adds thermo confidence bar to ecFBA simulation panel"""

OLD_AFTER_SUMMARY = """        }
        
        engineeringSimChart = new Chart(ctxSim, {"""

NEW_AFTER_SUMMARY = """        }

        // ── Thermodynamic Confidence Badge ──────────────────────────────────
        const thermoConf = document.getElementById('engineering-sim-thermo-confidence');
        if (thermoConf) {
            thermoConf.innerHTML = '<i class="fa-solid fa-spinner fa-spin"></i> Checking thermodynamic context…';
            fetchThermoContext(tfId).then(ctx => {
                if (!ctx) {
                    thermoConf.innerHTML = '<span style="color:var(--text-muted);font-size:10px;">No thermodynamic data available for this TF.</span>';
                    return;
                }
                const lvl = ctx.thermo_support_level;
                const conf = ctx.ko_thermo_confidence;
                const n = ctx.n_locked;
                const tot = ctx.total_reactions;
                const confPct = Math.round(conf * 100);
                const lvlConfig = {
                    'strong':   { color: '#16a34a', bg: 'rgba(22,163,74,0.08)',  label: '🔒 Strong' },
                    'moderate': { color: '#d97706', bg: 'rgba(217,119,6,0.08)',   label: '⚠️ Moderate' },
                    'weak':     { color: '#9ca3af', bg: 'rgba(156,163,175,0.06)', label: '〰️ Weak' },
                    'none':     { color: '#9ca3af', bg: 'rgba(156,163,175,0.04)', label: '❓ No data' },
                };
                const lc = lvlConfig[lvl] || lvlConfig['none'];
                thermoConf.innerHTML = `
                <div style="background:${lc.bg};border:1px solid ${lc.color}20;border-radius:8px;padding:8px 10px;margin-top:8px;">
                    <div style="display:flex;align-items:center;justify-content:space-between;margin-bottom:6px;">
                        <span style="font-size:11px;font-weight:700;color:${lc.color};">${lc.label} Thermodynamic Support</span>
                        <span style="font-size:10px;color:var(--text-secondary);">${n}/${tot} reactions direction-locked</span>
                    </div>
                    <div style="display:flex;align-items:center;gap:8px;">
                        <div style="flex:1;background:#e2e8f0;border-radius:4px;height:5px;overflow:hidden;">
                            <div style="background:${lc.color};height:100%;width:${confPct}%;transition:width 0.5s;"></div>
                        </div>
                        <span style="font-size:11px;font-weight:700;color:${lc.color};">${confPct}%</span>
                    </div>
                    <div style="font-size:9.5px;color:var(--text-muted);margin-top:4px;">
                        Prediction confidence based on thermodynamic direction constraints (Noor et al. 2013)
                    </div>
                </div>`;
            });
        }
        
        engineeringSimChart = new Chart(ctxSim, {"""

with open("web/app.js", "r", encoding="utf-8") as f:
    content = f.read()

# Count occurrences
count = content.count(OLD_AFTER_SUMMARY)
print(f"Occurrences of marker: {count}")

if count == 1:
    content = content.replace(OLD_AFTER_SUMMARY, NEW_AFTER_SUMMARY, 1)
    with open("web/app.js", "w", encoding="utf-8") as f:
        f.write(content)
    print("SUCCESS: Thermo confidence bar added to ecFBA simulation panel")
else:
    print("ERROR: marker not found or ambiguous")
