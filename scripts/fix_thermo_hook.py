"""scripts/fix_thermo_hook.py — fixes the thermo context hook by overriding fetchMetabolicImpact properly"""

# The MutationObserver approach doesn't work well.
# Instead, we patch fetchMetabolicImpact directly to call the thermo context.
# Find the fetchMetabolicImpact function and add a chained call.

# The original fetchMetabolicImpact calls renderMetabolicImpact(data, locusTag, nodeType) inside .then()
# We patch the .then() to also call fetchThermoContext

OLD_LOAD_IMPACT = """    loadImpact
        .then(data => {
            if (detailLocusTag.textContent !== locusTag) return;
            renderMetabolicImpact(data, locusTag, nodeType);
        })
        .catch(err => {
            console.error('Error fetching metabolic impact:', err);
            if (detailLocusTag.textContent === locusTag) {
                container.innerHTML = '<div class=\"metabolic-empty\">Failed to load metabolic model mapping.</div>'; 
            }
        });
}"""

NEW_LOAD_IMPACT = """    loadImpact
        .then(data => {
            if (detailLocusTag.textContent !== locusTag) return;
            renderMetabolicImpact(data, locusTag, nodeType);
            // ── Thermodynamic Context Card ────────────────────────────────
            // Append thermo card after metabolic impact renders
            if (typeof fetchThermoContext === 'function') {
                const section = document.getElementById('detail-metabolic-impact-section');
                if (section) {
                    // Create or reuse the thermo container
                    let thermoDiv = document.getElementById('thermo-context-card');
                    if (!thermoDiv) {
                        thermoDiv = document.createElement('div');
                        thermoDiv.id = 'thermo-context-card';
                        thermoDiv.style.cssText = 'margin:0 0 8px;padding:0;';
                        section.appendChild(thermoDiv);
                    }
                    thermoDiv.innerHTML = '<div style="padding:6px 0;color:var(--text-muted);font-size:10.5px;"><i class="fa-solid fa-spinner fa-spin"></i> Loading thermodynamic context\u2026</div>';
                    fetchThermoContext(locusTag).then(ctx => {
                        if (ctx && typeof renderThermoContextCard === 'function') {
                            renderThermoContextCard(ctx, thermoDiv);
                        } else if (thermoDiv) {
                            thermoDiv.innerHTML = '';
                        }
                    });
                }
            }
        })
        .catch(err => {
            console.error('Error fetching metabolic impact:', err);
            if (detailLocusTag.textContent === locusTag) {
                container.innerHTML = '<div class=\"metabolic-empty\">Failed to load metabolic model mapping.</div>'; 
            }
        });
}"""

with open("web/app.js", "r", encoding="utf-8") as f:
    content = f.read()

count = content.count(OLD_LOAD_IMPACT)
print(f"Marker occurrences: {count}")

if count == 1:
    # Also remove the broken MutationObserver approach from the thermo UI block
    content = content.replace(OLD_LOAD_IMPACT, NEW_LOAD_IMPACT, 1)
    
    # Remove the MutationObserver block since we have direct integration now
    OLD_OBSERVER = """// Listen for gene detail load events using MutationObserver on the metabolic section
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
            thermoContainer.innerHTML = `<div style="padding:8px 0;color:var(--text-muted);font-size:11px;"><i class="fa-solid fa-spinner fa-spin"></i> Loading thermodynamic context\u2026</div>`;
            fetchThermoContext(locus).then(ctx => {
                if (ctx) renderThermoContextCard(ctx, thermoContainer);
                else thermoContainer.innerHTML = '';
            });
        });
        obs.observe(metabolicSection, { childList: true, subtree: false, attributes: true, attributeFilter: ['style'] });
    }"""
    if OLD_OBSERVER in content:
        content = content.replace(OLD_OBSERVER, "// MutationObserver removed — thermo context now integrated in fetchMetabolicImpact", 1)
        print("OK: MutationObserver removed")
    
    with open("web/app.js", "w", encoding="utf-8") as f:
        f.write(content)
    print("SUCCESS: fetchMetabolicImpact patched with direct thermo context call")
else:
    print("ERROR: marker not found")
    idx = content.find("loadImpact")
    print(f"  Found loadImpact at char {idx}")
