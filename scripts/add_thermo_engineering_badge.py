"""scripts/add_thermo_engineering_badge.py — adds thermo badge to Engineering Targets rows"""

# Load the pre-computed thermo validation report to find thermo-essential genes
# In the Engineering Targets table, we want to add a 🔥 badge to TFs whose
# target genes include thermodynamically locked reactions

OLD_BADGE_SECTION = """        // Check if the TF itself is essential
        const locusLower = candidate.tfId.toLowerCase();
        const tfIsEssential = essentialGenes[locusLower] || (cgToCgl[locusLower] && essentialGenes[cgToCgl[locusLower].toLowerCase()]);
        const essentialBadge = tfIsEssential ? ` <span style="background:#fee2e2; color:#dc2626; font-size:8px; padding:1px 4px; border-radius:3px; font-weight:600; display:inline-block; vertical-align:middle; margin-left:4px;" title="This TF is essential for growth. Downregulation/knockout is lethal."><i class="fa-solid fa-triangle-exclamation"></i> Essential</span>` : '';

        // Check if the TF has Abasy systemic role warning
        const abasyRoleInfo = abasyRoles[locusLower] || (cgToCgl[locusLower] && abasyRoles[cgToCgl[locusLower].toLowerCase()]);
        let abasyBadge = '';
        if (abasyRoleInfo) {
            const role = abasyRoleInfo.role;
            if (role === 'Global Regulator' || role === 'Basal Machinery') {
                abasyBadge = ` <span style="background:#fef3c7; color:#d97706; font-size:8px; padding:1px 4px; border-radius:3px; font-weight:600; display:inline-block; vertical-align:middle; margin-left:4px;" title="Abasy role: ${role}. Modification of global hubs carries high pleiotropic risk of metabolic failure."><i class="fa-solid fa-circle-nodes"></i> Global Hub</span>`;
            }
        }"""

NEW_BADGE_SECTION = """        // Check if the TF itself is essential
        const locusLower = candidate.tfId.toLowerCase();
        const tfIsEssential = essentialGenes[locusLower] || (cgToCgl[locusLower] && essentialGenes[cgToCgl[locusLower].toLowerCase()]);
        const essentialBadge = tfIsEssential ? ` <span style="background:#fee2e2; color:#dc2626; font-size:8px; padding:1px 4px; border-radius:3px; font-weight:600; display:inline-block; vertical-align:middle; margin-left:4px;" title="This TF is essential for growth. Downregulation/knockout is lethal."><i class="fa-solid fa-triangle-exclamation"></i> Essential</span>` : '';

        // Check if the TF has Abasy systemic role warning
        const abasyRoleInfo = abasyRoles[locusLower] || (cgToCgl[locusLower] && abasyRoles[cgToCgl[locusLower].toLowerCase()]);
        let abasyBadge = '';
        if (abasyRoleInfo) {
            const role = abasyRoleInfo.role;
            if (role === 'Global Regulator' || role === 'Basal Machinery') {
                abasyBadge = ` <span style="background:#fef3c7; color:#d97706; font-size:8px; padding:1px 4px; border-radius:3px; font-weight:600; display:inline-block; vertical-align:middle; margin-left:4px;" title="Abasy role: ${role}. Modification of global hubs carries high pleiotropic risk of metabolic failure."><i class="fa-solid fa-circle-nodes"></i> Global Hub</span>`;
            }
        }

        // Thermodynamic context badge (from centrality data)
        const thermoCtx = _thermoContextCache.get(locusLower) || _thermoContextCache.get(candidate.tfId);
        let thermoBadge = '';
        if (thermoCtx) {
            const lvl = thermoCtx.thermo_support_level;
            const nLocked = thermoCtx.n_locked || 0;
            if (lvl === 'strong') {
                thermoBadge = ` <span style="background:#dcfce7; color:#16a34a; font-size:8px; padding:1px 5px; border-radius:3px; font-weight:700; display:inline-block; vertical-align:middle; margin-left:4px;" title="Thermodynamically constrained: ${nLocked} of ${thermoCtx.total_reactions} reactions are direction-locked. Predictions are thermodynamically supported."><i class="fa-solid fa-fire"></i> Thermo-constrained</span>`;
            } else if (lvl === 'moderate') {
                thermoBadge = ` <span style="background:#fef9c3; color:#d97706; font-size:8px; padding:1px 5px; border-radius:3px; font-weight:600; display:inline-block; vertical-align:middle; margin-left:4px;" title="Partial thermodynamic coverage: ${nLocked} locked reactions"><i class="fa-solid fa-bolt"></i> Thermo-partial</span>`;
            }
        }"""

OLD_TF_DISPLAY = "                <td><strong>${tfDisplay}</strong>${essentialBadge}${abasyBadge}</td>"
NEW_TF_DISPLAY = "                <td><strong>${tfDisplay}</strong>${essentialBadge}${abasyBadge}${thermoBadge}</td>"

with open("web/app.js", "r", encoding="utf-8") as f:
    content = f.read()

changed = 0
if OLD_BADGE_SECTION in content:
    content = content.replace(OLD_BADGE_SECTION, NEW_BADGE_SECTION, 1)
    changed += 1
    print("OK: Thermo badge code injected into Engineering Targets")
else:
    print("WARN: badge section not found")

if OLD_TF_DISPLAY in content:
    content = content.replace(OLD_TF_DISPLAY, NEW_TF_DISPLAY, 1)
    changed += 1
    print("OK: thermoBadge added to TF display cell")
else:
    print("WARN: TF display cell not found")

if changed > 0:
    with open("web/app.js", "w", encoding="utf-8") as f:
        f.write(content)
    print(f"SUCCESS: {changed} changes applied")
