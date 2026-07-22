/**
 * srnaTargetCompetitionView.js
 * ============================
 * Renders IntaRNA free energy ranking bar chart and seed binding alignment boxes.
 */

window.SrnaTargetCompetitionView = {
    async render(containerId = "srnaCompetitionContainer", srnaId = "scgl257.1") {
        const container = document.getElementById(containerId);
        if (!container) return;

        container.innerHTML = `<div class="p-4 text-center text-slate-400">Loading sRNA Target Competition Rankings for <strong>${srnaId}</strong>...</div>`;

        try {
            const res = await fetch(`/api/ncrna/targets?locus=${encodeURIComponent(srnaId)}`).then(r => r.json());
            const targets = res.targets || [];

            if (targets.length === 0) {
                container.innerHTML = `<div class="p-4 text-slate-400">No sRNA target competition data for ${srnaId}</div>`;
                return;
            }

            let html = `
                <div class="bg-slate-900 border border-slate-700/80 rounded-xl p-5 shadow-2xl text-slate-100 font-sans space-y-4">
                    <div class="flex items-center justify-between border-b border-slate-700/60 pb-3">
                        <h3 class="text-lg font-bold text-purple-400 flex items-center gap-2">
                            <span>📡</span>
                            <span>sRNA Target Competition & IntaRNA Binding Free Energy (${srnaId})</span>
                        </h3>
                        <span class="text-xs text-purple-300 bg-purple-950 px-2 py-0.5 rounded border border-purple-800">
                            ${targets.length} Target mRNAs
                        </span>
                    </div>

                    <!-- Target Energy Ranking Bars -->
                    <div class="space-y-2.5">
            `;

            const minE = Math.min(...targets.map(t => t.binding_energy_kcal || 0));

            targets.slice(0, 10).forEach((t, i) => {
                const energy = t.binding_energy_kcal || 0;
                const percent = Math.min(100, Math.max(15, (Math.abs(energy) / Math.abs(minE || 1)) * 100));

                html += `
                    <div class="bg-slate-950 p-3 rounded-lg border border-slate-800 space-y-1.5 text-xs">
                        <div class="flex items-center justify-between">
                            <span class="font-bold text-slate-200">#${i + 1} Target: <span class="text-purple-300">${t.target_name} (${t.target_locus})</span></span>
                            <span class="font-mono text-purple-400 font-bold">${energy.toFixed(2)} kcal/mol</span>
                        </div>
                        <!-- Energy Bar -->
                        <div class="w-full bg-slate-900 h-2.5 rounded-full overflow-hidden border border-slate-800">
                            <div class="bg-gradient-to-r from-purple-500 to-indigo-500 h-full rounded-full transition-all" style="width: ${percent}%"></div>
                        </div>
                        <div class="flex items-center justify-between text-[11px] text-slate-400">
                            <span>Mechanism: ${t.regulatory_mechanism}</span>
                            <span>Region: ${t.target_region_type}</span>
                        </div>
                    </div>
                `;
            });

            html += `</div></div>`;
            container.innerHTML = html;
        } catch (err) {
            console.error("Failed to render sRNA Competition View:", err);
            container.innerHTML = `<div class="p-4 text-rose-500">Error rendering sRNA competition view: ${err.message}</div>`;
        }
    }
};
