/**
 * allostericFeedbackView.js
 * =========================
 * Renders circular metabolite-TF-enzyme allosteric feedback loop diagrams.
 */

window.AllostericFeedbackView = {
    async render(containerId = "allostericFeedbackContainer", filterQuery = "cAMP") {
        const container = document.getElementById(containerId);
        if (!container) return;

        container.innerHTML = `<div class="p-4 text-center text-slate-400">Loading Allosteric Feedback Loops for <strong>${filterQuery}</strong>...</div>`;

        try {
            const res = await fetch(`/api/network/allosteric-feedback?query=${encodeURIComponent(filterQuery)}`).then(r => r.json());
            const loops = res.loops || [];

            if (loops.length === 0) {
                container.innerHTML = `<div class="p-4 text-slate-400">No allosteric feedback loops found for query: ${filterQuery}</div>`;
                return;
            }

            let html = `
                <div class="bg-slate-900 border border-slate-700/80 rounded-xl p-5 shadow-2xl text-slate-100 font-sans space-y-4">
                    <div class="flex items-center justify-between border-b border-slate-700/60 pb-3">
                        <h3 class="text-lg font-bold text-amber-400 flex items-center gap-2">
                            <span>🔄</span>
                            <span>Metabolite-TF-Enzyme Allosteric Feedback Loops (${loops.length} Loops)</span>
                        </h3>
                        <span class="text-xs text-slate-400">Signal-Driven Autoregulation</span>
                    </div>

                    <div class="grid grid-cols-1 md:grid-cols-2 gap-4">
            `;

            loops.forEach((item, idx) => {
                html += `
                    <div class="bg-slate-950 p-4 rounded-xl border border-slate-800 hover:border-amber-500/50 transition-all shadow-md space-y-3">
                        <div class="flex items-center justify-between">
                            <span class="text-xs font-semibold text-amber-400 bg-amber-950/80 px-2 py-0.5 rounded border border-amber-700/40">
                                💊 Effector: ${item.effector_molecule || 'Small Molecule'}
                            </span>
                            <span class="text-xs text-slate-400">${item.physiological_signal || 'Metabolic Sensor'}</span>
                        </div>

                        <!-- Closed Loop Node Flow -->
                        <div class="flex items-center justify-around py-3 bg-slate-900/80 rounded-lg border border-slate-800 text-xs">
                            <div class="text-center">
                                <span class="block text-slate-400">TF Node</span>
                                <span class="font-bold text-teal-300 text-sm">${item.tf_name || item.tf_locus}</span>
                            </div>
                            <span class="text-amber-400 text-lg">➔</span>
                            <div class="text-center">
                                <span class="block text-slate-400">Regulation</span>
                                <span class="font-semibold text-emerald-400">Target</span>
                            </div>
                            <span class="text-amber-400 text-lg">➔</span>
                            <div class="text-center">
                                <span class="block text-slate-400">Target Enzyme</span>
                                <span class="font-bold text-indigo-300 text-sm">${item.target_name || item.target_locus}</span>
                            </div>
                        </div>
                    </div>
                `;
            });

            html += `</div></div>`;
            container.innerHTML = html;
        } catch (err) {
            console.error("Failed to render Allosteric Feedback Loops:", err);
            container.innerHTML = `<div class="p-4 text-rose-500">Error rendering feedback loops: ${err.message}</div>`;
        }
    }
};
