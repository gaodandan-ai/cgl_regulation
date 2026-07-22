/**
 * icaConditionHeatmapView.js
 * ==========================
 * Renders 9-condition by 87-iModulon activity matrix heatmaps and F1-score overlap charts.
 */

window.IcaConditionHeatmapView = {
    async render(containerId = "icaHeatmapContainer", conditionName = "Glucose") {
        const container = document.getElementById(containerId);
        if (!container) return;

        container.innerHTML = `<div class="p-4 text-center text-slate-400">Loading ICA Condition Activity Matrix for <strong>${conditionName}</strong>...</div>`;

        try {
            const [actRes, overlapRes] = await Promise.all([
                fetch(`/api/imodulon/condition?condition=${encodeURIComponent(conditionName)}`).then(r => r.json()),
                fetch(`/api/imodulon/overlap`).then(r => r.json())
            ]);

            const activities = actRes.activities || [];
            const overlaps = overlapRes.overlaps || [];

            let html = `
                <div class="bg-slate-900 border border-slate-700/80 rounded-xl p-5 shadow-2xl text-slate-100 font-sans space-y-4">
                    <div class="flex items-center justify-between border-b border-slate-700/60 pb-3">
                        <h3 class="text-lg font-bold text-teal-400 flex items-center gap-2">
                            <span>🔥</span>
                            <span>ICA Condition Activity Heatmap (${activities.length} Active iModulons)</span>
                        </h3>
                        <span class="text-xs text-slate-400">Condition: ${conditionName}</span>
                    </div>

                    <!-- Heatmap Grid -->
                    <div class="grid grid-cols-1 md:grid-cols-3 gap-3 text-xs">
            `;

            activities.forEach(item => {
                const score = item.activity_score || 0;
                const isPos = score >= 0;
                const colorClass = isPos ? "bg-emerald-950/80 border-emerald-700/50 text-emerald-300" : "bg-sky-950/80 border-sky-700/50 text-sky-300";

                html += `
                    <div class="${colorClass} p-3 rounded-lg border flex items-center justify-between shadow-sm">
                        <div>
                            <span class="font-bold block">${item.imodulon_name}</span>
                            <span class="text-[11px] opacity-80">${item.linked_regulator} (${item.category})</span>
                        </div>
                        <span class="text-sm font-mono font-bold">${score > 0 ? '+' : ''}${score.toFixed(2)}</span>
                    </div>
                `;
            });

            html += `
                    </div>

                    <!-- iModulon-Regulon Overlap Scores -->
                    <div class="border-t border-slate-800 pt-3 space-y-2">
                        <h4 class="text-xs font-semibold text-slate-300">iModulon-Regulon F1-Score Overlap Alignments (${overlaps.length} Overlaps)</h4>
                        <div class="grid grid-cols-2 sm:grid-cols-4 gap-2 text-xs">
            `;

            overlaps.slice(0, 8).forEach(ov => {
                html += `
                    <div class="bg-slate-950 p-2 rounded border border-slate-800 flex items-center justify-between">
                        <span class="text-slate-300">${ov.imodulon_id} ➔ ${ov.tf_name}</span>
                        <span class="font-mono text-teal-400 font-bold">${(ov.f1_score * 100).toFixed(0)}%</span>
                    </div>
                `;
            });

            html += `</div></div></div>`;
            container.innerHTML = html;
        } catch (err) {
            console.error("Failed to render ICA Heatmap View:", err);
            container.innerHTML = `<div class="p-4 text-rose-500">Error rendering ICA heatmap view: ${err.message}</div>`;
        }
    }
};
