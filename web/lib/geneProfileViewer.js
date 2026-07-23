/**
 * geneProfileViewer.js
 * ====================
 * Renders HTML5 Canvas linear genomic tracks (+/- 20kb genomic window)
 * and 360-degree gene profile cards.
 */

window.GeneProfileViewer = {
    async loadGeneProfile(geneId, containerId = "geneProfileContainer") {
        const container = document.getElementById(containerId);
        if (!container) return;

        container.innerHTML = `<div class="p-4 text-center text-slate-400">Loading 360° Profile for <strong>${geneId}</strong>...</div>`;

        try {
            const [profileRes, neighborhoodRes] = await Promise.all([
                fetch(`/api/gene/profile/${encodeURIComponent(geneId)}`).then(r => r.ok ? r.json() : null),
                fetch(`/api/gene/neighborhood/${encodeURIComponent(geneId)}?window_bp=20000`).then(r => r.ok ? r.json() : null)
            ]);

            if (!profileRes) {
                container.innerHTML = `<div class="p-4 text-rose-500">Gene profile not found for ${geneId}</div>`;
                return;
            }

            this.renderProfileCard(container, profileRes, neighborhoodRes);
        } catch (err) {
            console.error("Failed to load gene profile:", err);
            container.innerHTML = `<div class="p-4 text-rose-500">Error loading gene profile: ${err.message}</div>`;
        }
    },

    renderProfileCard(container, profile, neighborhood) {
        const genes = neighborhood ? neighborhood.genes || [] : [];

        container.innerHTML = `
            <div class="bg-slate-900 border border-slate-700/80 rounded-xl p-5 shadow-2xl text-slate-100 font-sans space-y-4">
                <!-- Header -->
                <div class="flex items-center justify-between border-b border-slate-700/60 pb-3">
                    <div>
                        <h3 class="text-xl font-bold text-teal-400 flex items-center gap-2">
                            <span>🧬</span>
                            <span>${profile.gene_name || profile.cg_locus} (${profile.cg_locus.toUpperCase()})</span>
                        </h3>
                        <p class="text-xs text-slate-400 mt-1">${profile.product || 'Hypothetical protein'}</p>
                    </div>
                    <div class="text-right">
                        <span class="inline-block bg-teal-950/80 text-teal-300 text-xs px-2.5 py-1 rounded-full border border-teal-700/50">
                            ${profile.abasy_role || 'Metabolic Gene'}
                        </span>
                    </div>
                </div>

                <!-- Attributes Grid -->
                <div class="grid grid-cols-2 sm:grid-cols-4 gap-3 text-xs">
                    <div class="bg-slate-800/60 p-2.5 rounded-lg border border-slate-700/40">
                        <span class="text-slate-400 block">RefSeq Locus</span>
                        <span class="font-semibold text-slate-200">${profile.cgl_locus || 'N/A'}</span>
                    </div>
                    <div class="bg-slate-800/60 p-2.5 rounded-lg border border-slate-700/40">
                        <span class="text-slate-400 block">Coordinates</span>
                        <span class="font-semibold text-slate-200">${profile.start_pos || 'N/A'} .. ${profile.end_pos || 'N/A'} (${profile.strand || '+'})</span>
                    </div>
                    <div class="bg-slate-800/60 p-2.5 rounded-lg border border-slate-700/40">
                        <span class="text-slate-400 block">TF Family</span>
                        <span class="font-semibold text-emerald-400">${profile.tf_family || 'N/A'}</span>
                    </div>
                    <div class="bg-slate-800/60 p-2.5 rounded-lg border border-slate-700/40">
                        <span class="text-slate-400 block">Effector Molecule</span>
                        <span class="font-semibold text-amber-400">${profile.effector_molecule || 'None'}</span>
                    </div>
                </div>

                <!-- 5-Track Interactive Genomic Track Browser -->
                <div id="genomicTrack5Container" class="mt-3"></div>
            </div>
        `;

        if (window.GenomicTrackBrowser) {
            window.GenomicTrackBrowser.render('genomicTrack5Container', profile.cg_locus || profile.cgl_locus);
        } else {
            this.drawTrackCanvas(profile, genes);
        }
    },

    drawTrackCanvas(centerGene, genes) {
        const canvas = document.getElementById("genomicTrackCanvas");
        if (!canvas) return;
        const ctx = canvas.getContext("2d");
        const W = canvas.width;
        const H = canvas.height;

        ctx.clearRect(0, 0, W, H);

        if (!genes || genes.length === 0) {
            ctx.fillStyle = "#64748b";
            ctx.font = "12px sans-serif";
            ctx.fillText("No genomic coordinate data available", W / 2 - 100, H / 2);
            return;
        }

        const minPos = Math.min(...genes.map(g => g.start_pos));
        const maxPos = Math.max(...genes.map(g => g.end_pos));
        const range = maxPos - minPos || 1;

        // Baseline
        ctx.strokeStyle = "#334155";
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(10, H / 2);
        ctx.lineTo(W - 10, H / 2);
        ctx.stroke();

        // Draw genes
        genes.forEach(g => {
            const x1 = 10 + ((g.start_pos - minPos) / range) * (W - 20);
            const x2 = 10 + ((g.end_pos - minPos) / range) * (W - 20);
            const w = Math.max(x2 - x1, 8);

            const isCenter = g.locus_tag.toLowerCase() === centerGene.cg_locus.toLowerCase();
            const y = g.strand === "+" ? H / 2 - 18 : H / 2 + 4;
            const h = 14;

            ctx.fillStyle = isCenter ? "#14b8a6" : (g.strand === "+" ? "#3b82f6" : "#ef4444");
            ctx.fillRect(x1, y, w, h);

            if (isCenter) {
                ctx.strokeStyle = "#f59e0b";
                ctx.lineWidth = 2;
                ctx.strokeRect(x1 - 2, y - 2, w + 4, h + 4);
            }

            // Label
            ctx.fillStyle = isCenter ? "#fbbf24" : "#94a3b8";
            ctx.font = "10px sans-serif";
            ctx.fillText(g.gene_name || g.locus_tag, x1, g.strand === "+" ? y - 4 : y + h + 12);
        });
    }
};
