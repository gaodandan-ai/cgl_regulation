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

                <!-- ChIP-seq Peak & Promoter Inspector Card -->
                <div id="chipseqPeakInspectorCard" class="mt-4 bg-slate-950/90 border border-sky-800/60 rounded-xl p-4 space-y-3">
                    <div class="flex items-center justify-between border-b border-slate-800 pb-2">
                        <div class="flex items-center gap-2">
                            <span class="text-sky-400">🎯</span>
                            <h4 class="text-sm font-bold text-slate-100">ChIP-seq 结合峰与启动子调控全景 (Peak & Promoter Inspector)</h4>
                        </div>
                        <span class="text-[11px] bg-sky-950 text-sky-300 border border-sky-700 px-2 py-0.5 rounded font-mono" id="chipseq-peak-badge-count">Querying peaks...</span>
                    </div>
                    <div id="chipseq-peak-table-wrap" class="overflow-x-auto text-xs">
                        <div class="text-center py-4 text-slate-400 font-mono"><i class="fa-solid fa-spinner fa-spin"></i> Loading experimental ChIP-seq binding summits...</div>
                    </div>
                </div>
            </div>
        `;

        const locus = profile.cg_locus || profile.cgl_locus;
        if (window.GenomicTrackBrowser) {
            window.GenomicTrackBrowser.render('genomicTrack5Container', locus);
        } else {
            this.drawTrackCanvas(profile, genes);
        }

        this.loadChipSeqPeaksForProfile(locus);
    },

    async loadChipSeqPeaksForProfile(locus) {
        const wrap = document.getElementById("chipseq-peak-table-wrap");
        const countBadge = document.getElementById("chipseq-peak-badge-count");
        if (!wrap) return;

        try {
            const resp = await fetch(`/api/chipseq_peaks/${encodeURIComponent(locus)}`);
            if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
            const data = await resp.json();

            if (data.is_public_deployment) {
                if (countBadge) {
                    countBadge.textContent = "🔒 实验室内网保护";
                    countBadge.className = "text-[11px] bg-amber-950 text-amber-300 border border-amber-700 px-2 py-0.5 rounded font-mono";
                }
                wrap.innerHTML = `
                    <div class="p-3 bg-amber-950/40 border border-amber-800/60 rounded-lg text-amber-200 text-xs font-sans space-y-1">
                        <div class="font-bold flex items-center gap-1.5 text-amber-300">
                            <span>🔒</span> <span>实验室内网未发布数据保护 (Lab Intranet Server 172.16.2.105)</span>
                        </div>
                        <p class="text-[11px] text-amber-200/80">13,673 条实测高分辨率 ChIP-seq Binding Summits 与空间定位轨迹仅在课题组内网服务器 (172.16.2.105:8010) 上提供展示。</p>
                    </div>
                `;
                return;
            }

            const peaks = (data.as_target_peaks || []).concat(data.as_tf_peaks || []);
            if (countBadge) {
                countBadge.textContent = `${peaks.length} Peaks Mapped (🧪 内网版 172.16.2.105)`;
            }

            if (peaks.length === 0) {
                wrap.innerHTML = `<div class="text-center py-4 text-slate-500 font-mono">No direct experimental ChIP-seq peak summits mapped for ${locus}</div>`;
                return;
            }

            let html = `
                <table class="w-full text-left border-collapse font-sans">
                    <thead>
                        <tr class="border-b border-slate-800 text-slate-400 font-semibold text-[11px]">
                            <th class="py-1.5 px-2">TF Name</th>
                            <th class="py-1.5 px-2">Peak Summit (bp)</th>
                            <th class="py-1.5 px-2">Signal Enrichment</th>
                            <th class="py-1.5 px-2">-log10(q)</th>
                            <th class="py-1.5 px-2">TSS Offset</th>
                            <th class="py-1.5 px-2">Spatial Category</th>
                        </tr>
                    </thead>
                    <tbody class="divide-y divide-slate-800/60 font-mono text-[11px]">
            `;

            peaks.slice(0, 10).forEach(p => {
                const tf = p.tf_name || p.tf_id || 'TF';
                const center = p.peak_center || p.peak_start || 'N/A';
                const score = (p.peak_score || p.peak_signal || 1.0).toFixed(2);
                const negq = (p.neglog10q || 0.0).toFixed(1);
                const relTss = p.rel_pos_to_tss != null ? (p.rel_pos_to_tss >= 0 ? `+${p.rel_pos_to_tss}` : `${p.rel_pos_to_tss}`) + ' bp' : 'distal';
                const spatial = p.spatial_confidence || 'PROMOTER_DIRECT';

                const badgeClass = spatial === 'PROMOTER_DIRECT' ? 'bg-emerald-950 text-emerald-300 border-emerald-700'
                    : (spatial === 'INTERGENIC_PROMOTER' ? 'bg-teal-950 text-teal-300 border-teal-700' : 'bg-slate-800 text-slate-300 border-slate-700');

                html += `
                    <tr class="hover:bg-slate-800/40 transition-colors">
                        <td class="py-1.5 px-2 text-sky-400 font-bold font-sans">${tf}</td>
                        <td class="py-1.5 px-2 text-amber-300">${typeof center === 'number' ? center.toLocaleString() : center}</td>
                        <td class="py-1.5 px-2 text-slate-200 font-bold">${score}x</td>
                        <td class="py-1.5 px-2 text-indigo-300">${negq}</td>
                        <td class="py-1.5 px-2 text-emerald-400">${relTss}</td>
                        <td class="py-1.5 px-2"><span class="px-1.5 py-0.5 rounded border ${badgeClass} text-[10px]">${spatial}</span></td>
                    </tr>
                `;
            });

            html += `</tbody></table>`;
            wrap.innerHTML = html;
        } catch (err) {
            wrap.innerHTML = `<div class="text-center py-3 text-slate-500">Failed to load peak details: ${err.message}</div>`;
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
