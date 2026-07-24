/**
 * genomicTrackBrowser.js
 * =======================
 * Interactive 5-Track Genomic Track Inspector for C. glutamicum.
 * Visualizes:
 *  Track 1: Coordinate Ruler (bp) & Strand Direction
 *  Track 2: CDS Gene Models (+ / - strand arrows)
 *  Track 3: TSS & Promoter 70bp Region Highlights
 *  Track 4: ChIP-seq / TFBS Binding Peak Signal Density Curves
 *  Track 5: sRNA & ncRNA Annotations
 */

window.GenomicTrackBrowser = {
    state: {
        locusTag: null,
        data: null,
        windowBp: 10000,
        centerPos: 0,
        minPos: 0,
        maxPos: 0,
    },

    async render(containerId, locusTag) {
        const container = typeof containerId === 'string' ? document.getElementById(containerId) : containerId;
        if (!container) return;

        this.state.locusTag = locusTag;
        container.innerHTML = `
            <div class="bg-slate-900 border border-slate-700/80 rounded-xl p-4 shadow-2xl text-slate-100 font-sans space-y-3">
                <div class="flex items-center justify-between border-b border-slate-800 pb-2">
                    <div class="flex items-center gap-2">
                        <span class="text-teal-400 text-base">🧬</span>
                        <h4 class="font-bold text-sm text-slate-100">Genomic Track Browser</h4>
                        <span class="text-xs bg-slate-800 text-slate-300 px-2 py-0.5 rounded font-mono">${locusTag}</span>
                    </div>
                    <div class="flex items-center gap-1">
                        <button onclick="GenomicTrackBrowser.zoom(0.6)" class="px-2 py-1 text-xs bg-slate-800 hover:bg-slate-700 text-slate-200 rounded border border-slate-700">🔍 In</button>
                        <button onclick="GenomicTrackBrowser.zoom(1.5)" class="px-2 py-1 text-xs bg-slate-800 hover:bg-slate-700 text-slate-200 rounded border border-slate-700">🔍 Out</button>
                        <button onclick="GenomicTrackBrowser.shift(-3000)" class="px-2 py-1 text-xs bg-slate-800 hover:bg-slate-700 text-slate-200 rounded border border-slate-700">◀ Left</button>
                        <button onclick="GenomicTrackBrowser.shift(3000)" class="px-2 py-1 text-xs bg-slate-800 hover:bg-slate-700 text-slate-200 rounded border border-slate-700">Right ▶</button>
                    </div>
                </div>
                <div id="gtb-canvas-wrap" class="relative overflow-hidden bg-slate-950 rounded-lg p-2 border border-slate-800/80">
                    <div class="text-center text-xs text-slate-400 py-8"><i class="fa-solid fa-spinner fa-spin"></i> Loading genomic tracks for ${locusTag}...</div>
                </div>
                <div id="gtb-inspector" class="text-xs bg-slate-950 p-2.5 rounded border border-slate-800 text-slate-300 font-mono hidden">
                    Hover over a gene, promoter or peak to inspect sequence & coordinates.
                </div>
                <!-- AI Omics Copilot Quick Command Pills Bar -->
                <div class="flex items-center justify-between pt-2 border-t border-slate-800/80">
                    <span class="text-[11px] text-teal-400 font-semibold flex items-center gap-1">
                        <span>🤖</span>
                        <span>AI Omics Copilot:</span>
                    </span>
                    <div class="flex items-center gap-2">
                        <button onclick="window.currentSelectedLocus='${locusTag}'; window.sendEngineeringAiCommand('crispri')"
                                class="px-2 py-1 text-xs bg-indigo-950/80 hover:bg-indigo-900 text-indigo-200 border border-indigo-700/60 rounded flex items-center gap-1 transition-all cursor-pointer">
                            <span>🧪</span> CRISPRi Guide
                        </button>
                        <button onclick="window.currentSelectedLocus='${locusTag}'; window.sendEngineeringAiCommand('bottleneck')"
                                class="px-2 py-1 text-xs bg-emerald-950/80 hover:bg-emerald-900 text-emerald-200 border border-emerald-700/60 rounded flex items-center gap-1 transition-all cursor-pointer">
                            <span>⚡</span> FBA Bottleneck
                        </button>
                        <button onclick="window.currentSelectedLocus='${locusTag}'; window.sendEngineeringAiCommand('promoter')"
                                class="px-2 py-1 text-xs bg-purple-950/80 hover:bg-purple-900 text-purple-200 border border-purple-700/60 rounded flex items-center gap-1 transition-all cursor-pointer">
                            <span>🧬</span> Promoter Mutagenesis
                        </button>
                    </div>
                </div>
            </div>
        `;

        await this.loadData(locusTag);
    },

    async loadData(locusTag) {
        const wrap = document.getElementById("gtb-canvas-wrap");
        if (!wrap) return;

        try {
            const resp = await fetch(`/api/genomic_tracks/${encodeURIComponent(locusTag)}?window_bp=${this.state.windowBp}`);
            if (!resp.ok) throw new Error(`Track API error ${resp.status}`);
            const data = await resp.json();
            this.state.data = data;
            this.state.minPos = data.window.min_pos;
            this.state.maxPos = data.window.max_pos;
            this.state.centerPos = data.window.center_pos;
            this.drawTracks();
        } catch (err) {
            wrap.innerHTML = `<div class="text-rose-400 text-xs py-4 text-center">Failed to load genomic tracks: ${err.message}</div>`;
        }
    },

    zoom(factor) {
        this.state.windowBp = Math.max(2000, Math.min(50000, Math.round(this.state.windowBp * factor)));
        if (this.state.locusTag) this.loadData(this.state.locusTag);
    },

    shift(offsetBp) {
        if (!this.state.data) return;
        this.state.minPos += offsetBp;
        this.state.maxPos += offsetBp;
        this.drawTracks();
    },

    drawTracks() {
        const wrap = document.getElementById("gtb-canvas-wrap");
        if (!wrap || !this.state.data) return;

        const data = this.state.data;
        const minPos = this.state.minPos;
        const maxPos = this.state.maxPos;
        const totalBp = maxPos - minPos;
        const width = wrap.clientWidth || 800;
        const posToX = (pos) => Math.max(0, Math.min(width, ((pos - minPos) / totalBp) * width));

        let svg = `
            <svg width="100%" height="280" viewBox="0 0 ${width} 280" xmlns="http://www.w3.org/2000/svg" class="select-none font-sans">
                <!-- Track Background Grid -->
                <rect x="0" y="0" width="${width}" height="280" fill="#020617" />

                <!-- TRACK 1: Ruler & Coordinate Scale -->
                <g id="track-ruler">
                    <line x1="0" y1="30" x2="${width}" y2="30" stroke="#334155" stroke-width="1.5" />
        `;

        // Ticks every ~100px
        const tickStepBp = Math.pow(10, Math.floor(Math.log10(totalBp / 5)));
        const firstTick = Math.ceil(minPos / tickStepBp) * tickStepBp;

        for (let pos = firstTick; pos <= maxPos; pos += tickStepBp) {
            const x = posToX(pos);
            svg += `
                <line x1="${x}" y1="24" x2="${x}" y2="30" stroke="#64748b" stroke-width="1" />
                <text x="${x}" y="18" fill="#94a3b8" font-size="10" text-anchor="middle" font-family="monospace">${pos.toLocaleString()} bp</text>
            `;
        }

        svg += `
                </g>

                <!-- TRACK 2: CDS Gene Models (+/- strands) -->
                <g id="track-cds">
                    <text x="8" y="52" fill="#64748b" font-size="10" font-weight="bold">TRACK 2: CDS Gene Models</text>
                    <line x1="0" y1="75" x2="${width}" y2="75" stroke="#1e293b" stroke-width="1" />
        `;

        (data.genes || []).forEach(g => {
            const x1 = posToX(g.start_pos);
            const x2 = posToX(g.end_pos);
            const w = Math.max(12, x2 - x1);
            const isTarget = g.locus_tag.toLowerCase() === this.state.locusTag.toLowerCase();
            const isFwd = g.strand === '+' || g.strand === '1';

            const bg = isTarget ? "#f59e0b" : (isFwd ? "#10b981" : "#6366f1");
            const stroke = isTarget ? "#b45309" : (isFwd ? "#047857" : "#4338ca");

            const labelStr = g.gene_name ? `${g.gene_name} (${g.locus_tag})` : g.locus_tag;

            svg += `
                <g class="cursor-pointer hover:opacity-80 transition-all"
                   onmouseenter="GenomicTrackBrowser.inspectGene('${g.locus_tag}', '${g.gene_name || ''}', ${g.start_pos}, ${g.end_pos}, '${g.strand}', '${(g.promoter_70bp || '').replace(/'/g, '')}')"
                   onclick="window.querySingleGene && window.querySingleGene('${g.locus_tag}')">
                    <rect x="${x1}" y="60" width="${w}" height="24" rx="4" fill="${bg}" stroke="${stroke}" stroke-width="1.5" opacity="0.88" />
                    <text x="${x1 + w / 2}" y="76" fill="#ffffff" font-size="10" font-weight="bold" text-anchor="middle">${labelStr}</text>
                </g>
            `;
        });

        svg += `</g>

                <!-- TRACK 3: TSS & Promoter 70bp Highlights -->
                <g id="track-promoter">
                    <text x="8" y="115" fill="#64748b" font-size="10" font-weight="bold">TRACK 3: Promoter & TSS Sites</text>
                    <line x1="0" y1="135" x2="${width}" y2="135" stroke="#1e293b" stroke-width="1" />
        `;

        (data.genes || []).forEach(g => {
            if (g.tss_position) {
                const x = posToX(g.tss_position);
                svg += `
                    <g class="cursor-pointer" onmouseenter="GenomicTrackBrowser.inspectPromoter('${g.locus_tag}', ${g.tss_position}, '${(g.promoter_70bp || '').replace(/'/g, '')}')">
                        <line x1="${x}" y1="125" x2="${x}" y2="145" stroke="#22c55e" stroke-width="2" stroke-dasharray="2,2" />
                        <path d="M ${x} 125 L ${x + 8} 125 L ${x + 5} 122 M ${x + 8} 125 L ${x + 5} 128" stroke="#22c55e" stroke-width="2" fill="none" />
                        <text x="${x + 10}" y="128" fill="#4ade80" font-size="9" font-family="monospace">TSS: ${g.tss_position}</text>
                    </g>
                `;
            }
        });

        svg += `</g>

                <!-- TRACK 4: ChIP-seq / TFBS Binding Peak Signal Density -->
                <g id="track-peaks">
                    <text x="8" y="170" fill="#64748b" font-size="10" font-weight="bold">TRACK 4: Binding Peak Signal Density (ChIP-seq Peak Envelopes) <tspan fill="#38bdf8">🧪 实验室内网版 (172.16.2.105)</tspan></text>
                    <line x1="0" y1="210" x2="${width}" y2="210" stroke="#1e293b" stroke-width="1" />
        `;

        // Render real continuous ChIP-seq peak envelopes
        if ((data.peaks || []).length > 0) {
            data.peaks.forEach(p => {
                const center = p.peak_center || p.pos || (data.target ? data.target.start_pos : 0);
                const pxCenter = posToX(center);
                const pStart = p.peak_start ? posToX(p.peak_start) : pxCenter - 18;
                const pEnd = p.peak_end ? posToX(p.peak_end) : pxCenter + 18;
                const halfW = Math.max(12, (pEnd - pStart) / 2);

                const score = p.peak_score || p.score || 1.0;
                const negq = p.neglog10q || 0;
                const tier = p.strength_tier || 'moderate';
                const spatialConf = p.spatial_confidence || 'PROMOTER_DIRECT';
                const relTss = p.rel_pos_to_tss != null ? p.rel_pos_to_tss : null;

                let color = '#38bdf8';
                let stroke = '#0284c7';
                let fillGrad = '#0284c7';

                if (tier === 'very_strong') {
                    color = '#f87171'; stroke = '#dc2626'; fillGrad = '#ef4444';
                } else if (tier === 'strong') {
                    color = '#fbbf24'; stroke = '#d97706'; fillGrad = '#f59e0b';
                } else if (tier === 'moderate') {
                    color = '#38bdf8'; stroke = '#0284c7'; fillGrad = '#0284c7';
                } else {
                    color = '#94a3b8'; stroke = '#64748b'; fillGrad = '#64748b';
                }

                const peakH = Math.min(38, Math.max(14, Math.log2(score + 1) * 8 + (negq > 0 ? Math.log10(negq + 1) * 4 : 0)));

                svg += `
                    <g class="cursor-pointer hover:opacity-90 transition-all"
                       onmouseenter="GenomicTrackBrowser.inspectPeakDetails('${(p.tf_name||'TF').replace(/'/g,'')}', '${(p.peak_id||'').replace(/'/g,'')}', ${center}, ${score}, ${negq}, '${tier}', '${spatialConf}', ${relTss}, '${(p.nearest_gene_locus||'').replace(/'/g,'')}')">
                        <path d="M ${pxCenter - halfW} 210 Q ${pxCenter} ${210 - peakH} ${pxCenter + halfW} 210 Z" fill="${color}" opacity="0.65" stroke="${stroke}" stroke-width="1.5" />
                        <circle cx="${pxCenter}" cy="${210 - peakH}" r="3.5" fill="${fillGrad}" stroke="#ffffff" stroke-width="0.8" />
                        <text x="${pxCenter}" y="${194 - peakH}" fill="${color}" font-size="9" text-anchor="middle" font-weight="bold">${p.tf_name || 'TF'}</text>
                    </g>
                `;
            });
        } else {
            svg += `<text x="${width/2}" y="195" fill="#475569" font-size="10" text-anchor="middle">No experimentally verified binding peaks in this window</text>`;
        }

        svg += `</g>

                <!-- TRACK 5: sRNA & ncRNA Annotations -->
                <g id="track-srna">
                    <text x="8" y="235" fill="#64748b" font-size="10" font-weight="bold">TRACK 5: sRNA & ncRNA Tracks</text>
                    <line x1="0" y1="265" x2="${width}" y2="265" stroke="#1e293b" stroke-width="1" />
        `;

        (data.rnas || []).forEach(r => {
            const rx1 = posToX(r.start_pos);
            const rx2 = posToX(r.end_pos);
            const rw = Math.max(14, rx2 - rx1);
            svg += `
                <g class="cursor-pointer">
                    <rect x="${rx1}" y="248" width="${rw}" height="16" rx="8" fill="#a855f7" stroke="#7e22ce" stroke-width="1" opacity="0.85" />
                    <text x="${rx1 + rw / 2}" y="260" fill="#ffffff" font-size="9" text-anchor="middle" font-weight="bold">${r.rna_name || r.rna_id}</text>
                </g>
            `;
        });

        svg += `</g></svg>`;
        wrap.innerHTML = svg;
    },

    inspectGene(locus, name, start, end, strand, promoter) {
        const insp = document.getElementById("gtb-inspector");
        if (!insp) return;
        insp.classList.remove("hidden");
        insp.innerHTML = `
            <div><strong class="text-teal-400">CDS Gene:</strong> ${locus} ${name ? '('+name+')' : ''} | <strong>Coordinates:</strong> ${start.toLocaleString()} - ${end.toLocaleString()} bp (${(end-start).toLocaleString()} bp, Strand: ${strand})</div>
            ${promoter ? `<div class="mt-1 text-slate-400"><strong>Promoter 70bp:</strong> <span class="text-amber-300 font-mono">${promoter}</span></div>` : ''}
        `;
    },

    inspectPromoter(locus, tss, promoter) {
        const insp = document.getElementById("gtb-inspector");
        if (!insp) return;
        insp.classList.remove("hidden");
        insp.innerHTML = `
            <div><strong class="text-emerald-400">TSS Position:</strong> ${locus} @ ${tss.toLocaleString()} bp</div>
            <div class="mt-1 text-slate-400"><strong>Promoter (-35 / -10 box 70bp):</strong> <span class="text-emerald-300 font-mono">${promoter || 'N/A'}</span></div>
        `;
    },

    inspectPeak(tf, pos, score, seq) {
        this.inspectPeakDetails(tf, '', pos, score, 0, 'moderate', 'PROMOTER_DIRECT', 0, '');
    },

    inspectPeakDetails(tf, peakId, center, score, negq, tier, spatialConf, relTss, nearestGene) {
        const insp = document.getElementById("gtb-inspector");
        if (!insp) return;
        insp.classList.remove("hidden");

        const tierBadge = tier === 'very_strong' ? '<span class="bg-rose-950 text-rose-300 stroke-rose-700 px-1.5 py-0.5 rounded text-[10px] border border-rose-700">VERY STRONG</span>'
            : (tier === 'strong' ? '<span class="bg-amber-950 text-amber-300 border border-amber-700 px-1.5 py-0.5 rounded text-[10px]">STRONG</span>'
            : '<span class="bg-sky-950 text-sky-300 border border-sky-700 px-1.5 py-0.5 rounded text-[10px]">MODERATE</span>');

        const confBadge = spatialConf === 'PROMOTER_DIRECT' ? '<span class="bg-emerald-950 text-emerald-300 border border-emerald-700 px-1.5 py-0.5 rounded text-[10px]">🎯 PROMOTER DIRECT</span>'
            : (spatialConf === 'INTERGENIC_PROMOTER' ? '<span class="bg-teal-950 text-teal-300 border border-teal-700 px-1.5 py-0.5 rounded text-[10px]">🧬 INTERGENIC PROMOTER</span>'
            : '<span class="bg-slate-800 text-slate-300 border border-slate-700 px-1.5 py-0.5 rounded text-[10px]">GENE BODY INTERNAL</span>');

        const tssStr = relTss != null ? `${relTss > 0 ? '+' : ''}${relTss} bp to TSS` : 'Distal';

        insp.innerHTML = `
            <div class="flex items-center gap-2 flex-wrap">
                <strong class="text-sky-400">ChIP-seq Binding Peak:</strong> <span class="text-slate-100 font-bold">${tf}</span> (${peakId || 'Peak'}) @ <span class="text-amber-300 font-mono">${center.toLocaleString()} bp</span>
                ${tierBadge} ${confBadge}
            </div>
            <div class="mt-1 text-slate-300 flex items-center gap-4 text-[11px] flex-wrap font-mono">
                <span>Signal Enrichment: <strong class="text-sky-300">${score.toFixed(2)}x</strong></span>
                <span>-log10(q): <strong class="text-indigo-300">${negq.toFixed(1)}</strong></span>
                <span>TSS Offset: <strong class="text-emerald-300">${tssStr}</strong></span>
                <span>Target: <strong class="text-amber-300">${nearestGene}</strong></span>
            </div>
        `;
    }
};
