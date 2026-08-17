/**
 * genomicTrackBrowser.js
 * =======================
 * Premium Light-Themed 5-Track & ChIP-seq Peak Inspector for C. glutamicum.
 * Features:
 *  - Light theme styling matching the main app design system
 *  - Structured ChIP-seq Binding Peak summary table
 *  - Non-overlapping, high-legibility vector tracks (CDS, TSS, ChIP-seq, sRNA)
 *  - Interactive hover inspector with clear scientific metrics
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
            <div style="background: #ffffff; border: 1px solid var(--border-color, #e2e8f0); border-radius: 10px; padding: 14px; font-family: var(--font-primary, sans-serif); color: var(--text-primary, #0f172a); box-shadow: 0 2px 10px rgba(15, 23, 42, 0.04); display: flex; flex-direction: column; gap: 12px;">
                <!-- Header -->
                <div style="display: flex; align-items: center; justify-content: space-between; border-bottom: 1px solid var(--border-color, #e2e8f0); padding-bottom: 10px;">
                    <div style="display: flex; align-items: center; gap: 8px;">
                        <span style="font-size: 14px; color: #6366f1;">🧬</span>
                        <h4 style="margin: 0; font-size: 13px; font-weight: 700; color: var(--text-primary, #0f172a);">Genomic Track & Peak Inspector</h4>
                        <span style="font-size: 11px; background: rgba(99, 102, 241, 0.08); color: #6366f1; padding: 2px 8px; border-radius: 4px; font-family: monospace; font-weight: 600;">${locusTag}</span>
                    </div>
                    <div style="display: flex; align-items: center; gap: 6px;">
                        <button onclick="GenomicTrackBrowser.zoom(0.6)" class="secondary-btn" style="padding: 2px 8px; font-size: 10.5px; height: 24px; border-radius: 4px; border: 1px solid var(--border-color, #cbd5e1); background: #ffffff; cursor: pointer;" title="Zoom in">🔍 In</button>
                        <button onclick="GenomicTrackBrowser.zoom(1.5)" class="secondary-btn" style="padding: 2px 8px; font-size: 10.5px; height: 24px; border-radius: 4px; border: 1px solid var(--border-color, #cbd5e1); background: #ffffff; cursor: pointer;" title="Zoom out">🔍 Out</button>
                        <button onclick="GenomicTrackBrowser.shift(-3000)" class="secondary-btn" style="padding: 2px 8px; font-size: 10.5px; height: 24px; border-radius: 4px; border: 1px solid var(--border-color, #cbd5e1); background: #ffffff; cursor: pointer;" title="Pan left">◀ Left</button>
                        <button onclick="GenomicTrackBrowser.shift(3000)" class="secondary-btn" style="padding: 2px 8px; font-size: 10.5px; height: 24px; border-radius: 4px; border: 1px solid var(--border-color, #cbd5e1); background: #ffffff; cursor: pointer;" title="Pan right">Right ▶</button>
                    </div>
                </div>

                <!-- ChIP-seq Peak Summary Card Table -->
                <div id="gtb-peaks-summary-box" style="display:none; border: 1px solid #e2e8f0; border-radius: 8px; background: #f8fafc; padding: 10px; overflow-x: auto;">
                    <div style="font-size: 11px; font-weight: 700; color: #475569; margin-bottom: 6px; display: flex; align-items: center; justify-content: space-between;">
                        <span><i class="fa-solid fa-chart-line" style="color:#0284c7;"></i> Mapped ChIP-seq Binding Peaks</span>
                        <span id="gtb-peaks-count" style="font-size: 10px; font-weight: 600; color: #0284c7;">-</span>
                    </div>
                    <table style="width: 100%; border-collapse: collapse; text-align: left; font-size: 10.5px;">
                        <thead>
                            <tr style="border-bottom: 1px solid #cbd5e1; color: #64748b;">
                                <th style="padding: 4px 6px;">Regulator / Peak</th>
                                <th style="padding: 4px 6px;">Genomic Position</th>
                                <th style="padding: 4px 6px;">Enrichment</th>
                                <th style="padding: 4px 6px;">Confidence</th>
                                <th style="padding: 4px 6px;">Target</th>
                            </tr>
                        </thead>
                        <tbody id="gtb-peaks-table-body" style="font-family: monospace;">
                        </tbody>
                    </table>
                </div>

                <!-- 5-Track SVG Canvas Wrapper -->
                <div id="gtb-canvas-wrap" style="position: relative; overflow: hidden; background: #f8fafc; border-radius: 8px; border: 1px solid #e2e8f0; padding: 4px;">
                    <div style="text-align: center; font-size: 11px; color: #64748b; padding: 24px 0;"><i class="fa-solid fa-spinner fa-spin"></i> Loading genomic tracks for ${locusTag}...</div>
                </div>

                <!-- Hover Details Inspector -->
                <div id="gtb-inspector" style="font-size: 11px; background: #f1f5f9; border: 1px solid #cbd5e1; border-radius: 6px; padding: 8px 12px; color: #0f172a; font-family: monospace; display: none;">
                    Hover over a gene, promoter site, or binding peak to inspect details.
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
            this.drawPeaksTable();
            this.drawTracks();
        } catch (err) {
            wrap.innerHTML = `<div style="color: #ef4444; font-size: 11px; padding: 16px; text-align: center;">Failed to load genomic tracks: ${err.message}</div>`;
        }
    },

    drawPeaksTable() {
        const box = document.getElementById("gtb-peaks-summary-box");
        const tbody = document.getElementById("gtb-peaks-table-body");
        const countSpan = document.getElementById("gtb-peaks-count");
        if (!box || !tbody || !this.state.data) return;

        const peaks = this.state.data.peaks || [];
        if (peaks.length === 0) {
            box.style.display = "none";
            return;
        }

        box.style.display = "block";
        if (countSpan) countSpan.textContent = `${peaks.length} peak(s)`;

        let rowsHtml = "";
        peaks.slice(0, 5).forEach(p => {
            const center = p.peak_center || p.pos || 0;
            const score = p.peak_score || p.score || 1.0;
            const tier = p.strength_tier || "moderate";
            const relTss = p.rel_pos_to_tss != null ? p.rel_pos_to_tss : null;
            const target = p.nearest_gene_locus || p.locus_tag || "-";
            const tf = p.tf_name || p.tf_id || "TF";

            let tierBg = "#e0f2fe";
            let tierColor = "#0369a1";
            if (tier === "very_strong") { tierBg = "#ffe4e6"; tierColor = "#be123c"; }
            else if (tier === "strong") { tierBg = "#fef3c7"; tierColor = "#b45309"; }

            const tssBadge = relTss != null ? `${relTss > 0 ? '+' : ''}${relTss} bp` : 'Distal';

            rowsHtml += `
                <tr style="border-bottom: 1px solid #f1f5f9;">
                    <td style="padding: 4px 6px; font-weight: 700; color: #0284c7;">${tf}</td>
                    <td style="padding: 4px 6px; color: #334155;">${center.toLocaleString()} bp <span style="font-size:9px; color:#64748b;">(${tssBadge})</span></td>
                    <td style="padding: 4px 6px; font-weight: 700; color: #0f172a;">${score.toFixed(2)}x</td>
                    <td style="padding: 4px 6px;"><span style="font-size:9px; background:${tierBg}; color:${tierColor}; padding:1px 5px; border-radius:3px; font-weight:700;">${tier.toUpperCase()}</span></td>
                    <td style="padding: 4px 6px; color: #475569;">${target}</td>
                </tr>
            `;
        });
        tbody.innerHTML = rowsHtml;
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
        const width = wrap.clientWidth || 340;
        const posToX = (pos) => Math.max(10, Math.min(width - 10, ((pos - minPos) / totalBp) * (width - 20) + 10));

        let svg = `
            <svg width="100%" height="240" viewBox="0 0 ${width} 240" xmlns="http://www.w3.org/2000/svg" style="display:block; user-select:none; font-family: sans-serif;">
                <!-- Background -->
                <rect x="0" y="0" width="${width}" height="240" fill="#ffffff" rx="6" />

                <!-- TRACK 1: Ruler -->
                <g id="track-ruler">
                    <line x1="10" y1="26" x2="${width - 10}" y2="26" stroke="#cbd5e1" stroke-width="1.5" />
        `;

        // Calculate tick step size cleanly (at least 70px spacing)
        const minPixelSpacing = 70;
        const rawStepBp = (totalBp / width) * minPixelSpacing;
        const tickStepBp = Math.pow(10, Math.floor(Math.log10(rawStepBp))) * (rawStepBp / Math.pow(10, Math.floor(Math.log10(rawStepBp))) > 5 ? 5 : (rawStepBp / Math.pow(10, Math.floor(Math.log10(rawStepBp))) > 2 ? 2 : 1));
        const firstTick = Math.ceil(minPos / tickStepBp) * tickStepBp;

        for (let pos = firstTick; pos <= maxPos; pos += tickStepBp) {
            const x = posToX(pos);
            svg += `
                <line x1="${x}" y1="20" x2="${x}" y2="26" stroke="#94a3b8" stroke-width="1" />
                <text x="${x}" y="15" fill="#64748b" font-size="8.5" text-anchor="middle" font-family="monospace">${pos.toLocaleString()} bp</text>
            `;
        }

        svg += `
                </g>

                <!-- TRACK 2: CDS Gene Models -->
                <g id="track-cds">
                    <text x="10" y="44" fill="#475569" font-size="9" font-weight="700">TRACK 2: CDS Gene Models</text>
                    <line x1="10" y1="65" x2="${width - 10}" y2="65" stroke="#e2e8f0" stroke-width="1" />
        `;

        (data.genes || []).forEach(g => {
            const x1 = posToX(g.start_pos);
            const x2 = posToX(g.end_pos);
            const w = Math.max(14, x2 - x1);
            const isTarget = g.locus_tag.toLowerCase() === this.state.locusTag.toLowerCase();
            const isFwd = g.strand === '+' || g.strand === '1';

            const fill = isTarget ? "#f59e0b" : (isFwd ? "#10b981" : "#6366f1");
            const stroke = isTarget ? "#b45309" : (isFwd ? "#047857" : "#4338ca");

            const dispName = (g.gene_name && g.gene_name !== '--') ? g.gene_name : g.locus_tag;
            const truncatedName = dispName.length > 7 ? dispName.substring(0, 6) + '..' : dispName;

            svg += `
                <g style="cursor:pointer;"
                   onmouseenter="GenomicTrackBrowser.inspectGene('${g.locus_tag}', '${g.gene_name || ''}', ${g.start_pos}, ${g.end_pos}, '${g.strand}', '${(g.promoter_70bp || '').replace(/'/g, '')}')"
                   onclick="window.querySingleGene && window.querySingleGene('${g.locus_tag}')">
                    <rect x="${x1}" y="52" width="${w}" height="20" rx="3" fill="${fill}" stroke="${stroke}" stroke-width="1.2" opacity="0.9" />
                    <text x="${x1 + w / 2}" y="65" fill="#ffffff" font-size="8.5" font-weight="700" text-anchor="middle">${truncatedName}</text>
                </g>
            `;
        });

        svg += `</g>

                <!-- TRACK 3: TSS & Promoter 70bp Highlights -->
                <g id="track-promoter">
                    <text x="10" y="98" fill="#475569" font-size="9" font-weight="700">TRACK 3: Promoter & TSS Sites</text>
                    <line x1="10" y1="118" x2="${width - 10}" y2="118" stroke="#e2e8f0" stroke-width="1" />
        `;

        (data.genes || []).forEach(g => {
            if (g.tss_position) {
                const x = posToX(g.tss_position);
                svg += `
                    <g style="cursor:pointer;" onmouseenter="GenomicTrackBrowser.inspectPromoter('${g.locus_tag}', ${g.tss_position}, '${(g.promoter_70bp || '').replace(/'/g, '')}')">
                        <line x1="${x}" y1="108" x2="${x}" y2="128" stroke="#16a34a" stroke-width="1.5" stroke-dasharray="2,2" />
                        <path d="M ${x} 108 L ${x + 6} 108 L ${x + 4} 105 M ${x + 6} 108 L ${x + 4} 111" stroke="#16a34a" stroke-width="1.5" fill="none" />
                        <text x="${x + 8}" y="112" fill="#15803d" font-size="8" font-family="monospace">TSS: ${g.tss_position}</text>
                    </g>
                `;
            }
        });

        svg += `</g>

                <!-- TRACK 4: ChIP-seq Binding Peak Signal Density -->
                <g id="track-peaks">
                    <text x="10" y="148" fill="#475569" font-size="9" font-weight="700">TRACK 4: Binding Peak Signal Density (ChIP-seq)</text>
                    <line x1="10" y1="182" x2="${width - 10}" y2="182" stroke="#e2e8f0" stroke-width="1" />
        `;

        if ((data.peaks || []).length > 0) {
            data.peaks.forEach(p => {
                const center = p.peak_center || p.pos || (data.target ? data.target.start_pos : 0);
                const pxCenter = posToX(center);
                const pStart = p.peak_start ? posToX(p.peak_start) : pxCenter - 14;
                const pEnd = p.peak_end ? posToX(p.peak_end) : pxCenter + 14;
                const halfW = Math.max(10, (pEnd - pStart) / 2);

                const score = p.peak_score || p.score || 1.0;
                const negq = p.neglog10q || 0;
                const tier = p.strength_tier || 'moderate';
                const spatialConf = p.spatial_confidence || 'PROMOTER_DIRECT';
                const relTss = p.rel_pos_to_tss != null ? p.rel_pos_to_tss : null;

                let color = '#0284c7';
                let fillGrad = 'rgba(2, 132, 199, 0.25)';

                if (tier === 'very_strong') {
                    color = '#e11d48'; fillGrad = 'rgba(225, 29, 72, 0.25)';
                } else if (tier === 'strong') {
                    color = '#d97706'; fillGrad = 'rgba(217, 119, 6, 0.25)';
                }

                const peakH = Math.min(30, Math.max(12, Math.log2(score + 1) * 6 + (negq > 0 ? Math.log10(negq + 1) * 3 : 0)));

                svg += `
                    <g style="cursor:pointer;"
                       onmouseenter="GenomicTrackBrowser.inspectPeakDetails('${(p.tf_name||'TF').replace(/'/g,'')}', '${(p.peak_id||'').replace(/'/g,'')}', ${center}, ${score}, ${negq}, '${tier}', '${spatialConf}', ${relTss}, '${(p.nearest_gene_locus||'').replace(/'/g,'')}')">
                        <path d="M ${pxCenter - halfW} 182 Q ${pxCenter} ${182 - peakH} ${pxCenter + halfW} 182 Z" fill="${fillGrad}" stroke="${color}" stroke-width="1.3" />
                        <circle cx="${pxCenter}" cy="${182 - peakH}" r="3" fill="${color}" stroke="#ffffff" stroke-width="0.8" />
                        <text x="${pxCenter}" y="${170 - peakH}" fill="${color}" font-size="8.5" text-anchor="middle" font-weight="700">${p.tf_name || 'TF'}</text>
                    </g>
                `;
            });
        } else {
            svg += `<text x="${width/2}" y="170" fill="#94a3b8" font-size="9" text-anchor="middle">No binding peaks mapped in this window</text>`;
        }

        svg += `</g>

                <!-- TRACK 5: sRNA Tracks -->
                <g id="track-srna">
                    <text x="10" y="202" fill="#475569" font-size="9" font-weight="700">TRACK 5: sRNA & ncRNA Annotations</text>
                    <line x1="10" y1="228" x2="${width - 10}" y2="228" stroke="#e2e8f0" stroke-width="1" />
        `;

        (data.rnas || []).forEach(r => {
            const rx1 = posToX(r.start_pos);
            const rx2 = posToX(r.end_pos);
            const rw = Math.max(12, rx2 - rx1);
            svg += `
                <g style="cursor:pointer;">
                    <rect x="${rx1}" y="214" width="${rw}" height="13" rx="6" fill="#9333ea" stroke="#7e22ce" stroke-width="1" opacity="0.85" />
                    <text x="${rx1 + rw / 2}" y="224" fill="#ffffff" font-size="8" text-anchor="middle" font-weight="700">${r.rna_name || r.rna_id}</text>
                </g>
            `;
        });

        svg += `</g></svg>`;
        wrap.innerHTML = svg;
    },

    inspectGene(locus, name, start, end, strand, promoter) {
        const insp = document.getElementById("gtb-inspector");
        if (!insp) return;
        insp.style.display = "block";
        insp.innerHTML = `
            <div><strong style="color: #0d9488;">CDS Gene:</strong> ${locus} ${name ? '('+name+')' : ''} | <strong>Coordinates:</strong> ${start.toLocaleString()} - ${end.toLocaleString()} bp (${(end-start).toLocaleString()} bp, Strand: ${strand})</div>
            ${promoter ? `<div style="margin-top: 3px; color: #475569;"><strong>Promoter 70bp:</strong> <span style="color: #b45309; font-family: monospace;">${promoter}</span></div>` : ''}
        `;
    },

    inspectPromoter(locus, tss, promoter) {
        const insp = document.getElementById("gtb-inspector");
        if (!insp) return;
        insp.style.display = "block";
        insp.innerHTML = `
            <div><strong style="color: #16a34a;">TSS Position:</strong> ${locus} @ ${tss.toLocaleString()} bp</div>
            <div style="margin-top: 3px; color: #475569;"><strong>Promoter (-35 / -10 box 70bp):</strong> <span style="color: #15803d; font-family: monospace;">${promoter || 'N/A'}</span></div>
        `;
    },

    inspectPeak(tf, pos, score, seq) {
        this.inspectPeakDetails(tf, '', pos, score, 0, 'moderate', 'PROMOTER_DIRECT', 0, '');
    },

    inspectPeakDetails(tf, peakId, center, score, negq, tier, spatialConf, relTss, nearestGene) {
        const insp = document.getElementById("gtb-inspector");
        if (!insp) return;
        insp.style.display = "block";

        const tierBadge = tier === 'very_strong' ? '<span style="background:#ffe4e6; color:#be123c; padding:1px 5px; border-radius:3px; font-size:9px; font-weight:700;">VERY STRONG</span>'
            : (tier === 'strong' ? '<span style="background:#fef3c7; color:#b45309; padding:1px 5px; border-radius:3px; font-size:9px; font-weight:700;">STRONG</span>'
            : '<span style="background:#e0f2fe; color:#0369a1; padding:1px 5px; border-radius:3px; font-size:9px; font-weight:700;">MODERATE</span>');

        const confBadge = spatialConf === 'PROMOTER_DIRECT' ? '<span style="background:#dcfce7; color:#15803d; padding:1px 5px; border-radius:3px; font-size:9px; font-weight:700;">🎯 PROMOTER DIRECT</span>'
            : '<span style="background:#f1f5f9; color:#475569; padding:1px 5px; border-radius:3px; font-size:9px; font-weight:700;">INTERGENIC</span>';

        const tssStr = relTss != null ? `${relTss > 0 ? '+' : ''}${relTss} bp to TSS` : 'Distal';

        insp.innerHTML = `
            <div style="display:flex; align-items:center; gap:6px; flex-wrap:wrap; margin-bottom:4px;">
                <strong style="color: #0284c7;">ChIP-seq Peak:</strong> <span style="font-weight:700; color:#0f172a;">${tf}</span> (${peakId || 'Peak'}) @ <span style="color:#b45309;">${center.toLocaleString()} bp</span>
                ${tierBadge} ${confBadge}
            </div>
            <div style="color: #334155; display:flex; align-items:center; gap:10px; font-size:10.5px; flex-wrap:wrap;">
                <span>Signal Fold: <strong style="color:#0284c7;">${score.toFixed(2)}x</strong></span>
                <span>-log10(q): <strong style="color:#4f46e5;">${negq.toFixed(1)}</strong></span>
                <span>TSS Offset: <strong style="color:#16a34a;">${tssStr}</strong></span>
                <span>Target: <strong style="color:#b45309;">${nearestGene}</strong></span>
            </div>
        `;
    }
};
