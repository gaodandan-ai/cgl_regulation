/**
 * networkTopology.js
 * Network topology analysis and motif identification for C. glutamicum regulatory network.
 * Pure client-side JavaScript, no external dependencies.
 */
(function (global) {
    'use strict';

    // ── Adjacency Builder ────────────────────────────────────────────────────

    function buildAdjacency(edges) {
        const outAdj = {}, inAdj = {}, edgeMap = {}, names = {};
        edges.forEach(e => {
            const { src, tgt, role, evidence, tfName, tgName } = e;
            if (!src || !tgt) return;
            if (!outAdj[src]) outAdj[src] = new Set();
            if (!inAdj[tgt]) inAdj[tgt] = new Set();
            outAdj[src].add(tgt);
            inAdj[tgt].add(src);
            const key = `${src}|${tgt}`;
            if (!edgeMap[key]) edgeMap[key] = { role, evidence };
            if (tfName) names[src] = tfName;
            if (tgName) names[tgt] = tgName;
        });
        return { outAdj, inAdj, edgeMap, names };
    }

    // ── Degree Distribution ───────────────────────────────────────────────────

    function computeDegrees(edges, outAdj, inAdj) {
        const allNodes = new Set();
        edges.forEach(e => { if (e.src) allNodes.add(e.src); if (e.tgt) allNodes.add(e.tgt); });
        const outDeg = {}, inDeg = {};
        allNodes.forEach(n => {
            outDeg[n] = outAdj[n] ? outAdj[n].size : 0;
            inDeg[n] = inAdj[n] ? inAdj[n].size : 0;
        });
        const outDist = {}, inDist = {};
        Object.values(outDeg).forEach(d => { outDist[d] = (outDist[d] || 0) + 1; });
        Object.values(inDeg).forEach(d => { inDist[d] = (inDist[d] || 0) + 1; });
        return { allNodes: [...allNodes], outDeg, inDeg, outDist, inDist };
    }

    // ── Hub TF Ranking ────────────────────────────────────────────────────────

    function computeHubTFs(edges, outAdj, inAdj, edgeMap, names) {
        const tfs = new Set(edges.map(e => e.src).filter(Boolean));
        const result = [];
        tfs.forEach(tf => {
            const targets = outAdj[tf] || new Set();
            let activationCount = 0, repressionCount = 0, autoRole = null;
            targets.forEach(tgt => {
                const e = edgeMap[`${tf}|${tgt}`] || {};
                const r = (e.role || '').toUpperCase();
                if (r === 'A') activationCount++;
                else if (r === 'R') repressionCount++;
                if (tgt === tf) autoRole = r;
            });
            result.push({
                locus: tf,
                name: names[tf] || tf,
                outDegree: targets.size,
                inDegree: (inAdj[tf] || new Set()).size,
                activationCount,
                repressionCount,
                isAutoRegulated: targets.has(tf),
                autoRole
            });
        });
        result.sort((a, b) => b.outDegree - a.outDegree);
        return result;
    }

    // ── Betweenness Centrality (sampled BFS) ──────────────────────────────────

    function computeBetweennessCentrality(allNodes, outAdj, maxSources) {
        maxSources = maxSources || 150;
        const bc = {};
        allNodes.forEach(n => { bc[n] = 0; });
        const sources = allNodes.length > maxSources ? allNodes.slice(0, maxSources) : allNodes;

        sources.forEach(s => {
            const pred = {}, dist = {}, sigma = {}, delta = {};
            allNodes.forEach(n => { pred[n] = []; sigma[n] = 0; dist[n] = -1; delta[n] = 0; });
            dist[s] = 0; sigma[s] = 1;
            const queue = [s], stack = [];
            while (queue.length > 0) {
                const v = queue.shift();
                stack.push(v);
                (outAdj[v] ? [...outAdj[v]] : []).forEach(w => {
                    if (dist[w] < 0) { queue.push(w); dist[w] = dist[v] + 1; }
                    if (dist[w] === dist[v] + 1) { sigma[w] += sigma[v]; pred[w].push(v); }
                });
            }
            while (stack.length > 0) {
                const w = stack.pop();
                pred[w].forEach(v => { delta[v] += (sigma[v] / sigma[w]) * (1 + delta[w]); });
                if (w !== s) bc[w] += delta[w];
            }
        });
        const n = allNodes.length;
        const scale = n > 2 ? 1 / ((n - 1) * (n - 2)) : 1;
        Object.keys(bc).forEach(k => { bc[k] *= scale; });
        return bc;
    }

    // ── Autoregulation ────────────────────────────────────────────────────────

    function findAutoRegulations(edges, edgeMap, names, outAdj) {
        const seen = new Set(), result = [];
        edges.forEach(e => {
            if (e.src && e.src === e.tgt && !seen.has(e.src)) {
                seen.add(e.src);
                const info = edgeMap[`${e.src}|${e.src}`] || {};
                result.push({
                    locus: e.src, name: names[e.src] || e.src,
                    role: info.role || '', evidence: info.evidence || '',
                    outDegree: outAdj[e.src] ? outAdj[e.src].size : 0
                });
            }
        });
        result.sort((a, b) => b.outDegree - a.outDegree);
        return result;
    }

    // ── Feed-forward Loops ───────────────────────────────────────────────────

    function classifyFFL(roleAB, roleBC, roleAC) {
        const r = v => (v || '').toUpperCase() === 'A' ? 1 : -1;
        const isCoherent = (r(roleAB) * r(roleBC)) === r(roleAC);
        const prefix = isCoherent ? 'C' : 'I';
        const subtype = `${prefix}${(roleAB||'?')}${(roleBC||'?')}${(roleAC||'?')}`;
        return { isCoherent, subtype };
    }

    function findFFLs(edges, outAdj, inAdj, edgeMap, names) {
        const tfs = new Set(edges.map(e => e.src).filter(Boolean));
        const byMasterTF = {};

        tfs.forEach(a => {
            const aTargets = outAdj[a] || new Set();
            const ffls = [];
            aTargets.forEach(b => {
                if (b === a || !tfs.has(b)) return;
                const bTargets = outAdj[b] || new Set();
                aTargets.forEach(c => {
                    if (c === a || c === b || !bTargets.has(c)) return;
                    const roleAB = (edgeMap[`${a}|${b}`] || {}).role || '';
                    const roleBC = (edgeMap[`${b}|${c}`] || {}).role || '';
                    const roleAC = (edgeMap[`${a}|${c}`] || {}).role || '';
                    const { isCoherent, subtype } = classifyFFL(roleAB, roleBC, roleAC);
                    ffls.push({
                        masterTF: a, masterTFName: names[a] || a,
                        intermediateTF: b, intermediateTFName: names[b] || b,
                        target: c, targetName: names[c] || c,
                        roleAB, roleBC, roleAC,
                        type: isCoherent ? 'coherent' : 'incoherent',
                        subtype, isCoherent
                    });
                });
            });
            if (ffls.length > 0) {
                byMasterTF[a] = {
                    locus: a, name: names[a] || a,
                    outDegree: aTargets.size,
                    ffls,
                    coherentCount: ffls.filter(f => f.isCoherent).length,
                    incoherentCount: ffls.filter(f => !f.isCoherent).length
                };
            }
        });

        const sorted = Object.values(byMasterTF).sort((a, b) => b.ffls.length - a.ffls.length);
        const totalFFLs = sorted.reduce((s, g) => s + g.ffls.length, 0);
        const totalCoherent = sorted.reduce((s, g) => s + g.coherentCount, 0);
        return {
            byMasterTF: sorted, totalFFLs,
            totalCoherent, totalIncoherent: totalFFLs - totalCoherent
        };
    }

    // ── Mutual Regulation ─────────────────────────────────────────────────────

    function findMutualRegulation(edges, edgeMap, names) {
        const edgeSet = new Set(edges.map(e => `${e.src}|${e.tgt}`));
        const seen = new Set(), result = [];
        edges.forEach(e => {
            const { src: a, tgt: b } = e;
            if (!a || !b || a === b) return;
            const key = a < b ? `${a}::${b}` : `${b}::${a}`;
            if (seen.has(key)) return;
            if (edgeSet.has(`${b}|${a}`)) {
                seen.add(key);
                const eAB = edgeMap[`${a}|${b}`] || {};
                const eBA = edgeMap[`${b}|${a}`] || {};
                result.push({
                    nodeA: a, nameA: names[a] || a, roleAB: eAB.role || '', evidenceAB: eAB.evidence || '',
                    nodeB: b, nameB: names[b] || b, roleBA: eBA.role || '', evidenceBA: eBA.evidence || ''
                });
            }
        });
        return result;
    }

    // ── Multi-input Motifs ────────────────────────────────────────────────────

    function findMultiInputMotifs(edges, inAdj, edgeMap, names, minTF) {
        minTF = minTF || 3;
        const result = [];
        Object.entries(inAdj).forEach(([gene, regulators]) => {
            const regArr = [...regulators];
            if (regArr.length < minTF) return;
            const tfs = regArr.map(tf => {
                const e = edgeMap[`${tf}|${gene}`] || {};
                return { locus: tf, name: names[tf] || tf, role: e.role || '', evidence: e.evidence || '' };
            }).sort((a, b) => {
                if (a.role !== b.role) return a.role === 'A' ? -1 : 1;
                return (a.name || a.locus).localeCompare(b.name || b.locus);
            });
            result.push({ gene, geneName: names[gene] || gene, tfCount: regArr.length, tfs });
        });
        result.sort((a, b) => b.tfCount - a.tfCount);
        return result;
    }

    // ── Bi-fan Patterns ───────────────────────────────────────────────────────

    function findBiFans(edges, outAdj, names, topN) {
        topN = topN || 100;
        const tfs = [...new Set(edges.map(e => e.src).filter(Boolean))];
        const result = [];
        for (let i = 0; i < tfs.length; i++) {
            for (let j = i + 1; j < tfs.length; j++) {
                const a = tfs[i], b = tfs[j];
                const aT = outAdj[a] || new Set();
                const bT = outAdj[b] || new Set();
                const shared = [...aT].filter(t => bT.has(t) && t !== a && t !== b);
                if (shared.length >= 2) {
                    result.push({
                        tfA: a, nameA: names[a] || a,
                        tfB: b, nameB: names[b] || b,
                        sharedCount: shared.length,
                        sharedTargets: shared.slice(0, 5).map(t => ({ locus: t, name: names[t] || t }))
                    });
                }
            }
        }
        result.sort((a, b) => b.sharedCount - a.sharedCount);
        return result.slice(0, topN);
    }

    // ── Main Entry ────────────────────────────────────────────────────────────

    function normaliseEdges(rawEdges) {
        return rawEdges.map(e => ({
            src: (e.source || e.TF_locusTag || e.tf_locus || e.from || '').trim(),
            tgt: (e.target || e.TG_locusTag || e.tg_locus || e.to || '').trim(),
            role: (e.role || e.Role || e.regulationType || '').trim().toUpperCase(),
            evidence: (e.evidence || e.Evidence || '').trim(),
            tfName: (e.tf_name || e.TF_name || e.tfName || '').trim(),
            tgName: (e.tg_name || e.TG_name || e.tgName || '').trim()
        })).filter(e => e.src && e.tgt);
    }

    function getTopologyReport(rawEdges, options) {
        options = options || {};
        const useBetweenness = options.useBetweenness !== false; // default true
        const maxCentralitySources = options.maxCentralitySources || 150;

        const edges = normaliseEdges(rawEdges);
        const { outAdj, inAdj, edgeMap, names } = buildAdjacency(edges);
        const { allNodes, outDeg, inDeg, outDist, inDist } = computeDegrees(edges, outAdj, inAdj);
        const tfs = [...new Set(edges.map(e => e.src))];

        const activationEdges = edges.filter(e => e.role === 'A').length;
        const repressionEdges = edges.filter(e => e.role === 'R').length;
        const experimentalEdges = edges.filter(e => e.evidence.toLowerCase().includes('experimental')).length;
        const predictedEdges = edges.filter(e => e.evidence.toLowerCase().includes('predicted') && !e.evidence.toLowerCase().includes('experimental')).length;

        const hubTFs = computeHubTFs(edges, outAdj, inAdj, edgeMap, names);
        const autoRegs = findAutoRegulations(edges, edgeMap, names, outAdj);
        const fflResult = findFFLs(edges, outAdj, inAdj, edgeMap, names);
        const mutualRegs = findMutualRegulation(edges, edgeMap, names);
        const multiInputs = findMultiInputMotifs(edges, inAdj, edgeMap, names, 3);
        const biFans = findBiFans(edges, outAdj, names, 100);

        let betweenness = null;
        if (useBetweenness) {
            betweenness = computeBetweennessCentrality(allNodes, outAdj, maxCentralitySources);
        }

        return {
            summary: {
                nodeCount: allNodes.length,
                edgeCount: edges.length,
                tfCount: tfs.length,
                avgOutDegree: tfs.length > 0
                    ? (tfs.reduce((s, t) => s + (outAdj[t] ? outAdj[t].size : 0), 0) / tfs.length).toFixed(2)
                    : '0',
                activationEdges, repressionEdges,
                experimentalEdges, predictedEdges,
                autoRegCount: autoRegs.length,
                fflCount: fflResult.totalFFLs,
                coherentFFL: fflResult.totalCoherent,
                incoherentFFL: fflResult.totalIncoherent,
                mutualRegCount: mutualRegs.length,
                multiInputCount: multiInputs.length,
                biFanCount: biFans.length
            },
            outDeg, inDeg, outDist, inDist,
            hubTFs, autoRegs,
            fflsByMasterTF: fflResult.byMasterTF,
            mutualRegs, multiInputs, biFans,
            betweenness, names
        };
    }

    global.networkTopology = {
        getTopologyReport,
        buildAdjacency,
        computeDegrees,
        computeHubTFs,
        computeBetweennessCentrality,
        findAutoRegulations,
        findFFLs,
        findMutualRegulation,
        findMultiInputMotifs,
        findBiFans
    };
})(window);
