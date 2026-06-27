(function (global) {
    const adapter = global.metabolicModelAdapter || {};
    const normalizeGeneId = adapter.normalizeGeneId || function (id) {
        return String(id || '').trim().replace(/^gene[:_]/i, '').replace(/^G_/i, '').toLowerCase();
    };

    function clamp(value, min = 0, max = 1) {
        const parsed = Number(value);
        if (!Number.isFinite(parsed)) return min;
        return Math.max(min, Math.min(max, parsed));
    }

    function toArray(value) {
        if (!value) return [];
        if (Array.isArray(value)) return value;
        if (typeof value.toArray === 'function') return value.toArray();
        if (typeof value.forEach === 'function') {
            const items = [];
            value.forEach(item => items.push(item));
            return items;
        }
        return [];
    }

    function callOrRead(value, key) {
        if (!value) return undefined;
        if (typeof value.data === 'function') return value.data(key);
        if (value.data && key in value.data) return value.data[key];
        if (key in value) return value[key];
        return undefined;
    }

    function collectNodes(graph) {
        if (!graph) return [];
        if (typeof graph.nodes === 'function') return toArray(graph.nodes());
        if (graph.nodes) return toArray(graph.nodes);
        if (graph.elements?.nodes) return toArray(graph.elements.nodes);
        if (Array.isArray(graph.elements)) return graph.elements.filter(item => item.group === 'nodes' || item.data?.id);
        return [];
    }

    function collectEdges(graph) {
        if (!graph) return [];
        if (typeof graph.edges === 'function') return toArray(graph.edges());
        if (graph.edges) return toArray(graph.edges);
        if (graph.elements?.edges) return toArray(graph.elements.edges);
        if (Array.isArray(graph.elements)) return graph.elements.filter(item => item.group === 'edges' || item.data?.source);
        return [];
    }

    function edgeSource(edge) {
        return String(callOrRead(edge, 'source') || '');
    }

    function edgeTarget(edge) {
        return String(callOrRead(edge, 'target') || '');
    }

    function nodeLabel(node) {
        return String(callOrRead(node, 'label') || callOrRead(node, 'name') || callOrRead(node, 'id') || '');
    }

    function regulationType(edge) {
        const raw = callOrRead(edge, 'regulation') || callOrRead(edge, 'type') || callOrRead(edge, 'role') || 'unknown';
        if (typeof raw === 'object') return 'unknown';
        const lower = String(raw).toLowerCase();
        if (lower === 'a' || lower.includes('activ')) return 'activation';
        if (lower === 'r' || lower.includes('repress') || lower.includes('inhibit')) return 'repression';
        return lower || 'unknown';
    }

    function edgeConfidence(edge) {
        return clamp(callOrRead(edge, 'confidence') ?? callOrRead(edge, 'confidenceScore') ?? 0);
    }

    function getPathwayList() {
        return adapter.getAllMetabolicPathways ? adapter.getAllMetabolicPathways() : [];
    }

    function findPathwayMetadata(pathwayIdOrName) {
        const query = normalizeGeneId(pathwayIdOrName);
        const match = getPathwayList().find(pathway => {
            const id = normalizeGeneId(pathway.pathwayId || pathway.id);
            const name = normalizeGeneId(pathway.pathwayName || pathway.name);
            return id === query || name === query || id.includes(query) || name.includes(query);
        });
        return {
            pathwayId: match?.pathwayId || match?.id || pathwayIdOrName || '',
            pathwayName: match?.pathwayName || match?.name || pathwayIdOrName || ''
        };
    }

    function getGenesForPathway(pathwayIdOrName) {
        if (!pathwayIdOrName) return [];
        const byGene = new Map();
        const genes = adapter.getGenesForPathway ? adapter.getGenesForPathway(pathwayIdOrName) : [];
        genes.forEach(gene => {
            const geneId = normalizeGeneId(gene.geneId || gene.locus || gene.id);
            if (!geneId) return;
            const entry = byGene.get(geneId) || {
                geneId,
                geneLabel: gene.geneLabel || gene.name || geneId,
                reactions: []
            };
            const seenReactions = new Set(entry.reactions.map(reaction => reaction.reactionId));
            (gene.reactions || []).forEach(reaction => {
                const reactionId = reaction.reactionId || reaction.id;
                if (!reactionId || seenReactions.has(reactionId)) return;
                seenReactions.add(reactionId);
                entry.reactions.push({
                    reactionId,
                    reactionName: reaction.reactionName || reaction.label || reactionId
                });
            });
            byGene.set(geneId, entry);
        });
        return Array.from(byGene.values()).sort((a, b) => a.geneId.localeCompare(b.geneId));
    }

    function getRegulatorsForPathway(graph, pathwayIdOrName) {
        const genes = getGenesForPathway(pathwayIdOrName);
        const pathwayGenes = new Set(genes.map(gene => normalizeGeneId(gene.geneId)));
        if (pathwayGenes.size === 0) return [];

        const nodesById = new Map();
        collectNodes(graph).forEach(node => {
            const id = normalizeGeneId(callOrRead(node, 'id'));
            if (id) nodesById.set(id, node);
        });

        const grouped = new Map();
        collectEdges(graph).forEach(edge => {
            const target = normalizeGeneId(edgeTarget(edge));
            if (!pathwayGenes.has(target)) return;
            const source = normalizeGeneId(edgeSource(edge));
            if (!source) return;
            const bucket = grouped.get(source) || { genes: new Set(), types: new Set(), confidences: [] };
            bucket.genes.add(target);
            bucket.types.add(regulationType(edge));
            bucket.confidences.push(edgeConfidence(edge));
            grouped.set(source, bucket);
        });

        return Array.from(grouped.entries()).map(([tfId, bucket]) => {
            const averageConfidence = bucket.confidences.length > 0
                ? bucket.confidences.reduce((sum, value) => sum + value, 0) / bucket.confidences.length
                : 0;
            const regulationDiversityScore = clamp(bucket.types.size / 4);
            const regulatorScore = clamp(0.50 * clamp(bucket.genes.size / 20) + 0.30 * averageConfidence + 0.20 * regulationDiversityScore);
            const label = nodeLabel(nodesById.get(tfId)) || tfId;
            return {
                tfId,
                tfLabel: label,
                regulatedGenes: Array.from(bucket.genes).sort(),
                regulationTypes: Array.from(bucket.types).sort(),
                averageConfidence: Number(averageConfidence.toFixed(3)),
                regulatorScore: Number(regulatorScore.toFixed(3)),
                explanation: `${label} regulates ${bucket.genes.size} genes associated with this pathway.`
            };
        }).sort((a, b) => b.regulatorScore - a.regulatorScore || b.regulatedGenes.length - a.regulatedGenes.length || a.tfId.localeCompare(b.tfId));
    }

    function generatePathwayRegulatoryExplanation(summary) {
        if (!summary || summary.totalRegulators === 0) {
            return 'No upstream transcription factors were found for this pathway based on the current regulatory network.';
        }
        const top = (summary.regulators || []).slice(0, 3).map(regulator => regulator.tfLabel || regulator.tfId);
        return `The selected pathway contains ${summary.totalGenes} genes mapped to ${summary.totalReactions} reactions. It is potentially regulated by ${summary.totalRegulators} transcription factors. The top predicted regulators include ${top.join(', ')}.`;
    }

    function getPathwayRegulatorySummary(graph, pathwayIdOrName) {
        const metadata = findPathwayMetadata(pathwayIdOrName);
        const genes = getGenesForPathway(pathwayIdOrName);
        const regulators = getRegulatorsForPathway(graph, pathwayIdOrName);
        const reactionIds = new Set();
        genes.forEach(gene => (gene.reactions || []).forEach(reaction => {
            if (reaction.reactionId) reactionIds.add(reaction.reactionId);
        }));
        const summary = {
            pathwayId: metadata.pathwayId,
            pathwayName: metadata.pathwayName,
            totalGenes: genes.length,
            totalReactions: reactionIds.size,
            totalRegulators: regulators.length,
            genes,
            regulators,
            explanation: ''
        };
        summary.explanation = generatePathwayRegulatoryExplanation(summary);
        return summary;
    }

    async function getPathwayRegulatorySummaryAsync(graph, pathwayIdOrName) {
        if (adapter.loadMetabolicPathways) {
            await adapter.loadMetabolicPathways();
        }
        return getPathwayRegulatorySummary(graph, pathwayIdOrName);
    }

    async function loadPathwayOptions() {
        const payload = adapter.loadMetabolicPathways ? await adapter.loadMetabolicPathways() : { pathways: [] };
        return payload.pathways || [];
    }

    global.pathwayRegulatoryView = {
        getGenesForPathway,
        getRegulatorsForPathway,
        generatePathwayRegulatoryExplanation,
        getPathwayRegulatorySummary,
        getPathwayRegulatorySummaryAsync,
        loadPathwayOptions
    };
})(window);
