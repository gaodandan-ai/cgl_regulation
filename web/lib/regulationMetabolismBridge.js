(function (global) {
    const adapter = global.metabolicModelAdapter || {};
    const normalizeGeneId = adapter.normalizeGeneId || function (id) {
        return String(id || '').trim().replace(/^gene[:_]/i, '').replace(/^G_/i, '').toLowerCase();
    };

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
        if (Array.isArray(graph.elements)) {
            return graph.elements.filter(item => item.group === 'nodes' || item.data?.id);
        }
        return [];
    }

    function collectEdges(graph) {
        if (!graph) return [];
        if (typeof graph.edges === 'function') return toArray(graph.edges());
        if (graph.edges) return toArray(graph.edges);
        if (graph.elements?.edges) return toArray(graph.elements.edges);
        if (Array.isArray(graph.elements)) {
            return graph.elements.filter(item => item.group === 'edges' || item.data?.source);
        }
        return [];
    }

    function nodeId(node) {
        return String(callOrRead(node, 'id') || '');
    }

    function edgeSource(edge) {
        return String(callOrRead(edge, 'source') || '');
    }

    function edgeTarget(edge) {
        return String(callOrRead(edge, 'target') || '');
    }

    function isGeneLike(node, id) {
        const type = String(callOrRead(node, 'type') || callOrRead(node, 'nodeType') || '').toLowerCase();
        if (['gene', 'target', 'mrna', 'transcript', 'orf'].includes(type)) return true;
        return /^(cg|cgl)\w+/i.test(id || '');
    }

    function regulationType(edge) {
        const raw = callOrRead(edge, 'regulation') || callOrRead(edge, 'type') || callOrRead(edge, 'role') || 'unknown';
        if (typeof raw === 'object') return 'unknown';
        const normalized = String(raw).toLowerCase();
        if (normalized === 'a' || normalized.includes('activ')) return 'activation';
        if (normalized === 'r' || normalized.includes('repress') || normalized.includes('inhibit')) return 'repression';
        return normalized || 'unknown';
    }

    function edgeConfidence(edge) {
        const raw = callOrRead(edge, 'confidence') ?? callOrRead(edge, 'confidenceScore');
        const parsed = Number(raw);
        return Number.isFinite(parsed) ? parsed : undefined;
    }

    function getTargetGenesForTF(graph, tfId) {
        const normalizedTf = normalizeGeneId(tfId);
        if (!normalizedTf) return [];

        const nodesById = new Map();
        collectNodes(graph).forEach(node => {
            const id = nodeId(node);
            if (id) nodesById.set(normalizeGeneId(id), node);
        });

        const targets = new Map();
        collectEdges(graph).forEach(edge => {
            const source = normalizeGeneId(edgeSource(edge));
            const rawTarget = edgeTarget(edge);
            const target = normalizeGeneId(rawTarget);
            if (!target || source !== normalizedTf) return;

            const targetNode = nodesById.get(target);
            if (targetNode && !isGeneLike(targetNode, rawTarget)) return;
            if (!targetNode && !isGeneLike(null, rawTarget)) return;

            if (!targets.has(target)) {
                targets.set(target, {
                    geneId: target,
                    node: targetNode,
                    confidence: edgeConfidence(edge),
                    regulation: regulationType(edge)
                });
            }
        });

        return Array.from(targets.values());
    }

    function getMetabolicContextForGene(geneId) {
        const normalizedGene = normalizeGeneId(geneId);
        if (!normalizedGene) return { geneId: '', reactions: [], pathways: [] };
        return {
            geneId: normalizedGene,
            reactions: adapter.getReactionsForGene ? adapter.getReactionsForGene(normalizedGene) : [],
            pathways: adapter.getPathwaysForGene ? adapter.getPathwaysForGene(normalizedGene) : []
        };
    }

    function getMetabolicImpactForTF(graph, tfId) {
        const targets = getTargetGenesForTF(graph, tfId);
        const mappedGenes = new Set();
        const reactionIds = new Set();
        const pathwayStats = new Map();

        targets.forEach(target => {
            const reactions = adapter.getReactionsForGene ? adapter.getReactionsForGene(target.geneId) : [];
            if (reactions.length > 0) mappedGenes.add(target.geneId);

            reactions.forEach(reaction => {
                const reactionKey = `${reaction.model || 'model'}:${reaction.id || reaction.label || 'reaction'}`;
                reactionIds.add(reactionKey);

                const pathwayId = reaction.pathway_id || reaction.pathway_name || 'Unassigned pathway';
                const pathwayName = reaction.pathway_name || pathwayId;
                const stat = pathwayStats.get(pathwayId) || {
                    name: pathwayName,
                    genes: new Set(),
                    reactions: new Set()
                };
                stat.genes.add(target.geneId);
                stat.reactions.add(reaction.id || reaction.label || reactionKey);
                pathwayStats.set(pathwayId, stat);
            });
        });

        const pathwaySummary = Array.from(pathwayStats.entries())
            .map(([pathwayId, stat]) => ({
                pathwayId,
                pathwayName: stat.name,
                geneCount: stat.genes.size,
                reactionCount: stat.reactions.size,
                genes: Array.from(stat.genes).sort(),
                reactions: Array.from(stat.reactions).sort()
            }))
            .sort((a, b) =>
                b.geneCount - a.geneCount ||
                b.reactionCount - a.reactionCount ||
                a.pathwayName.localeCompare(b.pathwayName)
            );

        return {
            tfId: normalizeGeneId(tfId),
            totalTargetGenes: targets.length,
            mappedTargetGenes: mappedGenes.size,
            totalReactions: reactionIds.size,
            totalPathways: pathwaySummary.length,
            pathwaySummary
        };
    }

    function generateMetabolicImpactExplanation(impact) {
        if (!impact || !impact.mappedTargetGenes || !impact.totalReactions) {
            return 'No metabolic model mapping available for this node.';
        }

        const topPathways = (impact.pathwaySummary || [])
            .slice(0, 3)
            .map(pathway => pathway.pathwayName)
            .filter(Boolean);
        const pathwayText = topPathways.length > 0
            ? ` The most affected pathways include ${topPathways.join(', ')}.`
            : ' No pathway annotations are available for the mapped reactions.';

        return `This transcription factor regulates ${impact.totalTargetGenes || 0} target genes, ${impact.mappedTargetGenes || 0} of which are mapped to metabolic reactions.${pathwayText}`;
    }

    global.regulationMetabolismBridge = {
        getTargetGenesForTF,
        getMetabolicContextForGene,
        getMetabolicImpactForTF,
        generateMetabolicImpactExplanation
    };
})(window);
