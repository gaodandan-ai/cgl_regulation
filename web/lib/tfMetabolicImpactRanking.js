(function (global) {
    const adapter = global.metabolicModelAdapter || {};
    const bridge = global.regulationMetabolismBridge || {};
    const normalizeGeneId = adapter.normalizeGeneId || function (id) {
        return String(id || '').trim().replace(/^gene[:_]/i, '').replace(/^G_/i, '').toLowerCase();
    };

    const KEY_CGL_PATHWAY_KEYWORDS = [
        'glutamate',
        'glutamic acid',
        'amino acid',
        'lysine',
        'arginine',
        'tca',
        'citric acid cycle',
        'central carbon',
        'glycolysis',
        'transport'
    ];

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

    function nodeId(node) {
        return String(callOrRead(node, 'id') || '');
    }

    function nodeLabel(node) {
        return String(callOrRead(node, 'label') || callOrRead(node, 'name') || nodeId(node));
    }

    function isTFNode(node) {
        const type = String(callOrRead(node, 'type') || '').toLowerCase();
        return ['tf', 'transcription_factor', 'regulator'].includes(type);
    }

    function edgeSource(edge) {
        return String(callOrRead(edge, 'source') || '');
    }

    function edgeConfidence(edge) {
        return clamp(callOrRead(edge, 'confidence') ?? callOrRead(edge, 'confidenceScore') ?? 0);
    }

    function getAllTFNodes(graph) {
        const seen = new Set();
        const tfs = [];
        collectNodes(graph).forEach(node => {
            const id = nodeId(node);
            const normalized = normalizeGeneId(id);
            if (!id || seen.has(normalized) || !isTFNode(node)) return;
            seen.add(normalized);
            tfs.push({
                id,
                label: nodeLabel(node),
                type: callOrRead(node, 'type'),
                data: node.data
            });
        });
        return tfs;
    }

    function calculateTFMetabolicImpactScore(input) {
        const mappedGenes = clamp((Number(input?.mappedTargetGenes) || 0) / 30);
        const reactions = clamp((Number(input?.totalReactions) || 0) / 60);
        const pathways = clamp((Number(input?.totalPathways) || 0) / 15);
        const confidence = clamp(Number(input?.averageConfidence) || 0);
        const keyHits = Math.max(0, Number(input?.keyPathwayHits) || 0);
        const base = 0.35 * mappedGenes + 0.25 * reactions + 0.20 * pathways + 0.20 * confidence;
        const bonus = 0.05 * keyHits;
        return Number(clamp(base + bonus).toFixed(3));
    }

    function detectKeyPathways(pathwayNames) {
        const found = new Map();
        (pathwayNames || []).forEach(name => {
            const lower = String(name || '').toLowerCase();
            if (!lower) return;
            if (KEY_CGL_PATHWAY_KEYWORDS.some(keyword => lower.includes(keyword))) {
                found.set(name, name);
            }
        });
        return Array.from(found.values());
    }

    function generateTFImpactRankExplanation(rank) {
        if (!rank || rank.impactScore < 0.15 || rank.mappedTargetGenes === 0) {
            return 'This TF has limited metabolic model coverage based on current mappings.';
        }
        const pathwayText = rank.keyPathways.length > 0
            ? ', including ' + rank.keyPathways.slice(0, 3).join(' and ')
            : '';
        const level = rank.impactScore >= 0.7 ? 'high' : rank.impactScore >= 0.4 ? 'moderate' : 'limited';
        return `${rank.tfLabel || rank.tfId} has a ${level} predicted metabolic impact. It regulates ${rank.mappedTargetGenes} genes mapped to ${rank.totalReactions} metabolic reactions across ${rank.totalPathways} pathways${pathwayText}.`;
    }

    function rankTFsByMetabolicImpact(graph, options = {}) {
        const limit = options.limit ?? 20;
        const includeZeroImpact = options.includeZeroImpact ?? false;
        const edges = collectEdges(graph);
        return getAllTFNodes(graph).map(tf => {
            const impact = bridge.getMetabolicImpactForTF
                ? bridge.getMetabolicImpactForTF(graph, tf.id)
                : { totalTargetGenes: 0, mappedTargetGenes: 0, totalReactions: 0, totalPathways: 0, pathwaySummary: [] };
            const targets = bridge.getTargetGenesForTF ? bridge.getTargetGenesForTF(graph, tf.id) : [];
            const outgoingConfidences = edges
                .filter(edge => normalizeGeneId(edgeSource(edge)) === normalizeGeneId(tf.id))
                .map(edgeConfidence);
            const averageConfidence = outgoingConfidences.length > 0
                ? outgoingConfidences.reduce((sum, value) => sum + value, 0) / outgoingConfidences.length
                : 0;
            const keyPathways = detectKeyPathways((impact.pathwaySummary || []).map(pathway => pathway.pathwayName));
            const impactScore = calculateTFMetabolicImpactScore({
                mappedTargetGenes: impact.mappedTargetGenes,
                totalReactions: impact.totalReactions,
                totalPathways: impact.totalPathways,
                averageConfidence,
                keyPathwayHits: keyPathways.length
            });
            const rank = {
                tfId: tf.id,
                tfLabel: tf.label || tf.id,
                totalTargetGenes: impact.totalTargetGenes || targets.length,
                mappedTargetGenes: impact.mappedTargetGenes || 0,
                totalReactions: impact.totalReactions || 0,
                totalPathways: impact.totalPathways || 0,
                averageConfidence: Number(averageConfidence.toFixed(3)),
                keyPathways,
                impactScore,
                explanation: ''
            };
            rank.explanation = generateTFImpactRankExplanation(rank);
            return rank;
        })
            .filter(rank => includeZeroImpact || rank.impactScore > 0)
            .sort((a, b) => b.impactScore - a.impactScore || b.mappedTargetGenes - a.mappedTargetGenes || a.tfLabel.localeCompare(b.tfLabel))
            .slice(0, limit);
    }

    function buildRankFromApiPayload(tf, payload, graph) {
        const summary = payload?.summary || {};
        const pathways = payload?.pathways || [];
        const edges = collectEdges(graph).filter(edge => normalizeGeneId(edgeSource(edge)) === normalizeGeneId(tf.id));
        const averageConfidence = edges.length > 0
            ? edges.map(edgeConfidence).reduce((sum, value) => sum + value, 0) / edges.length
            : 0;
        const pathwaySummary = pathways.map(pathway => ({
            pathwayId: pathway.id || pathway.name || 'Unassigned pathway',
            pathwayName: pathway.name || pathway.id || 'Unassigned pathway',
            geneCount: Number(pathway.gene_count || 0),
            reactionCount: Number(pathway.reaction_count || 0),
            genes: pathway.genes || [],
            reactions: pathway.reactions || []
        }));
        const keyPathways = detectKeyPathways(pathwaySummary.map(pathway => pathway.pathwayName));
        const rank = {
            tfId: tf.id,
            tfLabel: tf.label || tf.id,
            totalTargetGenes: Number(summary.target_gene_count || 0),
            mappedTargetGenes: Number(summary.mapped_gene_count || 0),
            totalReactions: Number(summary.reaction_count || 0),
            totalPathways: Number(summary.pathway_count || pathways.length || 0),
            averageConfidence: Number(averageConfidence.toFixed(3)),
            keyPathways,
            pathwaySummary,
            impactScore: 0,
            explanation: ''
        };
        rank.impactScore = calculateTFMetabolicImpactScore({
            mappedTargetGenes: rank.mappedTargetGenes,
            totalReactions: rank.totalReactions,
            totalPathways: rank.totalPathways,
            averageConfidence: rank.averageConfidence,
            keyPathwayHits: keyPathways.length
        });
        rank.explanation = generateTFImpactRankExplanation(rank);
        return rank;
    }

    async function rankTFsByMetabolicImpactAsync(graph, options = {}) {
        const limit = options.limit ?? 20;
        const includeZeroImpact = options.includeZeroImpact ?? false;
        const tfs = getAllTFNodes(graph);
        const ranks = [];
        const batchSize = options.batchSize || 8;

        for (let i = 0; i < tfs.length; i += batchSize) {
            const batch = tfs.slice(i, i + batchSize);
            const partial = await Promise.all(batch.map(async tf => {
                try {
                    const payload = adapter.loadMetabolicImpact
                        ? await adapter.loadMetabolicImpact(tf.id)
                        : null;
                    return buildRankFromApiPayload(tf, payload, graph);
                } catch (err) {
                    console.warn('Failed to rank TF metabolic impact:', tf.id, err);
                    return null;
                }
            }));
            partial.forEach(rank => {
                if (rank) ranks.push(rank);
            });
        }

        return ranks
            .filter(rank => includeZeroImpact || rank.impactScore > 0)
            .sort((a, b) => b.impactScore - a.impactScore || b.mappedTargetGenes - a.mappedTargetGenes || a.tfLabel.localeCompare(b.tfLabel))
            .slice(0, limit);
    }

    global.tfMetabolicImpactRanking = {
        KEY_CGL_PATHWAY_KEYWORDS,
        getAllTFNodes,
        calculateTFMetabolicImpactScore,
        rankTFsByMetabolicImpact,
        rankTFsByMetabolicImpactAsync,
        generateTFImpactRankExplanation
    };
})(window);
