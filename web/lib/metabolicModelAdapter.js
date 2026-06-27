(function (global) {
    const impactCache = new Map();
    let pathwayMappingCache = null;

    function normalizeGeneId(id) {
        return String(id || '')
            .trim()
            .replace(/^gene[:_]/i, '')
            .replace(/^G_/i, '')
            .toLowerCase();
    }

    function cacheImpact(nodeId, payload) {
        const normalized = normalizeGeneId(nodeId || payload?.query);
        if (normalized) impactCache.set(normalized, payload);

        (payload?.affected_genes || []).forEach(gene => {
            const geneId = normalizeGeneId(gene.locus);
            if (!geneId) return;
            impactCache.set(geneId, {
                ...payload,
                mode: 'gene',
                affected_genes: [gene]
            });
        });
    }

    function getReactionsForGene(geneId) {
        const normalized = normalizeGeneId(geneId);
        const cached = impactCache.get(normalized);
        const gene = (cached?.affected_genes || []).find(item =>
            normalizeGeneId(item.locus) === normalized
        );
        return gene?.reactions || [];
    }

    function getPathwaysForGene(geneId) {
        const pathways = new Map();
        getReactionsForGene(geneId).forEach(reaction => {
            const id = reaction.pathway_id || reaction.pathway_name || 'Unassigned pathway';
            const pathway = pathways.get(id) || {
                id,
                name: reaction.pathway_name || id,
                model: reaction.model,
                gene_count: 1,
                reaction_count: 0,
                genes: [normalizeGeneId(geneId)],
                reactions: []
            };
            if (!pathway.reactions.includes(reaction.id)) {
                pathway.reactions.push(reaction.id);
            }
            pathway.reaction_count = pathway.reactions.length;
            pathways.set(id, pathway);
        });
        return Array.from(pathways.values());
    }

    function getEnzymeConstrainedReactionsForGene(geneId) {
        return getReactionsForGene(geneId).filter(reaction =>
            Boolean(reaction.enzyme_constraint || reaction.kcat || reaction.molecular_weight || reaction.kcat_MW)
        );
    }

    function getReactionVariantsForGene(geneId) {
        const grouped = {};
        getEnzymeConstrainedReactionsForGene(geneId).forEach(reaction => {
            const baseId = reaction.variant_of || reaction.id;
            if (!grouped[baseId]) grouped[baseId] = [];
            grouped[baseId].push(reaction);
        });
        return grouped;
    }

    function getMetabolicImpactForTF(tfId, graph) {
        const normalizedTf = normalizeGeneId(tfId);
        const targets = new Set();
        const graphEdges = graph?.edges || [];

        graphEdges.forEach(edge => {
            const source = normalizeGeneId(edge.source || edge.data?.source);
            const target = normalizeGeneId(edge.target || edge.data?.target);
            const type = edge.type || edge.data?.type || 'regulates';
            if (source === normalizedTf && target && type === 'regulates') {
                targets.add(target);
            }
        });

        const mappedGenes = [];
        const reactionIds = new Set();
        const pathwayStats = new Map();

        targets.forEach(geneId => {
            const reactions = getReactionsForGene(geneId);
            if (reactions.length === 0) return;
            mappedGenes.push({ geneId, reactions });

            reactions.forEach(reaction => {
                reactionIds.add(`${reaction.model || 'model'}:${reaction.id}`);
                const pathwayId = reaction.pathway_id || reaction.pathway_name || 'Unassigned pathway';
                const pathway = pathwayStats.get(pathwayId) || {
                    id: pathwayId,
                    name: reaction.pathway_name || pathwayId,
                    model: reaction.model,
                    gene_count: 0,
                    reaction_count: 0,
                    genes: [],
                    reactions: []
                };
                if (!pathway.genes.includes(geneId)) pathway.genes.push(geneId);
                if (!pathway.reactions.includes(reaction.id)) pathway.reactions.push(reaction.id);
                pathway.gene_count = pathway.genes.length;
                pathway.reaction_count = pathway.reactions.length;
                pathwayStats.set(pathwayId, pathway);
            });
        });

        const pathways = Array.from(pathwayStats.values()).sort((a, b) =>
            (b.gene_count || 0) - (a.gene_count || 0) ||
            (b.reaction_count || 0) - (a.reaction_count || 0) ||
            a.name.localeCompare(b.name)
        );

        return {
            tfId: normalizedTf,
            targetGeneCount: targets.size,
            mappedGeneCount: mappedGenes.length,
            reactionCount: reactionIds.size,
            pathwayCount: pathways.length,
            pathways,
            mappedGenes
        };
    }

    async function loadMetabolicImpact(nodeId) {
        const response = await fetch(`/api/metabolic_impact?gene=${encodeURIComponent(nodeId)}`);
        if (!response.ok) {
            throw new Error(`Metabolic impact request failed: ${response.status}`);
        }
        const payload = await response.json();
        cacheImpact(nodeId, payload);
        return payload;
    }

    async function loadMetabolicPathways(query = '') {
        if (!query && pathwayMappingCache) return pathwayMappingCache;
        const suffix = query ? `?query=${encodeURIComponent(query)}` : '';
        const response = await fetch(`/api/metabolic_pathways${suffix}`);
        if (!response.ok) {
            throw new Error(`Metabolic pathway request failed: ${response.status}`);
        }
        const payload = await response.json();
        if (!query) pathwayMappingCache = payload;
        return payload;
    }

    function getAllMetabolicPathways() {
        return pathwayMappingCache?.pathways || [];
    }

    function getGenesForPathway(pathwayIdOrName) {
        const query = normalizeGeneId(pathwayIdOrName);
        if (!query) return [];
        return (pathwayMappingCache?.pathways || [])
            .filter(pathway => {
                const id = normalizeGeneId(pathway.pathwayId || pathway.id);
                const name = normalizeGeneId(pathway.pathwayName || pathway.name);
                return id === query || name === query || id.includes(query) || name.includes(query);
            })
            .flatMap(pathway => pathway.genes || []);
    }

    global.metabolicModelAdapter = {
        normalizeGeneId,
        getReactionsForGene,
        getPathwaysForGene,
        getEnzymeConstrainedReactionsForGene,
        getReactionVariantsForGene,
        getMetabolicImpactForTF,
        loadMetabolicImpact,
        loadMetabolicPathways,
        getAllMetabolicPathways,
        getGenesForPathway
    };
})(window);
