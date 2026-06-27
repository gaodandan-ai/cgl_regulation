(function (global) {
    const adapter = global.metabolicModelAdapter || {};
    const ranking = global.tfMetabolicImpactRanking || {};
    const normalizeGeneId = adapter.normalizeGeneId || function (id) {
        return String(id || '').trim().replace(/^gene[:_]/i, '').replace(/^G_/i, '').toLowerCase();
    };

    const ENGINEERING_PATHWAY_KEYWORDS = [
        'glutamate',
        'glutamic acid',
        'amino acid',
        'lysine',
        'arginine',
        'tca',
        'citric acid cycle',
        'central carbon',
        'glycolysis',
        'pyruvate',
        'acetyl-coa',
        'transport',
        'nitrogen',
        'carbon metabolism'
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

    function edgeRegulation(edge) {
        const raw = callOrRead(edge, 'regulation') || callOrRead(edge, 'type') || callOrRead(edge, 'role') || 'unknown';
        if (typeof raw === 'object') return 'unknown';
        const lower = String(raw).toLowerCase();
        if (lower === 'a' || lower.includes('activ')) return 'activation';
        if (lower === 'r' || lower.includes('repress') || lower.includes('inhibit')) return 'repression';
        if (lower.includes('predict')) return 'predicted';
        return lower || 'unknown';
    }

    function calculateEngineeringCandidateScore(input) {
        const mappedGenes = clamp((Number(input?.mappedTargetGenes) || 0) / 30);
        const reactions = clamp((Number(input?.totalReactions) || 0) / 60);
        const pathways = clamp((Number(input?.totalPathways) || 0) / 15);
        const confidence = clamp(Number(input?.averageConfidence) || 0);
        const keyPathways = clamp((Number(input?.keyPathwayHits) || 0) / 8);
        const keyGenes = clamp((Number(input?.regulatedKeyGeneCount) || 0) / 20);
        const base = 0.25 * mappedGenes + 0.20 * reactions + 0.15 * pathways + 0.20 * confidence + 0.15 * keyPathways + 0.05 * keyGenes;
        return Number(clamp(base).toFixed(3));
    }

    function getRegulationProfile(edges) {
        const profile = { activationCount: 0, repressionCount: 0, predictedCount: 0, unknownCount: 0 };
        (edges || []).forEach(edge => {
            const type = edgeRegulation(edge);
            if (type === 'activation') profile.activationCount += 1;
            else if (type === 'repression') profile.repressionCount += 1;
            else if (type === 'predicted') profile.predictedCount += 1;
            else profile.unknownCount += 1;
        });
        return profile;
    }

    function recommendationLevel(score) {
        if (score >= 0.75) return 'high';
        if (score >= 0.45) return 'medium';
        return 'low';
    }

    function matchesEngineeringPathway(pathway, filter) {
        const text = String(pathway || '').toLowerCase();
        if (!text) return false;
        if (filter && !text.includes(String(filter).toLowerCase())) return false;
        return ENGINEERING_PATHWAY_KEYWORDS.some(keyword => text.includes(keyword));
    }

    function generateEngineeringTargetRationale(candidate) {
        if (!candidate || candidate.candidateScore < 0.45) {
            return 'This TF has limited evidence as an engineering target based on the current regulatory and metabolic mappings.';
        }
        const levelText = candidate.recommendationLevel === 'high' ? 'high-priority' : 'medium-priority';
        const pathwayText = candidate.keyPathways.length > 0
            ? ', including ' + candidate.keyPathways.slice(0, 3).join(' and ')
            : '';
        const keyGeneText = candidate.regulatedKeyGenes.length > 0
            ? ' Its targets include ' + candidate.regulatedKeyGenes.slice(0, 5).join(', ') + ', which are associated with key metabolic modules.'
            : '';
        return `${candidate.tfLabel || candidate.tfId} is a ${levelText} candidate engineering regulator. It regulates ${candidate.mappedTargetGenes} genes mapped to ${candidate.totalReactions} metabolic reactions across ${candidate.totalPathways} pathways${pathwayText}.${keyGeneText} This TF may influence metabolic phenotype and could be prioritized for further experimental or simulation-based evaluation.`;
    }

    function candidateFromRank(rank, graph) {
        const allEdges = collectEdges(graph);
        const tfEdges = allEdges.filter(edge => normalizeGeneId(edgeSource(edge)) === normalizeGeneId(rank.tfId));
        const regulationProfile = getRegulationProfile(tfEdges);
        const keyPathwayEntries = (rank.pathwaySummary || []).filter(pathway =>
            matchesEngineeringPathway(`${pathway.pathwayName || ''} ${pathway.pathwayId || ''}`)
        );
        const keyPathways = Array.from(new Set(keyPathwayEntries.map(pathway => pathway.pathwayName || pathway.pathwayId).filter(Boolean)));
        const regulatedKeyGenes = Array.from(new Set(keyPathwayEntries.flatMap(pathway => pathway.genes || []))).sort();
        const candidateScore = calculateEngineeringCandidateScore({
            mappedTargetGenes: rank.mappedTargetGenes,
            totalReactions: rank.totalReactions,
            totalPathways: rank.totalPathways,
            averageConfidence: rank.averageConfidence,
            keyPathwayHits: keyPathways.length,
            regulatedKeyGeneCount: regulatedKeyGenes.length,
            repressionCount: regulationProfile.repressionCount,
            activationCount: regulationProfile.activationCount
        });
        const candidate = {
            tfId: rank.tfId,
            tfLabel: rank.tfLabel || rank.tfId,
            candidateScore,
            totalTargetGenes: rank.totalTargetGenes || 0,
            mappedTargetGenes: rank.mappedTargetGenes || 0,
            totalReactions: rank.totalReactions || 0,
            totalPathways: rank.totalPathways || 0,
            keyPathways,
            regulatedKeyGenes,
            averageConfidence: rank.averageConfidence || 0,
            regulationProfile,
            recommendationLevel: recommendationLevel(candidateScore),
            rationale: ''
        };
        candidate.rationale = generateEngineeringTargetRationale(candidate);
        return candidate;
    }

    function applyCandidateFilters(candidates, options = {}) {
        const limit = options.limit ?? 20;
        const minCandidateScore = Number(options.minCandidateScore ?? 0);
        const pathwayFilter = String(options.pathwayKeywordFilter || '').trim().toLowerCase();
        return (candidates || [])
            .filter(candidate => candidate.candidateScore >= minCandidateScore)
            .filter(candidate => !options.recommendationLevel || candidate.recommendationLevel === options.recommendationLevel)
            .filter(candidate => !options.searchTf || `${candidate.tfId} ${candidate.tfLabel}`.toLowerCase().includes(String(options.searchTf).toLowerCase()))
            .filter(candidate => !pathwayFilter || candidate.keyPathways.some(pathway => pathway.toLowerCase().includes(pathwayFilter)))
            .filter(candidate => options.includeLowConfidence || candidate.averageConfidence > 0 || candidate.mappedTargetGenes > 0)
            .sort((a, b) => b.candidateScore - a.candidateScore || b.mappedTargetGenes - a.mappedTargetGenes || a.tfLabel.localeCompare(b.tfLabel))
            .slice(0, limit);
    }

    function findEngineeringTargetCandidates(graph, options = {}) {
        const ranks = ranking.rankTFsByMetabolicImpact
            ? ranking.rankTFsByMetabolicImpact(graph, { limit: 500, includeZeroImpact: true })
            : [];
        return applyCandidateFilters(ranks.map(rank => candidateFromRank(rank, graph)), options);
    }

    async function findEngineeringTargetCandidatesAsync(graph, options = {}) {
        const ranks = ranking.rankTFsByMetabolicImpactAsync
            ? await ranking.rankTFsByMetabolicImpactAsync(graph, { limit: 500, includeZeroImpact: true, batchSize: options.batchSize || 8 })
            : [];
        return applyCandidateFilters(ranks.map(rank => candidateFromRank(rank, graph)), options);
    }

    global.candidateEngineeringTargets = {
        ENGINEERING_PATHWAY_KEYWORDS,
        calculateEngineeringCandidateScore,
        getRegulationProfile,
        findEngineeringTargetCandidates,
        findEngineeringTargetCandidatesAsync,
        generateEngineeringTargetRationale
    };
})(window);
