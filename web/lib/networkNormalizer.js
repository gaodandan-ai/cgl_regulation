(function attachNetworkNormalizer(root, factory) {
    const scoring = root && root.CglEvidenceScoring
        ? root.CglEvidenceScoring
        : (typeof module === 'object' && module.exports ? require('./evidenceScoring.js') : null);
    const moduleApi = factory(scoring);
    if (typeof module === 'object' && module.exports) module.exports = moduleApi;
    if (root) root.CglNetworkNormalizer = moduleApi;
})(typeof window !== 'undefined' ? window : globalThis, function createNetworkNormalizerModule(scoring) {
    'use strict';

    const SCHEMA_VERSION = 'normalized-network-v1.0.0';
    const TYPE_RANK = Object.freeze({ query: 4, TF: 3, sRNA: 2, Target: 1 });

    function cleanString(value) {
        if (value === undefined || value === null || value === 'null') return '';
        return String(value).trim();
    }

    function mergeNode(nodesById, node) {
        if (!node) return null;
        const key = node.id.toLowerCase();
        const existing = nodesById[key];
        if (!existing) {
            nodesById[key] = node;
            return node;
        }
        const chosenType = (TYPE_RANK[node.type] || 0) > (TYPE_RANK[existing.type] || 0)
            ? node.type
            : existing.type;
        const merged = {
            ...existing,
            ...node,
            type: chosenType,
            aliases: { ...(existing.aliases || {}), ...(node.aliases || {}) },
        };
        nodesById[key] = merged;
        return merged;
    }

    function createNormalizer({ resolveLabel, getRfPrediction, scoringApi = scoring } = {}) {
        if (!scoringApi) throw new Error('CglEvidenceScoring is required');
        const labelResolver = typeof resolveLabel === 'function'
            ? resolveLabel
            : (_id, label) => label;
        const rfResolver = typeof getRfPrediction === 'function'
            ? getRfPrediction
            : () => null;

        function normalizeNodeRecord(id, label, type, aliases = {}) {
            const cleanId = cleanString(id);
            if (!cleanId) return null;
            return {
                id: cleanId,
                label: labelResolver(cleanId, label || cleanId),
                type,
                aliases,
                dataSource: 'local_csv',
                schemaVersion: SCHEMA_VERSION,
            };
        }

        function normalizeTfEdge(row = {}, index = 0) {
            const source = cleanString(row.TF_locusTag || row.TF_altLocusTag || row.TF_name);
            const target = cleanString(row.TG_locusTag || row.TG_altLocusTag || row.TG_name);
            if (!source || !target) return null;
            const regulationType = scoringApi.normalizeRegulationType(row.Role, 'TF-TG');
            const factors = {
                motif: scoringApi.confidenceFromMotif(row.Binding_site),
                chip: scoringApi.confidenceFromChip(row),
                expression: scoringApi.confidenceFromExpression(row, 'TF-TG'),
                database: scoringApi.confidenceFromEvidence(row.Evidence || row.Source),
            };
            const heuristicConfidenceScore = scoringApi.combineConfidenceScores(factors);
            const rfPrediction = rfResolver(source, target);
            const confidenceScore = rfPrediction?.predictedConfidence ?? heuristicConfidenceScore;
            if (rfPrediction) factors.randomForest = rfPrediction.predictedConfidence;
            const role = cleanString(row.Role);
            return {
                id: `edge_${source}_${target}_${index}`,
                source,
                target,
                sourceType: 'TF',
                targetType: 'Target',
                regulationType,
                role,
                legacyRole: role,
                interactionClass: 'TF-TG',
                schemaVersion: SCHEMA_VERSION,
                confidenceScore,
                heuristicConfidenceScore,
                predictedConfidence: rfPrediction?.predictedConfidence ?? null,
                confidenceModel: rfPrediction ? 'random_forest' : 'heuristic',
                rfConfidenceRank: rfPrediction?.confidenceRank || '',
                confidenceLevel: scoringApi.confidenceLevel(confidenceScore),
                confidenceFactors: factors,
                evidence: {
                    motifSequence: cleanString(row.Binding_site),
                    databaseEvidence: cleanString(row.Evidence),
                    source: cleanString(row.Source),
                    pmid: cleanString(row.PMID),
                    expressionCorrelation: cleanString(
                        row.expression_correlation ?? row.Expression_correlation ?? row.correlation ?? ''
                    ),
                    rfConfidenceRank: rfPrediction?.confidenceRank || '',
                    rfSampleType: rfPrediction?.sampleType || '',
                    rfLabel: rfPrediction?.label || '',
                    rfFeatureMissingCount: rfPrediction?.featureMissingCount || '',
                    rfExpressionFeatureAvailable: rfPrediction?.expressionFeatureAvailable || '',
                    rfTargetMappedReactionCount: rfPrediction?.targetMappedReactionCount || '',
                    rfTargetMappedPathwayCount: rfPrediction?.targetMappedPathwayCount || '',
                    rfTargetEnzymeConstrainedReactionCount: rfPrediction?.targetEnzymeConstrainedReactionCount || '',
                    rfTargetKcatMedian: rfPrediction?.targetKcatMedian || '',
                    rfTargetKcatMwMedian: rfPrediction?.targetKcatMwMedian || '',
                },
                original: row,
            };
        }

        function normalizeSrnaEdge(row = {}, index = 0) {
            const source = cleanString(row.srna);
            const target = cleanString(row.mrna);
            if (!source || !target) return null;
            const factors = {
                motif: 0,
                chip: 0,
                expression: scoringApi.confidenceFromExpression(row, 'sRNA-mRNA'),
                database: 0.45,
            };
            const confidenceScore = scoringApi.combineConfidenceScores(factors);
            return {
                id: `edge_srna_${source}_${target}_${index}`,
                source,
                target,
                sourceType: 'sRNA',
                targetType: 'Target',
                regulationType: 'post_transcriptional_repression',
                role: 'sRNA',
                legacyRole: 'sRNA',
                interactionClass: 'sRNA-mRNA',
                schemaVersion: SCHEMA_VERSION,
                confidenceScore,
                confidenceLevel: scoringApi.confidenceLevel(confidenceScore),
                confidenceFactors: factors,
                evidence: {
                    rank: row.rank,
                    energy: row.energy,
                    copraPvalue: row.copra_pvalue,
                    copraFdr: row.copra_fdr,
                    source: 'sRNA prediction',
                },
                original: row,
            };
        }

        function normalizeNetwork(tfRows = [], srnaRows = []) {
            const nodes = {};
            const edges = [];
            tfRows.forEach((row, index) => {
                const edge = normalizeTfEdge(row, index);
                if (!edge) return;
                edges.push(edge);
                mergeNode(nodes, normalizeNodeRecord(edge.source, cleanString(row.TF_name), 'TF', {
                    altLocus: cleanString(row.TF_altLocusTag),
                }));
                mergeNode(nodes, normalizeNodeRecord(edge.target, cleanString(row.TG_name), 'Target', {
                    altLocus: cleanString(row.TG_altLocusTag),
                    operon: cleanString(row.Operon),
                }));
            });
            srnaRows.forEach((row, index) => {
                const edge = normalizeSrnaEdge(row, index);
                if (!edge) return;
                edges.push(edge);
                mergeNode(nodes, normalizeNodeRecord(edge.source, edge.source, 'sRNA'));
                mergeNode(nodes, normalizeNodeRecord(edge.target, edge.target, 'Target'));
            });
            return { nodes, edges, schemaVersion: SCHEMA_VERSION };
        }

        return { normalizeNodeRecord, normalizeTfEdge, normalizeSrnaEdge, normalizeNetwork };
    }

    return { SCHEMA_VERSION, TYPE_RANK, mergeNode, createNormalizer };
});
