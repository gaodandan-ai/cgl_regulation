(function attachEvidenceScoring(root, factory) {
    const moduleApi = factory();
    if (typeof module === 'object' && module.exports) module.exports = moduleApi;
    if (root) root.CglEvidenceScoring = moduleApi;
})(typeof window !== 'undefined' ? window : globalThis, function createEvidenceScoringModule() {
    'use strict';

    const RULESET_VERSION = 'heuristic-confidence-v1.0.0';
    const FACTOR_WEIGHTS = Object.freeze({
        motif: 0.25,
        chip: 0.30,
        expression: 0.20,
        database: 0.25,
    });
    const LEVEL_THRESHOLDS = Object.freeze({ high: 0.75, medium: 0.50 });

    function cleanString(value) {
        if (value === undefined || value === null || value === 'null') return '';
        return String(value).trim();
    }

    function parseConfidenceScore(value) {
        const parsed = parseFloat(value);
        if (Number.isNaN(parsed)) return null;
        return Math.max(0, Math.min(1, parsed));
    }

    function normalizeRegulationType(role, sourceType = 'TF-TG') {
        const cleanRole = cleanString(role).toUpperCase();
        if (sourceType === 'sRNA-mRNA') return 'post_transcriptional_repression';
        if (cleanRole === 'A') return 'activation';
        if (cleanRole === 'R') return 'repression';
        if (cleanRole === 'DUAL') return 'dual';
        if (cleanRole === 'SIGMA') return 'sigma';
        return 'unknown';
    }

    function confidenceFromEvidence(evidence) {
        const text = cleanString(evidence).toLowerCase();
        if (text.includes('experimental') && text.includes('predicted')) return 0.78;
        if (text.includes('experimental')) return 0.86;
        if (text.includes('curated') || text.includes('literature')) return 0.74;
        if (text.includes('predicted')) return 0.42;
        return 0.32;
    }

    function confidenceFromMotif(bindingSite) {
        const site = cleanString(bindingSite);
        if (!site) return 0;
        const sites = site.split(';').map(value => value.trim()).filter(Boolean);
        if (sites.length >= 2) return 0.78;
        const longest = sites.reduce((max, value) => Math.max(max, value.length), 0);
        return longest >= 10 ? 0.66 : 0.48;
    }

    function confidenceFromChip(row = {}) {
        const evidence = [row.Evidence, row.Source, row.Method, row.Assay]
            .map(cleanString).join(' ').toLowerCase();
        if (evidence.includes('chip-exo')) return 0.95;
        if (evidence.includes('chip-seq') || evidence.includes('chip_seq') || evidence.includes('chip')) return 0.90;
        return 0;
    }

    function confidenceFromExpression(row = {}, sourceType = 'TF-TG') {
        if (sourceType === 'sRNA-mRNA') {
            const p = parseFloat(row.copra_pvalue);
            const fdr = parseFloat(row.copra_fdr);
            const energy = parseFloat(row.energy);
            let score = 0.35;
            if (!Number.isNaN(p)) score += p <= 0.001 ? 0.25 : p <= 0.01 ? 0.18 : p <= 0.05 ? 0.10 : 0;
            if (!Number.isNaN(fdr)) score += fdr <= 0.05 ? 0.20 : fdr <= 0.25 ? 0.12 : 0;
            if (!Number.isNaN(energy)) score += energy <= -20 ? 0.15 : energy <= -12 ? 0.08 : 0;
            return Math.min(0.90, score);
        }
        const correlation = parseFloat(
            row.expression_correlation ?? row.Expression_correlation ?? row.correlation ?? row.Correlation
        );
        return Number.isNaN(correlation) ? 0 : Math.min(0.95, Math.abs(correlation));
    }

    function combineConfidenceScores(factors = {}) {
        let weighted = 0;
        let usedWeight = 0;
        for (const [key, weight] of Object.entries(FACTOR_WEIGHTS)) {
            const value = Number(factors[key]) || 0;
            if (value > 0) {
                weighted += value * weight;
                usedWeight += weight;
            }
        }
        if (usedWeight === 0) return 0.25;
        const normalized = weighted / usedWeight;
        const evidenceCount = Object.values(factors).filter(value => Number(value) > 0.1).length;
        const multiEvidenceBonus = evidenceCount >= 2 ? 0.06 : 0;
        return Math.max(0.05, Math.min(0.99, normalized + multiEvidenceBonus));
    }

    function confidenceLevel(score) {
        if (score >= LEVEL_THRESHOLDS.high) return 'high';
        if (score >= LEVEL_THRESHOLDS.medium) return 'medium';
        return 'low';
    }

    function roleLabelFromType(role, regulationType) {
        if (regulationType === 'activation' || role === 'A') return 'Activation (+)';
        if (regulationType === 'repression' || role === 'R') return 'Repression (-)';
        if (regulationType === 'post_transcriptional_repression' || role === 'sRNA') {
            return 'sRNA / post-transcriptional repression';
        }
        if (regulationType === 'sigma') return 'Sigma factor';
        if (regulationType === 'dual' || role === 'Dual') return 'Dual regulation';
        return 'Unknown / pending';
    }

    function describeRules() {
        return {
            version: RULESET_VERSION,
            factorWeights: { ...FACTOR_WEIGHTS },
            levelThresholds: { ...LEVEL_THRESHOLDS },
            multiEvidenceBonus: 0.06,
            multiEvidenceMinimumFactors: 2,
        };
    }

    return {
        RULESET_VERSION,
        FACTOR_WEIGHTS,
        LEVEL_THRESHOLDS,
        parseConfidenceScore,
        normalizeRegulationType,
        confidenceFromEvidence,
        confidenceFromMotif,
        confidenceFromChip,
        confidenceFromExpression,
        combineConfidenceScores,
        confidenceLevel,
        roleLabelFromType,
        describeRules,
    };
});
