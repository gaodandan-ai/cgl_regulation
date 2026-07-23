'use strict';

const assert = require('node:assert/strict');
const versions = require('../web/method_versions.json');
const scoring = require('../web/lib/evidenceScoring.js');

function run() {
    assert.equal(scoring.RULESET_VERSION, versions.evidence_scoring);
    assert.deepEqual(scoring.describeRules().factorWeights, {
        motif: 0.25, chip: 0.30, expression: 0.20, database: 0.25,
    });

    assert.equal(scoring.parseConfidenceScore('bad'), null);
    assert.equal(scoring.parseConfidenceScore(-1), 0);
    assert.equal(scoring.parseConfidenceScore(2), 1);
    assert.equal(scoring.normalizeRegulationType(' a '), 'activation');
    assert.equal(scoring.normalizeRegulationType('R'), 'repression');
    assert.equal(scoring.normalizeRegulationType('anything', 'sRNA-mRNA'), 'post_transcriptional_repression');

    assert.equal(scoring.confidenceFromEvidence('experimental and predicted'), 0.78);
    assert.equal(scoring.confidenceFromEvidence('experimental'), 0.86);
    assert.equal(scoring.confidenceFromEvidence('curated literature'), 0.74);
    assert.equal(scoring.confidenceFromEvidence('predicted'), 0.42);
    assert.equal(scoring.confidenceFromEvidence(''), 0.32);

    assert.equal(scoring.confidenceFromMotif(''), 0);
    assert.equal(scoring.confidenceFromMotif('ATGC'), 0.48);
    assert.equal(scoring.confidenceFromMotif('ATGCGTACGT'), 0.66);
    assert.equal(scoring.confidenceFromMotif('ATGC; CGTA'), 0.78);
    assert.equal(scoring.confidenceFromChip({ Assay: 'ChIP-exo' }), 0.95);
    assert.equal(scoring.confidenceFromChip({ Method: 'ChIP-seq' }), 0.90);
    assert.equal(scoring.confidenceFromChip({}), 0);

    assert.equal(scoring.confidenceFromExpression({ correlation: -0.8 }), 0.8);
    assert.equal(scoring.confidenceFromExpression({}, 'sRNA-mRNA'), 0.35);
    assert.equal(scoring.confidenceFromExpression({
        copra_pvalue: 0.001, copra_fdr: 0.05, energy: -20,
    }, 'sRNA-mRNA'), 0.9);

    assert.equal(scoring.combineConfidenceScores({}), 0.25);
    assert.equal(scoring.combineConfidenceScores({ motif: 0.66 }), 0.66);
    assert.ok(Math.abs(
        scoring.combineConfidenceScores({ motif: 0.66, database: 0.86 }) - 0.82
    ) < 1e-12);
    assert.equal(scoring.confidenceLevel(0.75), 'high');
    assert.equal(scoring.confidenceLevel(0.50), 'medium');
    assert.equal(scoring.confidenceLevel(0.49), 'low');
    assert.equal(scoring.roleLabelFromType('sRNA', ''), 'sRNA / post-transcriptional repression');

    console.log('frontend evidence scoring tests passed');
}

run();
