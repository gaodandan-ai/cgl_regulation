'use strict';

const assert = require('node:assert/strict');
const versions = require('../web/method_versions.json');
const { SCHEMA_VERSION, mergeNode, createNormalizer } = require('../web/lib/networkNormalizer.js');

function run() {
    assert.equal(SCHEMA_VERSION, versions.network_normalization);
    const normalizer = createNormalizer({
        resolveLabel: (id, label) => `${label || id} (${id})`,
        getRfPrediction: (source, target) => source === 'cg0001' && target === 'cg0002'
            ? { predictedConfidence: 0.91, confidenceRank: '3', sampleType: 'held-out' }
            : null,
    });

    const row = {
        TF_locusTag: ' cg0001 ', TF_name: 'TF One', TF_altLocusTag: 'Cgl0001',
        TG_locusTag: 'cg0002', TG_name: 'Gene Two', TG_altLocusTag: 'Cgl0002',
        Role: 'A', Binding_site: 'ATGCGTACGT', Evidence: 'experimental',
        Source: 'curated', PMID: '123', expression_correlation: 0.7, Operon: 'op1',
    };
    const rfEdge = normalizer.normalizeTfEdge(row, 4);
    assert.equal(rfEdge.id, 'edge_cg0001_cg0002_4');
    assert.equal(rfEdge.regulationType, 'activation');
    assert.equal(rfEdge.confidenceScore, 0.91);
    assert.equal(rfEdge.confidenceModel, 'random_forest');
    assert.equal(rfEdge.confidenceFactors.randomForest, 0.91);
    assert.equal(rfEdge.evidence.rfSampleType, 'held-out');
    assert.equal(rfEdge.schemaVersion, SCHEMA_VERSION);
    assert.equal(rfEdge.original, row);
    assert.equal(normalizer.normalizeTfEdge({ TF_locusTag: '', TG_locusTag: 'cg1' }), null);

    const heuristicEdge = normalizer.normalizeTfEdge({
        TF_locusTag: 'cg0003', TG_locusTag: 'cg0004', Role: 'R', Evidence: 'predicted',
    });
    assert.equal(heuristicEdge.confidenceModel, 'heuristic');
    assert.equal(heuristicEdge.predictedConfidence, null);
    assert.equal(heuristicEdge.regulationType, 'repression');

    const srnaEdge = normalizer.normalizeSrnaEdge({
        srna: 'ncRNA1', mrna: 'cg0002', rank: 1,
        copra_pvalue: 0.001, copra_fdr: 0.05, energy: -20,
    });
    assert.equal(srnaEdge.interactionClass, 'sRNA-mRNA');
    assert.equal(srnaEdge.regulationType, 'post_transcriptional_repression');
    assert.equal(srnaEdge.confidenceLevel, 'medium');

    const nodes = {};
    mergeNode(nodes, { id: 'cg1', label: 'target', type: 'Target', aliases: { a: 1 } });
    mergeNode(nodes, { id: 'CG1', label: 'tf', type: 'TF', aliases: { b: 2 } });
    assert.equal(nodes.cg1.type, 'TF');
    assert.deepEqual(nodes.cg1.aliases, { a: 1, b: 2 });

    const graph = normalizer.normalizeNetwork([row], [{ srna: 'ncRNA1', mrna: 'cg0002' }]);
    assert.equal(graph.schemaVersion, SCHEMA_VERSION);
    assert.equal(graph.edges.length, 2);
    assert.equal(Object.keys(graph.nodes).length, 3);
    assert.equal(graph.nodes.cg0001.label, 'TF One (cg0001)');
    assert.equal(graph.nodes.cg0002.aliases.operon, 'op1');

    console.log('frontend network normalizer tests passed');
}

run();
