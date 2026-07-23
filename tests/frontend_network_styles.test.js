'use strict';

const assert = require('node:assert/strict');
const versions = require('../web/method_versions.json');
const networkStyles = require('../web/lib/networkStyles.js');

assert.equal(networkStyles.STYLE_VERSION, versions.network_styles);
const baseNodes = networkStyles.createBaseNodeStyles();
assert.deepEqual(baseNodes.map(rule => rule.selector), [
    'node', 'node[type="TF"]', 'node[type="sRNA"]', 'node[type="query"]', 'node.shared-target',
]);
assert.equal(baseNodes[0].style.label, 'data(name)');
assert.equal(baseNodes[1].style['border-color'], '#1976d2');
assert.equal(baseNodes[2].style.shape, 'hexagon');
assert.equal(baseNodes[3].style.width, '34px');
assert.equal(baseNodes[4].style['border-width'], '2.5px');

const [baseEdge] = networkStyles.createBaseEdgeStyles();
const confidenceEdge = value => ({ data: key => key === 'confidenceScore' ? value : undefined });
assert.equal(baseEdge.selector, 'edge');
assert.equal(baseEdge.style.width(confidenceEdge(0.5)), 2.8);
assert(Math.abs(baseEdge.style.opacity(confidenceEdge(0.5)) - 0.65) < 1e-12);
assert.equal(baseEdge.style.width(confidenceEdge(0)), 2.0, 'zero retains the historical default');

const regulationEdges = networkStyles.createRegulationEdgeStyles();
assert.deepEqual(regulationEdges.map(rule => rule.selector), [
    'edge[regulationType="ppi"]',
    'edge[regulationType="activation"]',
    'edge[role="A"]',
    'edge[regulationType="repression"]',
    'edge[role="R"]',
    'edge[regulationType="dual"], edge[regulationType="sigma"], edge[regulationType="unknown"]',
    'edge[role="Dual"]',
    'edge[regulationType="post_transcriptional_repression"]',
    'edge[role="sRNA"]',
    'edge.confidence-high',
    'edge.confidence-medium',
    'edge.confidence-low',
]);
const ppiEdge = value => ({ data: key => key === 'score' ? value : undefined });
assert.equal(regulationEdges[0].style.width(ppiEdge(700)), 2);
assert.equal(regulationEdges[0].style.width(ppiEdge(1000)), 4);
assert.equal(regulationEdges[3].style['target-arrow-shape'], 'tee');
assert.equal(regulationEdges[8].style['line-style'], 'dashed');
assert.equal(regulationEdges[11].style.opacity, 0.42);

const interactionStates = networkStyles.createInteractionStateStyles();
assert.deepEqual(interactionStates.map(rule => rule.selector), [
    '.dimmed', '.rnaseq-hidden', 'node.highlighted', 'edge.highlighted',
    'node.sim-up', 'node.sim-down', 'node.sim-dual',
]);
assert.equal(interactionStates[0].style.opacity, 0.15);
assert.equal(interactionStates[1].style.display, 'none');
assert.equal(interactionStates[4].style['shadow-color'], '#2e7d32');
assert.equal(interactionStates[5].style['shadow-color'], '#d32f2f');
assert.equal(interactionStates[6].style['shadow-color'], '#e65100');
assert.equal(networkStyles.numericThreshold('0.01', 0.05), 0.01);
assert.equal(networkStyles.numericThreshold('invalid', 0.05), 0.05);
assert.equal(networkStyles.baseNodeSize('query'), 34);
assert.equal(networkStyles.baseNodeSize('TF'), 26);
assert.equal(networkStyles.baseNodeSize('sRNA'), 26);
assert.equal(networkStyles.baseNodeSize('Target'), 22);
assert.equal(networkStyles.expressionNodeSize(undefined, 'Target'), 22);
assert.equal(networkStyles.expressionNodeSize(2, 'TF'), 34);
assert.equal(networkStyles.expressionNodeSize(100, 'Target'), 38, 'expression sizing must be capped');

const thresholds = {
    'rnaseq-p-threshold': '0.05',
    'rnaseq-lfc-threshold': '1.0',
};
const [rule] = networkStyles.createRnaSeqStyles({
    colorForLog2FoldChange: value => value > 0 ? 'up' : 'down',
    thresholdValue: (id, fallback) => thresholds[id] ?? fallback,
});
const node = values => ({ data: key => values[key] });
const significantUp = node({ rnaseq_log2fc: 2, rnaseq_pvalue: 0.01, type: 'query' });
const insignificant = node({ rnaseq_log2fc: -0.5, rnaseq_pvalue: 0.2, type: 'Target' });

assert.equal(rule.selector, 'node.rnaseq-node');
assert.equal(rule.style['background-color'](significantUp), 'up');
assert.equal(rule.style['border-width'](significantUp), '3.5px');
assert.equal(rule.style['border-color'](significantUp), '#0f172a');
assert.equal(rule.style.width(significantUp), 42);
assert.equal(rule.style.height(significantUp), 42);
assert.equal(rule.style['shadow-blur'](significantUp), '12px');
assert.equal(rule.style['shadow-color'](significantUp), '#ef4444');
assert.equal(rule.style['border-width'](insignificant), '2px');
assert.equal(rule.style['shadow-blur'](insignificant), '0px');
assert.equal(rule.style['shadow-color'](insignificant), '#2563eb');

console.log('frontend network styles tests passed');
