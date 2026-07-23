'use strict';

const assert = require('node:assert/strict');
const versions = require('../web/method_versions.json');
const networkGraph = require('../web/lib/networkGraph.js');

assert.equal(networkGraph.GRAPH_VERSION, versions.network_graph);
assert.equal(networkGraph.elementCount(null), 0);
assert.equal(networkGraph.elementCount([1, 2]), 2);
assert.equal(networkGraph.elementCount({ nodes: [1, 2], edges: [3, 4, 5] }), 5);

const smallOptions = networkGraph.createOptions({
    container: 'canvas',
    elements: { nodes: new Array(100), edges: new Array(150) },
    styles: ['style'],
    layoutName: 'cose',
});
assert.equal(smallOptions.hideEdgesOnViewport, false);
assert.equal(smallOptions.layout.name, 'cose');
assert.equal(smallOptions.style[0], 'style');

const largeOptions = networkGraph.createOptions({
    container: 'canvas',
    elements: { nodes: new Array(151), edges: new Array(100) },
    styles: [],
    layoutName: 'grid',
});
assert.equal(largeOptions.hideEdgesOnViewport, true, '251 total elements enable viewport optimization');

let receivedOptions = null;
const created = networkGraph.createGraph({
    cytoscapeImpl: options => { receivedOptions = options; return { graph: true }; },
    container: 'canvas',
    elements: [],
    styles: [],
    layoutName: 'circle',
});
assert.deepEqual(created, { graph: true });
assert.equal(receivedOptions.layout.name, 'circle');
assert.throws(() => networkGraph.createGraph({ cytoscapeImpl: null }), TypeError);

const existingIds = new Set(['ppi-cross-existing-other']);
let added = [];
const graph = {
    getElementById: id => ({ length: existingIds.has(id) ? 1 : 0 }),
    add: edges => { added = edges; edges.forEach(edge => existingIds.add(edge.data.id)); },
};
const addedCount = networkGraph.addPpiEdges(graph, [
    { source: 'cg1', target: 'cg2', score: 900 },
    { source: 'cg2', target: 'cg1', score: 900 },
    { source: 'other', target: 'existing', score: 800 },
    { source: '', target: 'cg3', score: 700 },
]);
assert.equal(addedCount, 1);
assert.equal(added.length, 1);
assert.deepEqual(added[0], {
    group: 'edges',
    data: {
        id: 'ppi-cross-cg1-cg2',
        source: 'cg1',
        target: 'cg2',
        role: 'protein-protein interaction',
        type: 'PPI',
        regulationType: 'ppi',
        score: 900,
        schemaVersion: 'unified-v1',
    },
});

console.log('frontend network graph tests passed');
