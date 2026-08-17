'use strict';

const assert = require('node:assert/strict');
const versions = require('../web/method_versions.json');
const binder = require('../web/lib/networkInteractionBinder.js');

assert.equal(binder.BINDER_VERSION, versions.network_interaction_binder);

const classes = [];
const allNodes = {
    length: 61,
    addClass: name => classes.push(`add:${name}`),
    removeClass: name => classes.push(`remove:${name}`),
};
const shared = { indegree: () => 2, addClass: name => classes.push(`shared:${name}`) };
const unshared = { indegree: () => 1, addClass: () => assert.fail('single-parent target marked shared') };
const handlers = {};
let zoom = 0.4;
const graph = {
    nodes: selector => selector ? [shared, unshared] : allNodes,
    zoom: () => zoom,
    on: (event, selectorOrHandler, maybeHandler) => {
        const key = maybeHandler ? `${event}:${selectorOrHandler}` : `${event}:canvas`;
        handlers[key] = maybeHandler || selectorOrHandler;
    },
};

const updateLabels = binder.bindLevelOfDetail(graph);
updateLabels();
assert(classes.includes('add:lod-hide-label'));
zoom = 0.8;
handlers['zoom:canvas']();
assert(classes.includes('remove:lod-hide-label'));
assert.equal(binder.markSharedTargets(graph), 1);
assert(classes.includes('shared:shared-target'));

const calls = [];
const times = [1000, 1200];
binder.bindInteractions(graph, {
    highlightSubnet: node => calls.push(`highlight:${node.id()}`),
    showNodeDetails: locus => calls.push(`details:${locus}`),
    querySingleGene: locus => calls.push(`query:${locus}`),
    toggleRightSidebar: open => calls.push(`sidebar:${open}`),
    now: () => times.shift(),
});
const node = { id: () => 'cg0001' };
handlers['tap:node']({ target: node });
handlers['tap:node']({ target: node });
handlers['tap:canvas']({ target: graph });
assert.deepEqual(calls, [
    'highlight:cg0001', 'details:cg0001', 'sidebar:true',
    'highlight:cg0001', 'details:cg0001', 'sidebar:true', 'query:cg0001',
    'sidebar:false',
]);

console.log('frontend network interaction binder tests passed');
