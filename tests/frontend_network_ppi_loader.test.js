'use strict';

const assert = require('node:assert/strict');
const versions = require('../web/method_versions.json');
const ppiLoader = require('../web/lib/networkPpiLoader.js');

async function run() {
    assert.equal(ppiLoader.LOADER_VERSION, versions.network_ppi_loader);
    assert.deepEqual(ppiLoader.normalizeQuery([' cg1 ', null, 'cg2']), ['cg1', 'cg2']);
    assert.equal(ppiLoader.partnerUrl('cg 1'), '/api/analysis/string_ppi?gene=cg%201');
    assert.equal(
        ppiLoader.visibleEdgesUrl(['cg1', 'cg2']),
        '/api/analysis/network_ppi_edges?genes=cg1%2Ccg2'
    );

    const calls = [];
    const warnings = [];
    const client = {
        getJson: async (url, options) => {
            calls.push({ url, options });
            if (url.includes('bad')) throw new Error('unavailable');
            const locus = url.includes('cg1') ? 'cg1' : 'cg2';
            return { partners: [{ partner: `${locus}-partner`, score: 900, type: 'experimental' }] };
        },
    };
    const signal = new AbortController().signal;
    const loading = ppiLoader.loadQueryInteractions({
        query: ['cg1', 'bad', 'cg2'],
        enabled: true,
        client,
        signal,
        onWarning: (locus, error) => warnings.push([locus, error.message]),
    });
    assert.equal(calls.length, 3, 'all batch requests should start before the first one settles');
    const interactions = await loading;
    assert.deepEqual(interactions.map(item => item.source), ['cg1', 'cg2']);
    assert.deepEqual(warnings, [['bad', 'unavailable']]);
    assert(calls.every(call => call.options.signal === signal));

    const disabled = await ppiLoader.loadQueryInteractions({
        query: ['cg1'], enabled: false, client, signal,
    });
    assert.deepEqual(disabled, []);
    assert.equal(calls.length, 3);

    const abortController = new AbortController();
    abortController.abort();
    await assert.rejects(() => ppiLoader.loadQueryInteractions({
        query: ['cg1'],
        enabled: true,
        signal: abortController.signal,
        client: { getJson: async () => { throw new Error('aborted'); } },
    }), /aborted/);

    let visibleCalls = 0;
    const oneNodeGraph = { nodes: () => [{ id: () => 'cg1' }] };
    assert.deepEqual(await ppiLoader.loadVisibleEdges({
        graph: oneNodeGraph,
        enabled: true,
        client: { getJson: async () => { visibleCalls += 1; } },
        signal,
    }), []);
    assert.equal(visibleCalls, 0);

    const twoNodeGraph = { nodes: () => [{ id: () => 'cg1' }, { id: () => 'cg2' }] };
    const visible = await ppiLoader.loadVisibleEdges({
        graph: twoNodeGraph,
        enabled: true,
        client: {
            getJson: async (url, options) => {
                visibleCalls += 1;
                assert.equal(url, '/api/analysis/network_ppi_edges?genes=cg1%2Ccg2');
                assert.equal(options.signal, signal);
                return { edges: [{ source: 'cg1', target: 'cg2' }] };
            },
        },
        signal,
    });
    assert.equal(visibleCalls, 1);
    assert.deepEqual(visible, [{ source: 'cg1', target: 'cg2' }]);

    console.log('frontend network PPI loader tests passed');
}

run().catch(error => {
    console.error(error);
    process.exitCode = 1;
});
