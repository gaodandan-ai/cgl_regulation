'use strict';

const assert = require('node:assert/strict');
const versions = require('../web/method_versions.json');
const navigation = require('../web/lib/queryNavigation.js');

function run() {
    assert.equal(navigation.STATE_VERSION, versions.query_navigation);
    assert.deepEqual(navigation.normalizeQueryList(null), []);
    assert.deepEqual(navigation.normalizeQueryList('cg1'), ['cg1']);
    assert.deepEqual(navigation.splitGeneQuery([' cg1,cg2 ', 'cg3']), ['cg1', 'cg2', 'cg3']);
    assert.equal(navigation.serializeGeneQuery([' cg1 ', '', 'cg2']), 'cg1,cg2');
    assert.equal(navigation.queriesEqual(['CG2', 'cg1'], ['cg1', 'cg2']), true);
    assert.equal(navigation.queriesEqual(['cg1'], ['cg1', 'cg1']), false);

    const parsed = navigation.parseUrlState('?workflow=gene&gene=cg1%2Ccg2');
    assert.equal(parsed.workflow, 'gene');
    assert.equal(parsed.gene, 'cg1,cg2');
    assert.deepEqual(parsed.genes, ['cg1', 'cg2']);
    assert.equal(navigation.parseUrlState('').workflow, 'gene');

    const built = navigation.buildUrlState(
        'https://example.test/app?old=1',
        { workflow: 'gene', gene: ['cg1', 'cg2'] }
    );
    assert.equal(built.href, 'https://example.test/app?old=1&workflow=gene&gene=cg1%2Ccg2');
    assert.deepEqual(built.state, {
        workflow: 'gene', gene: 'cg1,cg2', version: navigation.STATE_VERSION,
    });
    const withoutGene = navigation.buildUrlState(
        'https://example.test/app?workflow=gene&gene=cg1',
        { workflow: 'pathway', gene: null }
    );
    assert.equal(withoutGene.href, 'https://example.test/app?workflow=pathway');

    const history = navigation.createHistory();
    assert.equal(history.record(null, ['cg1']), false);
    assert.equal(history.record(['cg1'], ['CG1']), false);
    assert.equal(history.record(['cg1'], ['cg2']), true);
    assert.deepEqual(history.snapshot(), {
        canBack: true, canForward: false, backDepth: 1, forwardDepth: 0, suspended: false,
    });
    assert.deepEqual(history.go('back', ['cg2']), ['cg1']);
    assert.equal(history.snapshot().canForward, true);
    assert.deepEqual(history.go('forward', ['cg1']), ['cg2']);
    history.suspend();
    assert.equal(history.record(['cg2'], ['cg3']), false);
    assert.equal(history.snapshot().suspended, true);
    history.resume();
    assert.equal(history.record(['cg2'], ['cg3']), true);
    assert.equal(history.snapshot().canForward, false);

    console.log('frontend query navigation tests passed');
}

run();
