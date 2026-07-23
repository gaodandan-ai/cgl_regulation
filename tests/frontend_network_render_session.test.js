'use strict';

const assert = require('node:assert/strict');
const versions = require('../web/method_versions.json');
const renderSession = require('../web/lib/networkRenderSession.js');

assert.equal(renderSession.SESSION_VERSION, versions.network_render_session);
assert.deepEqual(renderSession.normalizeQuery([' cg1 ', null, '', 'cg2']), ['cg1', 'cg2']);

const session = renderSession.createSession();
assert.deepEqual(session.snapshot(), {
    id: 0, status: 'idle', query: [], error: null, meta: {},
});

const first = session.begin(['cg1']);
assert.equal(session.isActive(first.id), true);
assert.equal(first.signal.aborted, false);

const second = session.begin('cg2');
assert.equal(first.signal.aborted, true, 'a newer render must cancel the previous request');
assert.equal(session.isActive(first.id), false);
assert.equal(session.complete(first.id), false, 'a stale render cannot commit');
assert.equal(session.isActive(second.id), true);
assert.equal(session.complete(second.id, { nodes: 4, edges: 3 }), true);
assert.deepEqual(session.snapshot(), {
    id: second.id,
    status: 'ready',
    query: ['cg2'],
    error: null,
    meta: { nodes: 4, edges: 3 },
});

const failed = session.begin('cg3');
assert.equal(session.fail(failed.id, 'empty-network'), true);
assert.equal(session.snapshot().status, 'error');
assert.equal(session.snapshot().error, 'empty-network');

const resetTarget = session.begin('cg4');
session.reset();
assert.equal(resetTarget.signal.aborted, true);
assert.equal(session.snapshot().status, 'idle');
assert.deepEqual(session.snapshot().query, []);

console.log('frontend network render session tests passed');
