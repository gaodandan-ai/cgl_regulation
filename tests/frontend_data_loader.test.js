'use strict';

const assert = require('node:assert/strict');
const { RequiredAssetError, loadAssets } = require('../web/lib/dataLoader.js');

async function run() {
    const client = {
        async getJson(url) {
            if (url === '/missing') {
                const error = new Error('Not found');
                error.status = 404;
                throw error;
            }
            return { url };
        },
        async getText(url) { return `text:${url}`; },
    };
    const result = await loadAssets({
        json: { url: '/data.json', fallback: () => ({}) },
        text: { url: '/data.csv', type: 'text', fallback: '' },
        optional: { url: '/missing', fallback: () => [] },
        custom: { load: async () => 42, fallback: null },
    }, { client });
    assert.deepEqual(result.values.json, { url: '/data.json' });
    assert.equal(result.values.text, 'text:/data.csv');
    assert.deepEqual(result.values.optional, []);
    assert.equal(result.values.custom, 42);
    assert.deepEqual(result.failures, [{
        key: 'optional', url: '/missing', message: 'Not found', status: 404,
    }]);

    await assert.rejects(
        () => loadAssets({
            core: { url: '/missing', required: true, errorMessage: 'Core unavailable' },
        }, { client }),
        error => error instanceof RequiredAssetError
            && error.key === 'core'
            && error.message === 'Core unavailable'
    );

    const started = [];
    const resolvers = {};
    const concurrentClient = {
        getJson(url) {
            started.push(url);
            return new Promise(resolve => { resolvers[url] = resolve; });
        },
    };
    const pending = loadAssets({
        first: { url: '/first' },
        second: { url: '/second' },
    }, { client: concurrentClient });
    await Promise.resolve();
    assert.deepEqual(started, ['/first', '/second']);
    resolvers['/first'](1);
    resolvers['/second'](2);
    assert.deepEqual((await pending).values, { first: 1, second: 2 });

    console.log('frontend data loader tests passed');
}

run().catch(error => {
    console.error(error);
    process.exitCode = 1;
});
