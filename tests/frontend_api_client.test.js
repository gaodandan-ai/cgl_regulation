'use strict';

const assert = require('node:assert/strict');
const { ApiError, createClient } = require('../web/lib/apiClient.js');

function response({ ok = true, status = 200, json, text }) {
    return {
        ok,
        status,
        async json() {
            if (json instanceof Error) throw json;
            return json;
        },
        async text() { return text; },
    };
}

async function run() {
    let captured;
    const client = createClient({
        fetchImpl: async (url, options) => {
            captured = { url, options };
            return response({ json: { result: 'ok' } });
        },
    });
    assert.deepEqual(await client.postJson('/api/test', { gene: 'cg0350' }), { result: 'ok' });
    assert.equal(captured.options.method, 'POST');
    assert.equal(captured.options.headers['Content-Type'], 'application/json');
    assert.deepEqual(JSON.parse(captured.options.body), { gene: 'cg0350' });

    const failed = createClient({
        fetchImpl: async () => response({ ok: false, status: 422, json: { detail: 'Bad gene' } }),
    });
    await assert.rejects(
        () => failed.getJson('/api/test'),
        error => error instanceof ApiError && error.status === 422 && error.message === 'Bad gene'
    );

    const textClient = createClient({
        fetchImpl: async () => response({ text: 'a,b\n1,2\n' }),
    });
    assert.equal(await textClient.getText('/data/test.csv'), 'a,b\n1,2\n');

    const timeoutClient = createClient({
        defaultTimeoutMs: 5,
        fetchImpl: (_url, options) => new Promise((_resolve, reject) => {
            options.signal.addEventListener('abort', () => reject(options.signal.reason), { once: true });
        }),
    });
    await assert.rejects(
        () => timeoutClient.getJson('/api/slow'),
        error => error instanceof ApiError && error.message.includes('timed out')
    );

    console.log('frontend api client tests passed');
}

run().catch(error => {
    console.error(error);
    process.exitCode = 1;
});
