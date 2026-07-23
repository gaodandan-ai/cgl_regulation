'use strict';

const assert = require('node:assert/strict');
const { escapeCsvCell, toCsv, download } = require('../web/lib/exportUtils.js');

async function run() {
    assert.equal(escapeCsvCell('plain'), 'plain');
    assert.equal(escapeCsvCell('a,b'), '"a,b"');
    assert.equal(escapeCsvCell('a"b'), '"a""b"');
    assert.equal(escapeCsvCell('=HYPERLINK("bad")'), '"\'=HYPERLINK(""bad"")"');
    assert.equal(escapeCsvCell('+cmd'), "'+cmd");
    assert.equal(toCsv([['gene', 'score'], ['cg0350', 0.9]]), 'gene,score\ncg0350,0.9');
    assert.equal(toCsv([['a']], { bom: true }), '\uFEFFa');

    const events = [];
    const anchor = {
        style: {},
        click() { events.push('click'); },
        remove() { events.push('remove'); },
    };
    const documentRef = {
        body: { appendChild(item) { assert.equal(item, anchor); events.push('append'); } },
        createElement(name) { assert.equal(name, 'a'); return anchor; },
    };
    const urlApi = {
        createObjectURL(blob) { assert.ok(blob); events.push('create'); return 'blob:test'; },
        revokeObjectURL(url) { assert.equal(url, 'blob:test'); events.push('revoke'); },
    };
    class BlobImpl {
        constructor(parts, options) { this.parts = parts; this.options = options; }
    }
    download('content', 'report.csv', 'text/csv', { documentRef, urlApi, BlobImpl });
    await new Promise(resolve => setTimeout(resolve, 5));
    assert.equal(anchor.download, 'report.csv');
    assert.deepEqual(events, ['create', 'append', 'click', 'remove', 'revoke']);

    console.log('frontend export utility tests passed');
}

run().catch(error => {
    console.error(error);
    process.exitCode = 1;
});
