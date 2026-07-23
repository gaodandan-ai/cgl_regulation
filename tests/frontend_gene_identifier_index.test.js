const assert = require('assert');
const fs = require('fs');
const path = require('path');

const identifierModule = require('../web/lib/geneIdentifierIndex.js');
const methodVersions = JSON.parse(
    fs.readFileSync(path.join(__dirname, '..', 'web', 'method_versions.json'), 'utf8')
);

assert.strictEqual(
    identifierModule.INDEX_VERSION,
    methodVersions.gene_identifier_index,
    'the runtime identifier-index version must match provenance metadata'
);
assert.strictEqual(identifierModule.normalizeCgl('CGL1234'), 'cgl1234');

const index = identifierModule.createIndex({
    geneMappings: [
        { cgl_locus: 'CGL0001', cg_locus: 'cg0001', gene_name: 'glxR', product: 'regulator' },
        { cgl_locus: 'cgl0002', cg_locus: 'cg0002', gene_name: 'dupMap', product: 'enzyme' },
        { cgl_locus: 'cgl0003', cg_locus: 'cg0003', gene_name: 'dupMap', product: 'other' },
    ],
    regulations: [
        { TF_locusTag: 'cg0001', TF_name: 'GlxR', TG_locusTag: 'cg0002', TG_name: 'dup' },
        { TF_locusTag: 'cg0002', TF_name: 'DualRole', TG_locusTag: 'cg0004', TG_name: 'target' },
        { TF_locusTag: 'cg0003', TF_name: 'dup', TG_locusTag: 'cg0005', TG_name: 'other' },
    ],
    rnaRegulations: [
        { srna: 'srnaA', mrna: 'cg0006' },
    ],
});

assert.strictEqual(index.resolve('CG0001').locusTag, 'cg0001');
assert.strictEqual(index.resolve('CGL0001').locusTag, 'cg0001');
assert.strictEqual(index.resolve('GLXR').locusTag, 'cg0001');
assert.strictEqual(index.resolve('unknown'), null);
assert.strictEqual(index.resolve('cg0002').type, 'TF', 'a later TF role must upgrade a target role');
assert.strictEqual(index.resolve('dup').locusTag, 'cg0003', 'the higher-priority TF alias must win');
assert.strictEqual(index.nameToCg.dupmap, 'cg0002', 'mapping collisions retain the first canonical mapping');
assert.strictEqual(index.getPrioritizedLabel('cg0001', 'GlxR'), 'cgl0001');
assert.strictEqual(index.cgToProduct.cg0001, 'regulator');
assert.strictEqual(index.cgToProduct.cgl0001, 'regulator');
assert(index.suggestions.some(item => item.display === 'cgl0001'));
assert(index.conflicts.some(item => item.alias === 'dup' && item.selected === 'cg0003'));
assert(index.conflicts.some(item => item.alias === 'dupmap' && item.selected === 'cg0002'));
assert.strictEqual(index.searchSuggestions('cg', 2).length, 2);
assert.deepStrictEqual(index.searchSuggestions('', 15), []);

console.log('frontend_gene_identifier_index.test.js: ok');
