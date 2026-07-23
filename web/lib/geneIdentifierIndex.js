(function attachGeneIdentifierIndex(root, factory) {
    const moduleApi = factory();
    if (typeof module === 'object' && module.exports) module.exports = moduleApi;
    if (root) root.CglGeneIdentifierIndex = moduleApi;
})(typeof window !== 'undefined' ? window : globalThis, function createGeneIdentifierIndexModule() {
    'use strict';

    const INDEX_VERSION = 'gene-identifier-index-v1.0.0';
    const TYPE_RANK = Object.freeze({ TF: 3, sRNA: 2, Target: 1 });

    function cleanString(value) {
        if (value === undefined || value === null || value === 'null') return '';
        return String(value).trim();
    }

    function normalizeCgl(value) {
        const cgl = cleanString(value);
        return cgl.toLowerCase().startsWith('cgl') ? `cgl${cgl.substring(3)}` : cgl;
    }

    function createIndex({ geneMappings = [], regulations = [], rnaRegulations = [] } = {}) {
        const geneIndex = {};
        const cglToCg = {};
        const cgToCgl = {};
        const nameToCg = {};
        const cgToProduct = {};
        const conflicts = [];

        function register(alias, record, source) {
            const key = cleanString(alias).toLowerCase();
            if (!key || !record.locusTag) return;
            const existing = geneIndex[key];
            if (!existing) {
                geneIndex[key] = { ...record };
                return;
            }
            if (existing.locusTag.toLowerCase() !== record.locusTag.toLowerCase()) {
                const replace = (TYPE_RANK[record.type] || 0) > (TYPE_RANK[existing.type] || 0);
                conflicts.push({
                    alias: key,
                    selected: replace ? record.locusTag : existing.locusTag,
                    rejected: replace ? existing.locusTag : record.locusTag,
                    source,
                });
                if (!replace) return;
            }
            geneIndex[key] = {
                ...existing,
                ...record,
                type: (TYPE_RANK[record.type] || 0) > (TYPE_RANK[existing.type] || 0)
                    ? record.type
                    : existing.type,
            };
        }

        geneMappings.forEach(row => {
            const cgl = normalizeCgl(row.cgl_locus);
            const cg = cleanString(row.cg_locus);
            const name = cleanString(row.gene_name);
            const product = cleanString(row.product);
            if (cgl && cg) {
                cglToCg[cgl.toLowerCase()] = cg;
                cgToCgl[cg.toLowerCase()] = cgl;
            }
            if (name && name !== '--' && cg) {
                const nameKey = name.toLowerCase();
                if (!nameToCg[nameKey]) nameToCg[nameKey] = cg;
                else if (nameToCg[nameKey].toLowerCase() !== cg.toLowerCase()) {
                    conflicts.push({
                        alias: nameKey, selected: nameToCg[nameKey], rejected: cg,
                        source: 'gene_mapping_name',
                    });
                }
            }
            if (cg && product) cgToProduct[cg.toLowerCase()] = product;
            if (cgl && product) cgToProduct[cgl.toLowerCase()] = product;
        });

        regulations.forEach(row => {
            const tfTag = cleanString(row.TF_locusTag);
            const tfName = cleanString(row.TF_name);
            const tgTag = cleanString(row.TG_locusTag);
            const tgName = cleanString(row.TG_name);
            if (tfTag) {
                const record = { locusTag: tfTag, name: tfName || tfTag, type: 'TF' };
                register(tfTag, record, 'TF_locusTag');
                if (tfName) register(tfName, record, 'TF_name');
            }
            if (tgTag) {
                const record = { locusTag: tgTag, name: tgName || tgTag, type: 'Target' };
                register(tgTag, record, 'TG_locusTag');
                if (tgName) register(tgName, record, 'TG_name');
            }
        });

        rnaRegulations.forEach(row => {
            const srna = cleanString(row.srna);
            const mrna = cleanString(row.mrna);
            if (srna) register(srna, { locusTag: srna, name: srna, type: 'sRNA' }, 'sRNA');
            if (mrna) register(mrna, { locusTag: mrna, name: mrna, type: 'Target' }, 'mRNA');
        });

        const uniqueSuggestions = {};
        Object.values(geneIndex).forEach(item => {
            uniqueSuggestions[item.locusTag] = item;
            if (item.name && item.name !== item.locusTag) uniqueSuggestions[item.name] = item;
            const cgl = cgToCgl[item.locusTag.toLowerCase()];
            if (cgl) uniqueSuggestions[cgl] = item;
        });
        const suggestions = Object.entries(uniqueSuggestions).map(([display, item]) => ({
            display,
            locusTag: item.locusTag,
            type: item.type,
            cgl: cgToCgl[item.locusTag.toLowerCase()] || '',
        })).sort((left, right) => left.display.localeCompare(right.display));

        function resolve(query) {
            const key = cleanString(query).toLowerCase();
            if (!key) return null;
            const canonical = cleanString(cglToCg[key] || nameToCg[key] || key).toLowerCase();
            return geneIndex[canonical] || geneIndex[key] || null;
        }

        function getPrioritizedLabel(locusTag, commonName) {
            const locus = cleanString(locusTag);
            if (!locus) return cleanString(commonName);
            const cgl = cgToCgl[locus.toLowerCase()];
            if (cgl) return cgl;
            const name = cleanString(commonName);
            return name && name !== locus && name !== '--' ? name : locus;
        }

        function searchSuggestions(query, limit = 15) {
            const key = cleanString(query).toLowerCase();
            if (!key) return [];
            return suggestions.filter(item =>
                item.display.toLowerCase().includes(key)
                || item.locusTag.toLowerCase().includes(key)
                || item.cgl.toLowerCase().includes(key)
            ).slice(0, Math.max(0, limit));
        }

        return {
            version: INDEX_VERSION,
            geneIndex,
            cglToCg,
            cgToCgl,
            nameToCg,
            cgToProduct,
            suggestions,
            conflicts,
            resolve,
            getPrioritizedLabel,
            searchSuggestions,
        };
    }

    return { INDEX_VERSION, TYPE_RANK, normalizeCgl, createIndex };
});
