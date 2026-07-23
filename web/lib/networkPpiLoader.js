(function attachNetworkPpiLoader(root, factory) {
    const moduleApi = factory();
    if (typeof module === 'object' && module.exports) module.exports = moduleApi;
    if (root) root.CglNetworkPpiLoader = moduleApi;
})(typeof window !== 'undefined' ? window : globalThis, function createNetworkPpiLoaderModule() {
    'use strict';

    const LOADER_VERSION = 'network-ppi-loader-v1.0.0';

    function normalizeQuery(query) {
        return (Array.isArray(query) ? query : [query])
            .map(value => value === undefined || value === null ? '' : String(value).trim())
            .filter(Boolean);
    }

    function partnerUrl(locus) {
        return `/api/analysis/string_ppi?gene=${encodeURIComponent(locus)}`;
    }

    function visibleEdgesUrl(loci) {
        return `/api/analysis/network_ppi_edges?genes=${encodeURIComponent(loci.join(','))}`;
    }

    async function loadQueryInteractions({
        query,
        enabled,
        client,
        signal,
        onWarning = () => {},
    }) {
        if (!enabled) return [];
        const loci = normalizeQuery(query);
        const batches = await Promise.all(loci.map(async locus => {
            try {
                const data = await client.getJson(partnerUrl(locus), { signal });
                const partners = Array.isArray(data?.partners) ? data.partners : [];
                return partners
                    .filter(partner => partner && partner.partner)
                    .map(partner => ({
                        source: locus,
                        target: partner.partner,
                        score: partner.score,
                        type: partner.type,
                    }));
            } catch (error) {
                if (signal?.aborted) throw error;
                onWarning(locus, error);
                return [];
            }
        }));
        return batches.flat();
    }

    async function loadVisibleEdges({ graph, enabled, client, signal }) {
        if (!enabled) return [];
        const loci = graph.nodes().map(node => node.id()).filter(Boolean);
        if (loci.length <= 1) return [];
        const data = await client.getJson(visibleEdgesUrl(loci), { signal });
        return Array.isArray(data?.edges) ? data.edges : [];
    }

    return {
        LOADER_VERSION,
        normalizeQuery,
        partnerUrl,
        visibleEdgesUrl,
        loadQueryInteractions,
        loadVisibleEdges,
    };
});
