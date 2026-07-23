(function attachQueryNavigation(root, factory) {
    const moduleApi = factory();
    if (typeof module === 'object' && module.exports) module.exports = moduleApi;
    if (root) root.CglQueryNavigation = moduleApi;
})(typeof window !== 'undefined' ? window : globalThis, function createQueryNavigationModule() {
    'use strict';

    const STATE_VERSION = 'query-navigation-v1.0.0';

    function normalizeQueryList(value) {
        if (value === undefined || value === null || value === '') return [];
        return (Array.isArray(value) ? value : [value]).slice();
    }

    function splitGeneQuery(value) {
        return normalizeQueryList(value)
            .flatMap(item => String(item).split(','))
            .map(item => item.trim())
            .filter(Boolean);
    }

    function serializeGeneQuery(value) {
        const genes = splitGeneQuery(value);
        return genes.length ? genes.join(',') : null;
    }

    function queryKey(value) {
        return splitGeneQuery(value)
            .map(item => item.toLowerCase())
            .sort()
            .join('\u0000');
    }

    function queriesEqual(left, right) {
        return queryKey(left) === queryKey(right);
    }

    function parseUrlState(search) {
        const params = search instanceof URLSearchParams
            ? search
            : new URLSearchParams(String(search || '').replace(/^.*\?/, '?'));
        const workflow = params.get('workflow') || 'gene';
        const gene = serializeGeneQuery(params.get('gene'));
        return { workflow, gene, genes: splitGeneQuery(gene), version: STATE_VERSION };
    }

    function buildUrlState(currentHref, state = {}) {
        const url = new URL(currentHref);
        const workflow = state.workflow || 'gene';
        const gene = serializeGeneQuery(state.gene);
        if (state.workflow) url.searchParams.set('workflow', workflow);
        if (gene) url.searchParams.set('gene', gene);
        else if (workflow !== 'gene') url.searchParams.delete('gene');
        return {
            href: url.href,
            state: { workflow, gene, version: STATE_VERSION },
        };
    }

    function createHistory() {
        const backStack = [];
        const forwardStack = [];
        let suspensionDepth = 0;

        function record(current, next) {
            const currentList = normalizeQueryList(current);
            if (suspensionDepth || !currentList.length || queriesEqual(currentList, next)) {
                return false;
            }
            backStack.push(currentList);
            forwardStack.length = 0;
            return true;
        }

        function go(direction, current) {
            const source = direction === 'back' ? backStack : forwardStack;
            const destination = direction === 'back' ? forwardStack : backStack;
            if (!source.length) return null;
            const target = source.pop();
            const currentList = normalizeQueryList(current);
            if (currentList.length) destination.push(currentList);
            return target.slice();
        }

        function suspend() { suspensionDepth += 1; }
        function resume() { suspensionDepth = Math.max(0, suspensionDepth - 1); }
        function snapshot() {
            return {
                canBack: backStack.length > 0,
                canForward: forwardStack.length > 0,
                backDepth: backStack.length,
                forwardDepth: forwardStack.length,
                suspended: suspensionDepth > 0,
            };
        }

        return { record, go, suspend, resume, snapshot };
    }

    return {
        STATE_VERSION,
        normalizeQueryList,
        splitGeneQuery,
        serializeGeneQuery,
        queriesEqual,
        parseUrlState,
        buildUrlState,
        createHistory,
    };
});
