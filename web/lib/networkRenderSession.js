(function attachNetworkRenderSession(root, factory) {
    const moduleApi = factory();
    if (typeof module === 'object' && module.exports) module.exports = moduleApi;
    if (root) root.CglNetworkRenderSession = moduleApi;
})(typeof window !== 'undefined' ? window : globalThis, function createNetworkRenderSessionModule() {
    'use strict';

    const SESSION_VERSION = 'network-render-session-v1.0.0';

    function normalizeQuery(query) {
        const values = Array.isArray(query) ? query : [query];
        return values
            .map(value => value === undefined || value === null ? '' : String(value).trim())
            .filter(Boolean);
    }

    function createSession() {
        let sequence = 0;
        let active = null;
        let state = {
            id: 0,
            status: 'idle',
            query: [],
            error: null,
            meta: {},
        };

        function snapshot() {
            return {
                ...state,
                query: [...state.query],
                meta: { ...state.meta },
            };
        }

        function begin(query) {
            if (active) active.controller.abort('superseded');
            const controller = new AbortController();
            const transaction = {
                id: ++sequence,
                query: normalizeQuery(query),
                controller,
                signal: controller.signal,
            };
            active = transaction;
            state = {
                id: transaction.id,
                status: 'loading',
                query: [...transaction.query],
                error: null,
                meta: {},
            };
            return transaction;
        }

        function isActive(id) {
            return Boolean(active && active.id === id && !active.signal.aborted);
        }

        function settle(id, status, { error = null, meta = {} } = {}) {
            if (!isActive(id)) return false;
            state = {
                id,
                status,
                query: [...active.query],
                error: error ? String(error) : null,
                meta: { ...meta },
            };
            active = null;
            return true;
        }

        function complete(id, meta = {}) {
            return settle(id, 'ready', { meta });
        }

        function fail(id, error) {
            return settle(id, 'error', { error });
        }

        function reset() {
            if (active) active.controller.abort('reset');
            active = null;
            state = { id: sequence, status: 'idle', query: [], error: null, meta: {} };
        }

        return { begin, isActive, complete, fail, reset, snapshot };
    }

    return { SESSION_VERSION, normalizeQuery, createSession };
});
