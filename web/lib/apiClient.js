(function attachApiClient(root, factory) {
    const api = factory();
    if (typeof module === 'object' && module.exports) module.exports = api;
    if (root) root.CglApiClient = api;
})(typeof window !== 'undefined' ? window : globalThis, function createApiClientModule() {
    'use strict';

    class ApiError extends Error {
        constructor(message, { status = 0, url = '', payload = null, cause = null } = {}) {
            super(message);
            this.name = 'ApiError';
            this.status = status;
            this.url = url;
            this.payload = payload;
            if (cause) this.cause = cause;
        }
    }

    function errorMessage(payload, fallback) {
        if (payload && typeof payload === 'object') {
            return String(payload.detail || payload.error || payload.message || fallback);
        }
        if (typeof payload === 'string' && payload.trim()) return payload.trim();
        return fallback;
    }

    function createClient({ fetchImpl, defaultTimeoutMs = 15000 } = {}) {
        const performFetch = fetchImpl || (
            typeof fetch === 'function' ? fetch.bind(globalThis) : null
        );
        if (!performFetch) throw new Error('A fetch implementation is required');

        async function request(url, options = {}) {
            const {
                responseType = 'json',
                timeoutMs = defaultTimeoutMs,
                ...fetchOptions
            } = options;
            const controller = new AbortController();
            const externalSignal = fetchOptions.signal;
            const abortFromExternal = () => controller.abort(externalSignal?.reason);
            if (externalSignal) {
                if (externalSignal.aborted) abortFromExternal();
                else externalSignal.addEventListener('abort', abortFromExternal, { once: true });
            }
            fetchOptions.signal = controller.signal;
            const timer = timeoutMs > 0
                ? setTimeout(() => controller.abort(new Error('Request timed out')), timeoutMs)
                : null;

            let response;
            try {
                response = await performFetch(url, fetchOptions);
            } catch (cause) {
                const timedOut = controller.signal.aborted && !externalSignal?.aborted;
                throw new ApiError(
                    timedOut ? `Request timed out after ${timeoutMs} ms` : 'Network request failed',
                    { url: String(url), cause }
                );
            } finally {
                if (timer) clearTimeout(timer);
                if (externalSignal) externalSignal.removeEventListener('abort', abortFromExternal);
            }

            let payload = null;
            if (response.status !== 204) {
                try {
                    payload = responseType === 'text' ? await response.text() : await response.json();
                } catch (cause) {
                    if (response.ok) {
                        throw new ApiError(`Invalid ${responseType} response`, {
                            status: response.status, url: String(url), cause,
                        });
                    }
                }
            }
            if (!response.ok) {
                throw new ApiError(
                    errorMessage(payload, `API returned HTTP ${response.status}`),
                    { status: response.status, url: String(url), payload }
                );
            }
            return payload;
        }

        const getJson = (url, options = {}) => request(url, { ...options, responseType: 'json' });
        const getText = (url, options = {}) => request(url, { ...options, responseType: 'text' });
        const postJson = (url, body, options = {}) => request(url, {
            ...options,
            method: 'POST',
            headers: { 'Content-Type': 'application/json', ...(options.headers || {}) },
            body: JSON.stringify(body),
            responseType: 'json',
        });

        return { request, getJson, getText, postJson };
    }

    const defaultClient = createClient();
    return { ApiError, createClient, ...defaultClient };
});
