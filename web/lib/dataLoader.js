(function attachDataLoader(root, factory) {
    const apiClient = root && root.CglApiClient
        ? root.CglApiClient
        : (typeof module === 'object' && module.exports ? require('./apiClient.js') : null);
    const moduleApi = factory(apiClient);
    if (typeof module === 'object' && module.exports) module.exports = moduleApi;
    if (root) root.CglDataLoader = moduleApi;
})(typeof window !== 'undefined' ? window : globalThis, function createDataLoaderModule(apiClient) {
    'use strict';

    class RequiredAssetError extends Error {
        constructor(key, descriptor, cause) {
            super(descriptor.errorMessage || `Required data asset failed to load: ${key}`);
            this.name = 'RequiredAssetError';
            this.key = key;
            this.url = descriptor.url || '';
            this.cause = cause;
        }
    }

    function fallbackValue(descriptor) {
        return typeof descriptor.fallback === 'function'
            ? descriptor.fallback()
            : descriptor.fallback;
    }

    async function loadOne(key, descriptor, client) {
        try {
            let value;
            if (typeof descriptor.load === 'function') {
                value = await descriptor.load();
            } else if (descriptor.type === 'text') {
                value = await client.getText(descriptor.url, descriptor.options);
            } else {
                value = await client.getJson(descriptor.url, descriptor.options);
            }
            return { key, value, failure: null };
        } catch (cause) {
            if (descriptor.required) throw new RequiredAssetError(key, descriptor, cause);
            return {
                key,
                value: fallbackValue(descriptor),
                failure: {
                    key,
                    url: descriptor.url || '',
                    message: cause instanceof Error ? cause.message : String(cause),
                    status: Number(cause && cause.status) || 0,
                },
            };
        }
    }

    async function loadAssets(specification, { client = apiClient } = {}) {
        if (!client) throw new Error('CglApiClient is required');
        const entries = Object.entries(specification || {});
        const settled = await Promise.all(
            entries.map(([key, descriptor]) => loadOne(key, descriptor, client))
        );
        const values = {};
        const failures = [];
        for (const result of settled) {
            values[result.key] = result.value;
            if (result.failure) failures.push(result.failure);
        }
        return { values, failures };
    }

    return { RequiredAssetError, loadAssets };
});
