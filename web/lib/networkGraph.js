(function attachNetworkGraph(root, factory) {
    const moduleApi = factory();
    if (typeof module === 'object' && module.exports) module.exports = moduleApi;
    if (root) root.CglNetworkGraph = moduleApi;
})(typeof window !== 'undefined' ? window : globalThis, function createNetworkGraphModule() {
    'use strict';

    const GRAPH_VERSION = 'network-graph-v1.0.0';

    function elementCount(elements) {
        if (Array.isArray(elements)) return elements.length;
        if (!elements || typeof elements !== 'object') return 0;
        const nodes = Array.isArray(elements.nodes) ? elements.nodes.length : 0;
        const edges = Array.isArray(elements.edges) ? elements.edges.length : 0;
        return nodes + edges;
    }

    function createOptions({
        container,
        elements,
        styles,
        layoutName,
        largeNetworkThreshold = 250,
    }) {
        const count = elementCount(elements);
        return {
            container,
            elements,
            style: styles,
            textureOnViewport: true,
            hideEdgesOnViewport: count > largeNetworkThreshold,
            pixelRatio: 'auto',
            wheelSensitivity: 0.2,
            minZoom: 0.15,
            maxZoom: 2.0,
            layout: {
                name: layoutName,
                animate: true,
                animationDuration: 400,
                fit: count > 3,
                padding: count <= 3 ? 120 : 40,
            },
        };
    }

    function safeFit(graph, padding = 40) {
        if (!graph || typeof graph.nodes !== 'function') return;
        const nodes = graph.nodes();
        if (nodes.length <= 3) {
            graph.center();
            graph.zoom(1.0);
            graph.center();
        } else {
            graph.fit(undefined, padding);
            if (graph.zoom() > 2.0) {
                graph.zoom(2.0);
                graph.center();
            }
        }
    }

    function createGraph({ cytoscapeImpl, ...options }) {
        if (typeof cytoscapeImpl !== 'function') throw new TypeError('cytoscapeImpl must be a function');
        const graph = cytoscapeImpl(createOptions(options));

        const clampZoom = () => safeFit(graph);

        if (graph && typeof graph.on === 'function') {
            graph.on('layoutstop', clampZoom);
        }
        setTimeout(clampZoom, 50);
        setTimeout(clampZoom, 450);

        return graph;
    }

    function ppiEdge(record) {
        if (!record || !record.source || !record.target) return null;
        return {
            group: 'edges',
            data: {
                id: `ppi-cross-${record.source}-${record.target}`,
                source: record.source,
                target: record.target,
                role: 'protein-protein interaction',
                type: 'PPI',
                regulationType: 'ppi',
                score: record.score,
                schemaVersion: 'unified-v1',
            },
        };
    }

    function addPpiEdges(graph, records = []) {
        const pending = [];
        const reservedIds = new Set();
        records.forEach(record => {
            const edge = ppiEdge(record);
            if (!edge) return;
            const forwardId = edge.data.id;
            const reverseId = `ppi-cross-${edge.data.target}-${edge.data.source}`;
            const alreadyPresent = graph.getElementById(forwardId).length
                || graph.getElementById(reverseId).length
                || reservedIds.has(forwardId)
                || reservedIds.has(reverseId);
            if (alreadyPresent) return;
            reservedIds.add(forwardId);
            reservedIds.add(reverseId);
            pending.push(edge);
        });
        if (pending.length) graph.add(pending);
        return pending.length;
    }

    return { GRAPH_VERSION, elementCount, createOptions, createGraph, ppiEdge, addPpiEdges, safeFit };
});
