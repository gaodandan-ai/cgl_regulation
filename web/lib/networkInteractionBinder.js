(function attachNetworkInteractionBinder(root, factory) {
    const moduleApi = factory();
    if (typeof module === 'object' && module.exports) module.exports = moduleApi;
    if (root) root.CglNetworkInteractionBinder = moduleApi;
})(typeof window !== 'undefined' ? window : globalThis, function createNetworkInteractionBinderModule() {
    'use strict';

    const BINDER_VERSION = 'network-interaction-binder-v1.0.0';

    function bindLevelOfDetail(graph, { maxLabeledNodes = 60, minLabelZoom = 0.45 } = {}) {
        const updateLabels = () => {
            const nodes = graph.nodes();
            if (nodes.length <= maxLabeledNodes) return;
            if (graph.zoom() < minLabelZoom) nodes.addClass('lod-hide-label');
            else nodes.removeClass('lod-hide-label');
        };
        graph.on('zoom', updateLabels);
        return updateLabels;
    }

    function markSharedTargets(graph) {
        let marked = 0;
        graph.nodes('[type="Target"]').forEach(node => {
            if (node.indegree(false) > 1) {
                node.addClass('shared-target');
                marked += 1;
            }
        });
        return marked;
    }

    function bindInteractions(graph, {
        highlightSubnet,
        showNodeDetails,
        querySingleGene,
        toggleRightSidebar,
        doubleTapWindowMs = 350,
        now = Date.now,
    }) {
        let lastTapNode = null;
        let lastTapAt = null;

        const handleNodeTap = event => {
            const node = event.target;
            const tappedAt = now();
            highlightSubnet(node);
            showNodeDetails(node.id());
            if (typeof toggleRightSidebar === 'function') toggleRightSidebar(true);
            if (lastTapNode === node && lastTapAt !== null && tappedAt - lastTapAt < doubleTapWindowMs) {
                querySingleGene(node.id());
            } else {
                lastTapNode = node;
                lastTapAt = tappedAt;
            }
        };

        const handleCanvasTap = event => {
            if (event.target === graph) toggleRightSidebar(false);
        };

        graph.on('tap', 'node', handleNodeTap);
        graph.on('tap', handleCanvasTap);
        return { handleNodeTap, handleCanvasTap };
    }

    return { BINDER_VERSION, bindLevelOfDetail, markSharedTargets, bindInteractions };
});
