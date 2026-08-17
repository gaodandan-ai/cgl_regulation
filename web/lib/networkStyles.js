(function attachNetworkStyles(root, factory) {
    const moduleApi = factory();
    if (typeof module === 'object' && module.exports) module.exports = moduleApi;
    if (root) root.CglNetworkStyles = moduleApi;
})(typeof window !== 'undefined' ? window : globalThis, function createNetworkStylesModule() {
    'use strict';

    const STYLE_VERSION = 'network-styles-v2.0.0';

    function createBaseNodeStyles() {
        return [
            {
                selector: 'node',
                style: {
                    label: 'data(name)',
                    'font-size': '11px',
                    color: '#0f172a',
                    'background-color': '#f5f5f5',
                    'text-valign': 'bottom',
                    'text-margin-y': '6px',
                    width: '22px',
                    height: '22px',
                    'border-width': '2px',
                    'border-color': '#757575',
                    'transition-property': 'background-color, line-color, target-arrow-color, width, height, border-width',
                    'transition-duration': '0.2s',
                },
            },
            {
                selector: 'node[type="TF"]',
                style: {
                    'background-color': '#e3f2fd',
                    'border-color': '#1976d2',
                    width: '26px',
                    height: '26px',
                },
            },
            {
                selector: 'node[type="sRNA"]',
                style: {
                    'background-color': '#f3e5f5',
                    'border-color': '#8e24aa',
                    width: '26px',
                    height: '26px',
                    shape: 'hexagon',
                },
            },
            {
                selector: 'node[type="query"]',
                style: {
                    'background-color': '#ffe0b2',
                    'border-color': '#f57c00',
                    width: '34px',
                    height: '34px',
                    'border-width': '3px',
                    'font-weight': 'bold',
                    'font-size': '13px',
                },
            },
            {
                selector: 'node.shared-target',
                style: {
                    'background-color': '#e0f2f1',
                    'border-color': '#00897b',
                    'border-width': '2.5px',
                },
            },
        ];
    }

    function edgeConfidence(edge) {
        return edge.data('confidenceScore') || 0.25;
    }

    function createBaseEdgeStyles() {
        return [{
            selector: 'edge',
            style: {
                width: edge => 1.2 + edgeConfidence(edge) * 3.2,
                'line-color': '#e65100',
                'target-arrow-color': '#e65100',
                'target-arrow-shape': 'triangle',
                'curve-style': 'bezier',
                'arrow-scale': 1.1,
                opacity: edge => 0.35 + edgeConfidence(edge) * 0.6,
                'transition-property': 'line-color, target-arrow-color, opacity, width',
                'transition-duration': '0.2s',
            },
        }];
    }

    function createRegulationEdgeStyles() {
        return [
            {
                selector: 'edge[regulationType="ppi"]',
                style: {
                    'line-color': '#10b981',
                    'line-style': 'dashed',
                    'line-dash-pattern': [10, 5],
                    'target-arrow-shape': 'none',
                    width: edge => 2 + (((edge.data('score') || 700) - 700) / 300) * 2,
                    opacity: 0.85,
                    'curve-style': 'bezier',
                },
            },
            {
                selector: 'edge[regulationType="activation"]',
                style: {
                    'line-color': '#2e7d32',
                    'target-arrow-color': '#2e7d32',
                    'target-arrow-shape': 'triangle',
                },
            },
            {
                selector: 'edge[role="A"]',
                style: { 'line-color': '#2e7d32', 'target-arrow-color': '#2e7d32' },
            },
            {
                selector: 'edge[regulationType="repression"]',
                style: {
                    'line-color': '#d32f2f',
                    'target-arrow-color': '#d32f2f',
                    'target-arrow-shape': 'tee',
                },
            },
            {
                selector: 'edge[role="R"]',
                style: {
                    'line-color': '#d32f2f',
                    'target-arrow-color': '#d32f2f',
                    'target-arrow-shape': 'tee',
                },
            },
            {
                selector: 'edge[regulationType="dual"], edge[regulationType="sigma"], edge[regulationType="unknown"]',
                style: {
                    'line-color': '#e65100',
                    'target-arrow-color': '#e65100',
                    'target-arrow-shape': 'triangle',
                },
            },
            {
                selector: 'edge[role="Dual"]',
                style: { 'line-color': '#e65100', 'target-arrow-color': '#e65100' },
            },
            {
                selector: 'edge[regulationType="post_transcriptional_repression"]',
                style: {
                    'line-color': '#7b1fa2',
                    'target-arrow-color': '#7b1fa2',
                    'line-style': 'dashed',
                    'target-arrow-shape': 'triangle-tee',
                },
            },
            {
                selector: 'edge[role="sRNA"]',
                style: {
                    'line-color': '#7b1fa2',
                    'target-arrow-color': '#7b1fa2',
                    'line-style': 'dashed',
                    'target-arrow-shape': 'triangle-tee',
                },
            },
            { selector: 'edge.confidence-high', style: { 'line-style': 'solid' } },
            { selector: 'edge.confidence-medium', style: { 'line-style': 'solid' } },
            { selector: 'edge.confidence-low', style: { 'line-style': 'dotted', opacity: 0.42 } },
        ];
    }

    function createInteractionStateStyles() {
        return [
            { selector: '.dimmed', style: { opacity: 0.15 } },
            { selector: '.rnaseq-hidden', style: { display: 'none' } },
            {
                selector: 'node.highlighted',
                style: {
                    'border-width': '3px',
                    'border-color': '#0f172a',
                    width: '38px',
                    height: '38px',
                },
            },
            { selector: 'edge.highlighted', style: { width: 3.5, opacity: 1.0 } },
            {
                selector: 'node.sim-up',
                style: {
                    'border-color': '#2e7d32',
                    'border-width': '4px',
                    'background-color': '#e8f5e9',
                    'shadow-blur': '10px',
                    'shadow-color': '#2e7d32',
                    'shadow-opacity': 0.8,
                },
            },
            {
                selector: 'node.sim-down',
                style: {
                    'border-color': '#d32f2f',
                    'border-width': '4px',
                    'background-color': '#ffebee',
                    'shadow-blur': '10px',
                    'shadow-color': '#d32f2f',
                    'shadow-opacity': 0.8,
                },
            },
            {
                selector: 'node.sim-dual',
                style: {
                    'border-color': '#e65100',
                    'border-width': '4px',
                    'background-color': '#fff3e0',
                    'shadow-blur': '10px',
                    'shadow-color': '#e65100',
                    'shadow-opacity': 0.8,
                },
            },
        ];
    }

    function numericThreshold(value, fallback) {
        const parsed = Number.parseFloat(value);
        return Number.isFinite(parsed) ? parsed : fallback;
    }

    function baseNodeSize(type) {
        if (type === 'query') return 34;
        if (type === 'TF' || type === 'sRNA') return 26;
        return 22;
    }

    function expressionNodeSize(log2FoldChange, type) {
        const base = baseNodeSize(type);
        if (log2FoldChange === undefined || Number.isNaN(Number(log2FoldChange))) return base;
        return base + Math.min(16, Math.abs(Number(log2FoldChange)) * 4);
    }

    function createRnaSeqStyles({ colorForLog2FoldChange, thresholdValue }) {
        const threshold = (id, fallback) => numericThreshold(thresholdValue(id, fallback), fallback);
        return [{
            selector: 'node.rnaseq-node',
            style: {
                'background-color': node => colorForLog2FoldChange(node.data('rnaseq_log2fc')),
                'border-width': node => {
                    const pValue = node.data('rnaseq_pvalue');
                    return pValue !== undefined && pValue <= threshold('rnaseq-p-threshold', 0.05)
                        ? '3.5px' : '2px';
                },
                'border-color': node => {
                    const pValue = node.data('rnaseq_pvalue');
                    return pValue !== undefined && pValue <= threshold('rnaseq-p-threshold', 0.05)
                        ? '#0f172a' : '#94a3b8';
                },
                'width': node => expressionNodeSize(node.data('rnaseq_log2fc'), node.data('type')),
                'height': node => expressionNodeSize(node.data('rnaseq_log2fc'), node.data('type')),
                'shadow-blur': node => {
                    const pValue = node.data('rnaseq_pvalue');
                    const log2FoldChange = node.data('rnaseq_log2fc');
                    const significant = pValue !== undefined
                        && pValue <= threshold('rnaseq-p-threshold', 0.05)
                        && Math.abs(log2FoldChange) >= threshold('rnaseq-lfc-threshold', 1.0);
                    return significant ? '12px' : '0px';
                },
                'shadow-color': node => {
                    const log2FoldChange = node.data('rnaseq_log2fc');
                    if (log2FoldChange === undefined) return 'transparent';
                    return log2FoldChange > 0 ? '#ef4444' : '#2563eb';
                },
                'shadow-opacity': 0.85,
                'shadow-offset-x': '0px',
                'shadow-offset-y': '0px',
            },
        }];
    }

    function getStyles({ theme = 'light' } = {}) {
        return [
            ...createBaseNodeStyles(),
            ...createBaseEdgeStyles(),
            ...createRegulationEdgeStyles(),
            ...createInteractionStateStyles(),
        ];
    }

    return {
        STYLE_VERSION,
        createBaseNodeStyles,
        createBaseEdgeStyles,
        createRegulationEdgeStyles,
        createInteractionStateStyles,
        numericThreshold,
        baseNodeSize,
        expressionNodeSize,
        createRnaSeqStyles,
        getStyles,
    };
});
