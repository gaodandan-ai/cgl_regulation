(function attachDandelionLayout(root, factory) {
    const moduleApi = factory();
    if (typeof module === 'object' && module.exports) module.exports = moduleApi;
    if (root) root.CglDandelionLayout = moduleApi;
})(typeof window !== 'undefined' ? window : globalThis, function createDandelionLayoutModule() {
    'use strict';

    const LAYOUT_VERSION = 'dandelion-layout-v1.0.0';

    const DEFAULT_SECTORS = {
        activation: { startAngle: -Math.PI / 4, span: Math.PI / 2.2, label: 'Activation (+)' },
        repression: { startAngle: (3 * Math.PI) / 4, span: Math.PI / 2.2, label: 'Repression (-)' },
        dual: { startAngle: Math.PI / 4, span: Math.PI / 2.5, label: 'Dual / Sigma (±)' },
        srna: { startAngle: (-3 * Math.PI) / 4, span: Math.PI / 2.5, label: 'sRNA' },
        ppi: { startAngle: Math.PI / 2, span: Math.PI / 3, label: 'PPI' },
        other: { startAngle: -Math.PI / 2, span: Math.PI / 3, label: 'Target / Other' },
    };

    function getNodeCategory(node, querySet) {
        if (!node) return 'other';
        const id = (node.id() || '').toLowerCase();
        if (node.data('type') === 'query' || querySet.has(id)) return 'center';

        const edges = node.connectedEdges();
        if (!edges || edges.length === 0) return 'other';

        let hasActivation = false;
        let hasRepression = false;
        let hasSrna = false;
        let hasPpi = false;

        edges.forEach(edge => {
            const regType = (edge.data('regulationType') || edge.data('role') || '').toLowerCase();
            if (regType === 'activation' || regType === 'a') hasActivation = true;
            else if (regType === 'repression' || regType === 'r') hasRepression = true;
            else if (regType === 'dual' || regType === 'sigma') { hasActivation = true; hasRepression = true; }
            else if (regType === 'srna' || regType.includes('post_transcriptional')) hasSrna = true;
            else if (regType === 'ppi') hasPpi = true;
        });

        if (hasSrna) return 'srna';
        if (hasActivation && hasRepression) return 'dual';
        if (hasActivation) return 'activation';
        if (hasRepression) return 'repression';
        if (hasPpi) return 'ppi';
        return 'other';
    }

    function computeDandelionPositions(cy, queryLoci = []) {
        if (!cy || typeof cy.nodes !== 'function') return {};

        const queryList = Array.isArray(queryLoci) ? queryLoci : (queryLoci ? [queryLoci] : []);
        const querySet = new Set(queryList.map(q => (q || '').toLowerCase()));
        const nodes = cy.nodes();
        const positions = {};

        if (nodes.length === 0) return positions;

        const centerPos = { x: 0, y: 0 };
        const categoryGroups = {
            activation: [],
            repression: [],
            dual: [],
            srna: [],
            ppi: [],
            other: [],
        };

        nodes.forEach(node => {
            const cat = getNodeCategory(node, querySet);
            if (cat === 'center') {
                positions[node.id()] = { x: centerPos.x, y: centerPos.y };
            } else if (categoryGroups[cat]) {
                categoryGroups[cat].push(node);
            } else {
                categoryGroups.other.push(node);
            }
        });

        const totalOuterNodes = nodes.length - (querySet.size || 1);
        const baseRadiusMin = totalOuterNodes > 80 ? 280 : (totalOuterNodes > 30 ? 200 : 140);
        const ringStep = 65;

        Object.keys(categoryGroups).forEach(catKey => {
            const groupNodes = categoryGroups[catKey];
            if (groupNodes.length === 0) return;

            const sector = DEFAULT_SECTORS[catKey] || DEFAULT_SECTORS.other;
            const secStart = sector.startAngle;
            const effectiveSpan = sector.span;

            let currentRing = 0;
            let currentRingIndex = 0;

            groupNodes.forEach((node, idx) => {
                const radius = baseRadiusMin + currentRing * ringStep;
                const minArcLength = 32;
                const arcCapacity = Math.max(3, Math.floor((radius * effectiveSpan) / minArcLength));

                if (currentRingIndex >= arcCapacity) {
                    currentRing += 1;
                    currentRingIndex = 0;
                }

                const ringRadius = baseRadiusMin + currentRing * ringStep + (idx % 2) * 10;
                const ringCap = Math.max(3, Math.floor((ringRadius * effectiveSpan) / minArcLength));
                const angleFraction = (currentRingIndex + 0.5) / ringCap;
                const angle = secStart + angleFraction * effectiveSpan;

                positions[node.id()] = {
                    x: centerPos.x + ringRadius * Math.cos(angle),
                    y: centerPos.y + ringRadius * Math.sin(angle),
                };

                currentRingIndex += 1;
            });
        });

        return positions;
    }

    function applyDandelionLayout(cy, options = {}) {
        if (!cy || typeof cy.nodes !== 'function') return null;

        const {
            queryLoci = [],
            animate = true,
            bloomDuration = 550,
            fit = true,
            padding = 50,
            onComplete = null,
        } = options;

        const positions = computeDandelionPositions(cy, queryLoci);
        const nodes = cy.nodes();

        if (nodes.length === 0) return null;

        cy.batch(() => {
            nodes.forEach(node => {
                const targetPos = positions[node.id()] || { x: 0, y: 0 };
                node.position(targetPos);
                node.style({ opacity: 1 });
            });

            if (typeof cy.edges === 'function') {
                cy.edges().forEach(edge => {
                    edge.style({
                        'curve-style': 'unbundled-bezier',
                        'control-point-distances': '25',
                        'control-point-weights': '0.5',
                    });
                });
            }
        });

        if (fit && typeof cy.fit === 'function') {
            cy.fit(padding);
            if (typeof cy.center === 'function') cy.center();
        }

        if (!animate) {
            if (typeof onComplete === 'function') onComplete();
            return positions;
        }

        const queryList = Array.isArray(queryLoci) ? queryLoci : (queryLoci ? [queryLoci] : []);
        const querySet = new Set(queryList.map(q => (q || '').toLowerCase()));
        const isLargeNetwork = nodes.length > 60;

        if (isLargeNetwork) {
            cy.batch(() => {
                nodes.forEach(node => {
                    const targetPos = positions[node.id()];
                    if (!targetPos) return;
                    const id = (node.id() || '').toLowerCase();
                    const isQuery = node.data('type') === 'query' || querySet.has(id);
                    if (!isQuery) {
                        node.position({
                            x: targetPos.x * 0.35,
                            y: targetPos.y * 0.35,
                        });
                    }
                });
            });

            if (typeof cy.fit === 'function') cy.fit(padding);

            nodes.forEach(node => {
                const targetPos = positions[node.id()];
                if (!targetPos) return;
                const id = (node.id() || '').toLowerCase();
                const isQuery = node.data('type') === 'query' || querySet.has(id);
                if (!isQuery && typeof node.animate === 'function') {
                    node.animate({
                        position: targetPos,
                    }, {
                        duration: bloomDuration,
                        easing: 'ease-out',
                    });
                }
            });
        } else {
            nodes.forEach((node, index) => {
                const targetPos = positions[node.id()];
                if (!targetPos) return;

                const id = (node.id() || '').toLowerCase();
                const isQuery = node.data('type') === 'query' || querySet.has(id);

                if (!isQuery) {
                    node.position({
                        x: targetPos.x * 0.3,
                        y: targetPos.y * 0.3,
                    });

                    const delay = Math.min(160, index * 8);
                    setTimeout(() => {
                        if (typeof node.animate === 'function') {
                            node.animate({
                                position: targetPos,
                            }, {
                                duration: bloomDuration,
                                easing: 'ease-out-back',
                            });
                        } else {
                            node.position(targetPos);
                        }
                    }, delay);
                }
            });
        }

        if (typeof onComplete === 'function') {
            setTimeout(onComplete, bloomDuration + 150);
        }

        return positions;
    }

    function DandelionCytoscapeLayout(options) {
        this.options = Object.assign({
            name: 'dandelion',
            animate: true,
            padding: 50,
            fit: true,
        }, options);
    }

    DandelionCytoscapeLayout.prototype.run = function() {
        const options = this.options;
        const cy = options.cy;
        if (!cy) return this;

        if (typeof cy.trigger === 'function') cy.trigger('layoutstart');

        applyDandelionLayout(cy, {
            queryLoci: options.queryLoci || [],
            animate: options.animate !== false,
            fit: options.fit !== false,
            padding: options.padding || 50,
            onComplete: () => {
                if (typeof cy.trigger === 'function') cy.trigger('layoutstop');
            },
        });

        return this;
    };

    function registerCytoscapeExtension(cytoscapeImpl) {
        if (typeof cytoscapeImpl === 'function') {
            try {
                cytoscapeImpl('layout', 'dandelion', DandelionCytoscapeLayout);
            } catch (err) {
                console.warn('Cytoscape dandelion layout extension registration skipped:', err);
            }
        }
    }

    let ambientTimer = null;
    let breezePhase = 0;

    function startAmbientBreeze(cy, amplitude = 2.5) {
        stopAmbientBreeze();
        if (!cy || typeof cy.nodes !== 'function') return;
        if (cy.nodes().length > 80) return;

        function step() {
            breezePhase += 0.04;
            const nodes = cy.nodes();
            nodes.forEach((node, idx) => {
                if (node.data('type') === 'query') return;
                const pos = node.position();
                const offsetX = Math.sin(breezePhase + idx * 0.5) * amplitude * 0.06;
                const offsetY = Math.cos(breezePhase + idx * 0.3) * amplitude * 0.06;
                node.position({
                    x: pos.x + offsetX,
                    y: pos.y + offsetY,
                });
            });
            ambientTimer = requestAnimationFrame(step);
        }
        ambientTimer = requestAnimationFrame(step);
    }

    function stopAmbientBreeze() {
        if (ambientTimer) {
            cancelAnimationFrame(ambientTimer);
            ambientTimer = null;
        }
    }

    if (typeof window !== 'undefined' && window.cytoscape) {
        registerCytoscapeExtension(window.cytoscape);
    }

    return {
        LAYOUT_VERSION,
        DEFAULT_SECTORS,
        getNodeCategory,
        computeDandelionPositions,
        applyDandelionLayout,
        DandelionCytoscapeLayout,
        registerCytoscapeExtension,
        startAmbientBreeze,
        stopAmbientBreeze,
    };
});
