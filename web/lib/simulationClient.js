(function (global) {
    const BASE_URL = 'http://localhost:8001';

    async function getModelStatus() {
        try {
            const res = await fetch(`${BASE_URL}/api/model/status`);
            if (!res.ok) throw new Error(`HTTP error: ${res.status}`);
            return await res.json();
        } catch (err) {
            return {
                loaded: false,
                reaction_count: 0,
                gene_count: 0,
                metabolite_count: 0,
                error: err.message || "FastAPI backend offline"
            };
        }
    }

    async function searchReactions(query) {
        try {
            const res = await fetch(`${BASE_URL}/api/model/reactions/search?q=${encodeURIComponent(query)}`);
            if (!res.ok) throw new Error(`HTTP error: ${res.status}`);
            return await res.json();
        } catch (err) {
            return { query, matches: [], error: err.message || "FastAPI backend offline" };
        }
    }

    async function runBaselineSimulation() {
        try {
            const res = await fetch(`${BASE_URL}/api/simulation/baseline`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' }
            });
            if (!res.ok) throw new Error(`HTTP error: ${res.status}`);
            return await res.json();
        } catch (err) {
            return { status: "error", error: err.message || "FastAPI backend offline" };
        }
    }

    async function runGeneKnockout(geneId, objective, trackReactionIds) {
        try {
            const res = await fetch(`${BASE_URL}/api/simulation/gene-knockout`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ geneId, objective, trackReactionIds })
            });
            if (!res.ok) throw new Error(`HTTP error: ${res.status}`);
            return await res.json();
        } catch (err) {
            return {
                status: "error",
                baselineObjective: 0,
                perturbedObjective: 0,
                objectiveChange: 0,
                objectiveChangePercent: 0,
                trackedFluxes: [],
                error: err.message || "FastAPI backend offline"
            };
        }
    }

    async function runGeneSetKnockout(geneIds, objective, trackReactionIds) {
        try {
            const res = await fetch(`${BASE_URL}/api/simulation/gene-set-knockout`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ geneIds, objective, trackReactionIds })
            });
            if (!res.ok) throw new Error(`HTTP error: ${res.status}`);
            return await res.json();
        } catch (err) {
            return {
                status: "error",
                baselineObjective: 0,
                perturbedObjective: 0,
                objectiveChange: 0,
                objectiveChangePercent: 0,
                trackedFluxes: [],
                error: err.message || "FastAPI backend offline"
            };
        }
    }

    async function runTFPerturbation(tfId, targetGeneIds, objective, trackReactionIds) {
        try {
            const res = await fetch(`${BASE_URL}/api/simulation/tf-perturbation`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ tfId, targetGeneIds, mode: "knockout", objective, trackReactionIds })
            });
            if (!res.ok) throw new Error(`HTTP error: ${res.status}`);
            return await res.json();
        } catch (err) {
            return {
                tfId,
                status: "error",
                targetGeneCount: targetGeneIds.length,
                mappedGeneCount: 0,
                missingGenes: targetGeneIds,
                baselineObjective: 0,
                perturbedObjective: 0,
                objectiveChange: 0,
                objectiveChangePercent: 0,
                trackedFluxes: [],
                error: err.message || "FastAPI backend offline"
            };
        }
    }

    async function runFluxVariabilityAnalysis(mode, geneId, targetGeneIds, objective, trackReactionIds, fractionOfOptimum = 0.95) {
        try {
            const res = await fetch(`${BASE_URL}/api/simulation/flux-variability`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ mode, geneId, targetGeneIds, objective, trackReactionIds, fractionOfOptimum })
            });
            if (!res.ok) throw new Error(`HTTP error: ${res.status}`);
            return await res.json();
        } catch (err) {
            return {
                status: "error",
                fractionOfOptimum,
                fvaRanges: [],
                error: err.message || "FastAPI backend offline"
            };
        }
    }

    global.simulationClient = {
        getModelStatus,
        searchReactions,
        runBaselineSimulation,
        runGeneKnockout,
        runGeneSetKnockout,
        runTFPerturbation,
        runFluxVariabilityAnalysis
    };
})(window);
