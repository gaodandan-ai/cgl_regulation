export type ModelStatus = {
  loaded: boolean;
  model_id?: string;
  reaction_count: number;
  gene_count: number;
  metabolite_count: number;
  error?: string;
};

export type ObjectiveConfig = {
  objectiveType: "biomass" | "reaction";
  reactionId?: string | null;
};

export type TrackedFluxResult = {
  reactionId: string;
  baselineFlux: number | null;
  perturbedFlux: number | null;
  fluxChange: number | null;
  fluxChangePercent: number | null;
};

export type SimulationResult = {
  status: string;
  objective?: {
    objectiveType: string;
    reactionId?: string | null;
    label: string;
  };
  baselineObjective: number;
  perturbedObjective: number;
  objectiveChange: number;
  objectiveChangePercent: number;
  trackedFluxes?: TrackedFluxResult[];
  missingGenes?: string[];
  warnings?: string[];
  error?: string;
};

export type TFPerturbationResult = {
  tfId: string;
  status: string;
  targetGeneCount: number;
  mappedGeneCount: number;
  missingGenes: string[];
  objective?: {
    objectiveType: string;
    reactionId?: string | null;
    label: string;
  };
  baselineObjective: number;
  perturbedObjective: number;
  objectiveChange: number;
  objectiveChangePercent: number;
  trackedFluxes?: TrackedFluxResult[];
  warnings?: string[];
  error?: string;
};

const BASE_URL = (typeof window !== 'undefined' && window.location.protocol === 'file:') ? 'http://localhost:8000' : '';

export async function getModelStatus(): Promise<ModelStatus> {
  try {
    const res = await fetch(`${BASE_URL}/api/model/status`);
    if (!res.ok) throw new Error(`HTTP error: ${res.status}`);
    return await res.json();
  } catch (err: any) {
    return {
      loaded: false,
      reaction_count: 0,
      gene_count: 0,
      metabolite_count: 0,
      error: err.message || "FastAPI backend offline"
    };
  }
}

export async function searchReactions(query: string): Promise<any> {
  try {
    const res = await fetch(`${BASE_URL}/api/model/reactions/search?q=${encodeURIComponent(query)}`);
    if (!res.ok) throw new Error(`HTTP error: ${res.status}`);
    return await res.json();
  } catch (err: any) {
    return { query, matches: [], error: err.message || "FastAPI backend offline" };
  }
}

export async function runBaselineSimulation(): Promise<any> {
  try {
    const res = await fetch(`${BASE_URL}/api/simulation/baseline`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' }
    });
    if (!res.ok) throw new Error(`HTTP error: ${res.status}`);
    return await res.json();
  } catch (err: any) {
    return { status: "error", error: err.message || "FastAPI backend offline" };
  }
}

export async function runGeneKnockout(
  geneId: string, 
  objective?: ObjectiveConfig, 
  trackReactionIds?: string[]
): Promise<SimulationResult> {
  try {
    const res = await fetch(`${BASE_URL}/api/simulation/gene-knockout`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ geneId, objective, trackReactionIds })
    });
    if (!res.ok) throw new Error(`HTTP error: ${res.status}`);
    return await res.json();
  } catch (err: any) {
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

export async function runGeneSetKnockout(
  geneIds: string[], 
  objective?: ObjectiveConfig, 
  trackReactionIds?: string[]
): Promise<SimulationResult> {
  try {
    const res = await fetch(`${BASE_URL}/api/simulation/gene-set-knockout`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ geneIds, objective, trackReactionIds })
    });
    if (!res.ok) throw new Error(`HTTP error: ${res.status}`);
    return await res.json();
  } catch (err: any) {
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

export async function runTFPerturbation(
  tfId: string, 
  targetGeneIds: string[], 
  objective?: ObjectiveConfig, 
  trackReactionIds?: string[]
): Promise<TFPerturbationResult> {
  try {
    const res = await fetch(`${BASE_URL}/api/simulation/tf-perturbation`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ tfId, targetGeneIds, mode: "knockout", objective, trackReactionIds })
    });
    if (!res.ok) throw new Error(`HTTP error: ${res.status}`);
    return await res.json();
  } catch (err: any) {
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

export type RFBAResult = {
  status: string;
  time: number[];
  growth_rate: number[];
  glutamate_export: number[];
  glucose_uptake: number[];
  glucose_concentration: number[];
  biomass_concentration: number[];
  warnings?: string[];
  error?: string;
};

export async function runDynamicRFBA(
  tfPerturbations: Record<string, string>,
  initialGlucose: number,
  initialBiomass: number,
  timeSteps: number = 24
): Promise<RFBAResult> {
  try {
    const res = await fetch(`${BASE_URL}/api/simulation/rfba`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ tfPerturbations, initialGlucose, initialBiomass, timeSteps })
    });
    if (!res.ok) throw new Error(`HTTP error: ${res.status}`);
    return await res.json();
  } catch (err: any) {
    return {
      status: "error",
      time: [],
      growth_rate: [],
      glutamate_export: [],
      glucose_uptake: [],
      glucose_concentration: [],
      biomass_concentration: [],
      error: err.message || "FastAPI backend offline"
    };
  }
}

export type ECFBAResult = {
  status: string;
  flux: number;
  poolLimit: number;
  poolUsage: number;
  warnings?: string[];
  error?: string;
  calibratedPerturbations?: Record<string, number>;
};

export async function runECFBA(
  proteinPoolLimit: number,
  enzymePerturbations: Record<string, number>,
  targetProduct: string,
  temperature: number = 30.0,
  calibrateTimepoint: string | null = null
): Promise<ECFBAResult> {
  try {
    const res = await fetch(`${BASE_URL}/api/simulation/ecfba`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ proteinPoolLimit, enzymePerturbations, targetProduct, temperature, calibrateTimepoint })
    });
    if (!res.ok) throw new Error(`HTTP error: ${res.status}`);
    return await res.json();
  } catch (err: any) {
    return {
      status: "error",
      flux: 0.0,
      poolLimit: proteinPoolLimit,
      poolUsage: 0.0,
      error: err.message || "FastAPI backend offline"
    };
  }
}
