from pydantic import BaseModel
from typing import List, Optional, Dict

class ModelStatusResponse(BaseModel):
    loaded: bool
    model_id: Optional[str] = None
    reaction_count: int = 0
    gene_count: int = 0
    metabolite_count: int = 0
    error: Optional[str] = None

class ReactionMatchSchema(BaseModel):
    reactionId: str
    name: str
    equation: str
    lowerBound: float
    upperBound: float
    metabolites: List[str]

class ReactionSearchResponse(BaseModel):
    query: str
    matches: List[ReactionMatchSchema]

class ObjectiveSchema(BaseModel):
    objectiveType: str  # "biomass" or "reaction"
    reactionId: Optional[str] = None

class ObjectiveResponseSchema(BaseModel):
    objectiveType: str
    reactionId: Optional[str] = None
    label: str

class TrackedFluxSchema(BaseModel):
    reactionId: str
    baselineFlux: Optional[float] = None
    perturbedFlux: Optional[float] = None
    fluxChange: Optional[float] = None
    fluxChangePercent: Optional[float] = None

class BaselineSimulationResponse(BaseModel):
    status: str
    objective_value: Optional[float] = None
    objective_expression: Optional[str] = None
    warnings: List[str] = []

class GeneKnockoutRequest(BaseModel):
    geneId: str
    objective: Optional[ObjectiveSchema] = None
    trackReactionIds: Optional[List[str]] = None
    method: str = "fba"

class GeneKnockoutResponse(BaseModel):
    status: str
    objective: ObjectiveResponseSchema
    baselineObjective: float
    perturbedObjective: float
    objectiveChange: float
    objectiveChangePercent: float
    trackedFluxes: List[TrackedFluxSchema] = []
    warnings: List[str] = []

class GeneSetKnockoutRequest(BaseModel):
    geneIds: List[str]
    objective: Optional[ObjectiveSchema] = None
    trackReactionIds: Optional[List[str]] = None
    method: str = "fba"

class GeneSetKnockoutResponse(BaseModel):
    status: str
    objective: ObjectiveResponseSchema
    baselineObjective: float
    perturbedObjective: float
    objectiveChange: float
    objectiveChangePercent: float
    trackedFluxes: List[TrackedFluxSchema] = []
    missingGenes: List[str] = []
    warnings: List[str] = []

class TFPerturbationRequest(BaseModel):
    tfId: str
    targetGeneIds: List[str]
    mode: str = "knockout"
    objective: Optional[ObjectiveSchema] = None
    trackReactionIds: Optional[List[str]] = None
    method: str = "fba"

class TFPerturbationResponse(BaseModel):
    tfId: str
    status: str
    targetGeneCount: int
    mappedGeneCount: int
    missingGenes: List[str]
    objective: ObjectiveResponseSchema
    baselineObjective: float
    perturbedObjective: float
    objectiveChange: float
    objectiveChangePercent: float
    trackedFluxes: List[TrackedFluxSchema] = []
    warnings: List[str] = []

class GlutamateCandidateSchema(BaseModel):
    reactionId: str
    name: str
    equation: str
    lowerBound: float
    upperBound: float
    classification: str
    confidence: str
    reason: str

class GlutamateCandidatesResponse(BaseModel):
    candidates: List[GlutamateCandidateSchema]
    warnings: List[str] = []

class FVARangeSchema(BaseModel):
    reactionId: str
    baselineMin: float
    baselineMax: float
    perturbedMin: float
    perturbedMax: float

class FVARequest(BaseModel):
    geneId: Optional[str] = None
    targetGeneIds: Optional[List[str]] = None
    mode: str = "baseline"
    objective: Optional[ObjectiveSchema] = None
    trackReactionIds: Optional[List[str]] = None
    fractionOfOptimum: float = 0.95

class FVAResponse(BaseModel):
    status: str
    fractionOfOptimum: float
    fvaRanges: List[FVARangeSchema]
    warnings: List[str] = []

class RFBARequest(BaseModel):
    tfPerturbations: Dict[str, str]  # e.g., {"sigH": "knockout"}
    initialGlucose: float = 100.0
    initialBiomass: float = 0.1
    timeSteps: int = 24

class RFBAResponse(BaseModel):
    status: str
    time: List[float]
    growth_rate: List[float]
    glutamate_export: List[float]
    glucose_uptake: List[float]
    glucose_concentration: List[float]
    biomass_concentration: List[float]
    warnings: List[str] = []

class ECFBARequest(BaseModel):
    proteinPoolLimit: float = 0.129
    enzymePerturbations: Dict[str, float]  # e.g., {"gdh": 1.0, "lysC": 1.0}
    targetProduct: str = "growth"  # "growth", "glutamate", "lysine"
    temperature: float = 30.0
    calibrateTimepoint: Optional[str] = None

class ECFBABottleneckSchema(BaseModel):
    reaction_id: str
    reaction_name: str
    genes: str
    flux: float
    usage: float
    shadow_price: float

class ECFBAResponse(BaseModel):
    status: str
    flux: float
    poolLimit: float
    poolUsage: float
    warnings: List[str] = []
    calibratedPerturbations: Optional[Dict[str, float]] = None
    bottlenecks: List[ECFBABottleneckSchema] = []

class MFAComparisonItem(BaseModel):
    reaction_id: str
    reaction_name: str
    pathway: str
    mfa_flux: float
    mfa_std: float
    sim_flux: float
    deviation_pct: float
    matched_model_id: Optional[str] = None
    reference: str

class MFAComparisonResponse(BaseModel):
    status: str
    items: List[MFAComparisonItem] = []
    pearson_r: float = 0.0
    rmse: float = 0.0
    mean_deviation_pct: float = 0.0
    warnings: List[str] = []

class PathwayReactionsRequest(BaseModel):
    reactionIds: List[str]
