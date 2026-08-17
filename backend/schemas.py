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
    databaseLinks: Optional[Dict] = None
    metaboliteLinks: Optional[Dict[str, Dict]] = None

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
    isEssential: Optional[bool] = None
    essentialGenes: Optional[List[str]] = None
    hasGlobalRegulator: Optional[bool] = None
    globalRegulators: Optional[List[str]] = None

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
    tracked_fluxes: Optional[Dict[str, List[float]]] = None
    warnings: List[str] = []

class RECFBARequest(BaseModel):
    tfPerturbations: Dict[str, str]  # e.g., {"sigH": "knockout"}
    proteinPoolLimit: float = 0.129
    temperature: float = 30.0
    initialGlucose: float = 100.0
    initialBiomass: float = 0.1
    timeSteps: int = 24

class RECFBAResponse(BaseModel):
    status: str
    time: List[float]
    growth_rate: List[float]
    glutamate_export: List[float]
    glucose_uptake: List[float]
    glucose_concentration: List[float]
    biomass_concentration: List[float]
    tracked_fluxes: Optional[Dict[str, List[float]]] = None
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

class CascadeEdgeSchema(BaseModel):
    from_node: str
    to_node: str
    role: str = "A"
    score: float = 1.0

class CascadePathSchema(BaseModel):
    nodes: List[str]
    length: int
    edges: List[Dict]

class GraphCascadeResponse(BaseModel):
    source: str
    target: str
    paths: List[CascadePathSchema] = []

class GraphMotifResponse(BaseModel):
    motif_type: str
    count: int
    items: List[Dict] = []


# ── New Standardized Response Schemas ───────────────────────────────────────

class HTTPError(BaseModel):
    detail: str


class GeneCoordinatesResponse(BaseModel):
    locus_tag: str
    gene_name: Optional[str] = None
    start_pos: Optional[int] = None
    end_pos: Optional[int] = None
    strand: Optional[str] = None
    gene_length: Optional[int] = None
    tss_position: Optional[int] = None
    promoter_70bp: Optional[str] = None
    contig: Optional[str] = None


class GeneNeighborhoodItem(BaseModel):
    locus_tag: str
    gene_name: Optional[str] = None
    start_pos: Optional[int] = None
    end_pos: Optional[int] = None
    strand: Optional[str] = None
    product: Optional[str] = None


class GenomicNeighborhoodResponse(BaseModel):
    center_gene: str
    window_bp: int
    count: int
    genes: List[GeneNeighborhoodItem] = []


class ChIPSeqPeakItem(BaseModel):
    tf_locus: Optional[str] = None
    tf_name: Optional[str] = None
    target_locus: Optional[str] = None
    peak_start: Optional[int] = None
    peak_end: Optional[int] = None
    enrichment_fold: Optional[float] = None
    p_value: Optional[float] = None
    q_value: Optional[float] = None
    condition: Optional[str] = None
    source_dataset: Optional[str] = None


class ChIPSeqPeaksResponse(BaseModel):
    query: str
    as_target_count: int = 0
    as_target_peaks: List[ChIPSeqPeakItem] = []
    as_tf_count: int = 0
    as_tf_peaks: List[ChIPSeqPeakItem] = []
    is_public_deployment: bool = False
    message: Optional[str] = None


class ModuleEvidenceItem(BaseModel):
    module_run_id: str
    mean_score: Optional[float] = None
    condition_count: int = 0
    supported_count: int = 0
    significant_context_count: int = 0


class InterventionTargetSchema(BaseModel):
    target_locus: str
    target_name: Optional[str] = None
    product: Optional[str] = None
    evidence_grade: str = "D"
    module_count: int = 0
    modules: List[str] = []
    evidence_score: float = 0.0
    systems_impact_score: float = 0.0
    risk_score: float = 0.0
    engineering_tractability_score: float = 0.0
    priority_score: float = 0.0
    strategy_class: str = "context_specific_candidate"
    essentiality_status: str = "non_essential"
    proteomics_detected: bool = False
    pathway_count: int = 0
    module_evidence: List[ModuleEvidenceItem] = []


class InterventionTargetsResponse(BaseModel):
    total: int = 0
    limit: int = 150
    targets: List[InterventionTargetSchema] = []


class IModulonSummarySchema(BaseModel):
    imodulon_id: str
    name: str
    category: str
    regulator: Optional[str] = None
    explained_variance: float = 0.0
    gene_count: int = 0


class IModulonsListResponse(BaseModel):
    total: int = 0
    imodulons: List[IModulonSummarySchema] = []


class IModulonDetailResponse(BaseModel):
    imodulon_id: str
    name: str
    category: str
    description: Optional[str] = None
    explained_variance: float = 0.0
    gene_count: int = 0
    regulator: Optional[str] = None
    precision: Optional[float] = None
    recall: Optional[float] = None
    f1_score: Optional[float] = None
    top_genes: List[Dict] = []
    pathway_enrichment: List[Dict] = []


class RAGQueryRequest(BaseModel):
    query: str
    top_k: int = 5
    include_raw: bool = False


class RAGSourceSchema(BaseModel):
    title: str
    authors: Optional[str] = None
    journal: Optional[str] = None
    year: Optional[int] = None
    doi: Optional[str] = None
    pmid: Optional[str] = None
    score: float = 0.0
    snippet: Optional[str] = None


class RAGQueryResponse(BaseModel):
    query: str
    answer: str
    sources: List[RAGSourceSchema] = []

