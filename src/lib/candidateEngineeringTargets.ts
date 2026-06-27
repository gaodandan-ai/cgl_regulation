import { RegulatoryGraph, normalizeGeneId } from "./metabolicModelAdapter";
import { getTargetGenesForTF } from "./regulationMetabolismBridge";
import {
  TFMetabolicImpactRank,
  getAllTFNodes,
  rankTFsByMetabolicImpact
} from "./tfMetabolicImpactRanking";

export type EngineeringTargetCandidate = {
  tfId: string;
  tfLabel: string;
  candidateScore: number;
  totalTargetGenes: number;
  mappedTargetGenes: number;
  totalReactions: number;
  totalPathways: number;
  keyPathways: string[];
  regulatedKeyGenes: string[];
  averageConfidence: number;
  regulationProfile: {
    activationCount: number;
    repressionCount: number;
    predictedCount: number;
    unknownCount: number;
  };
  recommendationLevel: "high" | "medium" | "low";
  rationale: string;
};

export const ENGINEERING_PATHWAY_KEYWORDS = [
  "glutamate",
  "glutamic acid",
  "amino acid",
  "lysine",
  "arginine",
  "tca",
  "citric acid cycle",
  "central carbon",
  "glycolysis",
  "pyruvate",
  "acetyl-coa",
  "transport",
  "nitrogen",
  "carbon metabolism"
];

function clamp(value: number, min = 0, max = 1): number {
  if (!Number.isFinite(value)) return min;
  return Math.max(min, Math.min(max, value));
}

function toArray(value: any): any[] {
  if (!value) return [];
  if (Array.isArray(value)) return value;
  if (typeof value.toArray === "function") return value.toArray();
  if (typeof value.forEach === "function") {
    const items: any[] = [];
    value.forEach((item: any) => items.push(item));
    return items;
  }
  return [];
}

function callOrRead(value: any, key: string): any {
  if (!value) return undefined;
  if (typeof value.data === "function") return value.data(key);
  if (value.data && key in value.data) return value.data[key];
  if (key in value) return value[key];
  return undefined;
}

function collectEdges(graph: any): any[] {
  if (!graph) return [];
  if (typeof graph.edges === "function") return toArray(graph.edges());
  if (graph.edges) return toArray(graph.edges);
  if (graph.elements?.edges) return toArray(graph.elements.edges);
  if (Array.isArray(graph.elements)) return graph.elements.filter((item: any) => item.group === "edges" || item.data?.source);
  return [];
}

function edgeSource(edge: any): string {
  return String(callOrRead(edge, "source") || "");
}

function edgeRegulation(edge: any): string {
  const raw = callOrRead(edge, "regulation") || callOrRead(edge, "type") || callOrRead(edge, "role") || "unknown";
  if (typeof raw === "object") return "unknown";
  const lower = String(raw).toLowerCase();
  if (lower === "a" || lower.includes("activ")) return "activation";
  if (lower === "r" || lower.includes("repress") || lower.includes("inhibit")) return "repression";
  if (lower.includes("predict")) return "predicted";
  return lower || "unknown";
}

function pathwayMatchesEngineeringKeyword(pathway: string, keyword?: string): boolean {
  const text = String(pathway || "").toLowerCase();
  if (!text) return false;
  if (keyword && !text.includes(keyword.toLowerCase())) return false;
  return ENGINEERING_PATHWAY_KEYWORDS.some(key => text.includes(key));
}

function recommendationLevel(score: number): "high" | "medium" | "low" {
  if (score >= 0.75) return "high";
  if (score >= 0.45) return "medium";
  return "low";
}

export function calculateEngineeringCandidateScore(input: {
  mappedTargetGenes: number;
  totalReactions: number;
  totalPathways: number;
  averageConfidence: number;
  keyPathwayHits: number;
  regulatedKeyGeneCount: number;
  repressionCount?: number;
  activationCount?: number;
}): number {
  const mappedGenes = clamp((Number(input?.mappedTargetGenes) || 0) / 30);
  const reactions = clamp((Number(input?.totalReactions) || 0) / 60);
  const pathways = clamp((Number(input?.totalPathways) || 0) / 15);
  const confidence = clamp(Number(input?.averageConfidence) || 0);
  const keyPathways = clamp((Number(input?.keyPathwayHits) || 0) / 8);
  const keyGenes = clamp((Number(input?.regulatedKeyGeneCount) || 0) / 20);
  const base = 0.25 * mappedGenes + 0.20 * reactions + 0.15 * pathways + 0.20 * confidence + 0.15 * keyPathways + 0.05 * keyGenes;
  return Number(clamp(base).toFixed(3));
}

export function getRegulationProfile(edges: any[]): {
  activationCount: number;
  repressionCount: number;
  predictedCount: number;
  unknownCount: number;
} {
  const profile = { activationCount: 0, repressionCount: 0, predictedCount: 0, unknownCount: 0 };
  for (const edge of edges || []) {
    const type = edgeRegulation(edge);
    if (type === "activation") profile.activationCount += 1;
    else if (type === "repression") profile.repressionCount += 1;
    else if (type === "predicted") profile.predictedCount += 1;
    else profile.unknownCount += 1;
  }
  return profile;
}

export function generateEngineeringTargetRationale(
  candidate: Omit<EngineeringTargetCandidate, "rationale">
): string {
  if (!candidate || candidate.candidateScore < 0.45) {
    return "This TF has limited evidence as an engineering target based on the current regulatory and metabolic mappings.";
  }
  const levelText = candidate.recommendationLevel === "high" ? "high-priority" : "medium-priority";
  const pathwayText = candidate.keyPathways.length > 0
    ? `, including ${candidate.keyPathways.slice(0, 3).join(" and ")}`
    : "";
  const keyGeneText = candidate.regulatedKeyGenes.length > 0
    ? ` Its targets include ${candidate.regulatedKeyGenes.slice(0, 5).join(", ")}, which are associated with key metabolic modules.`
    : "";
  return `${candidate.tfLabel || candidate.tfId} is a ${levelText} candidate engineering regulator. It regulates ${candidate.mappedTargetGenes} genes mapped to ${candidate.totalReactions} metabolic reactions across ${candidate.totalPathways} pathways${pathwayText}.${keyGeneText} This TF may influence metabolic phenotype and could be prioritized for further experimental or simulation-based evaluation.`;
}

function candidateFromRank(rank: TFMetabolicImpactRank & { pathwaySummary?: any[] }, graph: RegulatoryGraph): EngineeringTargetCandidate {
  const allEdges = collectEdges(graph);
  const tfEdges = allEdges.filter(edge => normalizeGeneId(edgeSource(edge)) === normalizeGeneId(rank.tfId));
  const regulationProfile = getRegulationProfile(tfEdges);
  const keyPathways = (rank.pathwaySummary || [])
    .filter(pathway => pathwayMatchesEngineeringKeyword(`${pathway.pathwayName || ""} ${pathway.pathwayId || ""}`))
    .map(pathway => pathway.pathwayName || pathway.pathwayId)
    .filter(Boolean);
  const regulatedKeyGenes = Array.from(new Set((rank.pathwaySummary || [])
    .filter(pathway => pathwayMatchesEngineeringKeyword(`${pathway.pathwayName || ""} ${pathway.pathwayId || ""}`))
    .flatMap(pathway => pathway.genes || []))).sort();
  const candidateScore = calculateEngineeringCandidateScore({
    mappedTargetGenes: rank.mappedTargetGenes,
    totalReactions: rank.totalReactions,
    totalPathways: rank.totalPathways,
    averageConfidence: rank.averageConfidence,
    keyPathwayHits: keyPathways.length,
    regulatedKeyGeneCount: regulatedKeyGenes.length,
    repressionCount: regulationProfile.repressionCount,
    activationCount: regulationProfile.activationCount
  });
  const partial = {
    tfId: rank.tfId,
    tfLabel: rank.tfLabel,
    candidateScore,
    totalTargetGenes: rank.totalTargetGenes,
    mappedTargetGenes: rank.mappedTargetGenes,
    totalReactions: rank.totalReactions,
    totalPathways: rank.totalPathways,
    keyPathways,
    regulatedKeyGenes,
    averageConfidence: rank.averageConfidence,
    regulationProfile,
    recommendationLevel: recommendationLevel(candidateScore)
  };
  return {
    ...partial,
    rationale: generateEngineeringTargetRationale(partial)
  };
}

export function findEngineeringTargetCandidates(
  graph: RegulatoryGraph,
  options?: {
    limit?: number;
    pathwayKeywordFilter?: string;
    minCandidateScore?: number;
    includeLowConfidence?: boolean;
  }
): EngineeringTargetCandidate[] {
  const limit = options?.limit ?? 20;
  const minCandidateScore = options?.minCandidateScore ?? 0;
  const filter = (options?.pathwayKeywordFilter || "").trim().toLowerCase();
  getAllTFNodes(graph); // Keeps the sync API resilient if graph has no TF nodes.

  return rankTFsByMetabolicImpact(graph, { limit: 500, includeZeroImpact: true })
    .map(rank => candidateFromRank(rank as TFMetabolicImpactRank & { pathwaySummary?: any[] }, graph))
    .filter(candidate => candidate.candidateScore >= minCandidateScore)
    .filter(candidate => options?.includeLowConfidence || candidate.averageConfidence > 0 || candidate.mappedTargetGenes > 0)
    .filter(candidate => !filter || candidate.keyPathways.some(pathway => pathway.toLowerCase().includes(filter)))
    .sort((a, b) => b.candidateScore - a.candidateScore || b.mappedTargetGenes - a.mappedTargetGenes || a.tfLabel.localeCompare(b.tfLabel))
    .slice(0, limit);
}
