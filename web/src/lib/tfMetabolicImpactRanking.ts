import {
  RegulatoryGraph,
  normalizeGeneId
} from "./metabolicModelAdapter";
import {
  generateMetabolicImpactExplanation,
  getMetabolicImpactForTF,
  getTargetGenesForTF
} from "./regulationMetabolismBridge";

export type RegulatoryNode = {
  id: string;
  label?: string;
  type?: string;
  data?: {
    id?: string;
    label?: string;
    name?: string;
    type?: string;
  };
};

export type TFMetabolicImpactRank = {
  tfId: string;
  tfLabel: string;
  totalTargetGenes: number;
  mappedTargetGenes: number;
  totalReactions: number;
  totalPathways: number;
  averageConfidence: number;
  keyPathways: string[];
  impactScore: number;
  explanation: string;
};

export const KEY_CGL_PATHWAY_KEYWORDS = [
  "glutamate",
  "glutamic acid",
  "amino acid",
  "lysine",
  "arginine",
  "tca",
  "citric acid cycle",
  "central carbon",
  "glycolysis",
  "transport"
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

function collectNodes(graph: any): any[] {
  if (!graph) return [];
  if (typeof graph.nodes === "function") return toArray(graph.nodes());
  if ((graph as any).nodes) return toArray((graph as any).nodes);
  if ((graph as any).elements?.nodes) return toArray((graph as any).elements.nodes);
  if (Array.isArray((graph as any).elements)) {
    return (graph as any).elements.filter((item: any) => item.group === "nodes" || item.data?.id);
  }
  return [];
}

function collectEdges(graph: any): any[] {
  if (!graph) return [];
  if (typeof graph.edges === "function") return toArray(graph.edges());
  if ((graph as any).edges) return toArray((graph as any).edges);
  if ((graph as any).elements?.edges) return toArray((graph as any).elements.edges);
  if (Array.isArray((graph as any).elements)) {
    return (graph as any).elements.filter((item: any) => item.group === "edges" || item.data?.source);
  }
  return [];
}

function nodeId(node: any): string {
  return String(callOrRead(node, "id") || "");
}

function nodeLabel(node: any): string {
  return String(callOrRead(node, "label") || callOrRead(node, "name") || nodeId(node));
}

function isTFNode(node: any): boolean {
  const type = String(callOrRead(node, "type") || "").toLowerCase();
  return ["tf", "transcription_factor", "regulator"].includes(type);
}

function edgeSource(edge: any): string {
  return String(callOrRead(edge, "source") || "");
}

function edgeConfidence(edge: any): number {
  const raw = callOrRead(edge, "confidence") ?? callOrRead(edge, "confidenceScore");
  const parsed = Number(raw);
  return Number.isFinite(parsed) ? clamp(parsed) : 0;
}

function keyPathwayNames(pathwayNames: string[]): string[] {
  const found = new Map<string, string>();
  for (const pathwayName of pathwayNames) {
    const lower = String(pathwayName || "").toLowerCase();
    if (!lower) continue;
    for (const keyword of KEY_CGL_PATHWAY_KEYWORDS) {
      if (lower.includes(keyword)) {
        found.set(pathwayName, pathwayName);
        break;
      }
    }
  }
  return Array.from(found.values());
}

export function getAllTFNodes(graph: RegulatoryGraph): RegulatoryNode[] {
  const seen = new Set<string>();
  const tfs: RegulatoryNode[] = [];

  for (const node of collectNodes(graph)) {
    const id = nodeId(node);
    if (!id || seen.has(normalizeGeneId(id)) || !isTFNode(node)) continue;
    seen.add(normalizeGeneId(id));
    tfs.push({
      id,
      label: nodeLabel(node),
      type: callOrRead(node, "type"),
      data: (node as any).data
    });
  }

  return tfs;
}

export function calculateTFMetabolicImpactScore(input: {
  mappedTargetGenes: number;
  totalReactions: number;
  totalPathways: number;
  averageConfidence: number;
  keyPathwayHits?: number;
}): number {
  const mappedGenes = clamp((Number(input?.mappedTargetGenes) || 0) / 30);
  const reactions = clamp((Number(input?.totalReactions) || 0) / 60);
  const pathways = clamp((Number(input?.totalPathways) || 0) / 15);
  const confidence = clamp(Number(input?.averageConfidence) || 0);
  const keyHits = Math.max(0, Number(input?.keyPathwayHits) || 0);

  const base = 0.35 * mappedGenes + 0.25 * reactions + 0.20 * pathways + 0.20 * confidence;
  const bonus = 0.05 * keyHits;
  return Number(clamp(base + bonus).toFixed(3));
}

export function generateTFImpactRankExplanation(rank: TFMetabolicImpactRank): string {
  if (!rank || rank.impactScore < 0.15 || rank.mappedTargetGenes === 0) {
    return "This TF has limited metabolic model coverage based on current mappings.";
  }

  const pathwayText = rank.keyPathways.length > 0
    ? `, including ${rank.keyPathways.slice(0, 3).join(" and ")}`
    : "";
  const level = rank.impactScore >= 0.7 ? "high" : rank.impactScore >= 0.4 ? "moderate" : "limited";

  return `${rank.tfLabel || rank.tfId} has a ${level} predicted metabolic impact. It regulates ${rank.mappedTargetGenes} genes mapped to ${rank.totalReactions} metabolic reactions across ${rank.totalPathways} pathways${pathwayText}.`;
}

export function rankTFsByMetabolicImpact(
  graph: RegulatoryGraph,
  options?: {
    limit?: number;
    includeZeroImpact?: boolean;
  }
): TFMetabolicImpactRank[] {
  const limit = options?.limit ?? 20;
  const includeZeroImpact = options?.includeZeroImpact ?? false;
  const edges = collectEdges(graph);

  const ranks = getAllTFNodes(graph).map(tf => {
    const impact = getMetabolicImpactForTF(graph, tf.id);
    const targets = getTargetGenesForTF(graph, tf.id);
    const outgoingConfidences = edges
      .filter(edge => normalizeGeneId(edgeSource(edge)) === normalizeGeneId(tf.id))
      .map(edgeConfidence);
    const averageConfidence = outgoingConfidences.length > 0
      ? outgoingConfidences.reduce((sum, value) => sum + value, 0) / outgoingConfidences.length
      : 0;
    const pathwayNames = impact.pathwaySummary.map(pathway => pathway.pathwayName);
    const keyPathways = keyPathwayNames(pathwayNames);
    const impactScore = calculateTFMetabolicImpactScore({
      mappedTargetGenes: impact.mappedTargetGenes,
      totalReactions: impact.totalReactions,
      totalPathways: impact.totalPathways,
      averageConfidence,
      keyPathwayHits: keyPathways.length
    });
    const rank: TFMetabolicImpactRank = {
      tfId: tf.id,
      tfLabel: tf.label || tf.id,
      totalTargetGenes: impact.totalTargetGenes || targets.length,
      mappedTargetGenes: impact.mappedTargetGenes,
      totalReactions: impact.totalReactions,
      totalPathways: impact.totalPathways,
      averageConfidence: Number(averageConfidence.toFixed(3)),
      keyPathways,
      impactScore,
      explanation: ""
    };
    rank.explanation = generateTFImpactRankExplanation(rank) || generateMetabolicImpactExplanation(impact);
    return rank;
  });

  return ranks
    .filter(rank => includeZeroImpact || rank.impactScore > 0)
    .sort((a, b) => b.impactScore - a.impactScore || b.mappedTargetGenes - a.mappedTargetGenes || a.tfLabel.localeCompare(b.tfLabel))
    .slice(0, limit);
}
