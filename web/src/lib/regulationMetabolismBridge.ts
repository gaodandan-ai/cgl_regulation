import {
  MetabolicPathway,
  MetabolicReaction,
  getPathwaysForGene,
  getReactionsForGene,
  normalizeGeneId
} from "./metabolicModelAdapter";

export type RegulatedTargetGene = {
  geneId: string;
  node?: unknown;
  confidence?: number;
  regulation?: "activation" | "repression" | "unknown" | string;
};

export type MetabolicGeneContext = {
  geneId: string;
  reactions: MetabolicReaction[];
  pathways: MetabolicPathway[];
};

export type MetabolicImpactPathwaySummary = {
  pathwayId: string;
  pathwayName: string;
  geneCount: number;
  reactionCount: number;
  genes: string[];
  reactions: string[];
};

export type MetabolicImpact = {
  tfId: string;
  totalTargetGenes: number;
  mappedTargetGenes: number;
  totalReactions: number;
  totalPathways: number;
  pathwaySummary: MetabolicImpactPathwaySummary[];
};

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
  if (graph.nodes) return toArray(graph.nodes);
  if (graph.elements?.nodes) return toArray(graph.elements.nodes);
  if (Array.isArray(graph.elements)) {
    return graph.elements.filter((item: any) => item.group === "nodes" || item.data?.id);
  }
  return [];
}

function collectEdges(graph: any): any[] {
  if (!graph) return [];
  if (typeof graph.edges === "function") return toArray(graph.edges());
  if (graph.edges) return toArray(graph.edges);
  if (graph.elements?.edges) return toArray(graph.elements.edges);
  if (Array.isArray(graph.elements)) {
    return graph.elements.filter((item: any) => item.group === "edges" || item.data?.source);
  }
  return [];
}

function nodeId(node: any): string {
  return String(callOrRead(node, "id") || "");
}

function edgeSource(edge: any): string {
  return String(callOrRead(edge, "source") || "");
}

function edgeTarget(edge: any): string {
  return String(callOrRead(edge, "target") || "");
}

function isGeneLike(node: any, id: string): boolean {
  const type = String(callOrRead(node, "type") || callOrRead(node, "nodeType") || "").toLowerCase();
  if (["gene", "target", "mrna", "transcript", "orf"].includes(type)) return true;
  return /^(cg|cgl)\w+/i.test(id);
}

function regulationType(edge: any): string {
  const raw = callOrRead(edge, "regulation") || callOrRead(edge, "type") || callOrRead(edge, "role") || "unknown";
  if (typeof raw === "object") return "unknown";
  const normalized = String(raw).toLowerCase();
  if (normalized === "a" || normalized.includes("activ")) return "activation";
  if (normalized === "r" || normalized.includes("repress") || normalized.includes("inhibit")) return "repression";
  return normalized || "unknown";
}

function edgeConfidence(edge: any): number | undefined {
  const raw = callOrRead(edge, "confidence") ?? callOrRead(edge, "confidenceScore");
  const parsed = Number(raw);
  return Number.isFinite(parsed) ? parsed : undefined;
}

export function getTargetGenesForTF(graph: any, tfId: string): RegulatedTargetGene[] {
  const normalizedTf = normalizeGeneId(tfId);
  if (!normalizedTf) return [];

  const nodesById = new Map<string, any>();
  for (const node of collectNodes(graph)) {
    const id = nodeId(node);
    if (id) nodesById.set(normalizeGeneId(id), node);
  }

  const targets = new Map<string, RegulatedTargetGene>();
  for (const edge of collectEdges(graph)) {
    const source = normalizeGeneId(edgeSource(edge));
    const targetId = edgeTarget(edge);
    const target = normalizeGeneId(targetId);
    if (!target || source !== normalizedTf) continue;

    const targetNode = nodesById.get(target);
    if (targetNode && !isGeneLike(targetNode, targetId)) continue;
    if (!targetNode && !isGeneLike(null, targetId)) continue;

    if (!targets.has(target)) {
      targets.set(target, {
        geneId: target,
        node: targetNode,
        confidence: edgeConfidence(edge),
        regulation: regulationType(edge)
      });
    }
  }

  return Array.from(targets.values());
}

export function getMetabolicContextForGene(geneId: string): MetabolicGeneContext {
  const normalizedGene = normalizeGeneId(geneId);
  if (!normalizedGene) {
    return { geneId: "", reactions: [], pathways: [] };
  }

  return {
    geneId: normalizedGene,
    reactions: getReactionsForGene(normalizedGene) || [],
    pathways: getPathwaysForGene(normalizedGene) || []
  };
}

export function getMetabolicImpactForTF(graph: any, tfId: string): MetabolicImpact {
  const targets = getTargetGenesForTF(graph, tfId);
  const reactionIds = new Set<string>();
  const pathwayStats = new Map<string, { name: string; genes: Set<string>; reactions: Set<string> }>();
  const mappedGenes = new Set<string>();

  for (const target of targets) {
    const reactions = getReactionsForGene(target.geneId) || [];
    if (reactions.length > 0) mappedGenes.add(target.geneId);

    for (const reaction of reactions) {
      const reactionKey = `${reaction.model || "model"}:${reaction.id || reaction.label || "reaction"}`;
      reactionIds.add(reactionKey);

      const pathwayId = reaction.pathway_id || reaction.pathway_name || "Unassigned pathway";
      const pathwayName = reaction.pathway_name || pathwayId;
      const stat = pathwayStats.get(pathwayId) || {
        name: pathwayName,
        genes: new Set<string>(),
        reactions: new Set<string>()
      };
      stat.genes.add(target.geneId);
      stat.reactions.add(reaction.id || reaction.label || reactionKey);
      pathwayStats.set(pathwayId, stat);
    }
  }

  const pathwaySummary = Array.from(pathwayStats.entries())
    .map(([pathwayId, stat]) => ({
      pathwayId,
      pathwayName: stat.name,
      geneCount: stat.genes.size,
      reactionCount: stat.reactions.size,
      genes: Array.from(stat.genes).sort(),
      reactions: Array.from(stat.reactions).sort()
    }))
    .sort((a, b) =>
      b.geneCount - a.geneCount ||
      b.reactionCount - a.reactionCount ||
      a.pathwayName.localeCompare(b.pathwayName)
    );

  return {
    tfId: normalizeGeneId(tfId),
    totalTargetGenes: targets.length,
    mappedTargetGenes: mappedGenes.size,
    totalReactions: reactionIds.size,
    totalPathways: pathwaySummary.length,
    pathwaySummary
  };
}

export function generateMetabolicImpactExplanation(impact?: Partial<MetabolicImpact> | null): string {
  if (!impact || !impact.mappedTargetGenes || !impact.totalReactions) {
    return "No metabolic model mapping available for this node.";
  }

  const topPathways = (impact.pathwaySummary || [])
    .slice(0, 3)
    .map(pathway => pathway.pathwayName)
    .filter(Boolean);

  const pathwayText = topPathways.length > 0
    ? ` The most affected pathways include ${topPathways.join(", ")}.`
    : " No pathway annotations are available for the mapped reactions.";

  return `This transcription factor regulates ${impact.totalTargetGenes || 0} target genes, ${impact.mappedTargetGenes || 0} of which are mapped to metabolic reactions.${pathwayText}`;
}
