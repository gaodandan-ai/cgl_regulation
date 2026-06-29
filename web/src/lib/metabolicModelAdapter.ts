export type MetabolicGene = {
  id: string;
  label: string;
  name?: string;
};

export type MetabolicReaction = {
  id: string;
  label: string;
  model?: string;
  equation?: string;
  gpr_rule?: string;
  pathway_id?: string;
  pathway_name?: string;
  source_file?: string;
  ec_number?: string;
  kcat?: number;
  molecular_weight?: number;
  kcat_MW?: number;
  uniprot_ids?: string[];
  protein_masses?: number[];
  variant_of?: string;
  reaction_variant?: string;
  lower_bound?: number;
  upper_bound?: number;
  kcat_source_count?: number;
  kcat_species_sample?: string[];
  enzyme_constraint?: {
    model_variant?: string;
    kcat?: number;
    molecular_weight?: number;
    kcat_MW?: number;
    uniprot_ids?: string[];
    protein_masses?: number[];
    ec_number?: string;
    variant_of?: string;
    lower_bound?: number;
    upper_bound?: number;
    kcat_source_count?: number;
    kcat_species_sample?: string[];
    global_parameters?: Record<string, unknown>;
  };
};

export type MetabolicPathway = {
  id: string;
  name: string;
  model?: string;
  gene_count?: number;
  reaction_count?: number;
  genes?: string[];
  reactions?: string[];
};

export type GeneReactionMapping = {
  geneId: string;
  reactions: MetabolicReaction[];
};

export type RegulatoryGraphEdge = {
  source?: string;
  target?: string;
  type?: string;
  data?: {
    source?: string;
    target?: string;
    type?: string;
  };
};

export type RegulatoryGraph = {
  edges?: RegulatoryGraphEdge[];
};

export type MetabolicImpactSummary = {
  tfId: string;
  targetGeneCount: number;
  mappedGeneCount: number;
  reactionCount: number;
  pathwayCount: number;
  pathways: MetabolicPathway[];
  mappedGenes: GeneReactionMapping[];
};

const impactCache = new Map<string, any>();
let pathwayMappingCache: any = null;

export function normalizeGeneId(id: string): string {
  return String(id || "")
    .trim()
    .replace(/^gene[:_]/i, "")
    .replace(/^G_/i, "")
    .toLowerCase();
}

export function getReactionsForGene(geneId: string): MetabolicReaction[] {
  const normalized = normalizeGeneId(geneId);
  const cached = impactCache.get(normalized);
  const gene = (cached?.affected_genes || []).find(
    (item: any) => normalizeGeneId(item.locus) === normalized
  );
  return gene?.reactions || [];
}

export function getPathwaysForGene(geneId: string): MetabolicPathway[] {
  const pathways = new Map<string, MetabolicPathway>();

  for (const reaction of getReactionsForGene(geneId)) {
    const id = reaction.pathway_id || reaction.pathway_name || "Unassigned pathway";
    const pathway = pathways.get(id) || {
      id,
      name: reaction.pathway_name || id,
      model: reaction.model,
      gene_count: 1,
      reaction_count: 0,
      genes: [normalizeGeneId(geneId)],
      reactions: []
    };
    if (!pathway.reactions?.includes(reaction.id)) {
      pathway.reactions?.push(reaction.id);
    }
    pathway.reaction_count = pathway.reactions?.length || 0;
    pathways.set(id, pathway);
  }

  return Array.from(pathways.values());
}

export function getEnzymeConstrainedReactionsForGene(geneId: string): MetabolicReaction[] {
  return getReactionsForGene(geneId).filter(reaction =>
    Boolean(reaction.enzyme_constraint || reaction.kcat || reaction.molecular_weight || reaction.kcat_MW)
  );
}

export function getReactionVariantsForGene(geneId: string): Record<string, MetabolicReaction[]> {
  const grouped: Record<string, MetabolicReaction[]> = {};
  for (const reaction of getEnzymeConstrainedReactionsForGene(geneId)) {
    const baseId = reaction.variant_of || reaction.id;
    if (!grouped[baseId]) grouped[baseId] = [];
    grouped[baseId].push(reaction);
  }
  return grouped;
}

export function getMetabolicImpactForTF(tfId: string, graph?: RegulatoryGraph): MetabolicImpactSummary {
  const normalizedTf = normalizeGeneId(tfId);
  const targets = new Set<string>();

  for (const edge of graph?.edges || []) {
    const source = normalizeGeneId(edge.source || edge.data?.source || "");
    const target = normalizeGeneId(edge.target || edge.data?.target || "");
    const edgeType = edge.type || edge.data?.type || "regulates";
    if (source === normalizedTf && target && edgeType === "regulates") {
      targets.add(target);
    }
  }

  const mappedGenes: GeneReactionMapping[] = [];
  const reactionIds = new Set<string>();
  const pathwayStats = new Map<string, MetabolicPathway>();

  for (const geneId of targets) {
    const reactions = getReactionsForGene(geneId);
    if (reactions.length === 0) continue;
    mappedGenes.push({ geneId, reactions });

    for (const reaction of reactions) {
      reactionIds.add(`${reaction.model || "model"}:${reaction.id}`);
      const pathwayId = reaction.pathway_id || reaction.pathway_name || "Unassigned pathway";
      const pathway = pathwayStats.get(pathwayId) || {
        id: pathwayId,
        name: reaction.pathway_name || pathwayId,
        model: reaction.model,
        gene_count: 0,
        reaction_count: 0,
        genes: [],
        reactions: []
      };
      if (!pathway.genes?.includes(geneId)) pathway.genes?.push(geneId);
      if (!pathway.reactions?.includes(reaction.id)) pathway.reactions?.push(reaction.id);
      pathway.gene_count = pathway.genes?.length || 0;
      pathway.reaction_count = pathway.reactions?.length || 0;
      pathwayStats.set(pathwayId, pathway);
    }
  }

  const pathways = Array.from(pathwayStats.values()).sort((a, b) =>
    (b.gene_count || 0) - (a.gene_count || 0) ||
    (b.reaction_count || 0) - (a.reaction_count || 0) ||
    a.name.localeCompare(b.name)
  );

  return {
    tfId: normalizedTf,
    targetGeneCount: targets.size,
    mappedGeneCount: mappedGenes.length,
    reactionCount: reactionIds.size,
    pathwayCount: pathways.length,
    pathways,
    mappedGenes
  };
}

export async function loadMetabolicImpact(nodeId: string): Promise<any> {
  const response = await fetch(`/api/metabolic_impact?gene=${encodeURIComponent(nodeId)}`);
  if (!response.ok) {
    throw new Error(`Metabolic impact request failed: ${response.status}`);
  }
  const payload = await response.json();
  impactCache.set(normalizeGeneId(nodeId), payload);

  for (const gene of payload.affected_genes || []) {
    impactCache.set(normalizeGeneId(gene.locus), {
      ...payload,
      mode: "gene",
      affected_genes: [gene]
    });
  }

  return payload;
}

export async function loadMetabolicPathways(query = ""): Promise<any> {
  if (!query && pathwayMappingCache) return pathwayMappingCache;
  const suffix = query ? `?query=${encodeURIComponent(query)}` : "";
  const response = await fetch(`/api/metabolic_pathways${suffix}`);
  if (!response.ok) {
    throw new Error(`Metabolic pathway request failed: ${response.status}`);
  }
  const payload = await response.json();
  if (!query) {
    pathwayMappingCache = payload;
    // Pre-populate impactCache with all gene-reaction mappings from the pathways payload
    if (payload && Array.isArray(payload.pathways)) {
      payload.pathways.forEach(pathway => {
        if (Array.isArray(pathway.genes)) {
          pathway.genes.forEach(g => {
            const geneId = normalizeGeneId(g.geneId || g.locus);
            if (!geneId) return;
            
            // Map reactions to standard schema format
            const formattedReactions = (g.reactions || []).map((r: any) => ({
              id: r.reactionId || r.id,
              label: r.reactionName || r.label || r.reactionId || r.id,
              model: r.model || pathway.model,
              ec_number: r.ecNumber || r.ec_number,
              kcat: r.kcat,
              molecular_weight: r.molecularWeight || r.molecular_weight,
              kcat_MW: r.kcatMW || r.kcat_MW,
              uniprot_ids: r.uniprotIds || r.uniprot_ids || [],
              reaction_variant: r.reactionVariant || r.reaction_variant,
              enzyme_constraint: r.enzymeConstraint || r.enzyme_constraint || {},
              pathway_id: pathway.pathwayId,
              pathway_name: pathway.pathwayName
            }));

            const existing = impactCache.get(geneId);
            if (existing) {
              // Merge reactions avoiding duplicates
              const mergedRxns = [...(existing.affected_genes[0].reactions || [])];
              formattedReactions.forEach((r: any) => {
                if (!mergedRxns.some((er: any) => er.id === r.id && er.model === r.model)) {
                  mergedRxns.push(r);
                }
              });
              existing.affected_genes[0].reactions = mergedRxns;
            } else {
              impactCache.set(geneId, {
                mode: "gene",
                query: geneId,
                affected_genes: [{
                  locus: geneId,
                  name: g.geneLabel || g.name || geneId,
                  reactions: formattedReactions
                }]
              });
            }
          });
        }
      });
    }
  }
  return payload;
}

export function getAllMetabolicPathways(): any[] {
  return pathwayMappingCache?.pathways || [];
}

export function getGenesForPathway(pathwayIdOrName: string): any[] {
  const query = normalizeGeneId(pathwayIdOrName);
  if (!query) return [];
  return (pathwayMappingCache?.pathways || [])
    .filter((pathway: any) => {
      const id = normalizeGeneId(pathway.pathwayId || pathway.id);
      const name = normalizeGeneId(pathway.pathwayName || pathway.name);
      return id === query || name === query || id.includes(query) || name.includes(query);
    })
    .flatMap((pathway: any) => pathway.genes || []);
}
