import {
  normalizeGeneId,
  getReactionsForGene,
  getPathwaysForGene,
  getEnzymeConstrainedReactionsForGene
} from "./metabolicModelAdapter";

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
  if (value.data && typeof value.data === "object" && key in value.data) return value.data[key];
  if (key in value) return value[key];
  return undefined;
}

function collectNodes(graph: any): any[] {
  if (!graph) return [];
  if (typeof graph.nodes === "function") return toArray(graph.nodes());
  if (graph.nodes) return toArray(graph.nodes);
  if (graph.elements?.nodes) return toArray(graph.elements.nodes);
  if (Array.isArray(graph.elements)) {
    return graph.elements.filter((item: any) => item.group === "nodes" || (item.data && item.data.id));
  }
  return [];
}

// Ensure support for both flat elements and Cytoscape collections
function collectEdges(graph: any): any[] {
  if (!graph) return [];
  if (typeof graph.edges === "function") return toArray(graph.edges());
  if (graph.edges) return toArray(graph.edges);
  if (graph.elements?.edges) return toArray(graph.elements.edges);
  if (Array.isArray(graph.elements)) {
    return graph.elements.filter((item: any) => item.group === "edges" || (item.data && item.data.source));
  }
  return [];
}

function nodeId(node: any): string {
  return String(callOrRead(node, "id") || "");
}

// 1. Add regulatory network summary
export function getRegulatoryNetworkSummary(graph: any) {
  const nodes = collectNodes(graph);
  const edges = collectEdges(graph);

  let tfCount = 0;
  let srnaCount = 0;
  let operonCount = 0;
  let geneCount = 0;

  for (const node of nodes) {
    const id = nodeId(node);
    const type = String(callOrRead(node, "type") || callOrRead(node, "nodeType") || "").toLowerCase();
    
    if (["tf", "transcription_factor", "regulator"].includes(type)) {
      tfCount++;
    } else if (type === "srna") {
      srnaCount++;
    } else if (type === "operon") {
      operonCount++;
    } else if (["gene", "target", "mrna", "transcript", "orf"].includes(type) || /^(cg|cgl)\w+/i.test(id)) {
      geneCount++;
    }
  }

  let tfGeneEdgeCount = 0;
  let srnaEdgeCount = 0;
  let activationCount = 0;
  let repressionCount = 0;
  let predictedCount = 0;
  let unknownRegulationCount = 0;

  for (const edge of edges) {
    const interaction = String(callOrRead(edge, "interactionClass") || "").toLowerCase();
    const sType = String(callOrRead(edge, "sourceType") || "").toLowerCase();
    const role = String(callOrRead(edge, "role") || "").toLowerCase();
    const reg = String(callOrRead(edge, "regulationType") || callOrRead(edge, "regulation") || callOrRead(edge, "type") || "").toLowerCase();

    if (interaction === "srna-mrna" || sType === "srna" || role === "srna" || reg.includes("post_transcriptional")) {
      srnaEdgeCount++;
    } else {
      tfGeneEdgeCount++;
    }

    if (reg === "activation" || reg.includes("activ") || role === "a") {
      activationCount++;
    } else if (reg === "repression" || reg.includes("repress") || role === "r" || reg.includes("post_transcriptional")) {
      repressionCount++;
    } else if (reg.includes("predict") || role.includes("predict")) {
      predictedCount++;
    } else {
      unknownRegulationCount++;
    }
  }

  return {
    totalNodes: nodes.length,
    totalEdges: edges.length,
    tfCount,
    geneCount,
    srnaCount,
    operonCount,
    tfGeneEdgeCount,
    srnaEdgeCount,
    activationCount,
    repressionCount,
    predictedCount,
    unknownRegulationCount
  };
}

// 2. Add confidence score summary
export function getConfidenceScoreSummary(graph: any) {
  const edges = collectEdges(graph);
  
  let totalEdgesWithConfidence = 0;
  let sumConfidence = 0;
  const scores: number[] = [];
  
  let rfConfidenceAvailableCount = 0;
  let heuristicConfidenceAvailableCount = 0;
  let highConfidenceEdgeCount = 0;
  let mediumConfidenceEdgeCount = 0;
  let lowConfidenceEdgeCount = 0;

  const rfScores: number[] = [];
  const heuristicScores: number[] = [];
  let sumAbsoluteDifference = 0;
  let absoluteDifferenceCount = 0;

  for (const edge of edges) {
    const scoreVal = callOrRead(edge, "confidenceScore") ?? callOrRead(edge, "confidence");
    const heuristicVal = callOrRead(edge, "heuristicConfidenceScore") ?? callOrRead(edge, "heuristicConfidence");
    const rfVal = callOrRead(edge, "predictedConfidence") ?? callOrRead(edge, "rfConfidence") ?? (edge.confidenceFactors?.randomForest);

    const score = typeof scoreVal === "number" && !Number.isNaN(scoreVal) ? scoreVal : undefined;
    const heuristic = typeof heuristicVal === "number" && !Number.isNaN(heuristicVal) ? heuristicVal : undefined;
    const rf = typeof rfVal === "number" && !Number.isNaN(rfVal) ? rfVal : undefined;

    if (score !== undefined) {
      totalEdgesWithConfidence++;
      sumConfidence += score;
      scores.push(score);

      if (score >= 0.75) {
        highConfidenceEdgeCount++;
      } else if (score >= 0.45) {
        mediumConfidenceEdgeCount++;
      } else {
        lowConfidenceEdgeCount++;
      }
    }

    if (rf !== undefined && rf !== null) {
      rfConfidenceAvailableCount++;
      rfScores.push(rf);
    }
    if (heuristic !== undefined && heuristic !== null) {
      heuristicConfidenceAvailableCount++;
      heuristicScores.push(heuristic);
    }

    if (rf !== undefined && rf !== null && heuristic !== undefined && heuristic !== null) {
      sumAbsoluteDifference += Math.abs(rf - heuristic);
      absoluteDifferenceCount++;
    }
  }

  let averageConfidence = 0;
  if (totalEdgesWithConfidence > 0) {
    averageConfidence = sumConfidence / totalEdgesWithConfidence;
  }

  let medianConfidence = 0;
  if (scores.length > 0) {
    scores.sort((a, b) => a - b);
    const mid = Math.floor(scores.length / 2);
    medianConfidence = scores.length % 2 !== 0 ? scores[mid] : (scores[mid - 1] + scores[mid]) / 2;
  }

  const averageRfConfidence = rfScores.length > 0 ? rfScores.reduce((a, b) => a + b, 0) / rfScores.length : null;
  const averageHeuristicConfidence = heuristicScores.length > 0 ? heuristicScores.reduce((a, b) => a + b, 0) / heuristicScores.length : null;
  const averageAbsoluteDifference = absoluteDifferenceCount > 0 ? sumAbsoluteDifference / absoluteDifferenceCount : null;

  return {
    totalEdgesWithConfidence,
    averageConfidence,
    medianConfidence,
    highConfidenceEdgeCount,
    mediumConfidenceEdgeCount,
    lowConfidenceEdgeCount,
    rfConfidenceAvailableCount,
    heuristicConfidenceAvailableCount,
    averageRfConfidence,
    averageHeuristicConfidence,
    averageAbsoluteDifference
  };
}

// Helper to collect all gene-like IDs from the graph nodes
function getGeneIdsFromGraph(graph: any): Set<string> {
  const nodes = collectNodes(graph);
  const geneIds = new Set<string>();

  for (const node of nodes) {
    const id = nodeId(node);
    const type = String(callOrRead(node, "type") || callOrRead(node, "nodeType") || "").toLowerCase();
    
    if (["tf", "transcription_factor", "regulator", "gene", "target", "mrna", "transcript", "orf"].includes(type) || /^(cg|cgl)\w+/i.test(id)) {
      if (type !== "operon" && type !== "srna") {
        const norm = normalizeGeneId(id);
        if (norm) geneIds.add(norm);
      }
    }
  }

  return geneIds;
}

// 3. Add metabolic mapping coverage summary
export function getMetabolicMappingCoverage(graph: any) {
  const geneIds = getGeneIdsFromGraph(graph);

  let genesMappedToReactions = 0;
  let genesMappedToPathways = 0;
  const uniqueReactions = new Set<string>();
  const uniquePathways = new Set<string>();
  const unmappedGenes: string[] = [];

  for (const geneId of geneIds) {
    const reactions = getReactionsForGene(geneId) || [];
    const pathways = getPathwaysForGene(geneId) || [];

    let mapped = false;
    if (reactions.length > 0) {
      genesMappedToReactions++;
      mapped = true;
      for (const rxn of reactions) {
        if (rxn.id) uniqueReactions.add(rxn.id);
      }
    }

    if (pathways.length > 0) {
      genesMappedToPathways++;
      for (const path of pathways) {
        if (path.id) uniquePathways.add(path.id);
      }
    }

    if (!mapped) {
      unmappedGenes.push(geneId);
    }
  }

  unmappedGenes.sort();

  return {
    regulatoryGeneCount: geneIds.size,
    genesMappedToReactions,
    genesMappedToPathways,
    mappedReactionCount: uniqueReactions.size,
    mappedPathwayCount: uniquePathways.size,
    unmappedGeneCount: unmappedGenes.length,
    unmappedGenes
  };
}

// 4. Add ecCGL1 enzyme coverage summary
export function getEcCGL1CoverageSummary(graph: any) {
  const geneIds = getGeneIdsFromGraph(graph);

  let genesWithEnzymeMapping = 0;
  const enzymeReactions = new Set<string>();
  const rxnsWithKcat = new Set<string>();
  const rxnsWithMW = new Set<string>();
  const rxnsWithKcatMw = new Set<string>();
  const rxnsWithEC = new Set<string>();
  const rxnsWithUniProt = new Set<string>();
  const potentialEnzymeConstrainedReactions = new Set<string>();
  const unmappedEnzymeGenes: string[] = [];

  for (const geneId of geneIds) {
    const rxns = getEnzymeConstrainedReactionsForGene(geneId) || [];
    if (rxns.length > 0) {
      genesWithEnzymeMapping++;
      for (const rxn of rxns) {
        const enzyme = rxn.enzyme_constraint || {};
        const ecNumber = rxn.ec_number || enzyme.ec_number;
        const kcat = rxn.kcat ?? enzyme.kcat;
        const molecularWeight = rxn.molecular_weight ?? enzyme.molecular_weight;
        const kcatMw = rxn.kcat_MW ?? enzyme.kcat_MW;
        const uniprotIds = rxn.uniprot_ids || enzyme.uniprot_ids || [];

        const rxnId = rxn.id;
        if (!rxnId) continue;
        
        enzymeReactions.add(rxnId);

        if (kcat !== undefined && kcat !== null && !Number.isNaN(kcat)) rxnsWithKcat.add(rxnId);
        if (molecularWeight !== undefined && molecularWeight !== null && !Number.isNaN(molecularWeight)) rxnsWithMW.add(rxnId);
        if (kcatMw !== undefined && kcatMw !== null && !Number.isNaN(kcatMw)) rxnsWithKcatMw.add(rxnId);
        if (ecNumber) rxnsWithEC.add(rxnId);
        if (uniprotIds && uniprotIds.length > 0) rxnsWithUniProt.add(rxnId);
        
        if (rxn.enzyme_constraint || rxn.reaction_variant || rxn.variant_of) {
          potentialEnzymeConstrainedReactions.add(rxnId);
        }
      }
    } else {
      unmappedEnzymeGenes.push(geneId);
    }
  }

  unmappedEnzymeGenes.sort();

  return {
    genesWithEnzymeMapping,
    enzymeAssociatedReactionCount: enzymeReactions.size,
    reactionsWithKcat: rxnsWithKcat.size,
    reactionsWithMolecularWeight: rxnsWithMW.size,
    reactionsWithKcatPerMW: rxnsWithKcatMw.size,
    reactionsWithECNumber: rxnsWithEC.size,
    reactionsWithUniProtId: rxnsWithUniProt.size,
    potentialEnzymeConstrainedReactionCount: potentialEnzymeConstrainedReactions.size,
    unmappedEnzymeGenes
  };
}

// 5. Add global quality summary
export function getAnalysisQualityReport(graph: any) {
  return {
    regulatoryNetwork: getRegulatoryNetworkSummary(graph),
    confidenceScores: getConfidenceScoreSummary(graph),
    metabolicMapping: getMetabolicMappingCoverage(graph),
    enzymeConstraintCoverage: getEcCGL1CoverageSummary(graph),
    generatedAt: new Date().toISOString()
  };
}
