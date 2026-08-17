import {
  getPathwayRegulatorySummary
} from "./pathwayRegulatoryView";
import {
  findEngineeringTargetCandidates
} from "./candidateEngineeringTargets";
import {
  getMetabolicMappingCoverage,
  getEcCGL1CoverageSummary
} from "./analysisQuality";

export type CaseStudy = {
  id: string;
  title: string;
  subtitle: string;
  question: string;
  workflow: string[];
  entryMode: "gene" | "tf" | "pathway" | "engineering-targets";
  query?: string;
  pathwayKeyword?: string;
  tfId?: string;
  geneId?: string;
  description: string;
  expectedOutputs: string[];
  limitations: string[];
};

export const BUILT_IN_CASE_STUDIES: CaseStudy[] = [
  {
    id: "glutamate-regulation",
    title: "Glutamate-associated Regulatory Analysis",
    subtitle: "Identify upstream TFs regulating glutamate pathways",
    question: "Which transcription factors may regulate glutamate-associated metabolic modules?",
    workflow: [
      "Identify all genes and reactions mapped to glutamate metabolism",
      "Retrieve upstream transcription factors targeting these genes",
      "Evaluate regulatory confidence scores and ecCGL1 constraint coverage"
    ],
    entryMode: "pathway",
    pathwayKeyword: "glutamate",
    description: "This case study focuses on highlighting regulatory networks surrounding glutamate metabolism, mapping TF control over glutamate biosynthesis and transport.",
    expectedOutputs: [
      "Upstream transcription factors and confidence scores",
      "List of regulated genes and metabolic reactions in the glutamate pathway",
      "Enzyme constraints coverage details (kcat, molecular weight)"
    ],
    limitations: [
      "Does not simulate metabolic fluxes or growth rates",
      "Confidence scores are based on predictive models and evidence curation, not direct new wet-lab experiments"
    ]
  },
  {
    id: "tca-cycle-regulators",
    title: "TCA Cycle Upstream Regulator Discovery",
    subtitle: "Discover regulators of the TCA cycle reactions",
    question: "Which regulators may influence TCA cycle reactions?",
    workflow: [
      "Identify genes and reactions involved in the citric acid cycle (TCA cycle)",
      "Retrieve all transcription factors regulating these structural genes",
      "Sort regulators by target gene count and average link confidence"
    ],
    entryMode: "pathway",
    pathwayKeyword: "tca",
    description: "This case study focuses on the tricarboxylic acid (TCA) cycle, discovering transcription factors that target structural genes in citrate, succinate, or malate reactions.",
    expectedOutputs: [
      "TCA cycle structural genes and upstream regulators",
      "Average confidence and target counts for each regulator",
      "Heuristic vs RF confidence scoring comparisons"
    ],
    limitations: [
      "Regulatory links may have differing physiological significance under varying growth conditions",
      "Excludes dynamic kinetic interactions or feedforward loops"
    ]
  },
  {
    id: "amino-acid-engineering-targets",
    title: "Amino Acid Biosynthesis Engineering Target Ranking",
    subtitle: "Rank upstream regulators by amino acid pathway impact",
    question: "Which TFs are prioritized as candidate regulators for amino acid biosynthesis?",
    workflow: [
      "Query candidates prioritized by global metabolic impact",
      "Filter candidates by regulation of amino acid pathways",
      "Check recommendation levels, kcat coverage, and regulatory mode"
    ],
    entryMode: "engineering-targets",
    pathwayKeyword: "amino acid",
    description: "This case study ranks transcription factors based on their potential impact on amino acid biosynthesis pathways, prioritizing candidates that regulate multiple pathway reactions with high confidence.",
    expectedOutputs: [
      "Prioritized transcription factors ranking table",
      "Recommendation levels (high/medium/low)",
      "Regulated reaction and pathway counts"
    ],
    limitations: [
      "Does not account for downstream transcriptional feedforward loops",
      "Prioritization scores represent metabolic impact potential, not guaranteed metabolic outcomes"
    ]
  }
];

export function generateCaseStudyNarrative(input: { caseStudy: CaseStudy; results: any }): string {
  const { caseStudy, results } = input;
  if (!caseStudy) return "";

  if (caseStudy.id === "glutamate-regulation") {
    const totalRegulators = results.pathwaySummary?.totalRegulators || 0;
    const pathwayName = results.pathwaySummary?.pathwayName || "glutamate metabolism";
    return `This case study searches for glutamate-associated pathways and identifies upstream transcription factors regulating mapped pathway genes. In the current network, ${totalRegulators} transcription factors may regulate genes associated with ${pathwayName}. Candidate TFs are ranked using regulatory confidence, metabolic mapping coverage, and pathway relevance. These predictions are hypothesis-generating and require experimental validation.`;
  }
  
  if (caseStudy.id === "tca-cycle-regulators") {
    const totalRegulators = results.pathwaySummary?.totalRegulators || 0;
    const pathwayName = results.pathwaySummary?.pathwayName || "TCA cycle";
    return `This case study retrieves TCA cycle reactions and structural genes, identifying upstream transcription factors targeting citrate, succinate, or malate reactions. We identified ${totalRegulators} regulators that may influence TCA cycle genes. These regulatory linkages are predicted and associated with metabolic control, serving as hypothesis-generating insights that require wet-lab experimental validation.`;
  }

  if (caseStudy.id === "amino-acid-engineering-targets") {
    const candidateCount = results.engineeringCandidates?.length || 0;
    return `This case study ranks all transcription factors in the regulatory network based on their potential impact on amino acid biosynthesis pathways. We prioritized ${candidateCount} candidate regulators that regulate multiple reactions with high confidence. These target prioritizations represent predicted metabolic impact potential and are hypothesis-generating. They require experimental validation to confirm metabolic effects.`;
  }

  return "";
}

export function runCaseStudy(caseStudyId: string, graph: any) {
  const caseStudy = BUILT_IN_CASE_STUDIES.find(cs => cs.id === caseStudyId);
  if (!caseStudy) {
    throw new Error(`Case study with ID '${caseStudyId}' not found.`);
  }

  const results: {
    pathwaySummary?: any;
    tfRanking?: any[];
    engineeringCandidates?: any[];
    qualitySummary?: any;
    enzymeConstraintSummary?: any;
  } = {};

  const warnings: string[] = [];

  try {
    const qual = getMetabolicMappingCoverage(graph);
    const enz = getEcCGL1CoverageSummary(graph);
    results.qualitySummary = qual;
    results.enzymeConstraintSummary = enz;

    // Warnings for metabolic mapping and enzyme constraint coverage are disabled per user request

    if (caseStudy.id === "glutamate-regulation") {
      results.pathwaySummary = getPathwayRegulatorySummary(graph, "glutamate");
      if (results.pathwaySummary && results.pathwaySummary.regulators) {
        results.tfRanking = results.pathwaySummary.regulators;
      }
    } else if (caseStudy.id === "tca-cycle-regulators") {
      results.pathwaySummary = getPathwayRegulatorySummary(graph, "tca");
      if (results.pathwaySummary && results.pathwaySummary.regulators) {
        results.tfRanking = results.pathwaySummary.regulators;
      }
    } else if (caseStudy.id === "amino-acid-engineering-targets") {
      results.engineeringCandidates = findEngineeringTargetCandidates(graph, {
        pathwayKeywordFilter: "amino acid",
        limit: 10
      });
    }
  } catch (err) {
    console.warn(`Error running case study '${caseStudyId}':`, err);
  }

  const narrative = generateCaseStudyNarrative({ caseStudy, results });

  return {
    caseStudy,
    results,
    narrative,
    warnings
  };
}
