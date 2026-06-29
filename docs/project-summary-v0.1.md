# Cgl Regulation Explorer v0.1.1 Project Summary

**Cgl Regulation Explorer** is an interactive, reproducible regulatory-metabolic analysis platform tailored for *Corynebacterium glutamicum* DSM 20300 / ATCC 13032. It enables researchers to trace transcription factor (TF) and sRNA regulatory links and project them directly onto metabolic reaction models and pathways.

---

## 1. Core Platform Capabilities
The platform supports three core analysis workflows along with model quality auditing and reproducible example case studies:
- **Gene / TF Explorer**: Traces downstream targets and upstream regulators of a given locus tag or common name, visualizing evidence attributes (motif location, ChIP-seq support, expression correlation, operon grouping) and metabolic impacts.
- **Pathway View**: Centered on a metabolic pathway keyword (e.g., TCA cycle, glutamate, lysine), this workflow aggregates all active pathway genes/reactions and discovers upstream regulators, ranking them by pathway coverage.
- **Engineering Targets**: Prioritizes transcription factors based on their potential global metabolic impact. Priority is calculated using TF target count, regulation mode (repression/activation), mapped reaction count, average edge confidence score, and overlap with key pathways.
- **Data & Model Quality Dashboard**: Transparent audit of total regulatory nodes/edges, edge confidence score distributions, iCW773 metabolic model mapping coverage, and ecCGL1 enzyme-constrained support.
- **Built-in Case Studies**: Reproducible example workflows showcasing glutamate regulation analysis, TCA cycle upstream regulator discovery, and amino acid biosynthesis target ranking.

---

## 2. Integrated Datasets
- **Transcriptional Regulatory Network (TRN)**: Curated database of transcription factors, locus tags, and operons.
- **Post-transcriptional sRNA-mRNA network**: Curated and predicted sRNA connections.
- **Edge Confidence Scores**: Machine learning scores (Random Forest predicted probabilities) based on combined molecular and physical evidence features.
- **iCW773 Metabolic Model**: Mapping genes to reactions, stoichiometric equations, subsystems, and KEGG pathways.
- **ecCGL1 Model Annotations**: kcat, molecular weight (MW), EC numbers, and UniProt database IDs.

---

## 3. Current Limitations

> [!WARNING]
> - **Hypothesis Generation**: All prioritizing targets and scoring models are hypothesis-generating. They represent potential metabolic impact rather than experimentally verified rates or outputs.
> - **No Flux Simulations**: The platform maps network connections to model elements but does not run Flux Balance Analysis (FBA) or kinetic simulations.
> - **Score Interpretation**: Edge confidence scores represent prioritization markers and should not be interpreted as biophysical probabilities.
> - **Experimental Validation**: Wet-lab validation is required before drawing biological engineering conclusions.

---

## 4. References & Data Attribution

Original resources remain attributable to their respective authors and providers. Please cite the original resources when using derived platform results:

### Regulatory Network Resources
- **CoryneRegNet 7**: Parise, M.T.D. et al. *CoryneRegNet 7, the reference database and analysis platform for corynebacterial gene regulatory networks.* Scientific Data 7, 142 (2020).
- **Extended C. glutamicum TRN**: Parise, M.T.D. et al. *The transcriptional regulatory network of Corynebacterium glutamicum: an update.* BMC Genomics 12, 608 (2011).

### Genome-Scale Metabolic Models
- **iCW773 Model**: Zhang, Y. et al. *A new genome-scale metabolic model of Corynebacterium glutamicum and its application.* Biotechnology for Biofuels 10, 169 (2017).
- **iCGB21FR Model (Support Context)**: Feierabend, M. et al. *High-Quality Genome-Scale Reconstruction of Corynebacterium glutamicum ATCC 13032.* Frontiers in Microbiology 12, 750206 (2021).

### Enzyme-Constrained Model
- **ecCGL1 Model**: Niu, J. et al. *Construction and Analysis of an Enzyme-Constrained Metabolic Model of Corynebacterium glutamicum.* Biomolecules 12, 1499 (2022).

### Pathway & Annotation Resources
- **KEGG Database**: Kanehisa, M. and Goto, S. *KEGG: Kyoto Encyclopedia of Genes and Genomes.* Nucleic Acids Research 28, 27–30 (2000).
- **KEGG Pathways Update**: Kanehisa, M. et al. *KEGG for taxonomy-based analysis of pathways and genomes.* Nucleic Acids Research 51, D587–D592 (2023).

---

## 5. How to Cite & Usage Disclaimer

### Platform Citation
If you use Cgl Regulation Explorer in your research, please cite:
```text
Cgl Regulation Explorer: an interactive regulatory-metabolic analysis platform for Corynebacterium glutamicum. Version 0.1.1. Available at: https://cgl-regulation.vercel.app/
```
*Note: A future DOI / preprint citation will be updated here once published.*

### Usage Disclaimer
Users should cite both this platform and the original data/model resources (CoryneRegNet, iCW773, ecCGL1, and KEGG) when using exported results, figures, or derived hypotheses.
