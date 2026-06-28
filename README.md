# Cgl Regulation Explorer (v0.1.0)

An interactive, reproducible regulatory-metabolic analysis platform for **Corynebacterium glutamicum** DSM 20300 / ATCC 13032.

Cgl Regulation Explorer connects transcription factor (TF) and sRNA regulatory network evidence with genome-scale metabolic reaction equations and protein constraints to support hypothesis generation and strain engineering target prioritization.

---

## 1. Platform Workflow Pipeline
The platform processes data across 7 distinct analytical phases:
1. **Regulatory Evidence**: Integrates curated databases, binding motifs, ChIP-seq records, and sRNA-mRNA prediction profiles.
2. **RF Confidence Scoring**: Computes prioritization scores using a Random Forest machine learning model trained on evidence feature matrices.
3. **Metabolic Model Mapping**: Maps target genes to reaction equations and pathways in the *iCW773* genome-scale model.
4. **ecCGL1 Enzyme Annotations**: Enriches reactions with protein context, including kcat values, molecular weight (MW), EC numbers, and parent reactions.
5. **Pathway-centered Analysis**: Traces likely upstream regulators controlling active pathway modules.
6. **Prioritized Engineering Targets**: Ranks candidate regulators by metabolic reaction coverage and confidence.
7. **Case Study Reports**: Generates reproducible case reports and downloadable JSON summaries (e.g. Glutamate-associated regulation, TCA cycle upstream regulators, and Amino Acid biosyntheses).

---

## 2. Main Features
- **Gene / TF Explorer**: Query target networks, visual evidence, operon context, and metabolic mappings.
- **Pathway View**: Analyze pathway-specific upstream TFs and reaction lists.
- **Engineering Targets Page**: Rank TFs globally or by specific pathway keyword filters using candidate scores.
- **Data & Model Quality Dashboard**: Audits platform-wide node/edge counts, edge confidence ranges, and model mapping stats.
- **Built-in Case Studies**: Reproducible example workflows with cautious narrative summaries.
- **Data Exports**: Download network maps, table records, and audit metrics as PNG, CSV, or JSON.

---

## 3. Data Sources
- **`data/regulations.csv`**: Curated database of C. glutamicum TF-target regulatory relationships.
- **`data/rna_regulation.csv`**: Post-transcriptional sRNA-mRNA regulatory interactions.
- **`data/operons.csv`**: Operon groupings and annotations.
- **`data/model/iCW773.omex`**: Genome-scale metabolic model used for structural reaction mapping.
- **`data/model/ecCGL1-main/`**: Enzyme-constrained metabolic model resource containing reaction parameters.
- **`data/edge_confidence/`**: Feature matrices and trained Random Forest model files.

---

## 4. Current Limitations

> [!WARNING]
> - **Hypothesis Generation**: Prioritization scores and network impact metrics serve as computational screenings for high-impact hubs, rather than predicting exact in vivo rates or yields.
> - **Optional Local Simulations**: Standard Flux Balance Analysis (FBA), MOMA, and Flux Variability Analysis (FVA) can be run locally using the *iCW773* model by starting the optional FastAPI simulation backend (port 8001). Static web deployments do not run dynamic flux calculations.
> - **ecCGL1 Parameter Integration**: Enzyme kinetic parameters ($k_{\mathrm{cat}}$, molecular weight, EC numbers, and UniProt IDs) are integrated as static pathway-context annotations to aid bottleneck analysis. The platform does not run active enzyme-constrained simulations (ec-FBA) in v0.1.0.
> - **Validation Requirement**: Wet-lab validation is required before drawing biological engineering conclusions.

---

## 5. References & Data Attribution

Please cite the original resources when utilizing derived platform results:

- **CoryneRegNet 7 (TRNs)**: Parise, M.T.D. et al. *CoryneRegNet 7, the reference database and analysis platform for corynebacterial gene regulatory networks.* Scientific Data 7, 142 (2020).
- **Extended C. glutamicum TRN**: Parise, M.T.D. et al. *The transcriptional regulatory network of Corynebacterium glutamicum: an update.* BMC Genomics 12, 608 (2011).
- **iCW773 Metabolic Model**: Zhang, Y. et al. *A new genome-scale metabolic model of Corynebacterium glutamicum and its application.* Biotechnology for Biofuels 10, 169 (2017).
- **iCGB21FR Model (Support Context)**: Feierabend, M. et al. *High-Quality Genome-Scale Reconstruction of Corynebacterium glutamicum ATCC 13032.* Frontiers in Microbiology 12, 750206 (2021).
- **ecCGL1 Enzyme Constraint Model**: Niu, J. et al. *Construction and Analysis of an Enzyme-Constrained Metabolic Model of Corynebacterium glutamicum.* Biomolecules 12, 1499 (2022).
- **KEGG Pathways**: Kanehisa, M. and Goto, S. *KEGG: Kyoto Encyclopedia of Genes and Genomes.* Nucleic Acids Research 28, 27–30 (2000).

### How to Cite This Platform
If you use Cgl Regulation Explorer in your research, please cite:
```text
Cgl Regulation Explorer: an interactive regulatory-metabolic analysis platform for Corynebacterium glutamicum. Version 0.1.0. Available at: https://cgl-regulation.vercel.app/
```
*Note: Preprints and publication DOIs will be updated here once published.*  
Zenodo DOI: [10.5281/zenodo.placeholder](https://doi.org/10.5281/zenodo.placeholder)

### Usage Disclaimer
Users should cite both this platform and the original data/model resources (CoryneRegNet, iCW773, ecCGL1, and KEGG) when using exported results, figures, or derived hypotheses.

---

## 6. Local Development & Setup

### Install Dependencies
Ensure you have Python 3.8+ installed. Install the required libraries:
```bash
pip install -r requirements.txt
```

### Start the Local Server
Launch the static server and api controller:
```bash
python run_server.py
```
Then open your browser and navigate to:
[http://localhost:8000/index.html](http://localhost:8000/index.html)

### AI Summary Features (Optional)
To activate literature and functional AI-generated summary panels, configure your API keys (Gemini, OpenAI, DeepSeek, etc.) in the left sidebar configuration panel. Key parameters are stored locally in browser `localStorage`.

---

## 7. Project Directory Structure

```text
cgl_regulation/
├── api/                  # Vercel Serverless function entrypoints
│   └── index.py
├── backend/              # FastAPI FBA simulation backend (port 8001)
│   ├── app.py
│   ├── simulation.py
│   ├── model_loader.py
│   └── schemas.py
├── data/                 # Regulators, expressions, and metabolic datasets
│   ├── gene_mapping.csv
│   ├── regulations.csv
│   ├── rna_seq/
│   └── rna_seq_analysis_results.json
├── docs/                 # Platform documentation
├── outputs/              # Generated network visualization figures & HTMLs
├── scratch/              # Temporary and testing scripts
├── scripts/              # Data processing and network visualization scripts
│   ├── visualize_network.py
│   └── rna_seq_network_analysis.py
└── web/                  # Frontend static site files (port 8000)
    ├── app.js
    ├── index.html
    ├── src/              # Frontend TypeScript source code
    └── lib/              # Compiled JavaScript libraries & static assets
```
