Language: [ English ](README.md) | [ 中文版 ](README_zh.md)

# Cgl Regulation Explorer (v1.3.1)

**A Systems Biology and Synthetic Biology Platform for *Corynebacterium glutamicum***

Cgl Regulation Explorer is an integrated computational platform designed for multi-omics data exploration, regulatory network analysis, genome-scale metabolic modeling, and AI-assisted genetic engineering targeting *Corynebacterium glutamicum* (DSM 20300 / ATCC 13032).

---

## Table of Contents

- [System Overview](#system-overview)
- [Key Capabilities](#key-capabilities)
  - [1. Multi-Omics Regulatory Network Visualization](#1-multi-omics-regulatory-network-visualization)
  - [2. Interactive 5-Track Genomic Track Browser](#2-interactive-5-track-genomic-track-browser)
  - [3. Condition-Dependent Omics and iModulon Regulons](#3-condition-dependent-omics-and-imodulon-regulons)
  - [4. Genome-Scale Metabolic Modeling and FBA Simulation](#4-genome-scale-metabolic-modeling-and-fba-simulation)
  - [5. Grounded Multi-Omics AI Engineering Copilot](#5-grounded-multi-omics-ai-engineering-copilot)
- [Repository Structure](#repository-structure)
- [Installation and Deployment](#installation-and-deployment)
  - [Method 1: Desktop GUI Launcher (Windows)](#method-1-desktop-gui-launcher-windows)
  - [Method 2: Development Server (Python FastAPI)](#method-2-development-server-python-fastapi)
  - [Method 3: Containerized Deployment (Docker)](#method-3-containerized-deployment-docker)
- [REST API Reference](#rest-api-reference)
- [Primary Data Sources](#primary-data-sources)
- [Testing and Quality Assurance](#testing-and-quality-assurance)
- [Citation and License](#citation-and-license)

---

## System Overview

*Corynebacterium glutamicum* is an industrial chassis strain extensively utilized for amino acid production and metabolic engineering. Cgl Regulation Explorer bridges the gap between raw transcriptional regulatory datasets, genome-scale metabolic networks, and computational genetic design.

The platform integrates transcription factor-target gene (TF-TG) interactions, sRNA-mRNA post-transcriptional networks, operon structures, position weight matrices (PWMs), experimental ChIP-seq binding peaks, independent component analysis (ICA) iModulon activities, STRING v12 protein-protein interactions, and the iCW773 genome-scale metabolic model (GEM).

---

## Key Capabilities

### 1. Multi-Omics Regulatory Network Visualization
- **Multi-layer Network Topologies**: Renders combined transcription factor, target gene, sRNA, and protein-protein interaction (PPI) networks.
- **Viewport Texture Caching**: Utilizes `textureOnViewport` acceleration for smooth rendering of large-scale graphs (>250 nodes).
- **Level-of-Detail (LOD) Rendering**: Dynamically adjusts label visibility based on zoom scale, maintaining 60 FPS performance during viewport manipulation.

### 2. Interactive 5-Track Genomic Track Browser
- **Track 1 (Genomic Ruler)**: Displays absolute base-pair coordinates and strand orientation ticks.
- **Track 2 (CDS Gene Models)**: Renders directional block glyphs for forward (+) and reverse (-) strand CDS annotations.
- **Track 3 (Promoter and TSS Sites)**: Displays transcription start sites (TSS) and highlights 70bp promoter regions (-35 and -10 consensus boxes).
- **Track 4 (ChIP-seq Binding Density)**: Visualizes signal intensity curves derived from experimental ChIP-seq and CollectTF binding scores.
- **Track 5 (sRNA/ncRNA Annotations)**: Displays non-coding RNA loci.
- **Interactive Controls**: Supports viewport zooming, panning, and sequence inspection.

### 3. Condition-Dependent Omics and iModulon Regulons
- **iModulon Activity Matrix**: Visualizes 87 iModulon activities across 9 condition categories with F1-score overlap alignments.
- **Environmental Sub-networks**: Filters regulatory sub-networks under specific conditions including iron, oxygen, nitrogen, and stress responses.
- **Interactive Charting**: Enables scatter plot region selection with synchronized highlighting on the main network canvas.

### 4. Genome-Scale Metabolic Modeling and FBA Simulation
- **iCW773 Integration**: Implements Flux Balance Analysis (FBA) and Minimization of Metabolic Adjustment (MOMA) algorithms.
- **Perturbation Predictions**: Evaluates gene knockout and overexpression effects on growth rate and target metabolite production fluxes.
- **Thermodynamic Feasibility**: Cross-checks reaction free energy constraints via integrated thermodynamic databases.

### 5. Grounded Multi-Omics AI Engineering Copilot
- **Grounded Omics Context**: Automatically queries RefSeq coordinates, promoter 70bp sequences, TF families, effector molecules, and ChIP-seq peaks from SQLite to ground LLM prompts with empirical data.
- **Specialized Engineering Commands**:
  - `/design-crispri`: Recommends 20bp gRNA target windows and PAM sites for dCas9 knockdown while avoiding core promoter boxes.
  - `/solve-bottleneck`: Identifies rate-limiting metabolic steps and recommends synergistic gene overexpression/knockout combinations.
  - `/promoter-engineering`: Guides promoter mutagenesis library design based on effector-responsive elements.
- **Multi-Provider Support**: Compatible with OpenAI, DeepSeek, Google Gemini, and local offline Ollama models.

---

## Repository Structure

```text
f:\cgl_regulation\
├── backend\                   # FastAPI application core
│   ├── app.py                 # REST API endpoints
│   ├── db_manager.py          # Thread-safe SQLite database manager
│   ├── ai_handlers.py         # Grounded AI context and engineering commands
│   ├── bio_handlers.py        # Regulatory network and pathway analysis
│   ├── graph_engine.py       # Network motif analysis and centrality
│   ├── simulation.py          # FBA and MOMA simulation algorithms
│   ├── model_loader.py        # SBML / iCW773 model parser
│   └── security.py            # Security boundaries and validation
├── web\                       # Native web client
│   ├── index.html             # User interface
│   ├── app.js                 # Frontend application logic
│   ├── style.css              # Style definitions
│   └── lib\                   # Modular client libraries
│       ├── genomicTrackBrowser.js  # 5-track genomic browser
│       ├── geneProfileViewer.js    # 360-degree gene profile viewer
│       ├── networkTopology.js      # Network topology module
│       ├── icaConditionHeatmapView.js # iModulon activity matrix
│       ├── conditionRegulationView.js # Condition regulation module
│       ├── interventionPriorityView.js # Engineering target priorities
│       └── vendor\            # Local third-party vendor libraries
├── data\                      # Database assets
│   ├── public_database.db     # Main SQLite database
│   └── literature_cache.json  # RAG literature cache
├── data_pipeline\             # Data ETL pipeline
│   ├── cli.py                 # Pipeline CLI entry
│   └── scripts\               # Omics parsing scripts
├── tests\                     # Pytest API, data, security, and model tests
├── scripts\                   # Development and utility scripts
│   ├── pipeline\              # GEO, ChIP-seq, and thermodynamic processing
│   ├── analysis\              # Network hierarchy and cross-strain analysis
│   ├── build\                 # Build and packaging utilities
│   └── archive\           # Archived legacy scripts
├── launcher.pyw               # Desktop GUI launcher (PyWebView)
├── run_server.py              # Development server entry point
├── cgl_regulation.spec        # PyInstaller specification file
├── requirements.txt           # Python dependency specification
├── Dockerfile                 # Container specification
└── README.md                  # Project documentation
```

---

## Installation and Deployment

### Requirements
- **Python**: Version 3.10 or higher
- **Operating System**: Windows 10/11, macOS, Linux

---

### Method 1: Cross-Platform Desktop GUI Launcher (Windows / macOS)

- **Windows Users**: Execute the batch script from the project root:
  ```cmd
  启动客户端.bat
  ```
- **macOS Users**: Execute the macOS launcher script in Terminal:
  ```bash
  chmod +x start_mac.sh
  ./start_mac.sh
  ```
The launcher starts the FastAPI backend server in the background and opens the native desktop container (WKWebView on macOS; WebSockets/WebView2 on Windows).

---

### Method 2: Development Server (Python FastAPI)

1. Clone the repository and install dependencies:
   ```bash
   git clone https://github.com/gaodandan-ai/cgl_regulation.git
   cd cgl_regulation
   pip install -r requirements.txt
   ```

2. Start the development server:
   ```bash
   python run_server.py
   ```
   Alternatively, launch via Uvicorn:
   ```bash
   uvicorn backend.app:app --host 127.0.0.1 --port 8000 --reload
   ```

3. Access the interface:
   Open `http://127.0.0.1:8000` in a web browser.

---

### Method 3: Containerized Deployment (Docker)

1. Build the container image:
   ```bash
   docker build -t cgl-regulation-explorer .
   ```

2. Run the container:
   ```bash
   docker run -d -p 8000:8000 --name cgl-app cgl-regulation-explorer
   ```

---

## REST API Reference

| Endpoint | Method | Description |
| :--- | :--- | :--- |
| `/api/gene/profile/{gene_id}` | `GET` | Queries 360-degree profile data for a target gene |
| `/api/genomic_tracks/{gene_id}` | `GET` | Queries 5-track genomic coordinates and peak data |
| `/api/imodulon/condition` | `GET` | Queries iModulon activity matrix for a specified condition |
| `/api/intervention-targets` | `GET` | Returns ranked metabolic engineering target priorities |
| `/api/graph/motifs` | `GET` | Identifies feed-forward loop (FFL) network motifs |
| `/api/ai/engineering_command` | `POST` | Executes AI engineering commands (`/design-crispri`, `/solve-bottleneck`, `/promoter-engineering`) |
| `/api/check-update` | `GET` | Checks local application version status |
| `/api/provenance` | `GET` | Reports schema version, source hashes, release manifest, and interpretation limits |

---

## Primary Data Sources

1. **CollectTF & RegPrecise**: Experimentally validated transcription factor binding sites and PWM matrices.
2. **Abasy Atlas**: System-level regulatory roles and pleiotropy classifications.
3. **iModulonDB**: Independent component analysis compendium based on *C. glutamicum* RNA-seq datasets.
4. **STRING v12**: Protein-protein interaction scoring network.
5. **NCBI RefSeq (NC_003450.3)**: Genomic coordinates and strand annotations.
6. **iCW773 GEM**: Genome-scale metabolic model for *C. glutamicum*.

---

## Testing and Quality Assurance

The codebase includes automated test suites covering backend API endpoints, database operations, and algorithmic logic:

Run the test suite:
```bash
python -m pytest
```

Pull requests and pushes to `main` run the suite on Python 3.10, 3.11, and
3.13. Scientific output integrity can be audited separately with
`python scripts/analysis/validate_scientific_outputs.py`. See
[Data provenance and scientific interpretation](docs/DATA_PROVENANCE.md) for
the evidence classes, reproducibility metadata, and validation limits.

---

## Citation and License

If this platform contributes to research publications, please cite:

```bibtex
@software{cgl_regulation_2026,
  author = {Gao, Dandan and DeepMind Agentic Coding Team},
  title = {Cgl Regulation Explorer: A Systems Biology and Synthetic Biology Platform for Corynebacterium glutamicum},
  year = {2026},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/gaodandan-ai/cgl_regulation}},
  version = {v1.3.1}
}
```

This software is released under the **MIT License**. Refer to the [LICENSE](LICENSE) file for terms.
