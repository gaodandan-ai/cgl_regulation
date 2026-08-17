# Changelog

All notable changes to the Cgl Regulation Explorer project will be documented in this file.

## [v1.5.0] - 2026-08-17
### Added
- **16:9 Multi-Row Proportional Hierarchy View**: Balanced multi-row layout for dense regulatory tiers, eliminating horizontal stretch and fitting widescreen displays naturally.
- **Intelligent sRNA Density Control & Full Locus Formatting**: Fixed locus ID regex truncation (`cgb_20715`, `ncgl1747.1`), added interactive active/connected sRNA filter with one-click full library expansion.
- **3-Column Seamless Gene Detail Fullscreen Workspace**: Restructured right-sidebar fullscreen portal into a structured 3-column composite dashboard (Identity & Annotations, Regulatory Core, Genomics & Proteomics) with zero empty whitespace gaps.
- **Containerized Intranet & Production Deployment Pipelines**: Added one-click intranet deployment scripts and Docker compose orchestrations.

---

## [v1.4.0] - 2026-07-28
### Added
- **Enhanced iModulon Transcriptional Regulatory Network**: Interactive subnetwork view with layout selector (Concentric, COSE, Circular), weight threshold filter, and hover node detail card.
- **Visual Regulon Alignment & Mechanism Analysis**: Visual metrics for precision, recall, F1 score, and automatic transcriptional regulatory rationale summary.
- **Regulon & Mode Details in Gene Membership**: Added TRN membership status badges (Regulon Member vs Module Discovery) and explicit TF regulation mode badges (Activation `+` / Repression `-` / Dual `+/-`).

### Removed
- **Redundant Non-Regulatory Features**: Removed GEM FBA metabolic reaction map (`imodulon-tab-reactions`) and knockout flux simulation / multi-iModulon comparison (`imodulon-tab-simulation`) from the iModulon Explorer modal to streamline user focus on transcriptional regulation.

---

## [v1.3.0] - 2026-07-22
### Added
- **Condition-aware regulation analysis** across iron, carbon, nitrogen, oxygen, and stress response modules, with condition harmonization and edge-level support scoring.
- **Intervention target priorities** combining cross-module evidence, essentiality risk, and engineering strategy classes.
- **Centralized SQLite data layer** covering regulatory, genomic, transcriptomic, proteomic, non-coding RNA, pathway, iModulon, and literature resources.
- **New analysis views** for gene profiles, allosteric feedback, sRNA competition, condition-specific iModulons, and target prioritization.
- **Lightweight public data API** and generated read-only deployment database for Vercel.
- **Deployment, security, graph, database, omics, CollecTF, and ncRNA regression tests**.

### Changed
- Expanded the ATCC 13032 integration and evidence-confidence pipeline with additional public datasets and cross-strain analyses.
- Reduced the Vercel function bundle by separating public data queries from desktop-only metabolic simulation.
- Improved desktop startup diagnostics, loopback-only binding, application health checks, and packaging exclusions.

### Security
- Removed public debug and traceback exposure, restricted CORS, added security response headers, and validated outbound AI endpoints against SSRF.
- Excluded local credentials, full databases, rebuild backups, raw downloads, and temporary analysis files from releases.

## [v0.6.0] - 2026-07-21
### Added
- **Auto-update Notification**: App checks for new versions on startup and shows a sky-blue toast in the header when an update is available, with a direct link to GitHub Releases.
- **GitHub Releases Distribution**: Automated CI/CD pipeline builds and publishes Windows installer (`cgl_setup.exe`) and portable executable (`cgl_regulation.exe`) on every version tag push.
- **PPI Expression Overlay**: RNA-seq log2FC coloring on protein–protein interaction network nodes with red↔blue gradient, LFC threshold slider, dataset switcher, and floating legend.
- **Centrality Bubble Chart Enhancements**: Configurable X/Y axes and bubble-size dimensions; gene highlight search; click-to-inspect gene detail panel with "Open in Gene Explorer" jump; essentiality filter pills (All / Essential Only / Non-Essential Only).
- **Sigma Cascade Network Type**: New "Sigma Factor Cascade (Sigma→TF→Target)" network query mode.
- **PPI Score Threshold Slider**: Adjustable minimum STRING PPI score (150–900, default 400) in Advanced Network panel.

### Changed
- **Nav Hierarchy**: Primary workflow buttons (Gene/TF, Hierarchy, Network Topology, iModulon, PPI, Pathway) visually distinguished from secondary analysis tools (Engineering, Advanced Analytics, Simulation) via a vertical divider and subdued styling.
- **Nav Order**: Reordered to Gene/TF → Hierarchy View → Network Topology → iModulon Explorer → PPI Explorer → Pathway View.
- **App Icon**: New sky-blue radial-gradient squircle icon with white "Cgl" text (programmatically generated with Pillow).
- **Script Organization**: Data pipeline scripts moved to `data_pipeline/scripts/`; icon generator at `scripts/gen_icon.py`.

## [v0.5.0] - 2026-07-15
### Added
- **GPR Rule Parser for rFBA**: Implemented a recursive descent parser that evaluates standard COBRApy `gene_reaction_rule` expressions under formal `min` (AND) and `max` (OR) logic, ensuring rate-limiting subunit dynamics and isozyme compensation are accurately modeled in dynamic simulations.
- **Thermodynamic-Stoichiometric Conflict Reporting**: Exposed formerly silent Safety Guard lock rollbacks in the UI dashboard card as "Thermo-Stoichiometric Conflicts" to highlight gaps between thermodynamic energy bounds and model feasibility.
- **Bilingual RAG Tokenization**: Added overlapping character bigram tokenization alongside English word boundaries in the TF-IDF search, optimizing literature retrieval accuracy for Chinese and mixed-language publications.
- **Dynamic FBA Speedups**: Replaced per-step context instantiation in `run_dynamic_rfba` by wrapping the simulation loop in a single model context, avoiding 23 deep-copy operations of the solver and boosting speed by 10-20x.

## [v0.4.0] - 2026-07-05
### Added
- **Abasy Atlas Integration**: Integrated systemic topology roles (Global Hubs, Modular Regulators, and Pathway Genes) for networks.
- **Structured RAG AI Summary**: Enhanced literature synthesis using formatted Markdown tables and interactive locus tag anchors linkable to Cytoscape layout coordinates.

## [v0.1.1] - 2026-06-29
### Added
- **iModulon UI Panel**: Integrated iModulon transcriptional module badges to display gene memberships dynamically.
- **TCS Signal Chain Panel**: Integrated Two-Component System context displaying signal flows, stimuli, and targeted pathways.
- **Sigma Factor Context Panel**: Integrated ECF consensus promoter details and standard binding region mappings.

### Fixed
- **Motif Prediction File Paths**: Corrected local filesystem references in `run_server.py` to target relocated data files under `data/reference/`.
- **Fetch Request Safety**: Added `.then(res => { if (!res.ok) throw ... })` checks in frontend API wrappers to cleanly handle non-JSON responses and prevent JSON decoding failures.
- **Error Format Compatibility**: Adapted parser to recognize both FastAPI HTTP detail fields and default server error responses.

## [v0.1.0] - 2026-06-27
### Added
- **Workflow-based Navigation**: Compact horizontal pill-shaped tab system (Gene/TF Explorer, Pathway View, Engineering Targets, Data & Model Quality, Examples, Release Notes, References).
- **Default Example Network**: Automated loading of cg0350 / whiB4 / sigH regulatory neighborhood on first startup.
- **Unified Regulatory Input Schema**: Integration of transcriptional TRN evidence database and predicted sRNA-mRNA interactions.
- **ML Edge Confidence Model**: Random Forest model predicting edge probability priority scores from multiple evidence columns.
- **Metabolic Model Mapping**: iCW773 genome-scale model adapter linking genes to reactions and subsystem pathways.
- **Enzyme-constrained Annotations**: Support for ecCGL1 properties (MW, kcat, EC numbers, parent reactions, UniProt IDs) in details tables.
- **Engineering Prioritization Ranks**: Scoring system for candidate transcription factors based on target count, regulation effects, and mapped reactions.
- **Data & Model Quality Dashboard**: Quantitative counts, confidence histograms, and mapping coverage rates (e.g. 13.1% network-to-model ratio).
- **Built-in Examples & Case Studies**: Glutamate biosynthesis, TCA upstream regulators, and amino acid prioritizations case runners with narrative coverage alerts.
- **References & Data Attribution Page**: Grouped citations list, usage disclaimers, and DOI placeholders.
- **JSON & CSV Export Utilities**: Save active networks, audit stats, case summaries, and target lists.
