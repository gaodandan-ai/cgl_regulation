# Changelog

All notable changes to the Cgl Regulation Explorer project will be documented in this file.

---

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
