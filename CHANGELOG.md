# Changelog

All notable changes to the Cgl Regulation Explorer project will be documented in this file.

---

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
