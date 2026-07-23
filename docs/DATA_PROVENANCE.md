# Data provenance and scientific interpretation

## What is recorded

The full SQLite database records source artifacts in `dataset_metadata`. Each
entry contains the repository-relative path, SHA-256 digest, byte size, record
count when available, import time, source version, and notes. Database schema
changes are recorded in `schema_migrations` and `PRAGMA user_version`.

The public, query-focused database is accompanied by
`data/deploy/cgl_regulation_public.manifest.json`, which records its byte size,
included tables and views, and table row counts. The runtime endpoint
`GET /api/provenance` exposes the metadata available in the active deployment
without exposing local filesystem paths.

## Evidence classes

- **Experimental**: direct or manually curated experimental evidence.
- **Computed**: a deterministic or statistical result derived from source data.
- **Predicted**: a model-based hypothesis that requires independent validation.
- **Unknown**: the evidence type or source could not be established.

These labels describe provenance, not correctness. In particular, missing
evidence is not evidence that a regulatory relationship is absent.

## Rebuilding and auditing

Build scripts live under `data_pipeline/scripts`. After rebuilding the full
database, finalize its schema and provenance metadata with
`data_pipeline/scripts/finalize_database.py`. Build the serverless subset with
`scripts/pipeline/build_public_database.py`.

Run the machine-readable integrity audit with:

```bash
python scripts/analysis/validate_scientific_outputs.py \
  --output analysis_output/scientific_validation.json
```

This audit checks schema versioning, source hashes, required result tables, and
score domains. It does **not** establish biological validity.

An independent binary gold-standard dataset can be supplied as a CSV containing
`score` values in `[0, 1]` and `label` values of `0` or `1`:

```bash
python scripts/analysis/validate_scientific_outputs.py \
  --gold-standard validation/held_out_edges.csv \
  --output analysis_output/scientific_validation.json
```

The report then includes AUROC, average precision, Brier score, expected
calibration error, and per-bin calibration observations. Alternative column
names can be selected with `--score-column` and `--label-column`. The validation
CSV must be genuinely independent of model training and feature construction;
the command cannot detect biological or publication-level data leakage.

## Validation still required before experimental use

Before treating a predicted edge or engineering target as actionable, compare
the relevant method against an independent held-out dataset and report the
sampling unit, split strategy, baseline, precision/recall or calibration error,
uncertainty interval, and known sources of leakage. Proposed genetic
interventions require expert review and experimental confirmation.
