#!/usr/bin/env python3
"""Build AmtR/ArgR/LtbR nitrogen and amino-acid condition scores."""

from __future__ import annotations

from pathlib import Path

from data_pipeline.scripts.build_condition_regulatory_scores import DB_PATH
from data_pipeline.scripts.build_oxygen_regulatory_scores import build_projection_regulatory_module


RUN_ID = "nitrogen_regulon_v1"
SOURCE_RECORD_IDS = {
    "8.0",    # glutamine addition under ammonium/urea depletion
    "9.0",    # glutamine under ammonium/urea depletion vs urea depletion
    "10.0",   # oxoproline addition to glutamine
    "11.0",   # oxoproline vs glutamine
    "12.0",   # glutamate vs glutamine
    "15.0",   # ammonium/urea depletion in GABA medium
    "79.0",   # urea/MOPS removal in bioreactor
    "80.0",   # urea removal in BioLector
    "183.0",  # ammonium concentration in delta-odhI
    "184.0",  # ammonium concentration in delta-odhI with glucose/glutamine
}
TF_CONFIG = {
    "AmtR": {"locus": "cg0986"},
    "ArgR": {"locus": "cg1585"},
    "LtbR": {"locus": "cg1486"},
}
IMODULON_CONFIG = {
    "AmtR": {
        "imodulon_id": "iM_60_AmtR", "orientation": -1.0,
        "confidence": "medium", "regulon_f1": 0.4444444444444445,
        "orientation_basis": "negative correlation with signed clear-role regulon projection",
    },
    "ArgR": {
        "imodulon_id": "iM_41_ArgR", "orientation": -1.0,
        "confidence": "medium", "regulon_f1": 0.5454545454545454,
        "orientation_basis": "negative correlation with signed clear-role regulon projection",
    },
}


def build_nitrogen_regulatory_scores(db_path: Path = DB_PATH) -> None:
    build_projection_regulatory_module(
        db_path,
        run_id=RUN_ID,
        condition_terms=(),
        tf_config=TF_CONFIG,
        scope="AmtR/ArgR/LtbR nitrogen-source and amino-acid response",
        view_prefix="nitrogen",
        source_record_ids=SOURCE_RECORD_IDS,
        imodulon_config=IMODULON_CONFIG,
    )


if __name__ == "__main__":
    build_nitrogen_regulatory_scores()
