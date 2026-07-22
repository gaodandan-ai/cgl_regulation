#!/usr/bin/env python3
"""Build SigH/OxyR/MtrA oxidative, envelope and secretion-stress scores."""

from __future__ import annotations

from pathlib import Path

from data_pipeline.scripts.build_condition_regulatory_scores import DB_PATH
from data_pipeline.scripts.build_oxygen_regulatory_scores import build_projection_regulatory_module


RUN_ID = "stress_regulon_v1"
SOURCE_RECORD_IDS = {
    "61.0",   # cumene hydroperoxide
    "62.0",   # cinnamic acid
    "63.0",   # ferulic acid
    "64.0",   # coumaric acid without PCA
    "65.0",   # coumaric acid
    "66.0",   # coumaric acid with citrate
    "67.0",   # caffeic acid
    "68.0",   # hydroxyphenyl propionic acid
    "69.0",   # hydroxybenzoic acid
    "76.0",   # cutinase induction
    "77.0",   # cutinase secretion, 1 h
    "78.0",   # cutinase secretion, 4 h
    "151.0",  # ethambutol in delta-mtrAB
    "154.0",  # heat shock in delta-sigH
    "167.0",  # heat shock in delta-sigE
    "168.0",  # AmyE secretion in delta-sigE
    "211.0",  # AmyE secretion, delta-sigB vs WT
    "212.0",  # AmyE induction in delta-sigB
}
TF_CONFIG = {
    "SigH": {"locus": "cg0876"},
    "OxyR": {"locus": "cg2109"},
    "MtrA": {"locus": "cg0862"},
}


def build_stress_regulatory_scores(db_path: Path = DB_PATH) -> None:
    build_projection_regulatory_module(
        db_path,
        run_id=RUN_ID,
        condition_terms=(),
        tf_config=TF_CONFIG,
        scope="SigH/OxyR/MtrA oxidative, envelope, heat and secretion stress response",
        view_prefix="stress",
        source_record_ids=SOURCE_RECORD_IDS,
    )


if __name__ == "__main__":
    build_stress_regulatory_scores()
