#!/usr/bin/env python3
"""Build condition-specific central-carbon regulatory scores.

The selected contrasts change carbon-source composition rather than merely
sharing a glucose medium across genotypes. All activities are signed clear-role
target projections; the weakly annotated RamA iModulon remains available as
supporting data but is not promoted to an independent activity estimate.
"""

from __future__ import annotations

from pathlib import Path

from data_pipeline.scripts.build_condition_regulatory_scores import DB_PATH
from data_pipeline.scripts.build_oxygen_regulatory_scores import build_projection_regulatory_module


RUN_ID = "carbon_regulon_v1"
SOURCE_RECORD_IDS = {
    "13.0",   # GABA vs glucose
    "16.0",   # pyruvate vs lactate
    "17.0",   # lactate vs glucose
    "19.0",   # citrate vs glucose in LB
    "20.0",   # citrate vs glucose in CGXII
    "21.0",   # glucose+citrate vs glucose
    "22.0",   # sodium citrate vs glucose
    "29.0",   # glucose+methanol vs glucose
    "30.0",   # methanol addition without glucose
    "59.0",   # acetate vs 1% glucose
    "60.0",   # acetate vs 4% glucose
    "254.0",  # lactate vs glucose in delta-ramA
}
TF_CONFIG = {
    "GlxR": {"locus": "cg0350"},
    "RamA": {"locus": "cg2831"},
    "RamB": {"locus": "cg0444"},
    "SugR": {"locus": "cg2115"},
}


def build_carbon_regulatory_scores(db_path: Path = DB_PATH) -> None:
    build_projection_regulatory_module(
        db_path,
        run_id=RUN_ID,
        condition_terms=(),
        tf_config=TF_CONFIG,
        scope="GlxR/RamA/RamB/SugR carbon-source and central-metabolism response",
        view_prefix="carbon",
        source_record_ids=SOURCE_RECORD_IDS,
    )


if __name__ == "__main__":
    build_carbon_regulatory_scores()
