#!/usr/bin/env python3
"""Build the read-only, query-focused database shipped with the public API."""

from __future__ import annotations

import argparse
import json
import os
import sqlite3
from pathlib import Path


PUBLIC_TABLES = {
    "abasy_roles",
    "biocyc_kegg_pathways",
    "brenda_kcat",
    "canonical_locus_map",
    "chebi_mappings",
    "cog_annotations",
    "collectf_tfbs",
    "condition_analysis_runs",
    "condition_edge_scores",
    "condition_regulon_summary",
    "condition_tf_activity",
    "essential_genes",
    "gene_coordinates",
    "gene_mappings",
    "gene_pathway_mappings",
    "imodulon_condition_activities",
    "imodulon_gene_weights",
    "imodulon_regulon_overlaps",
    "imodulons",
    "intervention_target_module_evidence",
    "intervention_target_scores",
    "literature_fts",
    "ncrna_target_interactions",
    "network_centrality",
    "network_edges_extended",
    "network_rewired_edges",
    "regulations",
    "rfam_ncrnas",
    "rhea_mappings",
    "string_interactions",
    "tf_families_effectors",
    "tf_gene_rf_scores",
    "tf_hierarchy_rankings",
}

PUBLIC_VIEWS = {
    "v_condition_regulation_top_edges",
    "v_condition_regulon_response",
    "v_gene_full_profile",
    "v_intervention_target_priorities",
    "v_metabolite_tf_feedback",
    "v_srna_competition_ranking",
}


def _quote(identifier: str) -> str:
    return '"' + identifier.replace('"', '""') + '"'


def build(source: Path, destination: Path) -> dict:
    if not source.is_file():
        raise FileNotFoundError(f"Source database not found: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    if temporary.exists():
        temporary.unlink()

    connection = sqlite3.connect(temporary)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=MEMORY")
    connection.execute("ATTACH DATABASE ? AS source", (str(source.resolve()),))

    objects = {
        (kind, name): sql
        for kind, name, sql in connection.execute(
            "SELECT type, name, sql FROM source.sqlite_master WHERE sql IS NOT NULL"
        )
    }
    missing = sorted(name for name in PUBLIC_TABLES if ("table", name) not in objects)
    if missing:
        raise RuntimeError(f"Required public tables are missing: {', '.join(missing)}")

    counts = {}
    with connection:
        for table in sorted(PUBLIC_TABLES):
            connection.execute(objects[("table", table)])
            columns = [row[1] for row in connection.execute(f"PRAGMA source.table_info({_quote(table)})")]
            column_list = ", ".join(_quote(column) for column in columns)
            connection.execute(
                f"INSERT INTO {_quote(table)} ({column_list}) "
                f"SELECT {column_list} FROM source.{_quote(table)}"
            )
            counts[table] = connection.execute(f"SELECT COUNT(*) FROM {_quote(table)}").fetchone()[0]

        for index_name, table_name, sql in connection.execute(
            "SELECT name, tbl_name, sql FROM source.sqlite_master "
            "WHERE type='index' AND sql IS NOT NULL ORDER BY name"
        ):
            if table_name in PUBLIC_TABLES:
                connection.execute(sql)

        for view in sorted(PUBLIC_VIEWS):
            sql = objects.get(("view", view))
            if sql:
                connection.execute(sql)

        source_version = connection.execute("PRAGMA source.user_version").fetchone()[0]
        connection.execute(f"PRAGMA user_version={int(source_version)}")

    connection.execute("DETACH DATABASE source")
    connection.execute("PRAGMA optimize")
    connection.execute("VACUUM")
    connection.close()
    os.replace(temporary, destination)

    manifest = {
        "source": source.name,
        "destination": destination.as_posix(),
        "bytes": destination.stat().st_size,
        "tables": len(PUBLIC_TABLES),
        "views": len(PUBLIC_VIEWS),
        "rows": counts,
    }
    destination.with_suffix(".manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=Path("data/reference/cgl_regulation.db"))
    parser.add_argument("--output", type=Path, default=Path("data/deploy/cgl_regulation_public.db"))
    args = parser.parse_args()
    print(json.dumps(build(args.source, args.output), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
