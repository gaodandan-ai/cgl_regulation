#!/usr/bin/env python3
"""Build condition-specific ArnR/GlxR/HrrA scores for oxygen limitation.

No matching independent iModulon activity is available for these regulators in
the integrated dataset. Activities are therefore signed projections of curated
clear-role targets and are explicitly labelled as target-derived.
"""

from __future__ import annotations

import json
import math
import sqlite3
from collections import defaultdict
from pathlib import Path

import numpy as np

from data_pipeline.scripts.build_condition_regulatory_scores import (
    DB_PATH,
    METHOD_VERSION,
    add_bh_q_values,
    dynamic_score,
    empirical_enrichment,
    load_edges,
    signed_projection,
    utc_now,
)


RUN_ID = "oxygen_regulon_v1"
OXYGEN_TERMS = ("anaerob", "microaerob", "oxygen")
TF_CONFIG = {
    "ArnR": {"locus": "cg1340"},
    "GlxR": {"locus": "cg0350"},
    "HrrA": {"locus": "cg3247"},
}


def build_projection_regulatory_module(
    db_path: Path,
    *,
    run_id: str,
    condition_terms: tuple[str, ...],
    tf_config: dict[str, dict[str, str]],
    scope: str,
    view_prefix: str,
    source_record_ids: set[str] | None = None,
    imodulon_config: dict[str, dict[str, object]] | None = None,
) -> None:
    """Build one clear-role target-projection module in the shared schema."""
    print(f"Building {scope} scores in: {db_path}")
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys=ON")
    cursor = conn.cursor()

    module_rows = []
    for comparison_id, source_id, label in cursor.execute("""
        SELECT comparison_id, source_record_id, source_label
        FROM omics_comparison_map
        WHERE omics_type='imodulon_activity'
    """):
        lower = (label or "").lower()
        selected = (
            str(source_id) in source_record_ids if source_record_ids is not None
            else any(term in lower for term in condition_terms)
        )
        if selected:
            module_rows.append((comparison_id, str(source_id), label))
    comparison_ids = sorted({row[0] for row in module_rows})
    if not comparison_ids:
        conn.close()
        raise RuntimeError(f"No condition comparisons were found for {run_id}")

    placeholders = ",".join("?" for _ in comparison_ids)
    sample_ids_by_comparison: dict[str, list[str]] = defaultdict(list)
    for comparison_id, sample_id in cursor.execute(f"""
        SELECT comparison_id, source_record_id FROM omics_comparison_map
        WHERE omics_type='transcriptomics' AND comparison_id IN ({placeholders})
    """, comparison_ids):
        sample_ids_by_comparison[comparison_id].append(sample_id)

    gene_values: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for comparison_id, _sample_id, locus, value in cursor.execute(f"""
        SELECT m.comparison_id, e.sample_id, lower(p.resolved_locus), e.processed_value
        FROM omics_comparison_map m
        JOIN expression_values e ON e.sample_id=m.source_record_id
        JOIN expression_probes p ON p.probe_id=e.probe_id
        WHERE m.omics_type='transcriptomics'
          AND m.comparison_id IN ({placeholders})
          AND p.resolved_locus IS NOT NULL
    """, comparison_ids):
        gene_values[comparison_id][locus].append(float(value))
    gene_means = {
        comparison_id: {gene: float(np.mean(values)) for gene, values in genes.items()}
        for comparison_id, genes in gene_values.items()
    }

    label_by_comparison = {comparison_id: label for comparison_id, _, label in module_rows}
    source_by_comparison = {comparison_id: source_id for comparison_id, source_id, _ in module_rows}
    edges_by_tf = {name: load_edges(cursor, cfg["locus"]) for name, cfg in tf_config.items()}

    # Make standalone reruns idempotent without disturbing the iron module.
    cursor.execute("DELETE FROM condition_regulon_summary WHERE run_id=?", (run_id,))
    cursor.execute("DELETE FROM condition_edge_scores WHERE run_id=?", (run_id,))
    cursor.execute("DELETE FROM condition_tf_activity WHERE run_id=?", (run_id,))
    cursor.execute("DELETE FROM condition_analysis_runs WHERE run_id=?", (run_id,))
    cursor.execute("INSERT INTO condition_analysis_runs VALUES (?,?,?,?,?,?)", (
        run_id, METHOD_VERSION, scope, utc_now(),
        json.dumps({
            "condition_terms": condition_terms,
            "source_record_ids": sorted(source_record_ids) if source_record_ids is not None else None,
            "permutations": 500,
            "activity_method": "median signed clear-role target projection",
            "score_formula": "base_confidence * (0.5 + 0.5 * weighted_dynamic_support)",
        }, ensure_ascii=False),
        "All TF activities are target-derived projections, not independent validation. Missing values reduce completeness.",
    ))

    activity_rows = []
    edge_rows = []
    summary_rows = []
    for comparison_index, comparison_id in enumerate(comparison_ids):
        label = label_by_comparison[comparison_id]
        means = gene_means.get(comparison_id, {})
        sample_count = len(sample_ids_by_comparison.get(comparison_id, []))
        for tf_index, (tf_name, cfg) in enumerate(tf_config.items()):
            edges = edges_by_tf[tf_name]
            lower_label = label.lower()
            regulator_knockout = any(token in lower_label for token in (
                f"∆{tf_name.lower()}", f"delta-{tf_name.lower()}", f"delta {tf_name.lower()}",
            ))
            if regulator_knockout:
                activity, contributing = None, 0
                confidence = "unavailable"
                activity_method = "regulator_knockout_no_activity_inference"
                source_imodulon_id = None
                activity_details = {"regulator_knockout": True}
            elif imodulon_config and tf_name in imodulon_config:
                im_cfg = imodulon_config[tf_name]
                source_imodulon_id = str(im_cfg["imodulon_id"])
                raw_row = cursor.execute("""
                    SELECT activity_score FROM imodulon_condition_activities
                    WHERE imodulon_id=? AND sample_id=? LIMIT 1
                """, (source_imodulon_id, source_by_comparison[comparison_id])).fetchone()
                raw_activity = float(raw_row[0]) if raw_row else None
                orientation = float(im_cfg.get("orientation", 1.0))
                activity = raw_activity * orientation if raw_activity is not None else None
                contributing = 0
                confidence = str(im_cfg.get("confidence", "medium")) if activity is not None else "unavailable"
                activity_method = f"ICA_{source_imodulon_id}_regulon_oriented"
                activity_details = {
                    "raw_ica_score": raw_activity,
                    "orientation_multiplier": orientation,
                    "orientation_basis": im_cfg.get("orientation_basis", "configured"),
                    "regulon_f1": im_cfg.get("regulon_f1"),
                    "regulator_knockout": False,
                }
            else:
                activity, contributing = signed_projection(edges, means)
                confidence = "medium" if contributing >= 5 else "low" if contributing else "unavailable"
                activity_method = "median_signed_clear-role_targets"
                source_imodulon_id = None
                activity_details = {
                    "clear_role_targets_used": contributing,
                    "regulator_knockout": False,
                    "circularity_warning": "target-derived activity is not independent validation",
                }
            normalized_activity = math.tanh(activity / 2.0) if activity is not None else None
            direction = (
                "increased" if activity is not None and activity > 0
                else "decreased" if activity is not None and activity < 0
                else "neutral_or_unavailable"
            )
            activity_rows.append((
                run_id, comparison_id, cfg["locus"], tf_name, activity,
                normalized_activity, direction, activity_method, source_imodulon_id,
                contributing, sample_count, confidence,
                json.dumps(activity_details, ensure_ascii=False),
            ))

            scored_for_summary = []
            responsive_count = 0
            direction_consistent_count = 0
            for edge in edges:
                target = str(edge["target_locus"])
                values = gene_values.get(comparison_id, {}).get(target, [])
                expression = float(np.mean(values)) if values else None
                expression_sd = float(np.std(values, ddof=1)) if len(values) > 1 else None
                consistency = None
                if values:
                    positive = sum(value > 0 for value in values)
                    negative = sum(value < 0 for value in values)
                    consistency = max(positive, negative) / len(values)
                    consistency *= min(1.0, len(values) / 3.0)
                direction_consistency = None
                if activity is not None and expression is not None and activity != 0 and expression != 0:
                    if edge["role"] == "A":
                        direction_consistency = float(activity * expression > 0)
                    elif edge["role"] == "R":
                        direction_consistency = float(activity * expression < 0)
                support, completeness, condition_score = dynamic_score(
                    float(edge["base_confidence"]), activity, expression,
                    consistency, direction_consistency,
                )
                expression_strength = math.tanh(abs(expression)) if expression is not None else 0.0
                activity_strength = math.tanh(abs(activity) / 2.0) if activity is not None else 0.0
                if direction_consistency == 0.0 and expression_strength >= 0.3 and activity_strength >= 0.3:
                    state = "direction_conflict"
                    condition_score = min(condition_score, float(edge["base_confidence"]) * 0.55)
                elif condition_score >= 0.70 and expression_strength >= 0.2 and activity_strength >= 0.2:
                    state = "condition_supported"
                elif expression is None or activity is None:
                    state = "insufficient_dynamic_data"
                    condition_score = min(condition_score, float(edge["base_confidence"]) * 0.65)
                else:
                    state = "weak_context_support"
                responsive_count += int(expression_strength >= 0.3)
                direction_consistent_count += int(direction_consistency == 1.0)
                edge_rows.append((
                    run_id, comparison_id, cfg["locus"], tf_name, target,
                    edge["target_name"], edge["role"], edge["base_confidence"],
                    activity, expression, expression_sd, len(values), consistency,
                    direction_consistency, support, completeness, condition_score,
                    state, json.dumps(edge["source_edge_ids"], ensure_ascii=False),
                    json.dumps({
                        "base_sources": edge["sources"],
                        "activity_strength": activity_strength,
                        "expression_strength": expression_strength,
                        "unknown_direction_not_penalized": edge["role"] == "Dual_or_unknown",
                    }, ensure_ascii=False),
                ))
                scored_for_summary.append((condition_score, target, edge["target_name"], expression, state))

            target_set = {
                str(edge["target_locus"]) for edge in edges if edge["role"] in {"A", "R"}
            } or {str(edge["target_locus"]) for edge in edges}
            enrichment, p_value, tested_count, permutations = empirical_enrichment(
                target_set, means,
                seed=20260722 + comparison_index * 17 + tf_index,
            )
            top_targets = sorted(
                scored_for_summary,
                key=lambda item: (item[0], abs(item[3]) if item[3] is not None else -1),
                reverse=True,
            )[:15]
            mean_score = float(np.mean([item[0] for item in scored_for_summary])) if scored_for_summary else None
            summary_rows.append((
                run_id, comparison_id, label, cfg["locus"], tf_name, activity,
                len(scored_for_summary), responsive_count, direction_consistent_count,
                mean_score, enrichment, p_value, tested_count, permutations,
                json.dumps([
                    {"locus": item[1], "name": item[2], "expression": item[3], "score": item[0], "state": item[4]}
                    for item in top_targets
                ], ensure_ascii=False), "pending_fdr",
            ))

    cursor.executemany("INSERT INTO condition_tf_activity VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", activity_rows)
    cursor.executemany("INSERT INTO condition_edge_scores VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", edge_rows)
    summary_rows = add_bh_q_values(summary_rows)
    cursor.executemany("INSERT INTO condition_regulon_summary VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", summary_rows)
    cursor.executescript(f"""
        DROP VIEW IF EXISTS v_{view_prefix}_condition_top_edges;
        DROP VIEW IF EXISTS v_{view_prefix}_regulon_response;
        DROP VIEW IF EXISTS v_condition_regulation_top_edges;
        DROP VIEW IF EXISTS v_condition_regulon_response;
        CREATE VIEW v_condition_regulation_top_edges AS
        SELECT e.run_id, r.scope, s.condition_label, e.comparison_id, e.tf_name, e.tf_locus,
               e.target_name, e.target_locus, e.regulation_role,
               e.tf_activity, e.target_expression_mean, e.replicate_consistency,
               e.condition_score, e.support_state, e.evidence_completeness
        FROM condition_edge_scores e
        JOIN condition_regulon_summary s
          ON s.run_id=e.run_id AND s.comparison_id=e.comparison_id AND s.tf_locus=e.tf_locus
        JOIN condition_analysis_runs r ON r.run_id=e.run_id;
        CREATE VIEW v_condition_regulon_response AS
        SELECT r.scope, s.*, a.normalized_activity, a.activity_direction,
               a.activity_method, a.confidence_label AS activity_confidence
        FROM condition_regulon_summary s
        JOIN condition_tf_activity a
          ON a.run_id=s.run_id AND a.comparison_id=s.comparison_id AND a.tf_locus=s.tf_locus
        JOIN condition_analysis_runs r ON r.run_id=s.run_id;
        CREATE VIEW v_{view_prefix}_condition_top_edges AS
        SELECT * FROM v_condition_regulation_top_edges WHERE run_id='{run_id}';
        CREATE VIEW v_{view_prefix}_regulon_response AS
        SELECT * FROM v_condition_regulon_response WHERE run_id='{run_id}';
    """)
    conn.commit()
    print(
        f"SUCCESS: scored {len(comparison_ids)} contrasts for {run_id}, "
        f"{len(activity_rows)} TF activities and {len(edge_rows)} condition-specific edges"
    )
    conn.close()


def build_oxygen_regulatory_scores(db_path: Path = DB_PATH) -> None:
    build_projection_regulatory_module(
        db_path,
        run_id=RUN_ID,
        condition_terms=OXYGEN_TERMS,
        tf_config=TF_CONFIG,
        scope="ArnR/GlxR/HrrA oxygen limitation response",
        view_prefix="oxygen",
    )


if __name__ == "__main__":
    build_oxygen_regulatory_scores()
