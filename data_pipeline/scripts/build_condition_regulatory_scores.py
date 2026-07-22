#!/usr/bin/env python3
"""Build a validated iron-response prototype for condition-specific TF edges.

DtxR uses its ICA activity, oriented with the two delta-dtxR contrasts. HrrA
has no separate iModulon and therefore uses a signed projection of its curated
activation/repression targets. Missing measurements reduce completeness; they
are never interpreted as negative evidence.
"""

from __future__ import annotations

import json
import math
import sqlite3
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[2]
DB_PATH = ROOT_DIR / "data" / "reference" / "cgl_regulation.db"
RUN_ID = "iron_regulon_v1"
METHOD_VERSION = "condition_edge_score_v1"
IRON_TERMS = ("feso4", "iron", "heme", "haemin", "hrra", "hrrsa", "dtxr")
TF_CONFIG = {
    "DtxR": {"locus": "cg2103", "activity_method": "ICA_iM_59_DtxR"},
    "HrrA": {"locus": "cg3247", "activity_method": "signed_regulon_projection"},
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def role_from_values(values: set[str]) -> str:
    tokens = {token.strip().upper() for value in values for token in value.split(";") if token.strip()}
    has_a, has_r = "A" in tokens, "R" in tokens
    if has_a and not has_r:
        return "A"
    if has_r and not has_a:
        return "R"
    return "Dual_or_unknown"


def load_edges(cursor: sqlite3.Cursor, locus: str) -> list[dict[str, object]]:
    grouped: dict[str, dict[str, object]] = {}
    for target, name, role, confidence, label, source, edge_id in cursor.execute("""
        SELECT target_locus, target_name, regulation_role, confidence_score,
               confidence_label, primary_source, edge_id
        FROM regulatory_edges
        WHERE regulator_locus=?
    """, (locus,)):
        if not target:
            continue
        item = grouped.setdefault(target.lower(), {
            "target_locus": target.lower(), "target_name": name or target,
            "roles": set(), "base_confidence": 0.0, "confidence_label": "LOW",
            "sources": set(), "source_edge_ids": set(),
        })
        item["roles"].add(role or "")
        if float(confidence or 0.0) >= float(item["base_confidence"]):
            item["base_confidence"] = float(confidence or 0.0)
            item["confidence_label"] = label or "LOW"
        if source:
            item["sources"].add(source)
        if edge_id:
            item["source_edge_ids"].add(edge_id)
    output = []
    for item in grouped.values():
        item["role"] = role_from_values(item.pop("roles"))
        item["sources"] = "; ".join(sorted(item["sources"]))
        item["source_edge_ids"] = sorted(item["source_edge_ids"])
        output.append(item)
    return sorted(output, key=lambda item: str(item["target_locus"]))


def signed_projection(edges: list[dict[str, object]], gene_means: dict[str, float]) -> tuple[float | None, int]:
    signed = []
    for edge in edges:
        value = gene_means.get(str(edge["target_locus"]))
        if value is None:
            continue
        if edge["role"] == "A":
            signed.append(value)
        elif edge["role"] == "R":
            signed.append(-value)
    return (float(np.median(signed)), len(signed)) if signed else (None, 0)


def dynamic_score(
    base: float,
    activity: float | None,
    expression: float | None,
    replicate_consistency: float | None,
    direction_consistency: float | None,
) -> tuple[float, float, float]:
    components: list[tuple[float, float]] = []
    if activity is not None:
        components.append((0.30, math.tanh(abs(activity) / 2.0)))
    if expression is not None:
        components.append((0.25, math.tanh(abs(expression))))
    if replicate_consistency is not None:
        components.append((0.20, replicate_consistency))
    if direction_consistency is not None:
        components.append((0.25, direction_consistency))
    completeness = sum(weight for weight, _ in components)
    support = sum(weight * value for weight, value in components) / completeness if completeness else 0.0
    final = max(0.0, min(1.0, base * (0.5 + 0.5 * support)))
    return support, completeness, final


def empirical_enrichment(
    target_genes: set[str], gene_means: dict[str, float], seed: int, permutations: int = 500,
) -> tuple[float | None, float | None, int, int]:
    observed_values = [abs(gene_means[g]) for g in target_genes if g in gene_means]
    background = np.array([abs(value) for value in gene_means.values()], dtype=float)
    size = len(observed_values)
    if size < 5 or len(background) < size:
        return None, None, size, permutations
    observed = float(np.mean(observed_values))
    background_mean = float(np.mean(background))
    rng = np.random.default_rng(seed)
    null = np.array([
        float(np.mean(rng.choice(background, size=size, replace=False)))
        for _ in range(permutations)
    ])
    p_value = float((1 + np.sum(null >= observed)) / (permutations + 1))
    return observed - background_mean, p_value, size, permutations


def add_bh_q_values(rows: list[tuple]) -> list[tuple]:
    """Insert a global Benjamini-Hochberg q-value after the p-value column."""
    valid = [(index, float(row[11])) for index, row in enumerate(rows) if row[11] is not None]
    ordered = sorted(valid, key=lambda item: item[1])
    q_values: dict[int, float] = {}
    running = 1.0
    total = len(ordered)
    for rank_index in range(total - 1, -1, -1):
        row_index, p_value = ordered[rank_index]
        rank = rank_index + 1
        running = min(running, p_value * total / rank)
        q_values[row_index] = min(1.0, running)
    adjusted = []
    for index, row in enumerate(rows):
        q_value = q_values.get(index)
        enrichment = row[10]
        if q_value is None:
            validation = "insufficient_expression"
        elif q_value <= 0.10 and enrichment is not None and enrichment > 0:
            validation = "response_enriched_fdr10"
        else:
            validation = "not_enriched"
        adjusted.append(row[:12] + (q_value,) + row[12:15] + (validation,))
    return adjusted


def build_condition_regulatory_scores(db_path: Path = DB_PATH) -> None:
    print(f"Building HrrA/DtxR condition-specific regulatory scores in: {db_path}")
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys=ON")
    cursor = conn.cursor()
    cursor.executescript("""
        DROP VIEW IF EXISTS v_iron_condition_top_edges;
        DROP VIEW IF EXISTS v_iron_regulon_response;
        DROP TABLE IF EXISTS condition_regulon_summary;
        DROP TABLE IF EXISTS condition_edge_scores;
        DROP TABLE IF EXISTS condition_tf_activity;
        DROP TABLE IF EXISTS condition_analysis_runs;
        CREATE TABLE condition_analysis_runs (
            run_id TEXT PRIMARY KEY, method_version TEXT NOT NULL,
            scope TEXT NOT NULL, created_at TEXT NOT NULL,
            parameters_json TEXT NOT NULL, notes TEXT NOT NULL
        );
        CREATE TABLE condition_tf_activity (
            run_id TEXT NOT NULL REFERENCES condition_analysis_runs(run_id),
            comparison_id TEXT NOT NULL REFERENCES condition_comparisons(comparison_id),
            tf_locus TEXT NOT NULL, tf_name TEXT NOT NULL,
            activity_score REAL, normalized_activity REAL,
            activity_direction TEXT NOT NULL, activity_method TEXT NOT NULL,
            source_imodulon_id TEXT, contributing_target_count INTEGER NOT NULL,
            expression_sample_count INTEGER NOT NULL, confidence_label TEXT NOT NULL,
            details_json TEXT NOT NULL,
            PRIMARY KEY(run_id, comparison_id, tf_locus)
        );
        CREATE TABLE condition_edge_scores (
            run_id TEXT NOT NULL REFERENCES condition_analysis_runs(run_id),
            comparison_id TEXT NOT NULL REFERENCES condition_comparisons(comparison_id),
            tf_locus TEXT NOT NULL, tf_name TEXT NOT NULL,
            target_locus TEXT NOT NULL, target_name TEXT,
            regulation_role TEXT NOT NULL, base_confidence REAL NOT NULL,
            tf_activity REAL, target_expression_mean REAL,
            target_expression_sd REAL, expression_sample_count INTEGER NOT NULL,
            replicate_consistency REAL, direction_consistency REAL,
            dynamic_support REAL NOT NULL, evidence_completeness REAL NOT NULL,
            condition_score REAL NOT NULL, support_state TEXT NOT NULL,
            source_edge_ids_json TEXT NOT NULL, details_json TEXT NOT NULL,
            PRIMARY KEY(run_id, comparison_id, tf_locus, target_locus)
        );
        CREATE TABLE condition_regulon_summary (
            run_id TEXT NOT NULL REFERENCES condition_analysis_runs(run_id),
            comparison_id TEXT NOT NULL REFERENCES condition_comparisons(comparison_id),
            condition_label TEXT NOT NULL, tf_locus TEXT NOT NULL, tf_name TEXT NOT NULL,
            activity_score REAL, scored_edge_count INTEGER NOT NULL,
            responsive_edge_count INTEGER NOT NULL, direction_consistent_count INTEGER NOT NULL,
            mean_condition_score REAL, response_enrichment REAL,
            empirical_p_value REAL, empirical_q_value REAL,
            tested_target_count INTEGER NOT NULL,
            permutation_count INTEGER NOT NULL, top_targets_json TEXT NOT NULL,
            validation_status TEXT NOT NULL,
            PRIMARY KEY(run_id, comparison_id, tf_locus)
        );
    """)

    iron_rows = []
    for comparison_id, source_id, label in cursor.execute("""
        SELECT comparison_id, source_record_id, source_label
        FROM omics_comparison_map
        WHERE omics_type='imodulon_activity'
    """):
        lower = (label or "").lower()
        if any(term in lower for term in IRON_TERMS):
            iron_rows.append((comparison_id, str(source_id), label))
    comparison_ids = sorted({row[0] for row in iron_rows})
    if not comparison_ids:
        raise RuntimeError("No iron-related condition comparisons were found")

    placeholders = ",".join("?" for _ in comparison_ids)
    sample_ids_by_comparison: dict[str, list[str]] = defaultdict(list)
    for comparison_id, sample_id in cursor.execute(f"""
        SELECT comparison_id, source_record_id FROM omics_comparison_map
        WHERE omics_type='transcriptomics' AND comparison_id IN ({placeholders})
    """, comparison_ids):
        sample_ids_by_comparison[comparison_id].append(sample_id)

    gene_values: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for comparison_id, sample_id, locus, value in cursor.execute(f"""
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

    edges_by_tf = {name: load_edges(cursor, cfg["locus"]) for name, cfg in TF_CONFIG.items()}
    ica_by_comparison: dict[str, float] = {}
    label_by_comparison = {comparison_id: label for comparison_id, _, label in iron_rows}
    for comparison_id, source_id, _ in iron_rows:
        row = cursor.execute("""
            SELECT activity_score FROM imodulon_condition_activities
            WHERE imodulon_id='iM_59_DtxR' AND sample_id=? LIMIT 1
        """, (source_id,)).fetchone()
        if row:
            ica_by_comparison[comparison_id] = float(row[0])
    knockout_values = [
        score for comparison_id, score in ica_by_comparison.items()
        if "dtxr" in label_by_comparison[comparison_id].lower()
        and ("delta" in label_by_comparison[comparison_id].lower()
             or "∆" in label_by_comparison[comparison_id])
    ]
    orientation = -1.0 if knockout_values and float(np.mean(knockout_values)) > 0 else 1.0

    cursor.execute("INSERT INTO condition_analysis_runs VALUES (?,?,?,?,?,?)", (
        RUN_ID, METHOD_VERSION, "HrrA/DtxR iron and heme response prototype", utc_now(),
        json.dumps({
            "iron_terms": IRON_TERMS, "permutations": 500,
            "dtxr_orientation_multiplier": orientation,
            "score_formula": "base_confidence * (0.5 + 0.5 * weighted_dynamic_support)",
        }, ensure_ascii=False),
        "DtxR activity is ICA-derived and knockout-oriented; HrrA is a signed target projection. Missing values reduce completeness.",
    ))

    activity_rows = []
    edge_rows = []
    summary_rows = []
    for comparison_index, comparison_id in enumerate(comparison_ids):
        label = label_by_comparison[comparison_id]
        means = gene_means.get(comparison_id, {})
        sample_count = len(sample_ids_by_comparison.get(comparison_id, []))
        for tf_name, cfg in TF_CONFIG.items():
            edges = edges_by_tf[tf_name]
            if tf_name == "DtxR":
                raw_activity = ica_by_comparison.get(comparison_id)
                activity = raw_activity * orientation if raw_activity is not None else None
                contributing = 0
                method = "ICA_iM_59_DtxR_knockout_oriented"
                imodulon_id = "iM_59_DtxR"
                confidence = "high" if activity is not None else "unavailable"
                activity_details = {"raw_ica_score": raw_activity, "orientation_multiplier": orientation}
            else:
                activity, contributing = signed_projection(edges, means)
                method = "median_signed_clear-role_targets"
                imodulon_id = None
                confidence = "medium" if contributing >= 5 else "low" if contributing else "unavailable"
                activity_details = {"clear_role_targets_used": contributing, "circularity_warning": "target-derived activity is not independent validation"}
            normalized_activity = math.tanh(activity / 2.0) if activity is not None else None
            direction = "increased" if activity is not None and activity > 0 else "decreased" if activity is not None and activity < 0 else "neutral_or_unavailable"
            activity_rows.append((
                RUN_ID, comparison_id, cfg["locus"], tf_name, activity,
                normalized_activity, direction, method, imodulon_id,
                contributing, sample_count, confidence,
                json.dumps(activity_details, ensure_ascii=False),
            ))

            scored_for_summary = []
            direction_consistent_count = 0
            responsive_count = 0
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
                details = {
                    "base_sources": edge["sources"],
                    "activity_strength": activity_strength,
                    "expression_strength": expression_strength,
                    "unknown_direction_not_penalized": edge["role"] == "Dual_or_unknown",
                }
                edge_rows.append((
                    RUN_ID, comparison_id, cfg["locus"], tf_name, target,
                    edge["target_name"], edge["role"], edge["base_confidence"],
                    activity, expression, expression_sd, len(values), consistency,
                    direction_consistency, support, completeness, condition_score,
                    state, json.dumps(edge["source_edge_ids"], ensure_ascii=False),
                    json.dumps(details, ensure_ascii=False),
                ))
                scored_for_summary.append((condition_score, target, edge["target_name"], expression, state))

            target_set = {
                str(edge["target_locus"]) for edge in edges
                if edge["role"] in {"A", "R"}
            } or {str(edge["target_locus"]) for edge in edges}
            enrichment, p_value, tested_count, permutations = empirical_enrichment(
                target_set, means, seed=20260722 + comparison_index * 17 + (0 if tf_name == "DtxR" else 1),
            )
            top_targets = sorted(
                scored_for_summary,
                key=lambda item: (item[0], abs(item[3]) if item[3] is not None else -1),
                reverse=True,
            )[:15]
            mean_score = float(np.mean([item[0] for item in scored_for_summary])) if scored_for_summary else None
            summary_rows.append((
                RUN_ID, comparison_id, label, cfg["locus"], tf_name, activity,
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
    cursor.executescript("""
        CREATE INDEX idx_condition_tf_activity_tf ON condition_tf_activity(tf_locus, comparison_id);
        CREATE INDEX idx_condition_edge_tf ON condition_edge_scores(tf_locus, comparison_id, condition_score DESC);
        CREATE INDEX idx_condition_edge_target ON condition_edge_scores(target_locus, comparison_id);
        CREATE INDEX idx_condition_edge_state ON condition_edge_scores(support_state, condition_score DESC);
        CREATE VIEW v_iron_condition_top_edges AS
        SELECT s.condition_label, e.comparison_id, e.tf_name, e.tf_locus,
               e.target_name, e.target_locus, e.regulation_role,
               e.tf_activity, e.target_expression_mean, e.replicate_consistency,
               e.condition_score, e.support_state, e.evidence_completeness
        FROM condition_edge_scores e
        JOIN condition_regulon_summary s
          ON s.run_id=e.run_id AND s.comparison_id=e.comparison_id AND s.tf_locus=e.tf_locus
        WHERE e.run_id='iron_regulon_v1';
        CREATE VIEW v_iron_regulon_response AS
        SELECT s.*, a.normalized_activity, a.activity_direction,
               a.activity_method, a.confidence_label AS activity_confidence
        FROM condition_regulon_summary s
        JOIN condition_tf_activity a
          ON a.run_id=s.run_id AND a.comparison_id=s.comparison_id AND a.tf_locus=s.tf_locus
        WHERE s.run_id='iron_regulon_v1';
    """)
    conn.commit()
    print(
        f"SUCCESS: scored {len(comparison_ids)} iron-related contrasts, "
        f"{len(activity_rows)} TF activities and {len(edge_rows)} condition-specific edges"
    )
    conn.close()


if __name__ == "__main__":
    build_condition_regulatory_scores()
