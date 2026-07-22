#!/usr/bin/env python3
"""Rank cross-module engineering targets with auditable evidence components."""

from __future__ import annotations

import bisect
import json
import math
import re
import sqlite3
from collections import defaultdict
from pathlib import Path

from data_pipeline.scripts.build_condition_regulatory_scores import DB_PATH, utc_now


RUN_ID = "cross_module_priority_v1"
METHOD_VERSION = "intervention_priority_v1"


def normalize_name(value: object) -> str:
    text = re.sub(r"\s*\[Corynebacterium.*$", "", str(value or ""), flags=re.I)
    return re.sub(r"[^a-z0-9]+", " ", text.lower()).strip()


def percentile(sorted_values: list[float], value: float) -> float:
    if not sorted_values:
        return 0.0
    return bisect.bisect_right(sorted_values, value) / len(sorted_values)


def build_intervention_priorities(db_path: Path = DB_PATH) -> None:
    print(f"Building cross-module intervention priorities in: {db_path}")
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys=ON")
    cursor = conn.cursor()
    cursor.executescript("""
        DROP VIEW IF EXISTS v_intervention_target_priorities;
        DROP TABLE IF EXISTS intervention_target_module_evidence;
        DROP TABLE IF EXISTS intervention_target_scores;
        DROP TABLE IF EXISTS proteomics_gene_mapping;
        DROP TABLE IF EXISTS intervention_priority_runs;
        CREATE TABLE intervention_priority_runs (
            run_id TEXT PRIMARY KEY, method_version TEXT NOT NULL,
            created_at TEXT NOT NULL, parameters_json TEXT NOT NULL, notes TEXT NOT NULL
        );
        CREATE TABLE proteomics_gene_mapping (
            protein_accession TEXT PRIMARY KEY,
            gene_locus TEXT NOT NULL, gene_name TEXT,
            mapping_method TEXT NOT NULL, mapping_confidence TEXT NOT NULL,
            unique_peptide_count INTEGER NOT NULL, transition_count INTEGER NOT NULL,
            details_json TEXT NOT NULL
        );
        CREATE TABLE intervention_target_scores (
            run_id TEXT NOT NULL REFERENCES intervention_priority_runs(run_id),
            target_locus TEXT NOT NULL, target_name TEXT, product TEXT,
            module_count INTEGER NOT NULL, condition_count INTEGER NOT NULL,
            regulator_count INTEGER NOT NULL, edge_observation_count INTEGER NOT NULL,
            supported_count INTEGER NOT NULL, conflict_count INTEGER NOT NULL,
            significant_context_count INTEGER NOT NULL,
            mean_condition_score REAL NOT NULL, max_condition_score REAL NOT NULL,
            mean_evidence_completeness REAL NOT NULL, direction_agreement REAL,
            pathway_count INTEGER NOT NULL, proteomics_detected INTEGER NOT NULL,
            unique_peptide_count INTEGER NOT NULL, variant_count INTEGER NOT NULL,
            phenotype_context_available INTEGER NOT NULL,
            essentiality_status TEXT NOT NULL,
            pagerank_percentile REAL NOT NULL, betweenness_percentile REAL NOT NULL,
            evidence_score REAL NOT NULL, systems_impact_score REAL NOT NULL,
            engineering_tractability_score REAL NOT NULL, risk_score REAL NOT NULL,
            priority_score REAL NOT NULL, evidence_grade TEXT NOT NULL,
            strategy_class TEXT NOT NULL, rationale_json TEXT NOT NULL,
            PRIMARY KEY(run_id, target_locus)
        );
        CREATE TABLE intervention_target_module_evidence (
            run_id TEXT NOT NULL, target_locus TEXT NOT NULL, module_run_id TEXT NOT NULL,
            condition_count INTEGER NOT NULL, regulator_count INTEGER NOT NULL,
            supported_count INTEGER NOT NULL, conflict_count INTEGER NOT NULL,
            significant_context_count INTEGER NOT NULL,
            mean_score REAL NOT NULL, max_score REAL NOT NULL,
            PRIMARY KEY(run_id, target_locus, module_run_id),
            FOREIGN KEY(run_id, target_locus)
                REFERENCES intervention_target_scores(run_id, target_locus)
        );
    """)
    cursor.execute("INSERT INTO intervention_priority_runs VALUES (?,?,?,?,?)", (
        RUN_ID, METHOD_VERSION, utc_now(),
        json.dumps({
            "evidence_weights": {
                "module_breadth": 0.20, "condition_breadth": 0.15,
                "regulator_breadth": 0.10, "mean_condition_score": 0.20,
                "supported_fraction": 0.15, "evidence_completeness": 0.10,
                "fdr_contexts": 0.10,
            },
            "priority_formula": "0.50 evidence + 0.25 systems impact + 0.25 engineering tractability",
            "proteomics_mapping": "unique exact normalized protein product name only",
        }, ensure_ascii=False),
        "Absence from the curated essential set is not proof of non-essentiality. Strategy classes require experimental validation and a product-specific objective.",
    ))

    genes: dict[str, tuple[str | None, str | None]] = {}
    product_index: dict[str, list[str]] = defaultdict(list)
    for locus, name, product in cursor.execute("""
        SELECT lower(canonical_cg), gene_name, product FROM genes WHERE canonical_cg IS NOT NULL
    """):
        genes[locus] = (name, product)
        for candidate in (name, product):
            normalized = normalize_name(candidate)
            if len(normalized) >= 4:
                product_index[normalized].append(locus)

    proteomics: dict[str, dict[str, int]] = defaultdict(lambda: {"peptides": 0, "transitions": 0})
    mapping_rows = []
    for accession, protein_name, peptides, transitions in cursor.execute("""
        SELECT protein_accession, protein_name, unique_peptide_count, transition_count
        FROM pxd022622_proteins
    """):
        normalized = normalize_name(protein_name)
        loci = sorted(set(product_index.get(normalized, [])))
        if len(loci) != 1:
            continue
        locus = loci[0]
        proteomics[locus]["peptides"] += int(peptides or 0)
        proteomics[locus]["transitions"] += int(transitions or 0)
        mapping_rows.append((
            accession, locus, genes.get(locus, (None, None))[0],
            "unique_exact_normalized_product_name", "medium",
            int(peptides or 0), int(transitions or 0),
            json.dumps({"protein_name": protein_name, "normalized_name": normalized}, ensure_ascii=False),
        ))
    cursor.executemany("INSERT INTO proteomics_gene_mapping VALUES (?,?,?,?,?,?,?,?)", mapping_rows)

    essential = {str(row[0]).lower() for row in cursor.execute("SELECT locus_tag FROM essential_genes")}
    pathways = {
        str(locus).lower(): int(count)
        for locus, count in cursor.execute("""
            SELECT gene_locus, COUNT(DISTINCT pathway_id) FROM gene_pathway_mappings GROUP BY gene_locus
        """)
    }
    centrality = {
        str(locus).lower(): (float(pagerank or 0), float(betweenness or 0))
        for locus, pagerank, betweenness in cursor.execute(
            "SELECT locus_tag, pagerank, betweenness FROM network_centrality"
        )
    }
    pageranks = sorted(value[0] for value in centrality.values())
    betweennesses = sorted(value[1] for value in centrality.values())

    variant_data: dict[str, dict[str, object]] = defaultdict(lambda: {"count": 0, "strains": set()})
    phenotype_strains = {str(row[0]).lower() for row in cursor.execute("SELECT DISTINCT strain FROM pxd022622_phenotypes")}
    for affected, strain in cursor.execute("SELECT affected_locus, strain FROM pxd022622_variants"):
        for locus in str(affected or "").lower().split("|"):
            if not locus:
                continue
            variant_data[locus]["count"] = int(variant_data[locus]["count"]) + 1
            variant_data[locus]["strains"].add(str(strain or "").lower())

    aggregates = cursor.execute("""
        SELECT lower(target_locus), MAX(target_name),
               COUNT(DISTINCT run_id), COUNT(DISTINCT comparison_id),
               COUNT(DISTINCT tf_locus), COUNT(*),
               SUM(CASE WHEN support_state='condition_supported' THEN 1 ELSE 0 END),
               SUM(CASE WHEN support_state='direction_conflict' THEN 1 ELSE 0 END),
               AVG(condition_score), MAX(condition_score), AVG(evidence_completeness),
               AVG(CASE WHEN direction_consistency IS NOT NULL THEN direction_consistency END)
        FROM condition_edge_scores
        GROUP BY lower(target_locus)
    """).fetchall()
    significant_counts = {
        str(locus).lower(): int(count)
        for locus, count in cursor.execute("""
            SELECT lower(e.target_locus),
                   COUNT(DISTINCT e.run_id || '|' || e.comparison_id || '|' || e.tf_locus)
            FROM condition_edge_scores e
            JOIN condition_regulon_summary s
              ON s.run_id=e.run_id AND s.comparison_id=e.comparison_id AND s.tf_locus=e.tf_locus
            WHERE s.validation_status='response_enriched_fdr10'
            GROUP BY lower(e.target_locus)
        """)
    }

    score_rows = []
    for row in aggregates:
        (locus, edge_name, module_count, condition_count, regulator_count, observations,
         supported, conflicts, mean_score, max_score, completeness, direction_agreement) = row
        significant_count = significant_counts.get(locus, 0)
        supported_fraction = supported / observations if observations else 0.0
        pr, btw = centrality.get(locus, (0.0, 0.0))
        pr_pct = percentile(pageranks, pr)
        btw_pct = percentile(betweennesses, btw)
        pathway_count = pathways.get(locus, 0)
        protein = proteomics.get(locus, {"peptides": 0, "transitions": 0})
        known_essential = locus in essential
        variants = variant_data.get(locus, {"count": 0, "strains": set()})
        phenotype_context = int(any(
            strain in phenotype_strains or strain.removeprefix("evolved_") in phenotype_strains
            for strain in variants["strains"]
        ))

        evidence_score = (
            0.20 * min(1.0, module_count / 5.0)
            + 0.15 * min(1.0, condition_count / 20.0)
            + 0.10 * min(1.0, regulator_count / 6.0)
            + 0.20 * float(mean_score or 0.0)
            + 0.15 * supported_fraction
            + 0.10 * float(completeness or 0.0)
            + 0.10 * min(1.0, significant_count / 4.0)
        )
        systems_impact = 0.45 * pr_pct + 0.35 * btw_pct + 0.20 * min(1.0, pathway_count / 8.0)
        essential_penalty = 1.0 if known_essential else 0.0
        engineering_tractability = (
            0.35 * (1.0 - essential_penalty)
            + 0.25 * (1.0 - (0.55 * pr_pct + 0.45 * btw_pct))
            + 0.15 * int(protein["peptides"] > 0)
            + 0.15 * float(completeness or 0.0)
            + 0.10 * min(1.0, condition_count / 10.0)
        )
        risk_score = min(1.0, 0.65 * essential_penalty + 0.25 * systems_impact + 0.10 * min(1.0, pathway_count / 8.0))
        priority = 0.50 * evidence_score + 0.25 * systems_impact + 0.25 * engineering_tractability
        if known_essential:
            strategy = "dynamic_tuning_only"
        elif systems_impact >= 0.75:
            strategy = "careful_titration"
        elif module_count >= 3 and supported_fraction >= 0.20:
            strategy = "multi_stress_control_node"
        elif pathway_count >= 3:
            strategy = "metabolic_intervention_candidate"
        else:
            strategy = "context_specific_candidate"
        grade = "A" if evidence_score >= 0.70 and module_count >= 3 else "B" if evidence_score >= 0.55 else "C" if evidence_score >= 0.40 else "D"
        gene_name, product = genes.get(locus, (edge_name, None))
        rationale = {
            "supported_fraction": supported_fraction,
            "proteomics_mapping_is_conservative": True,
            "absence_from_essential_set_is_not_nonessential_proof": not known_essential,
            "variant_strains": sorted(variants["strains"]),
            "strategy_requires_product_objective": True,
        }
        score_rows.append((
            RUN_ID, locus, gene_name or edge_name, product,
            module_count, condition_count, regulator_count, observations,
            supported, conflicts, significant_count, float(mean_score or 0.0),
            float(max_score or 0.0), float(completeness or 0.0), direction_agreement,
            pathway_count, int(protein["peptides"] > 0), int(protein["peptides"]),
            int(variants["count"]), phenotype_context,
            "known_essential" if known_essential else "not_in_curated_essential_set",
            pr_pct, btw_pct, evidence_score, systems_impact,
            engineering_tractability, risk_score, priority, grade, strategy,
            json.dumps(rationale, ensure_ascii=False),
        ))
    cursor.executemany("INSERT INTO intervention_target_scores VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", score_rows)

    module_rows = cursor.execute("""
        SELECT ?, lower(e.target_locus), e.run_id,
               COUNT(DISTINCT e.comparison_id), COUNT(DISTINCT e.tf_locus),
               SUM(CASE WHEN e.support_state='condition_supported' THEN 1 ELSE 0 END),
               SUM(CASE WHEN e.support_state='direction_conflict' THEN 1 ELSE 0 END),
               COUNT(DISTINCT CASE WHEN s.validation_status='response_enriched_fdr10'
                    THEN e.comparison_id || '|' || e.tf_locus END),
               AVG(e.condition_score), MAX(e.condition_score)
        FROM condition_edge_scores e
        JOIN condition_regulon_summary s
          ON s.run_id=e.run_id AND s.comparison_id=e.comparison_id AND s.tf_locus=e.tf_locus
        GROUP BY lower(e.target_locus), e.run_id
    """, (RUN_ID,)).fetchall()
    cursor.executemany("INSERT INTO intervention_target_module_evidence VALUES (?,?,?,?,?,?,?,?,?,?)", module_rows)
    cursor.executescript("""
        CREATE INDEX idx_intervention_priority_score
            ON intervention_target_scores(priority_score DESC);
        CREATE INDEX idx_intervention_strategy
            ON intervention_target_scores(strategy_class, priority_score DESC);
        CREATE INDEX idx_intervention_module_target
            ON intervention_target_module_evidence(module_run_id, target_locus);
        CREATE VIEW v_intervention_target_priorities AS
        SELECT s.*,
               (SELECT GROUP_CONCAT(module_run_id, ',')
                FROM intervention_target_module_evidence m
                WHERE m.run_id=s.run_id AND m.target_locus=s.target_locus) AS modules
        FROM intervention_target_scores s
        WHERE s.run_id='cross_module_priority_v1';
    """)
    conn.commit()
    print(
        f"SUCCESS: ranked {len(score_rows)} cross-module targets and conservatively mapped "
        f"{len(mapping_rows)} proteomics accessions"
    )
    conn.close()


if __name__ == "__main__":
    build_intervention_priorities()
