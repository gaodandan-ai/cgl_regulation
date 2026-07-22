#!/usr/bin/env python3
"""Finalize the SQLite database with normalized, traceable core entities.

This compatibility migration keeps the legacy tables used by the web API while
adding normalized gene/ncRNA/regulatory-evidence/publication layers, explicit
dataset metadata, schema versioning, and the indexes used by common queries.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

import networkx as nx


ROOT_DIR = Path(__file__).resolve().parents[2]
REF_DIR = ROOT_DIR / "data" / "reference"
DB_PATH = REF_DIR / "cgl_regulation.db"
STRAIN_ID = "ATCC13032"
SCHEMA_VERSION = 12


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def normalize_locus(value: object) -> str:
    value = str(value or "").strip().lower()
    match = re.fullmatch(r"cg0*(\d+)", value)
    if match:
        return f"cg{int(match.group(1)):04d}"
    match = re.fullmatch(r"ncgl0*(\d+)", value)
    if match:
        return f"ncgl{int(match.group(1)):04d}"
    match = re.fullmatch(r"cgl0*(\d+)", value)
    if match:
        return f"cgl{int(match.group(1)):04d}"
    return value


def split_pmids(value: object) -> list[str]:
    return sorted(set(re.findall(r"\b\d{6,9}\b", str(value or ""))))


def add_column_if_missing(cursor: sqlite3.Cursor, table: str, definition: str) -> None:
    column = definition.split()[0]
    columns = {row[1] for row in cursor.execute(f"PRAGMA table_info({table})")}
    if column not in columns:
        cursor.execute(f"ALTER TABLE {table} ADD COLUMN {definition}")


def create_metadata_layer(cursor: sqlite3.Cursor) -> None:
    cursor.executescript("""
        CREATE TABLE IF NOT EXISTS organisms (
            strain_id TEXT PRIMARY KEY,
            species TEXT NOT NULL,
            strain_name TEXT NOT NULL,
            ncbi_taxonomy_id INTEGER,
            genome_accession TEXT,
            notes TEXT
        );
        CREATE TABLE IF NOT EXISTS genome_releases (
            release_id TEXT PRIMARY KEY,
            strain_id TEXT NOT NULL REFERENCES organisms(strain_id),
            accession TEXT NOT NULL,
            annotation_version TEXT,
            coordinate_method TEXT NOT NULL,
            is_estimated INTEGER NOT NULL CHECK(is_estimated IN (0, 1)),
            notes TEXT
        );
        CREATE TABLE IF NOT EXISTS schema_migrations (
            version INTEGER PRIMARY KEY,
            description TEXT NOT NULL,
            applied_at TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS dataset_metadata (
            dataset_id TEXT PRIMARY KEY,
            file_path TEXT NOT NULL,
            sha256 TEXT NOT NULL,
            byte_size INTEGER NOT NULL,
            record_count INTEGER,
            imported_at TEXT NOT NULL,
            source_version TEXT,
            notes TEXT
        );
    """)
    cursor.execute("""
        INSERT OR REPLACE INTO organisms VALUES (?, ?, ?, ?, ?, ?)
    """, (
        STRAIN_ID, "Corynebacterium glutamicum", "DSM 20300 / ATCC 13032",
        196627, "NC_003450.3", "Primary strain used by the regulatory platform",
    ))
    cursor.execute("""
        INSERT OR REPLACE INTO organisms VALUES (?, ?, ?, ?, ?, ?)
    """, (
        "STRAIN_R", "Corynebacterium glutamicum", "strain R", 196627,
        None, "Cross-strain ChIP-chip evidence retained separately from ATCC 13032",
    ))
    has_refseq_gff = (REF_DIR / "genome" / "NC_003450.3.gff3").exists()
    cursor.execute("""
        INSERT OR REPLACE INTO genome_releases VALUES (?, ?, ?, ?, ?, ?, ?)
    """, (
        "NC_003450.3-refseq-gff3" if has_refseq_gff else "NC_003450.3-estimated-layout",
        STRAIN_ID, "NC_003450.3", "NCBI RefSeq current annotation" if has_refseq_gff else "project-derived",
        "NCBI RefSeq GFF3 import" if has_refseq_gff else "ordinal estimate supplemented with local TSS annotations",
        0 if has_refseq_gff else 1,
        "Coordinates parsed from the NCBI RefSeq GFF3" if has_refseq_gff else "Coordinates are estimates, not parsed RefSeq features",
    ))
    cursor.execute("""
        INSERT OR REPLACE INTO schema_migrations(version, description, applied_at)
        VALUES (?, ?, ?)
    """, (SCHEMA_VERSION, "Add auditable cross-module intervention target priorities", utc_now()))


def create_gene_layer(cursor: sqlite3.Cursor) -> None:
    cursor.executescript("""
        CREATE TABLE IF NOT EXISTS genes (
            gene_id TEXT PRIMARY KEY,
            canonical_cg TEXT,
            canonical_cgl TEXT,
            gene_name TEXT,
            product TEXT,
            uniprot_id TEXT,
            strain_id TEXT NOT NULL REFERENCES organisms(strain_id),
            mapping_status TEXT NOT NULL,
            coordinate_status TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS gene_aliases (
            alias TEXT NOT NULL,
            gene_id TEXT NOT NULL REFERENCES genes(gene_id) ON DELETE CASCADE,
            alias_type TEXT NOT NULL,
            source TEXT NOT NULL,
            is_ambiguous INTEGER NOT NULL DEFAULT 0 CHECK(is_ambiguous IN (0, 1)),
            PRIMARY KEY(alias, gene_id)
        );
        CREATE INDEX IF NOT EXISTS idx_gene_cg ON genes(canonical_cg);
        CREATE INDEX IF NOT EXISTS idx_gene_cgl ON genes(canonical_cgl);
        CREATE INDEX IF NOT EXISTS idx_gene_name ON genes(gene_name);
        CREATE INDEX IF NOT EXISTS idx_gene_alias ON gene_aliases(alias);
        CREATE INDEX IF NOT EXISTS idx_gm_cg ON gene_mappings(cg_locus);
        CREATE INDEX IF NOT EXISTS idx_gm_cgl ON gene_mappings(cgl_locus);
        CREATE INDEX IF NOT EXISTS idx_gm_name ON gene_mappings(gene_name);
        CREATE INDEX IF NOT EXISTS idx_gm_uniprot ON gene_mappings(uniprot_id);
    """)
    cursor.execute("DELETE FROM gene_aliases")
    cursor.execute("DELETE FROM genes")

    mapping_rows = cursor.execute("""
        SELECT cgl_locus, cg_locus, gene_name, product, uniprot_id
        FROM gene_mappings ORDER BY rowid
    """).fetchall()
    genes: dict[str, dict[str, str]] = {}
    aliases: list[tuple[str, str, str, str, int]] = []
    for cgl, cg, name, product, uniprot in mapping_rows:
        cg_n, cgl_n = normalize_locus(cg), normalize_locus(cgl)
        gene_id = cgl_n or cg_n or normalize_locus(name)
        if not gene_id:
            continue
        record = genes.setdefault(gene_id, {
            "cg": cg_n, "cgl": cgl_n, "name": str(name or "").strip(),
            "product": str(product or "").strip(), "uniprot": str(uniprot or "").strip(),
            "mapping": "mapped", "coordinate": "missing",
        })
        for key, value in (("cg", cg_n), ("cgl", cgl_n), ("name", str(name or "").strip()),
                           ("product", str(product or "").strip()), ("uniprot", str(uniprot or "").strip())):
            if value and not record[key]:
                record[key] = value
        for alias, alias_type in ((cg_n, "cg_locus"), (cgl_n, "cgl_locus"), (normalize_locus(name), "gene_name")):
            if alias:
                aliases.append((alias, gene_id, alias_type, "gene_mapping.csv", 0))

    coordinate_columns = {row[1] for row in cursor.execute("PRAGMA table_info(gene_coordinates)")}
    optional = [name if name in coordinate_columns else f"NULL AS {name}"
                for name in ("coordinate_quality", "refseq_locus", "old_locus_tag", "product")]
    coordinate_rows = cursor.execute(
        f"SELECT locus_tag, gene_name, {', '.join(optional)} FROM gene_coordinates"
    ).fetchall()
    alias_to_gene = {alias: gene_id for alias, gene_id, *_ in aliases}
    for locus, name, quality, refseq_locus, old_locus, coordinate_product in coordinate_rows:
        locus_n = normalize_locus(locus)
        gene_id = alias_to_gene.get(locus_n, locus_n)
        if gene_id not in genes:
            genes[gene_id] = {
                "cg": locus_n if re.fullmatch(r"cg\d+", locus_n) else "", "cgl": "",
                "name": str(name or locus_n), "product": str(coordinate_product or ""), "uniprot": "",
                "mapping": "coordinate_only", "coordinate": quality or "unknown",
            }
            aliases.append((locus_n, gene_id, "cg_locus", "gene_coordinates", 0))
        else:
            genes[gene_id]["coordinate"] = quality or "unknown"
            if coordinate_product and not genes[gene_id]["product"]:
                genes[gene_id]["product"] = str(coordinate_product)
        for alias, alias_type in ((normalize_locus(refseq_locus), "refseq_locus"),
                                  (normalize_locus(old_locus), "old_locus_tag")):
            if alias:
                aliases.append((alias, gene_id, alias_type, "NCBI RefSeq NC_003450.3", 0))

    gene_rows = [(
        gene_id, data["cg"] or None, data["cgl"] or None, data["name"] or None,
        data["product"] or None, data["uniprot"] or None, STRAIN_ID,
        data["mapping"], data["coordinate"],
    ) for gene_id, data in genes.items()]
    cursor.executemany("INSERT INTO genes VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", gene_rows)

    alias_counts: dict[str, set[str]] = {}
    for alias, gene_id, *_ in aliases:
        alias_counts.setdefault(alias, set()).add(gene_id)
    alias_rows = {(alias, gene_id, kind, source, int(len(alias_counts[alias]) > 1))
                  for alias, gene_id, kind, source, _ in aliases}
    cursor.executemany("INSERT OR IGNORE INTO gene_aliases VALUES (?, ?, ?, ?, ?)", sorted(alias_rows))

    # Ensure every evidence endpoint has a resolvable placeholder entity.
    endpoints = set()
    for table, a, b in (("regulations", "TF_locusTag", "TG_locusTag"),
                        ("evidence_credibility", "tf", "tg")):
        for left, right in cursor.execute(f"SELECT {a}, {b} FROM {table}"):
            endpoints.update(filter(None, (normalize_locus(left), normalize_locus(right))))
    for locus in sorted(endpoints):
        if not cursor.execute("SELECT 1 FROM gene_aliases WHERE alias=?", (locus,)).fetchone():
            evidence_strain = "STRAIN_R" if locus == "gntr1" or locus.startswith("cgr") else STRAIN_ID
            cursor.execute("INSERT OR IGNORE INTO genes VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                           (locus, locus if re.fullmatch(r"cg\d+", locus) else None,
                            locus if re.fullmatch(r"cgl\d+", locus) else None, locus, None, None,
                            evidence_strain, "cross_strain_evidence" if evidence_strain == "STRAIN_R" else "evidence_only", "missing"))
            cursor.execute("INSERT OR IGNORE INTO gene_aliases VALUES (?, ?, ?, ?, ?)",
                           (locus, locus, "evidence_locus", "regulatory evidence", 0))

    add_column_if_missing(cursor, "gene_coordinates", "coordinate_source TEXT")
    add_column_if_missing(cursor, "gene_coordinates", "coordinate_quality TEXT")
    add_column_if_missing(cursor, "gene_coordinates", "strain_id TEXT")
    cursor.execute("""
        UPDATE gene_coordinates
        SET coordinate_source=coalesce(coordinate_source, 'NC_003450.3 ordinal estimate + local TSS annotations'),
            coordinate_quality=coalesce(coordinate_quality, 'estimated'), strain_id=coalesce(strain_id, ?)
    """, (STRAIN_ID,))


def create_ncrna_layer(cursor: sqlite3.Cursor) -> None:
    cursor.executescript("""
        CREATE TABLE IF NOT EXISTS ncrnas (
            ncrna_id TEXT PRIMARY KEY,
            ncrna_name TEXT NOT NULL,
            rna_type TEXT NOT NULL,
            annotation_class TEXT NOT NULL,
            rfam_acc TEXT,
            ena_acc TEXT,
            start_pos INTEGER,
            end_pos INTEGER,
            strand TEXT,
            strain_id TEXT NOT NULL REFERENCES organisms(strain_id),
            source TEXT NOT NULL,
            description TEXT
        );
        CREATE INDEX IF NOT EXISTS idx_ncrna_type ON ncrnas(rna_type);
        CREATE INDEX IF NOT EXISTS idx_ncrna_rfam ON ncrnas(rfam_acc);
    """)
    cursor.execute("DELETE FROM ncrnas")
    cursor.execute("""
        INSERT INTO ncrnas
        SELECT lower(ncrna_id), ncrna_name, rna_type, 'manual_reference',
               CASE WHEN rfam_acc LIKE '%_like' THEN NULL ELSE rfam_acc END,
               ena_acc, NULL, NULL, NULL, ?, 'Project-curated reference list', description
        FROM rfam_ncrnas
    """, (STRAIN_ID,))
    coordinate_columns = {row[1] for row in cursor.execute("PRAGMA table_info(gene_coordinates)")}
    if "gene_biotype" in coordinate_columns:
        for row in cursor.execute("""
            SELECT locus_tag, gene_name, gene_biotype, start_pos, end_pos, strand, refseq_locus, product
            FROM gene_coordinates
            WHERE lower(gene_biotype) IN ('rrna', 'trna', 'ncrna', 'tmrna')
        """).fetchall():
            cursor.execute("""
                INSERT OR REPLACE INTO ncrnas
                VALUES (?, ?, ?, 'refseq_annotation', NULL, NULL, ?, ?, ?, ?, 'NCBI RefSeq NC_003450.3 GFF3', ?)
            """, (row[0], row[1], row[2], row[3], row[4], row[5], STRAIN_ID, row[7]))
    sources = cursor.execute("""
        SELECT DISTINCT lower(source_locus) FROM network_edges_extended
        WHERE edge_type='srna_mrna' AND trim(source_locus)<>''
    """).fetchall()
    for (source,) in sources:
        cursor.execute("""
            INSERT OR IGNORE INTO ncrnas
            VALUES (?, ?, 'sRNA', 'computational_prediction', NULL, NULL, NULL, NULL, NULL, ?, ?, ?)
        """, (source, f"sRNA {source}", STRAIN_ID, "CopraRNA/IntaRNA", "Predicted sRNA with target interactions"))
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_ncrna_interaction_pair ON ncrna_target_interactions(srna_id, target_locus)")
    cursor.executescript("""
        DROP TABLE IF EXISTS ncrna_prediction_summary;
        CREATE TABLE ncrna_prediction_summary AS
        SELECT srna_id,
               count(*) AS predicted_target_count,
               sum(confidence_tier='HIGH') AS high_confidence_count,
               sum(confidence_tier='MEDIUM') AS medium_confidence_count,
               min(copra_fdr) AS best_copra_fdr,
               min(binding_energy_kcal) AS strongest_binding_energy_kcal,
               'computational predictions; not experimental validation' AS interpretation
        FROM ncrna_target_interactions GROUP BY srna_id;
        CREATE UNIQUE INDEX idx_ncrna_prediction_summary_id
        ON ncrna_prediction_summary(srna_id);
    """)


def enrich_tf_annotations(cursor: sqlite3.Cursor) -> None:
    """Add conservative RefSeq-derived families and expose mapping conflicts."""
    add_column_if_missing(cursor, "tf_families_effectors", "annotation_note TEXT")
    cursor.executescript("""
        DROP TABLE IF EXISTS tf_annotation_conflicts;
        CREATE TABLE tf_annotation_conflicts (
            tf_locus TEXT PRIMARY KEY,
            asserted_tf_name TEXT,
            refseq_gene_name TEXT,
            refseq_product TEXT,
            conflict_type TEXT NOT NULL,
            resolution TEXT NOT NULL
        );
    """)
    family_patterns = (
        (r"\b(tetr|acrr)\b", "TetR/AcrR family"),
        (r"\barsi?r|smtb\b", "ArsR/SmtB family"),
        (r"\blys?r\b", "LysR family"),
        (r"\bgnt?r\b", "GntR family"),
        (r"\bmar?r\b", "MarR family"),
        (r"\barac\b", "AraC family"),
        (r"\biclr\b", "IclR family"),
        (r"\bdeor\b", "DeoR family"),
        (r"\blaci\b", "LacI family"),
        (r"\bwhib\b", "WhiB family"),
        (r"sigma(?:-70)? factor", "Sigma factor"),
        (r"response regulator", "Two-component response regulator"),
    )
    rows = cursor.execute("""
        SELECT tf.tf_locus, tf.tf_name, tf.annotation_status,
               gc.gene_name, gc.product
        FROM tf_families_effectors tf
        LEFT JOIN gene_coordinates gc ON gc.locus_tag=tf.tf_locus
    """).fetchall()
    for locus, asserted_name, status, refseq_name, product in rows:
        asserted = str(asserted_name or "").strip()
        refseq = str(refseq_name or "").strip()
        product_text = str(product or "").strip()
        refseq_is_named = refseq and not refseq.lower().startswith("cgl_rs")
        asserted_is_named = asserted and asserted.lower() != str(locus).lower()
        if refseq_is_named and asserted_is_named and refseq.lower() != asserted.lower():
            cursor.execute("""
                INSERT INTO tf_annotation_conflicts VALUES (?, ?, ?, ?, ?, ?)
            """, (locus, asserted, refseq, product_text or None, "gene_name_mismatch",
                  "Do not transfer RefSeq product annotation to this regulator record"))
            cursor.execute("""
                UPDATE tf_families_effectors
                SET annotation_status='mapping_conflict', annotation_note=?
                WHERE tf_locus=?
            """, (f"RefSeq gene name at locus is {refseq}; asserted regulator name is {asserted}", locus))
            continue
        if status != "unclassified" or not product_text:
            continue
        lower_product = product_text.lower()
        if not any(token in lower_product for token in ("transcription", "response regulator", "sigma factor")):
            continue
        family = "Transcriptional regulator (family unspecified)"
        for pattern, inferred_family in family_patterns:
            if re.search(pattern, lower_product):
                family = inferred_family
                break
        cursor.execute("""
            UPDATE tf_families_effectors
            SET tf_family=?, annotation_status='refseq_product_inferred',
                annotation_source='NCBI RefSeq NC_003450.3 product annotation',
                annotation_note='Family inferred only when explicitly supported by product text'
            WHERE tf_locus=?
        """, (family, locus))


def create_regulatory_layer(cursor: sqlite3.Cursor) -> None:
    cursor.executescript("""
        DROP TABLE IF EXISTS edge_evidence;
        DROP TABLE IF EXISTS regulatory_edges;
        CREATE TABLE regulatory_edges (
            edge_id TEXT PRIMARY KEY,
            regulator_locus TEXT NOT NULL,
            target_locus TEXT NOT NULL,
            regulator_gene_id TEXT NOT NULL REFERENCES genes(gene_id),
            target_gene_id TEXT NOT NULL REFERENCES genes(gene_id),
            regulator_name TEXT,
            target_name TEXT,
            regulation_role TEXT,
            is_sigma_factor INTEGER NOT NULL DEFAULT 0 CHECK(is_sigma_factor IN (0, 1)),
            confidence_score REAL CHECK(confidence_score IS NULL OR (confidence_score >= 0 AND confidence_score <= 1)),
            confidence_label TEXT CHECK(confidence_label IN ('HIGH', 'MEDIUM', 'LOW')),
            primary_source TEXT,
            strain_id TEXT NOT NULL REFERENCES organisms(strain_id),
            UNIQUE(regulator_locus, target_locus, strain_id)
        );
        CREATE TABLE edge_evidence (
            evidence_id INTEGER PRIMARY KEY AUTOINCREMENT,
            edge_id TEXT NOT NULL REFERENCES regulatory_edges(edge_id) ON DELETE CASCADE,
            evidence_type TEXT NOT NULL,
            source_name TEXT NOT NULL,
            source_accession TEXT,
            pmid TEXT,
            binding_site TEXT,
            evidence_score REAL CHECK(evidence_score IS NULL OR (evidence_score >= 0 AND evidence_score <= 1)),
            strain_id TEXT REFERENCES organisms(strain_id),
            details TEXT,
            UNIQUE(edge_id, evidence_type, source_name, pmid, binding_site)
        );
        CREATE INDEX idx_reg_edge_regulator ON regulatory_edges(regulator_locus);
        CREATE INDEX idx_reg_edge_target ON regulatory_edges(target_locus);
        CREATE INDEX idx_edge_evidence_edge ON edge_evidence(edge_id);
        CREATE INDEX idx_edge_evidence_pmid ON edge_evidence(pmid);
    """)

    edges: dict[tuple[str, str], dict] = {}
    evidence_rows: set[tuple] = set()
    publication_ids: set[str] = set()

    def canonicalize_locus(value: object) -> str:
        alias = normalize_locus(value)
        if not alias:
            return ""
        row = cursor.execute("""
            SELECT coalesce(g.canonical_cg, g.canonical_cgl, g.gene_id)
            FROM gene_aliases a JOIN genes g ON g.gene_id=a.gene_id
            WHERE a.alias=? ORDER BY a.is_ambiguous, g.gene_id LIMIT 1
        """, (alias,)).fetchone()
        return row[0] if row else alias

    def add_edge(tf: object, tg: object, tf_name: object, tg_name: object,
                 role: object, sigma: object, score: object, label: object,
                 source: str, evidence_type: object, pmids: object = "",
                 binding_site: object = "", accession: object = "", details: object = None,
                 strain_id: str = STRAIN_ID) -> None:
        tf_n, tg_n = canonicalize_locus(tf), canonicalize_locus(tg)
        if not tf_n or not tg_n:
            return
        key = (tf_n, tg_n, strain_id)
        try:
            numeric_score = float(score) if str(score or "").strip() else None
        except ValueError:
            numeric_score = None
        numeric_score = min(1.0, max(0.0, numeric_score)) if numeric_score is not None else None
        item = edges.setdefault(key, {
            "tf_name": str(tf_name or tf_n).strip(), "tg_name": str(tg_name or tg_n).strip(),
            "roles": set(), "sigma": 0, "score": None, "label": "LOW", "sources": set(),
        })
        if str(role or "").strip():
            item["roles"].add(str(role).strip())
        item["sigma"] = max(item["sigma"], int(str(sigma).lower() in {"1", "yes", "true"}))
        if numeric_score is not None and (item["score"] is None or numeric_score > item["score"]):
            item["score"] = numeric_score
        label_n = str(label or "").upper()
        rank = {"LOW": 0, "MEDIUM": 1, "HIGH": 2}
        if rank.get(label_n, -1) > rank.get(item["label"], -1):
            item["label"] = label_n
        item["sources"].add(source)
        edge_id = "reg_" + hashlib.sha1(f"{tf_n}|{tg_n}|{strain_id}".encode()).hexdigest()[:20]
        pmid_list = split_pmids(pmids) or [None]
        for pmid in pmid_list:
            if pmid:
                publication_ids.add(pmid)
            evidence_rows.add((
                edge_id, str(evidence_type or "unspecified"), source,
                str(accession or "") or None, pmid, str(binding_site or "") or None,
                numeric_score, strain_id,
                json.dumps(details, ensure_ascii=False, sort_keys=True) if details is not None else None,
            ))

    for filename in ("regulations.csv", "chipseq_regulations.csv", "regprecise_regulations.csv"):
        with (REF_DIR / filename).open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                add_edge(
                    row.get("TF_locusTag"), row.get("TG_locusTag"), row.get("TF_name"), row.get("TG_name"),
                    row.get("Role"), row.get("Is_sigma_factor"), row.get("evidence_score"),
                    row.get("confidence_label"), row.get("Source") or filename,
                    row.get("Evidence"), row.get("PMID"), row.get("Binding_site"),
                    row.get("source_accession"), {"file": filename, "strain_note": row.get("strain_note")},
                    "STRAIN_R" if "strain_r" in str(row.get("strain_group", "")).lower() else STRAIN_ID,
                )

    with (REF_DIR / "tcs_regulations.csv").open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            targets = [normalize_locus(x) for x in row.get("target_genes", "").split(";") if x.strip()]
            names = [x.strip() for x in row.get("target_names", "").split(";")]
            for index, target in enumerate(targets):
                add_edge(
                    row.get("rr_locus"), target, row.get("rr_name"), names[index] if index < len(names) else target,
                    row.get("regulation_role"), False, 1.0 if row.get("evidence") == "experimental" else 0.6,
                    "HIGH" if row.get("evidence") == "experimental" else "MEDIUM", "TCS curated",
                    row.get("evidence"), row.get("pmid"), accession=row.get("system_name"), details=row,
                )

    collectf_path = REF_DIR / "collectf_tfbs.json"
    if collectf_path.exists():
        for row in json.loads(collectf_path.read_text(encoding="utf-8")):
            add_edge(
                row.get("tf_locus"), row.get("target_locus"), row.get("tf_name"), row.get("target_name"),
                "binding", False, row.get("credibility_score", 0.95), "HIGH", row.get("source_tag", "CollecTF"),
                row.get("technique", "TFBS"), row.get("pmid"), row.get("sequence"), details={"coordinates": [row.get("start_pos"), row.get("end_pos")]},
            )

    for row in cursor.execute("SELECT * FROM evidence_credibility").fetchall():
        data = dict(row)
        types = int(data.get("n_evidence_types") or 0)
        label = "HIGH" if types >= 3 else "MEDIUM" if types == 2 else "LOW"
        evidence_strain = "STRAIN_R" if (
            str(data.get("tf", "")).lower() == "gntr1"
            or str(data.get("tf", "")).lower().startswith("cgr")
            or str(data.get("tg", "")).lower().startswith("cgr")
        ) else STRAIN_ID
        add_edge(
            data.get("tf"), data.get("tg"), data.get("tf_name"), data.get("tg_name"), None, False,
            data.get("credibility_score"), label, "Evidence credibility synthesis", "aggregated_evidence",
            details={k: data[k] for k in data if k.startswith("n_") or k == "credibility_tier"},
            strain_id=evidence_strain,
        )

    # Add endpoint placeholders before enforcing the edge-to-gene relationships logically.
    for tf_n, tg_n, _edge_strain in edges:
        for locus in (tf_n, tg_n):
            if not cursor.execute("SELECT 1 FROM gene_aliases WHERE alias=?", (locus,)).fetchone():
                cursor.execute("INSERT OR IGNORE INTO genes VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                               (locus, locus if re.fullmatch(r"cg\d+", locus) else None,
                                locus if re.fullmatch(r"cgl\d+", locus) else None, locus, None, None,
                                _edge_strain,
                                "cross_strain_evidence" if _edge_strain == "STRAIN_R" else "evidence_only",
                                "missing"))
                cursor.execute("INSERT OR IGNORE INTO gene_aliases VALUES (?, ?, ?, ?, ?)",
                               (locus, locus, "evidence_locus", "regulatory evidence", 0))

    regulatory_rows = []
    for (tf_n, tg_n, edge_strain), item in sorted(edges.items()):
        score = item["score"] if item["score"] is not None else 0.2
        label = item["label"] if item["label"] in {"HIGH", "MEDIUM", "LOW"} else ("HIGH" if score >= 0.8 else "MEDIUM" if score >= 0.4 else "LOW")
        regulator_gene_id = cursor.execute(
            "SELECT gene_id FROM gene_aliases WHERE alias=? ORDER BY is_ambiguous, gene_id LIMIT 1", (tf_n,)
        ).fetchone()[0]
        target_gene_id = cursor.execute(
            "SELECT gene_id FROM gene_aliases WHERE alias=? ORDER BY is_ambiguous, gene_id LIMIT 1", (tg_n,)
        ).fetchone()[0]
        regulatory_rows.append((
            "reg_" + hashlib.sha1(f"{tf_n}|{tg_n}|{edge_strain}".encode()).hexdigest()[:20], tf_n, tg_n,
            regulator_gene_id, target_gene_id, item["tf_name"], item["tg_name"],
            ";".join(sorted(item["roles"])), item["sigma"], score, label,
            "; ".join(sorted(item["sources"])), edge_strain,
        ))
    cursor.executemany("INSERT INTO regulatory_edges VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", regulatory_rows)
    cursor.executemany("""
        INSERT OR IGNORE INTO edge_evidence
        (edge_id, evidence_type, source_name, source_accession, pmid, binding_site,
         evidence_score, strain_id, details) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, sorted(evidence_rows, key=lambda row: tuple("" if value is None else str(value) for value in row)))

    create_publication_layer(cursor, publication_ids)

    add_column_if_missing(cursor, "network_edges_extended", "confidence_score REAL")
    add_column_if_missing(cursor, "network_edges_extended", "binding_energy_kcal REAL")
    add_column_if_missing(cursor, "network_edges_extended", "strain_id TEXT")
    cursor.execute("UPDATE network_edges_extended SET binding_energy_kcal=score WHERE edge_type='srna_mrna' AND binding_energy_kcal IS NULL")
    cursor.execute("UPDATE network_edges_extended SET strain_id=? WHERE strain_id IS NULL", (STRAIN_ID,))
    cursor.execute("DELETE FROM network_edges_extended WHERE edge_type='tf_dna'")
    cursor.executemany("""
        INSERT INTO network_edges_extended
        (edge_id, source_locus, source_name, target_locus, target_name, edge_type,
         evidence_level, score, details, confidence_score, binding_energy_kcal, strain_id)
        VALUES (?, ?, ?, ?, ?, 'tf_dna', ?, ?, ?, ?, NULL, ?)
    """, [(
        "tf_dna_" + row[0][4:], row[1], row[5], row[2], row[6],
        "strong" if row[10] in {"HIGH", "MEDIUM"} else "all", row[9],
        json.dumps({"role": row[7], "label": row[10], "sources": row[11]}, ensure_ascii=False), row[9], row[12],
    ) for row in regulatory_rows])
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_edge_src_type ON network_edges_extended(source_locus, edge_type, evidence_level)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_edge_tgt_type ON network_edges_extended(target_locus, edge_type, evidence_level)")


def recompute_network_centrality(cursor: sqlite3.Cursor) -> None:
    graph = nx.DiGraph()
    for source, target, score in cursor.execute(
        "SELECT regulator_locus, target_locus, confidence_score FROM regulatory_edges WHERE strain_id=?",
        (STRAIN_ID,)
    ).fetchall():
        graph.add_edge(source, target, weight=float(score or 0.0))
    if not graph:
        return
    pagerank = nx.pagerank(graph, weight="weight")
    closeness = nx.closeness_centrality(graph)
    sample_size = min(300, graph.number_of_nodes())
    betweenness = nx.betweenness_centrality(graph, k=sample_size, normalized=True, seed=42)
    add_column_if_missing(cursor, "network_centrality", "method TEXT")
    add_column_if_missing(cursor, "network_centrality", "edge_count INTEGER")
    add_column_if_missing(cursor, "network_centrality", "computed_at TEXT")
    cursor.execute("DELETE FROM network_centrality")
    rows = []
    for node in graph.nodes:
        in_degree = graph.in_degree(node)
        out_degree = graph.out_degree(node)
        rows.append((node, in_degree + out_degree, in_degree, out_degree,
                     betweenness.get(node, 0.0), closeness.get(node, 0.0),
                     pagerank.get(node, 0.0), "normalized regulatory_edges; approximate betweenness k=300",
                     graph.number_of_edges(), utc_now()))
    cursor.executemany("""
        INSERT INTO network_centrality
        (locus_tag, degree, in_degree, out_degree, betweenness, closeness, pagerank,
         method, edge_count, computed_at) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, rows)


def create_publication_layer(cursor: sqlite3.Cursor, publication_ids: set[str]) -> None:
    cursor.executescript("""
        CREATE TABLE IF NOT EXISTS publications (
            pmid TEXT PRIMARY KEY,
            title TEXT,
            abstract TEXT,
            journal TEXT,
            publication_year INTEGER,
            doi TEXT,
            metadata_source TEXT NOT NULL,
            retrieved_at TEXT
        );
        CREATE TABLE IF NOT EXISTS literature_documents (
            doc_id TEXT PRIMARY KEY,
            pmid TEXT REFERENCES publications(pmid),
            gene_locus TEXT,
            title TEXT NOT NULL,
            abstract TEXT,
            source TEXT NOT NULL,
            retrieved_at TEXT
        );
        CREATE INDEX IF NOT EXISTS idx_lit_doc_pmid ON literature_documents(pmid);
        CREATE INDEX IF NOT EXISTS idx_lit_doc_gene ON literature_documents(gene_locus);
    """)
    cache_path = REF_DIR / "literature" / "pubmed_records.json"
    cache = {}
    if cache_path.exists():
        raw = json.loads(cache_path.read_text(encoding="utf-8"))
        cache = raw.get("records", raw) if isinstance(raw, dict) else {}
    for pmid in sorted(publication_ids):
        item = cache.get(pmid, {}) if isinstance(cache, dict) else {}
        cursor.execute("""
            INSERT OR REPLACE INTO publications
            (pmid, title, abstract, journal, publication_year, doi, metadata_source, retrieved_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            pmid, item.get("title"), item.get("abstract"), item.get("journal"), item.get("year"),
            item.get("doi"), "NCBI PubMed" if item else "local evidence PMID", item.get("retrieved_at") or utc_now(),
        ))

    cursor.execute("DELETE FROM literature_documents")
    publication_rows = cursor.execute(
        "SELECT pmid, title, abstract, retrieved_at FROM publications WHERE title IS NOT NULL"
    ).fetchall()
    for row in publication_rows:
        genes = [r[0] for r in cursor.execute("""
            SELECT DISTINCT re.regulator_locus FROM edge_evidence ee
            JOIN regulatory_edges re ON re.edge_id=ee.edge_id WHERE ee.pmid=?
            UNION
            SELECT DISTINCT re.target_locus FROM edge_evidence ee
            JOIN regulatory_edges re ON re.edge_id=ee.edge_id WHERE ee.pmid=?
        """, (row[0], row[0])).fetchall()]
        cursor.execute("INSERT INTO literature_documents VALUES (?, ?, ?, ?, ?, ?, ?)",
                       (f"pmid:{row[0]}", row[0], ",".join(genes), row[1], row[2], "NCBI PubMed", row[3]))

    cursor.execute("DELETE FROM literature_fts")
    cursor.execute("""
        INSERT INTO literature_fts(doc_id, gene_locus, title, abstract, pmid)
        SELECT doc_id, gene_locus, title, coalesce(abstract, ''), pmid FROM literature_documents
    """)


def record_dataset_metadata(cursor: sqlite3.Cursor) -> None:
    candidates = [
        "gene_mapping.csv", "regulations.csv", "chipseq_regulations.csv",
        "regprecise_regulations.csv", "tcs_regulations.csv", "rna_regulation.csv",
        "network_centrality.json", "thermo_dgr_data.json", "literature/pubmed_records.json",
        "genome/NC_003450.3.gff3", "kegg_cache/kegg_pathway_hierarchy.txt",
        "expression_compendium/raw/Filtered differential expression results.xlsx",
        "model/ecCGL1-main/ecCGL1-main/iCW773_get_data/reaction_kcat_MW.csv",
        "genome/GCF_000011325.1/feature_table.txt.gz",
        "genome/GCF_000011325.1/genomic.gbff.gz",
        "genome/GCF_000011325.1/protein.faa.gz",
        "expression_compendium/geo/GSE169361_series_matrix.txt.gz",
        "expression_compendium/geo/GSE171301_series_matrix.txt.gz",
        "expression_compendium/geo/GSE169361_family.xml.tgz",
        "expression_compendium/geo/GSE171301_family.xml.tgz",
        "transcript_structure/Table_S1.xlsx", "transcript_structure/Table_S3.xlsx",
        "transcript_structure/Table_S4.xlsx", "transcript_structure/Table_S5.xlsx",
        "transcript_structure/Table_S9.xlsx",
        "proteomics/PXD022622/CorynemanualDB200320.txt",
        "proteomics/PXD022622/Data_Sheet_1.PDF",
        "proteomics/PXD022622/PRJNA678589_runinfo.csv",
    ]
    now = utc_now()
    for relative in candidates:
        path = REF_DIR / relative
        if not path.exists():
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        count = None
        try:
            if path.suffix == ".csv":
                delimiter = "\t" if path.name == "rna_regulation.csv" else ","
                with path.open(encoding="utf-8-sig", newline="") as handle:
                    count = sum(1 for _ in csv.DictReader(handle, delimiter=delimiter))
            elif path.name == "network_centrality.json":
                count = len(json.loads(path.read_text(encoding="utf-8")).get("nodes", {}))
            elif path.name == "thermo_dgr_data.json":
                count = len(json.loads(path.read_text(encoding="utf-8")).get("reactions", {}))
            elif path.name == "pubmed_records.json":
                count = len(json.loads(path.read_text(encoding="utf-8")).get("records", {}))
        except (OSError, ValueError, TypeError):
            count = None
        cursor.execute("""
            INSERT OR REPLACE INTO dataset_metadata
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, (relative.replace("/", ":"), relative, digest, path.stat().st_size, count, now, None, None))


def recreate_views(cursor: sqlite3.Cursor) -> None:
    cursor.executescript("""
        DROP VIEW IF EXISTS v_gene_full_profile;
        CREATE VIEW v_gene_full_profile AS
        SELECT g.cg_locus, g.cgl_locus, g.gene_name, g.product,
               c.start_pos, c.end_pos, c.strand, c.gene_length,
               c.tss_position, c.promoter_70bp, c.coordinate_source, c.coordinate_quality,
               tf.tf_family, tf.hth_domain, tf.effector_molecule,
               tf.physiological_signal, tf.regulatory_role,
               a.systemic_role AS abasy_role,
               nc.degree, nc.in_degree, nc.out_degree, nc.betweenness,
               nc.closeness, nc.pagerank
        FROM gene_mappings g
        LEFT JOIN gene_coordinates c ON c.locus_tag IN (lower(g.cg_locus), lower(g.cgl_locus))
        LEFT JOIN tf_families_effectors tf ON lower(g.cg_locus)=tf.tf_locus
        LEFT JOIN abasy_roles a ON lower(g.cg_locus)=a.locus_tag
        LEFT JOIN network_centrality nc ON lower(g.cg_locus)=nc.locus_tag;

        DROP VIEW IF EXISTS v_metabolite_tf_feedback;
        CREATE VIEW v_metabolite_tf_feedback AS
        SELECT tf.tf_locus, tf.tf_name, tf.effector_molecule,
               tf.physiological_signal, e.target_locus, e.target_name,
               e.confidence_score AS score, e.details
        FROM tf_families_effectors tf
        JOIN network_edges_extended e ON tf.tf_locus=e.source_locus
        WHERE e.edge_type='tf_dna' AND tf.effector_molecule IS NOT NULL
              AND tf.effector_molecule<>'';

        DROP VIEW IF EXISTS v_srna_competition_ranking;
        CREATE VIEW v_srna_competition_ranking AS
        SELECT source_locus AS srna_id, target_locus AS mrna_id,
               target_name AS mrna_name, binding_energy_kcal AS binding_energy, details
        FROM network_edges_extended WHERE edge_type='srna_mrna'
        ORDER BY binding_energy_kcal ASC;
    """)


def finalize_database() -> None:
    if not DB_PATH.exists():
        raise FileNotFoundError(f"Database not found: {DB_PATH}")
    connection = sqlite3.connect(DB_PATH)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys=ON")
    try:
        with connection:
            cursor = connection.cursor()
            create_metadata_layer(cursor)
            # Normalized regulatory tables depend on genes and are rebuilt
            # below; remove them first so this finalizer is safely repeatable.
            cursor.executescript("DROP TABLE IF EXISTS edge_evidence; DROP TABLE IF EXISTS regulatory_edges;")
            create_gene_layer(cursor)
            enrich_tf_annotations(cursor)
            create_ncrna_layer(cursor)
            create_regulatory_layer(cursor)
            recompute_network_centrality(cursor)
            record_dataset_metadata(cursor)
            recreate_views(cursor)
            cursor.execute(f"PRAGMA user_version={SCHEMA_VERSION}")
            cursor.execute("PRAGMA application_id=1128741970")
        connection.execute("ANALYZE")
        issues = connection.execute("PRAGMA foreign_key_check").fetchall()
        if issues:
            raise RuntimeError(f"Foreign-key validation failed: {issues[:5]}")
    finally:
        connection.close()
    print(f"SUCCESS: finalized schema v{SCHEMA_VERSION} at {DB_PATH}")


if __name__ == "__main__":
    finalize_database()
