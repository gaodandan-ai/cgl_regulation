#!/usr/bin/env python3
"""Import authoritative NC_003450.3 gene coordinates from the NCBI GFF3."""

from __future__ import annotations

import os
import re
import sqlite3
import urllib.parse
from collections import Counter
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[2]
REF_DIR = ROOT_DIR / "data" / "reference"
DB_PATH = REF_DIR / "cgl_regulation.db"
GFF_PATH = REF_DIR / "genome" / "NC_003450.3.gff3"


def parse_attributes(value: str) -> dict[str, str]:
    result = {}
    for item in value.strip().split(";"):
        if "=" in item:
            key, raw = item.split("=", 1)
            result[key] = urllib.parse.unquote(raw)
    return result


def import_refseq_genome() -> None:
    if not GFF_PATH.exists():
        raise FileNotFoundError(f"NCBI GFF3 not found: {GFF_PATH}")

    features: dict[str, dict] = {}
    children: dict[str, dict] = {}
    with GFF_PATH.open(encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, raw_attrs = fields
            attrs = parse_attributes(raw_attrs)
            if feature_type == "gene":
                feature_id = attrs.get("ID", f"gene:{start}:{end}")
                features[feature_id] = {
                    "refseq_locus": attrs.get("locus_tag", ""),
                    "old_locus": attrs.get("old_locus_tag", "").split(",")[0],
                    "gene_name": attrs.get("gene") or attrs.get("Name") or attrs.get("locus_tag", ""),
                    "gene_biotype": attrs.get("gene_biotype", feature_type),
                    "ncbi_gene_id": next((x.split(":", 1)[1] for x in attrs.get("Dbxref", "").split(",") if x.startswith("GeneID:")), None),
                    "start": int(start), "end": int(end), "strand": strand,
                    "source": source, "product": None, "protein_id": None,
                }
            elif feature_type in {"CDS", "rRNA", "tRNA", "ncRNA", "tmRNA"}:
                parent = attrs.get("Parent")
                if parent:
                    children[parent] = {
                        "product": attrs.get("product"),
                        "protein_id": attrs.get("protein_id"),
                        "gene_biotype": feature_type,
                    }

    for feature_id, child in children.items():
        if feature_id in features:
            features[feature_id].update({k: v for k, v in child.items() if v})

    connection = sqlite3.connect(DB_PATH)
    cursor = connection.cursor()
    # Cgl#### (project legacy mapping) and NCgl#### (RefSeq old locus tag)
    # are distinct numbering systems.  They must not be joined by simply
    # stripping the leading "N".  Only an unambiguous gene-name match is safe
    # enough to attach a RefSeq feature directly to a canonical cg locus.
    names: dict[str, set[str]] = {}
    for gene_name, cg in cursor.execute("SELECT gene_name, cg_locus FROM gene_mappings"):
        name = str(gene_name or "").strip().lower()
        cg_locus = str(cg or "").strip().lower()
        if name and cg_locus and name not in {"-", "none", "unknown"}:
            names.setdefault(name, set()).add(cg_locus)
    refseq_name_counts = Counter(
        str(item["gene_name"] or "").strip().lower()
        for item in features.values()
        if item.get("gene_name") and not str(item["gene_name"]).lower().startswith("cgl_rs")
    )
    mapping_by_name = {
        name: next(iter(loci)) for name, loci in names.items()
        if len(loci) == 1 and refseq_name_counts.get(name) == 1
    }
    # Many RefSeq gene features have no `gene=` qualifier, while the child CDS
    # product ends with an explicit symbol (for example "transcriptional
    # regulator SugR").  Accept that exact terminal symbol only when it is
    # unique in both datasets.
    product_symbol_counts = Counter()
    feature_product_symbol: dict[str, str] = {}
    for feature_id, item in features.items():
        match = re.search(r"\b([A-Za-z][A-Za-z0-9]{2,})$", str(item.get("product") or "").strip())
        symbol = match.group(1).lower() if match else ""
        if symbol in names:
            feature_product_symbol[feature_id] = symbol
            product_symbol_counts[symbol] += 1
    mapping_by_product_symbol = {
        symbol: next(iter(names[symbol])) for symbol, count in product_symbol_counts.items()
        if count == 1 and len(names[symbol]) == 1
    }

    cursor.execute("DROP TABLE IF EXISTS gene_coordinates")
    cursor.execute("""
        CREATE TABLE gene_coordinates (
            locus_tag TEXT PRIMARY KEY,
            gene_name TEXT,
            start_pos INTEGER NOT NULL CHECK(start_pos > 0),
            end_pos INTEGER NOT NULL CHECK(end_pos >= start_pos),
            strand TEXT NOT NULL CHECK(strand IN ('+', '-')),
            gene_length INTEGER NOT NULL,
            tss_position REAL,
            promoter_70bp TEXT,
            refseq_locus TEXT NOT NULL,
            old_locus_tag TEXT,
            gene_biotype TEXT,
            ncbi_gene_id TEXT,
            protein_id TEXT,
            product TEXT,
            coordinate_source TEXT NOT NULL,
            coordinate_quality TEXT NOT NULL CHECK(coordinate_quality='reference'),
            strain_id TEXT NOT NULL
        )
    """)

    tss_path = REF_DIR / "tss_promoter_annotations.json"
    import json
    tss_data = json.loads(tss_path.read_text(encoding="utf-8")) if tss_path.exists() else {}
    rows = []
    for feature_id, item in features.items():
        old = item["old_locus"]
        old_alias = old.lower()
        refseq_locus = item["refseq_locus"].lower()
        gene_name = str(item["gene_name"] or "").strip().lower()
        product_symbol = feature_product_symbol.get(feature_id, "")
        locus = (mapping_by_name.get(gene_name)
                 or mapping_by_product_symbol.get(product_symbol)
                 or refseq_locus or old_alias)
        tss = tss_data.get(locus, {})
        rows.append((
            locus, item["gene_name"], item["start"], item["end"], item["strand"],
            item["end"] - item["start"] + 1, tss.get("tss_position"),
            tss.get("promoter_70bp_upstream", ""), refseq_locus,
            old_alias or None, item["gene_biotype"], item["ncbi_gene_id"],
            item["protein_id"], item["product"], "NCBI RefSeq NC_003450.3 GFF3",
            "reference", "ATCC13032",
        ))
    cursor.executemany("INSERT OR REPLACE INTO gene_coordinates VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", rows)
    cursor.execute("CREATE INDEX idx_coord_start ON gene_coordinates(start_pos)")
    cursor.execute("CREATE INDEX idx_coord_refseq ON gene_coordinates(refseq_locus)")
    cursor.execute("CREATE INDEX idx_coord_old_locus ON gene_coordinates(old_locus_tag)")
    cursor.execute("CREATE INDEX idx_coord_gene_name ON gene_coordinates(gene_name)")
    connection.commit()
    connection.close()
    print(f"SUCCESS: imported {len(rows)} authoritative RefSeq gene coordinates from {GFF_PATH}")


if __name__ == "__main__":
    import_refseq_genome()
