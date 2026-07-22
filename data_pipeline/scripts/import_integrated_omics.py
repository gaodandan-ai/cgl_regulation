#!/usr/bin/env python3
"""Import identifier, GEO expression, transcript-structure and PXD022622 layers."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


ROOT_DIR = Path(__file__).resolve().parents[2]
REF_DIR = ROOT_DIR / "data" / "reference"
DB_PATH = REF_DIR / "cgl_regulation.db"
GENOME_DIR = REF_DIR / "genome" / "GCF_000011325.1"
GEO_DIR = REF_DIR / "expression_compendium" / "geo"
TRANSCRIPT_DIR = REF_DIR / "transcript_structure"
PXD_DIR = REF_DIR / "proteomics" / "PXD022622"


def clean(value: object) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    return text or None


def normalize_identifier(value: object) -> str:
    text = str(value or "").strip().lower()
    for pattern, prefix in ((r"cg0*(\d+)", "cg"), (r"cgl0*(\d+)", "cgl"),
                            (r"ncgl0*(\d+)", "ncgl")):
        match = re.fullmatch(pattern, text)
        if match:
            return f"{prefix}{int(match.group(1)):04d}"
    return text


def import_identifier_layer(cursor: sqlite3.Cursor) -> None:
    cursor.executescript("""
        DROP TABLE IF EXISTS gene_identifier_crosswalk;
        DROP TABLE IF EXISTS genome_assembly_catalog;
        CREATE TABLE genome_assembly_catalog (
            assembly_id TEXT PRIMARY KEY,
            accession TEXT,
            assembly_name TEXT,
            sequence_accession TEXT,
            coordinate_system TEXT NOT NULL,
            genome_length INTEGER,
            status TEXT NOT NULL,
            source_url TEXT NOT NULL,
            notes TEXT
        );
        CREATE TABLE gene_identifier_crosswalk (
            feature_id TEXT PRIMARY KEY,
            assembly_id TEXT NOT NULL REFERENCES genome_assembly_catalog(assembly_id),
            refseq_locus TEXT,
            ncgl_locus TEXT,
            cg_locus TEXT,
            cgl_locus TEXT,
            gene_symbol TEXT,
            ncbi_gene_id TEXT,
            protein_accession TEXT,
            nonredundant_protein_accession TEXT,
            uniprot_id TEXT,
            mapping_method TEXT NOT NULL,
            mapping_confidence TEXT NOT NULL CHECK(mapping_confidence IN ('HIGH','MEDIUM','UNMAPPED')),
            source TEXT NOT NULL
        );
        CREATE INDEX idx_crosswalk_refseq ON gene_identifier_crosswalk(refseq_locus);
        CREATE INDEX idx_crosswalk_ncgl ON gene_identifier_crosswalk(ncgl_locus);
        CREATE INDEX idx_crosswalk_cg ON gene_identifier_crosswalk(cg_locus);
        CREATE INDEX idx_crosswalk_symbol ON gene_identifier_crosswalk(gene_symbol);
        CREATE INDEX idx_crosswalk_protein ON gene_identifier_crosswalk(protein_accession);
    """)
    cursor.executemany("INSERT INTO genome_assembly_catalog VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", [
        ("GCF_000011325.1", "GCF_000011325.1", "ASM1132v1", "NC_003450.3",
         "NCBI RefSeq", 3282708, "reference", "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000011325.1/",
         "Feature table, GenBank flatfile and protein FASTA downloaded from NCBI"),
        ("LEGACY_CG_TRANSCRIPTOME", "BX927147", "legacy cg annotation", "BX927147",
         "Study-reported cg coordinate system", None, "legacy",
         "https://pmc.ncbi.nlm.nih.gov/articles/PMC3890552/",
         "Coordinates used by the 2013 transcript-structure supplement; never mixed with NC_003450.3"),
        ("ATCC_HYBRID_2019", None, "ATCC 13032 hybrid assembly", None,
         "ATCC Genome Portal", 3310828, "alternate_complete",
         "https://genomes.atcc.org/genomes/a24105b785224292",
         "Illumina plus Oxford Nanopore circular assembly; separate coordinate release"),
    ])

    mapping_rows = cursor.execute("""
        SELECT lower(cg_locus), lower(cgl_locus), lower(gene_name), uniprot_id
        FROM gene_mappings
    """).fetchall()
    by_symbol: dict[str, set[str]] = defaultdict(set)
    mapping_details: dict[str, tuple[str | None, str | None]] = {}
    for cg, cgl, symbol, uniprot in mapping_rows:
        if cg:
            mapping_details[cg] = (cgl, uniprot)
        if symbol and cg:
            by_symbol[symbol].add(cg)
    unique_symbol = {symbol: next(iter(ids)) for symbol, ids in by_symbol.items() if len(ids) == 1}

    feature_path = GENOME_DIR / "feature_table.txt.gz"
    genes: dict[str, dict] = {}
    with gzip.open(feature_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        reader.fieldnames = [name.lstrip("# ") for name in reader.fieldnames or []]
        for row in reader:
            locus = str(row.get("locus_tag") or "").strip().lower()
            if not locus:
                continue
            item = genes.setdefault(locus, {
                "symbol": None, "gene_id": None, "ncgl": None, "protein": None,
                "nr_protein": None, "product": None,
            })
            symbol = str(row.get("symbol") or "").strip().lower()
            if symbol:
                item["symbol"] = symbol
            item["gene_id"] = item["gene_id"] or clean(row.get("GeneID"))
            attrs = str(row.get("attributes") or "")
            old_match = re.search(r"(?:^|;)old_locus_tag=([^;]+)", attrs)
            if old_match:
                item["ncgl"] = normalize_identifier(old_match.group(1).split(",")[0])
            if row.get("feature") == "CDS":
                item["protein"] = clean(row.get("product_accession"))
                item["nr_protein"] = clean(row.get("non-redundant_refseq"))
                item["product"] = clean(row.get("name"))

    symbol_counts = Counter(item["symbol"] for item in genes.values() if item["symbol"])
    rows = []
    for refseq_locus, item in sorted(genes.items()):
        symbol = item["symbol"]
        cg = unique_symbol.get(symbol) if symbol and symbol_counts[symbol] == 1 else None
        method = "unique_gene_symbol" if cg else "no_safe_legacy_match"
        confidence = "HIGH" if cg else "UNMAPPED"
        cgl, uniprot = mapping_details.get(cg, (None, None)) if cg else (None, None)
        rows.append((
            f"GCF_000011325.1:{refseq_locus}", "GCF_000011325.1", refseq_locus,
            item["ncgl"], cg, cgl, symbol, item["gene_id"], item["protein"],
            item["nr_protein"], uniprot, method, confidence,
            "NCBI GCF_000011325.1 feature_table + project gene_mapping.csv",
        ))
    cursor.executemany("INSERT INTO gene_identifier_crosswalk VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)", rows)


def parse_series_matrix(path: Path) -> tuple[dict[str, list[list[str]]], list[str], list[tuple[str, list[str]]]]:
    metadata: dict[str, list[list[str]]] = defaultdict(list)
    samples: list[str] = []
    matrix_rows: list[tuple[str, list[str]]] = []
    in_matrix = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_matrix = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            values = next(csv.reader([line], delimiter="\t"))
            if not values:
                continue
            if not in_matrix and values[0].startswith("!Sample_"):
                metadata[values[0][1:]].append(values[1:])
                if values[0] == "!Sample_geo_accession":
                    samples = values[1:]
            elif in_matrix:
                if values[0].strip('"') == "ID_REF":
                    if not samples:
                        samples = values[1:]
                    continue
                matrix_rows.append((values[0], values[1:]))
    return metadata, samples, matrix_rows


def import_expression_layer(cursor: sqlite3.Cursor) -> None:
    cursor.executescript("""
        DROP TABLE IF EXISTS expression_values;
        DROP TABLE IF EXISTS expression_sample_attributes;
        DROP TABLE IF EXISTS expression_samples;
        DROP TABLE IF EXISTS expression_probes;
        DROP TABLE IF EXISTS expression_datasets;
        CREATE TABLE expression_datasets (
            dataset_id TEXT PRIMARY KEY,
            parent_dataset_id TEXT,
            title TEXT NOT NULL,
            platform TEXT,
            sample_count INTEGER NOT NULL,
            value_semantics TEXT NOT NULL,
            raw_data_url TEXT,
            processed_source TEXT NOT NULL
        );
        CREATE TABLE expression_samples (
            sample_id TEXT PRIMARY KEY,
            dataset_id TEXT NOT NULL REFERENCES expression_datasets(dataset_id),
            title TEXT,
            source_ch1 TEXT,
            source_ch2 TEXT,
            genotype_ch1 TEXT,
            genotype_ch2 TEXT,
            condition_ch1 TEXT,
            condition_ch2 TEXT,
            replicate INTEGER,
            channel_count INTEGER,
            platform TEXT,
            organism TEXT
        );
        CREATE TABLE expression_sample_attributes (
            sample_id TEXT NOT NULL REFERENCES expression_samples(sample_id),
            attribute_name TEXT NOT NULL,
            attribute_ordinal INTEGER NOT NULL,
            attribute_value TEXT,
            PRIMARY KEY(sample_id, attribute_name, attribute_ordinal)
        );
        CREATE TABLE expression_probes (
            probe_id TEXT PRIMARY KEY,
            resolved_locus TEXT,
            resolution_method TEXT NOT NULL,
            source_platform TEXT NOT NULL
        );
        CREATE TABLE expression_values (
            sample_id TEXT NOT NULL REFERENCES expression_samples(sample_id),
            probe_id TEXT NOT NULL REFERENCES expression_probes(probe_id),
            processed_value REAL NOT NULL,
            dataset_id TEXT NOT NULL,
            PRIMARY KEY(sample_id, probe_id)
        ) WITHOUT ROWID;
        CREATE INDEX idx_expression_values_probe ON expression_values(probe_id);
        CREATE INDEX idx_expression_samples_dataset ON expression_samples(dataset_id);
    """)
    aliases: dict[str, set[str]] = defaultdict(set)
    for cg, cgl, name in cursor.execute("SELECT lower(cg_locus), lower(cgl_locus), lower(gene_name) FROM gene_mappings"):
        for alias in (cg, cgl, name):
            if alias and cg:
                aliases[alias].add(cg)
    unique_alias = {alias: next(iter(ids)) for alias, ids in aliases.items() if len(ids) == 1}

    dataset_paths = [
        ("GSE169361", GEO_DIR / "GSE169361_series_matrix.txt.gz", 927, "with raw GPR files"),
        ("GSE171301", GEO_DIR / "GSE171301_series_matrix.txt.gz", 287, "processed values only"),
    ]
    for dataset_id, path, expected_count, note in dataset_paths:
        metadata, samples, matrix_rows = parse_series_matrix(path)
        if len(samples) != expected_count:
            raise ValueError(f"{dataset_id}: expected {expected_count} samples, found {len(samples)}")
        cursor.execute("INSERT INTO expression_datasets VALUES (?,?,?,?,?,?,?,?)", (
            dataset_id, "GSE171302", "A compendium of expression profiles for C. glutamicum ATCC13032",
            "GPL29897", len(samples),
            "GEO-supplied processed two-channel comparison value; not absolute abundance",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171302/suppl/GSE171302_RAW.tar",
            str(path.relative_to(ROOT_DIR)).replace("\\", "/") + f" ({note})",
        ))

        def first(key: str, index: int) -> str | None:
            groups = metadata.get(key, [])
            return clean(groups[0][index]) if groups and index < len(groups[0]) else None

        characteristics: dict[tuple[str, int], dict[str, str]] = defaultdict(dict)
        for key in ("Sample_characteristics_ch1", "Sample_characteristics_ch2"):
            channel = key[-3:]
            for ordinal, group in enumerate(metadata.get(key, []), 1):
                for index, value in enumerate(group):
                    text = clean(value)
                    if text and ":" in text:
                        name, parsed = text.split(":", 1)
                        characteristics[(channel, index)][name.strip().lower()] = parsed.strip()

        sample_rows = []
        attribute_rows = []
        for index, sample_id in enumerate(samples):
            title = first("Sample_title", index)
            rep_match = re.search(r"(?:rep(?:licate)?[-_ ]?)(\d+)\s*$", title or "", re.I)
            sample_rows.append((
                sample_id, dataset_id, title, first("Sample_source_name_ch1", index),
                first("Sample_source_name_ch2", index),
                characteristics.get(("ch1", index), {}).get("genotype"),
                characteristics.get(("ch2", index), {}).get("genotype"),
                characteristics.get(("ch1", index), {}).get("condition"),
                characteristics.get(("ch2", index), {}).get("condition"),
                int(rep_match.group(1)) if rep_match else None,
                int(first("Sample_channel_count", index) or 0) or None,
                first("Sample_platform_id", index) or "GPL29897",
                first("Sample_organism_ch1", index),
            ))
            for attribute_name, groups in metadata.items():
                for ordinal, group in enumerate(groups, 1):
                    attribute_rows.append((
                        sample_id, attribute_name, ordinal,
                        clean(group[index]) if index < len(group) else None,
                    ))
        cursor.executemany("INSERT INTO expression_samples VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", sample_rows)
        cursor.executemany("INSERT INTO expression_sample_attributes VALUES (?,?,?,?)", attribute_rows)

        probe_rows = []
        for probe_id, _ in matrix_rows:
            probe_norm = normalize_identifier(probe_id)
            resolved = unique_alias.get(probe_norm)
            if resolved:
                method = "unique_project_alias"
            elif re.fullmatch(r"cg\d{4}", probe_norm):
                resolved, method = probe_norm, "cg_identifier"
            else:
                method = "unresolved_platform_label"
            probe_rows.append((probe_id, resolved, method, "GPL29897"))
        cursor.executemany("INSERT OR IGNORE INTO expression_probes VALUES (?,?,?,?)", probe_rows)

        batch: list[tuple[str, str, float, str]] = []
        for probe_id, values in matrix_rows:
            for index, raw in enumerate(values[:len(samples)]):
                raw = raw.strip()
                if not raw or raw.upper() in {"NA", "NAN", "NULL"}:
                    continue
                try:
                    batch.append((samples[index], probe_id, float(raw), dataset_id))
                except ValueError:
                    continue
                if len(batch) >= 100_000:
                    cursor.executemany("INSERT OR REPLACE INTO expression_values VALUES (?,?,?,?)", batch)
                    batch.clear()
        if batch:
            cursor.executemany("INSERT OR REPLACE INTO expression_values VALUES (?,?,?,?)", batch)


def read_supplement(number: int, header: int = 1) -> pd.DataFrame:
    return pd.read_excel(TRANSCRIPT_DIR / f"Table_S{number}.xlsx", header=header, engine="openpyxl")


def import_transcript_structure(cursor: sqlite3.Cursor) -> None:
    cursor.executescript("""
        DROP TABLE IF EXISTS transcription_unit_members;
        DROP TABLE IF EXISTS transcription_units;
        DROP TABLE IF EXISTS transcription_terminators;
        DROP TABLE IF EXISTS ribosome_binding_sites;
        DROP TABLE IF EXISTS five_prime_utrs;
        DROP TABLE IF EXISTS transcription_start_sites;
        DROP TABLE IF EXISTS novel_transcripts;
        CREATE TABLE transcription_start_sites (
            tss_id TEXT PRIMARY KEY, position INTEGER NOT NULL, strand TEXT NOT NULL,
            gene_locus TEXT, minus10_score REAL, minus35_score REAL,
            upstream_sequence TEXT, minus10_distance INTEGER, promoter_spacer_length INTEGER,
            assembly_id TEXT NOT NULL, source TEXT NOT NULL
        );
        CREATE TABLE five_prime_utrs (
            utr_id TEXT PRIMARY KEY, gene_locus TEXT NOT NULL, strand TEXT NOT NULL,
            tss_position INTEGER NOT NULL, gene_start INTEGER, start_codon TEXT,
            utr_length INTEGER, utr_sequence TEXT, assembly_id TEXT NOT NULL, source TEXT NOT NULL
        );
        CREATE TABLE ribosome_binding_sites (
            rbs_id TEXT PRIMARY KEY, gene_locus TEXT, strand TEXT, tss_position INTEGER,
            utr_length INTEGER, motif_score REAL, search_sequence TEXT, spacing INTEGER,
            assembly_id TEXT NOT NULL, source TEXT NOT NULL
        );
        CREATE TABLE transcription_units (
            transcription_unit_id TEXT PRIMARY KEY, genes_raw TEXT NOT NULL,
            classification TEXT NOT NULL, strand TEXT, gene_count INTEGER,
            tss_evidence TEXT, assembly_id TEXT NOT NULL, source TEXT NOT NULL
        );
        CREATE TABLE transcription_unit_members (
            transcription_unit_id TEXT NOT NULL REFERENCES transcription_units(transcription_unit_id),
            member_order INTEGER NOT NULL, gene_locus TEXT NOT NULL,
            PRIMARY KEY(transcription_unit_id, member_order)
        );
        CREATE TABLE transcription_terminators (
            terminator_id TEXT PRIMARY KEY, gene_locus TEXT, start_pos INTEGER, end_pos INTEGER,
            strand TEXT, hairpin_score REAL, tail_score REAL, confidence_value REAL,
            transcript_end_evidence TEXT, assembly_id TEXT NOT NULL, source TEXT NOT NULL
        );
        CREATE TABLE novel_transcripts (
            transcript_id TEXT PRIMARY KEY, transcript_class TEXT NOT NULL,
            start_pos INTEGER, end_pos INTEGER, strand TEXT, length INTEGER,
            related_locus TEXT, tss_evidence TEXT, minus10_region TEXT,
            sequence TEXT, assembly_id TEXT NOT NULL, source TEXT NOT NULL
        );
        CREATE INDEX idx_tss_gene ON transcription_start_sites(gene_locus);
        CREATE INDEX idx_tss_position ON transcription_start_sites(position);
        CREATE INDEX idx_tu_member_gene ON transcription_unit_members(gene_locus);
        CREATE INDEX idx_terminator_gene ON transcription_terminators(gene_locus);
    """)
    assembly = "LEGACY_CG_TRANSCRIPTOME"
    source = "PMC3890552 supplementary tables"

    tss = read_supplement(1)
    tss_rows = []
    for i, row in tss.iterrows():
        position = int(row.iloc[0])
        strand, gene = clean(row.iloc[1]), normalize_identifier(row.iloc[2])
        tss_rows.append((
            f"tss_{position}_{strand}_{i}", position, strand, gene,
            float(row.iloc[3]) if pd.notna(row.iloc[3]) else None,
            float(row.iloc[4]) if pd.notna(row.iloc[4]) else None,
            clean(row.iloc[5]), int(row.iloc[6]) if pd.notna(row.iloc[6]) else None,
            int(row.iloc[7]) if pd.notna(row.iloc[7]) else None, assembly, source,
        ))
    cursor.executemany("INSERT INTO transcription_start_sites VALUES (?,?,?,?,?,?,?,?,?,?,?)", tss_rows)

    utrs = read_supplement(3)
    utr_rows = []
    for i, row in utrs.iterrows():
        gene, strand, position = normalize_identifier(row.iloc[0]), clean(row.iloc[1]), int(row.iloc[2])
        utr_rows.append((
            f"utr_{gene}_{position}_{i}", gene, strand, position,
            int(row.iloc[3]) if pd.notna(row.iloc[3]) else None, clean(row.iloc[4]),
            int(row.iloc[5]) if pd.notna(row.iloc[5]) else None, clean(row.iloc[6]), assembly, source,
        ))
    cursor.executemany("INSERT INTO five_prime_utrs VALUES (?,?,?,?,?,?,?,?,?,?)", utr_rows)

    rbs = read_supplement(4)
    rbs_rows = []
    for i, row in rbs.iterrows():
        gene = normalize_identifier(row.iloc[0])
        position = int(row.iloc[2])
        rbs_rows.append((
            f"rbs_{gene}_{position}_{i}", gene, clean(row.iloc[1]), position,
            int(row.iloc[3]) if pd.notna(row.iloc[3]) else None,
            float(row.iloc[4]) if pd.notna(row.iloc[4]) else None, clean(row.iloc[5]),
            int(row.iloc[6]) if pd.notna(row.iloc[6]) else None, assembly, source,
        ))
    cursor.executemany("INSERT INTO ribosome_binding_sites VALUES (?,?,?,?,?,?,?,?,?,?)", rbs_rows)

    units = read_supplement(5)
    unit_rows, member_rows = [], []
    for i, row in units.iterrows():
        genes_raw = str(row.iloc[0]).strip()
        unit_id = "tu_" + hashlib.sha1(f"{i}|{genes_raw}|{row.iloc[1]}".encode()).hexdigest()[:16]
        unit_rows.append((
            unit_id, genes_raw, str(row.iloc[1]).strip(), clean(row.iloc[2]),
            int(row.iloc[3]) if pd.notna(row.iloc[3]) else None, clean(row.iloc[4]), assembly, source,
        ))
        members = re.findall(r"cg\d{4}|(?:4\.5S|6C|M1) RNA|tmRNA|tRNA\w*", genes_raw, re.I)
        for order, member in enumerate(members, 1):
            member_rows.append((unit_id, order, normalize_identifier(member)))
    cursor.executemany("INSERT INTO transcription_units VALUES (?,?,?,?,?,?,?,?)", unit_rows)
    cursor.executemany("INSERT INTO transcription_unit_members VALUES (?,?,?)", member_rows)

    terms = read_supplement(9)
    term_rows = []
    for i, row in terms.iterrows():
        gene = normalize_identifier(row.iloc[0])
        term_rows.append((
            f"term_{gene}_{i}", gene, int(row.iloc[1]), int(row.iloc[2]), clean(row.iloc[3]),
            float(row.iloc[4]) if pd.notna(row.iloc[4]) else None,
            float(row.iloc[5]) if pd.notna(row.iloc[5]) else None,
            float(row.iloc[6]) if pd.notna(row.iloc[6]) else None, clean(row.iloc[7]), assembly, source,
        ))
    cursor.executemany("INSERT INTO transcription_terminators VALUES (?,?,?,?,?,?,?,?,?,?,?)", term_rows)

    novel_rows = []
    intergenic = read_supplement(6, header=2)
    for i, row in intergenic.iterrows():
        if pd.isna(row.iloc[0]):
            continue
        start, end = int(row.iloc[0]), int(row.iloc[1])
        novel_rows.append((f"novel_intergenic_{i}", "intergenic", start, end, clean(row.iloc[3]),
                           int(row.iloc[2]), clean(row.iloc[5]), clean(row.iloc[6]), clean(row.iloc[7]),
                           clean(row.iloc[4]), assembly, source))
    antisense = read_supplement(7)
    for i, row in antisense.iterrows():
        start, end = int(row.iloc[0]), int(row.iloc[1])
        novel_rows.append((f"novel_antisense_{i}", "antisense", start, end, clean(row.iloc[3]),
                           int(row.iloc[2]), normalize_identifier(row.iloc[5]), clean(row.iloc[6]),
                           clean(row.iloc[7]), clean(row.iloc[4]), assembly, source))
    intragenic = read_supplement(8)
    for i, row in intragenic.iterrows():
        start = int(row.iloc[0])
        novel_rows.append((f"novel_intragenic_{i}", "intragenic", start, None, clean(row.iloc[1]),
                           None, normalize_identifier(row.iloc[2]), clean(row.iloc[3]), clean(row.iloc[4]),
                           None, assembly, source))
    cursor.executemany("INSERT INTO novel_transcripts VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", novel_rows)


PHENOTYPE_ROWS = [
    ("ATCC13032","R1",.45,.44,1.47,4.86,2.20,7.16,.81), ("ATCC13032","R2",.45,.45,1.41,4.61,1.92,6.41,.84),
    ("ATCC13032","R3",.47,.46,1.53,4.98,2.23,7.24,.82), ("delta_pck","R1",.36,.37,1.26,4.03,1.63,5.32,.80),
    ("delta_pck","R2",.38,.39,1.27,4.16,1.75,5.69,.81), ("delta_ppc","R1",.45,.46,1.57,5.08,2.43,7.60,.81),
    ("delta_ppc","R2",.43,.43,1.48,5.16,2.36,7.68,.77), ("delta_pyc","R1",.38,.40,1.44,4.60,2.24,7.29,.80),
    ("delta_pyc","R2",.37,.40,1.55,4.70,2.32,7.56,.79), ("delta_malE","R1",.48,.46,1.50,5.04,2.02,7.36,.82),
    ("delta_malE","R2",.49,.46,1.37,4.72,1.82,6.64,.85), ("delta_pck_delta_malE","R1",.45,.45,1.50,4.94,2.30,7.45,.82),
    ("delta_pck_delta_malE","R2",.47,.47,1.48,4.77,2.26,7.37,.87), ("delta_pck_delta_malE","R3",.43,.45,1.51,4.85,2.15,6.98,.81),
    ("delta_ppc_delta_malE","R1",.46,.43,1.46,5.03,2.14,6.96,.77), ("delta_ppc_delta_malE","R2",.42,.42,1.48,4.77,1.89,6.16,.77),
    ("delta_ppc_delta_malE","R3",.43,.43,1.40,4.73,2.23,7.25,.83), ("delta_pyc_delta_odx","R1",.39,.40,1.40,4.54,2.00,6.51,.78),
    ("delta_pyc_delta_odx","R2",.38,.38,1.31,4.15,1.87,6.07,.82), ("delta_pyc_delta_odx","R3",.43,.42,1.26,4.44,1.78,6.36,.83),
    ("delta_ppc_delta_pyc","R1",.27,.27,1.09,3.62,2.30,7.52,.81), ("delta_ppc_delta_pyc","R2",.26,.26,1.20,3.85,2.23,7.27,.74),
]

MUTATION_ROWS = [
    ("intergenic",3163180,"C","T",None,None,"cg3314|cg3315",181,1.0),
    ("intergenic",2030620,"C","A",None,None,"cg2136|cg2137",125,1.0),
    ("coding",None,"C","T","T306T","synonymous","cg1574",116,1.0),
    ("coding",None,"T","C","L328S","missense","cg1245",102,1.0),
    ("coding",None,"T","C","A22A","synonymous","cg3237",125,1.0),
    ("coding",None,"A","G","Y158H","missense","cg1676",106,.99),
    ("coding",None,"A","C","E324A","missense","cg1451",75,.35),
    ("coding",None,"C","T","G577S","missense","cg0766",148,.33),
    ("coding",None,"G","A","P583S","missense","cg0766",149,.14),
    ("coding",None,"A","T","L508*","stop_gained","cg2267",89,.11),
    ("coding",None,"A","T","V117E","missense","cg0237",76,.11),
]


def import_pxd022622(cursor: sqlite3.Cursor) -> None:
    cursor.executescript("""
        DROP TABLE IF EXISTS pxd022622_peptides;
        DROP TABLE IF EXISTS pxd022622_proteins;
        DROP TABLE IF EXISTS pxd022622_samples;
        DROP TABLE IF EXISTS pxd022622_phenotypes;
        DROP TABLE IF EXISTS pxd022622_variants;
        DROP TABLE IF EXISTS omics_data_availability;
        CREATE TABLE pxd022622_proteins (
            protein_accession TEXT PRIMARY KEY, protein_name TEXT,
            unique_peptide_count INTEGER NOT NULL, transition_count INTEGER NOT NULL,
            max_library_confidence REAL, source TEXT NOT NULL
        );
        CREATE TABLE pxd022622_peptides (
            protein_accession TEXT NOT NULL REFERENCES pxd022622_proteins(protein_accession),
            stripped_sequence TEXT NOT NULL, modification_sequence TEXT,
            precursor_charge INTEGER, retention_time REAL, irt REAL,
            max_relative_intensity REAL, transition_count INTEGER NOT NULL,
            shared INTEGER, decoy INTEGER, source TEXT NOT NULL,
            PRIMARY KEY(protein_accession, stripped_sequence, modification_sequence, precursor_charge)
        );
        CREATE TABLE pxd022622_samples (
            filename TEXT PRIMARY KEY, sample_index_1 INTEGER, sample_index_2 INTEGER,
            sample_index_3 INTEGER, file_type TEXT NOT NULL,
            condition_mapping_status TEXT NOT NULL, source_url TEXT NOT NULL
        );
        CREATE TABLE pxd022622_phenotypes (
            strain TEXT NOT NULL, replicate TEXT NOT NULL,
            growth_rate_bv_h REAL, growth_rate_cdw_h REAL,
            glucose_uptake_bv REAL, glucose_uptake_cdw REAL,
            co2_formation_bv REAL, co2_formation_cdw REAL,
            carbon_balance REAL, condition TEXT NOT NULL, source TEXT NOT NULL,
            PRIMARY KEY(strain, replicate)
        );
        CREATE TABLE pxd022622_variants (
            variant_id INTEGER PRIMARY KEY, region_type TEXT NOT NULL,
            genome_position INTEGER, reference_allele TEXT, alternate_allele TEXT,
            amino_acid_change TEXT, consequence TEXT, affected_locus TEXT,
            supporting_reads INTEGER, allele_frequency REAL,
            strain TEXT NOT NULL, source TEXT NOT NULL
        );
        CREATE TABLE omics_data_availability (
            dataset_id TEXT NOT NULL, data_type TEXT NOT NULL,
            access_status TEXT NOT NULL, local_import_status TEXT NOT NULL,
            source_url TEXT NOT NULL, notes TEXT,
            PRIMARY KEY(dataset_id, data_type)
        );
        CREATE INDEX idx_pxd_peptide_accession ON pxd022622_peptides(protein_accession);
    """)
    proteins: dict[str, dict] = {}
    peptides: dict[tuple, dict] = {}
    library_path = PXD_DIR / "CorynemanualDB200320.txt"
    with library_path.open(encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = str(row.get("uniprot_id") or "").strip()
            sequence = str(row.get("stripped_sequence") or "").strip()
            if not accession or not sequence:
                continue
            name = clean(row.get("protein_name"))
            try:
                confidence = float(row.get("confidence") or 0)
            except ValueError:
                confidence = 0.0
            protein = proteins.setdefault(accession, {"name": name, "peptides": set(), "transitions": 0, "confidence": 0.0})
            protein["peptides"].add(sequence)
            protein["transitions"] += 1
            protein["confidence"] = max(protein["confidence"], confidence)
            try:
                charge = int(float(row.get("prec_z") or 0)) or None
            except ValueError:
                charge = None
            modification = str(row.get("modification_sequence") or "").strip()
            key = (accession, sequence, modification, charge)
            peptide = peptides.setdefault(key, {
                "rt": None, "irt": None, "intensity": None, "transitions": 0,
                "shared": 0, "decoy": 0,
            })
            for field, target in (("RT_detected", "rt"), ("iRT", "irt"), ("relative_intensity", "intensity")):
                try:
                    value = float(row.get(field) or "")
                    peptide[target] = max(peptide[target], value) if peptide[target] is not None else value
                except ValueError:
                    pass
            peptide["transitions"] += 1
            peptide["shared"] = max(peptide["shared"], int(str(row.get("shared") or "").upper() == "TRUE"))
            peptide["decoy"] = max(peptide["decoy"], int(str(row.get("decoy") or "").upper() == "TRUE"))
    cursor.executemany("INSERT INTO pxd022622_proteins VALUES (?,?,?,?,?,?)", [
        (accession, item["name"], len(item["peptides"]), item["transitions"], item["confidence"],
         "PXD022622 CorynemanualDB200320 spectral library")
        for accession, item in proteins.items()
    ])
    cursor.executemany("INSERT INTO pxd022622_peptides VALUES (?,?,?,?,?,?,?,?,?,?,?)", [
        (accession, sequence, modification, charge, item["rt"], item["irt"], item["intensity"],
         item["transitions"], item["shared"], item["decoy"], "PXD022622 spectral library")
        for (accession, sequence, modification, charge), item in peptides.items()
    ])

    sample_rows = []
    with (PXD_DIR / "README.txt").open(encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            name = str(row.get("NAME") or "").strip()
            match = re.fullmatch(r"Proteom\s+(\d+)\.(\d+)\.(\d+)\.wiff", name, re.I)
            if match:
                sample_rows.append((name, *(int(v) for v in match.groups()), "RAW",
                                    "opaque_repository_code; no public code-to-strain mapping",
                                    str(row.get("URI") or "").replace("ftp://", "https://")))
    cursor.executemany("INSERT INTO pxd022622_samples VALUES (?,?,?,?,?,?,?)", sample_rows)
    cursor.executemany("INSERT INTO pxd022622_phenotypes VALUES (?,?,?,?,?,?,?,?,?,?,?)", [
        (*row, "CGXII with 1% D-glucose; exponential growth", "PMC7855459 Supplementary Table S1")
        for row in PHENOTYPE_ROWS
    ])
    cursor.executemany("""
        INSERT INTO pxd022622_variants
        (region_type, genome_position, reference_allele, alternate_allele,
         amino_acid_change, consequence, affected_locus, supporting_reads,
         allele_frequency, strain, source)
        VALUES (?,?,?,?,?,?,?,?,?,?,?)
    """, [(*row, "evolved_delta_ppc_delta_pyc", "PMC7855459 Supplementary Table S6") for row in MUTATION_ROWS])
    cursor.executemany("INSERT INTO omics_data_availability VALUES (?,?,?,?,?,?)", [
        ("PXD022622", "proteomics_raw", "public", "sample inventory imported; raw files not duplicated",
         "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD022622-1&test=no",
         "200 vendor RAW files; condition codes are kept opaque because a public mapping was not supplied"),
        ("PXD022622", "proteomics_spectral_library", "public", "protein and peptide evidence imported",
         "https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/09/PXD022622/",
         "Spectral evidence is not condition-level quantitative abundance"),
        ("PMC7855459", "quantitative_metabolomics", "available_on_request", "study metadata only",
         "https://pmc.ncbi.nlm.nih.gov/articles/PMC7855459/",
         "Article states that non-proteomics datasets are available from the corresponding author; plotted values were not digitized"),
        ("PRJNA678589", "whole_genome_sequencing", "public", "run metadata and reported variants imported",
         "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA678589",
         "SRR13069938; 11 reported variants imported from Supplementary Table S6"),
    ])


def import_integrated_omics() -> None:
    connection = sqlite3.connect(DB_PATH)
    connection.execute("PRAGMA foreign_keys=ON")
    connection.execute("PRAGMA temp_store=MEMORY")
    try:
        with connection:
            cursor = connection.cursor()
            print("[integrated 1/4] Identifier crosswalk and genome releases...")
            import_identifier_layer(cursor)
            print("[integrated 2/4] GEO GSE171302 expression compendium...")
            import_expression_layer(cursor)
            print("[integrated 3/4] TSS, UTR, RBS, transcription units and terminators...")
            import_transcript_structure(cursor)
            print("[integrated 4/4] PXD022622 spectral evidence, phenotypes and variants...")
            import_pxd022622(cursor)
        connection.execute("ANALYZE")
    finally:
        connection.close()
    print(f"SUCCESS: imported integrated public omics layers into {DB_PATH}")


if __name__ == "__main__":
    import_integrated_omics()
