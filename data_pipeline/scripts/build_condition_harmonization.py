#!/usr/bin/env python3
"""Build a conservative, provenance-aware cross-omics condition layer.

GEO and iModulon records are contrasts, PXD phenotypes are single-arm
measurements, and public PXD raw-file codes currently have no strain map. This
module preserves those distinctions instead of guessing missing links.
"""

from __future__ import annotations

import hashlib
import json
import re
import sqlite3
import unicodedata
from collections import defaultdict
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[2]
DB_PATH = ROOT_DIR / "data" / "reference" / "cgl_regulation.db"
PARSER_VERSION = "condition_parser_v1"
SUBSTRATES = (
    "glucose", "acetate", "pyruvate", "lactate", "fructose", "sucrose",
    "glycerol", "citrate", "gluconate", "ribose", "arabinose", "xylose",
    "methanol", "glutamate", "glutamine", "isoleucine", "leucine", "valine",
    "serine", "gaba", "propionate", "benzoate", "adipate", "oxoadipate",
)


def normalized_text(value: str | None) -> str:
    text = unicodedata.normalize("NFKC", value or "")
    text = text.replace("∆", "delta ").replace("Δ", "delta ")
    return re.sub(r"\s+", " ", text).strip().lower()


def normalized_contrast(value: str | None) -> str:
    text = normalized_text(value)
    text = re.sub(r"(?:[-_ ]?rep(?:licate)?[-_ ]?\d+)\s*$", "", text, flags=re.I)
    return re.sub(r"\s+", " ", text).strip()


def condition_signature(value: str | None) -> str:
    """Normalize formatting and dates while retaining biological qualifiers."""
    text = normalized_text(value)
    text = re.sub(r",?\s*date of experiment:\s*[0-9/\-]+(?:\s*/\s*[0-9/\-]+)*", "", text)
    text = text.replace("harvested during", "").replace("harvest during", "")
    text = re.sub(r"\s*,\s*", ", ", text)
    return re.sub(r"\s+", " ", text).strip(" ,")


def stable_id(prefix: str, *parts: object) -> str:
    payload = "\x1f".join("" if part is None else str(part) for part in parts)
    return f"{prefix}_{hashlib.sha1(payload.encode('utf-8')).hexdigest()[:16]}"


def normalize_genotype(value: str | None, source_text: str | None = None) -> str | None:
    raw = normalized_text(value)
    if not raw:
        raw = normalized_text(source_text).split(" harvested", 1)[0].split(" harvest", 1)[0]
        raw = raw.split(" (", 1)[0]
    if raw in {"wt", "wild type", "wild-type", "atcc13032", "atcc 13032"}:
        return "ATCC13032"
    raw = re.sub(r"\bwild[- ]?type\b", "ATCC13032", raw, flags=re.I)
    raw = re.sub(r"\bwt\b", "ATCC13032", raw, flags=re.I)
    raw = re.sub(r"\bdelta\s+", "delta_", raw)
    raw = re.sub(r"\s+", " ", raw).strip(" ,;()")
    return raw or None


def first_float(pattern: str, text: str) -> float | None:
    match = re.search(pattern, text, re.I)
    if not match:
        return None
    try:
        return float(match.group(1).replace(",", "."))
    except ValueError:
        return None


def parse_condition(source_text: str | None, genotype_hint: str | None = None) -> dict[str, object]:
    text = normalized_text(source_text)
    genotype = normalize_genotype(genotype_hint, source_text)
    medium = None
    for label, pattern in (
        ("CGXII", r"\bcgxii\b"), ("CGXIII", r"\bcgxiii\b"),
        ("BHIS", r"\bbhis\b"), ("BHI", r"\bbhi\b"), ("LB", r"\blb\b"),
    ):
        if re.search(pattern, text, re.I):
            medium = label
            break
    substrates = sorted({item for item in SUBSTRATES if re.search(rf"\b{re.escape(item)}\b", text)})
    growth_phase = None
    for label, pattern in (
        ("early_exponential", r"early exponential"),
        ("late_exponential", r"late exponential"),
        ("early_stationary", r"early stationary"),
        ("stationary", r"stationary"),
        ("production", r"production phase"),
        ("exponential", r"exponential|growth phase"),
    ):
        if re.search(pattern, text):
            growth_phase = label
            break
    oxygen = None
    if re.search(r"anaerob|oxygen deprivation", text):
        oxygen = "anaerobic"
    elif re.search(r"microaerob", text):
        oxygen = "microaerobic"
    elif re.search(r"\baerob", text):
        oxygen = "aerobic"
    temperature = first_float(r"(-?\d+(?:[.,]\d+)?)\s*°?c\b", text)
    ph = first_float(r"\bph\s*([0-9]+(?:[.,][0-9]+)?)", text)
    od_match = re.search(r"\bod\s*([0-9]+(?:[.,][0-9]+)?)", text)
    optical_density = od_match.group(1).replace(",", ".") if od_match else None
    time_min = first_float(r"(?:after|harvest after|,\s*)([0-9]+(?:[.,][0-9]+)?)\s*min", text)
    if time_min is None:
        hours = first_float(r"(?:after|harvest after|,\s*)([0-9]+(?:[.,][0-9]+)?)\s*h(?:ours?)?\b", text)
        time_min = hours * 60.0 if hours is not None else None
    perturbations: list[str] = []
    if genotype and ("delta_" in genotype or "::" in genotype):
        perturbations.append("genetic_perturbation")
    for label, pattern in (
        ("IPTG induction", r"\biptg\b"),
        ("phosphate limitation", r"phosphate limitation|w/o.*kh2po4"),
        ("nitrogen limitation", r"w/o.*(?:urea|nh4)"),
        ("iron limitation", r"iron limitation|\b1\s*μm\s*feso4"),
        ("oxygen perturbation", r"anaerob|microaerob|oxygen deprivation"),
        ("temperature perturbation", r"cold shock|heat shock|\d+\s*°?c"),
        ("antibiotic or chemical treatment", r"kasugamycin|ethambutol|bacitracin|mitomycin|glyphosate|hydroperoxide"),
    ):
        if re.search(pattern, text, re.I):
            perturbations.append(label)
    identified = sum(value is not None and value != [] for value in (
        genotype, medium, substrates, growth_phase, oxygen, temperature, ph,
        optical_density, perturbations, time_min,
    ))
    confidence = "high" if identified >= 4 else "medium" if identified >= 2 else "low"
    canonical = {
        "organism": "Corynebacterium glutamicum ATCC 13032",
        "genotype": genotype, "medium": medium, "substrates": substrates,
        "growth_phase": growth_phase, "oxygen_regime": oxygen,
        "temperature_c": temperature, "ph": ph, "optical_density": optical_density,
        "perturbations": perturbations, "treatment_time_min": time_min,
        "condition_signature": condition_signature(source_text),
    }
    condition_id = stable_id("cond", json.dumps(canonical, sort_keys=True, ensure_ascii=False))
    label_parts = [genotype or "unspecified genotype", medium]
    if substrates:
        label_parts.append("+".join(substrates))
    label_parts.extend([growth_phase, oxygen])
    canonical.update(
        condition_id=condition_id,
        canonical_label=" | ".join(str(item) for item in label_parts if item),
        mapping_confidence=confidence,
    )
    return canonical


def split_contrast(value: str) -> tuple[str, str | None]:
    parts = re.split(r"\s+vs\.?\s+", value, maxsplit=1, flags=re.I)
    return parts[0].strip(), parts[1].strip() if len(parts) == 2 else None


def build_condition_harmonization(db_path: Path = DB_PATH) -> None:
    print(f"Building cross-omics condition harmonization layer in: {db_path}")
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys=ON")
    cursor = conn.cursor()
    cursor.executescript("""
        DROP VIEW IF EXISTS v_condition_phenotype_summary;
        DROP VIEW IF EXISTS v_expression_imodulon_condition_link;
        DROP VIEW IF EXISTS v_condition_multiomics_coverage;
        DROP TABLE IF EXISTS omics_sample_condition_map;
        DROP TABLE IF EXISTS omics_comparison_map;
        DROP TABLE IF EXISTS condition_mapping_issues;
        DROP TABLE IF EXISTS condition_aliases;
        DROP TABLE IF EXISTS condition_comparisons;
        DROP TABLE IF EXISTS standard_conditions;
        CREATE TABLE standard_conditions (
            condition_id TEXT PRIMARY KEY, organism TEXT NOT NULL, genotype TEXT,
            medium TEXT, substrates_json TEXT NOT NULL, growth_phase TEXT,
            oxygen_regime TEXT, temperature_c REAL, ph REAL, optical_density TEXT,
            perturbations_json TEXT NOT NULL, treatment_time_min REAL,
            condition_signature TEXT NOT NULL, canonical_label TEXT NOT NULL, mapping_confidence TEXT NOT NULL,
            parser_version TEXT NOT NULL
        );
        CREATE TABLE condition_comparisons (
            comparison_id TEXT PRIMARY KEY,
            test_condition_id TEXT REFERENCES standard_conditions(condition_id),
            reference_condition_id TEXT REFERENCES standard_conditions(condition_id),
            canonical_label TEXT NOT NULL, comparison_kind TEXT NOT NULL
        );
        CREATE TABLE condition_aliases (
            alias_id TEXT PRIMARY KEY,
            condition_id TEXT NOT NULL REFERENCES standard_conditions(condition_id),
            source_system TEXT NOT NULL, source_text TEXT NOT NULL,
            normalized_text TEXT NOT NULL, mapping_method TEXT NOT NULL,
            mapping_confidence TEXT NOT NULL
        );
        CREATE TABLE omics_comparison_map (
            map_id TEXT PRIMARY KEY, dataset_id TEXT NOT NULL,
            source_record_id TEXT NOT NULL, omics_type TEXT NOT NULL,
            comparison_id TEXT REFERENCES condition_comparisons(comparison_id),
            source_label TEXT, mapping_method TEXT NOT NULL,
            mapping_confidence TEXT NOT NULL
        );
        CREATE TABLE omics_sample_condition_map (
            map_id TEXT PRIMARY KEY, dataset_id TEXT NOT NULL,
            source_record_id TEXT NOT NULL, omics_type TEXT NOT NULL,
            condition_role TEXT NOT NULL,
            condition_id TEXT REFERENCES standard_conditions(condition_id),
            comparison_id TEXT REFERENCES condition_comparisons(comparison_id),
            mapping_status TEXT NOT NULL, notes TEXT
        );
        CREATE TABLE condition_mapping_issues (
            issue_id TEXT PRIMARY KEY, issue_type TEXT NOT NULL,
            source_system TEXT NOT NULL, source_label TEXT,
            affected_records_json TEXT NOT NULL, details TEXT NOT NULL,
            resolution_status TEXT NOT NULL
        );
    """)
    conditions: dict[str, dict[str, object]] = {}
    comparisons: dict[str, tuple] = {}
    aliases: dict[str, tuple] = {}
    comparison_maps: dict[str, tuple] = {}
    sample_maps: dict[str, tuple] = {}
    issues: dict[str, tuple] = {}
    geo_contrast_index: dict[str, str] = {}
    geo_label_comparisons: dict[str, set[str]] = defaultdict(set)
    geo_label_records: dict[str, list[str]] = defaultdict(list)

    def register_condition(source_system: str, source_text: str | None, genotype: str | None = None) -> str:
        parsed = parse_condition(source_text, genotype)
        condition_id = str(parsed["condition_id"])
        conditions.setdefault(condition_id, parsed)
        if source_text:
            alias_id = stable_id("alias", source_system, normalized_text(source_text), parsed.get("genotype"))
            aliases[alias_id] = (
                alias_id, condition_id, source_system, source_text, normalized_text(source_text),
                PARSER_VERSION, parsed["mapping_confidence"],
            )
        return condition_id

    def register_comparison(test_id: str, ref_id: str | None, label: str, kind: str) -> str:
        comparison_id = stable_id("cmp", test_id, ref_id, normalized_contrast(label), kind)
        comparisons.setdefault(comparison_id, (comparison_id, test_id, ref_id, label, kind))
        return comparison_id

    for row in cursor.execute("""
        SELECT sample_id, dataset_id, title, genotype_ch1, genotype_ch2,
               condition_ch1, condition_ch2 FROM expression_samples
    """):
        sample_id, dataset_id, title, genotype1, genotype2, text1, text2 = row
        test_id = register_condition(dataset_id, text1 or title, genotype1)
        ref_id = register_condition(dataset_id, text2, genotype2) if text2 else None
        comparison_id = register_comparison(test_id, ref_id, normalized_contrast(title), "two_channel_contrast")
        contrast_key = normalized_contrast(title)
        geo_contrast_index.setdefault(contrast_key, comparison_id)
        geo_label_comparisons[contrast_key].add(comparison_id)
        geo_label_records[contrast_key].append(sample_id)
        map_id = stable_id("cmap", dataset_id, sample_id)
        comparison_maps[map_id] = (
            map_id, dataset_id, sample_id, "transcriptomics", comparison_id,
            title, "parsed_GEO_channels", "high",
        )
        for role, condition_id in (("test", test_id), ("reference", ref_id)):
            if condition_id:
                arm_id = stable_id("smap", dataset_id, sample_id, role)
                sample_maps[arm_id] = (
                    arm_id, dataset_id, sample_id, "transcriptomics", role,
                    condition_id, comparison_id, "mapped", None,
                )

    for key, comparison_ids in geo_label_comparisons.items():
        if len(comparison_ids) > 1:
            issue_id = stable_id("issue", "conflicting_GEO_arms", key)
            issues[issue_id] = (
                issue_id, "conflicting_GEO_arms", "GSE171302", key,
                json.dumps(geo_label_records[key], ensure_ascii=False),
                "Replicates share a contrast title but their channel metadata resolve to different conditions; kept separate.",
                "open_source_metadata_conflict",
            )

    for source_id, label in cursor.execute("""
        SELECT sample_id, condition_name FROM imodulon_condition_activities
        GROUP BY sample_id, condition_name
    """):
        key = normalized_contrast(label)
        comparison_id = geo_contrast_index.get(key)
        if comparison_id:
            method, confidence = "normalized_GEO_title_match", "high"
        else:
            test_text, reference_text = split_contrast(label)
            test_id = register_condition("ICA_GSE171302", test_text)
            ref_id = register_condition("ICA_GSE171302", reference_text) if reference_text else None
            comparison_id = register_comparison(test_id, ref_id, key, "inferred_contrast")
            method, confidence = "parsed_contrast_only", "medium"
            issue_id = stable_id("issue", "ICA_without_GEO_title_match", source_id)
            issues[issue_id] = (
                issue_id, "ICA_without_GEO_title_match", "ICA_GSE171302", label,
                json.dumps([str(source_id)], ensure_ascii=False),
                "No normalized GEO sample-title match; contrast was parsed but not linked to an expression replicate.",
                "retained_unlinked",
            )
        map_id = stable_id("cmap", "ICA_GSE171302", source_id)
        comparison_maps[map_id] = (
            map_id, "ICA_GSE171302", str(source_id), "imodulon_activity",
            comparison_id, label, method, confidence,
        )

    for strain, replicate, text in cursor.execute(
        "SELECT strain, replicate, condition FROM pxd022622_phenotypes"
    ):
        condition_id = register_condition("PXD022622_phenotype", text, strain)
        source_id = f"{strain}:{replicate}"
        map_id = stable_id("smap", "PXD022622_phenotype", source_id)
        sample_maps[map_id] = (
            map_id, "PXD022622", source_id, "phenotype", "single",
            condition_id, None, "mapped", "Supplementary Table S1",
        )
    for filename, status in cursor.execute(
        "SELECT filename, condition_mapping_status FROM pxd022622_samples"
    ):
        map_id = stable_id("smap", "PXD022622_raw", filename)
        sample_maps[map_id] = (
            map_id, "PXD022622", filename, "proteomics_raw", "single",
            None, None, "unmapped", status,
        )
    issue_id = stable_id("issue", "PXD_sample_sheet_missing", "PXD022622")
    issues[issue_id] = (
        issue_id, "PXD_sample_sheet_missing", "PXD022622", None,
        json.dumps(["200 public vendor raw files"], ensure_ascii=False),
        "Repository file codes cannot be assigned to strains or conditions without the requested author sample sheet.",
        "awaiting_external_metadata",
    )

    cursor.executemany("INSERT INTO standard_conditions VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", [(
        item["condition_id"], item["organism"], item["genotype"], item["medium"],
        json.dumps(item["substrates"], ensure_ascii=False), item["growth_phase"],
        item["oxygen_regime"], item["temperature_c"], item["ph"], item["optical_density"],
        json.dumps(item["perturbations"], ensure_ascii=False), item["treatment_time_min"],
        item["condition_signature"], item["canonical_label"], item["mapping_confidence"], PARSER_VERSION,
    ) for item in conditions.values()])
    cursor.executemany("INSERT INTO condition_comparisons VALUES (?,?,?,?,?)", comparisons.values())
    cursor.executemany("INSERT INTO condition_aliases VALUES (?,?,?,?,?,?,?)", aliases.values())
    cursor.executemany("INSERT INTO omics_comparison_map VALUES (?,?,?,?,?,?,?,?)", comparison_maps.values())
    cursor.executemany("INSERT INTO omics_sample_condition_map VALUES (?,?,?,?,?,?,?,?,?)", sample_maps.values())
    cursor.executemany("INSERT INTO condition_mapping_issues VALUES (?,?,?,?,?,?,?)", issues.values())
    cursor.executescript("""
        CREATE INDEX idx_condition_alias_condition ON condition_aliases(condition_id);
        CREATE INDEX idx_condition_alias_normalized ON condition_aliases(normalized_text);
        CREATE INDEX idx_comparison_pair ON condition_comparisons(test_condition_id, reference_condition_id);
        CREATE INDEX idx_omics_comparison ON omics_comparison_map(comparison_id, omics_type);
        CREATE INDEX idx_omics_condition ON omics_sample_condition_map(condition_id, omics_type);
        CREATE INDEX idx_omics_source ON omics_sample_condition_map(dataset_id, source_record_id);
        CREATE INDEX idx_condition_issues_type ON condition_mapping_issues(issue_type, resolution_status);
        CREATE VIEW v_condition_multiomics_coverage AS
        SELECT c.condition_id, c.canonical_label, c.genotype, c.medium,
               c.substrates_json, c.growth_phase, c.oxygen_regime,
               SUM(CASE WHEN m.omics_type='transcriptomics' THEN 1 ELSE 0 END) AS transcriptomic_arms,
               SUM(CASE WHEN m.omics_type='phenotype' THEN 1 ELSE 0 END) AS phenotype_samples,
               SUM(CASE WHEN m.omics_type='proteomics_raw' THEN 1 ELSE 0 END) AS proteomics_raw_files
        FROM standard_conditions c
        LEFT JOIN omics_sample_condition_map m ON m.condition_id=c.condition_id
        GROUP BY c.condition_id;
        CREATE VIEW v_expression_imodulon_condition_link AS
        SELECT geo.comparison_id, geo.dataset_id AS expression_dataset,
               geo.source_record_id AS expression_sample_id,
               ica.source_record_id AS imodulon_condition_id,
               ica.source_label AS condition_name,
               ica.mapping_method, ica.mapping_confidence
        FROM omics_comparison_map geo
        JOIN omics_comparison_map ica ON ica.comparison_id=geo.comparison_id
        WHERE geo.omics_type='transcriptomics' AND ica.omics_type='imodulon_activity';
        CREATE VIEW v_condition_phenotype_summary AS
        SELECT m.condition_id, c.canonical_label, p.strain,
               COUNT(*) AS replicate_count,
               AVG(p.growth_rate_cdw_h) AS mean_growth_rate_cdw_h,
               AVG(p.glucose_uptake_cdw) AS mean_glucose_uptake_cdw,
               AVG(p.co2_formation_cdw) AS mean_co2_formation_cdw,
               AVG(p.carbon_balance) AS mean_carbon_balance
        FROM pxd022622_phenotypes p
        JOIN omics_sample_condition_map m
          ON m.dataset_id='PXD022622' AND m.omics_type='phenotype'
         AND m.source_record_id=(p.strain || ':' || p.replicate)
        JOIN standard_conditions c ON c.condition_id=m.condition_id
        GROUP BY m.condition_id, p.strain;
    """)
    conn.commit()
    print(
        "SUCCESS: condition layer contains "
        f"{len(conditions)} standardized conditions, {len(comparisons)} contrasts, "
        f"{len(comparison_maps)} comparison mappings, {len(sample_maps)} sample-arm mappings, "
        f"and {len(issues)} auditable mapping issues"
    )
    conn.close()


if __name__ == "__main__":
    build_condition_harmonization()
