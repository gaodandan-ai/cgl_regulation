#!/usr/bin/env python3
"""Build TF-target edge evidence features for confidence modeling.

This script creates a flat feature table that can later be used to train a
RandomForestClassifier or another edge-confidence model. It intentionally does
not train a model yet.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import math
import random
import re
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
DEFAULT_OUTPUT = DATA_DIR / "reference" / "edge_confidence" / "tf_gene_edge_features.csv"

# Paths for new integration data
IMODULON_WEIGHTS_PATH = DATA_DIR / "reference" / "imodulon" / "imodulon_gene_weights.json"
IMODULON_BY_GENE_PATH = DATA_DIR / "reference" / "imodulon" / "imodulon_by_gene.json"
TCS_SYSTEMS_PATH      = DATA_DIR / "reference" / "tcs_systems.json"
SIGMA_ANNOT_PATH      = DATA_DIR / "reference" / "sigma_factor_annotations.json"
DEFAULT_COMPENDIUM_CORRELATIONS = DATA_DIR / "reference" / "expression_compendium" / "tf_target_compendium_correlations.csv"


FEATURE_COLUMNS = [
    "tf_locus",
    "tf_name",
    "target_locus",
    "target_name",
    "label",
    "sample_type",
    "regulation_role",
    "is_activation",
    "is_repression",
    "is_sigma_factor",
    "has_binding_site",
    "binding_site_length",
    "binding_site_gc_fraction",
    "evidence_text",
    "evidence_score",
    "has_experimental_evidence",
    "has_pmid",
    "pmid_count",
    "source",
    "has_database_source",
    "source_coryneregnet",
    "target_operon",
    "target_in_operon",
    "target_operon_size",
    "tf_target_same_operon",
    "target_has_srna_prediction",
    "target_best_srna_rank",
    "target_best_srna_copra_fdr",
    "target_best_srna_energy",
    "target_srna_count_rank_le_10",
    "target_mapped_reaction_count",
    "target_mapped_pathway_count",
    "target_enzyme_constrained_reaction_count",
    "target_has_enzyme_constraint",
    "target_ec_number_count",
    "target_uniprot_count",
    "target_kcat_median",
    "target_molecular_weight_median",
    "target_kcat_mw_median",
    "tf_mapped_reaction_count",
    "tf_target_share_metabolic_pathway",
    "tf_target_shared_pathway_count",
    "expression_feature_available",
    "expression_correlation",
    "expression_pvalue",
    "expression_abs_correlation",
    "expression_sample_count",
    "expression_source",
    "feature_missing_count",
    # --- iModulon / TCS / Sigma features (new) ---
    "target_imodulon_count",
    "tf_imodulon_count",
    "imodulon_coactivation_score",
    "tf_is_tcs_regulator",
    "tf_sigma_class",
]


def clean(value: Any) -> str:
    return str(value or "").strip()


def safe_float(value: Any) -> Optional[float]:
    try:
        if value is None or value == "":
            return None
        parsed = float(value)
        return parsed if math.isfinite(parsed) else None
    except Exception:
        return None


def safe_int(value: Any) -> Optional[int]:
    parsed = safe_float(value)
    return int(parsed) if parsed is not None else None


def format_optional(value: Optional[float], digits: int = 6) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}g}"


def gc_fraction(seq: str) -> str:
    bases = [b.upper() for b in seq if b.upper() in {"A", "C", "G", "T"}]
    if not bases:
        return ""
    return format_optional(sum(1 for b in bases if b in {"G", "C"}) / len(bases), 4)


def normalize_locus(value: str, cg_to_cgl: Dict[str, str], cgl_to_cg: Dict[str, str], name_to_cg: Dict[str, str]) -> str:
    raw = clean(value).lower()
    if not raw:
        return ""
    raw = raw.replace("gene:", "")
    if raw in cgl_to_cg:
        return cgl_to_cg[raw]
    if raw in name_to_cg:
        return name_to_cg[raw]
    return raw


def alias_set(value: str, cg_to_cgl: Dict[str, str], cgl_to_cg: Dict[str, str], name_to_cg: Dict[str, str]) -> set[str]:
    normalized = normalize_locus(value, cg_to_cgl, cgl_to_cg, name_to_cg)
    aliases = {clean(value).lower(), normalized}
    if normalized in cg_to_cgl:
        aliases.add(cg_to_cgl[normalized].lower())
    if normalized in cgl_to_cg:
        aliases.add(cgl_to_cg[normalized].lower())
    return {a for a in aliases if a}


def read_gene_mapping(path: Path) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str], Dict[str, str], Dict[str, str]]:
    cg_to_cgl: Dict[str, str] = {}
    cgl_to_cg: Dict[str, str] = {}
    name_to_cg: Dict[str, str] = {}
    gene_name: Dict[str, str] = {}
    product: Dict[str, str] = {}
    if not path.exists():
        return cg_to_cgl, cgl_to_cg, name_to_cg, gene_name, product

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            cg = clean(row.get("cg_locus")).lower()
            cgl = clean(row.get("cgl_locus")).lower()
            name = clean(row.get("gene_name"))
            prod = clean(row.get("product"))
            if cg and cgl:
                cg_to_cgl[cg] = cgl
                cgl_to_cg[cgl] = cg
            if cg and name:
                name_to_cg[name.lower()] = cg
                gene_name[cg] = name
            if cg and prod:
                product[cg] = prod
    return cg_to_cgl, cgl_to_cg, name_to_cg, gene_name, product


def read_regulations(path: Path, cg_to_cgl: Dict[str, str], cgl_to_cg: Dict[str, str], name_to_cg: Dict[str, str]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            tf = normalize_locus(row.get("TF_locusTag", ""), cg_to_cgl, cgl_to_cg, name_to_cg)
            target = normalize_locus(row.get("TG_locusTag", ""), cg_to_cgl, cgl_to_cg, name_to_cg)
            if not tf or not target:
                continue
            item = {k: clean(v) for k, v in row.items() if k is not None}
            item["tf_locus"] = tf
            item["target_locus"] = target
            rows.append(item)
    return rows


def read_operons(path: Path, cg_to_cgl: Dict[str, str], cgl_to_cg: Dict[str, str], name_to_cg: Dict[str, str]) -> Tuple[Dict[str, str], Dict[str, int]]:
    gene_to_operon: Dict[str, str] = {}
    operon_size: Dict[str, int] = {}
    if not path.exists():
        return gene_to_operon, operon_size
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        next(reader, None)
        for row in reader:
            if len(row) < 3:
                continue
            operon = clean(row[0]).lstrip(">")
            genes = [normalize_locus(g, cg_to_cgl, cgl_to_cg, name_to_cg) for g in row[2:] if clean(g)]
            genes = [g for g in genes if g]
            if not operon or not genes:
                continue
            operon_size[operon] = len(genes)
            for gene in genes:
                gene_to_operon[gene] = operon
    return gene_to_operon, operon_size


def read_srna_predictions(path: Path, cg_to_cgl: Dict[str, str], cgl_to_cg: Dict[str, str], name_to_cg: Dict[str, str]) -> Dict[str, Dict[str, Any]]:
    by_target: Dict[str, Dict[str, Any]] = {}
    if not path.exists():
        return by_target
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            target = normalize_locus(row.get("mrna", ""), cg_to_cgl, cgl_to_cg, name_to_cg)
            if not target:
                continue
            rank = safe_int(row.get("rank"))
            fdr = safe_float(row.get("copra_fdr"))
            energy = safe_float(row.get("energy"))
            stat = by_target.setdefault(target, {
                "count_rank_le_10": 0,
                "best_rank": None,
                "best_fdr": None,
                "best_energy": None,
            })
            if rank is not None and rank <= 10:
                stat["count_rank_le_10"] += 1
            if rank is not None and (stat["best_rank"] is None or rank < stat["best_rank"]):
                stat["best_rank"] = rank
                stat["best_fdr"] = fdr
                stat["best_energy"] = energy
    return by_target


def load_run_server_module() -> Optional[Any]:
    path = ROOT / "run_server.py"
    if not path.exists():
        return None
    root_string = str(ROOT)
    if root_string not in sys.path:
        sys.path.insert(0, root_string)
    spec = importlib.util.spec_from_file_location("cgl_run_server_for_features", path)
    if spec is None or spec.loader is None:
        return None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    old_cwd = Path.cwd()
    try:
        import os
        os.chdir(ROOT)
        spec.loader.exec_module(module)
        return module
    except Exception as exc:
        print(f"Warning: failed to import run_server.py for metabolic features: {exc}", file=sys.stderr)
        return None
    finally:
        import os
        os.chdir(old_cwd)


def load_metabolic_features(cg_to_cgl: Dict[str, str], cgl_to_cg: Dict[str, str], name_to_cg: Dict[str, str]) -> Dict[str, Dict[str, Any]]:
    module = load_run_server_module()
    if module is None:
        return {}
    old_cwd = Path.cwd()
    try:
        import os
        os.chdir(ROOT)
        module.METABOLIC_MODEL_CACHE = None
        mapping = module.load_metabolic_model_mappings()
    finally:
        import os
        os.chdir(old_cwd)

    features: Dict[str, Dict[str, Any]] = {}
    for alias, reactions in mapping.get("gene_to_reactions", {}).items():
        gene = normalize_locus(alias, cg_to_cgl, cgl_to_cg, name_to_cg)
        if not gene:
            continue
        stat = features.setdefault(gene, {
            "reaction_ids": set(),
            "pathway_ids": set(),
            "enzyme_reaction_ids": set(),
            "ec_numbers": set(),
            "uniprot_ids": set(),
            "kcat_values": [],
            "mw_values": [],
            "kcat_mw_values": [],
        })
        for reaction in reactions or []:
            model = clean(reaction.get("model")) or "model"
            reaction_id = clean(reaction.get("id"))
            if reaction_id:
                stat["reaction_ids"].add(f"{model}:{reaction_id}")
            pathway = clean(reaction.get("pathway_id")) or clean(reaction.get("pathway_name"))
            if pathway:
                stat["pathway_ids"].add(f"{model}:{pathway}")
            if reaction.get("enzyme_constraint") or reaction.get("kcat") is not None or reaction.get("kcat_MW") is not None:
                if reaction_id:
                    stat["enzyme_reaction_ids"].add(f"{model}:{reaction_id}")
            if reaction.get("ec_number"):
                stat["ec_numbers"].add(clean(reaction.get("ec_number")))
            for uniprot in reaction.get("uniprot_ids") or []:
                if uniprot:
                    stat["uniprot_ids"].add(clean(uniprot))
            for key, bucket in (("kcat", "kcat_values"), ("molecular_weight", "mw_values"), ("kcat_MW", "kcat_mw_values")):
                value = safe_float(reaction.get(key))
                if value is not None:
                    stat[bucket].append(value)
    return features


def read_expression_correlations(path: Optional[Path], cg_to_cgl: Dict[str, str], cgl_to_cg: Dict[str, str], name_to_cg: Dict[str, str]) -> Dict[Tuple[str, str], Dict[str, Optional[float]]]:
    correlations: Dict[Tuple[str, str], Dict[str, Optional[float]]] = {}
    if not path or not path.exists():
        return correlations
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            tf = normalize_locus(row.get("tf") or row.get("TF") or row.get("source") or row.get("regulator"), cg_to_cgl, cgl_to_cg, name_to_cg)
            target = normalize_locus(row.get("target") or row.get("gene") or row.get("TG") or row.get("target_gene"), cg_to_cgl, cgl_to_cg, name_to_cg)
            if not tf or not target:
                continue
            corr = safe_float(row.get("correlation") or row.get("pearson") or row.get("spearman") or row.get("r"))
            pvalue = safe_float(row.get("pvalue") or row.get("p_value") or row.get("padj"))
            abs_corr = safe_float(row.get("abs_correlation") or row.get("abs_pearson") or row.get("abs_r"))
            sample_count = safe_float(row.get("sample_count") or row.get("n_samples") or row.get("condition_count"))
            source = clean(row.get("source") or row.get("dataset"))
            correlations[(tf, target)] = {
                "correlation": corr,
                "pvalue": pvalue,
                "abs_correlation": abs_corr if abs_corr is not None else (abs(corr) if corr is not None else None),
                "sample_count": sample_count,
                "source": source,
            }
    return correlations


def evidence_score(text: str, source: str, pmid: str, binding_site: str) -> float:
    combined = f"{text} {source}".lower()
    score = 0.0
    if "experimental" in combined:
        score += 0.55
    if "chip" in combined:
        score += 0.2
    if "literature" in combined or pmid:
        score += 0.15
    if source:
        score += 0.05
    if binding_site:
        score += 0.05
    return min(score, 1.0)


def median_or_none(values: Sequence[float]) -> Optional[float]:
    return median(values) if values else None


def build_feature_row(
    tf: str,
    target: str,
    label: int,
    sample_type: str,
    regulation_row: Optional[Dict[str, str]],
    gene_name: Dict[str, str],
    gene_to_operon: Dict[str, str],
    operon_size: Dict[str, int],
    srna_features: Dict[str, Dict[str, Any]],
    metabolic_features: Dict[str, Dict[str, Any]],
    expression_correlations: Dict[Tuple[str, str], Dict[str, Optional[float]]],
) -> Dict[str, Any]:
    row = regulation_row or {}
    binding_site = clean(row.get("Binding_site"))
    role = clean(row.get("Role"))
    role_upper = role.upper()
    evidence = clean(row.get("Evidence"))
    pmid = clean(row.get("PMID"))
    source = clean(row.get("Source"))
    target_operon = clean(row.get("Operon")) or gene_to_operon.get(target, "")
    tf_operon = gene_to_operon.get(tf, "")
    srna = srna_features.get(target, {})
    target_met = metabolic_features.get(target, {})
    tf_met = metabolic_features.get(tf, {})
    target_pathways = target_met.get("pathway_ids", set())
    tf_pathways = tf_met.get("pathway_ids", set())
    shared_pathways = target_pathways & tf_pathways
    expr = expression_correlations.get((tf, target), {})

    values = {
        "tf_locus": tf,
        "tf_name": clean(row.get("TF_name")) or gene_name.get(tf, tf),
        "target_locus": target,
        "target_name": clean(row.get("TG_name")) or gene_name.get(target, target),
        "label": label,
        "sample_type": sample_type,
        "regulation_role": role or "unknown",
        "is_activation": 1 if role_upper == "A" else 0,
        "is_repression": 1 if role_upper == "R" else 0,
        "is_sigma_factor": 1 if clean(row.get("Is_sigma_factor")).lower() in {"yes", "true", "1"} else 0,
        "has_binding_site": 1 if binding_site else 0,
        "binding_site_length": len(re.sub(r"[^A-Za-z]", "", binding_site)) if binding_site else 0,
        "binding_site_gc_fraction": gc_fraction(binding_site),
        "evidence_text": evidence,
        "evidence_score": format_optional(evidence_score(evidence, source, pmid, binding_site), 4),
        "has_experimental_evidence": 1 if "experimental" in evidence.lower() else 0,
        "has_pmid": 1 if pmid else 0,
        "pmid_count": len([p for p in re.split(r"[;,| ]+", pmid) if p.strip()]),
        "source": source,
        "has_database_source": 1 if source else 0,
        "source_coryneregnet": 1 if "coryneregnet" in source.lower() else 0,
        "target_operon": target_operon,
        "target_in_operon": 1 if target_operon else 0,
        "target_operon_size": operon_size.get(target_operon, 0) if target_operon else 0,
        "tf_target_same_operon": 1 if tf_operon and target_operon and tf_operon == target_operon else 0,
        "target_has_srna_prediction": 1 if target in srna_features else 0,
        "target_best_srna_rank": srna.get("best_rank", ""),
        "target_best_srna_copra_fdr": format_optional(srna.get("best_fdr"), 6),
        "target_best_srna_energy": format_optional(srna.get("best_energy"), 6),
        "target_srna_count_rank_le_10": srna.get("count_rank_le_10", 0),
        "target_mapped_reaction_count": len(target_met.get("reaction_ids", set())),
        "target_mapped_pathway_count": len(target_pathways),
        "target_enzyme_constrained_reaction_count": len(target_met.get("enzyme_reaction_ids", set())),
        "target_has_enzyme_constraint": 1 if target_met.get("enzyme_reaction_ids") else 0,
        "target_ec_number_count": len(target_met.get("ec_numbers", set())),
        "target_uniprot_count": len(target_met.get("uniprot_ids", set())),
        "target_kcat_median": format_optional(median_or_none(target_met.get("kcat_values", [])), 6),
        "target_molecular_weight_median": format_optional(median_or_none(target_met.get("mw_values", [])), 6),
        "target_kcat_mw_median": format_optional(median_or_none(target_met.get("kcat_mw_values", [])), 6),
        "tf_mapped_reaction_count": len(tf_met.get("reaction_ids", set())),
        "tf_target_share_metabolic_pathway": 1 if shared_pathways else 0,
        "tf_target_shared_pathway_count": len(shared_pathways),
        "expression_feature_available": 1 if expr else 0,
        "expression_correlation": format_optional(expr.get("correlation"), 6),
        "expression_pvalue": format_optional(expr.get("pvalue"), 6),
        "expression_abs_correlation": format_optional(expr.get("abs_correlation"), 6),
        "expression_sample_count": format_optional(expr.get("sample_count"), 6),
        "expression_source": expr.get("source", ""),
    }

    missing = 0
    for key in (
        "binding_site_gc_fraction",
        "target_best_srna_rank",
        "target_best_srna_copra_fdr",
        "target_best_srna_energy",
        "target_kcat_median",
        "target_molecular_weight_median",
        "target_kcat_mw_median",
        "expression_correlation",
    ):
        if values.get(key) in ("", None):
            missing += 1
    values["feature_missing_count"] = missing

    # ── iModulon / TCS / Sigma extra features (injected externally if available) ──
    # These are filled in main() after build_feature_row() returns, so defaults here:
    values.setdefault("target_imodulon_count", "")
    values.setdefault("tf_imodulon_count", "")
    values.setdefault("imodulon_coactivation_score", "")
    values.setdefault("tf_is_tcs_regulator", "")
    values.setdefault("tf_sigma_class", "")

    return values


def stable_pair_hash(tf: str, target: str) -> str:
    return hashlib.sha1(f"{tf}->{target}".encode("utf-8")).hexdigest()[:12]


def make_negative_pairs(positives: set[Tuple[str, str]], tfs: Sequence[str], targets: Sequence[str], ratio: float, seed: int) -> List[Tuple[str, str]]:
    requested = int(len(positives) * max(ratio, 0))
    if requested <= 0:
        return []
    rng = random.Random(seed)
    tfs = sorted(set(tfs))
    targets = sorted(set(targets))
    if not tfs or not targets:
        return []
    negatives: set[Tuple[str, str]] = set()
    max_attempts = max(requested * 50, 1000)
    attempts = 0
    while len(negatives) < requested and attempts < max_attempts:
        attempts += 1
        tf = rng.choice(tfs)
        target = rng.choice(targets)
        if tf == target:
            continue
        pair = (tf, target)
        if pair in positives or pair in negatives:
            continue
        negatives.add(pair)
    return sorted(negatives, key=lambda pair: stable_pair_hash(*pair))


def write_csv(path: Path, rows: Iterable[Dict[str, Any]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FEATURE_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in FEATURE_COLUMNS})
            count += 1
    return count


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build TF-target edge evidence features for confidence modeling.")
    parser.add_argument("--regulations", type=Path, default=DATA_DIR / "regulations.csv")
    parser.add_argument("--gene-mapping", type=Path, default=DATA_DIR / "gene_mapping.csv")
    parser.add_argument("--operons", type=Path, default=DATA_DIR / "operons.csv")
    parser.add_argument("--srna-regulation", type=Path, default=DATA_DIR / "rna_regulation.csv")
    parser.add_argument(
        "--expression-correlations",
        type=Path,
        default=None,
        help=(
            "Optional CSV with tf,target,correlation,pvalue columns. If omitted, "
            "the expression compendium output is used when available."
        ),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--negative-ratio", type=float, default=1.0, help="Number of sampled negative TF-gene pairs per positive edge.")
    parser.add_argument("--seed", type=int, default=13)
    return parser.parse_args()


def load_integration_data() -> Dict[str, Any]:
    """Load iModulon, TCS, and sigma factor data for feature enrichment."""
    import json
    result: Dict[str, Any] = {
        "imodulon_by_gene": {},
        "imodulon_weights": {},
        "tcs_hk": set(),
        "tcs_rr": set(),
        "sigma_by_locus": {},
    }
    try:
        if IMODULON_BY_GENE_PATH.exists():
            with open(IMODULON_BY_GENE_PATH, encoding="utf-8") as f:
                result["imodulon_by_gene"] = json.load(f)
        if IMODULON_WEIGHTS_PATH.exists():
            with open(IMODULON_WEIGHTS_PATH, encoding="utf-8") as f:
                result["imodulon_weights"] = json.load(f)
        if TCS_SYSTEMS_PATH.exists():
            with open(TCS_SYSTEMS_PATH, encoding="utf-8") as f:
                for tcs in json.load(f):
                    if tcs.get("hk_locus"):
                        result["tcs_hk"].add(tcs["hk_locus"].lower())
                    if tcs.get("rr_locus"):
                        result["tcs_rr"].add(tcs["rr_locus"].lower())
        if SIGMA_ANNOT_PATH.exists():
            with open(SIGMA_ANNOT_PATH, encoding="utf-8") as f:
                for key, ann in json.load(f).items():
                    if ann.get("locus"):
                        result["sigma_by_locus"][ann["locus"].lower()] = ann.get("sigma_class", "ECF_sigma")
    except Exception as e:
        print(f"Warning: integration data load partial: {e}")
    return result


def main() -> int:
    args = parse_args()
    if args.expression_correlations is None and DEFAULT_COMPENDIUM_CORRELATIONS.exists():
        args.expression_correlations = DEFAULT_COMPENDIUM_CORRELATIONS
        print(f"Using expression compendium correlations: {args.expression_correlations}")
    cg_to_cgl, cgl_to_cg, name_to_cg, gene_name, _product = read_gene_mapping(args.gene_mapping)
    regulation_rows = read_regulations(args.regulations, cg_to_cgl, cgl_to_cg, name_to_cg)
    gene_to_operon, operon_size = read_operons(args.operons, cg_to_cgl, cgl_to_cg, name_to_cg)
    srna_features = read_srna_predictions(args.srna_regulation, cg_to_cgl, cgl_to_cg, name_to_cg)
    metabolic_features = load_metabolic_features(cg_to_cgl, cgl_to_cg, name_to_cg)
    expression_correlations = read_expression_correlations(args.expression_correlations, cg_to_cgl, cgl_to_cg, name_to_cg)

    # Load new integration data
    integ = load_integration_data()
    imod_by_gene  = integ["imodulon_by_gene"]
    imod_weights  = integ["imodulon_weights"]
    tcs_hk_set    = integ["tcs_hk"]
    tcs_rr_set    = integ["tcs_rr"]
    sigma_by_locus = integ["sigma_by_locus"]
    print(f"Integration data: {len(imod_by_gene)} iModulon gene entries, "
          f"{len(tcs_hk_set)+len(tcs_rr_set)} TCS loci, {len(sigma_by_locus)} sigma factor loci.")

    positives: Dict[Tuple[str, str], Dict[str, str]] = {}
    for row in regulation_rows:
        positives.setdefault((row["tf_locus"], row["target_locus"]), row)

    tfs = [tf for tf, _target in positives]
    targets = sorted(set(gene_name) | {target for _tf, target in positives} | set(metabolic_features))
    negative_pairs = make_negative_pairs(set(positives), tfs, targets, args.negative_ratio, args.seed)

    def enrich_with_integration(row: Dict[str, Any], tf: str, target: str) -> None:
        """Add iModulon/TCS/sigma features in-place."""
        tf_l = tf.lower()
        tgt_l = target.lower()

        tf_mods  = imod_by_gene.get(tf_l, [])
        tgt_mods = imod_by_gene.get(tgt_l, [])
        shared   = set(tf_mods) & set(tgt_mods)

        # Co-activation score: sum of weight products for shared iModulons
        coact = 0.0
        for im_id in shared:
            im = imod_weights.get(im_id, {})
            genes = im.get("genes", {})
            w_tf  = genes.get(tf_l, 0.0)
            w_tgt = genes.get(tgt_l, 0.0)
            coact += float(w_tf) * float(w_tgt)

        row["target_imodulon_count"] = len(tgt_mods)
        row["tf_imodulon_count"]     = len(tf_mods)
        row["imodulon_coactivation_score"] = format_optional(coact, 4) if coact else ""
        row["tf_is_tcs_regulator"]   = 1 if (tf_l in tcs_hk_set or tf_l in tcs_rr_set) else 0
        row["tf_sigma_class"]        = sigma_by_locus.get(tf_l, "")

    feature_rows: List[Dict[str, Any]] = []
    for (tf, target), row in sorted(positives.items()):
        frow = build_feature_row(
            tf, target, 1, "curated_positive", row, gene_name, gene_to_operon,
            operon_size, srna_features, metabolic_features, expression_correlations
        )
        enrich_with_integration(frow, tf, target)
        feature_rows.append(frow)
    for tf, target in negative_pairs:
        frow = build_feature_row(
            tf, target, 0, "sampled_negative", None, gene_name, gene_to_operon,
            operon_size, srna_features, metabolic_features, expression_correlations
        )
        enrich_with_integration(frow, tf, target)
        feature_rows.append(frow)

    count = write_csv(args.output, feature_rows)
    positives_count = len(positives)
    negatives_count = len(negative_pairs)
    enzyme_count = sum(1 for row in feature_rows if int(row.get("target_has_enzyme_constraint") or 0))
    imodulon_enriched = sum(1 for r in feature_rows if r.get("imodulon_coactivation_score"))
    print(f"Wrote {count} feature rows to {args.output}")
    print(f"Positive edges: {positives_count}")
    print(f"Sampled negatives: {negatives_count}")
    print(f"Rows with enzyme-constrained target features: {enzyme_count}")
    print(f"Rows with iModulon co-activation score: {imodulon_enriched}")
    if expression_correlations:
        print(f"Rows with expression compendium evidence: {sum(1 for r in feature_rows if r.get('expression_feature_available'))}")
    else:
        print("Note: expression correlation columns remain empty unless --expression-correlations is provided or compendium correlations are generated.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
