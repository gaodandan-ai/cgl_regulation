"""
backend/metabolic_mapper.py
===========================
Parse and cache gene→reaction mappings from:
  - CSV mapping files in data/reference/metabolic_models/
  - SBML (.xml / .omex) files
  - ecCGL1 JSON enzyme-constrained model files
"""
import os
import csv
import json
import zipfile
import xml.etree.ElementTree as ET

from gene_utils import (
    get_absolute_path,
    expand_gene_aliases,
    split_mapping_values,
    first_row_value,
    infer_model_name_from_file,
    safe_float,
    extract_genes_from_gpr_rule,
    reaction_equation_from_metabolites,
)

METABOLIC_MODEL_DIR   = get_absolute_path(os.path.join("data", "reference", "metabolic_models"))
METABOLIC_MODEL_CACHE = None   # populated on first call to load_metabolic_model_mappings()

# ── ecCGL1 helpers ────────────────────────────────────────────────────────────

def find_ecgl1_root() -> str:
    candidates = [
        get_absolute_path(os.path.join("data", "reference", "model", "ecCGL1-main", "ecCGL1-main")),
        get_absolute_path(os.path.join("data", "reference", "model", "ecCGL1-main")),
        get_absolute_path(os.path.join("data", "reference", "model", "ecCGL1")),
    ]
    for candidate in candidates:
        if os.path.isdir(candidate) and os.path.isdir(os.path.join(candidate, "model")):
            return candidate
    return ""


def load_ecgl1_metadata(ecgl1_root: str):
    gene_uniprot, gene_mass, reaction_metrics, reaction_kcat_sources = {}, {}, {}, {}
    if not ecgl1_root:
        return gene_uniprot, gene_mass, reaction_metrics, reaction_kcat_sources

    # Gene → UniProt mapping
    gene_json = os.path.join(ecgl1_root, "model", "iCW773_uniprot_modification.json")
    if not os.path.exists(gene_json):
        gene_json = os.path.join(ecgl1_root, "model", "iCW773_uniprot.json")
    try:
        with open(gene_json, "r", encoding="utf-8") as f:
            payload = json.load(f)
        for gene in payload.get("genes", []):
            gene_id = str(gene.get("id", "")).strip()
            uniprot = gene.get("annotation", {}).get("uniprot")
            if gene_id and uniprot:
                values = uniprot if isinstance(uniprot, list) else [uniprot]
                gene_uniprot[gene_id.lower()] = [str(v) for v in values if v]
    except Exception:
        pass

    # Gene → protein mass
    mass_json = os.path.join(ecgl1_root, "iCW773_get_data", "iCW773_mean_protein_id_mass_mapping.json")
    try:
        with open(mass_json, "r", encoding="utf-8") as f:
            payload = json.load(f)
        for gene_id, mass in payload.items():
            parsed = safe_float(mass)
            if parsed is not None:
                gene_mass[str(gene_id).lower()] = parsed
    except Exception:
        pass

    # Reaction kcat/MW
    metrics_csv = os.path.join(ecgl1_root, "iCW773_get_data", "reaction_kcat_MW.csv")
    try:
        with open(metrics_csv, "r", encoding="utf-8-sig", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                rid = (first_row_value(row, ["", "reaction", "reaction_id", "id"])
                       or next(iter(row.values()), "")).strip()
                if rid:
                    reaction_metrics[rid] = {
                        "kcat":             safe_float(first_row_value(row, ["kcat"])),
                        "molecular_weight": safe_float(first_row_value(row, ["MW", "mw", "molecular_weight"])),
                        "kcat_MW":          safe_float(first_row_value(row, ["kcat_MW", "kcat_mw"])),
                    }
    except Exception:
        pass

    # Reaction kcat species sources
    kcat_json = os.path.join(ecgl1_root, "iCW773_get_data",
                             "iCW773_mean_reactions_kcat_mapping_combined.json")
    try:
        with open(kcat_json, "r", encoding="utf-8") as f:
            payload = json.load(f)
        for reaction_id, data in payload.items():
            if not isinstance(data, dict):
                continue
            species = []
            for key in ("forward_species_list", "reverse_species_list"):
                for item in data.get(key) or []:
                    if item and item not in species:
                        species.append(str(item))
            reaction_kcat_sources[str(reaction_id)] = {
                "forward":       safe_float(data.get("forward")),
                "reverse":       safe_float(data.get("reverse")),
                "species_count": len(species),
                "species_sample": species[:8],
            }
    except Exception:
        pass

    return gene_uniprot, gene_mass, reaction_metrics, reaction_kcat_sources


def parse_ecgl1_json_mappings(ecgl1_root: str) -> list:
    if not ecgl1_root:
        return []

    gene_uniprot, gene_mass, reaction_metrics, reaction_kcat_sources = load_ecgl1_metadata(ecgl1_root)
    model_dir  = os.path.join(ecgl1_root, "model")
    json_names = [
        "iCW773_irr_enz_constraint.json",
        "iCW773_irr_enz_constraint_adj.json",
        "iCW773_irr_enz_constraint_adj_PDH.json",
        "iCW773_irr_enz_constraint_irreversible.json",
        "iCW773_uniprot_modification_del_irreversible.json",
    ]
    records = []

    for filename in json_names:
        path = os.path.join(model_dir, filename)
        if not os.path.exists(path):
            continue
        model_name   = "ecCGL1:" + os.path.splitext(filename)[0]
        variant_name = os.path.splitext(filename)[0]
        try:
            with open(path, "r", encoding="utf-8") as f:
                payload = json.load(f)
        except Exception:
            continue

        global_constraint = payload.get("enzyme_constraint") if isinstance(payload.get("enzyme_constraint"), dict) else {}

        for reaction in payload.get("reactions", []):
            reaction_id = str(reaction.get("id", "")).strip()
            if not reaction_id:
                continue
            genes = extract_genes_from_gpr_rule(reaction.get("gene_reaction_rule", ""))
            if not genes:
                continue

            metrics     = reaction_metrics.get(reaction_id, {})
            kcat        = safe_float(reaction.get("kcat")) or metrics.get("kcat")
            kcat_mw     = safe_float(reaction.get("kcat_MW")) or metrics.get("kcat_MW")
            mol_weight  = metrics.get("molecular_weight")
            kcat_source = reaction_kcat_sources.get(reaction_id, {})
            annotation  = reaction.get("annotation") if isinstance(reaction.get("annotation"), dict) else {}
            notes       = reaction.get("notes")       if isinstance(reaction.get("notes"), dict)       else {}

            uniprot_ids, protein_masses = [], []
            for gene in genes:
                key = gene.lower()
                for uid in gene_uniprot.get(key, []):
                    if uid not in uniprot_ids:
                        uniprot_ids.append(uid)
                mass = gene_mass.get(key)
                if mass is not None:
                    protein_masses.append(mass)

            if mol_weight is None and protein_masses:
                mol_weight = sum(protein_masses) / len(protein_masses)

            ec = annotation.get("ec-code") or annotation.get("ec_code") or ""

            for gene in genes:
                records.append({
                    "model":              model_name,
                    "gene":               gene,
                    "reaction_id":        reaction_id,
                    "reaction_name":      reaction.get("name") or reaction_id,
                    "equation":           reaction_equation_from_metabolites(reaction.get("metabolites")),
                    "gpr_rule":           reaction.get("gene_reaction_rule", ""),
                    "pathway_id":         "enzyme_constrained_model",
                    "pathway_name":       "Enzyme-constrained model reactions",
                    "source_file":        filename,
                    "ec_number":          ec,
                    "kcat":               kcat,
                    "molecular_weight":   mol_weight,
                    "kcat_MW":            kcat_mw,
                    "uniprot_ids":        uniprot_ids,
                    "protein_masses":     protein_masses,
                    "variant_of":         notes.get("reflection", ""),
                    "reaction_variant":   variant_name,
                    "lower_bound":        reaction.get("lower_bound"),
                    "upper_bound":        reaction.get("upper_bound"),
                    "kcat_source_count":  kcat_source.get("species_count", 0),
                    "kcat_species_sample": kcat_source.get("species_sample", []),
                    "enzyme_constraint": {
                        "model_variant":     variant_name,
                        "kcat":              kcat,
                        "molecular_weight":  mol_weight,
                        "kcat_MW":           kcat_mw,
                        "uniprot_ids":       uniprot_ids,
                        "protein_masses":    protein_masses,
                        "ec_number":         ec,
                        "variant_of":        notes.get("reflection", ""),
                        "lower_bound":       reaction.get("lower_bound"),
                        "upper_bound":       reaction.get("upper_bound"),
                        "kcat_source_count": kcat_source.get("species_count", 0),
                        "kcat_species_sample": kcat_source.get("species_sample", []),
                        "global_parameters": global_constraint,
                    },
                })
    return records

# ── SBML parser ───────────────────────────────────────────────────────────────

def _ns_attr(element, namespace: str, name: str) -> str:
    return element.attrib.get(f"{{{namespace}}}{name}") or element.attrib.get(name, "")


def parse_sbml_gene_reaction_mappings(xml_bytes: bytes, source_name: str):
    core_ns   = "http://www.sbml.org/sbml/level3/version1/core"
    fbc_ns    = "http://www.sbml.org/sbml/level3/version1/fbc/version2"
    groups_ns = "http://www.sbml.org/sbml/level3/version1/groups/version1"

    try:
        root = ET.fromstring(xml_bytes)
    except Exception as e:
        raise ValueError(f"SBML parse failed: {e}")

    model_el = root.find(f".//{{{core_ns}}}model")
    model    = (model_el.attrib.get("id", infer_model_name_from_file(source_name))
                if model_el is not None else infer_model_name_from_file(source_name))

    gene_products = {}
    for gp in root.findall(f".//{{{fbc_ns}}}geneProduct"):
        gid   = _ns_attr(gp, fbc_ns, "id")
        label = _ns_attr(gp, fbc_ns, "label") or gid
        if gid:
            gene_products[gid] = label.replace("G_", "")

    reaction_to_pathway = {}
    for group in root.findall(f".//{{{groups_ns}}}group"):
        pathway_id   = _ns_attr(group, groups_ns, "id")   or group.attrib.get("id", "")
        pathway_name = _ns_attr(group, groups_ns, "name") or pathway_id or "Unassigned pathway"
        for member in group.findall(f".//{{{groups_ns}}}member"):
            rid = _ns_attr(member, groups_ns, "idRef")
            if rid:
                reaction_to_pathway.setdefault(rid, {"pathway_id": pathway_id, "pathway_name": pathway_name})

    records = []
    for reaction in root.findall(f".//{{{core_ns}}}reaction"):
        reaction_id   = reaction.attrib.get("id", "")
        if not reaction_id:
            continue
        reaction_name = reaction.attrib.get("name", reaction_id)
        gene_refs = []
        for ref in reaction.findall(f".//{{{fbc_ns}}}geneProductRef"):
            gp_id = _ns_attr(ref, fbc_ns, "geneProduct")
            gene  = gene_products.get(gp_id, gp_id.replace("G_", "") if gp_id else "")
            if gene:
                gene_refs.append(gene)
        if not gene_refs:
            continue
        pathway = reaction_to_pathway.get(reaction_id, {"pathway_id": "", "pathway_name": "Unassigned pathway"})
        for gene in gene_refs:
            records.append({
                "model":        model,
                "gene":         gene,
                "reaction_id":  reaction_id,
                "reaction_name": reaction_name,
                "equation":     "",
                "gpr_rule":     " ".join(sorted(set(gene_refs))),
                "pathway_id":   pathway["pathway_id"],
                "pathway_name": pathway["pathway_name"],
                "source_file":  os.path.basename(source_name),
            })
    return model, records

# ── Main mapping loader ───────────────────────────────────────────────────────

def load_metabolic_model_mappings() -> dict:
    global METABOLIC_MODEL_CACHE
    if METABOLIC_MODEL_CACHE is not None:
        return METABOLIC_MODEL_CACHE

    from gene_utils import load_gene_mappings
    load_gene_mappings()

    gene_to_reactions, reaction_to_pathways = {}, {}
    files_loaded, warnings = [], []

    def add_mapping_record(record: dict) -> bool:
        model      = record.get("model") or "model"
        gene       = record.get("gene")  or ""
        reaction_id = record.get("reaction_id") or ""
        if not gene or not reaction_id:
            return False
        reaction = {
            "id":           reaction_id,
            "label":        record.get("reaction_name") or reaction_id,
            "model":        model,
            "equation":     record.get("equation", ""),
            "gpr_rule":     record.get("gpr_rule", ""),
            "pathway_id":   record.get("pathway_id", ""),
            "pathway_name": record.get("pathway_name") or "Unassigned pathway",
            "source_file":  record.get("source_file", ""),
        }
        for key in (
            "ec_number", "kcat", "molecular_weight", "kcat_MW", "uniprot_ids",
            "protein_masses", "variant_of", "reaction_variant", "lower_bound",
            "upper_bound", "kcat_source_count", "kcat_species_sample", "enzyme_constraint",
        ):
            if record.get(key) not in (None, "", []):
                reaction[key] = record[key]

        reaction_key = f"{model}:{reaction_id}"
        reaction_to_pathways[reaction_key] = {
            "id":    reaction["pathway_id"] or reaction["pathway_name"],
            "label": reaction["pathway_name"],
            "model": model,
        }
        for alias in expand_gene_aliases(gene):
            gene_to_reactions.setdefault(alias, [])
            if not any(r["model"] == model and r["id"] == reaction_id
                       for r in gene_to_reactions[alias]):
                gene_to_reactions[alias].append(reaction)
        return True

    if not os.path.isdir(METABOLIC_MODEL_DIR):
        METABOLIC_MODEL_CACHE = {
            "loaded": False, "files": [], "models": [],
            "gene_to_reactions": {}, "reaction_to_pathways": {},
            "warnings": [f"Missing mapping directory: {METABOLIC_MODEL_DIR}"],
        }
        return METABOLIC_MODEL_CACHE

    # Load CSV mapping files
    for filename in sorted(os.listdir(METABOLIC_MODEL_DIR)):
        lower = filename.lower()
        if "example" in lower or lower.endswith(".template.csv"):
            continue
        if not (lower.endswith(".csv") and ("reaction" in lower or "gpr" in lower or "mapping" in lower)):
            continue
        path       = os.path.join(METABOLIC_MODEL_DIR, filename)
        file_model = infer_model_name_from_file(path)
        try:
            with open(path, "r", encoding="utf-8-sig", newline="") as f:
                reader = csv.DictReader(f)
                if not reader.fieldnames:
                    warnings.append(f"{filename} has no header row")
                    continue
                row_count = 0
                for row in reader:
                    model      = first_row_value(row, ["model", "model_id", "source_model"]) or file_model
                    gene_field = first_row_value(row, [
                        "gene", "genes", "gene_id", "gene_locus", "locus", "locus_tag",
                        "cg_locus", "cgl_locus", "gene_reaction_rule_genes",
                    ])
                    reaction_id = first_row_value(row, ["reaction_id", "reaction", "rxn_id", "rxn", "id"])
                    if not gene_field or not reaction_id:
                        continue
                    for gene in (split_mapping_values(gene_field) or [gene_field]):
                        if add_mapping_record({
                            "model":        model,
                            "gene":         gene,
                            "reaction_id":  reaction_id,
                            "reaction_name": first_row_value(row, ["reaction_name", "rxn_name", "name", "description"]) or reaction_id,
                            "equation":     first_row_value(row, ["equation", "reaction_equation", "formula"]),
                            "gpr_rule":     first_row_value(row, ["gpr_rule", "gene_reaction_rule", "gpr", "grRule"]),
                            "pathway_id":   first_row_value(row, ["pathway_id", "subsystem_id", "pathway", "subsystem", "module_id"]),
                            "pathway_name": first_row_value(row, ["pathway_name", "subsystem_name", "pathway", "subsystem", "module", "category"])
                                            or first_row_value(row, ["pathway_id", "subsystem_id"]) or "Unassigned pathway",
                            "source_file":  filename,
                        }):
                            row_count += 1
            files_loaded.append({"file": filename, "model": file_model, "rows": row_count})
        except Exception as e:
            warnings.append(f"{filename}: {e}")

    # Load SBML / OMEX model files
    model_dirs = [METABOLIC_MODEL_DIR, get_absolute_path(os.path.join("data", "reference", "model"))]
    for model_dir in model_dirs:
        if not os.path.isdir(model_dir):
            continue
        for filename in sorted(os.listdir(model_dir)):
            if not filename.lower().endswith((".omex", ".xml", ".sbml")):
                continue
            path = os.path.join(model_dir, filename)
            try:
                xml_payloads = []
                if path.lower().endswith(".omex"):
                    with zipfile.ZipFile(path) as archive:
                        for name in archive.namelist():
                            if name.lower().endswith((".xml", ".sbml")) and "manifest" not in name.lower():
                                xml_payloads.append((name, archive.read(name)))
                else:
                    with open(path, "rb") as f:
                        xml_payloads.append((filename, f.read()))

                total_rows, parsed_model = 0, infer_model_name_from_file(path)
                for inner_name, xml_bytes in xml_payloads:
                    parsed_model, records = parse_sbml_gene_reaction_mappings(xml_bytes, inner_name or path)
                    for record in records:
                        record["source_file"] = filename
                        if add_mapping_record(record):
                            total_rows += 1
                if total_rows > 0:
                    files_loaded.append({"file": filename, "model": parsed_model, "rows": total_rows})
            except Exception as e:
                warnings.append(f"{filename}: {e}")

    # Load ecCGL1 JSON mappings
    ecgl1_root = find_ecgl1_root()
    if ecgl1_root:
        try:
            rows_by_model: dict = {}
            for record in parse_ecgl1_json_mappings(ecgl1_root):
                if add_mapping_record(record):
                    mn = record.get("model") or "ecCGL1"
                    rows_by_model[mn] = rows_by_model.get(mn, 0) + 1
            for mn, cnt in sorted(rows_by_model.items()):
                files_loaded.append({"file": "ecCGL1-main", "model": mn, "rows": cnt, "type": "enzyme_constrained"})
        except Exception as e:
            warnings.append(f"ecCGL1-main: {e}")

    models = sorted({f["model"] for f in files_loaded})
    METABOLIC_MODEL_CACHE = {
        "loaded":              len(files_loaded) > 0,
        "files":               files_loaded,
        "models":              models,
        "gene_to_reactions":   gene_to_reactions,
        "reaction_to_pathways": reaction_to_pathways,
        "warnings":            warnings,
    }
    return METABOLIC_MODEL_CACHE
