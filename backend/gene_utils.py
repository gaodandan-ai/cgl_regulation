"""
backend/gene_utils.py
=====================
Gene identifier utilities: normalization, alias expansion, and parsing helpers.

All global mapping dicts (CG_TO_CGL, CGL_TO_CG, etc.) live here so that every
module that imports from gene_utils shares the same in-process objects.
"""
import os
import re
import csv
import math

# ── Root directory resolution ─────────────────────────────────────────────────
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR  = os.path.dirname(_THIS_DIR)   # project root (one level up from backend/)

def get_absolute_path(relative_path: str) -> str:
    if not relative_path:
        return relative_path
    if os.path.isabs(relative_path):
        return relative_path
    return os.path.normpath(os.path.join(ROOT_DIR, relative_path))

# ── Shared gene mapping globals ───────────────────────────────────────────────
CG_TO_CGL:  dict = {}   # cg_locus (lower) -> cgl_locus
CGL_TO_CG:  dict = {}   # cgl_locus (lower) -> cg_locus
GENE_NAMES: dict = {}   # any locus (lower) -> display name
NAME_TO_CG: dict = {}   # gene name (lower) -> cg_locus
GENE_TO_UNIPROT: dict = {} # locus (lower) -> uniprot accession id

_gene_mappings_loaded = False

def load_gene_mappings() -> None:
    """Load gene_mapping.csv into the global dicts (idempotent)."""
    global _gene_mappings_loaded
    if _gene_mappings_loaded:
        return
    _gene_mappings_loaded = True

    path = get_absolute_path("data/reference/gene_mapping.csv")
    if not os.path.exists(path):
        return
    try:
        with open(path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                cg   = row.get("cg_locus",  "").strip()
                cgl  = row.get("cgl_locus", "").strip()
                name = row.get("gene_name", "").strip()
                uniprot = row.get("uniprot_id", "").strip()
                if cg and cgl:
                    CG_TO_CGL[cg.lower()]  = cgl
                    CGL_TO_CG[cgl.lower()] = cg
                if cg and name:
                    GENE_NAMES[cg.lower()] = name
                    NAME_TO_CG.setdefault(name.lower(), cg)
                if cgl and name:
                    GENE_NAMES[cgl.lower()] = name
                if uniprot:
                    if cg:
                        GENE_TO_UNIPROT[cg.lower()] = uniprot
                    if cgl:
                        GENE_TO_UNIPROT[cgl.lower()] = uniprot
    except Exception as e:
        print("Error loading gene mapping CSV:", e)

# ── Normalization & alias expansion ──────────────────────────────────────────

def normalize_gene_locus(locus: str) -> str:
    locus = (locus or "").strip()
    if not locus:
        return ""
    lower = locus.lower()
    if lower in CGL_TO_CG:
        return CGL_TO_CG[lower].lower()
    if lower in CG_TO_CGL:
        return lower
    if lower in NAME_TO_CG:
        return NAME_TO_CG[lower].lower()
    return lower


def expand_gene_aliases(locus: str) -> set:
    """Return a set of all known aliases for a gene locus tag."""
    load_gene_mappings()
    aliases: set = set()
    lower = (locus or "").strip().lower()
    if not lower:
        return aliases
    aliases.add(lower)

    # Support NCBI lcl prot gene identifiers like "lcl_123"
    if "_" in lower:
        parts = lower.split("_")
        if parts[-1].isdigit():
            num = int(parts[-1])
            aliases.add(f"cg{num:04d}")

    canonical = normalize_gene_locus(lower)
    if canonical:
        aliases.add(canonical.lower())

    for alias in list(aliases):
        if alias in CG_TO_CGL:
            aliases.add(CG_TO_CGL[alias].lower())
        if alias in CGL_TO_CG:
            aliases.add(CGL_TO_CG[alias].lower())

    return aliases

# ── Parsing helpers ───────────────────────────────────────────────────────────

def split_mapping_values(value: str) -> list:
    text = (value or "").strip()
    if not text:
        return []
    parts = re.split(r"[;,|]+|\s+and\s+|\s+or\s+", text, flags=re.IGNORECASE)
    return [p.strip() for p in parts if p.strip()]


def first_row_value(row: dict, names: list) -> str:
    for name in names:
        if name in row and row.get(name):
            return str(row.get(name, "")).strip()
    lower_map = {k.lower(): v for k, v in row.items()}
    for name in names:
        val = lower_map.get(name.lower())
        if val:
            return str(val).strip()
    return ""


def infer_model_name_from_file(path: str) -> str:
    name = os.path.splitext(os.path.basename(path))[0]
    for suffix in ("_gene_reaction_mapping", "_gene_reaction_map", "_reaction_mapping"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return name


def safe_float(value) -> float | None:
    try:
        if value is None or value == "":
            return None
        parsed = float(value)
        return parsed if math.isfinite(parsed) else None
    except Exception:
        return None


def extract_genes_from_gpr_rule(rule: str) -> list:
    text = str(rule or "")
    if not text:
        return []
    genes = re.findall(r"\b(?:Cgl|cg|cgb|NCgl|cgR)[:_ -]?\d{3,5}\b", text, flags=re.IGNORECASE)
    cleaned = []
    for gene in genes:
        gene = re.sub(r"[:_ -]+", "", gene).strip()
        if gene and gene not in cleaned:
            cleaned.append(gene)
    if cleaned:
        return cleaned
    tokens = re.split(r"\s+and\s+|\s+or\s+|[(),;|]+", text, flags=re.IGNORECASE)
    return [t.strip() for t in tokens if t.strip() and t.strip().lower() not in {"and", "or"}]


def reaction_equation_from_metabolites(metabolites: dict) -> str:
    if not isinstance(metabolites, dict) or not metabolites:
        return ""
    reactants, products = [], []
    for metabolite, coefficient in metabolites.items():
        value = safe_float(coefficient)
        if value is None or value == 0:
            continue
        label  = str(metabolite)
        amount = abs(value)
        term   = label if amount == 1 else f"{amount:g} {label}"
        if value < 0:
            reactants.append(term)
        else:
            products.append(term)
    if not reactants and not products:
        return ""
    return f"{' + '.join(reactants) or '0'} -> {' + '.join(products) or '0'}"


def evaluate_gpr_rule(rule_str: str, gene_values: dict, default_val: float = 1.0) -> float:
    """
    Evaluates a COBRApy gene_reaction_rule string under min/max semantics:
      - 'and' -> min
      - 'or' -> max
    """
    if not rule_str or not rule_str.strip():
        return default_val

    # Tokenize: find parenthesis, words, operators
    tokens = re.findall(r"\(|\)|\w+", rule_str)
    index = 0

    def parse_expression() -> float:
        nonlocal index
        val = parse_term()
        while index < len(tokens) and tokens[index].lower() == 'or':
            index += 1
            right = parse_term()
            val = max(val, right)
        return val

    def parse_term() -> float:
        nonlocal index
        val = parse_factor()
        while index < len(tokens) and tokens[index].lower() == 'and':
            index += 1
            right = parse_factor()
            val = min(val, right)
        return val

    def parse_factor() -> float:
        nonlocal index
        if index >= len(tokens):
            return default_val

        token = tokens[index]
        if token == '(':
            index += 1
            val = parse_expression()
            if index < len(tokens) and tokens[index] == ')':
                index += 1
            return val
        else:
            index += 1
            norm_token = normalize_gene_locus(token)
            
            val = None
            if token in gene_values:
                val = gene_values[token]
            elif norm_token in gene_values:
                val = gene_values[norm_token]
            else:
                for k, v in gene_values.items():
                    if normalize_gene_locus(k) == norm_token:
                        val = v
                        break
            if val is None:
                return default_val
            return val

    try:
        return parse_expression()
    except Exception:
        return default_val
