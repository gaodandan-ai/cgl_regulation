import os
import json
import logging
import sys

try:
    from db_manager import get_db_manager
    _DB_MANAGER_AVAILABLE = True
except ImportError:
    try:
        from backend.db_manager import get_db_manager
        _DB_MANAGER_AVAILABLE = True
    except ImportError:
        get_db_manager = None
        _DB_MANAGER_AVAILABLE = False

try:
    import run_server
except ImportError:
    import backend.run_server as run_server

logger = logging.getLogger("app.services.reference_data")

ESSENTIAL_GENES = {}
PRODORIC_PWMS = {}
BRENDA_KCAT_MAPPINGS = {}
STRING_INTERACTIONS = {}
ABASY_ROLES = {}
RHEA_MAPPINGS = {}
CHEBI_MAPPINGS = {}
COG_ANNOTATIONS = {}


def load_json_reference(file_name: str, parent_dir: str):
    root_dir = parent_dir
    ref_path = os.path.join(root_dir, "data", "reference", file_name)
    if os.path.exists(ref_path):
        with open(ref_path, "r", encoding="utf-8") as f:
            data = json.load(f)
            logger.info(f"Loaded {len(data)} items from {file_name}")
            return data
    else:
        logger.warning(f"{file_name} not found at {ref_path}")
        return {}


def load_all_reference_data(parent_dir: str):
    ESSENTIAL_GENES.clear()
    ESSENTIAL_GENES.update(load_json_reference("essential_genes.json", parent_dir))

    PRODORIC_PWMS.clear()
    PRODORIC_PWMS.update(load_json_reference("prodoric_pwms.json", parent_dir))

    BRENDA_KCAT_MAPPINGS.clear()
    BRENDA_KCAT_MAPPINGS.update(load_json_reference("brenda_kcat_mappings.json", parent_dir))

    STRING_INTERACTIONS.clear()
    STRING_INTERACTIONS.update(load_json_reference("string_interactions.json", parent_dir))

    ABASY_ROLES.clear()
    ABASY_ROLES.update(load_json_reference("abasy_roles.json", parent_dir))

    RHEA_MAPPINGS.clear()
    RHEA_MAPPINGS.update(load_json_reference("rhea_mappings.json", parent_dir))

    CHEBI_MAPPINGS.clear()
    CHEBI_MAPPINGS.update(load_json_reference("chebi_mappings.json", parent_dir))

    COG_ANNOTATIONS.clear()
    COG_ANNOTATIONS.update(load_json_reference("cog_annotations.json", parent_dir))


def check_essentiality(gene_id: str):
    """
    Check if a gene locus tag (e.g. cg0001) or its aliases are classified as essential.
    """
    if not gene_id:
        return None
    if _DB_MANAGER_AVAILABLE:
        db = get_db_manager()
        res = db.get_essential_gene(gene_id)
        if res:
            return res["details"]
        aliases = run_server.expand_gene_aliases(gene_id)
        for alias in aliases:
            res = db.get_essential_gene(alias)
            if res:
                return res["details"]

    # Fallback to memory dict
    g_lower = gene_id.strip().lower()
    if g_lower in ESSENTIAL_GENES:
        return ESSENTIAL_GENES[g_lower]

    aliases = run_server.expand_gene_aliases(g_lower)
    for alias in aliases:
        a_lower = alias.lower()
        if a_lower in ESSENTIAL_GENES:
            return ESSENTIAL_GENES[a_lower]

    return None


def check_abasy_role(gene_id: str):
    """
    Check if a gene locus tag has an Abasy role classification.
    """
    if not gene_id:
        return None
    if _DB_MANAGER_AVAILABLE:
        db = get_db_manager()
        role = db.get_abasy_role(gene_id)
        if role:
            return role
        aliases = run_server.expand_gene_aliases(gene_id)
        for alias in aliases:
            role = db.get_abasy_role(alias)
            if role:
                return role

    # Fallback to memory dict
    g_lower = gene_id.strip().lower()
    if g_lower in ABASY_ROLES:
        return ABASY_ROLES[g_lower]

    aliases = run_server.expand_gene_aliases(g_lower)
    for alias in aliases:
        a_lower = alias.lower()
        if a_lower in ABASY_ROLES:
            return ABASY_ROLES[a_lower]

    return None
