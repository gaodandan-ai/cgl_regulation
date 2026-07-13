"""scripts/add_thermo_gene_api.py — injects /api/thermo/gene_context endpoint into app.py"""

ENDPOINT_CODE = r'''
# ── Thermodynamic Gene Context API ────────────────────────────────────────────
_THERMO_DATA_CACHE = None

def _load_thermo_gene_data():
    """Load thermo data, cached after first call."""
    global _THERMO_DATA_CACHE
    if _THERMO_DATA_CACHE is not None:
        return _THERMO_DATA_CACHE
    path = os.path.join(os.path.dirname(BACKEND_DIR), "data", "reference", "thermo_dgr_data.json")
    if not os.path.exists(path):
        return {}
    with open(path, "r", encoding="utf-8") as f:
        _THERMO_DATA_CACHE = json.load(f)
    return _THERMO_DATA_CACHE


@app.get("/api/thermo/gene_context")
def get_thermo_gene_context(gene: str = ""):
    """
    Return thermodynamic status for all reactions associated with a given gene.
    Used to annotate ecFBA results with thermodynamic context.
    
    Returns:
      - reactions: list of {id, name, direction_locked, dgr0, dgr_min, dgr_max, confidence}
      - thermo_support_level: 'strong' | 'moderate' | 'weak' | 'none'
      - n_locked: number of this gene's reactions that are direction-locked
      - ko_thermo_confidence: 0.0-1.0 score
    """
    if not gene:
        raise HTTPException(status_code=400, detail="Missing gene parameter")

    thermo_data = _load_thermo_gene_data()
    thermo_rxns = thermo_data.get("reactions", {})

    # Get the loaded model to find gene→reaction associations
    from model_loader import load_model_if_needed
    model = load_model_if_needed()
    if model is None:
        raise HTTPException(status_code=503, detail="Model not loaded")

    # Expand aliases
    gene_lower = gene.strip().lower()
    aliases = set()
    try:
        aliases = run_server.expand_gene_aliases(gene_lower)
    except Exception:
        aliases = {gene_lower, gene.strip()}

    # Find reactions associated with this gene
    gene_reactions = []
    for rxn in model.reactions:
        rxn_gene_ids = {g.id.lower() for g in rxn.genes}
        if aliases & rxn_gene_ids:
            gene_reactions.append(rxn.id)

    if not gene_reactions:
        return {
            "gene": gene,
            "gene_reactions": [],
            "thermo_annotated": [],
            "thermo_support_level": "none",
            "n_locked": 0,
            "ko_thermo_confidence": 0.0,
            "message": "No reactions found for this gene in the metabolic model"
        }

    # Annotate each reaction with thermo data
    annotated = []
    n_locked = 0
    n_confirmed = 0
    confidence_sum = 0.0
    conf_weights = {"HIGH": 1.0, "MED": 0.6, "LOW": 0.3}

    for rxn_id in gene_reactions:
        rxn = model.reactions.get_by_id(rxn_id)
        entry = thermo_rxns.get(rxn_id, {})

        direction_locked = entry.get("direction_locked", "none")
        dgr0 = entry.get("dgr_prime_0")
        dgr_min = entry.get("dgr_prime_min")
        dgr_max = entry.get("dgr_prime_max")
        conf = entry.get("confidence", "NONE")
        conf_weight = conf_weights.get(conf, 0.0)

        is_locked = direction_locked in ("forward", "reverse")
        if is_locked:
            n_locked += 1
            confidence_sum += conf_weight
        elif dgr0 is not None:
            n_confirmed += 1
            confidence_sum += conf_weight * 0.5

        annotated.append({
            "reaction_id": rxn_id,
            "reaction_name": rxn.name or rxn_id,
            "direction_locked": direction_locked,
            "current_lb": rxn.lower_bound,
            "current_ub": rxn.upper_bound,
            "has_thermo_data": dgr0 is not None,
            "dgr_prime_0": dgr0,
            "dgr_prime_min": dgr_min,
            "dgr_prime_max": dgr_max,
            "confidence": conf,
            "thermo_label": (
                f"ΔG'°={dgr0:.1f} kJ/mol" if dgr0 is not None else "No ΔG' data"
            ),
        })

    # Compute overall thermo support level
    n_rxns = len(gene_reactions)
    lock_fraction = n_locked / n_rxns if n_rxns > 0 else 0
    data_fraction = (n_locked + n_confirmed) / n_rxns if n_rxns > 0 else 0
    ko_confidence = min(confidence_sum / max(n_rxns, 1), 1.0)

    if lock_fraction >= 0.5:
        support_level = "strong"
    elif lock_fraction > 0 or data_fraction >= 0.5:
        support_level = "moderate"
    elif data_fraction > 0:
        support_level = "weak"
    else:
        support_level = "none"

    return {
        "gene": gene,
        "gene_reactions": gene_reactions,
        "thermo_annotated": sorted(annotated, key=lambda x: x["direction_locked"] != "none", reverse=True),
        "thermo_support_level": support_level,
        "n_locked": n_locked,
        "n_confirmed": n_confirmed,
        "n_no_data": n_rxns - n_locked - n_confirmed,
        "total_reactions": n_rxns,
        "lock_fraction": round(lock_fraction, 3),
        "ko_thermo_confidence": round(ko_confidence, 3),
    }

'''

MARKER = '@app.get("/api/metabolic_impact")\ndef metabolic_impact'
REPLACEMENT = ENDPOINT_CODE + '@app.get("/api/metabolic_impact")\ndef metabolic_impact'

with open("backend/app.py", "r", encoding="utf-8") as f:
    content = f.read()

if MARKER in content:
    content = content.replace(MARKER, REPLACEMENT, 1)
    with open("backend/app.py", "w", encoding="utf-8") as f:
        f.write(content)
    print("SUCCESS: /api/thermo/gene_context endpoint added")
else:
    print("ERROR: marker not found")
