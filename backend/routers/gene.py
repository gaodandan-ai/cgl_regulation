from fastapi import APIRouter, HTTPException, Response
import os
import sys
import logging

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

from services.reference_data import (
    ESSENTIAL_GENES,
    BRENDA_KCAT_MAPPINGS,
    ABASY_ROLES,
    COG_ANNOTATIONS,
    STRING_INTERACTIONS
)

try:
    from schemas import (
        GeneCoordinatesResponse,
        GenomicNeighborhoodResponse,
        ChIPSeqPeaksResponse,
        HTTPError
    )
except ImportError:
    from backend.schemas import (
        GeneCoordinatesResponse,
        GenomicNeighborhoodResponse,
        ChIPSeqPeaksResponse,
        HTTPError
    )

router = APIRouter(tags=["Gene & Quality"])
logger = logging.getLogger("app.routers.gene")


@router.get(
    "/api/gene/coordinates/{gene_id}",
    response_model=GeneCoordinatesResponse,
    responses={404: {"model": HTTPError}}
)
def get_gene_coordinates_api(gene_id: str):
    db = get_db_manager()
    res = db.get_gene_coordinates(gene_id)
    if not res:
        raise HTTPException(status_code=404, detail=f"Coordinates not found for {gene_id}")
    return res


@router.get(
    "/api/gene/profile/{gene_id}",
    responses={404: {"model": HTTPError}}
)
def get_full_gene_profile_api(gene_id: str):
    db = get_db_manager()
    res = db.get_full_gene_profile(gene_id)
    if not res:
        raise HTTPException(status_code=404, detail=f"Profile not found for {gene_id}")
    return res


@router.get(
    "/api/gene/neighborhood/{gene_id}",
    response_model=GenomicNeighborhoodResponse
)
def get_genomic_neighborhood_api(gene_id: str, window_bp: int = 20000):
    db = get_db_manager()
    genes = db.get_genomic_neighborhood(gene_id, window_bp=window_bp)
    return {"center_gene": gene_id, "window_bp": window_bp, "count": len(genes), "genes": genes}


@router.get("/api/genomic_tracks/{gene_id}")
def get_genomic_tracks_api(gene_id: str, window_bp: int = 10000):
    db = get_db_manager()
    data = db.get_genomic_track_data(gene_id, window_bp=window_bp)
    if not data:
        raise HTTPException(status_code=404, detail=f"Genomic track data unavailable for {gene_id}")
    return data


try:
    from security import PUBLIC_DEPLOYMENT
except ImportError:
    try:
        from backend.security import PUBLIC_DEPLOYMENT
    except ImportError:
        PUBLIC_DEPLOYMENT = False


@router.get("/api/chipseq_peaks/{gene_id}", response_model=ChIPSeqPeaksResponse)
def get_chipseq_peaks_api(gene_id: str):
    """Retrieve all experimental ChIP-seq binding peaks associated with target gene or TF."""
    if PUBLIC_DEPLOYMENT:
        return {
            "query": gene_id,
            "as_target_count": 0,
            "as_target_peaks": [],
            "as_tf_count": 0,
            "as_tf_peaks": [],
            "is_public_deployment": True,
            "message": "Full experimental ChIP-seq peak datasets are available exclusively on the laboratory intranet server (172.16.2.105:8010)."
        }

    db = get_db_manager()
    gene_peaks = db.get_gene_chipseq_peaks(gene_id)
    tf_peaks = db.get_tf_chipseq_peaks(gene_id)
    return {
        "query": gene_id,
        "as_target_count": len(gene_peaks),
        "as_target_peaks": gene_peaks,
        "as_tf_count": len(tf_peaks),
        "as_tf_peaks": tf_peaks,
        "is_public_deployment": False
    }


@router.get("/api/quality/essential")
def get_essential_genes(response: Response):
    """Return the database of C. glutamicum essential genes."""
    response.headers["Cache-Control"] = "public, max-age=3600"
    return ESSENTIAL_GENES


@router.get("/api/quality/brenda")
def get_brenda_mappings(response: Response):
    """Return the database of C. glutamicum BRENDA kcat mappings."""
    response.headers["Cache-Control"] = "public, max-age=3600"
    return BRENDA_KCAT_MAPPINGS


@router.get("/api/quality/abasy")
def get_abasy_roles(response: Response):
    """Return the database of C. glutamicum Abasy roles."""
    response.headers["Cache-Control"] = "public, max-age=3600"
    return ABASY_ROLES


@router.get("/api/quality/cog")
def get_cog_annotations(response: Response):
    """Return the database of C. glutamicum COG functional annotations."""
    response.headers["Cache-Control"] = "public, max-age=3600"
    return COG_ANNOTATIONS


@router.get("/api/quality/ppi")
def quality_ppi():
    try:
        genes = [k for k in STRING_INTERACTIONS.keys() if k != "_meta"]
        total_proteins = len(genes)

        edges = set()
        scores = []
        very_high = 0
        high = 0
        medium = 0
        low = 0

        channels = {
            "experimental": 0,
            "database": 0,
            "coexpression": 0,
            "textmining": 0,
            "neighborhood": 0,
            "cooccurrence": 0,
            "fusion": 0
        }

        for g in genes:
            for p in STRING_INTERACTIONS[g]:
                partner = p.get("partner")
                if not partner:
                    continue
                edge_key = tuple(sorted([g, partner.lower()]))
                if edge_key not in edges:
                    edges.add(edge_key)
                    score = p.get("score", 0)
                    scores.append(score)

                    if score >= 900:
                        very_high += 1
                    elif score >= 700:
                        high += 1
                    elif score >= 400:
                        medium += 1
                    else:
                        low += 1

                    for ch in channels.keys():
                        if p.get(ch, 0) > 0:
                            channels[ch] += 1

        total_edges = len(edges)
        avg_score = sum(scores) / len(scores) if scores else 0

        return {
            "total_proteins": total_proteins,
            "total_interactions": total_edges,
            "avg_partners": round(total_edges * 2 / total_proteins, 1) if total_proteins > 0 else 0,
            "avg_score": round(avg_score, 1),
            "score_distribution": {
                "very_high": very_high,
                "high": high,
                "medium": medium,
                "low": low
            },
            "channel_support": channels
        }
    except Exception as e:
        logger.error(f"Error computing PPI quality: {e}")
        return {"error": str(e)}


@router.get("/api/quality/icgb21fr")
def quality_icgb21fr():
    """Compute regulatory gene coverage statistics against the iCGB21FR model."""
    import cobra, csv as csv_mod
    try:
        backend_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        root_dir = os.path.dirname(backend_dir)
        model_path = os.path.join(root_dir, "data", "reference", "model", "iCGB21FR.xml")
        if not os.path.exists(model_path):
            raise HTTPException(status_code=404, detail=f"iCGB21FR.xml not found at {model_path}")

        model = cobra.io.read_sbml_model(model_path)

        # Build gene-to-reactions and gene-to-subsystems maps
        gene_to_rxns: dict = {}
        gene_to_paths: dict = {}
        for rxn in model.reactions:
            subsystem = (rxn.subsystem or "").strip()
            for gene in rxn.genes:
                g_id = gene.id.strip().lower()
                all_ids = run_server.expand_gene_aliases(g_id)
                all_ids.add(g_id)
                for aid in all_ids:
                    gene_to_rxns.setdefault(aid.lower(), set()).add(rxn.id)
                    if subsystem:
                        gene_to_paths.setdefault(aid.lower(), set()).add(subsystem)

        # Load regulatory genes from regulations.csv
        reg_path = os.path.join(root_dir, "data", "reference", "regulations.csv")
        reg_genes_raw: set = set()
        with open(reg_path, "r", encoding="utf-8") as csvf:
            reader = csv_mod.DictReader(csvf)
            for row in reader:
                for field in ("TF_locusTag", "TF_altLocusTag", "TG_locusTag", "TG_altLocusTag", "TF", "Target", "tf_locus", "target_locus"):
                    val = row.get(field, "").strip().lower()
                    if val:
                        reg_genes_raw.add(val)

        # Expand aliases for all regulatory genes
        unique_reg_genes: set = set()
        for rg in reg_genes_raw:
            unique_reg_genes.add(rg)
            for alias in run_server.expand_gene_aliases(rg):
                unique_reg_genes.add(alias.lower())

        # Compute coverage
        mapped_rxn_genes = 0
        mapped_path_genes = 0
        unique_rxns: set = set()
        unique_paths: set = set()
        unmapped = []

        for gene in unique_reg_genes:
            rxns = gene_to_rxns.get(gene, set())
            paths = gene_to_paths.get(gene, set())
            if rxns:
                mapped_rxn_genes += 1
                unique_rxns.update(rxns)
            if paths:
                mapped_path_genes += 1
                unique_paths.update(paths)
            if not rxns:
                unmapped.append(gene)

        unmapped = sorted(set(unmapped))

        return {
            "model_id": model.id,
            "model_genes": len(model.genes),
            "regulatory_gene_count": len(unique_reg_genes),
            "genes_mapped_to_reactions": mapped_rxn_genes,
            "genes_mapped_to_pathways": mapped_path_genes,
            "unique_mapped_reactions": len(unique_rxns),
            "unique_mapped_pathways": len(unique_paths),
            "unmapped_gene_count": len(unmapped),
            "unmapped_genes": unmapped[:100]
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"iCGB21FR quality endpoint failed: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/homolog_alignment")
def homolog_alignment(gene_name: str = "", accession: str = ""):
    try:
        result = run_server.handle_homolog_alignment(gene_name, accession)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/list_organisms")
def list_organisms():
    try:
        organisms = []
        backend_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        root_dir = os.path.dirname(backend_dir)
        folder = os.path.join(root_dir, 'data', 'reference', 'AllOrganismsFiles')
        if os.path.exists(folder):
            for filename in os.listdir(folder):
                if filename.endswith('_regulations.csv'):
                    org_id = filename[:-16]
                    if not org_id:
                        continue
                    name = org_id
                    parts = org_id.split('_', 2)
                    if len(parts) >= 2:
                        key = f"{parts[0]}_{parts[1]}"
                        rest = parts[2] if len(parts) > 2 else ""
                        if key in run_server.SPECIES_MAP:
                            clean_rest = rest.replace('_', ' ').strip()
                            name = f"{run_server.SPECIES_MAP[key]} {clean_rest}".strip()
                        else:
                            name = org_id.replace('_', ' ')
                    else:
                        name = org_id.replace('_', ' ')
                    rna_file = f"{org_id}_rna_regulation.csv"
                    has_rna = os.path.exists(os.path.join(folder, rna_file))
                    organisms.append({
                        "id": org_id,
                        "name": name,
                        "has_rna": has_rna
                    })
        return {"organisms": organisms}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/ncrna/list")
def get_ncrna_list_api(rna_type: str = None):
    db = get_db_manager()
    rnas = db.get_ncrnas(rna_type=rna_type)
    return {"rna_type": rna_type, "count": len(rnas), "ncrnas": rnas}


@router.get("/api/ncrna/targets")
def get_srna_targets_api(locus: str = None):
    db = get_db_manager()
    targets = db.get_srna_targets(locus=locus)
    return {"locus": locus, "count": len(targets), "targets": targets}
