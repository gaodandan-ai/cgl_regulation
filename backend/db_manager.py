#!/usr/bin/env python3
"""
db_manager.py
=============
Thread-safe database manager for querying cgl_regulation.db with connection pooling
and prepared SQL statements. Provides fast, indexed lookups for all backend handlers.
"""

import os
import sys
import json
import sqlite3
import threading
import logging

logger = logging.getLogger("cgl_db")

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_PUBLIC_DB_PATH = os.path.join(ROOT_DIR, "data", "deploy", "cgl_regulation_public.db")
DB_PATH = os.environ.get("CGL_DATABASE_PATH") or (
    _PUBLIC_DB_PATH
    if (os.environ.get("VERCEL") == "1" and os.path.isfile(_PUBLIC_DB_PATH))
    else os.path.join(ROOT_DIR, "data", "reference", "cgl_regulation.db")
)

class CglDatabaseManager:
    """
    Singleton Manager for querying the SQLite cgl_regulation.db.
    Uses thread-local storage for safe concurrent access in FastAPI/Uvicorn workers.
    """
    _instance = None
    _lock = threading.Lock()

    def __new__(cls, db_path=None):
        with cls._lock:
            if cls._instance is None:
                cls._instance = super(CglDatabaseManager, cls).__new__(cls)
                cls._instance._db_path = db_path or DB_PATH
                cls._instance._local = threading.local()
                cls._instance._initialized = True
            return cls._instance

    def get_connection(self):
        if not hasattr(self._local, "conn") or self._local.conn is None:
            if not os.path.exists(self._db_path):
                logger.warning(f"Database file not found at {self._db_path}. Fallbacks will be used.")
                return None
            db_uri = f"file:{os.path.abspath(self._db_path).replace(os.sep, '/')}?mode=ro"
            self._local.conn = sqlite3.connect(db_uri, uri=True, check_same_thread=False, timeout=5)
            self._local.conn.row_factory = sqlite3.Row
            self._local.conn.execute("PRAGMA foreign_keys = ON")
            self._local.conn.execute("PRAGMA query_only = ON")
            self._local.conn.execute("PRAGMA busy_timeout = 5000")
        return self._local.conn

    def close(self):
        if hasattr(self._local, "conn") and self._local.conn is not None:
            self._local.conn.close()
            self._local.conn = None

    def get_canonical_locus(self, alias: str):
        if not alias:
            return None
        conn = self.get_connection()
        if not conn:
            return None
        cursor = conn.cursor()
        cursor.execute(
            "SELECT alias, canonical_cg, canonical_cgl, gene_name, product FROM canonical_locus_map WHERE alias = ? OR alias = ?",
            (alias, alias.lower())
        )
        row = cursor.fetchone()
        return dict(row) if row else None

    def get_essential_gene(self, gene_id: str):
        if not gene_id:
            return None
        conn = self.get_connection()
        if not conn:
            return None
        cursor = conn.cursor()
        cursor.execute(
            "SELECT locus_tag, symbol, essentiality, details FROM essential_genes WHERE locus_tag = ?",
            (gene_id.lower(),)
        )
        row = cursor.fetchone()
        if row:
            res = dict(row)
            try:
                res["details"] = json.loads(res["details"]) if res["details"] else {}
            except Exception:
                pass
            return res
        return None

    def get_abasy_role(self, gene_id: str):
        if not gene_id:
            return None
        conn = self.get_connection()
        if not conn:
            return None
        cursor = conn.cursor()
        cursor.execute(
            "SELECT systemic_role, details FROM abasy_roles WHERE locus_tag = ?",
            (gene_id.lower(),)
        )
        row = cursor.fetchone()
        if row:
            return row["systemic_role"]
        return None

    def get_all_abasy_roles(self):
        conn = self.get_connection()
        if not conn:
            return {}
        cursor = conn.cursor()
        cursor.execute("SELECT locus_tag, systemic_role FROM abasy_roles")
        return {row["locus_tag"]: row["systemic_role"] for row in cursor.fetchall()}

    def get_string_interactions(self, gene_id: str, min_score: float = 400.0):
        if not gene_id:
            return []
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        g_lower = gene_id.lower()
        cursor.execute(
            "SELECT gene_b as partner, score FROM string_interactions WHERE gene_a = ? AND score >= ? ORDER BY score DESC",
            (g_lower, min_score)
        )
        return [dict(r) for r in cursor.fetchall()]

    def get_regulations_for_tf(self, tf_id: str):
        if not tf_id:
            return []
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        cursor.execute(
            "SELECT * FROM regulations WHERE TF_locusTag = ? OR TF_name = ? OR LOWER(TF_locusTag) = ? OR LOWER(TF_name) = ?",
            (tf_id, tf_id, tf_id.lower(), tf_id.lower())
        )
        return [dict(r) for r in cursor.fetchall()]

    def get_regulations_for_tg(self, tg_id: str):
        if not tg_id:
            return []
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        cursor.execute(
            "SELECT * FROM regulations WHERE TG_locusTag = ? OR TG_name = ? OR LOWER(TG_locusTag) = ? OR LOWER(TG_name) = ?",
            (tg_id, tg_id, tg_id.lower(), tg_id.lower())
        )
        return [dict(r) for r in cursor.fetchall()]

    def get_brenda_kcat(self, rxn_id: str):
        if not rxn_id:
            return None
        conn = self.get_connection()
        if not conn:
            return None
        cursor = conn.cursor()
        cursor.execute(
            "SELECT reaction_id, ec_number, kcat_val, substrate, details FROM brenda_kcat WHERE reaction_id = ?",
            (rxn_id,)
        )
        row = cursor.fetchone()
        if row:
            res = dict(row)
            try:
                res["details"] = json.loads(res["details"]) if res["details"] else {}
            except Exception:
                pass
            return res
        return None

    def get_rhea_mapping(self, rxn_id: str):
        if not rxn_id:
            return None
        conn = self.get_connection()
        if not conn:
            return None
        cursor = conn.cursor()
        cursor.execute("SELECT links FROM rhea_mappings WHERE rxn_id = ?", (rxn_id,))
        row = cursor.fetchone()
        if row and row["links"]:
            try:
                return json.loads(row["links"])
            except Exception:
                pass
        return None

    def get_chebi_mapping(self, met_id: str):
        if not met_id:
            return None
        conn = self.get_connection()
        if not conn:
            return None
        cursor = conn.cursor()
        cursor.execute("SELECT details FROM chebi_mappings WHERE met_id = ?", (met_id,))
        row = cursor.fetchone()
        if row and row["details"]:
            try:
                return json.loads(row["details"])
            except Exception:
                pass
        return None

    def get_cog_annotation(self, locus: str):
        if not locus:
            return None
        conn = self.get_connection()
        if not conn:
            return None
        cursor = conn.cursor()
        cursor.execute(
            "SELECT cog_id, category, description FROM cog_annotations WHERE locus_tag = ?",
            (locus.lower(),)
        )
        row = cursor.fetchone()
        return dict(row) if row else None

    def get_network_centrality(self, locus: str):
        if not locus:
            return None
        conn = self.get_connection()
        if not conn:
            return None
        cursor = conn.cursor()
        cursor.execute(
            "SELECT degree, in_degree, out_degree, betweenness, closeness, pagerank FROM network_centrality WHERE locus_tag = ?",
            (locus.lower(),)
        )
        row = cursor.fetchone()
        return dict(row) if row else None

    def search_literature_fts(self, query: str, limit: int = 20):
        if not query:
            return []
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        # Clean query for FTS5
        clean_query = query.replace('"', '').replace("'", "").strip()
        if not clean_query:
            return []
        try:
            cursor.execute(
                """
                SELECT gene_locus, title, abstract, pmid, rank
                FROM literature_fts
                WHERE literature_fts MATCH ?
                ORDER BY rank
                LIMIT ?
                """,
                (f'"{clean_query}" OR {clean_query}*', limit)
            )
            return [dict(r) for r in cursor.fetchall()]
        except Exception as e:
            logger.warning(f"FTS search error: {e}")
            return []

    def get_extended_edges(self, locus_tag: str, mode: str = "all", edge_type: str = None):
        """
        Query network_edges_extended by source or target locus tag.
        mode: 'strong', 'all'
        edge_type: 'tf_dna', 'srna_mrna', or None for both
        """
        if not locus_tag:
            return []
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        g_lower = locus_tag.lower().strip()

        query = "SELECT * FROM network_edges_extended WHERE (source_locus = ? OR target_locus = ?)"
        params = [g_lower, g_lower]

        if mode == "strong":
            query += " AND evidence_level = 'strong'"
        if edge_type:
            query += " AND edge_type = ?"
            params.append(edge_type)

        query += " ORDER BY CASE WHEN edge_type = 'srna_mrna' THEN binding_energy_kcal END ASC, CASE WHEN edge_type = 'tf_dna' THEN confidence_score END DESC LIMIT 100"
        cursor.execute(query, params)
        rows = cursor.fetchall()
        res = []
        for r in rows:
            item = dict(r)
            if item.get("details"):
                try:
                    item["details"] = json.loads(item["details"])
                except Exception:
                    pass
            res.append(item)
        return res

    def get_gene_coordinates(self, locus_tag: str):
        """
        Query NCBI RefSeq genomic coordinates for a given gene locus tag.
        """
        if not locus_tag:
            return None
        conn = self.get_connection()
        if not conn:
            return None
        cursor = conn.cursor()
        cursor.execute(
            "SELECT locus_tag, gene_name, start_pos, end_pos, strand, gene_length, tss_position, promoter_70bp FROM gene_coordinates WHERE locus_tag = ?",
            (locus_tag.lower().strip(),)
        )
        row = cursor.fetchone()
        return dict(row) if row else None

    def get_tf_effector_info(self, tf_id: str):
        """
        Query UniProt TF structural family and small-molecule effector annotations.
        """
        if not tf_id:
            return None
        conn = self.get_connection()
        if not conn:
            return None
        cursor = conn.cursor()
        g_lower = tf_id.lower().strip()
        cursor.execute(
            "SELECT tf_locus, tf_name, tf_family, hth_domain, effector_molecule, physiological_signal, regulatory_role FROM tf_families_effectors WHERE tf_locus = ? OR LOWER(tf_name) = ?",
            (g_lower, g_lower)
        )
        row = cursor.fetchone()
        return dict(row) if row else None

    def get_full_gene_profile(self, locus_tag: str):
        """
        Query the v_gene_full_profile SQL View joining gene mappings, RefSeq coordinates,
        TF families & effectors, Abasy roles, and Network Centrality.
        """
        if not locus_tag:
            return None
        conn = self.get_connection()
        if not conn:
            return None
        cursor = conn.cursor()
        g_lower = locus_tag.lower().strip()
        cursor.execute(
            "SELECT * FROM v_gene_full_profile WHERE LOWER(cg_locus) = ? OR LOWER(cgl_locus) = ? OR LOWER(gene_name) = ?",
            (g_lower, g_lower, g_lower)
        )
        row = cursor.fetchone()
        return dict(row) if row else None

    def get_genomic_neighborhood(self, locus_tag: str, window_bp: int = 20000):
        """
        Spatial genomic query fetching all genes within +/- window_bp of locus_tag.
        """
        target_coords = self.get_gene_coordinates(locus_tag)
        if not target_coords or not target_coords.get("start_pos"):
            return []
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        center = target_coords["start_pos"]
        min_pos = max(1, center - window_bp)
        max_pos = center + window_bp

        cursor.execute(
            """
            SELECT locus_tag, gene_name, start_pos, end_pos, strand, gene_length, tss_position
            FROM gene_coordinates
            WHERE start_pos BETWEEN ? AND ?
            ORDER BY start_pos ASC
            """,
            (min_pos, max_pos)
        )
        return [dict(r) for r in cursor.fetchall()]

    def get_genomic_track_data(self, locus_tag: str, window_bp: int = 10000):
        """
        Query 5-track genomic data (CDS genes, TSS promoter sites, TFBS/ChIP-seq binding peaks, sRNAs)
        within +/- window_bp of locus_tag.
        """
        target_coords = self.get_gene_coordinates(locus_tag)
        if not target_coords or not target_coords.get("start_pos"):
            return None

        conn = self.get_connection()
        if not conn:
            return None
        cursor = conn.cursor()

        center = target_coords["start_pos"]
        min_pos = max(1, center - window_bp)
        max_pos = center + window_bp

        # 1. Neighboring CDS genes
        cursor.execute(
            """
            SELECT locus_tag, gene_name, start_pos, end_pos, strand, gene_length, tss_position, promoter_70bp
            FROM gene_coordinates
            WHERE (start_pos BETWEEN ? AND ?) OR (end_pos BETWEEN ? AND ?)
            ORDER BY start_pos ASC
            """,
            (min_pos, max_pos, min_pos, max_pos)
        )
        genes = [dict(r) for r in cursor.fetchall()]

        # 2. CollectTF / TFBS & Experimental ChIP-seq binding peaks in window
        peaks = []
        try:
            cursor.execute(
                """
                SELECT tf_name, locus_tag, score, site_seq, rel_pos
                FROM collectf_tfbs
                WHERE locus_tag = ? OR LOWER(tf_name) = ?
                """,
                (locus_tag.lower().strip(), locus_tag.lower().strip())
            )
            for r in cursor.fetchall():
                d = dict(r)
                d["pos"] = center + (d.get("rel_pos") or 0)
                d["source"] = "CollectTF"
                peaks.append(d)
        except Exception:
            pass

        try:
            cursor.execute(
                """
                SELECT peak_id, tf_id, tf_name, peak_start, peak_end, peak_center,
                       peak_score, peak_signal, neglog10q, strength_tier, nearest_gene_locus,
                       rel_pos_to_tss, spatial_confidence, genomic_region_chip
                FROM chipseq_peaks
                WHERE (peak_start BETWEEN ? AND ?) OR (peak_end BETWEEN ? AND ?)
                   OR nearest_gene_locus = ? OR LOWER(tf_name) = ?
                ORDER BY peak_start ASC
                """,
                (min_pos, max_pos, min_pos, max_pos, locus_tag, locus_tag.lower())
            )
            for r in cursor.fetchall():
                d = dict(r)
                d["pos"] = d.get("peak_center") or d.get("peak_start")
                d["score"] = d.get("peak_score") or d.get("peak_signal") or 1.0
                d["source"] = "ChIP-seq (Internal/Experimental)"
                peaks.append(d)
        except Exception:
            pass

        # 3. sRNA / ncRNA annotations in window
        rnas = []
        try:
            cursor.execute(
                """
                SELECT rna_id, rna_name, rna_type, start_pos, end_pos, strand
                FROM rfam_ncrna
                WHERE start_pos BETWEEN ? AND ?
                """,
                (min_pos, max_pos)
            )
            rnas = [dict(r) for r in cursor.fetchall()]
        except Exception:
            pass

        return {
            "query_locus": locus_tag,
            "window": {"min_pos": min_pos, "max_pos": max_pos, "center_pos": center},
            "target": target_coords,
            "genes": genes,
            "peaks": peaks,
            "rnas": rnas
        }

    def get_gene_chipseq_peaks(self, locus_tag: str):
        """Query all ChIP-seq binding peaks mapped to target locus_tag."""
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        try:
            cursor.execute(
                """
                SELECT * FROM chipseq_peaks
                WHERE nearest_gene_locus = ? OR gene_list LIKE ?
                ORDER BY peak_score DESC
                """,
                (locus_tag, f"%{locus_tag}%")
            )
            return [dict(r) for r in cursor.fetchall()]
        except Exception:
            return []

    def get_tf_chipseq_peaks(self, tf_identifier: str):
        """Query all ChIP-seq binding peaks for a given TF (locus_tag or name)."""
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        try:
            cursor.execute(
                """
                SELECT * FROM chipseq_peaks
                WHERE LOWER(tf_id) = ? OR LOWER(tf_name) = ?
                ORDER BY peak_score DESC
                """,
                (tf_identifier.lower().strip(), tf_identifier.lower().strip())
            )
            return [dict(r) for r in cursor.fetchall()]
        except Exception:
            return []

    def get_allosteric_feedback_loops(self, tf_or_metabolite: str = None):
        """
        Query v_metabolite_tf_feedback for closed-loop metabolite-TF allosteric regulation.
        """
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        if tf_or_metabolite:
            term = f"%{tf_or_metabolite.lower().strip()}%"
            cursor.execute(
                """
                SELECT * FROM v_metabolite_tf_feedback
                WHERE LOWER(tf_locus) LIKE ? OR LOWER(tf_name) LIKE ? OR LOWER(effector_molecule) LIKE ?
                ORDER BY score DESC LIMIT 50
                """,
                (term, term, term)
            )
        else:
            cursor.execute("SELECT * FROM v_metabolite_tf_feedback ORDER BY score DESC LIMIT 50")

        rows = cursor.fetchall()
        res = []
        for r in rows:
            item = dict(r)
            if item.get("details"):
                try:
                    item["details"] = json.loads(item["details"])
                except Exception:
                    pass
            res.append(item)
        return res

    def get_srna_target_competition(self, srna_id: str = None):
        """
        Query v_srna_competition_ranking for sRNA target binding energy and competition ranks.
        """
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        if srna_id:
            cursor.execute(
                "SELECT * FROM v_srna_competition_ranking WHERE LOWER(srna_id) = ? ORDER BY binding_energy ASC LIMIT 50",
                (srna_id.lower().strip(),)
            )
        else:
            cursor.execute("SELECT * FROM v_srna_competition_ranking ORDER BY binding_energy ASC LIMIT 50")

        rows = cursor.fetchall()
        res = []
        for r in rows:
            item = dict(r)
            if item.get("details"):
                try:
                    item["details"] = json.loads(item["details"])
                except Exception:
                    pass
            res.append(item)
        return res

    def get_imodulons_for_gene(self, gene_locus: str):
        """
        Query iModulon co-expression modules and weights for a given gene.
        """
        if not gene_locus:
            return []
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        cursor.execute(
            """
            SELECT w.imodulon_id, w.weight, m.name, m.explanation
            FROM imodulon_gene_weights w
            LEFT JOIN imodulons m ON w.imodulon_id = m.imodulon_id
            WHERE w.gene_locus = ?
            ORDER BY ABS(w.weight) DESC
            """,
            (gene_locus.lower().strip(),)
        )
        return [dict(r) for r in cursor.fetchall()]

    def get_rf_edge_scores(self, locus: str, min_confidence: float = 0.3):
        """
        Query Random Forest machine learning edge confidence scores for a given TF or Target locus.
        """
        if not locus:
            return []
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        g_lower = locus.lower().strip()
        cursor.execute(
            """
            SELECT * FROM tf_gene_rf_scores
            WHERE (LOWER(tf_locus) = ? OR LOWER(target_locus) = ?) AND predicted_confidence >= ?
            ORDER BY predicted_confidence DESC LIMIT 50
            """,
            (g_lower, g_lower, min_confidence)
        )
        return [dict(r) for r in cursor.fetchall()]

    def get_tf_hierarchy_rankings(self):
        """
        Query TF 3-tier pyramid hierarchy rankings (Master/Top, Middle, Bottom).
        """
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM tf_hierarchy_rankings ORDER BY rank ASC")
        return [dict(r) for r in cursor.fetchall()]

    def get_rewired_edges(self, locus: str = None):
        """
        Query cross-strain evolutionary rewired regulation edges.
        """
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        if locus:
            g_lower = locus.lower().strip()
            cursor.execute(
                "SELECT * FROM network_rewired_edges WHERE LOWER(tf_cgl) = ? OR LOWER(tg_cgl) = ? OR LOWER(tf_name) = ? OR LOWER(tg_name) = ?",
                (g_lower, g_lower, g_lower, g_lower)
            )
        else:
            cursor.execute("SELECT * FROM network_rewired_edges LIMIT 100")
        return [dict(r) for r in cursor.fetchall()]

    def get_collectf_tfbs(self, locus: str = None):
        """
        Query CollecTF experimentally validated TFBS records deduplicated against RegPrecise.
        """
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        if locus:
            g_lower = locus.lower().strip()
            cursor.execute(
                """
                SELECT * FROM collectf_tfbs
                WHERE LOWER(tf_locus) = ? OR LOWER(target_locus) = ? OR LOWER(tf_name) = ? OR LOWER(target_name) = ?
                ORDER BY site_id ASC
                """,
                (g_lower, g_lower, g_lower, g_lower)
            )
        else:
            cursor.execute("SELECT * FROM collectf_tfbs ORDER BY site_id ASC LIMIT 100")
        return [dict(r) for r in cursor.fetchall()]

    def get_pathways_for_gene(self, locus: str):
        """
        Query all BioCyc/KEGG pathways associated with a given gene.
        """
        if not locus:
            return []
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        g_lower = locus.lower().strip()
        cursor.execute(
            """
            SELECT m.gene_locus, m.pathway_id, m.pathway_name, p.category, p.gene_list
            FROM gene_pathway_mappings m
            LEFT JOIN biocyc_kegg_pathways p ON m.pathway_id = p.pathway_id
            WHERE LOWER(m.gene_locus) = ?
            ORDER BY m.pathway_name ASC
            """,
            (g_lower,)
        )
        return [dict(r) for r in cursor.fetchall()]

    def get_genes_in_pathway(self, pathway_id: str):
        """
        Query all genes belonging to a given BioCyc/KEGG pathway ID or name.
        """
        if not pathway_id:
            return []
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        p_lower = pathway_id.lower().strip()
        cursor.execute(
            """
            SELECT m.gene_locus, g.cgl_locus, g.gene_name, g.product, m.pathway_id, m.pathway_name
            FROM gene_pathway_mappings m
            LEFT JOIN gene_mappings g ON LOWER(m.gene_locus) = LOWER(g.cg_locus)
            WHERE LOWER(m.pathway_id) = ? OR LOWER(m.pathway_name) LIKE ?
            ORDER BY m.gene_locus ASC
            """,
            (p_lower, f"%{p_lower}%")
        )
        return [dict(r) for r in cursor.fetchall()]

    def get_ncrnas(self, rna_type: str = None):
        """
        Query Rfam / ENA non-coding RNA annotations (sRNA, Riboswitch, 6S RNA, tmRNA, RNase P).
        """
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        if rna_type:
            cursor.execute(
                "SELECT * FROM rfam_ncrnas WHERE LOWER(rna_type) LIKE ? ORDER BY start_pos ASC",
                (f"%{rna_type.lower().strip()}%",)
            )
        else:
            cursor.execute("SELECT * FROM rfam_ncrnas ORDER BY start_pos ASC")
        return [dict(r) for r in cursor.fetchall()]

    def get_srna_targets(self, locus: str = None):
        """
        Query sRNA-mRNA target interaction mechanics, IntaRNA binding energy, and regulatory mechanisms.
        """
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        if locus:
            g_lower = locus.lower().strip()
            cursor.execute(
                """
                SELECT * FROM ncrna_target_interactions
                WHERE LOWER(srna_id) = ? OR LOWER(target_locus) = ? OR LOWER(target_name) = ?
                ORDER BY binding_energy_kcal ASC LIMIT 50
                """,
                (g_lower, g_lower, g_lower)
            )
        else:
            cursor.execute("SELECT * FROM ncrna_target_interactions ORDER BY binding_energy_kcal ASC LIMIT 50")
        return [dict(r) for r in cursor.fetchall()]

    def get_condition_specific_regulons(self, condition_name: str = None):
        """
        Query condition-specific active iModulons and regulons across growth conditions.
        """
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        if condition_name:
            c_lower = condition_name.lower().strip()
            cursor.execute(
                """
                SELECT * FROM imodulon_condition_activities
                WHERE LOWER(condition_name) LIKE ? AND is_significant = 1
                ORDER BY ABS(activity_score) DESC
                """,
                (f"%{c_lower}%",)
            )
        else:
            cursor.execute("SELECT * FROM imodulon_condition_activities WHERE is_significant = 1 ORDER BY ABS(activity_score) DESC LIMIT 50")
        return [dict(r) for r in cursor.fetchall()]

    def get_imodulon_regulon_overlap(self, imodulon_id: str = None):
        """
        Query iModulon-regulon structural overlap matrix, F1-scores, and TF alignments.
        """
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        if imodulon_id:
            i_lower = imodulon_id.lower().strip()
            cursor.execute(
                "SELECT * FROM imodulon_regulon_overlaps WHERE LOWER(imodulon_id) = ? OR LOWER(tf_locus) = ? OR LOWER(tf_name) = ?",
                (i_lower, i_lower, i_lower)
            )
        else:
            cursor.execute("SELECT * FROM imodulon_regulon_overlaps ORDER BY f1_score DESC LIMIT 50")
        return [dict(r) for r in cursor.fetchall()]

    def get_condition_regulation_runs(self):
        """List the available condition-specific regulation analysis modules."""
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        cursor.execute("""
            SELECT r.run_id, r.method_version, r.scope, r.created_at, r.notes,
                   COUNT(DISTINCT s.comparison_id) AS comparison_count,
                   COUNT(DISTINCT s.tf_locus) AS transcription_factor_count
            FROM condition_analysis_runs r
            LEFT JOIN condition_regulon_summary s ON s.run_id=r.run_id
            GROUP BY r.run_id, r.method_version, r.scope, r.created_at, r.notes
            ORDER BY r.run_id
        """)
        return [dict(row) for row in cursor.fetchall()]

    def get_condition_regulation_conditions(self, run_id: str = "iron_regulon_v1"):
        """List scored contrasts and available TF analyses for one module."""
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        cursor.execute("""
            SELECT comparison_id, condition_label,
                   GROUP_CONCAT(tf_name, ',') AS transcription_factors,
                   MAX(CASE WHEN validation_status='response_enriched_fdr10' THEN 1 ELSE 0 END) AS has_fdr_signal,
                   MAX(scored_edge_count) AS max_scored_edges
            FROM v_condition_regulon_response
            WHERE run_id = ?
            GROUP BY comparison_id, condition_label
            ORDER BY condition_label COLLATE NOCASE
        """, (run_id,))
        rows = []
        for row in cursor.fetchall():
            item = dict(row)
            item["transcription_factors"] = [
                value for value in (item.get("transcription_factors") or "").split(",") if value
            ]
            item["has_fdr_signal"] = bool(item.get("has_fdr_signal"))
            rows.append(item)
        return rows

    def get_condition_regulation_summary(
        self, comparison_id: str = None, tf_name: str = None,
        run_id: str = "iron_regulon_v1",
    ):
        """Return TF activity, enrichment and regulon-level summaries."""
        conn = self.get_connection()
        if not conn:
            return []
        cursor = conn.cursor()
        clauses = ["run_id = ?"]
        params = [run_id]
        if comparison_id:
            clauses.append("comparison_id = ?")
            params.append(comparison_id)
        if tf_name:
            clauses.append("LOWER(tf_name) = ?")
            params.append(tf_name.lower().strip())
        where = " WHERE " + " AND ".join(clauses) if clauses else ""
        cursor.execute(
            "SELECT * FROM v_condition_regulon_response" + where +
            " ORDER BY condition_label COLLATE NOCASE, tf_name COLLATE NOCASE",
            params,
        )
        rows = []
        for row in cursor.fetchall():
            item = dict(row)
            try:
                item["top_targets"] = json.loads(item.pop("top_targets_json") or "[]")
            except (TypeError, ValueError):
                item["top_targets"] = []
            rows.append(item)
        return rows

    def get_condition_regulation_edges(
        self,
        comparison_id: str,
        tf_name: str = None,
        support_state: str = None,
        min_score: float = 0.0,
        limit: int = 100,
        offset: int = 0,
        run_id: str = "iron_regulon_v1",
    ):
        """Query paginated, condition-specific regulatory edge scores."""
        conn = self.get_connection()
        if not conn or not comparison_id:
            return {"total": 0, "edges": []}
        cursor = conn.cursor()
        allowed_states = {
            "condition_supported", "direction_conflict", "weak_context_support",
            "insufficient_dynamic_data",
        }
        clauses = ["run_id = ?", "comparison_id = ?", "condition_score >= ?"]
        params = [run_id, comparison_id, max(0.0, min(1.0, float(min_score)))]
        if tf_name:
            clauses.append("LOWER(tf_name) = ?")
            params.append(tf_name.lower().strip())
        if support_state in allowed_states:
            clauses.append("support_state = ?")
            params.append(support_state)
        where = " WHERE " + " AND ".join(clauses)
        cursor.execute("SELECT COUNT(*) FROM v_condition_regulation_top_edges" + where, params)
        total = cursor.fetchone()[0]
        safe_limit = max(1, min(int(limit), 500))
        safe_offset = max(0, int(offset))
        cursor.execute(
            "SELECT * FROM v_condition_regulation_top_edges" + where +
            " ORDER BY condition_score DESC, ABS(COALESCE(target_expression_mean, 0)) DESC LIMIT ? OFFSET ?",
            params + [safe_limit, safe_offset],
        )
        return {"total": total, "edges": [dict(row) for row in cursor.fetchall()]}

    def get_intervention_targets(
        self, query: str = None, strategy: str = None, min_modules: int = 1,
        max_risk: float = 1.0, evidence_grade: str = None,
        include_known_essential: bool = True, limit: int = 100, offset: int = 0,
    ):
        """Return paginated cross-module engineering target priorities."""
        conn = self.get_connection()
        if not conn:
            return {"total": 0, "targets": []}
        clauses = ["module_count >= ?", "risk_score <= ?"]
        params = [max(1, int(min_modules)), max(0.0, min(1.0, float(max_risk)))]
        allowed_strategies = {
            "dynamic_tuning_only", "careful_titration", "multi_stress_control_node",
            "metabolic_intervention_candidate", "context_specific_candidate",
        }
        if strategy in allowed_strategies:
            clauses.append("strategy_class = ?")
            params.append(strategy)
        if evidence_grade and evidence_grade.upper() in {"A", "B", "C", "D"}:
            clauses.append("evidence_grade = ?")
            params.append(evidence_grade.upper())
        if not include_known_essential:
            clauses.append("essentiality_status <> 'known_essential'")
        if query:
            clauses.append("(LOWER(target_locus) LIKE ? OR LOWER(COALESCE(target_name,'')) LIKE ? OR LOWER(COALESCE(product,'')) LIKE ?)")
            token = f"%{query.lower().strip()}%"
            params.extend([token, token, token])
        where = " WHERE " + " AND ".join(clauses)
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM v_intervention_target_priorities" + where, params)
        total = cursor.fetchone()[0]
        safe_limit = max(1, min(int(limit), 500))
        safe_offset = max(0, int(offset))
        cursor.execute(
            "SELECT * FROM v_intervention_target_priorities" + where
            + " ORDER BY priority_score DESC, target_locus LIMIT ? OFFSET ?",
            params + [safe_limit, safe_offset],
        )
        targets = []
        for row in cursor.fetchall():
            item = dict(row)
            item["modules"] = [value for value in (item.get("modules") or "").split(",") if value]
            try:
                item["rationale"] = json.loads(item.pop("rationale_json") or "{}")
            except (TypeError, ValueError):
                item["rationale"] = {}
            targets.append(item)
        return {"total": total, "targets": targets}

    def get_intervention_target_detail(self, locus: str):
        """Return one priority score and its per-module evidence decomposition."""
        conn = self.get_connection()
        if not conn or not locus:
            return None
        cursor = conn.cursor()
        cursor.execute(
            "SELECT * FROM v_intervention_target_priorities WHERE LOWER(target_locus)=LOWER(?) LIMIT 1",
            (locus,),
        )
        row = cursor.fetchone()
        if not row:
            return None
        target = dict(row)
        target["modules"] = [value for value in (target.get("modules") or "").split(",") if value]
        try:
            target["rationale"] = json.loads(target.pop("rationale_json") or "{}")
        except (TypeError, ValueError):
            target["rationale"] = {}
        cursor.execute("""
            SELECT module_run_id, condition_count, regulator_count, supported_count,
                   conflict_count, significant_context_count, mean_score, max_score
            FROM intervention_target_module_evidence
            WHERE run_id='cross_module_priority_v1' AND LOWER(target_locus)=LOWER(?)
            ORDER BY mean_score DESC
        """, (locus,))
        target["module_evidence"] = [dict(item) for item in cursor.fetchall()]
        return target

# Export singleton instance getter
def get_db_manager(db_path=None):
    return CglDatabaseManager(db_path=db_path)
