"""
scripts/build_string_data.py
============================
Download full STRING PPI network for Corynebacterium glutamicum ATCC 13032
(taxon ID 196627) and convert to the project's string_interactions.json format.

STRING API used:
  /api/json/network  (no API key required)
  IDs format: "196627.cgXXXX"

Output: data/reference/string_interactions.json
  {
    "cgXXXX": [
      {
        "partner":      "cgYYYY",
        "score":        875,          # 0-1000
        "neighborhood": 0,            # nscore * 1000
        "fusion":       0,            # fscore * 1000
        "cooccurrence": 210,          # pscore * 1000
        "coexpression": 320,          # ascore * 1000
        "experimental": 650,          # escore * 1000
        "database":     0,            # dscore * 1000
        "textmining":   180,          # tscore * 1000
        "type":         "experimental" # dominant channel
      }
    ]
  }

Usage:
  python scripts/build_string_data.py                # score >= 400 (medium)
  python scripts/build_string_data.py --min-score 700  # high confidence only
"""

import os
import json
import math
import time
import logging
import argparse
import urllib.request
import urllib.parse
import urllib.error
import csv as csv_mod

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("build_string")

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR     = os.path.dirname(SCRIPT_DIR)
MAPPING_CSV  = os.path.join(ROOT_DIR, "data", "reference", "gene_mapping.csv")
OUTPUT_PATH  = os.path.join(ROOT_DIR, "data", "reference", "string_interactions.json")

TAXON_ID     = 196627          # C. glutamicum ATCC 13032
CALLER_ID    = "cgl_regulation_app"
BATCH_SIZE   = 500             # proteins per API call (STRING limit)
API_DELAY    = 1.0             # seconds between batches (rate limiting)
STRING_BASE  = "https://string-db.org/api/json"


# ── Channel name mapping ──────────────────────────────────────────────────────
CHANNEL_KEYS = {
    "nscore": "neighborhood",
    "fscore": "fusion",
    "pscore": "cooccurrence",
    "ascore": "coexpression",
    "escore": "experimental",
    "dscore": "database",
    "tscore": "textmining",
}

CHANNEL_PRIORITY = [
    "experimental", "database", "coexpression",
    "neighborhood", "cooccurrence", "textmining", "fusion",
]


def _dominant_type(row_scores: dict) -> str:
    """Pick the highest-scoring channel as the interaction type label."""
    best = max(CHANNEL_PRIORITY, key=lambda ch: row_scores.get(ch, 0))
    if row_scores.get(best, 0) > 50:
        return best
    return "functional association"


def load_cg_loci() -> list:
    """Load all cg locus tags from gene_mapping.csv."""
    loci = set()
    try:
        with open(MAPPING_CSV, encoding="utf-8") as f:
            reader = csv_mod.DictReader(f)
            for row in reader:
                cg = row.get("cg_locus", "").strip().lower()
                if cg and cg.startswith("cg") and cg[2:].isdigit():
                    loci.add(cg)
    except Exception as e:
        logger.error(f"Cannot read {MAPPING_CSV}: {e}")
        raise
    loci = sorted(loci)
    logger.info(f"Loaded {len(loci)} CG locus tags from gene_mapping.csv")
    return loci


def fetch_string_batch(loci_batch: list, min_score_int: int) -> list:
    """
    Query STRING /api/json/network for a batch of loci.
    Returns list of edge dicts.
    """
    identifiers = "%0d".join(f"{TAXON_ID}.{locus}" for locus in loci_batch)
    params = urllib.parse.urlencode({
        "species":         TAXON_ID,
        "required_score":  min_score_int,
        "caller_identity": CALLER_ID,
        "network_type":    "functional",
    })
    url = f"{STRING_BASE}/network?identifiers={identifiers}&{params}"

    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:
            logger.warning(f"HTTP {e.code} on attempt {attempt+1}; retrying...")
            time.sleep(2 ** attempt)
        except Exception as e:
            logger.warning(f"Request error on attempt {attempt+1}: {e}; retrying...")
            time.sleep(2 ** attempt)
    logger.error("All retries exhausted for batch.")
    return []


def build_string_data(min_score: int = 400):
    """
    Download full STRING network and write string_interactions.json.

    min_score : 0-1000 combined score threshold
    """
    logger.info(f"Building STRING data (taxon={TAXON_ID}, min_score={min_score})...")

    loci = load_cg_loci()
    n_batches = math.ceil(len(loci) / BATCH_SIZE)

    # Collect all edges (deduplicated)
    # Store as set of frozensets to avoid counting both A→B and B→A twice
    edge_map: dict = {}   # (locusA, locusB) -> merged score entry (keep highest score)

    for batch_idx in range(n_batches):
        batch = loci[batch_idx * BATCH_SIZE:(batch_idx + 1) * BATCH_SIZE]
        logger.info(f"  Batch {batch_idx + 1}/{n_batches}: {batch[0]}..{batch[-1]}")

        edges = fetch_string_batch(batch, min_score)
        logger.info(f"    → {len(edges)} edges returned")

        for edge in edges:
            # Parse locus tags from stringId: "196627.cg0001" → "cg0001"
            id_a = edge.get("stringId_A", "").split(".")[-1].lower()
            id_b = edge.get("stringId_B", "").split(".")[-1].lower()

            if not id_a.startswith("cg") or not id_b.startswith("cg"):
                continue
            if id_a == id_b:
                continue  # skip self-loops

            raw_score = edge.get("score", 0)
            score_int = round(raw_score * 1000) if raw_score <= 1 else int(raw_score)

            # Build sub-scores dict
            sub = {}
            for api_key, ch_name in CHANNEL_KEYS.items():
                val = edge.get(api_key, 0)
                sub[ch_name] = round(val * 1000) if val <= 1 else int(val)

            # Canonical edge key (alphabetical order)
            key = (min(id_a, id_b), max(id_a, id_b))
            existing = edge_map.get(key)
            if existing is None or score_int > existing["score"]:
                edge_map[key] = {
                    "a": id_a, "b": id_b,
                    "score": score_int,
                    **sub,
                }

        if batch_idx < n_batches - 1:
            time.sleep(API_DELAY)

    logger.info(f"Total unique edges collected: {len(edge_map)}")

    # Build adjacency list
    adjacency: dict = {}
    for key, entry in edge_map.items():
        a, b = entry["a"], entry["b"]
        score_int = entry["score"]
        sub_scores = {ch: entry.get(ch, 0) for ch in CHANNEL_PRIORITY}
        dom_type = _dominant_type(sub_scores)

        record = {
            "partner":      b,
            "score":        score_int,
            **{ch: sub_scores[ch] for ch in CHANNEL_PRIORITY},
            "type":         dom_type,
        }
        adjacency.setdefault(a, []).append(record)

        # Reciprocal
        record_rev = {
            "partner":      a,
            "score":        score_int,
            **{ch: sub_scores[ch] for ch in CHANNEL_PRIORITY},
            "type":         dom_type,
        }
        adjacency.setdefault(b, []).append(record_rev)

    # Sort each gene's partner list by score descending
    for locus in adjacency:
        adjacency[locus].sort(key=lambda x: -x["score"])

    n_genes    = len(adjacency)
    n_edges    = len(edge_map)
    n_high     = sum(1 for e in edge_map.values() if e["score"] >= 700)
    n_medium   = n_edges - n_high

    logger.info("=" * 60)
    logger.info(f"Genes with ≥1 interaction : {n_genes}")
    logger.info(f"Unique edges (score≥{min_score}) : {n_edges}")
    logger.info(f"  High-confidence (≥700)  : {n_high}")
    logger.info(f"  Medium (400-699)         : {n_medium}")
    logger.info("=" * 60)

    output = {
        "_meta": {
            "description": f"STRING v12 PPI network for C. glutamicum ATCC 13032 (taxon {TAXON_ID})",
            "source":       "https://string-db.org",
            "version":      "12.0",
            "taxon_id":     TAXON_ID,
            "min_score":    min_score,
            "n_genes":      n_genes,
            "n_edges":      n_edges,
            "n_high_conf":  n_high,
            "n_medium_conf": n_medium,
            "generated":    "2026-07-13",
            "score_channels": {
                "neighborhood": "gene neighborhood (nscore)",
                "fusion":       "gene fusion (fscore)",
                "cooccurrence": "phylogenetic co-occurrence (pscore)",
                "coexpression": "co-expression (ascore)",
                "experimental": "experimental evidence (escore)",
                "database":     "curated databases (dscore)",
                "textmining":   "text mining (tscore)",
            },
        },
        **adjacency,
    }

    with open(OUTPUT_PATH, "w", encoding="utf-8") as f:
        json.dump(output, f, ensure_ascii=False, separators=(",", ":"))
    size_mb = os.path.getsize(OUTPUT_PATH) / 1e6
    logger.info(f"Written: {OUTPUT_PATH} ({size_mb:.1f} MB)")
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build STRING interaction data for C. glutamicum")
    parser.add_argument("--min-score", type=int, default=400,
                        help="Minimum combined score (0-1000). Default: 400")
    args = parser.parse_args()
    build_string_data(min_score=args.min_score)
