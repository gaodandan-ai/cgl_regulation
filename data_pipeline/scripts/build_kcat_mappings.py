"""
scripts/build_kcat_mappings.py
================================
Build brenda_kcat_mappings.json by merging:
  Layer 1 (Priority): existing brenda_kcat_mappings.json (11 curated)
  Layer 2 (Fill):     dlkcat_predicted_kcat.json (1850 predicted, source='experimental')

Output format (per reaction):
  {
    "kcat": float,           # 1/s
    "source": str,           # "BRENDA" | "DLKcat"
    "confidence": str,       # "HIGH" | "MEDIUM"
    "reference": str         # citation or "DLKcat ML prediction"
  }
"""
import os, json

ROOT       = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BRENDA_IN  = os.path.join(ROOT, "data", "reference", "brenda_kcat_mappings.json")
DLKCAT_IN  = os.path.join(ROOT, "data", "reference", "dlkcat_predicted_kcat.json")
OUT        = BRENDA_IN   # overwrite brenda file with merged result

def main():
    # Load Layer 1: existing BRENDA (priority)
    with open(BRENDA_IN, encoding="utf-8") as f:
        brenda = json.load(f)
    print(f"Loaded {len(brenda)} BRENDA entries (priority layer)")

    # Normalize BRENDA entries
    brenda_clean = {}
    for rxn_id, entry in brenda.items():
        brenda_clean[rxn_id] = {
            "kcat":       float(entry.get("kcat", 0)),
            "source":     "BRENDA",
            "confidence": "HIGH",
            "reference":  entry.get("reference", "BRENDA database"),
        }

    # Load Layer 2: DLKcat predictions
    with open(DLKCAT_IN, encoding="utf-8") as f:
        dlkcat = json.load(f)
    print(f"Loaded {len(dlkcat)} DLKcat entries (fill layer)")

    # Merge: BRENDA takes priority
    merged = dict(brenda_clean)
    n_added = n_skipped = 0
    for rxn_id, entry in dlkcat.items():
        if rxn_id in merged:
            n_skipped += 1
            continue
        kcat_val = entry.get("kcat") if isinstance(entry, dict) else float(entry)
        src      = entry.get("source", "DLKcat") if isinstance(entry, dict) else "DLKcat"
        if kcat_val is None or kcat_val <= 0:
            continue
        merged[rxn_id] = {
            "kcat":       round(float(kcat_val), 4),
            "source":     "DLKcat",
            "confidence": "MEDIUM",   # ML prediction
            "reference":  "DLKcat ML prediction (Li et al. 2022, Nat Commun)",
        }
        n_added += 1

    print(f"Merged: {len(merged)} total ({len(brenda_clean)} BRENDA + {n_added} DLKcat, {n_skipped} DLKcat skipped - BRENDA has priority)")

    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(merged, f, indent=2, ensure_ascii=False)
    print(f"Written: {OUT}")

    # Stats
    sources = {}
    for e in merged.values():
        sources[e["source"]] = sources.get(e["source"], 0) + 1
    for src, cnt in sorted(sources.items()):
        print(f"  {src}: {cnt}")

if __name__ == "__main__":
    main()
