"""
scripts/enrich_regulations.py
================================
Add evidence_score column to regulations.csv based on Evidence field.

Mapping:
  'experimental'              -> 1.0  (HIGH)
  'experimental + predicted'  -> 0.8  (HIGH)
  'predicted'                 -> 0.4  (MEDIUM)
  ''  (unknown)               -> 0.2  (LOW)

Also adds a 'confidence_label' column: HIGH / MEDIUM / LOW
"""
import os, csv

ROOT    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN_CSV  = os.path.join(ROOT, "data", "reference", "regulations.csv")
OUT_CSV = IN_CSV  # overwrite in-place

EVIDENCE_MAP = {
    "experimental":             (1.0, "HIGH"),
    "experimental + predicted": (0.8, "HIGH"),
    "predicted":                (0.4, "MEDIUM"),
}

def main():
    with open(IN_CSV, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    
    fieldnames = list(rows[0].keys())
    if "evidence_score" not in fieldnames:
        fieldnames.append("evidence_score")
    if "confidence_label" not in fieldnames:
        fieldnames.append("confidence_label")

    n_high = n_medium = n_low = 0
    for row in rows:
        ev = row.get("Evidence", "").strip().lower()
        score, label = EVIDENCE_MAP.get(ev, (0.2, "LOW"))
        row["evidence_score"]   = score
        row["confidence_label"] = label
        if label == "HIGH":    n_high   += 1
        elif label == "MEDIUM": n_medium += 1
        else:                   n_low    += 1

    with open(OUT_CSV, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Written {len(rows)} rows to {OUT_CSV}")
    print(f"  HIGH   (experimental)       : {n_high}")
    print(f"  MEDIUM (predicted)          : {n_medium}")
    print(f"  LOW    (unknown/empty)      : {n_low}")

if __name__ == "__main__":
    main()
