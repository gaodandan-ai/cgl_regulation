"""
process_regprecise_sites.py
===========================
Converts RegPrecise ATCC 13032 regulatory sites FASTA to:
  1. regprecise_binding_sites.tsv  -- one row per binding site instance
  2. regprecise_pwm.json           -- per-TF consensus PWM
  3. Updates regprecise_regulations.csv with binding_site, evidence_level, source_accession

Run from the project root:
    python data_pipeline/scripts/process_regprecise_sites.py
"""

import re, csv, json, math
from collections import defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent.parent   # f:/cgl_regulation
DATA = ROOT / "data" / "reference"

FASTA_PATH   = Path(r"C:\Users\Tsuki\Documents\Codex\2026-07-21\interconnected-iron-and-heme-regulatory-networks\outputs\other_tf_resources\RegPrecise_Cglutamicum\RegPrecise_ATCC13032_regulatory_sites_official.fasta")
CHIPSEQ_PATH = DATA / "chipseq_regulations.csv"
TSS_PATH     = DATA / "tss_promoter_annotations.json"
GENE_MAP     = DATA / "gene_mapping.csv"

OUT_TSV  = DATA / "regprecise_binding_sites.tsv"
OUT_PWM  = DATA / "regprecise_pwm.json"
OUT_REGS = DATA / "regprecise_regulations.csv"   # updated in-place


# ---------------------------------------------------------------------------
# 1. Parse FASTA
# ---------------------------------------------------------------------------
HEADER_RE = re.compile(
    r'^>([\w-]+)(?:\(([\w-]+)\))?\s+Score=([\d.]+)\s+Pos=(-?\d+)\s+\[([^\]]+)\]'
)
TF_RE = re.compile(r'^#\s*TF\s*-\s*(\w+):\s*([\w,\s]+)')


def parse_fasta(fasta_path: Path) -> list:
    records = []
    current_tf_name = ""
    current_tf_locus = ""
    pending_header = None

    with open(fasta_path, 'r', encoding='utf-8') as f:
        for raw in f:
            line = raw.rstrip()
            if not line:
                continue

            m = TF_RE.match(line)
            if m:
                current_tf_name  = m.group(1).strip()
                current_tf_locus = m.group(2).strip().split(',')[0].strip()
                pending_header   = None
                continue

            if line.startswith('>'):
                hm = HEADER_RE.match(line)
                if hm:
                    pending_header = {
                        'tf_name':   current_tf_name,
                        'tf_locus':  current_tf_locus,
                        'tg_locus':  hm.group(1).lower().replace('-', '_'),
                        'tg_gene':   (hm.group(2) or hm.group(1)),
                        'score':     float(hm.group(3)),
                        'pos_bp':    int(hm.group(4)),
                        'genome':    hm.group(5).strip(),
                    }
                else:
                    print(f"  WARN unmatched header: {line}")
                    pending_header = None
                continue

            if pending_header is not None and re.match(r'^[ACGTacgtNn]+$', line):
                rec = dict(pending_header)
                rec['sequence']  = line.upper()
                rec['site_width'] = len(line)
                records.append(rec)
                pending_header = None

    return records


# ---------------------------------------------------------------------------
# 2. Gene-name alias map
# ---------------------------------------------------------------------------
def load_gene_name_map(gene_map_path: Path) -> dict:
    gmap = {}
    if not gene_map_path.exists():
        return gmap
    with open(gene_map_path, encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):
            locus = (row.get('cg_locus') or row.get('locusTag') or '').lower().strip()
            name  = (row.get('gene_name') or row.get('name') or '').strip()
            if locus and name:
                gmap[locus] = name
    return gmap


# ---------------------------------------------------------------------------
# 3. TSS annotations
# ---------------------------------------------------------------------------
def load_tss(tss_path: Path) -> dict:
    with open(tss_path, encoding='utf-8') as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# 4. ChIP-seq evidence
# ---------------------------------------------------------------------------
def load_chipseq(chip_path: Path):
    chip_pairs = {}         # (tf_locus_lower, tg_locus_lower) -> [hits]
    tf_name_map = {}        # tf_name_lower -> tf_locus_lower

    with open(chip_path, encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):
            tf_l = (row.get('TF_locusTag') or '').lower().strip()
            tf_n = (row.get('TF_name') or '').lower().strip()
            tg_l = (row.get('TG_locusTag') or '').lower().strip()
            if not tf_l or not tg_l:
                continue
            key = (tf_l, tg_l)
            chip_pairs.setdefault(key, []).append({
                'evidence':     row.get('Evidence', ''),
                'source':       row.get('Source', ''),
                'pmid':         row.get('PMID', ''),
                'strain_group': row.get('strain_group', ''),
                'confidence':   row.get('confidence_label', ''),
            })
            if tf_n:
                tf_name_map[tf_n] = tf_l

    return chip_pairs, tf_name_map


# ---------------------------------------------------------------------------
# 5. PWM computation
# ---------------------------------------------------------------------------
def compute_pwm(sequences: list):
    if not sequences:
        return {}, '', 0.0

    # Use the most common width to avoid multi-width contamination (e.g. NrtR)
    from collections import Counter
    width_counts = Counter(len(s) for s in sequences)
    majority_width = width_counts.most_common(1)[0][0]
    seqs = [s for s in sequences if len(s) == majority_width]
    if len(seqs) < len(sequences):
        print(f"    [PWM] Dropped {len(sequences)-len(seqs)} off-width sequences (kept w={majority_width})")

    counts = [{b: 0 for b in 'ACGT'} for _ in range(majority_width)]
    for seq in seqs:
        for i, base in enumerate(seq.upper()):
            if base in 'ACGT':
                counts[i][base] += 1

    PC = 0.5   # pseudocount
    pwm       = {}
    consensus = []
    ic_total  = 0.0

    for i, col in enumerate(counts):
        n     = sum(col.values()) + 4 * PC
        freqs = {b: (col[b] + PC) / n for b in 'ACGT'}
        pwm[i] = freqs

        ic = 2.0 + sum(f * math.log2(f) if f > 0 else 0 for f in freqs.values())
        ic_total += ic

        top = sorted(freqs.items(), key=lambda x: -x[1])
        if top[0][1] >= 0.5 or (top[0][1] - top[1][1]) >= 0.15:
            consensus.append(top[0][0])
        else:
            consensus.append('N')

    return pwm, ''.join(consensus), round(ic_total, 3)


# ---------------------------------------------------------------------------
# 6. Evidence level classification
# ---------------------------------------------------------------------------
def classify_evidence(chip_hits: list) -> str:
    if not chip_hits:
        return 'PREDICTED'
    strains = [h['strain_group'] for h in chip_hits]
    if 'ATCC13032' in strains:
        return 'EXPERIMENTAL_ATCC13032'
    return 'EXPERIMENTAL_OTHER_STRAIN'


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    print("Loading data sources...")
    records             = parse_fasta(FASTA_PATH)
    gene_map            = load_gene_name_map(GENE_MAP)
    tss                 = load_tss(TSS_PATH)
    chip_pairs, tf_name_map = load_chipseq(CHIPSEQ_PATH)
    print(f"  FASTA records     : {len(records)}")
    print(f"  TSS genes         : {len(tss)}")
    print(f"  ChIP edge pairs   : {len(chip_pairs)}")

    # ── Enrich records ──────────────────────────────────────────────────────
    enriched = []
    for rec in records:
        tg_l = rec['tg_locus']

        # Gene name
        rec['tg_gene_name'] = gene_map.get(tg_l, rec.get('tg_gene', tg_l))

        # TSS alignment
        tss_info = tss.get(tg_l, {})
        rec['tss_position'] = tss_info.get('tss_position', '')
        rec['strand']       = tss_info.get('strand', '')
        rec['tss_type']     = tss_info.get('tss_type', '')
        rec['promoter_70bp'] = tss_info.get('promoter_70bp_upstream', '')

        # Estimated absolute genomic position of binding site
        tss_pos = tss_info.get('tss_position')
        if tss_pos:
            delta = rec['pos_bp']
            strand = rec['strand']
            rec['genomic_pos_est'] = int(tss_pos) + (delta if strand == 'Fwd' else -delta)
        else:
            rec['genomic_pos_est'] = ''

        # ChIP support lookup
        tf_l = rec['tf_locus'].lower()
        tf_n = rec['tf_name'].lower()
        chip_hits = (
            chip_pairs.get((tf_l, tg_l)) or
            chip_pairs.get((tf_name_map.get(tf_n, ''), tg_l)) or
            []
        )

        rec['has_chip_support']   = bool(chip_hits)
        rec['chip_evidence_type'] = '; '.join(sorted({h['evidence'] for h in chip_hits}))
        rec['chip_source']        = '; '.join(sorted({h['source'] for h in chip_hits}))
        rec['chip_pmid']          = '; '.join(sorted({h['pmid'] for h in chip_hits if h['pmid']}))
        rec['chip_strain']        = '; '.join(sorted({h['strain_group'] for h in chip_hits}))
        rec['evidence_level']     = classify_evidence(chip_hits)
        rec['source_accession']   = 'RegPrecise_3.2_ATCC13032'

        enriched.append(rec)

    # ── Write binding sites TSV ──────────────────────────────────────────────
    TSV_COLS = [
        'tf_name', 'tf_locus', 'tg_locus', 'tg_gene_name',
        'score', 'pos_bp', 'genomic_pos_est', 'strand', 'tss_type', 'tss_position',
        'sequence', 'site_width',
        'has_chip_support', 'chip_evidence_type', 'chip_source', 'chip_pmid', 'chip_strain',
        'evidence_level', 'source_accession', 'genome',
    ]
    with open(OUT_TSV, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=TSV_COLS, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(enriched)
    print(f"\nWrote binding sites: {OUT_TSV}  ({len(enriched)} rows)")

    # ── Compute per-TF PWMs ──────────────────────────────────────────────────
    print("\nComputing PWMs per TF:")
    tf_groups = defaultdict(list)
    for rec in enriched:
        tf_groups[rec['tf_name']].append(rec)

    pwm_output = {}
    for tf in sorted(tf_groups):
        group = tf_groups[tf]
        seqs  = [r['sequence'] for r in group]
        chip_n = sum(1 for r in group if r['has_chip_support'])
        exp_n  = sum(1 for r in group if r['evidence_level'] == 'EXPERIMENTAL_ATCC13032')

        pwm, consensus, ic = compute_pwm(seqs)
        pwm_output[tf] = {
            'tf_name':                tf,
            'tf_locus':               group[0]['tf_locus'],
            'n_sites':                len(seqs),
            'site_width':             len(seqs[0]),
            'consensus':              consensus,
            'ic_bits_total':          ic,
            'chip_supported_sites':   chip_n,
            'experimental_atcc_sites': exp_n,
            'source_accession':       'RegPrecise_3.2_ATCC13032',
            'pwm': {str(i): {b: round(v, 4) for b, v in col.items()}
                    for i, col in pwm.items()},
        }
        chip_tag = f" [ChIP={chip_n}/{len(seqs)}]" if chip_n else ""
        print(f"  {tf:15s}: {len(seqs):3d} sites | {consensus} | IC={ic:5.1f}{chip_tag}")

    with open(OUT_PWM, 'w', encoding='utf-8') as f:
        json.dump(pwm_output, f, indent=2, ensure_ascii=False)
    print(f"\nWrote PWMs: {OUT_PWM}  ({len(pwm_output)} TFs)")

    # ── Patch regprecise_regulations.csv with best binding site per TF-TG ───
    # Best = highest Score
    best_site = {}
    for rec in enriched:
        key = (rec['tf_name'].lower(), rec['tg_locus'])
        if key not in best_site or rec['score'] > best_site[key]['score']:
            best_site[key] = rec

    with open(OUT_REGS, encoding='utf-8-sig') as f:
        regs_rows = list(csv.DictReader(f))

    # Ensure columns exist
    fieldnames = list(regs_rows[0].keys()) if regs_rows else []
    for col in ('Binding_site', 'evidence_level', 'source_accession'):
        if col not in fieldnames:
            fieldnames.append(col)

    updated = 0
    for row in regs_rows:
        tf_n = (row.get('TF_name') or '').lower()
        tg_l = (row.get('TG_locusTag') or '').lower()
        site = best_site.get((tf_n, tg_l))
        if site:
            row['Binding_site']     = site['sequence']
            row['evidence_level']   = site['evidence_level']
            row['source_accession'] = site['source_accession']
            updated += 1
        else:
            row.setdefault('Binding_site', '')
            row.setdefault('evidence_level', 'PREDICTED')
            row.setdefault('source_accession', 'RegPrecise_3.2_ATCC13032')

    with open(OUT_REGS, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        w.writeheader()
        w.writerows(regs_rows)
    print(f"\nUpdated {updated}/{len(regs_rows)} rows in {OUT_REGS.name}")

    # ── Final summary ────────────────────────────────────────────────────────
    ev_counts = defaultdict(int)
    for r in enriched:
        ev_counts[r['evidence_level']] += 1

    tss_aligned = sum(1 for r in enriched if r['tss_position'])
    gpos_est    = sum(1 for r in enriched if r['genomic_pos_est'])

    print("\n=== Final Summary ===")
    print(f"  Total TFBS records         : {len(enriched)}")
    print(f"  TSS-aligned records        : {tss_aligned}")
    print(f"  Genomic position estimated : {gpos_est}")
    for ev, cnt in sorted(ev_counts.items()):
        print(f"  {ev:35s}: {cnt}")
    print(f"  TFs with PWM computed      : {len(pwm_output)}")


if __name__ == '__main__':
    main()
