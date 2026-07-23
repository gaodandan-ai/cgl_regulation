#!/usr/bin/env python3
"""
analyze_novel_discoveries.py
=============================
新发现 1: 孤儿 iModulon 调控因子预测
新发现 2: 代谢通路"调控盲区"识别 + 工程靶点预测

Approach 1 – Orphan iModulon TF prediction:
  Use expression correlation between known TFs and genes in each
  uncharacterized iModulon to predict candidate regulators.
  Logic: if TF X co-varies with 70% of genes in an orphan module,
  X is likely the undiscovered regulator of that module.

Approach 2 – Metabolic regulation gap:
  Map regulatory network onto metabolic reaction graph.
  Identify amino acid biosynthesis steps with ZERO known TF regulator
  (regulatory orphan reactions) → priority targets for strain engineering.
"""

import json, os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from collections import defaultdict

warnings.filterwarnings('ignore')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data', 'reference')
OUT  = os.path.join(ROOT, 'analysis_output', 'novel_discoveries')
os.makedirs(OUT, exist_ok=True)

print("=" * 65)
print("  Novel Discovery Analysis – C. glutamicum")
print("=" * 65)

# ─── load data ───────────────────────────────────────────────────────────────
print("\n[Loading data]")

with open(os.path.join(DATA, 'imodulon', 'imodulon_gene_weights.json')) as f:
    imod_raw = json.load(f)
with open(os.path.join(DATA, 'imodulon', 'imodulon_metadata.json')) as f:
    imod_meta = json.load(f)

reg  = pd.read_csv(os.path.join(DATA, 'regulations.csv'))
chip = pd.read_csv(os.path.join(DATA, 'chipseq_regulations.csv'))
corr = pd.read_csv(os.path.join(DATA, 'expression_compendium',
                                'tf_target_compendium_correlations.csv'))
ec   = pd.read_csv(os.path.join(DATA, 'edge_confidence', 'tf_gene_edge_scores.csv'))
gmap = pd.read_csv(os.path.join(DATA, 'gene_mapping.csv'))
met  = pd.read_csv(os.path.join(DATA, 'metabolic_models', 'gene_reaction_mapping.csv'))

name_map = {}
for _, r in gmap.iterrows():
    nm = r['gene_name'] if pd.notna(r['gene_name']) else r['cgl_locus']
    name_map[r['cgl_locus']] = nm
    if pd.notna(r['cg_locus']):
        name_map[r['cg_locus']] = nm
def gn(l): return name_map.get(l, l)

# All curated + ChIP TF -> target relationships
all_edges = pd.concat([
    reg[['TF_locusTag','TG_locusTag']].rename(columns={'TF_locusTag':'tf','TG_locusTag':'tg'}),
    chip[['TF_locusTag','TG_locusTag']].rename(columns={'TF_locusTag':'tf','TG_locusTag':'tg'})
]).drop_duplicates()

# TF -> target set
tf_regulon = defaultdict(set)
for _, r in all_edges.iterrows():
    tf_regulon[r['tf']].add(r['tg'])

all_tfs = set(all_edges['tf'])
print(f"  {len(imod_raw)} iModulons loaded")
print(f"  {len(corr)} expression correlation pairs")
print(f"  {len(met)} gene-reaction mapping rows")
print(f"  {len(all_tfs)} TFs with regulatory data")


# ═══════════════════════════════════════════════════════════════════════════
# DISCOVERY 1: Orphan iModulon TF prediction
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*65)
print("  DISCOVERY 1: Orphan iModulon → Candidate TF Predictor")
print("="*65)

# Collect orphan iModulons with enough genes to be meaningful
orphan_imods = []
for imod_key, imod_data in imod_raw.items():
    if imod_data.get('linked_regulator') is None and imod_data.get('category') == 'uncharacterized':
        genes = list(imod_data.get('genes', {}).keys())
        genes = [g for g in genes if g and isinstance(g, str)]
        ve = imod_data.get('variance_explained', 0)
        if len(genes) >= 5:   # only modules with meaningful gene count
            orphan_imods.append({
                'imod_key': imod_key,
                'name': imod_data.get('name'),
                'genes': genes,
                'gene_count': len(genes),
                'variance_explained': ve,
                'pathways': [p['pathway_name'] for p in imod_data.get('enriched_pathways', [])[:3]],
            })

orphan_imods.sort(key=lambda x: -x['variance_explained'])
print(f"\n  Uncharacterized iModulons with ≥5 genes: {len(orphan_imods)}")

# For each orphan iModulon, find the best candidate TF
# Method: expression correlation overlap score
# Score = fraction of orphan iModulon genes that correlate with TF (|r|≥0.3)

# Build TF -> correlated gene set from expression data
TF_CORR_THRESHOLD = 0.30
tf_corr_genes = defaultdict(set)
for _, r in corr.iterrows():
    if abs(r['correlation']) >= TF_CORR_THRESHOLD:
        tf_corr_genes[r['tf']].add(r['target'])

# Also build reverse: target -> correlated TFs
target_corr_tfs = defaultdict(set)
for _, r in corr.iterrows():
    if abs(r['correlation']) >= TF_CORR_THRESHOLD:
        target_corr_tfs[r['target']].add(r['tf'])

predictions = []
for imod in orphan_imods:
    imod_genes = set(imod['genes'])

    # Score each TF by how many orphan genes it correlates with
    tf_scores = {}
    for tf in all_tfs | set(corr['tf']):
        corr_with_imod = imod_genes & tf_corr_genes.get(tf, set())
        if len(corr_with_imod) == 0:
            continue
        overlap_frac  = len(corr_with_imod) / len(imod_genes)
        # Bonus: if TF already regulates some iModulon genes (curated)
        curated_overlap = imod_genes & tf_regulon.get(tf, set())
        curated_bonus   = len(curated_overlap) / max(len(imod_genes), 1) * 0.5
        score = overlap_frac + curated_bonus
        tf_scores[tf] = {
            'score': round(score, 4),
            'n_corr': len(corr_with_imod),
            'pct_corr': round(len(corr_with_imod)/len(imod_genes)*100, 1),
            'n_curated': len(curated_overlap),
            'corr_genes': list(corr_with_imod)[:5],
        }

    if not tf_scores:
        continue

    top3 = sorted(tf_scores.items(), key=lambda x: -x[1]['score'])[:3]
    predictions.append({
        'imod': imod,
        'top_candidates': top3,
    })

print(f"\n  {'iModulon':<25} {'Genes':>6} {'VarExp%':>8}  {'Top Candidate TF':<14} {'Score':>7}  {'% correlated':>13}  {'Curated support'}")
print("  " + "-" * 90)
for p in predictions:
    im = p['imod']
    if not p['top_candidates']:
        continue
    tf, stats_d = p['top_candidates'][0]
    tf_nm = gn(tf)
    paths = ', '.join(im['pathways'][:2])
    print(f"  {im['name']:<25} {im['gene_count']:>6} {im['variance_explained']*100:>7.2f}%  "
          f"{tf_nm:<14} {stats_d['score']:>7.4f}  "
          f"{stats_d['pct_corr']:>10.1f}%   "
          f"{stats_d['n_curated']} known edges")


# Best prediction: uncharacterized-10 (highest var explained + most genes)
best = predictions[0]
best_imod = best['imod']
print(f"\n  >>> BEST PREDICTION <<<")
print(f"  iModulon: {best_imod['name']}  ({best_imod['gene_count']} genes, {best_imod['variance_explained']*100:.2f}% variance)")
print(f"  Enriched pathways: {', '.join(best_imod['pathways'])}")
print()
print(f"  {'Rank':<5} {'TF Locus':<12} {'TF Name':<14} {'Score':>7} {'%corr':>7} {'Curated':>8}")
print("  " + "-" * 55)
for i, (tf, s) in enumerate(best['top_candidates'], 1):
    print(f"  {i:<5} {tf:<12} {gn(tf):<14} {s['score']:>7.4f} {s['pct_corr']:>6.1f}% {s['n_curated']:>8}")
    if i == 1:
        print(f"       Example correlated genes: {', '.join(gn(g) for g in s['corr_genes'])}")


# ═══════════════════════════════════════════════════════════════════════════
# DISCOVERY 2: Metabolic regulatory gaps
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*65)
print("  DISCOVERY 2: Metabolic 'Regulatory Blind Spots'")
print("="*65)

print("\n  Metabolic model columns:", met.columns.tolist())
print(met.head(3).to_string())

# Map gene -> reactions
gene_col = next((c for c in met.columns if c.lower() == 'gene'), None)
rxn_col  = next((c for c in met.columns if 'reaction' in c.lower() and 'name' not in c.lower()), None)
# Prefer human-readable pathway_name over pathway_id
path_col = 'pathway_name' if 'pathway_name' in met.columns else \
           next((c for c in met.columns if 'pathway' in c.lower() or 'subsystem' in c.lower()), None)

print(f"\n  Using columns: gene={gene_col}, reaction={rxn_col}, pathway={path_col}")

# All known regulated genes (has at least 1 TF in curated or ChIP)
tg_with_regulator = set(all_edges['tg'])

# For each gene in metabolic model, check if it has a known regulator
met_clean = met.dropna(subset=[gene_col])
gene_has_reg = {}
for g in met_clean[gene_col].unique():
    gene_has_reg[g] = g in tg_with_regulator

total_met_genes = met_clean[gene_col].nunique()
unregulated_met = sum(1 for v in gene_has_reg.values() if not v)
regulated_met   = sum(1 for v in gene_has_reg.values() if v)
print(f"\n  Total metabolic genes: {total_met_genes}")
print(f"  With known TF regulator: {regulated_met} ({regulated_met/total_met_genes*100:.1f}%)")
print(f"  Without any TF regulator: {unregulated_met} ({unregulated_met/total_met_genes*100:.1f}%)")

# Per-pathway regulatory coverage
if path_col:
    pathway_stats = []
    for pathway, grp in met_clean.groupby(path_col):
        genes_in_path = grp[gene_col].unique()
        n_reg   = sum(1 for g in genes_in_path if g in tg_with_regulator)
        n_unreg = len(genes_in_path) - n_reg
        pathway_stats.append({
            'pathway': pathway,
            'n_genes': len(genes_in_path),
            'n_regulated': n_reg,
            'n_unregulated': n_unreg,
            'pct_regulated': n_reg / len(genes_in_path) * 100 if len(genes_in_path) > 0 else 0,
        })
    path_df = pd.DataFrame(pathway_stats).sort_values('pct_regulated')
    path_df.to_csv(os.path.join(OUT, 'pathway_regulatory_coverage.csv'), index=False)

    print(f"\n  === Pathways most in need of regulatory discovery ===")
    print(f"  (sorted by % unregulated, ≥5 genes)")
    print(f"  {'Pathway':<45} {'Genes':>6} {'Reg':>5} {'Unreg':>6} {'Reg%':>7}")
    print("  " + "-" * 75)
    target_paths = ['amino acid', 'lysine', 'glutamate', 'nitrogen', 'carbon']
    aa_paths = path_df[(path_df['n_genes']>=5)]
    # highlight amino acid pathways
    for _, row in aa_paths.head(20).iterrows():
        flag = ''
        for tp in target_paths:
            if tp in str(row['pathway']).lower():
                flag = ' ***'
                break
        print(f"  {str(row['pathway']):<45} {int(row['n_genes']):>6} "
              f"{int(row['n_regulated']):>5} {int(row['n_unregulated']):>6} "
              f"{row['pct_regulated']:>6.1f}%{flag}")

    # Find which TF best covers each unregulated pathway
    print(f"\n  === Best TF candidate for each under-regulated pathway ===")
    aa_related = path_df[path_df['pathway'].apply(
        lambda x: any(t in str(x).lower() for t in ['amino','lysin','glutam','nitr','carbon','amin'])
    ) & (path_df['n_genes']>=5)]

    for _, row in aa_related.iterrows():
        path_genes = set(met_clean[met_clean[path_col]==row['pathway']][gene_col])
        unreg_genes = path_genes - tg_with_regulator
        if len(unreg_genes) == 0:
            continue
        # Find TF that regulates most OTHER genes in this pathway (best known regulator)
        best_tf = None
        best_cnt = 0
        for tf, targets in tf_regulon.items():
            cnt = len(targets & path_genes)
            if cnt > best_cnt:
                best_cnt = cnt
                best_tf = tf
        unreg_names = ', '.join(gn(g) for g in list(unreg_genes)[:5])
        print(f"\n  {row['pathway']}")
        print(f"    Unregulated genes ({len(unreg_genes)}): {unreg_names}...")
        if best_tf:
            print(f"    Best known regulator: {gn(best_tf)} (covers {best_cnt}/{len(path_genes)} pathway genes)")
            print(f"    => {len(unreg_genes)} genes have NO known TF — candidate for extension study")


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE
# ═══════════════════════════════════════════════════════════════════════════
print("\n[Generating figures...]")

fig = plt.figure(figsize=(18, 12))
fig.suptitle('Novel Regulatory Discoveries in C. glutamicum\n'
             'Orphan iModulon TF Prediction & Metabolic Blind Spots',
             fontsize=14, fontweight='bold')

# ── Panel A: iModulon overview (size vs variance, colored by status) ────────
ax_a = fig.add_subplot(2, 3, 1)
for im_key, im_data in imod_raw.items():
    genes = list(im_data.get('genes', {}).keys())
    genes = [g for g in genes if g and isinstance(g, str)]
    ve = im_data.get('variance_explained', 0) * 100
    n = len(genes)
    linked = im_data.get('linked_regulator')
    cat = im_data.get('category', '')
    if 'uncharacterized' in cat.lower():
        color, zorder = '#ef4444', 3
    elif linked:
        color, zorder = '#10b981', 2
    else:
        color, zorder = '#94a3b8', 1
    ax_a.scatter(n, ve, c=color, alpha=0.75, s=60, zorder=zorder)

leg_patches = [
    mpatches.Patch(color='#10b981', label='Linked to known TF'),
    mpatches.Patch(color='#ef4444', label='Uncharacterized (orphan)'),
    mpatches.Patch(color='#94a3b8', label='Other (no TF assigned)'),
]
# Annotate top orphans
for p in predictions[:4]:
    im = p['imod']
    n = im['gene_count']
    ve = im['variance_explained'] * 100
    ax_a.annotate(im['name'].replace('uncharacterized-','U-'),
                  (n, ve), fontsize=7.5, ha='left',
                  xytext=(5, 2), textcoords='offset points')
ax_a.set_xlabel('Number of Genes in iModulon', fontsize=10)
ax_a.set_ylabel('Variance Explained (%)', fontsize=10)
ax_a.set_title('(A) iModulon Landscape\n(red = no known regulator)', fontsize=11)
ax_a.legend(handles=leg_patches, fontsize=8, loc='upper right')
ax_a.grid(alpha=0.3)

# ── Panel B: Top candidate TFs for each orphan iModulon ────────────────────
ax_b = fig.add_subplot(2, 3, 2)
valid_preds = [(p, p['top_candidates'][0]) for p in predictions if p['top_candidates']][:8]
y_pos   = range(len(valid_preds))
pct_vals = [s['pct_corr'] for _, (_, s) in valid_preds]
labels  = [p['imod']['name'].replace('uncharacterized-','U-') + '\n→ ' + gn(tf)
           for p, (tf, s) in valid_preds]
colors_b = ['#ef4444' if s['pct_corr']>=30 else '#f59e0b' if s['pct_corr']>=15 else '#94a3b8'
            for _, (_, s) in valid_preds]
ax_b.barh(list(y_pos), pct_vals, color=colors_b, alpha=0.85)
ax_b.set_yticks(list(y_pos))
ax_b.set_yticklabels(labels, fontsize=8)
ax_b.invert_yaxis()
ax_b.set_xlabel('% of iModulon Genes Co-varying\nwith Candidate TF', fontsize=10)
ax_b.set_title('(B) Predicted TF Regulators\nfor Orphan iModulons', fontsize=11)
ax_b.axvline(20, color='#64748b', ls='--', alpha=0.5, label='20% threshold')
ax_b.legend(fontsize=8)
ax_b.grid(axis='x', alpha=0.3)

# ── Panel C: Best prediction detail (uncharacterized-10 if available) ──────
ax_c = fig.add_subplot(2, 3, 3)
if predictions:
    best_p = predictions[0]
    best_im = best_p['imod']
    tfs_plot = best_p['top_candidates'][:6]
    tf_names = [gn(tf) for tf, _ in tfs_plot]
    tf_scores_plt = [s['score'] for _, s in tfs_plot]
    tf_pcts   = [s['pct_corr'] for _, s in tfs_plot]
    x = range(len(tfs_plot))
    bars = ax_c.bar(x, tf_scores_plt, color='#6366f1', alpha=0.85)
    ax_c.set_xticks(list(x))
    ax_c.set_xticklabels(tf_names, rotation=30, ha='right', fontsize=9)
    ax_c.set_ylabel('Candidate Score', fontsize=10)
    ax_c.set_title(f'(C) Top TF Candidates for\n{best_im["name"]}'
                   f' ({best_im["gene_count"]} genes)', fontsize=11)
    for bar, pct in zip(bars, tf_pcts):
        ax_c.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.002,
                  f'{pct:.0f}%', ha='center', fontsize=8)
    ax_c.grid(axis='y', alpha=0.3)

# ── Panel D: Metabolic regulatory coverage per pathway ─────────────────────
if path_col and 'path_df' in dir():
    ax_d = fig.add_subplot(2, 1, 2)
    # Show top 20 pathways by gene count, colored by regulation %
    top_paths = path_df[path_df['n_genes']>=5].sort_values('n_genes', ascending=False).head(20)
    x = range(len(top_paths))

    # Stacked bar: regulated (green) + unregulated (red)
    ax_d.bar(list(x), top_paths['n_regulated'],   label='Has known TF regulator', color='#10b981', alpha=0.85)
    ax_d.bar(list(x), top_paths['n_unregulated'], bottom=top_paths['n_regulated'],
             label='NO known TF regulator', color='#ef4444', alpha=0.85)

    ax_d.set_xticks(list(x))
    ax_d.set_xticklabels(
        [str(p)[:30] for p in top_paths['pathway']],
        rotation=40, ha='right', fontsize=8
    )
    ax_d.set_ylabel('Number of Metabolic Genes', fontsize=11)
    ax_d.set_title('(D) Metabolic Pathway Regulatory Coverage\n'
                   '(red = genes with no known TF regulator = potential discovery targets)',
                   fontsize=12)
    ax_d.legend(fontsize=10)
    ax_d.grid(axis='y', alpha=0.3)

plt.tight_layout(rect=[0, 0, 1, 0.96])
p_out = os.path.join(OUT, 'fig_novel_discoveries.png')
plt.savefig(p_out, dpi=200, bbox_inches='tight')
plt.close()
print(f"  Saved: {p_out}")

# Save prediction table
pred_rows = []
for p in predictions:
    for rank, (tf, s) in enumerate(p['top_candidates'], 1):
        pred_rows.append({
            'imodulon': p['imod']['name'],
            'gene_count': p['imod']['gene_count'],
            'variance_explained': p['imod']['variance_explained'],
            'candidate_rank': rank,
            'candidate_tf': tf,
            'candidate_tf_name': gn(tf),
            'score': s['score'],
            'pct_genes_correlated': s['pct_corr'],
            'n_curated_edges': s['n_curated'],
            'corr_gene_examples': ', '.join(gn(g) for g in s['corr_genes']),
        })
pd.DataFrame(pred_rows).to_csv(os.path.join(OUT, 'orphan_imodulon_predictions.csv'), index=False)

print(f"\n  Saved: orphan_imodulon_predictions.csv")
print(f"  Saved: pathway_regulatory_coverage.csv")
print("\n  Analysis complete!")
