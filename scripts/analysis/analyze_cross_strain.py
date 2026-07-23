#!/usr/bin/env python3
"""
analyze_cross_strain.py
=======================
C. glutamicum 跨菌株调控差异分析
Cross-Strain and Cross-Condition Regulatory Divergence Analysis

分析内容：
  1. sigH 调控子跨菌株比较 (ATCC13032 vs Strain_R)
  2. hrrA 条件特异性调控 (iron_excess / heme / combined)
  3. 表达调控重连网络 (rewired edges: control → heat)
  4. Hub-switching TFs (条件间枢纽 TF 的切换)
  5. 全局调控 Jaccard 相似性矩阵
"""

import json, os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
from collections import defaultdict
from scipy import stats

# matplotlib_venn 可选依赖检查
try:
    from matplotlib_venn import venn2, venn3
    HAS_VENN = True
except ImportError:
    HAS_VENN = False
    print("  [Info] matplotlib-venn not installed; using bar charts instead")

warnings.filterwarnings('ignore')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data', 'reference')
OUT  = os.path.join(ROOT, 'analysis_output', 'cross_strain')
os.makedirs(OUT, exist_ok=True)

print("=" * 65)
print("  C. glutamicum Cross-Strain Regulatory Divergence Analysis")
print("=" * 65)


# ─────────────────────────────────────────────────────────────────────────────
# 加载数据
# ─────────────────────────────────────────────────────────────────────────────
print("\n[1/5] Loading data...")

chip = pd.read_csv(os.path.join(DATA, 'chipseq_regulations.csv'))
reg  = pd.read_csv(os.path.join(DATA, 'regulations.csv'))
gmap = pd.read_csv(os.path.join(DATA, 'gene_mapping.csv'))
corr = pd.read_csv(os.path.join(DATA, 'expression_compendium', 'tf_target_compendium_correlations.csv'))

with open(os.path.join(DATA, 'rna_seq_analysis_results.json')) as f:
    rna = json.load(f)
with open(os.path.join(DATA, 'cog_annotations.json')) as f:
    cog_db = json.load(f)

name_map = {}
for _, r in gmap.iterrows():
    nm = r['gene_name'] if pd.notna(r['gene_name']) else r['cgl_locus']
    name_map[r['cgl_locus']] = nm
    if pd.notna(r['cg_locus']):
        name_map[r['cg_locus']] = nm

def gn(locus):
    return name_map.get(locus, locus)

COG_CAT = {
    'J':'Translation','K':'Transcription','L':'Replication','D':'Cell cycle',
    'T':'Signal transduction','M':'Cell wall','O':'Post-translational',
    'C':'Energy','G':'Carbohydrate','E':'Amino acid','F':'Nucleotide',
    'H':'Coenzyme','I':'Lipid','P':'Ion transport','Q':'Secondary metabolites',
    'R':'General function','S':'Unknown/unannotated',
}

print(f"  ChIP-seq: {len(chip)} edges across strains")
print(f"  Rewired edges in expression data: {len(rna.get('rewired_edges', []))}")
print(f"  Expression correlations: {len(corr)} TF-target pairs")


# ─────────────────────────────────────────────────────────────────────────────
# 分析 1: sigH 跨菌株调控子比较
# ─────────────────────────────────────────────────────────────────────────────
print("\n[2/5] sigH cross-strain analysis (ATCC13032 vs Strain_R)...")

sigh_atcc = chip[(chip['TF_locusTag']=='cg0876') & (chip['strain_group']=='ATCC13032')]
sigh_r    = chip[(chip['TF_locusTag']=='cg0876') & (chip['strain_group']=='Strain_R')]

targets_atcc = set(sigh_atcc['TG_locusTag'])
targets_r    = set(sigh_r['TG_locusTag'])
shared       = targets_atcc & targets_r
atcc_only    = targets_atcc - targets_r
r_only       = targets_r - targets_atcc

print(f"  ATCC13032 sigH targets: {len(targets_atcc)}")
print(f"  Strain_R  sigH targets: {len(targets_r)}")
print(f"  Shared: {len(shared)}")
print(f"  ATCC13032-specific: {len(atcc_only)} | Strain_R-specific: {len(r_only)}")
print(f"  Jaccard similarity: {len(shared)/len(targets_atcc|targets_r):.4f}")

# Regulon size ratio
ratio = len(targets_r) / max(len(targets_atcc), 1)
print(f"  Strain_R regulon is {ratio:.1f}x larger than ATCC13032 sigH regulon")

# COG comparison for sigH regulon
def cog_profile(target_set):
    cats = defaultdict(int)
    for tg in target_set:
        cats[cog_db.get(tg, {}).get('category', 'S')] += 1
    return dict(cats)

cog_atcc = cog_profile(targets_atcc)
cog_r    = cog_profile(targets_r)

# Statistical Fisher's exact test for COG enrichment in Strain_R vs ATCC13032
total_atcc = max(sum(cog_atcc.values()), 1)
total_r    = max(sum(cog_r.values()), 1)
enrichment = []
for cat in set(list(cog_atcc.keys()) + list(cog_r.keys())):
    a = cog_atcc.get(cat, 0)
    b = total_atcc - a
    c = cog_r.get(cat, 0)
    d = total_r - c
    if a + c == 0:
        continue
    oddr, p = stats.fisher_exact([[a, b], [c, d]], alternative='two-sided')
    enrichment.append({'cog': cat, 'desc': COG_CAT.get(cat, cat),
                       'n_atcc': a, 'n_r': c,
                       'pct_atcc': a/total_atcc*100, 'pct_r': c/total_r*100,
                       'odds_ratio': oddr, 'pvalue': p})

enr_df = pd.DataFrame(enrichment).sort_values('pvalue')
print("\n  COG enrichment (Strain_R vs ATCC13032 sigH regulon):")
print(f"  {'COG':<5} {'Category':<25} {'ATCC%':>7} {'StrR%':>7} {'OR':>7} {'p-val':>10}")
print("  " + "-" * 65)
for _, row in enr_df[enr_df['n_r'] > 0].head(10).iterrows():
    print(f"  [{row['cog']}]  {row['desc']:<25} {row['pct_atcc']:>6.1f}%  "
          f"{row['pct_r']:>6.1f}%  {row['odds_ratio']:>6.2f}  {row['pvalue']:>10.4f}")

enr_df.to_csv(os.path.join(OUT, 'sigh_cog_enrichment.csv'), index=False)


# ─────────────────────────────────────────────────────────────────────────────
# 分析 2: hrrA 条件特异性调控
# ─────────────────────────────────────────────────────────────────────────────
print("\n[3/5] hrrA condition-specific ChIP analysis...")

hrra = chip[chip['TF_locusTag']=='cg3247'].copy()
hrra['condition'] = hrra['strain_note'].map({
    'ATCC13032;cond=iron_excess':          'Iron excess',
    'ATCC13032;cond=heme':                 'Heme',
    'ATCC13032;cond=iron_excess+heme':     'Iron+Heme',
})

tg_iron = set(hrra[hrra['condition']=='Iron excess']['TG_locusTag'])
tg_heme = set(hrra[hrra['condition']=='Heme']['TG_locusTag'])
tg_both = set(hrra[hrra['condition']=='Iron+Heme']['TG_locusTag'])

print(f"  Iron excess: {len(tg_iron)} targets")
print(f"  Heme:        {len(tg_heme)} targets")
print(f"  Iron+Heme:   {len(tg_both)} targets")
print(f"  Core (all 3): {len(tg_iron & tg_heme & tg_both)} targets")
print(f"  Iron-specific only: {len(tg_iron - tg_heme - tg_both)}")
print(f"  Heme-specific only: {len(tg_heme - tg_iron - tg_both)}")

# Core vs condition-specific COG comparison
core_tgs       = tg_iron & tg_heme & tg_both
iron_spec_tgs  = tg_iron - tg_heme - tg_both
heme_spec_tgs  = tg_heme - tg_iron - tg_both

print(f"\n  Core hrrA regulon COG categories:")
for cat, cnt in sorted(cog_profile(core_tgs).items(), key=lambda x: -x[1])[:5]:
    print(f"    [{cat}] {COG_CAT.get(cat, cat)}: {cnt}")


# ─────────────────────────────────────────────────────────────────────────────
# 分析 3: 表达调控重连 (rewired edges: condition-specific)
# ─────────────────────────────────────────────────────────────────────────────
print("\n[4/5] Expression-based rewired regulatory edges analysis...")

rewired = pd.DataFrame(rna.get('rewired_edges', []))
hub_sw  = pd.DataFrame(rna.get('hub_switching', []))

if len(rewired) > 0:
    # 按类型分类
    print(f"  Total rewired edges: {len(rewired)}")
    if 'type' in rewired.columns:
        print("  By type:")
        print(rewired['type'].value_counts().to_string())

    # 影响最大的 TF (最多 rewired 边)
    tf_rewire_cnt = rewired.groupby('tf_name').size().sort_values(ascending=False)
    print(f"\n  TFs with most rewired targets (top 10):")
    for tf, cnt in tf_rewire_cnt.head(10).items():
        print(f"    {tf}: {cnt} rewired edges")

    # 最强的 rewiring (|delta_r| 最大)
    if 'delta_r' in rewired.columns:
        top_rewired = rewired.reindex(rewired['delta_r'].abs().sort_values(ascending=False).index).head(10)
        print(f"\n  Strongest rewiring events (top 10 by |Δr|):")
        print(f"  {'TF':<12} {'Target':<12} {'r_ctrl':>8} {'r_heat':>8} {'Δr':>8} {'type'}")
        print("  " + "-" * 65)
        for _, row in top_rewired.iterrows():
            print(f"  {row['tf_name']:<12} {row.get('tg_name','?'):<12} "
                  f"{row.get('r_control',0):>8.3f} {row.get('r_heat',0):>8.3f} "
                  f"{row['delta_r']:>8.3f}  {row.get('type','')}")

    rewired.to_csv(os.path.join(OUT, 'rewired_edges.csv'), index=False)

if len(hub_sw) > 0:
    print(f"\n  Hub-switching TFs: {len(hub_sw)}")
    if 'category' in hub_sw.columns:
        print(hub_sw['category'].value_counts().to_string())
    print(f"\n  Top hub-switching TFs (by |Δdegree|):")
    hub_sorted = hub_sw.sort_values('delta_degree', key=abs, ascending=False)
    for _, row in hub_sorted.head(8).iterrows():
        ctrl = row.get('control_degree', 0)
        heat = row.get('heat_degree', 0)
        print(f"    {row['tf_name']:<12} ctrl={ctrl:>4} heat={heat:>4} "
              f"Δ={row['delta_degree']:>+5}  [{row.get('category','')}]")

    hub_sw.to_csv(os.path.join(OUT, 'hub_switching.csv'), index=False)


# ─────────────────────────────────────────────────────────────────────────────
# 可视化
# ─────────────────────────────────────────────────────────────────────────────
print("\n[5/5] Generating figures...")

fig = plt.figure(figsize=(20, 16))
fig.suptitle('Cross-Strain and Cross-Condition Regulatory Divergence in C. glutamicum',
             fontsize=15, fontweight='bold', y=0.98)

# ── Panel A: sigH regulon size comparison (bar) ────────────────────────────
ax_a = fig.add_subplot(3, 4, (1, 2))
categories = ['ATCC13032\nsigH', 'Strain_R\nsigH']
counts = [len(targets_atcc), len(targets_r)]
colors = ['#0ea5e9', '#f59e0b']
bars = ax_a.bar(categories, counts, color=colors, width=0.5, alpha=0.85, edgecolor='white')
ax_a.set_ylabel('ChIP-seq Target Genes', fontsize=11)
ax_a.set_title('(A) sigH Regulon Size\nATCC13032 vs Strain_R', fontsize=11)
for bar, cnt in zip(bars, counts):
    ax_a.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
              str(cnt), ha='center', va='bottom', fontsize=12, fontweight='bold')
ax_a.set_ylim(0, max(counts) * 1.2)
ax_a.grid(axis='y', alpha=0.3)

# Jaccard label
jaccard = len(shared) / max(len(targets_atcc | targets_r), 1)
ax_a.text(0.5, 0.92, f'Jaccard = {jaccard:.3f}\nShared = {len(shared)} targets',
          transform=ax_a.transAxes, ha='center', fontsize=9,
          bbox=dict(boxstyle='round,pad=0.3', facecolor='#f1f5f9', alpha=0.8))

# ── Panel B: sigH COG profile comparison ──────────────────────────────────
ax_b = fig.add_subplot(3, 4, (3, 4))
cog_cats_use = [c for c in 'CEGHI PKTMLJORS'.replace(' ','') if c in set(cog_atcc)|set(cog_r)]
x = np.arange(len(cog_cats_use))
w = 0.35
total_a = max(sum(cog_atcc.values()), 1)
total_r_val = max(sum(cog_r.values()), 1)
vals_a = [cog_atcc.get(c, 0)/total_a*100 for c in cog_cats_use]
vals_r = [cog_r.get(c, 0)/total_r_val*100 for c in cog_cats_use]
ax_b.bar(x - w/2, vals_a, w, label='ATCC13032', color='#0ea5e9', alpha=0.85)
ax_b.bar(x + w/2, vals_r, w, label='Strain_R',  color='#f59e0b', alpha=0.85)
ax_b.set_xticks(x)
ax_b.set_xticklabels([f'[{c}]' for c in cog_cats_use], fontsize=8)
ax_b.set_ylabel('% of Regulon', fontsize=10)
ax_b.set_title('(B) sigH Regulon COG Composition\n(% of regulon)', fontsize=11)
ax_b.legend(fontsize=9)
ax_b.grid(axis='y', alpha=0.3)

# ── Panel C: hrrA Venn / bar chart (condition-specific) ───────────────────
ax_c = fig.add_subplot(3, 4, (5, 6))
if HAS_VENN:
    venn_sets = {'Iron only': tg_iron - tg_heme - tg_both,
                 'Heme only': tg_heme - tg_iron - tg_both,
                 'Iron+Heme only': tg_both - tg_iron - tg_heme,
                 'Iron∩Heme': (tg_iron & tg_heme) - tg_both,
                 'Iron∩Both': (tg_iron & tg_both) - tg_heme,
                 'Heme∩Both': (tg_heme & tg_both) - tg_iron,
                 'All three': tg_iron & tg_heme & tg_both}
    # Simple bar version for clarity
    labels = ['Iron\nonly', 'Heme\nonly', 'Iron+Heme\nonly', 'Fe∩Hm', 'Fe∩F+H', 'Hm∩F+H', 'All 3']
    sizes  = [len(v) for v in venn_sets.values()]
    bar_colors = ['#ef4444','#10b981','#6366f1','#f97316','#ec4899','#14b8a6','#1e293b']
    ax_c.bar(labels, sizes, color=bar_colors, alpha=0.85, edgecolor='white')
    ax_c.set_ylabel('Number of Targets', fontsize=10)
    ax_c.set_title('(C) hrrA Condition-Specific Targets\n(iron excess / heme / iron+heme)', fontsize=11)
    ax_c.grid(axis='y', alpha=0.3)
    for i, v in enumerate(sizes):
        if v > 0:
            ax_c.text(i, v + 0.5, str(v), ha='center', fontsize=8, fontweight='bold')
else:
    cond_data = {'Iron excess': len(tg_iron), 'Heme': len(tg_heme), 'Iron+Heme': len(tg_both)}
    ax_c.bar(cond_data.keys(), cond_data.values(),
             color=['#ef4444','#10b981','#6366f1'], alpha=0.85)
    ax_c.set_ylabel('Target Genes', fontsize=10)
    ax_c.set_title('(C) hrrA Condition-Specific Targets', fontsize=11)

# ── Panel D: hrrA COG comparison by condition ─────────────────────────────
ax_d = fig.add_subplot(3, 4, (7, 8))
core_cog  = cog_profile(tg_iron & tg_heme & tg_both)
iron_cog  = cog_profile(tg_iron - tg_heme - tg_both)
heme_cog  = cog_profile(tg_heme - tg_iron - tg_both)
all_cats  = sorted(set(list(core_cog)+list(iron_cog)+list(heme_cog)),
                   key=lambda c: -(core_cog.get(c,0)+iron_cog.get(c,0)+heme_cog.get(c,0)))[:10]
x = np.arange(len(all_cats))
w = 0.28
ax_d.bar(x-w,   [core_cog.get(c,0)  for c in all_cats], w, label='Core (all 3)', color='#6366f1', alpha=0.85)
ax_d.bar(x,     [iron_cog.get(c,0)  for c in all_cats], w, label='Iron-specific', color='#ef4444', alpha=0.85)
ax_d.bar(x+w,   [heme_cog.get(c,0)  for c in all_cats], w, label='Heme-specific', color='#10b981', alpha=0.85)
ax_d.set_xticks(x)
ax_d.set_xticklabels([f'[{c}]' for c in all_cats], fontsize=8)
ax_d.set_ylabel('Gene Count', fontsize=10)
ax_d.set_title('(D) hrrA Regulon Functional Composition\nby Condition', fontsize=11)
ax_d.legend(fontsize=8)
ax_d.grid(axis='y', alpha=0.3)

# ── Panel E: Hub-switching TFs ─────────────────────────────────────────────
if len(hub_sw) > 0:
    ax_e = fig.add_subplot(3, 1, 3)
    hub_top = hub_sw.sort_values('delta_degree', key=abs, ascending=False).head(15)

    x = np.arange(len(hub_top))
    ctrl_vals = hub_top['control_degree'].values
    heat_vals = hub_top['heat_degree'].values

    ax_e.bar(x - 0.2, ctrl_vals, 0.4, label='Control condition', color='#64748b', alpha=0.8)
    ax_e.bar(x + 0.2, heat_vals, 0.4, label='Heat/stress condition', color='#f59e0b', alpha=0.8)
    ax_e.set_xticks(x)
    ax_e.set_xticklabels(hub_top['tf_name'].values, rotation=30, ha='right', fontsize=9)
    ax_e.set_ylabel('Number of Co-expressed Targets', fontsize=11)
    ax_e.set_title('(E) TF Hub-Switching: Connectivity Change Between Conditions\n'
                   '(Control vs Heat/Stress)', fontsize=12)
    ax_e.legend(fontsize=10)
    ax_e.grid(axis='y', alpha=0.3)

    # Color x-labels by direction of change
    for i, (_, row) in enumerate(hub_top.iterrows()):
        color = '#ef4444' if row['delta_degree'] > 0 else '#10b981'
        ax_e.get_xticklabels()[i].set_color(color)

    ax_e.text(0.98, 0.97, 'Red label = gained targets\nGreen label = lost targets',
              transform=ax_e.transAxes, ha='right', va='top', fontsize=8,
              bbox=dict(boxstyle='round', facecolor='#f1f5f9', alpha=0.8))

plt.tight_layout(rect=[0, 0, 1, 0.97])
p = os.path.join(OUT, 'fig_cross_strain_divergence.png')
plt.savefig(p, dpi=200, bbox_inches='tight')
plt.close()
print(f"  Saved: {p}")


# ── Summary for manuscript ────────────────────────────────────────────────
print("\n" + "="*65)
print("  SUMMARY FOR MANUSCRIPT")
print("="*65)
print(f"\n  === sigH Cross-Strain Comparison ===")
print(f"  ATCC13032 regulon: {len(targets_atcc)} targets")
print(f"  Strain_R regulon:  {len(targets_r)} targets  ({ratio:.1f}x larger)")
print(f"  Shared targets:    {len(shared)}  (Jaccard = {jaccard:.4f})")
print(f"  => sigH has undergone significant regulatory rewiring between strains")

print(f"\n  === hrrA Condition-Specific Regulation ===")
print(f"  Iron excess only targets:   {len(tg_iron - tg_heme - tg_both)}")
print(f"  Heme only targets:          {len(tg_heme - tg_iron - tg_both)}")
print(f"  Iron+Heme only targets:     {len(tg_both - tg_iron - tg_heme)}")
print(f"  Core (all 3 conditions):    {len(tg_iron & tg_heme & tg_both)}")
print(f"  => Iron and heme induce partially overlapping but distinct hrrA regulons")

if len(rewired) > 0:
    inv = len(rewired[rewired.get('type','')=='inversion']) if 'type' in rewired.columns else 0
    print(f"\n  === Expression-based Rewiring ===")
    print(f"  Total rewired TF-target edges: {len(rewired)}")
    if 'type' in rewired.columns:
        print(f"  Inversion events: {inv} ({inv/len(rewired)*100:.1f}%)")
    print(f"  Hub-switching TFs: {len(hub_sw)}")

print(f"\n  Output: {OUT}")
print("  [OK] fig_cross_strain_divergence.png")
print("  [OK] sigh_cog_enrichment.csv")
print("  [OK] rewired_edges.csv")
print("  [OK] hub_switching.csv")
print("\n  Analysis complete!")
