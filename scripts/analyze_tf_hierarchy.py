#!/usr/bin/env python3
"""
analyze_tf_hierarchy.py
=======================
C. glutamicum 主调控因子层级分析
Hierarchical Classification of Master Regulators in C. glutamicum

分析内容：
  1. 基于多维网络指标将 TF 分为 Tier 1/2/3 层级
  2. 每个 TF 调控子（regulon）的 COG 功能组成分析
  3. TF→TF 层级调控级联网络（谁调控谁）
  4. 代谢网络连通性（TF 对代谢通路的覆盖度）
  5. iModulon 交叉验证（转录组独立成分分析支持）
  6. 输出：排名表、层级图、regulon 功能热图
"""

import json, os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from collections import defaultdict

warnings.filterwarnings('ignore')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data', 'reference')
OUT  = os.path.join(ROOT, 'analysis_output', 'tf_hierarchy')
os.makedirs(OUT, exist_ok=True)

print("=" * 65)
print("  C. glutamicum TF Hierarchy Analysis")
print("=" * 65)


# ─────────────────────────────────────────────────────────────────────────────
# 1. 加载所有数据
# ─────────────────────────────────────────────────────────────────────────────
print("\n[1/7] Loading data...")

reg    = pd.read_csv(os.path.join(DATA, 'regulations.csv'))
chip   = pd.read_csv(os.path.join(DATA, 'chipseq_regulations.csv'))
gmap   = pd.read_csv(os.path.join(DATA, 'gene_mapping.csv'))
edges  = pd.read_csv(os.path.join(DATA, 'edge_confidence', 'tf_gene_edge_scores.csv'))

with open(os.path.join(DATA, 'network_centrality.json')) as f:
    centrality = json.load(f)
with open(os.path.join(DATA, 'cog_annotations.json')) as f:
    cog_db = json.load(f)
with open(os.path.join(DATA, 'imodulon', 'imodulon_metadata.json')) as f:
    imod_meta = json.load(f)

# Gene name lookup
name_map = {}
for _, r in gmap.iterrows():
    name_map[r['cgl_locus']] = r['gene_name'] if pd.notna(r['gene_name']) else r['cgl_locus']
    if pd.notna(r['cg_locus']):
        name_map[r['cg_locus']] = r['gene_name'] if pd.notna(r['gene_name']) else r['cg_locus']

def get_name(locus):
    return name_map.get(locus, locus)

print(f"  Loaded {len(reg)} curated regulations, {len(chip)} ChIP-seq regulations")
print(f"  {edges['tf_locus'].nunique()} TFs with confidence scores, {len(edges)} edge records")


# ─────────────────────────────────────────────────────────────────────────────
# 2. 构建 TF 综合指标矩阵
# ─────────────────────────────────────────────────────────────────────────────
print("\n[2/7] Building TF metric matrix...")

cent_nodes = centrality['nodes']
cent_df_rows = []
for locus, nd in cent_nodes.items():
    if nd.get('is_tf'):
        cent_df_rows.append({
            'locus': locus,
            'name': get_name(locus) if get_name(locus) != locus else nd.get('name', locus),
            'out_degree': nd.get('out_degree', 0),
            'in_degree': nd.get('in_degree', 0),
            'betweenness': nd.get('betweenness', 0),
            'pagerank': nd.get('pagerank', 0),
            'hub_score': nd.get('hub_score', 0),
            'authority_score': nd.get('auth_score', 0),
            'closeness': nd.get('closeness', 0),
            'importance': nd.get('importance', 0),
            'activation_ratio': nd.get('activation_ratio', np.nan),
            'is_sigma': nd.get('is_sigma', False),
            'n_activations': nd.get('n_activations', 0),
            'n_repressions': nd.get('n_repressions', 0),
        })

tf_df = pd.DataFrame(cent_df_rows)
print(f"  {len(tf_df)} TFs in centrality data")

# ChIP targets per TF
chip_counts = chip.groupby('TF_locusTag').size().rename('chipseq_targets')
tf_df = tf_df.merge(chip_counts, left_on='locus', right_index=True, how='left')
tf_df['chipseq_targets'] = tf_df['chipseq_targets'].fillna(0).astype(int)

# Curated targets per TF
reg_counts = reg.groupby('TF_locusTag').size().rename('curated_targets')
tf_df = tf_df.merge(reg_counts, left_on='locus', right_index=True, how='left')
tf_df['curated_targets'] = tf_df['curated_targets'].fillna(0).astype(int)
tf_df['chipseq_coverage'] = np.where(
    tf_df['curated_targets'] > 0,
    tf_df['chipseq_targets'] / tf_df['curated_targets'],
    0
)

# ── 来源多样性：计算每个 TF 的独立文献数和来源数（用于偏差修正）─────────────────
chip_diversity = chip.groupby('TF_locusTag').agg(
    n_chip_pmids  =('PMID',   'nunique'),
    n_chip_sources=('Source', 'nunique'),
).reset_index()
tf_df = tf_df.merge(chip_diversity, left_on='locus', right_on='TF_locusTag', how='left').drop(columns=['TF_locusTag'])
tf_df['n_chip_pmids']   = tf_df['n_chip_pmids'].fillna(0).astype(int)
tf_df['n_chip_sources'] = tf_df['n_chip_sources'].fillna(0).astype(int)

# diversity_weight: 单篇文章覆盖 >80% ChIP 靶点 → 0.5 惩罚
# 多来源（≥2 papers）→ 不惩罚
def _diversity_weight(row):
    if row['chipseq_targets'] == 0:
        return 1.0      # 无 ChIP 证据，不参与 ChIP 评分，不惩罚
    if row['n_chip_pmids'] <= 1 and row['chipseq_targets'] >= 100:
        return 0.5      # 大量靶点全来自单篇 → 明显偏差，降权
    if row['n_chip_pmids'] <= 1 and row['chipseq_targets'] >= 30:
        return 0.75     # 中等规模单篇 → 轻微降权
    return 1.0          # 多来源 or 小规模 → 不惩罚

tf_df['diversity_weight'] = tf_df.apply(_diversity_weight, axis=1)
tf_df['chipseq_targets_weighted'] = (
    tf_df['chipseq_targets'] * tf_df['diversity_weight']
).astype(int)

# Edge confidence stats per TF
edge_tf = edges.groupby('tf_locus').agg(
    mean_confidence=('predicted_confidence', 'mean'),
    n_metabolic_targets=('target_mapped_reaction_count', lambda x: (x > 0).sum()),
    n_shared_pathway=('tf_target_share_metabolic_pathway', 'sum'),
    mean_expr_corr=('expression_correlation', lambda x: x.abs().mean()),
).reset_index()
tf_df = tf_df.merge(edge_tf, left_on='locus', right_on='tf_locus', how='left').drop(columns=['tf_locus'])

# TF-to-TF regulation
tf_set = set(tf_df['locus'])
tf_tg_reg = reg[reg['TG_locusTag'].isin(tf_set)]
tf2tf_out = tf_tg_reg.groupby('TF_locusTag').size().rename('tf_targets_count')
tf2tf_in  = tf_tg_reg.groupby('TG_locusTag').size().rename('tf_regulators_count')

tf_df = tf_df.merge(tf2tf_out, left_on='locus', right_index=True, how='left')
tf_df = tf_df.merge(tf2tf_in,  left_on='locus', right_index=True, how='left')
tf_df['tf_targets_count']    = tf_df['tf_targets_count'].fillna(0).astype(int)
tf_df['tf_regulators_count'] = tf_df['tf_regulators_count'].fillna(0).astype(int)

print(f"  TF-to-TF edges: {len(tf_tg_reg)}")


# ─────────────────────────────────────────────────────────────────────────────
# 3. Tier 分类
# ─────────────────────────────────────────────────────────────────────────────
print("\n[3/7] Classifying TF tiers...")

def norm(s):
    mn, mx = s.min(), s.max()
    return (s - mn) / (mx - mn) if mx > mn else s * 0

tf_df['norm_out']        = norm(tf_df['out_degree'])
tf_df['norm_between']    = norm(tf_df['betweenness'])
tf_df['norm_tf_targets'] = norm(tf_df['tf_targets_count'])
tf_df['norm_chip']       = norm(tf_df['chipseq_targets_weighted'])  # ← 已做偏差修正

tf_df['master_score'] = (
    0.35 * tf_df['norm_out'] +
    0.25 * tf_df['norm_between'] +
    0.25 * tf_df['norm_tf_targets'] +
    0.15 * tf_df['norm_chip']   # ChIP 权重用修正后的值
)

def assign_tier(row):
    if row['master_score'] >= 0.30 or row['out_degree'] >= 50:
        return 'Tier 1 – Master Regulator'
    elif row['master_score'] >= 0.10 or row['out_degree'] >= 10:
        return 'Tier 2 – Intermediate Regulator'
    else:
        return 'Tier 3 – Local Regulator'

tf_df['tier'] = tf_df.apply(assign_tier, axis=1)
tf_df = tf_df.sort_values('master_score', ascending=False).reset_index(drop=True)
tf_df['rank'] = tf_df.index + 1

for t, n in sorted(tf_df['tier'].value_counts().items()):
    print(f"  {t}: {n} TFs")


# ─────────────────────────────────────────────────────────────────────────────
# 4. COG 功能分析
# ─────────────────────────────────────────────────────────────────────────────
print("\n[4/7] COG functional analysis per TF regulon...")

COG_CATEGORIES = {
    'J': 'Translation', 'K': 'Transcription', 'L': 'Replication/repair',
    'D': 'Cell cycle', 'T': 'Signal transduction', 'M': 'Cell wall',
    'O': 'Post-translational mod.', 'C': 'Energy production',
    'G': 'Carbohydrate metabolism', 'E': 'Amino acid metabolism',
    'F': 'Nucleotide metabolism', 'H': 'Coenzyme metabolism',
    'I': 'Lipid metabolism', 'P': 'Inorganic ion transport',
    'Q': 'Secondary metabolites', 'R': 'General function', 'S': 'Unknown',
}

all_regs = pd.concat([
    reg[['TF_locusTag', 'TG_locusTag']],
    chip[['TF_locusTag', 'TG_locusTag']]
]).drop_duplicates()

regulon_cog = {}
for tf_locus in tf_df['locus']:
    targets = all_regs[all_regs['TF_locusTag'] == tf_locus]['TG_locusTag'].unique()
    cats = defaultdict(int)
    for tg in targets:
        cat = cog_db.get(tg, {}).get('category', 'S')
        cats[cat] += 1
    regulon_cog[tf_locus] = dict(cats)

tier_cogs = {'T1': defaultdict(int), 'T2': defaultdict(int), 'T3': defaultdict(int)}
tier_keys  = {'Tier 1 – Master Regulator': 'T1',
              'Tier 2 – Intermediate Regulator': 'T2',
              'Tier 3 – Local Regulator': 'T3'}

for _, row in tf_df.iterrows():
    tk = tier_keys[row['tier']]
    for cat, cnt in regulon_cog.get(row['locus'], {}).items():
        tier_cogs[tk][cat] += cnt

print("  Top COG categories in Tier 1 regulons:")
for cat, cnt in sorted(tier_cogs['T1'].items(), key=lambda x: -x[1])[:5]:
    print(f"    [{cat}] {COG_CATEGORIES.get(cat, cat)}: {cnt} genes")


# ─────────────────────────────────────────────────────────────────────────────
# 5. iModulon 交叉验证
# ─────────────────────────────────────────────────────────────────────────────
print("\n[5/7] iModulon cross-validation...")

imod_tf_support = {}
for meta in imod_meta:                        # list of iModulon objects
    imod_name = meta.get('name', meta.get('id', ''))
    tf_assoc  = meta.get('linked_regulator', '') or ''
    if isinstance(tf_assoc, list):
        tf_assoc = ','.join(tf_assoc)
    for tf_locus, tf_name in name_map.items():
        if tf_name.lower() in tf_assoc.lower() or tf_locus.lower() in tf_assoc.lower():
            imod_tf_support.setdefault(tf_locus, []).append(imod_name)

tf_df['imodulon_count'] = tf_df['locus'].map(lambda x: len(imod_tf_support.get(x, [])))
tf_df['imodulon_names'] = tf_df['locus'].map(lambda x: '; '.join(imod_tf_support.get(x, [])))
print(f"  {(tf_df['imodulon_count']>0).sum()} TFs have iModulon support")


# ─────────────────────────────────────────────────────────────────────────────
# 6. 保存结果表格
# ─────────────────────────────────────────────────────────────────────────────
print("\n[6/7] Saving results tables...")

output_cols = ['rank','locus','name','tier','master_score','out_degree','betweenness',
               'pagerank','tf_targets_count','tf_regulators_count','chipseq_targets',
               'chipseq_coverage','curated_targets','n_metabolic_targets',
               'mean_confidence','activation_ratio','is_sigma','imodulon_count','imodulon_names']
output_cols = [c for c in output_cols if c in tf_df.columns]
result_df = tf_df[output_cols].copy()
for col in ['master_score','chipseq_coverage','mean_confidence']:
    if col in result_df: result_df[col] = result_df[col].round(4)
for col in ['betweenness','pagerank']:
    if col in result_df: result_df[col] = result_df[col].round(6)

result_path = os.path.join(OUT, 'tf_hierarchy_rankings.csv')
result_df.to_csv(result_path, index=False)
print(f"  Saved: {result_path}")

# Console top-20 table
print("\n  ── Top 20 Regulators ──────────────────────────────────────────────")
print(f"  {'Rank':<5} {'Locus':<10} {'Name':<12} {'Tier':<8} {'Score':<7} {'OutDeg':<8} {'ChIP':<7} {'TF→TF'}")
print("  " + "-"*68)
for _, row in result_df.head(20).iterrows():
    ts = row['tier'].split('–')[0].strip()
    print(f"  {int(row['rank']):<5} {row['locus']:<10} {str(row['name']):<12} "
          f"{ts:<8} {row['master_score']:<7.4f} {int(row['out_degree']):<8} "
          f"{int(row['chipseq_targets']):<7} {int(row['tf_targets_count'])}")


# ─────────────────────────────────────────────────────────────────────────────
# 7. 可视化（3 张图）
# ─────────────────────────────────────────────────────────────────────────────
print("\n[7/7] Generating figures...")

TIER_COLORS = {
    'Tier 1 – Master Regulator':       '#e11d48',
    'Tier 2 – Intermediate Regulator': '#f59e0b',
    'Tier 3 – Local Regulator':        '#64748b',
}

# ── Fig 1: Overview (bubble + bar + COG heatmap) ───────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle('Hierarchical Classification of C. glutamicum Transcription Factors',
             fontsize=14, fontweight='bold', y=1.01)

# Panel A: scatter
ax = axes[0]
for tier, grp in tf_df.groupby('tier'):
    sizes = (grp['chipseq_targets'] + 1) * 8
    ax.scatter(grp['out_degree'], grp['betweenness'] * 1000,
               c=TIER_COLORS[tier], s=sizes, alpha=0.75,
               label=tier.split('–')[1].strip(), edgecolors='white', linewidth=0.5)
for _, row in tf_df.head(10).iterrows():
    ax.annotate(row['name'], (row['out_degree'], row['betweenness'] * 1000),
                fontsize=7, ha='left', va='bottom', xytext=(3,3), textcoords='offset points')
ax.set_xlabel('Out-degree (target genes)', fontsize=11)
ax.set_ylabel('Betweenness Centrality (×10⁻³)', fontsize=11)
ax.set_title('(A) Network Position', fontsize=12)
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# Panel B: horizontal bar
ax = axes[1]
top25 = tf_df.head(25)
colors = [TIER_COLORS[t] for t in top25['tier']]
ax.barh(range(len(top25)), top25['master_score'], color=colors, alpha=0.85)
ax.set_yticks(range(len(top25)))
ax.set_yticklabels(top25['name'], fontsize=8)
ax.invert_yaxis()
ax.set_xlabel('Master Regulator Score', fontsize=11)
ax.set_title('(B) Top 25 Ranked TFs', fontsize=12)
for i, (_, row) in enumerate(top25.iterrows()):
    ax.text(row['master_score'] + 0.002, i, f"n={int(row['out_degree'])}",
            va='center', fontsize=6.5, color='#374151')
ax.grid(axis='x', alpha=0.3)

# Panel C: COG heatmap
ax = axes[2]
cog_cats_use = ['C','E','G','H','I','P','K','T','M','J','L','O','R','S']
tier_cog_list = [tier_cogs['T1'], tier_cogs['T2'], tier_cogs['T3']]
heatmap_data = np.zeros((len(cog_cats_use), 3))
for j, tcog in enumerate(tier_cog_list):
    total = sum(tcog.values()) or 1
    for i, cat in enumerate(cog_cats_use):
        heatmap_data[i, j] = tcog.get(cat, 0) / total * 100

cmap = LinearSegmentedColormap.from_list('tier', ['#f8fafc', '#0ea5e9', '#1e3a5f'])
im = ax.imshow(heatmap_data, cmap=cmap, aspect='auto', vmin=0)
ax.set_xticks(range(3))
ax.set_xticklabels(['Tier 1\nMaster', 'Tier 2\nIntermediate', 'Tier 3\nLocal'], fontsize=9)
ax.set_yticks(range(len(cog_cats_use)))
ax.set_yticklabels([f'[{c}] {COG_CATEGORIES.get(c,c)[:22]}' for c in cog_cats_use], fontsize=7.5)
ax.set_title('(C) Regulon COG Composition (%)', fontsize=12)
plt.colorbar(im, ax=ax, label='% of regulon', shrink=0.8)
for i in range(len(cog_cats_use)):
    for j in range(3):
        v = heatmap_data[i,j]
        if v > 3:
            ax.text(j, i, f'{v:.1f}', ha='center', va='center',
                    fontsize=7, color='white' if v > 10 else 'black')

plt.tight_layout()
p = os.path.join(OUT, 'fig1_tf_hierarchy_overview.png')
plt.savefig(p, dpi=200, bbox_inches='tight'); plt.close()
print(f"  Saved: {p}")

# ── Fig 2: TF-to-TF cascade network ────────────────────────────────────────
try:
    import networkx as nx
    tier1_loci = set(tf_df[tf_df['tier']=='Tier 1 – Master Regulator']['locus'])
    tier2_loci = set(tf_df[tf_df['tier']=='Tier 2 – Intermediate Regulator']['locus'])
    tf_all     = tier1_loci | tier2_loci

    G = nx.DiGraph()
    for _, row in tf_tg_reg.iterrows():
        if row['TF_locusTag'] in tf_all and row['TG_locusTag'] in tf_all:
            G.add_edge(row['TF_locusTag'], row['TG_locusTag'],
                       role=row.get('Role','A'))

    if len(G.nodes) > 0:
        fig2, ax2 = plt.subplots(figsize=(12, 10))
        ax2.set_facecolor('#0f172a'); fig2.patch.set_facecolor('#0f172a')
        pos = nx.spring_layout(G, k=2.5, seed=42, iterations=80)

        nc = ['#e11d48' if n in tier1_loci else '#f59e0b' for n in G.nodes()]
        ns = [1200 if n in tier1_loci else 500 for n in G.nodes()]
        ec = ['#10b981' if G[u][v].get('role','A')=='A' else '#f43f5e' for u,v in G.edges()]

        nx.draw_networkx_nodes(G, pos, node_color=nc, node_size=ns, alpha=0.9, ax=ax2)
        nx.draw_networkx_edges(G, pos, edge_color=ec, alpha=0.65, arrows=True,
                               arrowsize=15, width=1.5,
                               connectionstyle='arc3,rad=0.1', ax=ax2)
        nx.draw_networkx_labels(G, pos, {n: get_name(n) for n in G.nodes()},
                                font_size=9, font_color='white', ax=ax2)
        ax2.set_title('TF–TF Regulatory Cascade Network\n(red=Master | amber=Intermediate | '
                      'green=activation | rose=repression)',
                      color='white', fontsize=12, pad=15)
        ax2.axis('off')
        legend_els = [mpatches.Patch(color=c, label=l) for c, l in [
            ('#e11d48','Tier 1: Master'), ('#f59e0b','Tier 2: Intermediate'),
            ('#10b981','Activation'), ('#f43f5e','Repression')]]
        ax2.legend(handles=legend_els, loc='lower left',
                   facecolor='#1e293b', labelcolor='white', fontsize=9)
        plt.tight_layout()
        p = os.path.join(OUT, 'fig2_tf_cascade_network.png')
        plt.savefig(p, dpi=200, bbox_inches='tight', facecolor='#0f172a'); plt.close()
        print(f"  Saved: {p}  ({len(G.nodes())} nodes, {len(G.edges())} edges)")
except ImportError:
    print("  networkx not available – skipping cascade figure")

# ── Fig 3: Regulon evidence composition ────────────────────────────────────
fig3, ax3 = plt.subplots(figsize=(12, 5))
top15 = tf_df.head(15)
x = np.arange(len(top15))
chip_n    = top15['chipseq_targets'].values
curated_x = np.maximum(top15['curated_targets'].values - chip_n, 0)
ax3.bar(x, chip_n,    0.55, label='ChIP-seq confirmed', color='#0ea5e9', alpha=0.9)
ax3.bar(x, curated_x, 0.55, bottom=chip_n, label='Curated (non-ChIP)', color='#6366f1', alpha=0.8)
ax3.set_xticks(x)
ax3.set_xticklabels(top15['name'], rotation=35, ha='right', fontsize=10)
ax3.set_ylabel('Number of Regulated Genes', fontsize=11)
ax3.set_title('Evidence Composition of Top 15 TF Regulons in C. glutamicum', fontsize=13)
ax3.legend(fontsize=10)
ax3.grid(axis='y', alpha=0.3)
for i, (_, row) in enumerate(top15.iterrows()):
    ax3.get_xticklabels()[i].set_color(TIER_COLORS[row['tier']])
plt.tight_layout()
p = os.path.join(OUT, 'fig3_regulon_evidence_composition.png')
plt.savefig(p, dpi=200, bbox_inches='tight'); plt.close()
print(f"  Saved: {p}")


# ── 论文用摘要统计 ────────────────────────────────────────────────────────────
t1 = tf_df[tf_df['tier']=='Tier 1 – Master Regulator']
t2 = tf_df[tf_df['tier']=='Tier 2 – Intermediate Regulator']
t3 = tf_df[tf_df['tier']=='Tier 3 – Local Regulator']

print("\n" + "="*65)
print("  SUMMARY FOR MANUSCRIPT")
print("="*65)
print(f"  Network: {centrality['_meta']['n_nodes']} genes, "
      f"{centrality['_meta']['n_edges']} edges, {centrality['_meta']['n_tfs']} TFs")
print(f"\n  Tier 1 (Master):        {len(t1):>3} TFs | "
      f"avg targets {t1['out_degree'].mean():.1f} | "
      f"avg ChIP {t1['chipseq_targets'].mean():.1f} | "
      f"regulate TFs: {(t1['tf_targets_count']>0).sum()}")
print(f"  Tier 2 (Intermediate):  {len(t2):>3} TFs | "
      f"avg targets {t2['out_degree'].mean():.1f} | "
      f"avg ChIP {t2['chipseq_targets'].mean():.1f}")
print(f"  Tier 3 (Local):         {len(t3):>3} TFs | "
      f"avg targets {t3['out_degree'].mean():.1f}")
print(f"  TF-TF cascade edges:    {len(tf_tg_reg)}")
print(f"  TFs with iModulon:      {(tf_df['imodulon_count']>0).sum()}")
print(f"\n  Output: {OUT}")
print("  [OK] tf_hierarchy_rankings.csv")
print("  [OK] fig1_tf_hierarchy_overview.png")
print("  [OK] fig2_tf_cascade_network.png")
print("  [OK] fig3_regulon_evidence_composition.png")
print("\n  Analysis complete!")
