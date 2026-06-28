import os, sys
import pandas as pd
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
prot_dir = os.path.join(ROOT_DIR, "data", "raw", "protemic")
out_dir = os.path.join(ROOT_DIR, "analysis", "outputs", "heat_stress")
os.makedirs(out_dir, exist_ok=True)
tps = ['1h', '4h', '24h']
tp_lb = ['1h (Acute)', '4h (Adaptation)', '24h (Chronic)']
tp_x = [1, 2, 3]

# Load proteomics tables
prot_data = {}
for tp in tps:
    path = os.path.join(prot_dir, f"{tp}_all_results.csv")
    if os.path.exists(path):
        df = pd.read_csv(path)
        prot_data[tp] = df.set_index('cgl_id')

# Pathway definitions
PATHWAYS = {
    'HSPs (Chaperones)': {
        'genes': ['Cgl0597', 'Cgl2716', 'Cgl0598', 'Cgl2800', 'Cgl2780', 'Cgl2799', 'Cgl2798'],
        'color': '#de2d26', 'marker': 'o', 'name': 'groES, groEL, groEL\', dnaK, clpB, grpE, dnaJ'
    },
    'Glycolysis': {
        'genes': ['Cgl0851', 'Cgl1250', 'Cgl2770', 'Cgl1586', 'Cgl2089', 'Cgl1588', 'Cgl0974'],
        'color': '#3182bd', 'marker': 's', 'name': 'pgi, pfkA, fda, tpi, pyk, gap, eno'
    },
    'TCA Cycle': {
        'genes': ['Cgl0829', 'Cgl0664', 'Cgl1129', 'Cgl2380', 'Cgl0371', 'Cgl0372', 'Cgl0370'],
        'color': '#31a354', 'marker': '^', 'name': 'gltA, icd, odhA, mdh, sdhA, sdhB, sdhCD'
    },
    'Lysine Synthesis': {
        'genes': ['Cgl0251', 'Cgl0252', 'Cgl1971', 'Cgl1973', 'Cgl1106', 'Cgl1109', 'Cgl1943', 'Cgl1180'],
        'color': '#e6550d', 'marker': 'd', 'name': 'lysC, asd, dapA, dapB, dapD, dapE, dapF, lysA'
    },
    'Glutamate/N-assim': {
        'genes': ['Cgl2079', 'Cgl0184', 'Cgl0185'],
        'color': '#756bb1', 'marker': 'v', 'name': 'gdh, gltB, gltD'
    }
}

# Calculate timeseries for each pathway
pathway_ts = {}
for pw, info in PATHWAYS.items():
    ys, es = [], []
    for tp in tps:
        df = prot_data[tp]
        present = [g for g in info['genes'] if g in df.index]
        vals = []
        for g in present:
            v = df.loc[g, 'logFC']
            if isinstance(v, pd.Series): v = v.mean()
            vals.append(v)
        ys.append(np.mean(vals) if vals else np.nan)
        es.append(np.std(vals) / np.sqrt(len(vals)) if vals else np.nan)
    pathway_ts[pw] = (ys, es)

# Setup plotting style
plt.rcParams.update({
    'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
    'font.size':9.5,'axes.labelsize':10,'axes.titlesize':10.5,'axes.titleweight':'bold',
    'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':8.5,
    'figure.dpi':300,'axes.linewidth':0.8,
    'axes.spines.top':False,'axes.spines.right':False,
})

fig = plt.figure(figsize=(15, 5))
fig.patch.set_facecolor('white')
gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.3)

# ── Panel A: Chaperone Induction ──
ax_A = fig.add_subplot(gs[0,0])
ax_A.set_facecolor('white')
ax_A.axhline(0, color='#999', lw=0.8, ls='--', zorder=1)

# Plot individual chaperone curves faintly
for g in PATHWAYS['HSPs (Chaperones)']['genes']:
    g_ys = []
    for tp in tps:
        df = prot_data[tp]
        if g in df.index:
            v = df.loc[g, 'logFC']
            if isinstance(v, pd.Series): v = v.mean()
            g_ys.append(v)
        else:
            g_ys.append(np.nan)
    gene_name = prot_data['1h'].loc[g, 'gene.name'] if g in prot_data['1h'].index else g
    ax_A.plot(tp_x, g_ys, color='#de2d26', alpha=0.3, lw=1, ls=':', zorder=2)
    ax_A.text(3.05, g_ys[-1], gene_name, color='#de2d26', alpha=0.6, fontsize=7.5, va='center')

# Plot mean chaperone curve
ys, es = pathway_ts['HSPs (Chaperones)']
ax_A.errorbar(tp_x, ys, yerr=es, color='#de2d26', fmt='-o', lw=2.5, ms=8, capsize=4, label='HSPs (Mean)', zorder=3)
ax_A.set_title('(A) Heat Shock Protein (HSP) Induction')
ax_A.set_xticks(tp_x); ax_A.set_xticklabels(tp_lb)
ax_A.set_ylabel('Protein Abundance log2FC (40°C vs 30°C)')
ax_A.set_ylim(-1, 6)
ax_A.legend(loc='upper left', frameon=True, edgecolor='#ccc', fancybox=False)

# ── Panel B: Metabolic Suppression ──
ax_B = fig.add_subplot(gs[0,1])
ax_B.set_facecolor('white')
ax_B.axhline(0, color='#999', lw=0.8, ls='--', zorder=1)

for pw in ['Glycolysis', 'TCA Cycle', 'Lysine Synthesis', 'Glutamate/N-assim']:
    ys, es = pathway_ts[pw]
    info = PATHWAYS[pw]
    ax_B.errorbar(tp_x, ys, yerr=es, color=info['color'], fmt='-'+info['marker'],
                  lw=2, ms=6.5, capsize=3, label=pw, zorder=3)

ax_B.set_title('(B) Metabolic Machinery Suppression')
ax_B.set_xticks(tp_x); ax_B.set_xticklabels(tp_lb)
ax_B.set_ylabel('Protein Abundance log2FC (40°C vs 30°C)')
ax_B.set_ylim(-1.2, 0.4)
ax_B.legend(loc='lower left', frameon=True, edgecolor='#ccc', fancybox=False)

# ── Panel C: Tug-of-War Contrast at 24h ──
ax_C = fig.add_subplot(gs[0,2])
ax_C.set_facecolor('white')
ax_C.axvline(0, color='#555', lw=0.8, zorder=1)

pws_to_compare = ['TCA Cycle', 'Glycolysis', 'Glutamate/N-assim', 'Lysine Synthesis', 'HSPs (Chaperones)']
labels_compare = ['TCA Cycle', 'Glycolysis', 'Glutamate/N-assim', 'Lysine Synthesis', 'HSPs']

fc_24h = [pathway_ts[pw][0][2] for pw in pws_to_compare]
err_24h = [pathway_ts[pw][1][2] for pw in pws_to_compare]
colors_24h = [PATHWAYS[pw]['color'] for pw in pws_to_compare]

bars = ax_C.barh(labels_compare, fc_24h, xerr=err_24h, height=0.45, color=colors_24h, edgecolor='black', lw=0.8, alpha=0.85, capsize=4, zorder=2)
# Add value labels
for bar in bars:
    width = bar.get_width()
    ax_C.text(width + (0.15 if width >= 0 else -0.55), bar.get_y() + bar.get_height()/2.,
              f"{width:+.2f}", va='center', ha='left' if width >= 0 else 'right', fontsize=8, fontweight='bold')

ax_C.set_title('(C) Proteome Reallocation at 24h')
ax_C.set_xlabel('log2FC (40°C vs 30°C)')
ax_C.set_xlim(-1.5, 5.8)

# Overall Title
fig.suptitle('Proteome Reallocation & Tug-of-War Effect\nExperimental Proteomics Validation of ecFBA Capacity Allocation',
             fontsize=12, fontweight='bold', y=1.02)

out_file = os.path.join(out_dir, "proteome_reallocation.png")
fig.savefig(out_file, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close(fig)

print(f"Successfully generated: {out_file}")
