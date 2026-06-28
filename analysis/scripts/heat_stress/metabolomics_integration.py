"""
Metabolomics Integration Analysis — Heat Stress (40°C vs 30°C)
==============================================================
1. Key metabolite time-series (1h / 4h / 24h)
2. AKG/Glu/Lys bifurcation evidence panel
3. Model-flux vs metabolomics FC correlation
4. TCA / redox / nitrogen overview heatmap
"""
import sys, os, csv, math
import warnings; warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding='utf-8')

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
DATA_DIR = os.path.join(ROOT_DIR, "data", "raw", "Metabolome")
OUT_DIR  = os.path.join(ROOT_DIR, "analysis", "outputs", "heat_stress")
os.makedirs(OUT_DIR, exist_ok=True)
TPS      = ['1h', '4h', '24h']
MODES    = ['NEG']

# ── Load all metabolomics data ─────────────────────────────────────────────────
def load_metabo():
    """Returns dict: {timepoint: {mode: [(compound, fc, lfc, pval), ...]}}"""
    data = {}
    for tp in TPS:
        data[tp] = {}
        for mode in MODES:
            path = os.path.join(DATA_DIR, tp, f"{tp}_{mode}_volcano.csv")
            rows = []
            if os.path.exists(path):
                with open(path, encoding='utf-8') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        cmp  = row.get('Compounds','').strip().strip('"')
                        fc   = row.get('FC','').strip().strip('"')
                        lfc  = row.get('log2(FC)','').strip().strip('"')
                        pval = row.get('p.value','').strip().strip('"')
                        if cmp:
                            try: rows.append((cmp, float(fc), float(lfc), float(pval)))
                            except: pass
            data[tp][mode] = rows
    return data

def find_compound(data, keywords, tp, prefer_mode='NEG'):
    """Find best match for a metabolite at a given timepoint."""
    kw_lower = [k.lower() for k in keywords]
    for mode in ([prefer_mode] + [m for m in MODES if m != prefer_mode]):
        for (cmp, fc, lfc, pval) in data[tp].get(mode, []):
            if any(k in cmp.lower() for k in kw_lower):
                return {'compound':cmp, 'fc':fc, 'lfc':lfc, 'pval':pval, 'mode':mode}
    return None

raw = load_metabo()
print("Data loaded.")

# ── Key metabolites of interest ────────────────────────────────────────────────
TARGETS = {
    'α-Ketoglutarate\n(AKG)':        ['alpha-ketoglutaric acid', '2-oxoglutarate'],
    'Glutamate\n(Glu)':              ['glutamate', 'l-glutamic acid', 'dl-glutamic acid'],
    'N-Acetylaspartate\n(Asp proxy)':['n-acetylaspartate', 'n-acetyl-aspartate'],
    'N6-Acetyl-Lys\n(Lys proxy)':    ['n6-acetyl-l-lysine', 'n6-acetyl-lysine'],
    'Glutamine\n(Gln)':              ['glutamine', 'd-glutamine', 'l-glutamine'],
    'Citrate /\nIsocitrate':         ['citric acid', 'isocitrate'],
    'Malate':                         ['malic acid', 'd-(+)-malic acid', '(s)-malate'],
    'Succinate':                      ['succinic acid', 'succinate'],
    'NAD':                            ['nad'],
    'NADP':                           ['nadp'],
}

# Build time-series matrix
ts_data = {}
for label, kws in TARGETS.items():
    ts_data[label] = {}
    for tp in TPS:
        hit = find_compound(raw, kws, tp)
        ts_data[label][tp] = hit

print("\n=== Key metabolite time-series (log2FC, 40°C vs 30°C) ===")
print(f"{'Metabolite':<32s}  {'1h':>8}  {'4h':>8}  {'24h':>8}")
print('-'*60)
for label, tpts in ts_data.items():
    row_vals = []
    for tp in TPS:
        h = tpts[tp]
        if h:
            sig = '**' if h['pval']<0.05 else ('*' if h['pval']<0.1 else '  ')
            row_vals.append(f"{h['lfc']:+.2f}{sig}")
        else:
            row_vals.append('  ---  ')
    print(f"{label.replace(chr(10),' '):<32s}  {row_vals[0]:>8}  {row_vals[1]:>8}  {row_vals[2]:>8}")

# ── Model predictions at 40°C vs 30°C ─────────────────────────────────────────
# From flux_partition_analysis output
MODEL_PRED = {
    'Growth rate (FBA)':      math.log2(0.0585/0.0657),   # -0.17
    'ASPK flux (FBA)':        math.log2(0.0492/0.0553),   # -0.17
    'GDH flux (FBA)':         math.log2(0.4012/0.4502),   # -0.17
    'AKGDH flux (FBA)':       math.log2(0.0039/0.0055),   # -0.49
    'Lys_FBA':                math.log2(0.0875/0.2096),   # -1.26
    'Glu_FBA':                math.log2(0.2056/0.2224),   # -0.11
    'α_LysC(40°C)':           math.log2(0.184/1.0),       # -2.44
    'α_GDH(40°C)':            math.log2(0.60/1.0),        # -0.74
}

print("\n=== Model predictions at 40°C (log2FC vs 30°C) ===")
for k,v in MODEL_PRED.items():
    print(f"  {k:<30s}  {v:+.3f}")

# ── Figure ─────────────────────────────────────────────────────────────────────
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

plt.rcParams.update({
    'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
    'font.size':9,'axes.labelsize':9.5,'axes.titlesize':10,'axes.titleweight':'bold',
    'xtick.labelsize':8.5,'ytick.labelsize':8.5,'legend.fontsize':8,
    'figure.dpi':300,'axes.linewidth':0.8,
    'axes.spines.top':False,'axes.spines.right':False,
})

def styled(ax):
    ax.set_facecolor('white')
    ax.tick_params(direction='out',length=3,width=0.8,colors='black')
    for sp in ['left','bottom']:
        ax.spines[sp].set_color('black'); ax.spines[sp].set_linewidth(0.8)
    return ax

tp_x  = [1, 2, 3]
tp_lb = ['1h', '4h', '24h']

# Colour scheme
CLR = {
    'AKG':  '#e6550d',   # orange
    'Glu':  '#3182bd',   # blue
    'Asp':  '#31a354',   # green
    'Lys':  '#de2d26',   # red
    'Gln':  '#756bb1',   # purple
    'TCA':  '#636363',   # grey
    'NAD':  '#fd8d3c',   # amber
}

fig = plt.figure(figsize=(18, 14))
fig.patch.set_facecolor('white')
gs  = gridspec.GridSpec(3, 4, figure=fig, hspace=0.55, wspace=0.45)

def plot_ts(ax, labels_kws, title, ylabel='log₂FC (40°C/30°C)',
            ref_line=True, ylim=None):
    styled(ax)
    ax.axhline(0, color='#999', lw=0.8, ls='--', zorder=0)
    for label, kws, clr, marker in labels_kws:
        ys, errs, sigs = [], [], []
        for tp in TPS:
            hit = find_compound(raw, kws, tp)
            if hit:
                ys.append(hit['lfc'])
                sigs.append('**' if hit['pval']<0.05 else ('*' if hit['pval']<0.1 else ''))
            else:
                ys.append(np.nan)
                sigs.append('')
        ax.plot(tp_x, ys, color=clr, lw=2, marker=marker, ms=7, label=label, zorder=3)
        for i,(y,s) in enumerate(zip(ys,sigs)):
            if not np.isnan(y) and s:
                ax.text(tp_x[i], y + (0.12 if y>=0 else -0.22), s, ha='center',
                        fontsize=8, color=clr, fontweight='bold')
    ax.set_xticks(tp_x); ax.set_xticklabels(tp_lb)
    ax.set_xlabel('Time at 40°C'); ax.set_ylabel(ylabel)
    ax.set_title(title)
    if ylim: ax.set_ylim(ylim)
    ax.legend(frameon=True, edgecolor='#ccc', fancybox=False, fontsize=7.5,
              loc='best')
    return ax

# Panel A: AKG accumulation (core finding)
ax_A = fig.add_subplot(gs[0,0])
plot_ts(ax_A,
    [('α-KG (AKG)', ['alpha-ketoglutaric acid'], CLR['AKG'], 'o')],
    '(A)  AKG accumulation\n(AKGDH+GDH suppressed)', ylim=(-1, 5))
ax_A.axhline(0, color='#999', lw=0.8, ls='--')
ax_A.fill_between([0.8,3.2],[0,0],[3.5,3.5], alpha=0.05, color=CLR['AKG'])
ax_A.text(2, 3.3, 'AKG ↑ → AKGDH+GDH\nheat inactivation', ha='center',
          fontsize=7.5, color=CLR['AKG'])

# Panel B: Glutamate decrease
ax_B = fig.add_subplot(gs[0,1])
plot_ts(ax_B,
    [('Glutamate', ['glutamic acid','l-glutamic acid','dl-glutamic acid'],
      CLR['Glu'], 'o'),
     ('Glutamate dipeptide (Glu-Gln)', ['glu-gln', 'glu gln'], '#08519c', 's')],
    '(B)  Glutamate dynamics\n(GDH heat inactivation)', ylim=(-2.5, 3))

# Panel C: Aspartate/Lysine collapse
ax_C = fig.add_subplot(gs[0,2])
plot_ts(ax_C,
    [('N-Acetylaspartate (Asp)', ['n-acetylaspartate'], CLR['Asp'], 'o'),
     ('N6-Acetyl-Lys (Lys)', ['n6-acetyl-l-lysine'], CLR['Lys'], 's')],
    '(C)  Asp/Lys pathway collapse\n(LysC heat inactivation)', ylim=(-7, 4))

# Panel D: Glutamine (N-metabolism)
ax_D = fig.add_subplot(gs[0,3])
plot_ts(ax_D,
    [('Glutamine', ['glutamine'], CLR['Gln'], 'o'),
     ('α-N-Acetyl-Gln', ['alpha-n-acetyl-glutamine'], '#bcbddc', 's')],
    '(D)  Glutamine / N-cycle\ndynamics', ylim=(-6, 2))

# Panel E: TCA intermediates
ax_E = fig.add_subplot(gs[1,0])
plot_ts(ax_E,
    [('Citrate/Isocitrate', ['citric acid','isocitrate'], '#636363', 'o'),
     ('Malate', ['malic acid','d-(+)-malic acid'], '#969696', 's'),
     ('Succinate', ['succinic acid'], '#bdbdbd', '^')],
    '(E)  TCA intermediates\n(carbon backbone)', ylim=(-2.5, 2))

# Panel F: Cofactors
ax_F = fig.add_subplot(gs[1,1])
plot_ts(ax_F,
    [('NAD', ['nad'], CLR['NAD'], 'o'),
     ('NADP+', ['nadp'], '#fdae6b', 's')],
    '(F)  Redox cofactors\n(NAD/NADP)', ylim=(-2, 4))

# Panel G: Model vs metabolomics comparison at 24h (correlation bar chart)
ax_G = styled(fig.add_subplot(gs[1,2:4]))
# Paired data: {metabolite: (metabo_lfc_24h, model_pred_log2, label)}
pairs = {
    'AKG': (None, None, 'alpha-ketoglutaric acid', 'AKGDH flux (FBA)'),
    'Glu': (None, None, 'l-glutamic acid', 'Glu_FBA'),
    'NAAsp': (None, None, 'n-acetylaspartate', 'Lys_FBA'),
    'N6ALys': (None, None, 'n6-acetyl-l-lysine', 'Lys_FBA'),
    'Gln': (None, None, 'glutamine', 'Growth rate (FBA)'),
}
comp_names, metabo_vals, model_vals, pvals_comp = [], [], [], []
metabo_targets_24h = {
    'AKG (24h)':           (['alpha-ketoglutaric acid'], 'AKGDH flux (FBA)', 1.0),
    'Glu (24h)':           (['l-glutamic acid','glutamic acid','dl-glutamic acid'], 'Glu_FBA', 1.0),
    'N-Acetylaspartate\n(24h)': (['n-acetylaspartate'], 'Lys_FBA', -1.0),
    'N6-Acetyl-Lys (24h)': (['n6-acetyl-l-lysine'], 'Lys_FBA', -1.0),
    'Glutamine (24h)':     (['glutamine'], 'Growth rate (FBA)', 1.0),
    'Citrate (1h)':        (['citric acid'], 'AKGDH flux (FBA)', -1.0),
    'NAD (1h)':            (['nad'], 'α_LysC(40°C)', -1.0),
}

bar_labels, bar_metabo, bar_model, bar_pvals = [], [], [], []
for label, (kws, model_key, _sign) in metabo_targets_24h.items():
    # Use 24h if available, else 4h
    tp = '24h' if '24h' in label else ('1h' if '1h' in label else '4h')
    hit = find_compound(raw, kws, tp)
    if hit:
        bar_labels.append(label.replace('\n',' '))
        bar_metabo.append(hit['lfc'])
        bar_model.append(MODEL_PRED.get(model_key, 0))
        bar_pvals.append(hit['pval'])

x  = np.arange(len(bar_labels))
w  = 0.35
b1 = ax_G.bar(x - w/2, bar_metabo, w, color='#4292c6', alpha=0.85, label='Metabolomics log₂FC')
b2 = ax_G.bar(x + w/2, bar_model,  w, color='#ef6548', alpha=0.85, label='Model prediction log₂')
ax_G.axhline(0, color='#555', lw=0.8)
ax_G.set_xticks(x); ax_G.set_xticklabels(bar_labels, fontsize=7.5, rotation=15, ha='right')
ax_G.set_ylabel('log₂FC / log₂(predicted change)')
ax_G.set_title('(G)  Metabolomics vs Model — directional concordance')
ax_G.legend(frameon=True, edgecolor='#ccc', fancybox=False)
# Concordance check
for i,(m,p,pv) in enumerate(zip(bar_metabo, bar_model, bar_pvals)):
    agree = (m*p > 0) or (abs(p) < 0.1)
    sym = '✓' if agree else '✗'
    color = '#2ca02c' if agree else '#d62728'
    ax_G.text(i, max(abs(m),abs(p))*1.05 + 0.2, sym, ha='center', fontsize=10, color=color)

# Panel H: Comprehensive heatmap of all significant metabolites (top hits by |lfc|)
ax_H = styled(fig.add_subplot(gs[2,:]))
# Build matrix for selected metabolites (all time points)
selected_for_heatmap = [
    ('α-Ketoglutaric acid',       ['alpha-ketoglutaric acid'],   'AKG'),
    ('Glutamic acid',             ['l-glutamic acid','dl-glutamic acid'], 'Glu'),
    ('Glutamine',                 ['glutamine'],                   'Gln'),
    ('N-Acetylaspartate',         ['n-acetylaspartate'],          'N-AcAsp'),
    ('N6-Acetyl-Lys',             ['n6-acetyl-l-lysine'],         'N6AcLys'),
    ('Citric acid',               ['citric acid'],                 'Citrate'),
    ('Malic acid',                ['malic acid','d-(+)-malic acid'], 'Malate'),
    ('Succinic acid',             ['succinic acid'],               'Succinate'),
    ('NAD',                       ['nad'],                         'NAD'),
    ('3-Phosphoglyceric acid',    ['3-phosphoglyceric acid'],      '3PG'),
    ('N-Acetyl-Glu',              ['n-acetyl-glutamic acid'],      'N-AcGlu'),
    ('Citrulline',                ['citrulline'],                  'Citrulline'),
    ('N-Formyl-Glu',              ['n-formyl-l-glutamic acid'],    'N-FmGlu'),
    ('N-MeAsp',                   ['n-methylaspartate'],           'N-MeAsp'),
    ('Pyruvate',                  ['pyruvic acid'],                'Pyruvate'),
]

hm_rows, hm_cols, hm_mat, hm_pval = [], [], [], []
hm_cols = ['1h', '4h', '24h']
for (fullname, kws, short) in selected_for_heatmap:
    row_lfc, row_pv = [], []
    for tp in TPS:
        for mode in MODES:
            hit = find_compound(raw, kws, tp, prefer_mode=mode)
            for (cmp, fc, lfc, pval) in raw[tp].get(mode, []):
                if any(k in cmp.lower() for k in [x.lower() for x in kws]):
                    row_lfc.append(lfc)
                    row_pv.append(pval)
                    break
            else:
                if len(row_lfc) < (TPS.index(tp)*len(MODES) + MODES.index(mode) + 1):
                    row_lfc.append(np.nan)
                    row_pv.append(1.0)
    # Pad if needed
    while len(row_lfc) < 3: row_lfc.append(np.nan)
    while len(row_pv) < 3: row_pv.append(1.0)
    hm_rows.append(short)
    hm_mat.append(row_lfc[:3])
    hm_pval.append(row_pv[:3])

hm_arr  = np.array(hm_mat, dtype=float)
pv_arr  = np.array(hm_pval, dtype=float)

from matplotlib.colors import TwoSlopeNorm
vmax = np.nanmax(np.abs(hm_arr))
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
im = ax_H.imshow(hm_arr, aspect='auto', cmap='RdBu_r', norm=norm)
ax_H.set_xticks(range(3)); ax_H.set_xticklabels(hm_cols, fontsize=8.5)
ax_H.set_yticks(range(len(hm_rows))); ax_H.set_yticklabels(hm_rows, fontsize=8.5)
ax_H.set_title('(H)  Metabolite log₂FC heatmap — 40°C vs 30°C across time points (NEG mode)', fontsize=10, fontweight='bold')
plt.colorbar(im, ax=ax_H, shrink=0.6, label='log₂FC (40°C/30°C)')
# Asterisks for significance
for i in range(len(hm_rows)):
    for j in range(3):
        v = hm_arr[i,j]
        p = pv_arr[i,j]
        if not np.isnan(v):
            sig = '**' if p<0.05 else ('*' if p<0.1 else '')
            txt_col = 'white' if abs(v) > vmax*0.6 else 'black'
            ax_H.text(j, i, f'{v:+.1f}{sig}', ha='center', va='center',
                      fontsize=7, color=txt_col)

# Draw category boxes
for i,y_sep in enumerate([-0.5, 0.5, 1.5, 2.5, 4.5, 7.5, 8.5, 9.5, 10.5]):
    ax_H.axhline(y_sep, color='#aaa', lw=0.5, ls=':')

fig.suptitle('Metabolomics Integration — C. glutamicum Heat Stress (40°C vs 30°C)\n'
             'Validation of ecFBA Model Predictions', fontsize=11, fontweight='bold', y=1.01)

out = os.path.join(OUT_DIR, 'metabolomics_integration.png')
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close(fig)
print(f'\nSaved: {out}')
print('Done.')
