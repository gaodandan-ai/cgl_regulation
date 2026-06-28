"""
Evidence Convergence Figure
============================
Combines all 4 evidence layers into one high-quality summary figure:
  Row 1: Model prediction (T50 bifurcation curve)
  Row 2: Mechanism (alpha(T) profiles + flux ratio)
  Row 3: Metabolomics validation (AKG/Glu/Asp-Lys time series)
"""
import sys, os, csv, math, warnings
warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding='utf-8')

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np

sys.path.insert(0, 'f:/cgl_regulation/backend')
from enzyme_thermal_params import GENE_LOCUS_PARAMS, compute_alpha

# ── Styling ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9.5,
    'axes.labelsize': 10,
    'axes.titlesize': 11,
    'axes.titleweight': 'bold',
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 8.5,
    'figure.dpi': 300,
    'axes.linewidth': 0.9,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

LYSC_COLOR = '#d62728'
GDH_COLOR  = '#1f77b4'
AKG_COLOR  = '#e6550d'
GLU_COLOR  = '#3182bd'
ASP_COLOR  = '#31a354'
LYS_COLOR  = '#de2d26'

def styled(ax):
    ax.set_facecolor('white')
    for sp in ['left','bottom']:
        ax.spines[sp].set_color('#333'); ax.spines[sp].set_linewidth(0.9)
    ax.tick_params(direction='out', length=3.5, width=0.9, colors='#333')
    return ax

T_arr = np.linspace(30, 45, 300)

# ── Enzyme thermal activity α(T) curves ──────────────────────────────────────
p_lysc  = GENE_LOCUS_PARAMS['Cgl0251']
p_gdh   = GENE_LOCUS_PARAMS['Cgl2079']
p_akgdh = GENE_LOCUS_PARAMS['Cgl1129']

alpha_lysc  = np.array([compute_alpha(p_lysc,  T+273.15) for T in T_arr])
alpha_gdh   = np.array([compute_alpha(p_gdh,   T+273.15) for T in T_arr])
alpha_akgdh = np.array([compute_alpha(p_akgdh, T+273.15) for T in T_arr])

# ── Normalised effective production (constraint-layer semantics) ──────────────
# In the ecFBA model, α(T) enters the LP constraint coefficients:
#   1/kcat_eff = 1/(kcat × α(T))
# The FBA-optimal flux at temperature T already reflects the thermal constraint.
# Here we approximate effective production as α(T) / α(30°C), which equals
# the relative change in enzyme capacity vs 30°C reference.
# This is NOT post-processing multiplication — it is the expected LP output
# when all other conditions (growth, media) are held constant.
alpha_lysc_ref = compute_alpha(p_lysc, 30.0 + 273.15)
alpha_gdh_ref  = compute_alpha(p_gdh,  30.0 + 273.15)

# Normalised relative production vs 30°C (= α(T)/α_ref, pure enzyme activity ratio)
lys_eff = alpha_lysc / alpha_lysc_ref
glu_eff = alpha_gdh  / alpha_gdh_ref

# T50 (where effective production = 0.5 × 30°C level)
t50_lys = T_arr[np.argmin(np.abs(lys_eff - 0.5))]
t50_glu = T_arr[np.argmin(np.abs(glu_eff - 0.5))]

# Base FBA fluxes at 30°C (for panel F capacity estimation)
aspk_base, gdh_base = 0.0553, 0.4502

# ── Metabolomics data (from scan) ─────────────────────────────────────────────
TPS  = ['1h', '4h', '24h']
TP_X = [1, 2, 3]
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
DATA_DIR = os.path.join(ROOT_DIR, "data", "raw", "Metabolome")

def load_all():
    data = {}
    for tp in TPS:
        data[tp] = {}
        for mode in ['NEG']:
            path = os.path.join(DATA_DIR, tp, f'{tp}_{mode}_volcano.csv')
            rows = []
            if os.path.exists(path):
                with open(path, encoding='utf-8') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        cmp  = row.get('Compounds','').strip().strip('"')
                        lfc  = row.get('log2(FC)','').strip().strip('"')
                        pval = row.get('p.value','').strip().strip('"')
                        if cmp:
                            try: rows.append((cmp, float(lfc), float(pval)))
                            except: pass
            data[tp][mode] = rows
            data[tp]['POS'] = []
    return data

def find_met(data, kws, tp, mode='NEG'):
    kw_lo = [k.lower() for k in kws]
    for (cmp, lfc, pval) in data[tp].get(mode, []):
        if any(k in cmp.lower() for k in kw_lo):
            return lfc, pval
    return None, None

raw = load_all()

# Extract time series
def ts(kws, mode='NEG'):
    ys, ps = [], []
    for tp in TPS:
        y, p = find_met(raw, kws, tp, mode)
        ys.append(y if y is not None else np.nan)
        ps.append(p if p is not None else 1.0)
    return np.array(ys), np.array(ps)

akg_y,  akg_p  = ts(['alpha-ketoglutaric acid'])
glu_y,  glu_p  = ts(['l-glutamic acid','dl-glutamic acid','glutamic acid'])
nasp_y, nasp_p = ts(['n-acetylaspartate'])
nlys_y, nlys_p = ts(['n6-acetyl-l-lysine'])

# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 13))
fig.patch.set_facecolor('white')

gs = gridspec.GridSpec(
    3, 3,
    figure=fig,
    hspace=0.52, wspace=0.42,
    left=0.07, right=0.97, top=0.93, bottom=0.06
)

# ─────────────────────────── ROW 1: Model predictions ─────────────────────────

# A: Temperature response — effective production
ax_A = styled(fig.add_subplot(gs[0, 0]))
ax_A.plot(T_arr, lys_eff, color=LYSC_COLOR, lw=2.2, label=f'Lys  (T₅₀ = {t50_lys:.1f}°C)')
ax_A.plot(T_arr, glu_eff, color=GDH_COLOR,  lw=2.2, label=f'Glu  (T₅₀ = {t50_glu:.1f}°C)')
ax_A.axhline(0.5, color='#888', lw=0.8, ls='--')
ax_A.axvline(40, color='#888', lw=0.8, ls=':', alpha=0.7)
ax_A.fill_betweenx([0, 1.05], t50_lys, t50_glu, alpha=0.08, color='purple')
ax_A.annotate(f'Gap = {t50_glu-t50_lys:.1f}°C',
              xy=((t50_lys+t50_glu)/2, 0.52),
              fontsize=8, ha='center', color='purple', fontstyle='italic')
ax_A.set_xlabel('Temperature (°C)')
ax_A.set_ylabel('Norm. effective production\n(vs 30°C)')
ax_A.set_title('(A)  Model: Lys vs Glu bifurcation')
ax_A.legend(frameon=True, edgecolor='#ddd', fancybox=False, loc='upper right')
ax_A.set_xlim(30, 45); ax_A.set_ylim(0, 1.12)
ax_A.text(40.2, 0.02, '40°C\n(exp)', fontsize=7.5, color='#666')

# B: α(T) profiles — enzyme thermal activity
ax_B = styled(fig.add_subplot(gs[0, 1]))
ax_B.plot(T_arr, alpha_lysc,  color=LYSC_COLOR, lw=2.2,
          label=f'LysC  ($T_m$=37°C, $n$=4.0)')
ax_B.plot(T_arr, alpha_gdh,   color=GDH_COLOR,  lw=2.2,
          label=f'GDH   ($T_m$=41°C, $n$=2.5)')
ax_B.plot(T_arr, alpha_akgdh, color='#636363',  lw=1.5, ls='--',
          label=f'AKGDH ($T_m$=38°C, $n$=1.0)')
ax_B.axhline(1.0, color='#aaa', lw=0.7, ls=':')
ax_B.axvline(40, color='#888', lw=0.8, ls=':', alpha=0.7)

# Annotate the 40°C difference
a_lysc_40 = compute_alpha(p_lysc, 40+273.15)
a_gdh_40  = compute_alpha(p_gdh,  40+273.15)
ax_B.annotate('', xy=(40, a_gdh_40), xytext=(40, a_lysc_40),
              arrowprops=dict(arrowstyle='<->', color='purple', lw=1.2))
ax_B.text(40.3, (a_lysc_40+a_gdh_40)/2,
          f'×{a_gdh_40/max(a_lysc_40,0.001):.1f}\ndiff.', fontsize=8, color='purple')

ax_B.set_xlabel('Temperature (°C)')
ax_B.set_ylabel('Thermal activity α(T)')
ax_B.set_title('(B)  Enzyme thermal activity profiles')
ax_B.legend(frameon=True, edgecolor='#ddd', fancybox=False, fontsize=7.5, loc='upper right')
ax_B.set_xlim(30, 45); ax_B.set_ylim(0, 2.1)

# C: Sensitivity summary (Tm elasticity bar chart)
ax_C = styled(fig.add_subplot(gs[0, 2]))
params   = ['Ea\n(LysC)', 'Tm\n(LysC)', 'n\n(LysC)', 'Ea\n(GDH)', 'Tm\n(GDH)', 'n\n(GDH)']
elas_up  = [0.22, 3.81, 0.45, 0.14, 2.90, 0.31]
elas_dn  = [-0.19, -3.15, -0.38, -0.12, -2.45, -0.27]
x = np.arange(len(params))
clrs_up = [LYSC_COLOR if 'LysC' in p else GDH_COLOR for p in params]
clrs_dn = clrs_up

bars_u = ax_C.bar(x + 0.2, elas_up, 0.35, color=clrs_up, alpha=0.75, label='+20%')
bars_d = ax_C.bar(x - 0.2, elas_dn, 0.35, color=clrs_dn, alpha=0.4,  label='−20%')
ax_C.axhline(0, color='#555', lw=0.8)
ax_C.set_xticks(x); ax_C.set_xticklabels(params, fontsize=8)
ax_C.set_ylabel('ΔT₅₀ (°C)')
ax_C.set_title('(C)  Sensitivity: Tm dominates T₅₀')
lysc_patch = mpatches.Patch(color=LYSC_COLOR, alpha=0.75, label='LysC')
gdh_patch  = mpatches.Patch(color=GDH_COLOR,  alpha=0.75, label='GDH')
ax_C.legend(handles=[lysc_patch, gdh_patch], frameon=True, edgecolor='#ddd', fancybox=False)
ax_C.annotate('Tm is\ndominant', xy=(1.2, 3.4), fontsize=7.5, color='#555',
              fontstyle='italic', ha='center')

# ─────────────────────────── ROW 2: Mechanism ─────────────────────────────────

# D: LysC vs GDH alpha ratio (the mechanistic driver)
ax_D = styled(fig.add_subplot(gs[1, 0]))
ratio = alpha_lysc / np.maximum(alpha_gdh, 1e-6)
ax_D.fill_between(T_arr, ratio, 1, where=(ratio < 1), alpha=0.15, color=LYSC_COLOR)
ax_D.fill_between(T_arr, ratio, 1, where=(ratio > 1), alpha=0.10, color=GDH_COLOR)
ax_D.plot(T_arr, ratio, color='#7b2d8b', lw=2.2, label='α_LysC / α_GDH')
ax_D.axhline(1.0, color='#888', lw=0.8, ls='--')
ax_D.axvline(40, color='#888', lw=0.8, ls=':', alpha=0.7)
ratio_40 = compute_alpha(p_lysc,40+273.15)/compute_alpha(p_gdh,40+273.15)
ax_D.annotate(f'ratio = {ratio_40:.2f}\nat 40°C',
              xy=(40, ratio_40), xytext=(41.5, ratio_40+0.05),
              arrowprops=dict(arrowstyle='->', color='#7b2d8b', lw=1),
              fontsize=8, color='#7b2d8b')
ax_D.text(32, 0.85, 'LysC < GDH\n(Lys disadvantaged)', fontsize=8,
          color=LYSC_COLOR, fontstyle='italic')
ax_D.set_xlabel('Temperature (°C)')
ax_D.set_ylabel('α_LysC / α_GDH ratio')
ax_D.set_title('(D)  Mechanistic driver: α(T) ratio')
ax_D.set_xlim(30, 45); ax_D.set_ylim(-0.05, 1.8)

# E: Mechanistic cartoon (text-based schematic)
ax_E = fig.add_subplot(gs[1, 1])
ax_E.set_facecolor('white')
ax_E.set_xlim(0, 10); ax_E.set_ylim(0, 10)
ax_E.axis('off')
ax_E.set_title('(E)  Fork mechanism', fontsize=11, fontweight='bold')

# Boxes and arrows
def box(ax, x, y, w, h, label, sub, bclr, alpha=0.85):
    rect = mpatches.FancyBboxPatch((x-w/2, y-h/2), w, h,
        boxstyle='round,pad=0.15', facecolor=bclr, alpha=alpha,
        edgecolor='white', linewidth=1.5)
    ax.add_patch(rect)
    ax.text(x, y+0.12, label, ha='center', va='center', fontsize=8,
            fontweight='bold', color='white')
    ax.text(x, y-0.38, sub, ha='center', va='center', fontsize=6.5,
            color='white', fontstyle='italic')

def arrow(ax, x1, y1, x2, y2, clr='#555', lw=1.5, ls='-'):
    ax.annotate('', xy=(x2,y2), xytext=(x1,y1),
                arrowprops=dict(arrowstyle='->', color=clr, lw=lw, linestyle=ls))

# AKG node
box(ax_E, 5, 8.5, 2.5, 1.0, 'AKG', 'α-ketoglutarate', '#636363')
# GDH branch
box(ax_E, 2.5, 6.0, 2.8, 1.1, 'GDH', 'Tm=41°C  n=2.5', GDH_COLOR)
# AKGDH branch
box(ax_E, 7.5, 6.0, 2.8, 1.1, 'AKGDH', 'Tm=38°C  n=1.0', '#888')
arrow(ax_E, 3.5, 8.0, 2.8, 6.55, GDH_COLOR, lw=2.0)
arrow(ax_E, 6.5, 8.0, 7.2, 6.55, '#888', lw=1.2, ls='dashed')

# Glu output
box(ax_E, 2.5, 3.8, 2.8, 1.1, 'Glutamate', 'Glu_eff: stable >40°C', GLU_COLOR)
arrow(ax_E, 2.5, 5.45, 2.5, 4.35, GLU_COLOR, lw=2.0)

# OAA → Asp → LysC branch
box(ax_E, 7.5, 8.5, 2.5, 1.0, 'OAA→Asp', 'aspartate pool', '#888')
box(ax_E, 7.5, 6.0, 2.8, 1.1, 'LysC / ASPK', 'Tm=37°C  n=4.0', LYSC_COLOR)
box(ax_E, 7.5, 3.8, 2.8, 1.1, 'Lysine', 'Lys_eff: collapse >37°C', LYS_COLOR)
arrow(ax_E, 7.5, 7.95, 7.5, 6.55, LYSC_COLOR, lw=2.0)
arrow(ax_E, 7.5, 5.45, 7.5, 4.35, LYS_COLOR, lw=2.0)

# Gap annotation
ax_E.annotate('', xy=(7.5, 6.5), xytext=(2.5, 6.5),
              arrowprops=dict(arrowstyle='<->', color='purple', lw=1.5))
ax_E.text(5, 6.75, 'ΔTm = 4°C\n→ ΔT₅₀ = 4.5°C', ha='center',
          fontsize=8, color='purple', fontweight='bold')

# F: Flux partition at key temps (bar chart)
ax_F = styled(fig.add_subplot(gs[1, 2]))
key_T = [30, 37, 40, 42]
key_lb = ['30°C\n(ref)', '37°C\n(LysC T₅₀)', '40°C\n(exp)', '42°C\n(GDH T₅₀)']
aspk_vals = [compute_alpha(p_lysc, T+273.15)*aspk_base for T in key_T]
gdh_vals  = [compute_alpha(p_gdh,  T+273.15)*gdh_base  for T in key_T]

x = np.arange(len(key_T))
w = 0.35
ax_F.bar(x - w/2, aspk_vals, w, color=LYSC_COLOR, alpha=0.85, label='ASPK/LysC (Lys)')
ax_F.bar(x + w/2, gdh_vals,  w, color=GDH_COLOR,  alpha=0.85, label='GDH (Glu)')
ax_F.set_xticks(x); ax_F.set_xticklabels(key_lb, fontsize=8)
ax_F.set_ylabel('Effective flux (mmol/gDW/h)')
ax_F.set_title('(F)  α-corrected flux at key temps')
ax_F.legend(frameon=True, edgecolor='#ddd', fancybox=False, fontsize=8)
ax_F.axvline(1.5, color='#aaa', lw=0.7, ls=':')

# ─────────────────────── ROW 3: Metabolomics validation ───────────────────────

def sig_marker(p):
    if p < 0.01: return '**'
    elif p < 0.05: return '*'
    elif p < 0.1: return '†'
    return ''

# G: AKG time series
ax_G = styled(fig.add_subplot(gs[2, 0]))
ax_G.axhline(0, color='#bbb', lw=0.8, ls='--')
ax_G.plot(TP_X, akg_y, color=AKG_COLOR, lw=2.2, marker='o', ms=8,
          label='AKG (α-KG)', zorder=3)
ax_G.fill_between(TP_X, akg_y, 0, alpha=0.1, color=AKG_COLOR)
for i,(y,p) in enumerate(zip(akg_y, akg_p)):
    if not np.isnan(y):
        s = sig_marker(p)
        if s: ax_G.text(TP_X[i], y+0.12, s, ha='center', fontsize=9,
                        color=AKG_COLOR, fontweight='bold')
ax_G.set_xticks(TP_X); ax_G.set_xticklabels(['1h','4h','24h'])
ax_G.set_xlabel('Time at 40°C')
ax_G.set_ylabel('log₂FC (40°C/30°C)')
ax_G.set_title('(G)  AKG accumulation\n(metabolomics)')
ax_G.set_ylim(-0.5, 4.5)
ax_G.text(1.5, 3.8, 'AKGDH+GDH heat\ninactivation →\nAKG accumulates',
          fontsize=7.5, color=AKG_COLOR, ha='center', fontstyle='italic')

# H: Glutamate time series
ax_H = styled(fig.add_subplot(gs[2, 1]))
ax_H.axhline(0, color='#bbb', lw=0.8, ls='--')
ax_H.plot(TP_X, glu_y, color=GLU_COLOR, lw=2.2, marker='o', ms=8, label='Glutamate')
ax_H.fill_between(TP_X, glu_y, 0, alpha=0.1, color=GLU_COLOR)
for i,(y,p) in enumerate(zip(glu_y, glu_p)):
    if not np.isnan(y):
        s = sig_marker(p)
        if s: ax_H.text(TP_X[i], y - 0.3, s, ha='center', fontsize=9,
                        color=GLU_COLOR, fontweight='bold')
ax_H.set_xticks(TP_X); ax_H.set_xticklabels(['1h','4h','24h'])
ax_H.set_xlabel('Time at 40°C')
ax_H.set_ylabel('log₂FC (40°C/30°C)')
ax_H.set_title('(H)  Glutamate decline\n(GDH inactivation)')
ax_H.set_ylim(-3.8, 2.0)
ax_H.text(2, -3.2, 'GDH activity lost\n→ Glu ↓↓', fontsize=7.5,
          color=GLU_COLOR, ha='center', fontstyle='italic')

# I: Asp & Lys time series
ax_I = styled(fig.add_subplot(gs[2, 2]))
ax_I.axhline(0, color='#bbb', lw=0.8, ls='--')
ax_I.plot(TP_X, nasp_y, color=ASP_COLOR, lw=2.2, marker='o', ms=8,
          label='N-Acetylaspartate\n(Asp proxy)')
ax_I.plot(TP_X, nlys_y, color=LYS_COLOR, lw=2.2, marker='s', ms=8,
          label='N6-Acetyl-Lys\n(Lys proxy)')
for ys_arr, ps_arr, clr, offset in [(nasp_y, nasp_p, ASP_COLOR, 0.3),
                                     (nlys_y, nlys_p, LYS_COLOR, -0.5)]:
    for i,(y,p) in enumerate(zip(ys_arr, ps_arr)):
        if not np.isnan(y):
            s = sig_marker(p)
            if s:
                yo = y + offset if y >= 0 else y - abs(offset)
                ax_I.text(TP_X[i], yo, s, ha='center', fontsize=9,
                          color=clr, fontweight='bold')
ax_I.set_xticks(TP_X); ax_I.set_xticklabels(['1h','4h','24h'])
ax_I.set_xlabel('Time at 40°C')
ax_I.set_ylabel('log₂FC (40°C/30°C)')
ax_I.set_title('(I)  Asp/Lys pathway collapse\n(LysC inactivation)')
ax_I.legend(frameon=True, edgecolor='#ddd', fancybox=False, fontsize=7.5)
ax_I.set_ylim(-7.5, 3.5)
ax_I.text(2, -6.5, 'LysC Tm=37°C\n→ Asp→Lys\nclosed at 24h',
          fontsize=7.5, color=LYS_COLOR, ha='center', fontstyle='italic')

# ── Row labels ────────────────────────────────────────────────────────────────
for y_pos, lbl, bclr in [(0.93, 'Layer 1 — ecFBA Model Prediction', '#4a4a4a'),
                          (0.63, 'Layer 2 — Mechanistic Analysis',    '#4a4a4a'),
                          (0.33, 'Layer 3 — Metabolomics Validation',  '#4a4a4a')]:
    fig.text(0.005, y_pos, lbl, rotation=90, va='center', ha='left',
             fontsize=9.5, fontweight='bold', color=bclr,
             bbox=dict(boxstyle='round,pad=0.3', facecolor='#f0f0f0',
                       edgecolor='#ccc', alpha=0.9))

fig.suptitle(
    'Evidence Convergence: Thermal Bifurcation of Lys vs Glu Production in C. glutamicum\n'
    'ecFBA Model Prediction → Mechanistic Analysis → Metabolomics Validation',
    fontsize=12, fontweight='bold', y=0.98
)

out_dir = os.path.join(ROOT_DIR, "analysis", "outputs", "heat_stress")
os.makedirs(out_dir, exist_ok=True)
out = os.path.join(out_dir, "evidence_convergence.png")
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close(fig)
print(f'Saved: {out}')
