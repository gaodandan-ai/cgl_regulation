"""
Flux Partition Analysis — Mechanistic Fork Diagram
===================================================
AKG fork:  AKGDH (TCA, Cgl1129 Tm=38.1C) vs GDH (Glu, Cgl2079 Tm=41.0C)
Asp fork:  ASPTA (Asp pool) vs ASPK/LysC (Lys, Cgl0251 Tm=37.0C)

For each temperature, runs pFBA with biomass objective to extract internal
metabolic fluxes. Then applies alpha(T) correction to get thermal-adjusted
flux estimates. Generates mechanistic bifurcation figure.
"""
import sys, os, copy, math, json
import warnings; warnings.filterwarnings("ignore")
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(ROOT_DIR, "backend"))
sys.path.insert(0, ROOT_DIR)
import cobra
from cobra.flux_analysis import pfba
from enzyme_thermal_params import get_params, compute_alpha, GENE_LOCUS_PARAMS

MODEL_PATH   = os.path.join(ROOT_DIR, "data", "reference", "model", "ecCGL1-main", "ecCGL1-main", "model", "iCW773_irr_enz_constraint.json")
OUT_DIR      = os.path.join(ROOT_DIR, "analysis", "outputs", "heat_stress")
os.makedirs(OUT_DIR, exist_ok=True)
PROTEIN_POOL = 0.129
R_GAS        = 8.314
T_REF_K      = 303.15

print("Loading model...", end=" ", flush=True)
BASE_MODEL = cobra.io.load_json_model(MODEL_PATH)
with open(MODEL_PATH) as f:
    d = json.load(f)

# Load DLKcat predictions
predictions_path = os.path.join(ROOT_DIR, "data", "reference", "dlkcat_predicted_kcat.json")
predictions = {}
if os.path.exists(predictions_path):
    with open(predictions_path, "r", encoding="utf-8") as pf:
        predictions = json.load(pf)

KCAT_MAP = {}
DEFAULT_KCAT = 7398.8133918117555

for r in d["reactions"]:
    rxn_id = r["id"]
    kcat_mw = float(r.get("kcat_MW", 0)) if r.get("kcat_MW") else 0.0
    kcat = float(r.get("kcat", 0)) if r.get("kcat") else 0.0
    
    # Calculate MW of the enzyme for this reaction
    mw = 50000.0 # fallback default
    if kcat > 0 and kcat_mw > 0:
        mw = (kcat * 3600 * 1000) / kcat_mw
        
    # Check if we have a prediction to override
    if rxn_id in predictions:
        pred_info = predictions[rxn_id]
        if pred_info["source"] == "dlkcat_prediction":
            # Only override if original kcat is missing, 0, or default fallback
            if kcat <= 0 or abs(kcat - DEFAULT_KCAT) < 1e-3:
                kcat = float(pred_info["kcat"])
                kcat_mw = (kcat * 3600 * 1000) / mw
                
    if kcat_mw > 0:
        KCAT_MAP[rxn_id] = kcat_mw

print(f"Loaded {len(KCAT_MAP)} enzyme constraints (including DLKcat predictions). Done.\n")

BIO_ID  = None
LYS_ID  = "EX_lys_L_e"
GLU_ID  = "EX_glu_L_e"
for bid in ["CG_biomass_cgl_ATCC13032","BIOMASS_Cgl_ATCC13032","Growth"]:
    if bid in BASE_MODEL.reactions: BIO_ID = bid; break

# Fork reaction IDs
ASPK_ID  = "ASPK"    # LysC (Cgl0251) - aspartate kinase (Lys pathway entry)
DHDPS_ID = "DHDPS"   # dihydrodipicolinate synthase (Lys pathway step 3)
GDH_ID   = "GLUDy"   # GDH (Cgl2079) - glutamate dehydrogenase
AKGDH_ID = "AKGDH"   # AKGDH (Cgl1129) - 2-oxoglutarate dehydrogenase (TCA)
ASPTA_ID = "ASPTA"   # aspartate transaminase (Asp pool source)

def hsp_frac(temp):
    dt = max(0.0, temp - 30.0)
    return 0.20 * (dt**3) / (11.0**3 + dt**3)

def build_model(temp):
    model = BASE_MODEL.copy()
    T_K   = temp + 273.15
    for rxn in model.exchanges:
        if rxn.lower_bound < 0: rxn.lower_bound = 0.0
    for nid in ["EX_glc__D_e","EX_glc_e","EX_glucose_e"]:
        if nid in model.reactions: model.reactions.get_by_id(nid).lower_bound = -10.0
    for nid in ["EX_o2_e","EX_nh4_e","EX_pi_e","EX_so4_e","EX_k_e",
                "EX_mg2_e","EX_h2o_e","EX_h_e","EX_fe2_e","EX_fe3_e"]:
        if nid in model.reactions: model.reactions.get_by_id(nid).lower_bound = -1000.0
    for eid in [LYS_ID, GLU_ID]:
        if eid in model.reactions:
            model.reactions.get_by_id(eid).lower_bound = 0.0
            model.reactions.get_by_id(eid).upper_bound = 1000.0

    f_hsp      = hsp_frac(temp)
    pool_avail = PROTEIN_POOL * (1.0 - f_hsp)
    coefs = {}
    # FIX Issue 5: use van't Hoff formula for η(ΔG) — same as simulation.py
    DH_RXN_MAP = {"ICDHyr": 8400.0, "MDH": 29000.0, "PGI": -2500.0}
    for rxn in model.reactions:
        km = KCAT_MAP.get(rxn.id)
        if not km: continue
        gene_loci = [g.id.replace("g_","").replace("gene_","") for g in rxn.genes]
        p = get_params(rxn.id, gene_loci)
        try:    alpha = compute_alpha(p, T_K)
        except: alpha = 1.0
        dG_ref = 3000.0
        if "ICDHyr" in rxn.id: dG_ref = 1000.0
        if "MDH"   in rxn.id and "UAMDH" not in rxn.id: dG_ref = 800.0
        if "PGI"   in rxn.id: dG_ref = 1500.0
        dH_rxn = next((DH_RXN_MAP[k] for k in DH_RXN_MAP if k in rxn.id), -10000.0)
        dG_T = max(100.0, dG_ref + dH_rxn * (1.0 - T_K / T_REF_K))
        eta  = max(0.01, math.tanh(dG_T / (2.0 * R_GAS * T_K)))
        coefs[rxn.forward_variable] = 1.0 / (km * alpha * eta)

    pool_con = model.problem.Constraint(0, lb=0, ub=pool_avail, name="pool")
    model.add_cons_vars(pool_con)
    model.solver.update()
    pool_con.set_linear_coefficients(coefs)
    return model

def run_temp(temp):
    model = build_model(temp)
    T_K   = temp + 273.15

    # pFBA with growth objective -> internal fluxes
    model.objective = BIO_ID
    try:
        sol = pfba(model, fraction_of_optimum=1.0)
        if sol.status != "optimal": return None
        growth  = float(sol.fluxes.get(BIO_ID, 0.0))
        aspk    = abs(float(sol.fluxes.get(ASPK_ID,  0.0)))
        dhdps   = abs(float(sol.fluxes.get(DHDPS_ID, 0.0)))
        gdh     = abs(float(sol.fluxes.get(GDH_ID,   0.0)))
        akgdh   = abs(float(sol.fluxes.get(AKGDH_ID, 0.0)))
        aspta   = abs(float(sol.fluxes.get(ASPTA_ID, 0.0)))
    except Exception as e:
        print(f"  pFBA failed at {temp}C: {e}")
        return None

    # Max Lys (95% growth floor)
    m2 = model.copy()
    if BIO_ID in m2.reactions: m2.reactions.get_by_id(BIO_ID).lower_bound = 0.95*growth
    m2.objective = LYS_ID
    try:
        s2 = m2.optimize()
        lys_max = max(0.0, float(s2.fluxes.get(LYS_ID,0.0))) if s2.status=="optimal" else 0.0
    except: lys_max = 0.0

    # Max Glu (95% growth floor)
    m3 = model.copy()
    if BIO_ID in m3.reactions: m3.reactions.get_by_id(BIO_ID).lower_bound = 0.95*growth
    m3.objective = GLU_ID
    try:
        s3 = m3.optimize()
        glu_max = max(0.0, float(s3.fluxes.get(GLU_ID,0.0))) if s3.status=="optimal" else 0.0
    except: glu_max = 0.0

    # α(T) is already embedded in the LP constraint coefficients (build_model L80):
    #   coefs[rxn.forward_variable] = 1.0 / (km * alpha * eta)
    # The LP-optimal fluxes lys_max / glu_max already reflect thermal enzyme capacity.
    # Report α values here for diagnostics / figure annotation only.
    # DO NOT multiply lys_max/glu_max by α again — that would cause α² distortion.
    a_lysc  = compute_alpha(GENE_LOCUS_PARAMS["Cgl0251"], T_K)
    a_gdh   = compute_alpha(GENE_LOCUS_PARAMS["Cgl2079"], T_K)
    a_akgdh = compute_alpha(GENE_LOCUS_PARAMS["Cgl1129"], T_K)

    return {
        "growth":      growth,
        "aspk":        aspk,
        "dhdps":       dhdps,
        "gdh":         gdh,
        "akgdh":       akgdh,
        "aspta":       aspta,
        # Constraint-layer FBA fluxes — α(T) already embedded in LP, no post-processing
        "lys_max":     lys_max,
        "glu_max":     glu_max,
        "lys_fba":     lys_max,   # identical: kept for backward compatibility
        "glu_fba":     glu_max,
        # α values for transparency / figure annotation (diagnostic only)
        "alpha_lysc":  a_lysc,
        "alpha_gdh":   a_gdh,
        "alpha_akgdh": a_akgdh,
        "hsp":         hsp_frac(temp),
    }

TEMPS = [30 + i*0.5 for i in range(31)]
print(f"\nRunning fork analysis sweep ({len(TEMPS)} temps)...")
print(f"{'Temp':>6}  {'Growth':>7}  {'ASPK':>7}  {'DHDPS':>7}  {'GDH':>7}  {'AKGDH':>7}  {'Lys_eff':>8}  {'Glu_eff':>8}")
print("-"*75)
results = {}
for T in TEMPS:
    r = run_temp(T)
    if r:
        results[T] = r
        print(f"  {T:5.1f}C  {r['growth']:7.4f}  {r['aspk']:7.4f}  {r['dhdps']:7.4f}  "
              f"{r['gdh']:7.4f}  {r['akgdh']:7.4f}  {r['lys_max']:8.4f}  {r['glu_max']:8.4f}")
    else:
        results[T] = None
        print(f"  {T:5.1f}C  INFEASIBLE")

# Interpolate gaps
Tl = sorted(results)
for i,T in enumerate(Tl):
    if results[T] is None:
        pT=next((Tl[j] for j in range(i-1,-1,-1) if results.get(Tl[j])),None)
        nT=next((Tl[j] for j in range(i+1,len(Tl)) if results.get(Tl[j])),None)
        if pT and nT:
            f=(T-pT)/(nT-pT)
            results[T]={k:results[pT][k]+f*(results[nT][k]-results[pT][k]) for k in results[pT]}
            print(f"  {T:5.1f}C  interpolated")

Tv = [T for T in TEMPS if results.get(T)]

# Reference at 30°C
ref = results[30.0]
aspk_ref  = max(ref["aspk"],  1e-9)
dhdps_ref = max(ref["dhdps"], 1e-9)
gdh_ref   = max(ref["gdh"],   1e-9)
akgdh_ref = max(ref["akgdh"], 1e-9)
lys_ref   = max(ref["lys_max"],1e-9)
glu_ref   = max(ref["glu_max"],1e-9)

print(f"\nRef (30C): ASPK={aspk_ref:.4f}  DHDPS={dhdps_ref:.4f}  GDH={gdh_ref:.4f}  AKGDH={akgdh_ref:.4f}")

# ── Figures ────────────────────────────────────────────────────────────────────
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

plt.rcParams.update({
    "font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "font.size":9.5,"axes.labelsize":10,"axes.titlesize":10,"axes.titleweight":"bold",
    "xtick.labelsize":8.5,"ytick.labelsize":8.5,"legend.fontsize":8,
    "figure.dpi":300,"axes.linewidth":0.8,
    "axes.spines.top":False,"axes.spines.right":False,
})

def styled(ax):
    ax.set_facecolor("white")
    ax.tick_params(direction="out",length=3,width=0.8,colors="black")
    for sp in ["left","bottom"]:
        ax.spines[sp].set_color("black"); ax.spines[sp].set_linewidth(0.8)
    return ax

CLR_LYS   = "#d6604d"   # warm red  -> Lys/ASPK/LysC
CLR_GLU   = "#4393c3"   # blue      -> Glu/GDH
CLR_AKGDH = "#2ca02c"   # green     -> AKGDH (TCA, competes with GDH)
CLR_DHDPS = "#ff7f0e"   # orange    -> DHDPS (downstream Lys)
CLR_SHADE = "#fef0e7"   # light salmon shading for 37-41°C window

# Derived arrays
aspk_n   = [results[T]["aspk"]/aspk_ref     for T in Tv]
dhdps_n  = [results[T]["dhdps"]/dhdps_ref   for T in Tv]
gdh_n    = [results[T]["gdh"]/gdh_ref       for T in Tv]
akgdh_n  = [results[T]["akgdh"]/akgdh_ref   for T in Tv]
lys_n    = [results[T]["lys_max"]/lys_ref   for T in Tv]
glu_n    = [results[T]["glu_max"]/glu_ref   for T in Tv]
al_lysc  = [results[T]["alpha_lysc"]        for T in Tv]
al_gdh   = [results[T]["alpha_gdh"]         for T in Tv]
al_akgdh = [results[T]["alpha_akgdh"]       for T in Tv]
# AKG partitioning: GDH/(GDH+AKGDH)
akg_part = []
for T in Tv:
    r = results[T]
    tot = r["gdh"] + r["akgdh"]
    akg_part.append(r["gdh"]/tot if tot>1e-9 else 0.5)
# Asp partitioning: ASPK/(ASPK+ASPTA/2)
asp_part = []
for T in Tv:
    r = results[T]
    tot = r["aspk"] + r["aspta"]
    asp_part.append(r["aspk"]/tot if tot>1e-9 else 0.5)

fig = plt.figure(figsize=(16, 12))
fig.patch.set_facecolor("white")
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.50, wspace=0.40)

def shade_window(ax, lo=37.0, hi=41.0):
    ax.axvspan(lo, hi, alpha=0.10, color="#f4a582", zorder=0)
    ax.axvline(37.0, color=CLR_LYS,   lw=0.9, ls=":", alpha=0.7)
    ax.axvline(41.0, color=CLR_GLU,   lw=0.9, ls=":", alpha=0.7)

# Panel A: ASPK vs GDH flux (normalized)
ax_A = styled(fig.add_subplot(gs[0,0]))
ax_A.plot(Tv, aspk_n,  color=CLR_LYS,   lw=2.0, label="ASPK (LysC) flux")
ax_A.plot(Tv, gdh_n,   color=CLR_GLU,   lw=2.0, label="GDH flux")
ax_A.axhline(1.0, color="#999", lw=0.6, ls="--")
shade_window(ax_A)
ax_A.set_xlim(29.5,45.5); ax_A.set_ylim(-0.05, 2.0)
ax_A.set_xlabel("Temperature (°C)"); ax_A.set_ylabel("Norm. flux (vs 30°C)")
ax_A.set_title("(A)  ASPK vs GDH — flux divergence")
ax_A.legend(frameon=True, edgecolor="#ccc", fancybox=False)
# annotation: arrow showing gap
ax_A.annotate("", xy=(41.0,0.5), xytext=(37.0,0.5),
              arrowprops=dict(arrowstyle="<->",color="#555",lw=1.2))
ax_A.text(39.0, 0.55, "4.0°C $T_m$ gap", ha="center", fontsize=8, color="#555")

# Panel B: DHDPS vs AKGDH flux
ax_B = styled(fig.add_subplot(gs[0,1]))
ax_B.plot(Tv, dhdps_n,  color=CLR_DHDPS, lw=2.0, label="DHDPS (Lys path step 3)")
ax_B.plot(Tv, akgdh_n,  color=CLR_AKGDH, lw=2.0, label="AKGDH (TCA, competes Glu)")
ax_B.axhline(1.0, color="#999", lw=0.6, ls="--")
shade_window(ax_B)
ax_B.set_xlim(29.5,45.5); ax_B.set_ylim(-0.05, 2.0)
ax_B.set_xlabel("Temperature (°C)"); ax_B.set_ylabel("Norm. flux (vs 30°C)")
ax_B.set_title("(B)  DHDPS & AKGDH — downstream fork")
ax_B.legend(frameon=True, edgecolor="#ccc", fancybox=False)

# Panel C: alpha(T) for LysC, GDH, AKGDH
ax_C = styled(fig.add_subplot(gs[0,2]))
ax_C.plot(Tv, al_lysc,  color=CLR_LYS,   lw=2.0, ls="-",  label="$\\alpha_{LysC}$ ($T_m$=37°C, $n$=4.0)")
ax_C.plot(Tv, al_gdh,   color=CLR_GLU,   lw=2.0, ls="-",  label="$\\alpha_{GDH}$ ($T_m$=41°C, $n$=2.5)")
ax_C.plot(Tv, al_akgdh, color=CLR_AKGDH, lw=1.5, ls="--", label="$\\alpha_{AKGDH}$ ($T_m$=38.1°C, $n$=1)")
ax_C.axhline(1.0, color="#999", lw=0.6, ls="--")
ax_C.axhline(0.5, color="#ccc", lw=0.6, ls=":")
shade_window(ax_C)
ax_C.set_xlim(29.5,45.5); ax_C.set_ylim(-0.05, 2.1)
ax_C.set_xlabel("Temperature (°C)"); ax_C.set_ylabel("Thermal activity $\\alpha$(T)")
ax_C.set_title("(C)  Enzyme thermal activity profiles")
ax_C.legend(frameon=True, edgecolor="#ccc", fancybox=False)

# Panel D: AKG partitioning GDH/(GDH+AKGDH)
ax_D = styled(fig.add_subplot(gs[1,0]))
ax_D.fill_between(Tv, akg_part, 0.5, where=[v>0.5 for v in akg_part],
                  alpha=0.25, color=CLR_GLU, label="GDH dominates")
ax_D.fill_between(Tv, akg_part, 0.5, where=[v<=0.5 for v in akg_part],
                  alpha=0.25, color=CLR_AKGDH, label="AKGDH dominates")
ax_D.plot(Tv, akg_part, color=CLR_GLU, lw=2.0)
ax_D.axhline(0.5, color="#888", lw=0.8, ls="--")
shade_window(ax_D)
ax_D.set_xlim(29.5,45.5); ax_D.set_ylim(0,1.05)
ax_D.set_xlabel("Temperature (°C)")
ax_D.set_ylabel("GDH/(GDH+AKGDH)")
ax_D.set_title("(D)  AKG partitioning — Glu vs TCA")
ax_D.legend(frameon=True, edgecolor="#ccc", fancybox=False, loc="lower left")
ax_D.text(40.5, 0.85, "More AKG\nto Glu", color=CLR_GLU, fontsize=8, ha="center")
ax_D.text(33.0, 0.15, "More AKG\nto TCA", color=CLR_AKGDH, fontsize=8, ha="center")

# Panel E: Asp partitioning ASPK/(ASPK+ASPTA)
ax_E = styled(fig.add_subplot(gs[1,1]))
ax_E.fill_between(Tv, asp_part, 0, alpha=0.20, color=CLR_LYS)
ax_E.plot(Tv, asp_part, color=CLR_LYS, lw=2.0, label="ASPK/(ASPK+ASPTA)")
ax_E.axhline(0.5, color="#888", lw=0.8, ls="--")
shade_window(ax_E)
ax_E.set_xlim(29.5,45.5); ax_E.set_ylim(0,1.05)
ax_E.set_xlabel("Temperature (°C)")
ax_E.set_ylabel("Fraction to Lys pathway")
ax_E.set_title("(E)  Asp partitioning — Lys vs other")
ax_E.legend(frameon=True, edgecolor="#ccc", fancybox=False)

# Panel F: Lys_eff vs Glu_eff (main output, matches bifurcation analysis)
ax_F = styled(fig.add_subplot(gs[1,2]))
ax_F.plot(Tv, lys_n, color=CLR_LYS, lw=2.0, label="Lys$_{eff}$")
ax_F.plot(Tv, glu_n, color=CLR_GLU, lw=2.0, label="Glu$_{eff}$")
ax_F.axhline(0.5, color="#aaa", lw=0.7, ls="--")
shade_window(ax_F)
t50_lys = t50_glu = None
for i,T in enumerate(Tv[1:],1):
    if lys_n[i-1]>=0.5 and lys_n[i]<0.5 and t50_lys is None:
        f=(lys_n[i-1]-0.5)/max(lys_n[i-1]-lys_n[i],1e-9); t50_lys=Tv[i-1]+f*(Tv[i]-Tv[i-1])
    if glu_n[i-1]>=0.5 and glu_n[i]<0.5 and t50_glu is None:
        f=(glu_n[i-1]-0.5)/max(glu_n[i-1]-glu_n[i],1e-9); t50_glu=Tv[i-1]+f*(Tv[i]-Tv[i-1])
if t50_lys: ax_F.axvline(t50_lys, color=CLR_LYS, lw=1.2, ls="--", alpha=0.7)
if t50_glu: ax_F.axvline(t50_glu, color=CLR_GLU, lw=1.2, ls="--", alpha=0.7)
ax_F.set_xlim(29.5,45.5); ax_F.set_ylim(-0.05,2.0)
ax_F.set_xlabel("Temperature (°C)"); ax_F.set_ylabel("Norm. production (vs 30°C)")
gap_str = f"{t50_glu-t50_lys:.1f}°C" if t50_lys and t50_glu else "N/A"
ax_F.set_title(f"(F)  Lys vs Glu  (Gap = {gap_str})")
ax_F.legend(frameon=True, edgecolor="#ccc", fancybox=False)

# Panel G-I: Stacked bar at key temperatures (30/37/40/42°C)
key_Ts = [30.0, 37.0, 40.0, 42.0]
labels = ["30°C\n(baseline)", "37°C\n(Lys T50)", "40°C\n(heat)", "42°C\n(Glu T50)"]

# Panel G: flux magnitudes at key Ts
ax_G = styled(fig.add_subplot(gs[2,0]))
x  = np.arange(len(key_Ts))
w  = 0.2
for i_r, (rname, clr, label) in enumerate([("aspk",CLR_LYS,"ASPK"),("dhdps",CLR_DHDPS,"DHDPS"),
                                            ("gdh",CLR_GLU,"GDH"),("akgdh",CLR_AKGDH,"AKGDH")]):
    vals = []
    for T in key_Ts:
        r = results.get(T) or {}
        vals.append(r.get(rname,0))
    ax_G.bar(x + i_r*w - 1.5*w, vals, w, color=clr, alpha=0.85, label=label)
ax_G.set_xticks(x); ax_G.set_xticklabels(labels, fontsize=8)
ax_G.set_ylabel("Flux (mmol/gDW/h)")
ax_G.set_title("(G)  Fork fluxes at key temperatures")
ax_G.legend(frameon=True, edgecolor="#ccc", fancybox=False, fontsize=7.5)

# Panel H: Ratio ASPK/GDH shows divergence
ax_H = styled(fig.add_subplot(gs[2,1]))
ratio_ag = []
for T in Tv:
    r = results[T]
    g = max(r["gdh"],1e-9)
    ratio_ag.append(r["aspk"]/g)
ax_H.plot(Tv, ratio_ag, color="#7b2d8b", lw=2.0)
ax_H.axhline(ratio_ag[0], color="#ccc", lw=0.7, ls="--", label=f"30°C baseline ({ratio_ag[0]:.2f})")
ax_H.axhline(1.0, color="#aaa", lw=0.6, ls=":")
shade_window(ax_H)
ax_H.set_xlim(29.5,45.5)
ax_H.set_xlabel("Temperature (°C)")
ax_H.set_ylabel("ASPK flux / GDH flux")
ax_H.set_title("(H)  Lys/Glu enzyme flux ratio")
ax_H.legend(frameon=True, edgecolor="#ccc", fancybox=False)
# Annotation: where ratio crosses 1
for i in range(1,len(Tv)):
    if ratio_ag[i-1]>1 and ratio_ag[i]<=1:
        ax_H.axvline(Tv[i], color="#7b2d8b", lw=0.8, ls=":", alpha=0.7)
        ax_H.text(Tv[i]+0.3, max(ratio_ag)*0.8, f"{Tv[i]:.1f}°C", fontsize=7.5, color="#7b2d8b")

# Panel I: Mechanistic schematic (text + arrows via annotations)
ax_I = fig.add_subplot(gs[2,2])
ax_I.set_xlim(0,10); ax_I.set_ylim(0,10)
ax_I.set_aspect("equal"); ax_I.axis("off")
ax_I.set_title("(I)  Mechanistic fork diagram", fontsize=10, fontweight="bold", pad=6)

# Draw metabolic fork schematic
style_box = dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#aaa", lw=0.8)
style_hot = dict(boxstyle="round,pad=0.3", facecolor="#fee8d6", edgecolor=CLR_LYS, lw=1.2)
style_cool= dict(boxstyle="round,pad=0.3", facecolor="#deeaf7", edgecolor=CLR_GLU, lw=1.2)
style_tca = dict(boxstyle="round,pad=0.3", facecolor="#e8f5e9", edgecolor=CLR_AKGDH, lw=1.2)

# OAA -> Asp -> [fork]
ax_I.text(5, 9.0, "OAA", ha="center", va="center", fontsize=8.5, fontweight="bold", bbox=style_box)
ax_I.annotate("", xy=(5,8.0), xytext=(5,8.7), arrowprops=dict(arrowstyle="->",color="#555",lw=1))
ax_I.text(5, 7.6, "Asp (ASPTA)", ha="center", va="center", fontsize=8.5, bbox=style_box)
ax_I.text(1, 9.0, "AKG", ha="center", va="center", fontsize=8.5, fontweight="bold", bbox=style_box)
# AKG fork
ax_I.annotate("", xy=(0.5,7.0), xytext=(1,8.6), arrowprops=dict(arrowstyle="->",color=CLR_GLU,lw=1.5))
ax_I.annotate("", xy=(2.5,7.2), xytext=(1.5,8.6), arrowprops=dict(arrowstyle="->",color=CLR_AKGDH,lw=1.5))
ax_I.text(0, 6.6, "GDH\n$T_m$=41°C", ha="center", va="center", fontsize=8, color=CLR_GLU, fontweight="bold",
          bbox=style_cool)
ax_I.text(3.0, 6.8, "AKGDH\n$T_m$=38.1°C", ha="center", va="center", fontsize=8, color=CLR_AKGDH,
          bbox=style_tca)
ax_I.annotate("", xy=(0,5.5), xytext=(0,6.2), arrowprops=dict(arrowstyle="->",color=CLR_GLU,lw=1.2))
ax_I.text(0, 5.1, "Glutamate", ha="center", va="center", fontsize=8, bbox=style_cool)
ax_I.annotate("", xy=(3,5.8), xytext=(3,6.4), arrowprops=dict(arrowstyle="->",color=CLR_AKGDH,lw=1.2))
ax_I.text(3, 5.4, "Succinyl-CoA\n(TCA)", ha="center", va="center", fontsize=7.5, bbox=style_tca)
# Asp -> LysC fork
ax_I.annotate("", xy=(5,6.6), xytext=(5,7.2), arrowprops=dict(arrowstyle="->",color=CLR_LYS,lw=1.5))
ax_I.text(5, 6.2, "ASPK/LysC\n$T_m$=37°C, $n$=4", ha="center", va="center", fontsize=8,
          color=CLR_LYS, fontweight="bold", bbox=style_hot)
ax_I.annotate("", xy=(5,5.2), xytext=(5,5.8), arrowprops=dict(arrowstyle="->",color=CLR_LYS,lw=1.2))
ax_I.text(5, 4.8, "Asp-phosphate", ha="center", va="center", fontsize=7.5, bbox=style_box)
ax_I.annotate("", xy=(5,4.0), xytext=(5,4.4), arrowprops=dict(arrowstyle="->",color=CLR_DHDPS,lw=1.2))
ax_I.text(5, 3.6, "DHDPS", ha="center", va="center", fontsize=8, color=CLR_DHDPS, bbox=style_box)
ax_I.annotate("", xy=(5,2.8), xytext=(5,3.2), arrowprops=dict(arrowstyle="->",color=CLR_DHDPS,lw=1.2))
ax_I.text(5, 2.4, "Lysine", ha="center", va="center", fontsize=8.5, fontweight="bold", bbox=style_hot)
# Delta Tm annotation
ax_I.annotate("", xy=(3.5,6.2), xytext=(3.5,5.6),
              arrowprops=dict(arrowstyle="<->",color="#333",lw=1.2))
ax_I.text(8.5, 5.0, "$\\Delta T_m$=4.0°C\nexplains 4.5°C\nbifurcation gap", ha="center",
          va="center", fontsize=8, color="#333",
          bbox=dict(boxstyle="round,pad=0.4",facecolor="#fffff0",edgecolor="#bbb",lw=0.8))
ax_I.annotate("", xy=(6.5,5.8), xytext=(8.2,5.8),
              arrowprops=dict(arrowstyle="->",color="#555",lw=0.8))

fig.suptitle("Mechanistic flux partitioning — $\\it{C. glutamicum}$ LysC ($T_m$=37.0°C) vs GDH ($T_m$=41.0°C)",
             fontsize=11, fontweight="bold", y=1.005)

out = os.path.join(OUT_DIR, "flux_partition_mechanism.png")
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
plt.close(fig)
print(f"\nSaved: {out}")
print("Done.")
