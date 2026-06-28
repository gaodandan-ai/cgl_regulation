"""
Sensitivity Analysis for ecCGL1 heat stress model
Ea / Tm / n +/- 20% perturbations on LysC and GDH
Reuses the full model infrastructure from temp_sweep_bifurcation.py
"""

import sys, os, copy, math
import warnings
warnings.filterwarnings("ignore")

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(ROOT_DIR, "backend"))
sys.path.insert(0, ROOT_DIR)
import cobra
from cobra.flux_analysis import pfba
import json
from enzyme_thermal_params import get_params, compute_alpha, GENE_LOCUS_PARAMS

MODEL_PATH   = os.path.join(ROOT_DIR, "data", "reference", "model", "ecCGL1-main", "ecCGL1-main", "model", "iCW773_irr_enz_constraint.json")
OUT_DIR      = os.path.join(ROOT_DIR, "analysis", "outputs", "heat_stress")
os.makedirs(OUT_DIR, exist_ok=True)
PROTEIN_POOL = 0.129
R_GAS        = 8.314
T_REF_K      = 303.15
PERTURB      = 0.20

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
print(f"BIO_ID: {BIO_ID}")

LYSC_LOCUS = "Cgl0251"
GDH_LOCUS  = "Cgl2079"

# Base params
p_lysc_base = copy.deepcopy(GENE_LOCUS_PARAMS.get(LYSC_LOCUS))
p_gdh_base  = copy.deepcopy(GENE_LOCUS_PARAMS.get(GDH_LOCUS))
n_lysc_base = p_lysc_base.get("hill_n", 1.0)
n_gdh_base  = p_gdh_base.get("hill_n", 1.0)
Tm_lysc_C   = p_lysc_base["H_d"]/p_lysc_base["S_d"]-273.15
Tm_gdh_C    = p_gdh_base["H_d"]/p_gdh_base["S_d"]-273.15

print(f"LysC: Ea={p_lysc_base['E_a']/1000:.1f}kJ/mol  Tm={Tm_lysc_C:.1f}C  n={n_lysc_base}")
print(f"GDH:  Ea={p_gdh_base['E_a']/1000:.1f}kJ/mol   Tm={Tm_gdh_C:.1f}C  n={n_gdh_base}")

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
    # FIX Issue 5: use van't Hoff formula for η(ΔG) — consistent with simulation.py
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
        if "MDH" in rxn.id and "UAMDH" not in rxn.id: dG_ref = 800.0
        if "PGI"  in rxn.id: dG_ref = 1500.0
        dH_rxn = next((DH_RXN_MAP[k] for k in DH_RXN_MAP if k in rxn.id), -10000.0)
        dG_T = max(100.0, dG_ref + dH_rxn * (1.0 - T_K / T_REF_K))
        eta  = max(0.01, math.tanh(dG_T / (2.0 * R_GAS * T_K)))
        coefs[rxn.forward_variable] = 1.0 / (km * alpha * eta)

    pool_con = model.problem.Constraint(0, lb=0, ub=pool_avail, name="pool")
    model.add_cons_vars(pool_con)
    model.solver.update()
    pool_con.set_linear_coefficients(coefs)
    return model, f_hsp

def run_fba_at(temp):
    """Run 3-objective FBA at temp, return raw FBA values (no alpha correction)."""
    model, _ = build_model(temp)
    model.objective = BIO_ID
    try:
        sol = pfba(model, fraction_of_optimum=1.0)
        if sol.status != "optimal": return None
        growth = float(sol.fluxes.get(BIO_ID, 0.0))
    except: return None
    m2 = model.copy()
    if BIO_ID in m2.reactions: m2.reactions.get_by_id(BIO_ID).lower_bound = 0.95*growth
    m2.objective = LYS_ID
    try:
        s2 = m2.optimize()
        lys = max(0.0, float(s2.fluxes.get(LYS_ID,0.0))) if s2.status=="optimal" else 0.0
    except: lys = 0.0
    m3 = model.copy()
    if BIO_ID in m3.reactions: m3.reactions.get_by_id(BIO_ID).lower_bound = 0.95*growth
    m3.objective = GLU_ID
    try:
        s3 = m3.optimize()
        glu = max(0.0, float(s3.fluxes.get(GLU_ID,0.0))) if s3.status=="optimal" else 0.0
    except: glu = 0.0
    return {"lys_fba": lys, "glu_fba": glu, "growth": growth}

# ── Baseline FBA sweep (1-degree steps for speed) ─────────────────────────────
TEMPS = [30+i*0.5 for i in range(31)]  # 30 to 45 in 0.5-step
print(f"\nRunning baseline FBA sweep ({len(TEMPS)} temperatures)...")
fba_raw = {}
for T in TEMPS:
    r = run_fba_at(T)
    if r:
        fba_raw[T] = r
        print(f"  {T:5.1f}C  Lys_FBA={r['lys_fba']:.4f}  Glu_FBA={r['glu_fba']:.4f}")
    else:
        fba_raw[T] = None
        print(f"  {T:5.1f}C  INFEASIBLE")

# Interpolate gaps
Tlist = sorted(fba_raw)
for i,T in enumerate(Tlist):
    if fba_raw[T] is None:
        pT=next((Tlist[j] for j in range(i-1,-1,-1) if fba_raw.get(Tlist[j])),None)
        nT=next((Tlist[j] for j in range(i+1,len(Tlist)) if fba_raw.get(Tlist[j])),None)
        if pT and nT:
            f=(T-pT)/(nT-pT)
            fba_raw[T]={k:fba_raw[pT][k]+f*(fba_raw[nT][k]-fba_raw[pT][k]) for k in fba_raw[pT]}
            print(f"  {T:5.1f}C  interpolated")

def compute_alpha_custom(p, T_K, hill_n=None):
    """Same formula as enzyme_thermal_params.compute_alpha."""
    E_a = p.get("E_a", 50000.0)
    H_d = p["H_d"]; S_d = p["S_d"]
    n   = hill_n if hill_n is not None else p.get("hill_n", 1.0)
    def fd(T):
        exp = max(-500, min(500, -n*(H_d - T*S_d)/(R_GAS*T)))
        return 1.0/(1.0+math.exp(exp))
    arr = math.exp(-E_a/R_GAS*(1.0/T_K - 1.0/T_REF_K))
    fd_ref = fd(T_REF_K)
    return 0.0 if fd_ref < 1e-12 else arr*fd(T_K)/fd_ref

def t50_from_alpha_curve(p, hill_n=None, T_scan=None):
    """
    FIX Issue 2: Compute T50 analytically from α(T)/α(30°C) curve.
    T50 = temperature where enzyme activity drops to 50% of 30°C value.
    This is a pure thermal-parameter property — no FBA fluxes involved.
    (FBA fluxes already have α embedded in LP coefficients; using them
    again would cause α² distortion.)
    """
    if T_scan is None:
        T_scan = [T_REF_K-273.15 + i*0.05 for i in range(800)]  # 30→70°C, 0.05°C steps
    a_ref = compute_alpha_custom(p, T_REF_K, hill_n)
    if a_ref < 1e-12: return None
    target = 0.5 * a_ref
    prev_a = a_ref
    for T in T_scan[1:]:
        a = compute_alpha_custom(p, T+273.15, hill_n)
        if prev_a >= target > a:
            frac = (prev_a - target) / max(prev_a - a, 1e-9)
            return (T - 0.05) + frac * 0.05
        prev_a = a
    return None

def metrics_from(p_lysc, p_gdh, n_lysc=None, n_gdh=None):
    """
    FIX Issue 2: Compute sensitivity metrics from α(T) curves directly.
    T50 = analytical from thermal model; Enzyme@40C = α(40°C)/α(30°C).
    No FBA required. No double-counting.
    """
    t50l = t50_from_alpha_curve(p_lysc, n_lysc)
    t50g = t50_from_alpha_curve(p_gdh,  n_gdh)
    T40K = 40.0 + 273.15
    lys_40 = compute_alpha_custom(p_lysc, T40K, n_lysc)   # relative to 30°C (α(30°C)=1.0)
    glu_40 = compute_alpha_custom(p_gdh,  T40K, n_gdh)
    
    # Calculate '_eff' trajectory for plotting
    eff = {}
    lys_ref = fba_raw[30.0]["lys_fba"]
    glu_ref = fba_raw[30.0]["glu_fba"]
    for T in TEMPS:
        T_K = T + 273.15
        al = compute_alpha_custom(p_lysc, T_K, n_lysc)
        ag = compute_alpha_custom(p_gdh,  T_K, n_gdh)
        eff[T] = {"lys_eff": lys_ref * al, "glu_eff": glu_ref * ag}
        
    return {
        "T50_Lys":  t50l,
        "T50_Glu":  t50g,
        "Gap":      (t50g - t50l) if (t50l and t50g) else None,
        "Lys@40C":  lys_40,
        "Glu@40C":  glu_40,
        "_eff":     eff,
    }



# Baseline
base = metrics_from(p_lysc_base, p_gdh_base, n_lysc_base, n_gdh_base)
print(f"\nBaseline: T50_Lys={base['T50_Lys']:.2f}C  T50_Glu={base['T50_Glu']:.2f}C  Gap={base['Gap']:.2f}C")

def shift_Tm(p_base, delta_C):
    p = copy.deepcopy(p_base)
    new_Tm_K = p["H_d"]/p["S_d"] + delta_C
    p["H_d"]  = new_Tm_K * p["S_d"]
    return p

params = ["Ea_LysC","Tm_LysC","n_LysC","Ea_GDH","Tm_GDH","n_GDH"]
all_m  = {"Baseline": base}

print(f"\n{'Case':30s}  {'DT50_Lys':>10}  {'DT50_Glu':>10}  {'DGap':>8}")
print("-"*65)
for sign, sfx in [(+1,"+20%"),(-1,"-20%")]:
    # LysC Ea
    p=copy.deepcopy(p_lysc_base); p["E_a"]*=(1+sign*PERTURB)
    m=metrics_from(p, p_gdh_base, n_lysc_base, n_gdh_base); name=f"Ea_LysC_{sfx}"
    all_m[name]=m; dl=(m['T50_Lys']-base['T50_Lys']) if m['T50_Lys'] else 0; dg=(m['T50_Glu']-base['T50_Glu']) if m['T50_Glu'] else 0; dG=(m['Gap']-base['Gap']) if m['Gap'] else 0
    print(f"  {name:28s}  {dl:+.2f}C  {dg:+.2f}C  {dG:+.2f}C")
    # LysC Tm
    p=shift_Tm(p_lysc_base, sign*PERTURB*Tm_lysc_C)
    m=metrics_from(p, p_gdh_base, n_lysc_base, n_gdh_base); name=f"Tm_LysC_{sfx}"
    all_m[name]=m; dl=(m['T50_Lys']-base['T50_Lys']) if m['T50_Lys'] else 0; dg=(m['T50_Glu']-base['T50_Glu']) if m['T50_Glu'] else 0; dG=(m['Gap']-base['Gap']) if m['Gap'] else 0
    print(f"  {name:28s}  {dl:+.2f}C  {dg:+.2f}C  {dG:+.2f}C")
    # LysC n
    n_new=n_lysc_base*(1+sign*PERTURB)
    m=metrics_from(p_lysc_base, p_gdh_base, n_new, n_gdh_base); name=f"n_LysC_{sfx}"
    all_m[name]=m; dl=(m['T50_Lys']-base['T50_Lys']) if m['T50_Lys'] else 0; dg=(m['T50_Glu']-base['T50_Glu']) if m['T50_Glu'] else 0; dG=(m['Gap']-base['Gap']) if m['Gap'] else 0
    print(f"  {name:28s}  {dl:+.2f}C  {dg:+.2f}C  {dG:+.2f}C")
    # GDH Ea
    p=copy.deepcopy(p_gdh_base); p["E_a"]*=(1+sign*PERTURB)
    m=metrics_from(p_lysc_base, p, n_lysc_base, n_gdh_base); name=f"Ea_GDH_{sfx}"
    all_m[name]=m; dl=(m['T50_Lys']-base['T50_Lys']) if m['T50_Lys'] else 0; dg=(m['T50_Glu']-base['T50_Glu']) if m['T50_Glu'] else 0; dG=(m['Gap']-base['Gap']) if m['Gap'] else 0
    print(f"  {name:28s}  {dl:+.2f}C  {dg:+.2f}C  {dG:+.2f}C")
    # GDH Tm
    p=shift_Tm(p_gdh_base, sign*PERTURB*Tm_gdh_C)
    m=metrics_from(p_lysc_base, p, n_lysc_base, n_gdh_base); name=f"Tm_GDH_{sfx}"
    all_m[name]=m; dl=(m['T50_Lys']-base['T50_Lys']) if m['T50_Lys'] else 0; dg=(m['T50_Glu']-base['T50_Glu']) if m['T50_Glu'] else 0; dG=(m['Gap']-base['Gap']) if m['Gap'] else 0
    print(f"  {name:28s}  {dl:+.2f}C  {dg:+.2f}C  {dG:+.2f}C")
    # GDH n
    n_new=n_gdh_base*(1+sign*PERTURB)
    m=metrics_from(p_lysc_base, p_gdh_base, n_lysc_base, n_new); name=f"n_GDH_{sfx}"
    all_m[name]=m; dl=(m['T50_Lys']-base['T50_Lys']) if m['T50_Lys'] else 0; dg=(m['T50_Glu']-base['T50_Glu']) if m['T50_Glu'] else 0; dG=(m['Gap']-base['Gap']) if m['Gap'] else 0
    print(f"  {name:28s}  {dl:+.2f}C  {dg:+.2f}C  {dG:+.2f}C")

# ── Figures ───────────────────────────────────────────────────────────────────
try:
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
        "font.size":10,"axes.labelsize":11,"axes.titlesize":11,"axes.titleweight":"bold",
        "xtick.labelsize":9,"ytick.labelsize":9,"legend.fontsize":8.5,"figure.dpi":300,
        "axes.linewidth":0.8,"axes.spines.top":False,"axes.spines.right":False})

    def styled(ax):
        ax.set_facecolor("white")
        ax.tick_params(direction="out",length=3,width=0.8,colors="black")
        for sp in ["left","bottom"]: ax.spines[sp].set_color("black"); ax.spines[sp].set_linewidth(0.8)
        return ax

    CLR_P="#d6604d"; CLR_N="#4393c3"
    metrics_keys = ["T50_Lys","T50_Glu","Gap","Lys@40C","Glu@40C"]
    metrics_disp = ["$T_{50}$ Lys (deg C)","$T_{50}$ Glu (deg C)","Gap (deg C)","Lys@40C","Glu@40C"]

    # ── Tornado ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1,3,figsize=(15,5)); fig.patch.set_facecolor("white")
    fig.subplots_adjust(wspace=0.42)
    for col,(mkey,mlabel) in enumerate(zip(metrics_keys[:3],metrics_disp[:3])):
        ax=styled(axes[col]); bv=base.get(mkey)
        if bv is None: ax.set_title(mlabel); continue
        bp=[]; bm=[]
        for param in params:
            vp=all_m.get(param+"_+20%",{}).get(mkey); vm=all_m.get(param+"_-20%",{}).get(mkey)
            bp.append((vp-bv) if vp is not None else 0.0)
            bm.append((vm-bv) if vm is not None else 0.0)
        order=sorted(range(len(params)),key=lambda i:-max(abs(bp[i]),abs(bm[i])))
        ps=[params[i] for i in order]; bp_=[bp[i] for i in order]; bm_=[bm[i] for i in order]
        y=np.arange(len(ps))
        ax.barh(y,    bp_, height=0.35, color=CLR_P, alpha=0.85, label="+20%")
        ax.barh(y-0.38,bm_,height=0.35, color=CLR_N, alpha=0.85, label="-20%")
        ax.axvline(0,color="black",lw=0.8,alpha=0.7)
        ax.set_yticks(y-0.19); ax.set_yticklabels([p.replace("_"," ") for p in ps],fontsize=9)
        ax.set_xlabel("Delta "+mlabel); ax.set_title(f"({chr(65+col)})  {mlabel}")
        if col==0: ax.legend(frameon=True,edgecolor="#cccccc",fancybox=False,fontsize=8.5)
        ax.text(0.01,-0.14,f"Baseline = {bv:.2f}",transform=ax.transAxes,fontsize=8,color="#555")
    fig.suptitle("Parameter sensitivity — Tornado chart (+/-20% perturbation)",fontsize=11,fontweight="bold",y=1.01)
    p1=os.path.join(OUT_DIR,"sensitivity_tornado.png")
    fig.savefig(p1,dpi=300,bbox_inches="tight",facecolor="white",edgecolor="none"); plt.close(fig)
    print(f"\nSaved: {p1}")

    # ── Perturbed curves ───────────────────────────────────────────────────────
    base_eff=base["_eff"]
    ref30=base_eff.get(30.0) or {}
    lr=max(ref30.get("lys_eff",1e-9),1e-9); gr=max(ref30.get("glu_eff",1e-9),1e-9)
    param_disp={"Ea_LysC":"LysC Ea","Tm_LysC":"LysC Tm","n_LysC":"LysC Hill n",
                "Ea_GDH":"GDH Ea","Tm_GDH":"GDH Tm","n_GDH":"GDH Hill n"}
    fig2,ax2s=plt.subplots(2,3,figsize=(15,9)); fig2.patch.set_facecolor("white")
    fig2.subplots_adjust(hspace=0.44,wspace=0.36)
    for idx,param in enumerate(params):
        row,ci=divmod(idx,3); ax=styled(ax2s[row,ci])
        bl=[base_eff.get(T,{}).get("lys_eff",0)/lr if base_eff.get(T) else 0 for T in TEMPS]
        bg=[base_eff.get(T,{}).get("glu_eff",0)/gr if base_eff.get(T) else 0 for T in TEMPS]
        ax.plot(TEMPS,bl,color=CLR_P,lw=1.8,label="Lys (baseline)")
        ax.plot(TEMPS,bg,color=CLR_N,lw=1.8,label="Glu (baseline)")
        for sign,ls,lbl in [(+1,"--","+20%"),(-1,":","-20%")]:
            ck=param+("_+20%" if sign>0 else "_-20%"); eff=all_m.get(ck,{}).get("_eff",{})
            lp=[eff.get(T,{}).get("lys_eff",0)/lr if eff.get(T) else 0 for T in TEMPS]
            gp=[eff.get(T,{}).get("glu_eff",0)/gr if eff.get(T) else 0 for T in TEMPS]
            ax.plot(TEMPS,lp,color=CLR_P,lw=1.1,ls=ls,alpha=0.8,label=f"Lys {lbl}")
            ax.plot(TEMPS,gp,color=CLR_N,lw=1.1,ls=ls,alpha=0.8,label=f"Glu {lbl}")
        ax.axhline(0.5,color="#888",lw=0.7,ls="--",alpha=0.7)
        ax.axvspan(37,41,alpha=0.07,color="#f4a582",zorder=0)
        ax.set_xlim(29.5,45.5); ax.set_ylim(-0.05,2.0)
        ax.set_xlabel("Temperature (deg C)"); ax.set_ylabel("Rel. production")
        ax.set_title(f"({chr(65+idx)})  {param_disp[param]}")
        ax.legend(frameon=True,edgecolor="#cccccc",fancybox=False,fontsize=7.5)
    fig2.suptitle("Production curves under +/-20% perturbations",fontsize=11,fontweight="bold",y=1.01)
    p2=os.path.join(OUT_DIR,"sensitivity_curves.png")
    fig2.savefig(p2,dpi=300,bbox_inches="tight",facecolor="white",edgecolor="none"); plt.close(fig2)
    print(f"Saved: {p2}")

    # ── Elasticity heatmap ─────────────────────────────────────────────────────
    fig3,ax3=plt.subplots(figsize=(9,5)); fig3.patch.set_facecolor("white"); styled(ax3)
    SC=np.zeros((len(params),len(metrics_keys)))
    for i,param in enumerate(params):
        for j,mkey in enumerate(metrics_keys):
            bv=base.get(mkey)
            if bv is None or abs(bv)<1e-9: continue
            vp=all_m.get(param+"_+20%",{}).get(mkey)
            vm=all_m.get(param+"_-20%",{}).get(mkey)
            if vp is None or vm is None: continue
            SC[i,j]=((vp-vm)/2.0/bv)/PERTURB
    im=ax3.imshow(SC,cmap="RdBu_r",vmin=-2,vmax=2,aspect="auto")
    ax3.set_xticks(range(len(metrics_keys))); ax3.set_xticklabels(metrics_disp,fontsize=9.5)
    ax3.set_yticks(range(len(params))); ax3.set_yticklabels([p.replace("_"," ") for p in params],fontsize=9.5)
    ax3.set_title("Sensitivity elasticity  (dOutput/Output_base) / (dParam/Param_base)",fontweight="bold")
    cb=plt.colorbar(im,ax=ax3,fraction=0.035,pad=0.03); cb.set_label("Elasticity",fontsize=9); cb.ax.tick_params(labelsize=8)
    for i in range(len(params)):
        for j in range(len(metrics_keys)):
            v=SC[i,j]
            ax3.text(j,i,f"{v:+.2f}",ha="center",va="center",fontsize=8.5,
                     color="white" if abs(v)>1.0 else "black",fontweight="bold" if abs(v)>1.0 else "normal")
    for sp in ax3.spines.values(): sp.set_visible(False)
    ax3.tick_params(length=0)
    fig3.suptitle("Parameter sensitivity elasticity — ecCGL1 heat stress model",fontsize=11,fontweight="bold",y=1.02)
    p3=os.path.join(OUT_DIR,"sensitivity_heatmap.png")
    fig3.savefig(p3,dpi=300,bbox_inches="tight",facecolor="white",edgecolor="none"); plt.close(fig3)
    print(f"Saved: {p3}")
    print("\nAll done.")

except Exception as e:
    import traceback; print("Figure error:",e); traceback.print_exc()
