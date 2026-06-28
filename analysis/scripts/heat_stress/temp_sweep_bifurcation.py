"""
Temperature Response Curves + Lys/Glu Bifurcation Analysis (FIXED)
====================================================================
Strategy per temperature:
  1. Maximize biomass -> growth*
  2. Fix growth >= 0.95 * growth*, maximize EX_lys_L_e -> max Lys production
  3. Fix growth >= 0.95 * growth*, maximize EX_glu_L_e -> max Glu production
  4. Extract LysC (ASPK) and GDH (GLUDy) internal fluxes from pFBA run
"""

import sys, math, json, os
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(ROOT_DIR, "backend"))

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

# All three timepoints are sampled at 40°C (heat stress time course at constant 40°C)
# 1h = acute response (1 hour after shift to 40°C)
# 4h = early adaptation (4 hours at 40°C)
# 24h = chronic response (24 hours at 40°C)
# For the temperature sweep, we use:
#   - At 40°C: proteomics-averaged calibration (best available data)
#   - All other temperatures: alpha(T) thermal factors only
PROT_PERTURBS = {
    "1h":  {"pgi": 0.96, "gdh": 0.98, "pyk": 1.42, "mdh": 1.01, "icdh": 1.10},
    "4h":  {"pgi": 0.86, "gdh": 1.14, "pyk": 1.09, "mdh": 1.14, "icdh": 0.77},
    "24h": {"pgi": 0.61, "gdh": 0.76, "pyk": 0.45, "mdh": 0.41, "icdh": 0.57},
}

# Average across all three 40°C timepoints for a representative 40°C calibration
PROT_AVG_40C = {
    k: (PROT_PERTURBS["1h"].get(k, 1.0) +
        PROT_PERTURBS["4h"].get(k, 1.0) +
        PROT_PERTURBS["24h"].get(k, 1.0)) / 3.0
    for k in PROT_PERTURBS["1h"]
}

def get_perturbs(temp):
    """Return proteomics perturbations: only available AT 40°C.
    All three experimental timepoints (1h/4h/24h) were sampled at 40°C.
    For other temperatures, rely on alpha(T) thermal correction only.
    """
    if 39.5 <= temp <= 40.5:
        # At 40°C: use averaged proteomics calibration
        return PROT_AVG_40C
    else:
        # No proteomics data at other temperatures — alpha(T) only
        return {}

BIO_ID  = None
LYS_ID  = "EX_lys_L_e"
GLU_ID  = "EX_glu_L_e"
ASPK_ID = "ASPK"
GDH_ID  = "GLUDy"

for bid in ["CG_biomass_cgl_ATCC13032", "BIOMASS_Cgl_ATCC13032", "Growth"]:
    if bid in BASE_MODEL.reactions: BIO_ID = bid; break

# ── Pre-compute alpha reference at 30°C ────────────────────────────────────
p_lysc = GENE_LOCUS_PARAMS.get("Cgl0251")
p_gdh  = GENE_LOCUS_PARAMS.get("Cgl2079")
alpha_lysc_ref = compute_alpha(p_lysc, T_REF_K) if p_lysc else 1.0
alpha_gdh_ref  = compute_alpha(p_gdh,  T_REF_K) if p_gdh  else 1.0

def build_base_model():
    """Build base model with media and protein pool constraint at 30°C for reference."""
    model = BASE_MODEL.copy()
    
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

    coefs = {}
    for rxn in model.reactions:
        km = KCAT_MAP.get(rxn.id)
        if not km: continue
        gene_loci = [g.id.replace("g_","").replace("gene_","") for g in rxn.genes]
        p = get_params(rxn.id, gene_loci)
        try:    alpha = compute_alpha(p, T_REF_K)
        except: alpha = 1.0
        dG_ref = 3000.0
        if "ICDHyr" in rxn.id:   dG_ref = 1000.0
        if "MDH" in rxn.id and "UAMDH" not in rxn.id: dG_ref = 800.0
        if "PGI"  in rxn.id:     dG_ref = 1500.0
        dG_T = max(100.0, dG_ref)
        eta  = max(0.01, math.tanh(dG_T / (2.0 * R_GAS * T_REF_K)))
        coefs[rxn.forward_variable] = 1.0 / (km * alpha * eta)

    pool_con = model.problem.Constraint(0, lb=0, ub=PROTEIN_POOL, name="pool")
    model.add_cons_vars(pool_con)
    model.solver.update()
    pool_con.set_linear_coefficients(coefs)

    model.objective = BIO_ID
    try:
        sol = pfba(model, fraction_of_optimum=1.0)
        aspk_ref = abs(float(sol.fluxes.get(ASPK_ID, 0.0)))
        gdh_ref  = abs(float(sol.fluxes.get(GDH_ID,  0.0)))
        return aspk_ref, gdh_ref
    except:
        return 0.1, 1.0

print("Computing 30C reference fluxes...", end=" ", flush=True)
ASPK_REF, GDH_REF = build_base_model()
# Use generous reference Vmax (scale up to prevent over-tightening at low T)
ASPK_VMAX = max(ASPK_REF * 1.1, 0.05)
GDH_VMAX  = max(GDH_REF  * 1.1, 0.5)
print(f"ASPK_ref={ASPK_REF:.4f}  GDH_ref={GDH_REF:.4f}")

def build_constrained_model(temp):
    """Build model for given temperature with full thermal constraints."""
    model  = BASE_MODEL.copy()
    T_K    = temp + 273.15
    t_diff = max(0.0, temp - 30.0)

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

    hsp_frac   = 0.20 * (t_diff**3) / (11.0**3 + t_diff**3)
    pool_avail = PROTEIN_POOL * (1.0 - hsp_frac)

    # FIX Issue 7: proteomics FC → constraint coefficient modification (ecFBA-correct)
    # Gene locus → enzyme name → FC mapping for constraint layer update
    PROT_LOCUS_MAP = {
        "Cgl0851": "pgi",   # PGI
        "Cgl2079": "gdh",   # GDH
        "Cgl2089": "pyk",   # PYK1
        "Cgl2910": "pyk",   # PYK2
        "Cgl2380": "mdh",   # MDH
        "Cgl0664": "icdh",  # ICDH
    }
    perturbs = get_perturbs(temp)

    coefs = {}
    for rxn in model.reactions:
        km = KCAT_MAP.get(rxn.id)
        if not km: continue
        gene_loci = [g.id.replace("g_","").replace("gene_","") for g in rxn.genes]
        p = get_params(rxn.id, gene_loci)
        try:    alpha = compute_alpha(p, T_K)
        except: alpha = 1.0
        dG_ref = 3000.0
        if "ICDHyr" in rxn.id:   dG_ref = 1000.0
        if "MDH" in rxn.id and "UAMDH" not in rxn.id: dG_ref = 800.0
        if "PGI"  in rxn.id:     dG_ref = 1500.0
        # FIX Issue 5 (η): use van't Hoff with per-reaction dH_rxn where known
        DH_RXN = {"ICDHyr": 8400.0, "MDH": 29000.0, "PGI": -2500.0}
        dH_rxn = next((DH_RXN[k] for k in DH_RXN if k in rxn.id), -10000.0)
        dG_T = max(100.0, dG_ref + dH_rxn * (1.0 - T_K / T_REF_K))
        eta  = max(0.01, math.tanh(dG_T / (2.0 * R_GAS * T_K)))
        base_coef = 1.0 / (km * alpha * eta)
        # FIX Issue 7: apply proteomics FC to constraint coefficient
        # More enzyme (FC > 1) → lower constraint cost → divide coef by FC
        prot_fc = 1.0
        for locus in gene_loci:
            enz_name = PROT_LOCUS_MAP.get(locus)
            if enz_name and enz_name in perturbs:
                prot_fc = perturbs[enz_name]
                break
        coefs[rxn.forward_variable] = base_coef / prot_fc

    pool_con = model.problem.Constraint(0, lb=0, ub=pool_avail, name="pool")
    model.add_cons_vars(pool_con)
    model.solver.update()
    pool_con.set_linear_coefficients(coefs)

    return model, hsp_frac

def run_temperature(temp):
    model, hsp_frac = build_constrained_model(temp)
    T_K = temp + 273.15

    # Step 1: Maximize growth
    model.objective = BIO_ID
    try:
        sol = pfba(model, fraction_of_optimum=1.0)
        if sol.status != "optimal": return None
        growth_opt = float(sol.fluxes.get(BIO_ID, 0.0))
        lysc_g     = abs(float(sol.fluxes.get(ASPK_ID, 0.0)))
        gdh_g      = abs(float(sol.fluxes.get(GDH_ID,  0.0)))
    except: return None

    # Step 2: Maximize Lys (with 95% growth floor)
    model2 = model.copy()
    if BIO_ID in model2.reactions:
        model2.reactions.get_by_id(BIO_ID).lower_bound = 0.95 * growth_opt
    model2.objective = LYS_ID
    try:
        sol2 = model2.optimize()
        lys_fba = max(0.0, float(sol2.fluxes.get(LYS_ID, 0.0))) if sol2.status=="optimal" else 0.0
    except: lys_fba = 0.0

    # Step 3: Maximize Glu (with 95% growth floor)
    model3 = model.copy()
    if BIO_ID in model3.reactions:
        model3.reactions.get_by_id(BIO_ID).lower_bound = 0.95 * growth_opt
    model3.objective = GLU_ID
    try:
        sol3 = model3.optimize()
        glu_fba = max(0.0, float(sol3.fluxes.get(GLU_ID, 0.0))) if sol3.status=="optimal" else 0.0
    except: glu_fba = 0.0

    # FIX Issue 1: α(T) is already embedded in LP constraint coefficients.
    # build_constrained_model() sets coef = 1/(km × α(T) × η) for every enzyme.
    # The LP-optimal fluxes lys_fba/glu_fba already reflect thermal enzyme capacity.
    # DO NOT multiply again — that was an α² distortion.
    # Keep α values for diagnostic reporting only.
    alpha_lysc_T = compute_alpha(p_lysc, T_K) if p_lysc else 1.0
    alpha_gdh_T  = compute_alpha(p_gdh,  T_K) if p_gdh  else 1.0

    # lys / glu are the constraint-layer FBA solutions (no post-processing multiply)
    lys_eff = lys_fba
    glu_eff = glu_fba

    return {
        "growth": max(0.0, growth_opt),
        "lys":    max(0.0, lys_eff),
        "glu":    max(0.0, glu_eff),
        "lys_fba": lys_fba,
        "glu_fba": glu_fba,
        "lysc":   lysc_g,
        "gdh":    gdh_g,
        "hsp":    hsp_frac,
        "pool":   PROTEIN_POOL * (1.0 - hsp_frac),
        "alpha_lysc": alpha_lysc_T,
        "alpha_gdh":  alpha_gdh_T,
    }

# ── Sweep ───────────────────────────────────────────────────────────────────
temps   = [30.0 + i*0.5 for i in range(31)]
results = {}

print(f"{'Temp':>6}  {'Growth':>8}  {'MaxLys':>8}  {'MaxGlu':>8}  {'LysC':>8}  {'GDH':>8}  {'HSP%':>7}")
print("-" * 68)
for T in temps:
    r = run_temperature(T)
    if r:
        results[T] = r
        print(f"  {T:>4.1f}C  {r['growth']:>8.4f}  {r['lys']:>8.4f}  {r['glu']:>8.4f}  "
              f"{r['lysc']:>8.4f}  {r['gdh']:>8.4f}  {r['hsp']*100:>6.1f}%")
    else:
        results[T] = None
        print(f"  {T:>4.1f}C  INFEASIBLE")

# FIX Issue 6: INFEASIBLE points get zero production (not interpolated).
# Linear interpolation across the T50 transition zone would fabricate smooth
# data precisely where the bifurcation signal is strongest, distorting T50.
infeasible_filled = []
for i, T in enumerate(temps):
    if results.get(T) is None:
        prev_T = next((temps[j] for j in range(i-1, -1, -1) if results.get(temps[j])), None)
        if prev_T is not None:
            # Carry growth forward; set production to 0 (LP infeasible = constraint violated)
            prev = results[prev_T]
            results[T] = {k: 0.0 if k in ("lys","glu","lys_fba","glu_fba") else prev.get(k, 0.0)
                          for k in prev}
        else:
            results[T] = None
        infeasible_filled.append(T)
        print(f"  {T:>4.1f}C  INFEASIBLE → production set to 0.0")

# ── Alpha curves ────────────────────────────────────────────────────────────
lysc_p = GENE_LOCUS_PARAMS.get("Cgl0251")
gdh_p  = GENE_LOCUS_PARAMS.get("Cgl2079")
valid_T = [T for T in temps if results.get(T)]

# ── Bifurcation ─────────────────────────────────────────────────────────────
print("\n--- Bifurcation Analysis ---")
r0 = results.get(30.0, {}) or {}
lys0 = r0.get("lys", 0.0); glu0 = r0.get("glu", 0.0)

lys_50 = glu_50 = None
prev_lr = prev_gr = 1.0
for T in valid_T[1:]:
    r = results[T]
    lr = r["lys"]/lys0 if lys0 > 1e-6 else 0.0
    gr = r["glu"]/glu0 if glu0 > 1e-6 else 0.0
    if lys_50 is None and lr < 0.5 and prev_lr >= 0.5: lys_50 = T
    if glu_50 is None and gr < 0.5 and prev_gr >= 0.5: glu_50 = T
    prev_lr, prev_gr = lr, gr

print(f"Lys 50% threshold temperature: {lys_50}°C")
print(f"Glu 50% threshold temperature: {glu_50}°C")
if lys_50 and glu_50:
    gap = glu_50 - lys_50
    print(f"Bifurcation gap: {gap:.1f}°C  (experimental: ~5.5°C)")

# ── Figures ─────────────────────────────────────────────────────────────────
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import numpy as np

    growth_v = [results[T]["growth"] for T in valid_T]
    lys_v    = [results[T]["lys"]    for T in valid_T]
    glu_v    = [results[T]["glu"]    for T in valid_T]
    lysc_v   = [results[T]["lysc"]   for T in valid_T]
    gdh_v    = [results[T]["gdh"]    for T in valid_T]
    hsp_v    = [results[T]["hsp"]*100 for T in valid_T]
    a_lysc   = [compute_alpha(lysc_p, T+273.15) for T in valid_T] if lysc_p else []
    a_gdh    = [compute_alpha(gdh_p,  T+273.15) for T in valid_T] if gdh_p  else []

    lys0n  = lys_v[0]  if lys_v[0]  > 1e-6 else 1.0
    glu0n  = glu_v[0]  if glu_v[0]  > 1e-6 else 1.0
    lysc0n = lysc_v[0] if lysc_v[0] > 1e-6 else 1.0
    gdh0n  = gdh_v[0]  if gdh_v[0]  > 1e-6 else 1.0

    lys_n  = [v/lys0n  for v in lys_v]
    glu_n  = [v/glu0n  for v in glu_v]
    lysc_n = [v/lysc0n for v in lysc_v]
    gdh_n  = [v/gdh0n  for v in gdh_v]

    # ── Academic color palette (colorblind-friendly, Nature/Science style) ────
    CLR = {
        "growth": "#2166ac",    # steel blue
        "lys":    "#d6604d",    # muted red
        "glu":    "#4393c3",    # medium blue
        "hsp":    "#f4a582",    # salmon
        "lysc":   "#d6604d",    # same as lys (LysC drives Lys)
        "gdh":    "#4393c3",    # same as glu (GDH drives Glu)
        "ref":    "#888888",    # medium gray
        "shade":  "#f0f0f0",    # light gray shading
    }

    # ── Font settings ─────────────────────────────────────────────────────────
    plt.rcParams.update({
        "font.family":       "sans-serif",
        "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size":         10,
        "axes.labelsize":    11,
        "axes.titlesize":    11,
        "axes.titleweight":  "bold",
        "xtick.labelsize":   9,
        "ytick.labelsize":   9,
        "legend.fontsize":   8.5,
        "figure.dpi":        300,
        "axes.linewidth":    0.8,
        "lines.linewidth":   1.5,
        "axes.spines.top":   False,
        "axes.spines.right": False,
    })

    def styled(ax):
        """Academic white style: clean spines, no background color."""
        ax.set_facecolor("white")
        ax.tick_params(direction="out", length=3, width=0.8, colors="black")
        for spine in ["left", "bottom"]:
            ax.spines[spine].set_color("black")
            ax.spines[spine].set_linewidth(0.8)
        ax.xaxis.label.set_color("black")
        ax.yaxis.label.set_color("black")
        ax.title.set_color("black")
        return ax

    # ─────────────────────── Figure 1: Temperature Response ──────────────────
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(hspace=0.38, wspace=0.35)

    # ── A: Growth rate ────────────────────────────────────────────────────────
    ax = styled(axes[0, 0])
    ax.plot(valid_T, growth_v, color=CLR["growth"], lw=1.8,
            marker="o", ms=3.5, mfc="white", mew=1.2, clip_on=False)
    if infeasible_filled:
        interp_g = [results[T]["growth"] for T in infeasible_filled if T in results]
        ax.scatter(infeasible_filled, interp_g, s=50, facecolors="none",
                   edgecolors=CLR["growth"], lw=1.2, zorder=5, label="Interpolated")
        ax.legend(frameon=True, edgecolor="#cccccc", fancybox=False)
    ax.axvspan(37, 41, alpha=0.08, color="#f4a582", zorder=0)
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Growth rate (gDW gDW$^{-1}$ h$^{-1}$)")
    ax.set_title("(A)  Growth rate vs. temperature")
    ax.set_xlim(29.5, 45.5)

    # ── B: Lys + Glu normalized production ───────────────────────────────────
    ax = styled(axes[0, 1])
    ax.plot(valid_T, lys_n, color=CLR["lys"], lw=1.8,
            marker="o", ms=3.5, mfc="white", mew=1.2,
            label=f"Lys × $\\alpha_{{LysC}}(T)$  [$T_m$=37.0°C, $n$=4]")
    ax.plot(valid_T, glu_n, color=CLR["glu"], lw=1.8,
            marker="s", ms=3.5, mfc="white", mew=1.2,
            label=f"Glu × $\\alpha_{{GDH}}(T)$   [$T_m$=41.0°C, $n$=2.5]")
    ax.axhline(0.5, color="#888888", lw=0.8, ls="--", alpha=0.9,
               label="50% threshold")
    if lys_50:
        ax.axvline(lys_50, color=CLR["lys"], lw=1.0, ls=":", alpha=0.9,
                   label=f"Lys $T_{{50}}$ = {lys_50}°C")
    if glu_50:
        ax.axvline(glu_50, color=CLR["glu"], lw=1.0, ls=":", alpha=0.9,
                   label=f"Glu $T_{{50}}$ = {glu_50}°C")
    if infeasible_filled:
        interp_l  = [results[T]["lys"]/lys0n for T in infeasible_filled if T in results]
        interp_g2 = [results[T]["glu"]/glu0n for T in infeasible_filled if T in results]
        ax.scatter(infeasible_filled, interp_l,  s=60, facecolors="none",
                   edgecolors=CLR["lys"], lw=1.2, zorder=5)
        ax.scatter(infeasible_filled, interp_g2, s=60, facecolors="none",
                   edgecolors=CLR["glu"], lw=1.2, zorder=5)
    ax.axvspan(37, 41, alpha=0.08, color="#f4a582", zorder=0)
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Relative production capacity (vs. 30°C)")
    ax.set_title("(B)  Lys/Glu production capacity")
    ax.set_xlim(29.5, 45.5)
    ax.legend(frameon=True, edgecolor="#cccccc", fancybox=False,
              handlelength=1.8, loc="upper right")

    # ── C: LysC vs GDH alpha(T) ───────────────────────────────────────────────
    ax = styled(axes[1, 0])
    if a_lysc:
        ax.plot(valid_T, a_lysc, color=CLR["lysc"], lw=1.8,
                label="$\\alpha_{LysC}(T)$  $T_m$=37.0°C, $n$=4  (α₂β₂ tetramer)")
    if a_gdh:
        ax.plot(valid_T, a_gdh, color=CLR["gdh"], lw=1.8, ls="--",
                label="$\\alpha_{GDH}(T)$   $T_m$=41.0°C, $n$=2.5 (hexamer)")
    ax.axvspan(35.5, 38.5, alpha=0.07, color=CLR["lysc"], zorder=0)
    ax.axvspan(39.5, 42.5, alpha=0.07, color=CLR["gdh"],  zorder=0)
    ax.axvline(37.0, color="#888888", lw=0.8, ls="--", alpha=0.7)
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Relative enzyme activity $\\alpha(T)$")
    ax.set_title("(C)  LysC vs. GDH thermal activity")
    ax.set_xlim(29.5, 45.5)
    ax.legend(frameon=True, edgecolor="#cccccc", fancybox=False, handlelength=2.0)
    if lysc_p and gdh_p:
        tm_l = lysc_p["H_d"] / lysc_p["S_d"] - 273.15
        tm_g = gdh_p["H_d"]  / gdh_p["S_d"]  - 273.15
        ax.annotate("", xy=(tm_l, 0.50), xytext=(tm_g, 0.50),
                    arrowprops=dict(arrowstyle="<->", color="black", lw=1.2))
        ax.text((tm_l+tm_g)/2, 0.54, f"$\\Delta T_m$ = {tm_g-tm_l:.1f}°C",
                color="black", ha="center", fontsize=9)

    # ── D: HSP cost + internal enzyme fluxes ──────────────────────────────────
    ax = styled(axes[1, 1])
    ax.fill_between(valid_T, hsp_v, color=CLR["hsp"], alpha=0.3, zorder=0)
    ax.plot(valid_T, hsp_v, color=CLR["hsp"], lw=1.8,
            label="HSP proteome fraction (%)")
    ax2 = ax.twinx()
    ax2.set_facecolor("white")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_color("black")
    ax2.spines["right"].set_linewidth(0.8)
    ax2.plot(valid_T, lysc_n, color=CLR["lysc"], lw=1.5, ls="--",
             label="LysC flux (norm.)")
    ax2.plot(valid_T, gdh_n,  color=CLR["gdh"],  lw=1.5, ls="-.",
             label="GDH flux (norm.)")
    ax2.tick_params(direction="out", length=3, width=0.8)
    ax2.set_ylabel("Normalized enzyme flux", color="black")
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("HSP proteome fraction (%)")
    ax.set_title("(D)  HSP cost and enzyme fluxes")
    ax.set_xlim(29.5, 45.5)
    lines1, lab1 = ax.get_legend_handles_labels()
    lines2, lab2 = ax2.get_legend_handles_labels()
    ax.legend(lines1+lines2, lab1+lab2, frameon=True, edgecolor="#cccccc",
              fancybox=False, fontsize=8)

    fig.suptitle("ecCGL1 temperature response curve — "
                 "$\\it{C. glutamicum}$ (30–45°C)",
                 fontsize=12, fontweight="bold", y=1.01, color="black")
    out1 = os.path.join(OUT_DIR, "temperature_response.png")
    fig.savefig(out1, dpi=300, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print("Saved: " + out1)

    # ─────────────────────── Figure 2: Bifurcation Analysis ──────────────────
    fig2 = plt.figure(figsize=(13, 8))
    fig2.patch.set_facecolor("white")
    gs = gridspec.GridSpec(2, 3, figure=fig2, hspace=0.46, wspace=0.40)

    # 2A: Lys/Glu ratio over temperature ─────────────────────────────────────
    ax2a = styled(fig2.add_subplot(gs[0, :2]))
    if lys0n > 1e-6 and glu0n > 1e-6:
        ratio = [l/g if g > 1e-6 else 0.0 for l, g in zip(lys_n, glu_n)]
        ax2a.plot(valid_T, ratio, color="#b35806", lw=1.8,
                  marker="D", ms=3.5, mfc="white", mew=1.2,
                  label="Lys / Glu production ratio (norm.)")
        ax2a.fill_between(valid_T, ratio, alpha=0.12, color="#b35806")
    ax2a.axhline(1.0, color="#888888", lw=0.8, ls="--", alpha=0.8,
                 label="Equal production")
    if lys_50:
        ax2a.axvline(lys_50, color=CLR["lys"], lw=1.0, ls=":",
                     label=f"Lys $T_{{50}}$ = {lys_50}°C")
    if glu_50:
        ax2a.axvline(glu_50, color=CLR["glu"], lw=1.0, ls=":",
                     label=f"Glu $T_{{50}}$ = {glu_50}°C")
    if lys_50 and glu_50:
        y_ann = 0.55
        ax2a.annotate("", xy=(lys_50, y_ann), xytext=(glu_50, y_ann),
                      arrowprops=dict(arrowstyle="<->", color="black", lw=1.2))
        ax2a.text((lys_50+glu_50)/2, y_ann+0.05,
                  f"Gap = {glu_50-lys_50:.1f}°C (exp. ~5.5°C)",
                  color="black", ha="center", fontsize=9)
    ax2a.set_xlabel("Temperature (°C)")
    ax2a.set_ylabel("Lys/Glu ratio (norm.)")
    ax2a.set_title("(A)  Lys/Glu production bifurcation")
    ax2a.set_xlim(29.5, 45.5)
    ax2a.legend(frameon=True, edgecolor="#cccccc", fancybox=False, fontsize=8.5)

    # 2B: Enzyme Tm ranking ────────────────────────────────────────────────────
    ax2b = styled(fig2.add_subplot(gs[0, 2]))
    tm_data = [
        ("LysC\n(Cgl0251)", 37.0, "#d6604d"),
        ("PFK\n(Cgl1250)",  36.9, "#f4a582"),
        ("AKGDH\n(Cgl1129)",38.1, "#fddbc7"),
        ("PDH E1\n(Cgl0766)",38.5,"#d1e5f0"),
        ("GDH\n(Cgl2079)",  41.0, "#4393c3"),
        ("MDH\n(Cgl2380)",  43.8, "#2166ac"),
        ("CS\n(Cgl0696)",   48.3, "#053061"),
        ("ICDH\n(Cgl0664)", 50.0, "#313695"),
    ]
    names = [x[0] for x in tm_data]
    tms   = [x[1] for x in tm_data]
    clrs  = [x[2] for x in tm_data]
    bars  = ax2b.barh(range(len(tms)), tms, color=clrs, height=0.6,
                      edgecolor="white", linewidth=0.3)
    ax2b.set_yticks(range(len(names)))
    ax2b.set_yticklabels(names, fontsize=7.5)
    ax2b.axvline(37, color="#d6604d", lw=1.2, ls="--", alpha=0.9, label="37°C")
    ax2b.axvline(40, color="#4393c3", lw=1.2, ls="--", alpha=0.9, label="40°C")
    ax2b.set_xlabel("$T_m$ (°C)")
    ax2b.set_title("(B)  Enzyme $T_m$ ranking")
    ax2b.legend(frameon=True, edgecolor="#cccccc", fancybox=False, fontsize=7.5)
    ax2b.annotate(f"$\\Delta T_m$ = {41.0-37.0:.1f}°C",
                  xy=(37.0, 4), xytext=(44, 3.4),
                  color="black", fontsize=8,
                  arrowprops=dict(arrowstyle="->", color="black", lw=1.0))
    ax2b.set_xlim(0, 55)

    # 2C: Alpha differential ───────────────────────────────────────────────────
    ax2c = styled(fig2.add_subplot(gs[1, :2]))
    if a_lysc and a_gdh:
        diff = [l - g for l, g in zip(a_lysc, a_gdh)]
        pos_c = "#d6604d"; neg_c = "#4393c3"
        ax2c.bar(valid_T, diff,
                 color=[pos_c if v > 0 else neg_c for v in diff],
                 width=0.42, alpha=0.75, label="$\\alpha_{LysC} - \\alpha_{GDH}$")
        ax2c.plot(valid_T, a_lysc, color=CLR["lysc"], lw=1.4, ls="--",
                  alpha=0.85, label="$\\alpha_{LysC}$")
        ax2c.plot(valid_T, a_gdh,  color=CLR["gdh"],  lw=1.4, ls="-.",
                  alpha=0.85, label="$\\alpha_{GDH}$")
        ax2c.axhline(0, color="black", lw=0.8, alpha=0.6)
        ax2c.set_xlabel("Temperature (°C)")
        ax2c.set_ylabel("$\\Delta\\alpha = \\alpha_{LysC} - \\alpha_{GDH}$")
        ax2c.set_title("(C)  LysC vs. GDH activity differential")
        ax2c.set_xlim(29.5, 45.5)
        ax2c.legend(frameon=True, edgecolor="#cccccc", fancybox=False, fontsize=8.5)

    # 2D: Flux heatmap ─────────────────────────────────────────────────────────
    ax2d = styled(fig2.add_subplot(gs[1, 2]))
    key_temps_hm = [30, 33, 35, 37, 38, 39, 40, 41, 42, 43, 44, 45]
    enz_labels = ["Lys\n(norm.)", "Glu\n(norm.)",
                  "LysC\n(norm.)", "GDH\n(norm.)", "Growth\n(norm.)"]
    g0 = growth_v[0] if growth_v[0] > 1e-6 else 1.0
    mat_rows = []
    for T in key_temps_hm:
        r = results.get(float(T))
        if not r:
            r = {"lys": 0, "glu": 0, "lysc": 0, "gdh": 0, "growth": 0}
        mat_rows.append([
            r["lys"]    / lys0n   if lys0n   > 1e-6 else 0,
            r["glu"]    / glu0n   if glu0n   > 1e-6 else 0,
            r["lysc"]   / lysc0n  if lysc0n  > 1e-6 else 0,
            r["gdh"]    / gdh0n   if gdh0n   > 1e-6 else 0,
            r["growth"] / g0,
        ])
    mat = np.array(mat_rows).T
    im = ax2d.imshow(mat, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1.5,
                     interpolation="nearest")
    ax2d.set_xticks(range(len(key_temps_hm)))
    ax2d.set_xticklabels([str(t) for t in key_temps_hm],
                          fontsize=6.5, rotation=45, ha="right")
    ax2d.set_yticks(range(len(enz_labels)))
    ax2d.set_yticklabels(enz_labels, fontsize=8)
    ax2d.set_xlabel("Temperature (°C)", fontsize=9)
    ax2d.set_title("(D)  Flux heatmap\n(norm. to 30°C)")
    cb = plt.colorbar(im, ax=ax2d, fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=7)
    # Remove imshow spines except bottom/left
    for sp in ax2d.spines.values(): sp.set_visible(False)

    fig2.suptitle("Lys/Glu bifurcation analysis — "
                  "$\\it{C. glutamicum}$  "
                  "LysC ($T_m$=37.0°C) vs. GDH ($T_m$=41.0°C)",
                  fontsize=11, fontweight="bold", y=1.01, color="black")

    out2 = os.path.join(OUT_DIR, "bifurcation_analysis.png")
    fig2.savefig(out2, dpi=300, bbox_inches="tight",
                 facecolor="white", edgecolor="none")
    plt.close(fig2)
    print("Saved: " + out2)
    print("Figures done.")

except Exception as e:
    import traceback
    print("Figure error: " + str(e))
    traceback.print_exc()

