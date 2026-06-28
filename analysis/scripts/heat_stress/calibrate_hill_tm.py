"""
Hill 系数 + 体内有效 Tm 校准
==============================
目标：
  Lys 产量 T50 = 37-38°C  (实验: Takeno 2010, Xu 2022)
  Glu 产量 T50 = 42-43°C  (实验: Takeno 2010)
  分叉温差   = ~5.5°C

方法：
  1. 在 compute_active_fraction 中引入 Hill 协同系数 n:
     K_eq = exp(-n * (H_d - T*S_d) / R*T)   →  更陡的 sigmoid
  2. 对 LysC 的体内有效 Tm (Tm_eff_C) 做网格搜索
  3. 对 n_lysc (LysC 四聚体协同性) 和 n_gdh (GDH 六聚体协同性) 做网格搜索
  4. 找到 RMSE(T50_lys - 37.5, T50_glu - 42.5) 最小的参数组合

物理依据:
  LysC: 四聚体 (α2β2), 别构调节 → n ≈ 2-4
  GDH:  六聚体, 协同解折叠   → n ≈ 2-3
  体内 Tm vs 体外 T50: 差值通常 3-6°C (底物竞争, 分子拥挤效应)
"""
import sys, os, math, json
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(ROOT_DIR, "backend"))
import numpy as np
from enzyme_thermal_params import GENE_LOCUS_PARAMS, compute_alpha

R = 8.314
T_REF_K = 303.15

def compute_active_hill(H_d, S_d, T_K, n=1.0):
    """Two-state unfolding with Hill cooperativity coefficient n."""
    try:
        exponent = -n * (H_d - T_K * S_d) / (R * T_K)
        exponent = max(-500.0, min(500.0, exponent))
        K_eq = math.exp(exponent)
        return 1.0 / (1.0 + K_eq)
    except:
        return 1e-6

def compute_alpha_hill(E_a, H_d, S_d, T_K, n=1.0):
    """Arrhenius × two-state-Hill combined correction factor."""
    arr = math.exp(-(E_a / R) * (1.0/T_K - 1.0/T_REF_K))
    f_ref = compute_active_hill(H_d, S_d, T_REF_K, n)
    f_t   = compute_active_hill(H_d, S_d, T_K,    n)
    alpha = arr * (f_t / f_ref) if f_ref > 1e-10 else 1e-6
    return max(1e-6, alpha)

def Tm_to_Hd_Sd(Tm_C, S_d_ref=None):
    """
    Convert a target Tm (°C) back to (H_d, S_d) keeping S_d fixed
    from the reference (body of literature); adjust H_d = Tm_K * S_d.
    """
    Tm_K = Tm_C + 273.15
    # Use reference S_d values (these encode the width of the transition)
    p_lysc = GENE_LOCUS_PARAMS["Cgl0251"]
    p_gdh  = GENE_LOCUS_PARAMS["Cgl2079"]
    return Tm_K * p_lysc["S_d"], Tm_K * p_gdh["S_d"]

# ── Current reference values ─────────────────────────────────────────────────
p_lysc = GENE_LOCUS_PARAMS["Cgl0251"]
p_gdh  = GENE_LOCUS_PARAMS["Cgl2079"]

print(f"Reference LysC: Tm={p_lysc['H_d']/p_lysc['S_d']-273.15:.1f}°C  Ea={p_lysc['E_a']/1000:.0f}kJ/mol")
print(f"Reference GDH:  Tm={p_gdh['H_d']/p_gdh['S_d']-273.15:.1f}°C  Ea={p_gdh['E_a']/1000:.0f}kJ/mol")
print()

# ── Simplified production model (without full FBA) ───────────────────────────
# Use alpha(T) as a proxy for relative production capacity
# T50 = temperature at which alpha falls to 50% of max alpha in 30-35°C range
def find_t50_alpha(Ea, Hd, Sd, n):
    """Find temperature where alpha(T) drops to 50% of peak in 30-35C range."""
    temps = np.arange(30, 48, 0.1)
    alphas = [compute_alpha_hill(Ea, Hd, Sd, T+273.15, n) for T in temps]
    peak_alpha = max(alphas[:51])   # max in 30-35C window
    thresh = 0.5 * peak_alpha
    for i, (T, a) in enumerate(zip(temps, alphas)):
        if i > 0 and a < thresh:
            return T
    return None

# ── Grid search ──────────────────────────────────────────────────────────────
# Targets: Lys T50 ≈ 37.5°C, Glu T50 ≈ 42.5°C, gap ≈ 5.0°C

TARGET_LYS_T50 = 37.5
TARGET_GLU_T50 = 42.5

Tm_lysc_range = np.arange(31.0, 38.0, 0.5)    # in-vivo Tm scan
n_lysc_range  = np.arange(1.0, 5.5, 0.5)       # Hill n for LysC
n_gdh_range   = np.arange(1.0, 4.5, 0.5)       # Hill n for GDH
Tm_gdh_range  = np.arange(39.0, 44.0, 0.5)     # GDH Tm scan (keep near 42)

best_rmse = 1e9
best_params = None
results = []

print("Grid searching (Tm_lysc, n_lysc, Tm_gdh, n_gdh)...")
total = len(Tm_lysc_range)*len(n_lysc_range)*len(Tm_gdh_range)*len(n_gdh_range)
count = 0
for Tm_l in Tm_lysc_range:
    for n_l in n_lysc_range:
        for Tm_g in Tm_gdh_range:
            for n_g in n_gdh_range:
                count += 1
                # Build params
                Hd_l = (Tm_l + 273.15) * p_lysc["S_d"]
                Hd_g = (Tm_g + 273.15) * p_gdh["S_d"]

                t50_l = find_t50_alpha(p_lysc["E_a"], Hd_l, p_lysc["S_d"], n_l)
                t50_g = find_t50_alpha(p_gdh["E_a"],  Hd_g, p_gdh["S_d"],  n_g)

                if t50_l is None or t50_g is None:
                    continue

                gap = t50_g - t50_l
                rmse = math.sqrt((t50_l - TARGET_LYS_T50)**2 + (t50_g - TARGET_GLU_T50)**2)

                results.append({
                    "Tm_lysc": Tm_l, "n_lysc": n_l, "Tm_gdh": Tm_g, "n_gdh": n_g,
                    "t50_lys": t50_l, "t50_glu": t50_g, "gap": gap, "rmse": rmse
                })

                if rmse < best_rmse:
                    best_rmse = rmse
                    best_params = {
                        "Tm_lysc": Tm_l, "n_lysc": n_l, "Tm_gdh": Tm_g, "n_gdh": n_g,
                        "t50_lys": t50_l, "t50_glu": t50_g, "gap": gap, "rmse": rmse
                    }

print(f"Searched {count} combinations, RMSE best = {best_rmse:.3f}°C\n")

# Sort and show top 10
results.sort(key=lambda x: x["rmse"])
print(f"{'Tm_lysc':>8} {'n_lysc':>7} {'Tm_gdh':>8} {'n_gdh':>7} {'T50_Lys':>9} {'T50_Glu':>9} {'Gap':>6} {'RMSE':>7}")
print("-"*72)
for r in results[:10]:
    print(f"  {r['Tm_lysc']:>6.1f}  {r['n_lysc']:>6.1f}  {r['Tm_gdh']:>7.1f}  {r['n_gdh']:>6.1f}  "
          f"{r['t50_lys']:>8.1f}  {r['t50_glu']:>8.1f}  {r['gap']:>5.1f}  {r['rmse']:>6.3f}")

print(f"\n{'='*72}")
print(f"BEST PARAMETERS:")
bp = best_params
print(f"  LysC: Tm_eff = {bp['Tm_lysc']:.1f}°C  (original: 36.9°C, reduced by "
      f"{36.9-bp['Tm_lysc']:.1f}°C for in-vivo effect)")
print(f"  LysC: Hill n = {bp['n_lysc']:.1f}  (n=1: no cooperativity; tetrameric ~2-4)")
print(f"  GDH:  Tm_eff = {bp['Tm_gdh']:.1f}°C  (original: 41.9°C)")
print(f"  GDH:  Hill n = {bp['n_gdh']:.1f}  (hexameric enzyme ~2-3)")
print(f"  Predicted T50_Lys = {bp['t50_lys']:.1f}°C  (target: {TARGET_LYS_T50}°C)")
print(f"  Predicted T50_Glu = {bp['t50_glu']:.1f}°C  (target: {TARGET_GLU_T50}°C)")
print(f"  Bifurcation gap   = {bp['gap']:.1f}°C  (target: ~5.5°C)")
print(f"  RMSE              = {bp['rmse']:.3f}°C")

# ── Update enzyme_thermal_params with best values ────────────────────────────
print(f"\n{'='*72}")
print("Computing new H_d values for GENE_LOCUS_PARAMS:")
bp = best_params
Hd_lysc_new = (bp['Tm_lysc'] + 273.15) * p_lysc["S_d"]
Hd_gdh_new  = (bp['Tm_gdh']  + 273.15) * p_gdh["S_d"]
print(f"  LysC Cgl0251: H_d = {Hd_lysc_new:.0f}  S_d = {p_lysc['S_d']}  n = {bp['n_lysc']}")
print(f"  GDH  Cgl2079: H_d = {Hd_gdh_new:.0f}   S_d = {p_gdh['S_d']}  n = {bp['n_gdh']}")
print(f"  Verify LysC Tm = {Hd_lysc_new/p_lysc['S_d']-273.15:.2f}°C")
print(f"  Verify GDH  Tm = {Hd_gdh_new /p_gdh['S_d'] -273.15:.2f}°C")
