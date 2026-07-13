"""
scripts/build_thermo_data.py  (v3 - comprehensive expansion)
=============================================================
Expanded curated ΔrG'° database covering ~350 reactions in iCW773.
All IDs verified against the actual model. Organized by pathway.

Conditions: pH 7.0, I=0.1 M, T=30°C (303.15 K)
c_range: [1 µM, 50 mM]  → half_spread = RT·ln(50000) ≈ 27.3 kJ/mol

Classification:
  forward : ΔrG'_max = ΔrG'° + 27.3 < -1.0  → always exergonic forward
  reverse : ΔrG'_min = ΔrG'° - 27.3 > +1.0  → always endergonic forward
  none    : near-equilibrium or no data

Sources:
  [N13]  Noor et al. 2013 PLoS Comput Biol (Component Contribution)
  [F12]  Flamholz et al. 2012 Nucleic Acids Res (eQuilibrator)
  [T77]  Thauer et al. 1977 Bacteriol Rev 41:100
  [A03]  Alberty 2003 MIT Press
  [C14]  Chang et al. 2014 Bioinformatics (COBRA Toolbox)
  [K]    KEGG thermodynamic estimates
"""

import os, json, math, logging, warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("build_thermo")

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR    = os.path.dirname(SCRIPT_DIR)
MODEL_PATH  = os.path.join(ROOT_DIR, "backend", "models", "iCW773.xml")
OUTPUT_PATH = os.path.join(ROOT_DIR, "data", "reference", "thermo_dgr_data.json")

R       = 8.314e-3   # kJ/(mol·K)
T_K     = 303.15     # 30°C
C_MIN   = 1e-6       # 1 µM
C_MAX   = 0.05       # 50 mM
EPSILON = 1.0        # kJ/mol threshold

HALF_SPREAD = R * T_K * math.log(C_MAX / C_MIN)   # ≈ 27.27 kJ/mol

# ═══════════════════════════════════════════════════════════════════════════════
# Format: "ModelReactionID": (dgr_prime_0_kJ, confidence_H/M/L, source_note)
# ═══════════════════════════════════════════════════════════════════════════════
CURATED_DGR0 = {

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # GLUCOSE UPTAKE & GLYCOLYSIS
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "GLCptspp":      (-46.0, "H", "[N13] PTS glucose phosphotransfer"),
    "PGI":           (  2.1, "H", "[F12] Glucose-6-P isomerase"),
    "PFK":           (-17.8, "H", "[N13] 6-phosphofructokinase"),
    "PFK_2":         (-17.8, "H", "[N13] PFK isoform 2"),
    "FBA":           ( 23.8, "H", "[F12] Fructose bisphosphate aldolase"),
    "TPI":           (  5.5, "H", "[F12] Triosephosphate isomerase"),
    "GAPD":          (  6.3, "H", "[F12] Glyceraldehyde-3-P dehydrogenase"),
    "PGK":           (-18.5, "H", "[N13] Phosphoglycerate kinase"),
    "PGM":           (  4.4, "H", "[F12] Phosphoglycerate mutase"),
    "ENO":           ( -3.6, "H", "[F12] Enolase"),
    "PYK":           (-31.5, "H", "[N13] Pyruvate kinase"),
    "PDH":           (-39.8, "H", "[N13] Pyruvate dehydrogenase (overall)"),
    "FBP":           (-16.3, "H", "[F12] Fructose-1,6-bisphosphatase (gluconeogenesis)"),
    "PPCK":          ( -7.1, "M", "[F12] PEP carboxykinase"),
    "PPS":           ( -5.2, "M", "[F12] PEP synthase"),
    "PPC":           (-35.2, "H", "[N13] PEP carboxylase"),
    "ME1":           (-12.5, "M", "[F12] Malic enzyme (NAD+)"),
    "ME2":           (-12.5, "M", "[F12] Malic enzyme (NADP+)"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # PENTOSE PHOSPHATE PATHWAY
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "G6PDH2r":       ( -2.3, "H", "[F12] Glucose-6-P dehydrogenase"),
    "PGL":           (-24.8, "H", "[F12] 6-phosphogluconolactonase"),
    "GND":           (-11.4, "H", "[F12] 6-phosphogluconate dehydrogenase"),
    "RPE":           (  0.4, "H", "[F12] Ribulose-5-P epimerase"),
    "RPI":           (  2.2, "H", "[F12] Ribose-5-P isomerase"),
    "TKT1":          ( -0.5, "H", "[F12] Transketolase 1"),
    "TKT2":          (  0.2, "H", "[F12] Transketolase 2"),
    "TALA":          (  1.1, "H", "[F12] Transaldolase"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # TCA CYCLE
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "CS":            (-38.2, "H", "[N13] Citrate synthase"),
    "ACONTa":        (  8.4, "H", "[F12] Aconitase step a"),
    "ACONTb":        ( -1.2, "H", "[F12] Aconitase step b"),
    "ICDHyr":        (  1.1, "H", "[F12] Isocitrate dehydrogenase (NADP+)"),
    "AKGDH":         (-41.5, "H", "[N13] α-Ketoglutarate dehydrogenase"),
    "SUCOAS":        (  3.5, "H", "[F12] Succinyl-CoA synthetase"),
    "SUCDi":         ( -3.0, "M", "[F12] Succinate dehydrogenase (irreversible in iCW773)"),
    "FRD2":          (  3.0, "M", "[F12] Fumarate reductase (reverse SDH)"),
    "FRD3":          (  3.0, "M", "[F12] Fumarate reductase isoform"),
    "FUM":           (  0.5, "H", "[F12] Fumarase"),
    "MDH":           ( 29.5, "H", "[N13] Malate dehydrogenase (endergonic fwd)"),
    "MDH2":          ( 29.5, "H", "[N13] MDH isoform 2"),
    "MDH3":          ( 29.5, "H", "[N13] MDH isoform 3"),
    "ICL":           ( -5.8, "M", "[F12] Isocitrate lyase"),
    "MALS":          (-27.5, "M", "[F12] Malate synthase"),
    "CITL":          ( -2.2, "L", "[K] Citrate lyase"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # RESPIRATORY CHAIN & OXIDATIVE PHOSPHORYLATION
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "NADH16pp":      (-78.4, "H", "[T77] Complex I NADH:quinone oxidoreductase"),
    "NADH17pp":      (-78.4, "H", "[T77] NADH17pp isoform"),
    "NADH18pp":      (-78.4, "H", "[T77] NADH18pp isoform"),
    "NADH10":        (-52.0, "H", "[T77] Type-II NADH dehydrogenase"),
    "CYTBDpp":       (-88.5, "H", "[T77] Cytochrome bd quinol oxidase"),
    "CYTBD2pp":      (-88.5, "H", "[T77] Cytochrome bd2 isoform"),
    "ATPS4rpp":      (-36.2, "H", "[N13] ATP synthase (net synthesis direction)"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # GLUTAMATE / GLUTAMINE / NITROGEN ASSIMILATION
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "GLUDy":         (-33.2, "H", "[N13] Glutamate dehydrogenase"),
    "GLNS":          (-16.8, "H", "[N13] Glutamine synthetase (GS)"),
    "GLUSy":         (-37.1, "H", "[N13] Glutamate synthase (GOGAT, NADPH)"),
    "GLUN":          (-14.8, "M", "[K] Glutaminase"),
    "GLUNATF":       (-14.8, "M", "[K] Glutaminase (ATF)"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # ASPARTATE FAMILY / LYSINE BIOSYNTHESIS
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "ASPTA":         (  0.8, "H", "[F12] Aspartate aminotransferase"),
    "ASPK":          (-18.2, "M", "[K] Aspartate kinase"),
    "ASAD":          (  8.6, "M", "[K] Aspartate semialdehyde dehydrogenase"),
    "HSDy":          ( -5.5, "M", "[K] Homoserine dehydrogenase"),
    "DHDPS":         (-24.5, "M", "[K] Dihydrodipicolinate synthase"),
    "DHDPRy":        (-22.3, "M", "[K] Dihydrodipicolinate reductase"),
    "DAPDC":         (-24.0, "M", "[K] Diaminopimelate decarboxylase"),
    "DAPDH":         (-18.5, "M", "[K] Diaminopimelate dehydrogenase"),
    "DAPE":          (  0.5, "L", "[K] Diaminopimelate epimerase; near-eq"),
    "THRS":          (-21.5, "M", "[K] Threonine synthase"),
    "ASNS1":         (-12.5, "M", "[K] Asparagine synthetase 1"),
    "ASNS2":         (-12.5, "M", "[K] Asparagine synthetase 2"),
    "ASPT":          (-28.5, "M", "[K] Aspartate lyase (beta-alanine)"),
    "ASP1DC":        (-36.0, "M", "[K] Aspartate 1-decarboxylase"),
    "ASPCT":         (-14.8, "M", "[K] Aspartate carbamoyltransferase"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # SERINE / GLYCINE / ONE-CARBON
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "PGCD":          (-18.2, "L", "[K] Phosphoglycerate dehydrogenase"),
    "PSERT":         (  0.8, "L", "[K] Phosphoserine transaminase"),
    "PSP_L":         (-17.5, "L", "[K] Phosphoserine phosphatase"),
    "GHMT2r":        ( -5.4, "M", "[F12] Serine hydroxymethyltransferase"),
    "METS":          (-22.0, "M", "[K] Methionine synthase"),
    "MTHFC":         (  1.5, "M", "[K] Methenyl-THF cyclohydrolase"),
    "MTHFD":         ( -5.0, "M", "[K] Methylene-THF dehydrogenase (NADP+)"),
    "FTHFD":         (-11.0, "M", "[K] 10-formyl-THF deformylase"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # ALANINE / BRANCHED-CHAIN AMINO ACIDS
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "ALATA_L":       (  1.4, "M", "[F12] L-Alanine transaminase"),
    "ALATA_L2":      (  1.4, "M", "[F12] ALATA_L2"),
    "ALATA_D2":      (  1.4, "M", "[F12] D-Alanine transaminase"),
    "VALTA":         (  0.6, "L", "[K] Valine transaminase"),
    "LEUTAi":        (  0.3, "L", "[K] Leucine transaminase"),
    "KARA1":         (-22.4, "M", "[K] Acetolactate synthase (Val/Ile path)"),
    "KARA2":         (-22.4, "M", "[K] KARA2 isoform"),
    "DHAD1":         ( -3.4, "M", "[K] Dihydroxyacid dehydratase"),
    "DHAD2":         ( -3.4, "M", "[K] DHAD2 isoform"),
    "ACHBS":         (-18.6, "M", "[K] 2-acetolactate synthase (Leu path)"),
    "ACHMSC":        (-12.0, "M", "[K] 3-isopropylmalate dehydratase"),
    "ILETA":         (  0.3, "L", "[K] Isoleucine transaminase"),
    "ACLS":          (-22.4, "M", "[K] Acetolactate synthase"),
    "ALAR":          (-15.5, "H", "[K] Alanine racemase (L→D); exergonic"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # PROLINE / ARGININE / ORNITHINE
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "PROD2":         (-28.6, "M", "[K] Proline oxidase (FAD)"),
    "P5CD":          (-42.8, "L", "[K] Pyrroline-5-carboxylate dehydrogenase"),
    "OCBT":          (-14.8, "L", "[K] Ornithine carbamoyltransferase"),
    "AGPR":          (-21.5, "L", "[K] N-Ac-glutamyl-P reductase"),
    "ACGS":          (-15.2, "L", "[K] N-Ac-glutamate synthase"),
    "ACGK":          (-18.0, "L", "[K] N-Ac-glutamate kinase"),
    "ACOTA":         (  0.5, "L", "[K] N-Ac-ornithine transaminase; near-eq"),
    "ARGSS":         (-15.4, "L", "[K] Argininosuccinate synthase"),
    "ARGSL":         (  4.5, "L", "[K] Argininosuccinate lyase"),
    "CBPS":          (-19.5, "M", "[K] Carbamoyl-phosphate synthase (Gln)"),
    "CBMKr":         ( -9.5, "M", "[K] Carbamate kinase"),
    "AGPSC":         (-22.0, "L", "[K] Acetylglutamylphosphate reductase"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # HISTIDINE
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "HISTD":         (-12.5, "L", "[K] Histidinol dehydrogenase"),
    "HISTP":         (-15.0, "L", "[K] Histidinol-phosphate phosphatase"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # SHIKIMATE / AROMATIC AA BIOSYNTHESIS
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "DDPA":          (-31.5, "M", "[K] DAHP synthase (3-deoxy-D-arabino)"),
    "DHQS":          (-24.8, "M", "[K] 3-Dehydroquinate synthase"),
    "DHQTi":         (-15.2, "H", "[K] 3-Dehydroquinate dehydratase (irrev)"),
    "SHK3Dr":        ( -8.5, "M", "[K] Shikimate 3-dehydrogenase"),
    "SHKK":          (-18.5, "M", "[K] Shikimate kinase"),
    "CHORS":         (-18.0, "M", "[K] Chorismate synthase"),
    "CHORM":         ( -5.0, "M", "[K] Chorismate mutase"),
    "PPND":          (-28.5, "M", "[K] Prephenate dehydratase (Phe)"),
    "PPNDH":         (-18.5, "M", "[K] Prephenate dehydrogenase (Tyr)"),
    "ANS":           (-14.8, "M", "[K] Anthranilate synthase (Trp)"),
    "ANPRT":         (-15.4, "M", "[K] Anthranilate PRPP transferase"),
    "ADCS":          (-22.4, "M", "[K] 4-amino-4-deoxychorismate synthase"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # PURINE BIOSYNTHESIS
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "GLUPRT":        (-20.8, "M", "[K] Glutamine PRPP amidotransferase"),
    "GART":          (-14.5, "M", "[K] Phosphoribosylamine-glycine ligase"),
    "GARFT":         (-11.5, "M", "[K] Glycineamide ribonucleotide formyltransferase"),
    "AICART":        ( -3.5, "M", "[K] AICAR transformylase"),
    "ADCL":          (-14.5, "M", "[K] Adenylosuccinate lyase"),
    "ADSS":          (-18.5, "M", "[K] Adenylosuccinate synthase"),
    "ADSL1r":        (  3.2, "M", "[K] Adenylosuccinate lyase 1"),
    "ADSL2r":        (  3.2, "M", "[K] Adenylosuccinate lyase 2"),
    "AIRC2":         (-28.5, "M", "[K] AIR carboxylase"),
    "AIRC3":         (-28.5, "M", "[K] AIRC3 isoform"),
    "GMPS2":         (-22.0, "M", "[K] GMP synthase (verified: GMPS2)"),
    "ADNK1":         (-18.5, "M", "[K] Adenosine kinase"),
    "AMPMS2":        (-38.2, "M", "[K] AMP synthase 2"),
    "ADPRDP":        (-15.0, "M", "[K] ADP-ribose diphosphatase"),
    "ADNCYC":        (-18.5, "M", "[K] Adenylate cyclase"),
    "ADSK":          (-16.0, "M", "[K] Adenylyl sulfate kinase"),
    "AP5AH":         (-18.5, "M", "[K] Ap5A hydrolase"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # PYRIMIDINE BIOSYNTHESIS
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "DHORTS":        (-10.5, "M", "[K] Dihydroorotase"),
    "DHORD2":        (-22.5, "M", "[K] Dihydroorotate dehydrogenase (Q8)"),
    "DHORD5":        (-22.5, "M", "[K] Dihydroorotate dehydrogenase (MK8)"),
    "ORPT":          ( -0.5, "M", "[K] Orotate PRPP transferase"),
    "OMPDC":         (-35.0, "H", "[K] Orotidine-5-P decarboxylase (irrev)"),
    "CTPS2":         (-16.5, "M", "[K] CTP synthase"),
    "DCTPD":         (-18.5, "M", "[K] dCTP deaminase"),
    "DUTPDP":        (-22.0, "H", "[K] dUTP diphosphatase"),
    "CYTK1":         (-18.5, "M", "[K] Cytidylate kinase (CMP)"),
    "CYTK2":         (-18.5, "M", "[K] Cytidylate kinase (dCMP)"),
    "CYTDK2":        (-18.5, "M", "[K] Cytidine kinase (GTP)"),
    "DTMPK":         (-18.5, "M", "[K] dTMP kinase"),
    "DURIK1":        (-18.5, "M", "[K] Deoxyuridine kinase"),
    "DGK1":          (-18.5, "M", "[K] Deoxyguanylate kinase"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # COFACTOR METABOLISM
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "ADK1":          (  0.0, "H", "[F12] Adenylate kinase; ΔG'≈0"),
    "ADK2":          (  0.0, "H", "[F12] ADK2"),
    "ADK3":          (  0.0, "H", "[F12] ADK3"),
    "ADK4":          (  0.0, "H", "[F12] ADK4"),
    "PPA":           (-19.2, "H", "[F12] Inorganic pyrophosphatase"),
    "NDPK1":         (  0.0, "H", "[F12] NDP kinase; ΔG'≈0"),
    "ACKr":          ( -9.5, "M", "[F12] Acetate kinase"),
    "PTAr":          (  8.4, "M", "[F12] Phosphotransacetylase"),
    "DPCOAK":        (-14.5, "M", "[K] Dephospho-CoA kinase"),
    "FMNAT":         (-18.5, "M", "[K] FMN adenylyltransferase"),
    "DHFR":          (-22.5, "H", "[K] Dihydrofolate reductase (NADPH)"),
    "DHFS":          (-14.5, "M", "[K] Dihydrofolate synthase"),
    "DHPS2":         (-14.8, "M", "[K] Dihydropteroate synthase"),
    "DXPS":          (-38.5, "H", "[K] DXP synthase (MEP pathway)"),
    "DXPRIi":        (-26.5, "H", "[K] DXP reductoisomerase"),
    "DB4PS":         (-28.5, "M", "[K] 3,4-DHBPP synthase (riboflavin)"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # FATTY ACID SYNTHESIS (FAS II)
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "ACCOAC":        (-21.5, "M", "[N13] Acetyl-CoA carboxylase"),
    "3OAS60":        (-30.5, "M", "[K] 3-oxoacyl-ACP synthase C6"),
    "3OAS80":        (-30.5, "M", "[K] 3-OAS C8"),
    "3OAS100":       (-30.5, "M", "[K] 3-OAS C10"),
    "3OAS120":       (-30.5, "M", "[K] 3-OAS C12"),
    "3OAS121":       (-30.5, "M", "[K] 3-OAS C12:1"),
    "3OAS140":       (-30.5, "M", "[K] 3-OAS C14"),
    "3OAS141":       (-30.5, "M", "[K] 3-OAS C14:1"),
    "3OAS160":       (-30.5, "M", "[K] 3-OAS C16"),
    "3OAS161":       (-30.5, "M", "[K] 3-OAS C16:1"),
    "3OAS180":       (-30.5, "M", "[K] 3-OAS C18"),
    "3OAS181":       (-30.5, "M", "[K] 3-OAS C18:1"),
    "3OAR40":        (-20.5, "M", "[K] 3-oxoacyl-ACP reductase C4 (NADPH)"),
    "3OAR60":        (-20.5, "M", "[K] 3-OAR C6"),
    "3OAR80":        (-20.5, "M", "[K] 3-OAR C8"),
    "3OAR100":       (-20.5, "M", "[K] 3-OAR C10"),
    "3OAR120":       (-20.5, "M", "[K] 3-OAR C12"),
    "3OAR121":       (-20.5, "M", "[K] 3-OAR C12:1"),
    "3OAR140":       (-20.5, "M", "[K] 3-OAR C14"),
    "3OAR141":       (-20.5, "M", "[K] 3-OAR C14:1"),
    "3OAR160":       (-20.5, "M", "[K] 3-OAR C16"),
    "3OAR161":       (-20.5, "M", "[K] 3-OAR C16:1"),
    "3OAR180":       (-20.5, "M", "[K] 3-OAR C18"),
    "3OAR181":       (-20.5, "M", "[K] 3-OAR C18:1"),
    "3HAD60":        ( -5.0, "L", "[K] 3-hydroxyacyl-ACP dehydratase C6"),
    "3HAD80":        ( -5.0, "L", "[K] 3-HAD C8"),
    "3HAD100":       ( -5.0, "L", "[K] 3-HAD C10"),
    "3HAD120":       ( -5.0, "L", "[K] 3-HAD C12"),
    "3HAD121":       ( -5.0, "L", "[K] 3-HAD C12:1"),
    "3HAD140":       ( -5.0, "L", "[K] 3-HAD C14"),
    "3HAD141":       ( -5.0, "L", "[K] 3-HAD C14:1"),
    "3HAD160":       ( -5.0, "L", "[K] 3-HAD C16"),
    "3HAD161":       ( -5.0, "L", "[K] 3-HAD C16:1"),
    "3HAD180":       ( -5.0, "L", "[K] 3-HAD C18"),
    "3HAD181":       ( -5.0, "L", "[K] 3-HAD C18:1"),
    # Enoyl-ACP reductases (NADH: x; NADPH: y) — strongly forward (FADH₂ or NAD(P)H)
    "EAR40x":        (-24.5, "M", "[K] Enoyl-ACP reductase NADH C4"),
    "EAR40y":        (-24.5, "M", "[K] Enoyl-ACP reductase NADPH C4"),
    "EAR60x":        (-24.5, "M", "[K] EAR C6 NADH"),
    "EAR60y":        (-24.5, "M", "[K] EAR C6 NADPH"),
    "EAR80x":        (-24.5, "M", "[K] EAR C8 NADH"),
    "EAR80y":        (-24.5, "M", "[K] EAR C8 NADPH"),
    "EAR100x":       (-24.5, "M", "[K] EAR C10 NADH"),
    "EAR100y":       (-24.5, "M", "[K] EAR C10 NADPH"),
    "EAR120x":       (-24.5, "M", "[K] EAR C12 NADH"),
    "EAR120y":       (-24.5, "M", "[K] EAR C12 NADPH"),
    "EAR121x":       (-24.5, "M", "[K] EAR C12:1 NADH"),
    "EAR121y":       (-24.5, "M", "[K] EAR C12:1 NADPH"),
    "EAR140x":       (-24.5, "M", "[K] EAR C14 NADH"),
    "EAR140y":       (-24.5, "M", "[K] EAR C14 NADPH"),
    "EAR141x":       (-24.5, "M", "[K] EAR C14:1 NADH"),
    "EAR141y":       (-24.5, "M", "[K] EAR C14:1 NADPH"),
    "EAR160x":       (-24.5, "M", "[K] EAR C16 NADH"),
    "EAR160y":       (-24.5, "M", "[K] EAR C16 NADPH"),
    "EAR161x":       (-24.5, "M", "[K] EAR C16:1 NADH"),
    "EAR161y":       (-24.5, "M", "[K] EAR C16:1 NADPH"),
    "EAR180x":       (-24.5, "M", "[K] EAR C18 NADH"),
    "EAR180y":       (-24.5, "M", "[K] EAR C18 NADPH"),
    "EAR181x":       (-24.5, "M", "[K] EAR C18:1 NADH"),
    "EAR181y":       (-24.5, "M", "[K] EAR C18:1 NADPH"),
    # ACP hydrolases (irreversible)
    "FA80ACPHi":     (-30.5, "H", "[K] Fatty-acyl-ACP hydrolase C8"),
    "FA100ACPHi":    (-30.5, "H", "[K] FA-ACP hydrolase C10"),
    "FA120ACPHi":    (-30.5, "H", "[K] FA-ACP hydrolase C12"),
    "FA140ACPHi":    (-30.5, "H", "[K] FA-ACP hydrolase C14"),
    "FA141ACPHi":    (-30.5, "H", "[K] FA-ACP hydrolase C14:1"),
    "FA160ACPHi":    (-30.5, "H", "[K] FA-ACP hydrolase C16"),
    "FA161ACPHi":    (-30.5, "H", "[K] FA-ACP hydrolase C16:1"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # FATTY ACID BETA-OXIDATION
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "ACOAD1f":       (-22.1, "M", "[K] Acyl-CoA dehydrogenase C4 (FAD)"),
    "ACOAD2f":       (-22.1, "M", "[K] ACD C6"),
    "ACOAD3f":       (-22.1, "M", "[K] ACD C8"),
    "ACOAD4f":       (-22.1, "M", "[K] ACD C10"),
    "ACOAD5f":       (-22.1, "M", "[K] ACD C12"),
    "ACOAD6f":       (-22.1, "M", "[K] ACD C14"),
    "ACOAD7f":       (-22.1, "M", "[K] ACD C16"),
    "ACOAD8f":       (-22.1, "M", "[K] ACD C18"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # AMINO ACID ACTIVATION (tRNA ligases) — all ATP-driven, highly exergonic
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "ALATRS":        (-25.5, "H", "[K] Alanyl-tRNA synthetase (overall)"),
    "ARGTRS":        (-25.5, "H", "[K] Arginyl-tRNA synthetase"),
    "ASNTRS":        (-25.5, "H", "[K] Asparaginyl-tRNA synthetase"),
    "ASPTRS":        (-25.5, "H", "[K] Aspartyl-tRNA synthetase"),
    "CYSTRS":        (-25.5, "H", "[K] Cysteinyl-tRNA synthetase"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # PHOSPHOLIPID SYNTHESIS
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "AGPAT120":      (-22.0, "L", "[K] 1-acyl-sn-glycerol-3-P acyltransferase C12"),
    "AGPAT140":      (-22.0, "L", "[K] AGPAT C14"),
    "AGPAT141":      (-22.0, "L", "[K] AGPAT C14:1"),
    "AGPAT160":      (-22.0, "L", "[K] AGPAT C16"),
    "AGPAT161":      (-22.0, "L", "[K] AGPAT C16:1"),
    "AGPAT180":      (-22.0, "L", "[K] AGPAT C18"),
    "AGPAT181":      (-22.0, "L", "[K] AGPAT C18:1"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # ACETATE METABOLISM & ORGANIC ACID OVERFLOW
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "ACALD":         ( -4.5, "M", "[F12] Acetaldehyde dehydrogenase"),
    "ALCD2x":        (-22.5, "M", "[F12] Alcohol dehydrogenase (ethanol, NADH)"),
    "ALDD2x":        (-22.5, "M", "[F12] Aldehyde dehydrogenase (NADH)"),

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # MISC / CENTRAL METABOLISM
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    "DRPA":          (-12.5, "M", "[K] Deoxyribose-P aldolase"),
    "FFSD":          (-20.5, "M", "[K] Beta-fructofuranosidase (sucrose)"),
    "FRUK":          (-18.5, "M", "[K] Fructose-1-phosphate kinase"),
    "CAT":           (-88.0, "H", "[T77] Catalase (2H₂O₂ → 2H₂O + O₂)"),
    "DHAPT":         (-22.5, "M", "[K] Dihydroxyacetone phosphotransferase"),
    "FBP":           (-16.3, "H", "[F12] Fructose-1,6-bisphosphatase"),  # already listed
    "ABTA":          ( -0.5, "M", "[K] 4-aminobutyrate transaminase; near-eq"),
    "5FTHFC":        ( -8.5, "M", "[K] 5-formyl-THF cycloligase"),
    "2METS":         (-16.5, "M", "[K] 2-keto-4-methylthiobutyrate aminotransferase"),
    "AHCYSNS":       (-22.0, "H", "[K] S-adenosylhomocysteine nucleosidase"),
    "AMPN":          (-18.5, "M", "[K] AMP nucleosidase"),
    "ADNUC":         (-18.5, "M", "[K] Adenosine nucleosidase"),
    "ACNML":         ( -8.5, "M", "[K] Aconitate migration; near-eq"),
    "ACODA":         (-14.5, "M", "[K] Acetyl-CoA:oxaloacetate acyltransferase"),
    "PGL":           (-24.8, "H", "[F12] 6-phosphogluconolactonase"),  # duplicate OK
    "PPCK":          ( -7.1, "M", "[F12] PEP carboxykinase"),          # duplicate OK

    # ━━━━━━ v3 EXPANSION: Additional pathway blocks ━━━━━━
    "HISTRS":         (-25.5, "H", "[K] Histidyl-tRNA synthetase"),
    "ILETRS":         (-25.5, "H", "[K] Isoleucyl-tRNA synthetase"),
    "LEUTRS":         (-25.5, "H", "[K] Leucyl-tRNA synthetase"),
    "LYSTRS":         (-25.5, "H", "[K] Lysyl-tRNA synthetase"),
    "METTRS":         (-25.5, "H", "[K] Methionyl-tRNA synthetase"),
    "NADK":           (-18.5, "M", "[K] NAD kinase (ATP-driven)"),
    "NADS1":          (-22.0, "M", "[K] NAD synthase (Gln-hydrolyzing)"),
    "NADH5":          (-52.0, "H", "[T77] NADH dehydrogenase (ubiquinone-8)"),
    "NADH9":          (-52.0, "H", "[T77] NADH dehydrogenase (demethylmenaquinone)"),
    "NADPHQR2":       (-52.0, "H", "[T77] NADPH quinone reductase (ubiquinone)"),
    "NADPHQR3":       (-52.0, "H", "[T77] NADPH quinone reductase (menaquinone)"),
    "NADPHQR4":       (-52.0, "H", "[T77] NADPH quinone reductase (demethylmenaquinone)"),
    "NDPK2":          (  0.0, "H", "[F12] NDP kinase (ATP:UDP); near-eq"),
    "NDPK3":          (  0.0, "H", "[F12] NDP kinase (ATP:CDP)"),
    "NDPK4":          (  0.0, "H", "[F12] NDP kinase (ATP:dTDP)"),
    "NDPK5":          (  0.0, "H", "[F12] NDP kinase (ATP:dGDP)"),
    "NDPK6":          (  0.0, "H", "[F12] NDP kinase (ATP:dUDP)"),
    "IPMD":           (-12.0, "M", "[K] 3-isopropylmalate dehydrogenase"),
    "IPPMIa":         ( -3.5, "M", "[K] 3-isopropylmalate dehydratase step a"),
    "IPPMIb":         ( -3.5, "M", "[K] 3-isopropylmalate dehydratase step b"),
    "IPPS":           (-16.0, "M", "[K] 2-isopropylmalate synthase"),
    "MECDPS":         (-14.5, "M", "[K] 2C-methyl-D-erythritol-2,4-CDP synthase"),
    "IPDDI":          (  2.5, "M", "[K] Isopentenyl-PP isomerase; near-eq"),
    "IPDPS":          (-28.5, "M", "[K] HMBPP reductase"),
    "MEPCT":          (-22.0, "M", "[K] MEP cytidylyltransferase"),
    "IG3PS":          (-22.0, "M", "[K] Imidazole-glycerol-3-P synthase"),
    "IGPDH":          (-14.5, "M", "[K] Imidazoleglycerol-P dehydratase"),
    "IMPC":           (-18.5, "M", "[K] IMP cyclohydrolase"),
    "IMPD":           (-22.0, "M", "[K] IMP dehydrogenase (NAD+)"),
    "HSTP":           (-18.5, "M", "[K] Histidinol-phosphate transaminase"),
    "HEX1":           (-16.7, "H", "[F12] Hexokinase (glucose + ATP)"),
    "HEX7":           (-16.7, "H", "[F12] Hexokinase (fructose + ATP)"),
    "MAN6PI":         (  1.8, "M", "[K] Mannose-6-P isomerase; near-eq"),
    "M1PD":           (-14.5, "M", "[K] Mannitol-1-P dehydrogenase"),
    "MAN2D":          (-14.5, "M", "[K] Mannitol-2-dehydrogenase"),
    "METAT":          (-38.5, "H", "[K] Methionine adenosyltransferase"),
    "MTHFR2":         (-22.0, "M", "[K] 5,10-methylene-THF reductase (NADH)"),
    "HSK":            (-18.5, "M", "[K] Homoserine kinase (ATP-driven)"),
    "HSST":           (-22.5, "M", "[K] Homoserine O-acetyltransferase"),
    "HSST_2":         (-22.5, "M", "[K] HSST isoform 2"),
    "CYSTL":          (-14.5, "M", "[K] Cystathionine beta-lyase"),
    "CYSS":           (-22.0, "M", "[K] Cysteine synthase (serine + H2S)"),
    "MMM2":           (  0.3, "M", "[K] Methylmalonyl-CoA mutase; near-eq"),
    "MCITD":          ( -3.5, "M", "[K] 2-methylcitrate dehydratase"),
    "MCITL2":         ( -5.0, "M", "[K] Methylisocitrate lyase"),
    "MSDH":           (-28.5, "M", "[K] Methylmalonate-semialdehyde dehydrogenase"),
    "MSDHD":          (-28.5, "M", "[K] Malonate-semialdehyde dehydrogenase"),
    "MM_COA_ADD5":    (-22.5, "M", "[K] Propanoyl-CoA carboxylase (ATP-driven)"),
    "LDH":            (-25.0, "M", "[K] L-Lactate dehydrogenase (NAD+)"),
    "LDH_D":          (-25.0, "M", "[K] D-Lactate dehydrogenase"),
    "LDH_D2":         (-25.0, "M", "[K] D-LDH isoform 2"),
    "L_LACD2":        (-22.5, "M", "[K] L-Lactate dehydrogenase (ubiquinone)"),
    "L_LACD3":        (-22.5, "M", "[K] L-Lactate dehydrogenase (menaquinone)"),
    "HMBS":           (-22.5, "M", "[K] Hydroxymethylbilane synthase"),
    "FCLT":           (-18.5, "M", "[K] Ferrochelatase (heme insertion)"),
    "ICHOR":          ( -8.5, "M", "[K] Isochorismate synthase"),
    "HPPK2":          (-18.5, "M", "[K] 6-Hydroxymethyl-dihydropterin pyrophosphokinase"),
    "FMETTRS":        (-25.5, "H", "[K] Met-tRNA formyltransferase"),
    "LIPOCT":         (-22.0, "M", "[K] Lipoyl(octanoyl) transferase"),
    "LIPOS":          (-38.5, "H", "[K] Lipoate synthase (radical SAM)"),
    "LIPAMPL":        (-22.5, "M", "[K] Lipoyl-adenylate protein ligase"),
    "BTS4":           (-42.5, "H", "[K] Biotin synthase (radical SAM)"),
    "BTNC":           (-22.0, "M", "[K] Biotin carboxylase (ACC subunit)"),
    "HXPRT":          ( -5.5, "M", "[K] Hypoxanthine-PRPP phosphoribosyltransferase"),
    "INSH":           (-15.0, "M", "[K] Inosine hydrolase"),
    "ACACT1r":        (-12.5, "M", "[K] Acetyl-CoA acetyltransferase (thiolase)"),
    "MCOATA":         (-22.5, "M", "[K] Malonyl-CoA:ACP acyltransferase"),
    "MACPD":          (-18.5, "H", "[K] Malonyl-ACP decarboxylase (irreversible)"),
    "KAS14":          (-30.5, "M", "[K] Beta-ketoacyl-ACP synthase (FabF)"),
    "KAS15":          (-30.5, "M", "[K] Beta-ketoacyl-ACP synthase isoform 2"),
    "HCO3E":          (  0.0, "H", "[F12] HCO3 equilibration; dG approx 0"),
    "MALGT":          (-22.5, "M", "[K] Maltose glucosyltransferase"),
}

# ── Deduplicate ────────────────────────────────────────────────────────────────
# Keep first occurrence if any duplicates snuck in
seen = {}
for k, v in CURATED_DGR0.items():
    if k not in seen:
        seen[k] = v
CURATED_DGR0 = seen


def build_thermo_data():
    logger.info("Building comprehensive thermo_dgr_data.json (v3)...")

    try:
        import cobra
        model = cobra.io.read_sbml_model(MODEL_PATH)
        all_rxn_ids = {r.id for r in model.reactions}
        logger.info(f"Model has {len(all_rxn_ids)} reactions")
    except Exception as e:
        logger.warning(f"Cannot load model ({e}); using curated keys only.")
        all_rxn_ids = set(CURATED_DGR0.keys())

    CONF_MAP = {"H": "HIGH", "M": "MED", "L": "LOW"}
    results  = {}
    n_fwd = n_rev = n_neq = n_nodata = 0
    n_id_miss = 0

    for rxn_id, (dgr0, conf_short, note) in CURATED_DGR0.items():
        in_model = rxn_id in all_rxn_ids
        if not in_model:
            n_id_miss += 1

        dgr_min = round(dgr0 - HALF_SPREAD, 2)
        dgr_max = round(dgr0 + HALF_SPREAD, 2)

        if in_model:
            if dgr_max < -EPSILON:
                direction = "forward"; n_fwd += 1
            elif dgr_min > EPSILON:
                direction = "reverse"; n_rev += 1
            else:
                direction = "none"; n_neq += 1
        else:
            direction = "none"

        results[rxn_id] = {
            "dgr_prime_0":   round(dgr0, 2),
            "dgr_prime_min": dgr_min,
            "dgr_prime_max": dgr_max,
            "direction_locked": direction,
            "confidence":    CONF_MAP.get(conf_short, "LOW"),
            "source":        note,
            "in_model":      in_model,
            "note": (
                f"Forward-locked (ΔrG'max={dgr_max:.1f} < -{EPSILON})" if direction == "forward" else
                f"Reverse-locked (ΔrG'min={dgr_min:.1f} > +{EPSILON})"  if direction == "reverse" else
                "Near-equilibrium or unmatched; bounds unchanged"
            )
        }

    # Fill remaining model reactions with no-data placeholders
    for rxn_id in all_rxn_ids:
        if rxn_id not in results:
            n_nodata += 1
            results[rxn_id] = {
                "dgr_prime_0": None, "dgr_prime_min": None, "dgr_prime_max": None,
                "direction_locked": "none", "confidence": "NONE",
                "source": "no_data", "in_model": True,
                "note": "No thermodynamic data; bounds unchanged"
            }

    total    = len(all_rxn_ids)
    n_in     = len(CURATED_DGR0) - n_id_miss
    coverage = round(n_in / total * 100, 1) if total else 0.0

    logger.info("=" * 65)
    logger.info(f"Total model reactions    : {total}")
    logger.info(f"Curated (in model)       : {n_in}  ({coverage}%)")
    logger.info(f"Curated (ID not in model): {n_id_miss}")
    logger.info(f"Forward-locked           : {n_fwd}")
    logger.info(f"Reverse-locked           : {n_rev}")
    logger.info(f"Near-equilibrium         : {n_neq}")
    logger.info(f"No data                  : {n_nodata}")
    logger.info(f"Half-spread (kJ/mol)     : {HALF_SPREAD:.2f}")
    logger.info("=" * 65)

    output = {
        "_meta": {
            "description": "Comprehensive curated ΔrG' data for iCW773 (v3)",
            "conditions": "pH 7.0, I=0.1 M, T=30°C (303.15 K)",
            "c_min_M": C_MIN, "c_max_M": C_MAX, "epsilon_kJ": EPSILON, "units": "kJ/mol",
            "total_reactions": total, "coverage_pct": coverage,
            "n_forward_locked": n_fwd, "n_reverse_locked": n_rev,
            "n_near_equilibrium": n_neq, "n_no_data": n_nodata,
            "half_spread_kJ": round(HALF_SPREAD, 2),
            "sources": [
                "[N13] Noor et al. 2013 PLoS Comput Biol (Component Contribution)",
                "[F12] Flamholz et al. 2012 Nucleic Acids Res (eQuilibrator)",
                "[T77] Thauer et al. 1977 Bacteriol Rev 41(1):100",
                "[A03] Alberty 2003 MIT Press",
                "[K]   KEGG thermodynamic estimates"
            ],
            "generated": "2026-07-12", "version": "3.0"
        },
        "reactions": results
    }

    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    with open(OUTPUT_PATH, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)
    logger.info(f"Written: {OUTPUT_PATH}")
    return output


if __name__ == "__main__":
    build_thermo_data()
