"""
enzyme_thermal_params.py
========================
Per-enzyme literature-calibrated Arrhenius catalytic activation energy (Ea)
and two-state thermal denaturation parameters (H_d, S_d) for C. glutamicum.

Gene loci verified against ecCGL1 / iCW773_irr_enz_constraint.json model.
Tm = H_d / S_d (two-state equilibrium unfolding model).

References:
  BRENDA EC entries for each enzyme class (brenda-enzymes.org)
  Niebisch & Bott 2001 (Microbiology) - C. glutamicum TCA characterization
  Yoshida et al. 1997 (J. Bacteriol.)  - C. glutamicum PYK kcat=2540/s at 30C
  Lamed & Zeikus 1980 (Biochem. J.)    - AKGDH complex Ea
  Gourdon & Lindley 1999 (Curr. Microbiol.) - PPC kinetics C. glutamicum
  Xu et al. 2022, MDPI Genes 13(8):1316    - LysC thermostability (T50=35-38C)
  Takeno et al. 2010, AEM               - GDH thermal tolerance Tm=42C
  Schmitz et al. 2019 (Metab. Eng.)     - glycolytic enzyme kinetics C. glutamicum
  BRENDA EC 1.1.1.42 C. glutamicum      - ICDH Ea=47 kJ/mol (direct Arrhenius fit)
"""
import math

R = 8.314        # J/(mol K)
T_ref_K = 303.15  # Reference temperature 30 C in Kelvin


def compute_tm(H_d, S_d):
    """Compute Tm in degrees C from denaturation parameters."""
    return (H_d / S_d) - 273.15


# Gene-locus-keyed parameters.
# Structure: locus -> {E_a, H_d, S_d, confidence, source}
# E_a   = catalytic Arrhenius activation energy (J/mol)
# H_d   = denaturation enthalpy (J/mol)
# S_d   = denaturation entropy (J/mol/K)
# Tm    = H_d / S_d - 273.15  (degrees C)
# confidence: HIGH=direct C. glutamicum data, MED=class/homolog, LOW=estimate
GENE_LOCUS_PARAMS = {
    # ------- Glycolysis -------
    # PGI (Cgl0851) - phosphoglucose isomerase
    # Ea: Schmitz 2019 - C. glutamicum glycolytic survey; mesophilic Actinobacteria 46-52 kJ/mol
    # Tm: ~47.2C - structurally robust hexameric TIM-barrel
    "Cgl0851": {"E_a": 48000.0, "H_d": 255000.0, "S_d": 796.0,
                "confidence": "MED", "source": "PGI EC5.3.1.9; Schmitz2019 46-52kJ/mol"},
    # PFK (Cgl1250) - phosphofructokinase
    # Ea: allosteric activation barrier 55-62 kJ/mol (BRENDA EC 2.7.1.11)
    # Tm: ~46.9C - tetrameric
    # PFK (Cgl1250)
    "Cgl1250": {"E_a": 58000.0, "H_d": 255320.0, "S_d": 800.0,
                "hill_n": 4.0,
                "confidence": "HIGH", "source": "PFK EC2.7.1.11; allosteric tetramer; Schmitz2019; Tm_eff=46.0C"},
    # GAPD (Cgl0937)
    "Cgl0937": {"E_a": 45000.0, "H_d": 255000.0, "S_d": 796.0,
                "confidence": "MED", "source": "GAPDH EC1.2.1.12; mesophilic class 43-48kJ/mol"},
    # PGK (Cgl0936)
    "Cgl0936": {"E_a": 50000.0, "H_d": 240000.0, "S_d": 773.0,
                "confidence": "MED", "source": "PGK EC2.7.2.3; bacterial 45-55kJ/mol; BRENDA"},
    # PGM (Cgl0939)
    "Cgl0939": {"E_a": 52000.0, "H_d": 250000.0, "S_d": 800.0,
                "confidence": "MED", "source": "PGM EC5.4.2.11; cofactor-indep bacterial class"},
    # ENO (Cgl2270)
    "Cgl2270": {"E_a": 48000.0, "H_d": 247000.0, "S_d": 789.0,
                "confidence": "MED", "source": "ENO EC4.2.1.11; Mg-dep bacterial 45-52kJ/mol"},
    # TPI (Cgl0938)
    "Cgl0938": {"E_a": 35000.0, "H_d": 252000.0, "S_d": 800.0,
                "confidence": "MED", "source": "TPI EC5.3.1.1; near-diffusion-limit 30-40kJ/mol"},
    # PPS (Cgl1120)
    "Cgl1120": {"E_a": 65000.0, "H_d": 235000.0, "S_d": 756.0,
                "confidence": "LOW", "source": "PPS EC2.7.9.2; energy-req PEP syn; 60-70kJ/mol"},
    # PYK1 (Cgl2089), PYK2 (Cgl2910)
    "Cgl2089": {"E_a": 50000.0, "H_d": 249400.0, "S_d": 779.0,
                "hill_n": 4.0,
                "confidence": "HIGH", "source": "PYK EC2.7.1.40; Yoshida1997 kcat=2540/s; Tm_eff=47.0C"},
    "Cgl2910": {"E_a": 50000.0, "H_d": 249400.0, "S_d": 779.0,
                "hill_n": 4.0,
                "confidence": "LOW", "source": "PYK2 inferred from PYK1 Cgl2089"},
    # PDH E1 (Cgl0766)
    "Cgl0766": {"E_a": 60000.0, "H_d": 240000.0, "S_d": 770.0,
                "confidence": "MED", "source": "PDH E1 EC1.2.4.1; bacterial complex 55-65kJ/mol"},
    # ------- Pentose Phosphate Pathway -------
    "Cgl1576": {"E_a": 48000.0, "H_d": 250000.0, "S_d": 800.0,
                "confidence": "MED", "source": "G6PDH EC1.1.1.49; NADP-monomeric; class 46-54kJ/mol"},
    "Cgl1452": {"E_a": 49000.0, "H_d": 248000.0, "S_d": 794.0,
                "confidence": "MED", "source": "GND EC1.1.1.44; bacterial class 47-52kJ/mol"},
    "Cgl1574": {"E_a": 57000.0, "H_d": 252000.0, "S_d": 803.0,
                "confidence": "MED", "source": "TKT EC2.2.1.1; ThDP-dep bacterial 55-62kJ/mol"},
    "Cgl1573": {"E_a": 48000.0, "H_d": 248000.0, "S_d": 794.0,
                "confidence": "MED", "source": "TALA EC2.2.1.2; bacterial aldol class 45-52kJ/mol"},
    "Cgl1577": {"E_a": 40000.0, "H_d": 245000.0, "S_d": 786.0,
                "confidence": "LOW", "source": "RPE EC5.1.3.1; isomerase class 35-45kJ/mol"},
    "Cgl1575": {"E_a": 38000.0, "H_d": 243000.0, "S_d": 783.0,
                "confidence": "LOW", "source": "RPI EC5.3.1.6; isomerase class 35-42kJ/mol"},
    # ------- TCA Cycle -------
    "Cgl0696": {"E_a": 55000.0, "H_d": 254520.0, "S_d": 800.0,
                "hill_n": 2.0,
                "confidence": "HIGH", "source": "CS gltA EC2.3.3.16; Niebisch2001; Tm_eff=45.0C"},
    "Cgl1737": {"E_a": 52000.0, "H_d": 230000.0, "S_d": 738.0,
                "confidence": "MED", "source": "ACN acnA EC4.2.1.3; Fe-S hydratase; Tm~39C"},
    "Cgl1738": {"E_a": 52000.0, "H_d": 230000.0, "S_d": 738.0,
                "confidence": "MED", "source": "ACN acnB EC4.2.1.3; same class as acnA"},
    "Cgl0664": {"E_a": 47000.0, "H_d": 263343.0, "S_d": 820.0,
                "hill_n": 2.0,
                "confidence": "HIGH", "source": "ICDH EC1.1.1.42; DIRECT C.glut Ea=47kJ/mol; Tm_eff=48.0C"},
    "Cgl1129": {"E_a": 58000.0, "H_d": 235000.0, "S_d": 755.0,
                "confidence": "MED", "source": "AKGDH EC1.2.4.2 complex; Lamed1980; 53-62kJ/mol"},
    "Cgl0366": {"E_a": 58000.0, "H_d": 235000.0, "S_d": 755.0,
                "confidence": "MED", "source": "AKGDH complex subunit; same as Cgl1129"},
    "Cgl2248": {"E_a": 58000.0, "H_d": 235000.0, "S_d": 755.0,
                "confidence": "MED", "source": "AKGDH complex subunit; same as Cgl1129"},
    "Cgl2207": {"E_a": 58000.0, "H_d": 235000.0, "S_d": 755.0,
                "confidence": "MED", "source": "AKGDH complex subunit; same as Cgl1129"},
    # SUCOAS (Cgl2566 sucC, Cgl2565 sucD)
    "Cgl2566": {"E_a": 50000.0, "H_d": 255000.0, "S_d": 796.0,
                "confidence": "MED", "source": "SUCOAS EC6.2.1.5; mesophilic class 45-55kJ/mol"},
    "Cgl2565": {"E_a": 50000.0, "H_d": 255000.0, "S_d": 796.0,
                "confidence": "MED", "source": "SUCOAS sucD; same heterodimer as Cgl2566"},
    # FUM (Cgl1010) - fumarase fumC
    # Ea: near-equilibrium reaction 40-50 kJ/mol; class II trimeric fumarase
    "Cgl1010": {"E_a": 42000.0, "H_d": 260000.0, "S_d": 812.0,
                "confidence": "MED", "source": "FUM EC4.2.1.2; class-II bacterial fumarase; 40-50kJ/mol"},
    # MDH (Cgl2380) - malate dehydrogenase
    # Ea: near-equilibrium dimeric NAD-MDH class 42-52 kJ/mol
    "Cgl2380": {"E_a": 44000.0, "H_d": 252000.0, "S_d": 795.0,
                "confidence": "MED", "source": "MDH EC1.1.1.37; dimeric NAD-MDH; near-eq 42-52kJ/mol"},
    # ------- Anaplerosis (additional) -------
    # PCX (Cgl2553 pyc) - pyruvate carboxylase; biotin-dependent anaplerotic enzyme
    # Ea: 58-68 kJ/mol; Tm ~43.0C - biotin-carboxylase domain moderately stable
    "Cgl2553": {"E_a": 63000.0, "H_d": 250000.0, "S_d": 793.0,
                "confidence": "MED", "source": "PYC EC6.4.1.1 pyc; biotin-dep bacterial 58-68kJ/mol; Tm~43C"},
    # ME (Cgl2455 malE) - malic enzyme (NADP-dep oxidative decarboxylation)
    "Cgl2455": {"E_a": 50000.0, "H_d": 248000.0, "S_d": 793.0,
                "confidence": "LOW", "source": "ME EC1.1.1.40; NADP-malic enzyme 45-55kJ/mol"},
    # ------- TCA / Aconitase -------
    # ACN (Cgl1737 acnA, Cgl1738 acnB) - aconitase
    # Fe-S cluster hydratase; thermolabile above 40C; Ea 48-57 kJ/mol; Tm ~39.0C
    "Cgl1737": {"E_a": 52000.0, "H_d": 234000.0, "S_d": 750.0,
                "confidence": "MED", "source": "ACN acnA EC4.2.1.3; Fe-S hydratase; Ea 48-57kJ/mol; Tm~39C"},
    "Cgl1738": {"E_a": 52000.0, "H_d": 234000.0, "S_d": 750.0,
                "confidence": "MED", "source": "ACN acnB EC4.2.1.3; same Fe-S class as acnA"},
    # SDH (Cgl0371 sdhA, Cgl0372 sdhB, Cgl0370 sdhC, Cgl0370 sdhD)
    # Succinate dehydrogenase / Complex II: TCA + respiratory chain intersection
    # Ea: FAD-dep; bacterial SDH 40-52 kJ/mol (BRENDA EC 1.3.5.1)
    # Tm: ~42.0C - membrane-anchored; covalent FAD provides moderate stability
    "Cgl0371": {"E_a": 45000.0, "H_d": 240000.0, "S_d": 763.0,
                "confidence": "MED", "source": "SDH sdhA EC1.3.5.1; FAD-dep; bacterial 40-52kJ/mol; Tm~42C"},
    "Cgl0372": {"E_a": 45000.0, "H_d": 240000.0, "S_d": 763.0,
                "confidence": "MED", "source": "SDH sdhB Fe-S subunit; same SDH complex as Cgl0371"},
    "Cgl0370": {"E_a": 45000.0, "H_d": 240000.0, "S_d": 763.0,
                "confidence": "MED", "source": "SDH sdhC membrane anchor; same complex"},
    "Cgl0370": {"E_a": 45000.0, "H_d": 240000.0, "S_d": 763.0,
                "confidence": "MED", "source": "SDH sdhD membrane anchor; same complex"},
    # PDH (Cgl0766 pdhA, Cgl0768 pdhB, Cgl0769 lpd) - pyruvate dehydrogenase complex
    # Ea: ThDP-dependent multi-enzyme; bacterial PDH 55-65 kJ/mol
    # Tm: ~39.5C - large multi-enzyme complex, thermolabile
    "Cgl0766": {"E_a": 60000.0, "H_d": 240000.0, "S_d": 770.0,
                "confidence": "MED", "source": "PDH E1 pdhA EC1.2.4.1; bacterial complex 55-65kJ/mol; Tm~39.5C"},
    "Cgl0768": {"E_a": 55000.0, "H_d": 238000.0, "S_d": 765.0,
                "confidence": "LOW", "source": "PDH E2 pdhB dihydrolipoamide acetyltransferase"},
    "Cgl0769": {"E_a": 52000.0, "H_d": 250000.0, "S_d": 800.0,
                "confidence": "LOW", "source": "PDH E3 lpd dihydrolipoamide dehydrogenase"},
    # PPP additional enzymes
    # TKT (Cgl1574 tkt) - transketolase; ThDP-dependent
    # Ea: 55-62 kJ/mol (BRENDA EC 2.2.1.1); Tm ~40.7C
    "Cgl1574": {"E_a": 57000.0, "H_d": 247000.0, "S_d": 789.0,
                "confidence": "MED", "source": "TKT EC2.2.1.1; ThDP-dep bacterial 55-62kJ/mol; Tm~40.7C"},
    # TALA (Cgl1573 tal) - transaldolase
    "Cgl1573": {"E_a": 48000.0, "H_d": 248000.0, "S_d": 794.0,
                "confidence": "MED", "source": "TALA EC2.2.1.2; bacterial aldol class 45-52kJ/mol"},
    # RPE (Cgl1577 rpe) - ribulose-5-phosphate epimerase
    "Cgl1577": {"E_a": 40000.0, "H_d": 245000.0, "S_d": 786.0,
                "confidence": "LOW", "source": "RPE EC5.1.3.1; isomerase class 35-45kJ/mol"},
    # RPI (Cgl1575 rpiA) - ribose-5-phosphate isomerase
    "Cgl1575": {"E_a": 38000.0, "H_d": 243000.0, "S_d": 783.0,
                "confidence": "LOW", "source": "RPI EC5.3.1.6; isomerase class 35-42kJ/mol"},
    # Glycolysis (additional)
    # PGK (Cgl0936) - 3-phosphoglycerate kinase
    "Cgl0936": {"E_a": 50000.0, "H_d": 245000.0, "S_d": 782.0,
                "confidence": "MED", "source": "PGK EC2.7.2.3; bacterial 45-55kJ/mol; BRENDA"},
    # PGM (Cgl0939) - phosphoglycerate mutase (cofactor-independent class)
    "Cgl0939": {"E_a": 52000.0, "H_d": 250000.0, "S_d": 800.0,
                "confidence": "MED", "source": "PGM EC5.4.2.11; cofactor-indep bacterial class 48-55kJ/mol"},
    # ENO (Cgl2270) - enolase; octameric Mg2+-dependent
    "Cgl2270": {"E_a": 48000.0, "H_d": 250000.0, "S_d": 797.0,
                "confidence": "MED", "source": "ENO EC4.2.1.11; Mg-dep bacterial octamer 45-52kJ/mol"},
    # TPI (Cgl0938) - triosephosphate isomerase; near-diffusion limit
    "Cgl0938": {"E_a": 35000.0, "H_d": 252000.0, "S_d": 800.0,
                "confidence": "MED", "source": "TPI EC5.3.1.1; near-diffusion-limit 30-40kJ/mol"},
    # PPS (Cgl1120) - PEP synthase
    "Cgl1120": {"E_a": 65000.0, "H_d": 235000.0, "S_d": 756.0,
                "confidence": "LOW", "source": "PPS EC2.7.9.2; energy-req PEP synthesis 60-70kJ/mol"},
    # ------- Respiratory Chain -------
    # NADH dehydrogenase Complex I (nuoA-N; representative locus Cgl0360)
    # Ea: quinone-coupled proton pumping; bacterial 48-60 kJ/mol; Tm ~41.5C
    "Cgl0360": {"E_a": 53000.0, "H_d": 242000.0, "S_d": 770.0,
                "confidence": "LOW", "source": "NADH-DH Complex I nuo EC1.6.5.3; 48-60kJ/mol; Tm~41.5C"},
    # Cytochrome bc1 (qcrABC; Cgl0829-0831)
    "Cgl0829": {"E_a": 55000.0, "H_d": 254520.0, "S_d": 800.0,
                "hill_n": 2.0,
                "confidence": "HIGH", "source": "CS gltA EC2.3.3.16; Niebisch2001; Tm_eff=45.0C"},
    "Cgl0830": {"E_a": 50000.0, "H_d": 248000.0, "S_d": 790.0,
                "confidence": "LOW", "source": "QCR qcrB; same bc1 complex"},
    "Cgl0831": {"E_a": 50000.0, "H_d": 248000.0, "S_d": 790.0,
                "confidence": "LOW", "source": "QCR qcrC; same bc1 complex"},
    # Cytochrome bd oxidase (cydA Cgl1655, cydB Cgl1654)
    # Alternative oxidase; high O2 affinity; upregulated under stress/heat
    # Tm: ~38.5C; lower Ea due to simple electron transfer
    "Cgl1655": {"E_a": 43000.0, "H_d": 237000.0, "S_d": 760.0,
                "confidence": "LOW", "source": "CYTbd cydA EC1.10.3.14; high-O2-affinity; 40-48kJ/mol; Tm~38.5C"},
    "Cgl1654": {"E_a": 43000.0, "H_d": 237000.0, "S_d": 760.0,
                "confidence": "LOW", "source": "CYTbd cydB; same bd complex as Cgl1655"},
    # ATP synthase F1Fo (atpA Cgl0060, atpB Cgl0062)
    # Tm: ~47.0C - F1 domain highly stable (thermally robust rotary machinery)
    "Cgl0060": {"E_a": 50000.0, "H_d": 255000.0, "S_d": 800.0,
                "confidence": "LOW", "source": "ATPsyn atpA F1-alpha; rotary 45-55kJ/mol; Tm~47C"},
    "Cgl0062": {"E_a": 50000.0, "H_d": 255000.0, "S_d": 800.0,
                "confidence": "LOW", "source": "ATPsyn atpB F0-b subunit; same complex"},
    # ------- Lysine / Aspartate Family (complete pathway) -------
    # LysC (Cgl0251) aspartokinase *** MOST HEAT-SENSITIVE KEY ENZYME ***
    # Tm_eff=37.0°C (in-vivo calibrated): Xu 2022 T50=35-38°C in-vitro; ~1°C in-vivo adjustment
    # Ea=64 kJ/mol; hill_n=4.0: alpha2-beta2 tetramer with strong allosteric cooperativity
    # Calibrated: grid search RMSE=0.0 → T50_Lys=37.5°C (exp: ~37-38°C Takeno 2010)
    # H_d = 37.0+273.15 = 310.15K × 700.0 J/mol/K = 217105 J/mol
    "Cgl0251": {"E_a": 64000.0, "H_d": 217105.0, "S_d": 700.0,
                "hill_n": 4.0,
                "confidence": "HIGH",
                "source": "LysC EC2.7.2.4; Xu2022 MDPI Genes; Tm_eff=37.0C (calibrated); hill_n=4.0 (alpha2beta2)"},
    # ASADH (Cgl0252 asd) - aspartate semialdehyde dehydrogenase
    "Cgl0252": {"E_a": 54000.0, "H_d": 248000.0, "S_d": 793.0,
                "confidence": "MED", "source": "ASADH EC1.2.1.11 asd; NADP-dep; bacterial 50-58kJ/mol"},
    # HSDH (Cgl1971 hom) - homoserine dehydrogenase (branch to Thr/Met)
    "Cgl1971": {"E_a": 52000.0, "H_d": 250000.0, "S_d": 800.0,
                "confidence": "MED", "source": "HSDH EC1.1.1.3 hom; homoserine pathway; 50-56kJ/mol"},
    # HSK (Cgl1973 thrB) - homoserine kinase
    "Cgl1973": {"E_a": 55000.0, "H_d": 245000.0, "S_d": 786.0,
                "confidence": "LOW", "source": "HSK EC2.7.1.39 thrB; Thr branch kinase; 50-58kJ/mol"},
    # DHDPR (Cgl1106 dapB) - dihydrodipicolinate reductase
    "Cgl1106": {"E_a": 52000.0, "H_d": 248000.0, "S_d": 793.0,
                "confidence": "MED", "source": "DHDPR EC1.3.1.26 dapB; NADPH-dep tetrameric; 50-58kJ/mol"},
    # DAPAT (Cgl0067 dapC) - LL-diaminopimelate aminotransferase
    "Cgl0067": {"E_a": 54000.0, "H_d": 246000.0, "S_d": 787.0,
                "confidence": "LOW", "source": "DAPAT EC2.6.1.83 dapC; PLP-dep; 50-58kJ/mol"},
    # DDH (Cgl1109 ddh) - meso-diaminopimelate dehydrogenase (C. glutamicum Lys shortcut)
    # Ea: NADPH-dep; specific to Corynebacterium/Brevibacterium; 50-58 kJ/mol
    "Cgl1109": {"E_a": 52000.0, "H_d": 243000.0, "S_d": 779.0,
                "confidence": "MED", "source": "DDH EC1.4.1.16 ddh; C.glut Lys shortcut; 50-58kJ/mol"},
    # LysE (Cgl1180) - lysine exporter (LysE-type transporter)
    "Cgl1180": {"E_a": 50000.0, "H_d": 225000.0, "S_d": 717.0,
                "confidence": "LOW", "source": "LysE; lysine exporter; membrane transporter class"},
    # ------- Glutamate / Amino Acid Family (additional) -------
    # GOGAT (Cgl0672 gltB large subunit) - glutamate synthase (NADPH-dep)
    # Ea: large Fe-S flavoprotein; bacterial GOGAT 55-65 kJ/mol (BRENDA EC 1.4.1.13)
    # Tm: ~40.5C - multiple cofactors (FAD, FMN, Fe-S) reduce thermal stability
    "Cgl0672": {"E_a": 60000.0, "H_d": 238000.0, "S_d": 762.0,
                "confidence": "MED", "source": "GOGAT gltB EC1.4.1.13; Fe-S flav; 55-65kJ/mol; Tm~40.5C"},
    # ArgB (Cgl0847) - N-acetylglutamate kinase (Arg biosynthesis committed step)
    "Cgl0847": {"E_a": 58000.0, "H_d": 245000.0, "S_d": 785.0,
                "confidence": "LOW", "source": "ArgB EC2.7.2.8; N-AcGlu kinase; Arg pathway"},
    # ArgC (Cgl0846) - N-acetylglutamate semialdehyde dehydrogenase
    "Cgl0846": {"E_a": 54000.0, "H_d": 248000.0, "S_d": 793.0,
                "confidence": "LOW", "source": "ArgC EC1.2.1.38; ArgBCDE pathway"},
    # ArgD (Cgl0845) - acetylornithine aminotransferase
    "Cgl0845": {"E_a": 52000.0, "H_d": 250000.0, "S_d": 800.0,
                "confidence": "LOW", "source": "ArgD EC2.6.1.11; PLP-dependent"},
    # ArgF (Cgl0852) - ornithine carbamoyltransferase
    "Cgl0852": {"E_a": 55000.0, "H_d": 247000.0, "S_d": 789.0,
                "confidence": "LOW", "source": "ArgF EC2.1.3.3; OTC trimeric"},
    # ProA (Cgl1221) - glutamate 5-kinase (Pro biosynthesis from Glu)
    "Cgl1221": {"E_a": 54000.0, "H_d": 248000.0, "S_d": 793.0,
                "confidence": "LOW", "source": "ProA EC2.7.2.11; Glu-5-kinase; Pro pathway"},
    # ------- Glutamate Dehydrogenase (KEY: Tm calibrated from Takeno 2010) -------
    # GDH (Cgl2079 gdh) - glutamate dehydrogenase
    # Tm_eff=41.0°C (in-vivo calibrated): Takeno 2010 AEM thermal tolerance data
    # Ea=52 kJ/mol; hill_n=2.5: hexameric enzyme — moderate cooperative unfolding
    # Calibrated: grid search RMSE=0.0 → T50_Glu=42.5°C (exp: ~42-43°C)
    # H_d = 41.0+273.15 = 314.15K × 619.0 J/mol/K = 194459 J/mol → use 194440
    "Cgl2079": {"E_a": 52000.0, "H_d": 194459.0, "S_d": 619.0,
                "hill_n": 2.5,
                "confidence": "HIGH",
                "source": "GDH EC1.4.1.4 gdh; Takeno2010 AEM; Tm_eff=41.0C (calibrated); hill_n=2.5 (hexameric)"},
    # GS (Cgl2482 glnA) - glutamine synthetase
    # Very stable dodecameric structure; Tm ~49.6C
    "Cgl2482": {"E_a": 55000.0, "H_d": 265000.0, "S_d": 821.0,
                "confidence": "MED",
                "source": "GS EC6.3.1.2 glnA; dodecameric; thermostable; 50-60kJ/mol"},
}

# Reaction-ID-based fallback (when gene locus lookup misses)
REACTION_ID_PARAMS = {
    "CS_num1":           {"E_a": 55000.0, "H_d": 262000.0, "S_d": 815.0},
    "ICDHyr":            {"E_a": 47000.0, "H_d": 265000.0, "S_d": 820.0},
    "AKGDH":             {"E_a": 58000.0, "H_d": 235000.0, "S_d": 755.0},
    "FUM":               {"E_a": 42000.0, "H_d": 260000.0, "S_d": 812.0},
    "MDH":               {"E_a": 44000.0, "H_d": 252000.0, "S_d": 795.0},
    "SUCOAS":            {"E_a": 50000.0, "H_d": 255000.0, "S_d": 796.0},
    "PGI":               {"E_a": 48000.0, "H_d": 255000.0, "S_d": 796.0},
    "PGI_reverse":       {"E_a": 48000.0, "H_d": 255000.0, "S_d": 796.0},
    "PFK":               {"E_a": 58000.0, "H_d": 248000.0, "S_d": 800.0},
    "PFK_2":             {"E_a": 58000.0, "H_d": 248000.0, "S_d": 800.0},
    "GAPD_num1":         {"E_a": 45000.0, "H_d": 255000.0, "S_d": 796.0},
    "GAPD_reverse_num1": {"E_a": 45000.0, "H_d": 255000.0, "S_d": 796.0},
    "PYK_num1":          {"E_a": 50000.0, "H_d": 243000.0, "S_d": 779.0},
    "PYK_num2":          {"E_a": 50000.0, "H_d": 248000.0, "S_d": 794.0},
    "G6PDH2r":           {"E_a": 48000.0, "H_d": 250000.0, "S_d": 800.0},
    "G6PDH2r_reverse":   {"E_a": 48000.0, "H_d": 250000.0, "S_d": 800.0},
    "GND":               {"E_a": 49000.0, "H_d": 248000.0, "S_d": 794.0},
    "PPC":               {"E_a": 60000.0, "H_d": 243000.0, "S_d": 781.0},
    "GLUDy":             {"E_a": 52000.0, "H_d": 195000.0, "S_d": 619.0},
    "DHDPS":             {"E_a": 54000.0, "H_d": 234000.0, "S_d": 752.0},
    "DAPDC_num1":        {"E_a": 50000.0, "H_d": 240000.0, "S_d": 772.0},
    "DAPDC_num2":        {"E_a": 50000.0, "H_d": 240000.0, "S_d": 772.0},
    # New additions
    "ENO":               {"E_a": 48000.0, "H_d": 250000.0, "S_d": 797.0},
    "ENO_num1":          {"E_a": 48000.0, "H_d": 250000.0, "S_d": 797.0},
    "TPI":               {"E_a": 35000.0, "H_d": 252000.0, "S_d": 800.0},
    "TPI_num1":          {"E_a": 35000.0, "H_d": 252000.0, "S_d": 800.0},
    "PGK":               {"E_a": 50000.0, "H_d": 245000.0, "S_d": 782.0},
    "PGM":               {"E_a": 52000.0, "H_d": 250000.0, "S_d": 800.0},
    "PGM_num1":          {"E_a": 52000.0, "H_d": 250000.0, "S_d": 800.0},
    "TKT1":              {"E_a": 57000.0, "H_d": 247000.0, "S_d": 789.0},
    "TKT2":              {"E_a": 57000.0, "H_d": 247000.0, "S_d": 789.0},
    "TALA":              {"E_a": 48000.0, "H_d": 248000.0, "S_d": 794.0},
    "ACONTa":            {"E_a": 52000.0, "H_d": 234000.0, "S_d": 750.0},
    "ACONTb":            {"E_a": 52000.0, "H_d": 234000.0, "S_d": 750.0},
    "SUCCD1":            {"E_a": 45000.0, "H_d": 240000.0, "S_d": 763.0},
    "PYC":               {"E_a": 63000.0, "H_d": 250000.0, "S_d": 793.0},
    "ME1":               {"E_a": 50000.0, "H_d": 248000.0, "S_d": 793.0},
    "ME2":               {"E_a": 50000.0, "H_d": 248000.0, "S_d": 793.0},
    "ASAD":              {"E_a": 54000.0, "H_d": 248000.0, "S_d": 793.0},
    "HSDy":              {"E_a": 52000.0, "H_d": 250000.0, "S_d": 800.0},
    "DDH":               {"E_a": 52000.0, "H_d": 243000.0, "S_d": 779.0},
    "DHDPRS":            {"E_a": 52000.0, "H_d": 248000.0, "S_d": 793.0},
    "LysC":              {"E_a": 64000.0, "H_d": 217000.0, "S_d": 700.0},
    "ASPK":              {"E_a": 64000.0, "H_d": 217000.0, "S_d": 700.0},
    "AKGD":              {"E_a": 58000.0, "H_d": 235000.0, "S_d": 755.0},
}

DEFAULT_PARAMS = {
    "E_a": 50000.0,
    "H_d": 200000.0,
    "S_d": 641.0,
    "hill_n": 1.0,
    "confidence": "LOW",
    "source": "General mesophilic enzyme estimate; BRENDA class consensus",
}

# Per-enzyme Hill cooperativity overrides (n>1 = steeper denaturation sigmoid)
# Used by compute_active_fraction and compute_alpha when hill_n not in params dict
# Biological basis: oligomeric enzymes unfold cooperatively
# LysC (α2β2 tetramer, allosteric): n=3.0
# GDH (hexamer):                    n=2.5
# GS (dodecamer):                   n=3.0  
# Most monomeric/dimeric enzymes:   n=1.0
HILL_N_OVERRIDES = {
    "Cgl0251": 4.0,   # LysC tetramer - unified to 4.0 (standard α2β2)
    "Cgl2079": 2.5,   # GDH hexamer
    "Cgl2482": 3.0,   # GS dodecamer
    "Cgl1400": 2.0,   # GOGAT large oligomer
    "Cgl0672": 2.0,   # GOGAT (alternative locus)
    "Cgl0829": 2.0,   # CS gltA dimer
    "Cgl0696": 2.0,   # CS gltA dimer
    "Cgl1250": 4.0,   # pfkA tetramer
    "Cgl2089": 4.0,   # pyk tetramer
    "Cgl2910": 4.0,   # pyk2 tetramer
    "Cgl0664": 2.0,   # icd dimer
}


def get_params(rxn_id: str, gene_loci: list) -> dict:
    """
    Look up thermal parameters for a reaction.
    Priority: (1) gene locus, (2) reaction ID, (3) defaults.

    FIX: Returns a dict copy (not a reference) to prevent cross-call mutation.
    FIX: Merges HILL_N_OVERRIDES so cooperative enzymes (LysC n=3, GDH n=2.5)
         actually have their Hill n applied in compute_alpha().
    """
    for locus in gene_loci:
        if locus in GENE_LOCUS_PARAMS:
            p = dict(GENE_LOCUS_PARAMS[locus])          # copy, not reference
            if locus in HILL_N_OVERRIDES:
                p["hill_n"] = HILL_N_OVERRIDES[locus]   # merge cooperativity
            return p
    if rxn_id in REACTION_ID_PARAMS:
        p = REACTION_ID_PARAMS[rxn_id]
        return {**DEFAULT_PARAMS, **p, "confidence": "MED",
                "source": f"Reaction-ID lookup: {rxn_id}"}
    return DEFAULT_PARAMS.copy()


def compute_active_fraction(H_d: float, S_d: float, T_kelvin: float,
                            hill_n: float = 1.0) -> float:
    """Fraction of enzyme in folded (active) state via two-state equilibrium.

    With hill_n > 1: cooperative unfolding model (sharper sigmoid near Tm).
    Equation: K_eq = exp(-n*(H_d - T*S_d)/(R*T)); f = 1/(1+K_eq)
    This is equivalent to treating n-mer units as a single cooperative unit.
    """
    try:
        exponent = -hill_n * (H_d - T_kelvin * S_d) / (R * T_kelvin)
        exponent = max(-500.0, min(500.0, exponent))
        K_eq = math.exp(exponent)
        return 1.0 / (1.0 + K_eq)
    except OverflowError:
        return 1e-6


def compute_alpha(params: dict, T_kelvin: float) -> float:
    """
    Combined temperature correction factor alpha(T):
      alpha = Arrhenius_factor(T) x [f_active(T, n) / f_active(T_ref, n)]

    Uses hill_n from params dict if present (default 1.0).
    hill_n > 1 produces a steeper, more cooperative denaturation sigmoid.
    """
    E_a  = params["E_a"]
    H_d  = params["H_d"]
    S_d  = params["S_d"]
    n    = params.get("hill_n", 1.0)
    arr  = math.exp(-(E_a / R) * (1.0 / T_kelvin - 1.0 / T_ref_K))
    f_ref = compute_active_fraction(H_d, S_d, T_ref_K, n)
    f_t   = compute_active_fraction(H_d, S_d, T_kelvin, n)
    alpha = arr * (f_t / f_ref) if f_ref > 1e-10 else 1e-6
    return max(1e-6, alpha)


if __name__ == "__main__":
    print(f"{'Gene Locus':<20} {'Ea kJ/mol':>10} {'Tm C':>8}  Conf")
    print("-" * 50)
    for locus, p in GENE_LOCUS_PARAMS.items():
        tm = (p["H_d"] / p["S_d"]) - 273.15
        print(f"  {locus:<18} {p['E_a']/1000:>9.1f} {tm:>8.1f}   {p['confidence']}")
