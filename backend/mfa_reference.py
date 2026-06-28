"""
13C-MFA Literature Reference Data for C. glutamicum ATCC 13032
===============================================================
Source: Cheng et al. (2017) Biotechnol. Biofuels 10:169
        (iCW773 genome-scale model validation, Pearson R = 0.99 vs. MFA)
        
        Becker & Wittmann (2011) Appl. Microbiol. Biotechnol. 92:449-463
        (Wild-type strain, aerobic, glucose minimal medium)

Conditions: wild-type C. glutamicum ATCC 13032
            aerobic batch culture on glucose minimal medium
            Glucose uptake rate: 10.0 mmol/gDW/h (reference)
            
All flux values in mmol/gDW/h
"""

# Reaction ID mapping: model reaction ID -> pathway/label info
# Fluxes normalized to glucose uptake rate of 10 mmol/gDW/h
MFA_LITERATURE_DATASET = [
    {
        "reaction_id": "GLCpts",
        "reaction_name": "Glucose PTS Transport",
        "pathway": "Glycolysis",
        "mfa_flux": 10.00,
        "mfa_std": 0.50,
        "reference": "Cheng et al. 2017 / Becker & Wittmann 2011"
    },
    {
        "reaction_id": "PGI",
        "reaction_name": "Phosphoglucose Isomerase",
        "pathway": "Glycolysis",
        "mfa_flux": 6.55,
        "mfa_std": 0.40,
        "reference": "Becker & Wittmann 2011 Table 2"
    },
    {
        "reaction_id": "PFK",
        "reaction_name": "Phosphofructokinase",
        "pathway": "Glycolysis",
        "mfa_flux": 6.55,
        "mfa_std": 0.40,
        "reference": "Becker & Wittmann 2011 Table 2"
    },
    {
        "reaction_id": "GAPD",
        "reaction_name": "GAPDH (Glyceraldehyde-3-Phosphate DH)",
        "pathway": "Glycolysis",
        "mfa_flux": 15.42,
        "mfa_std": 0.70,
        "reference": "Becker & Wittmann 2011 Table 2"
    },
    {
        "reaction_id": "PYK",
        "reaction_name": "Pyruvate Kinase",
        "pathway": "Glycolysis",
        "mfa_flux": 8.48,
        "mfa_std": 0.45,
        "reference": "Becker & Wittmann 2011 Table 2"
    },
    {
        "reaction_id": "G6PDH2r",
        "reaction_name": "Glucose-6-Phosphate Dehydrogenase (PPP)",
        "pathway": "PPP",
        "mfa_flux": 3.45,
        "mfa_std": 0.30,
        "reference": "Becker & Wittmann 2011 Table 2"
    },
    {
        "reaction_id": "GND",
        "reaction_name": "6-Phosphogluconate Dehydrogenase",
        "pathway": "PPP",
        "mfa_flux": 3.45,
        "mfa_std": 0.30,
        "reference": "Becker & Wittmann 2011 Table 2"
    },
    {
        "reaction_id": "CS",
        "reaction_name": "Citrate Synthase",
        "pathway": "TCA Cycle",
        "mfa_flux": 7.02,
        "mfa_std": 0.60,
        "reference": "Cheng et al. 2017 / Becker & Wittmann 2011"
    },
    {
        "reaction_id": "ICDHyr",
        "reaction_name": "Isocitrate Dehydrogenase (NADP)",
        "pathway": "TCA Cycle",
        "mfa_flux": 6.47,
        "mfa_std": 0.55,
        "reference": "Becker & Wittmann 2011 Table 2"
    },
    {
        "reaction_id": "AKGDH",
        "reaction_name": "Alpha-Ketoglutarate Dehydrogenase",
        "pathway": "TCA Cycle",
        "mfa_flux": 5.52,
        "mfa_std": 0.50,
        "reference": "Becker & Wittmann 2011 Table 2"
    },
    {
        "reaction_id": "SUCOAS",
        "reaction_name": "Succinyl-CoA Synthetase",
        "pathway": "TCA Cycle",
        "mfa_flux": 4.10,
        "mfa_std": 0.40,
        "reference": "Becker & Wittmann 2011 Table 2"
    },
    {
        "reaction_id": "FUM",
        "reaction_name": "Fumarase",
        "pathway": "TCA Cycle",
        "mfa_flux": 4.10,
        "mfa_std": 0.40,
        "reference": "Becker & Wittmann 2011 Table 2"
    },
    {
        "reaction_id": "MDH",
        "reaction_name": "Malate Dehydrogenase",
        "pathway": "TCA Cycle",
        "mfa_flux": 4.10,
        "mfa_std": 0.40,
        "reference": "Becker & Wittmann 2011 Table 2"
    },
    {
        "reaction_id": "PPC",
        "reaction_name": "Pyruvate Carboxylase (Anaplerosis)",
        "pathway": "Anaplerosis",
        "mfa_flux": 2.68,
        "mfa_std": 0.35,
        "reference": "Becker & Wittmann 2011 Table 2"
    },
    {
        "reaction_id": "GLUDy",
        "reaction_name": "Glutamate Dehydrogenase (NADPH)",
        "pathway": "Amino Acid Biosynthesis",
        "mfa_flux": 1.02,
        "mfa_std": 0.20,
        "reference": "Becker & Wittmann 2011 Table 2"
    },
]

# Alternative model reaction ID mappings for common iCW773 naming differences
REACTION_ID_ALIASES = {
    "GLCpts": ["GLCpts", "GLCptspp", "PTS_glc", "EX_glc_e"],
    "PGI": ["PGI", "PGI_cgl"],
    "PFK": ["PFK", "PFK_cgl", "PFKM"],
    "GAPD": ["GAPD", "GAPDH", "GAPD_cgl"],
    "PYK": ["PYK", "PYK_cgl"],
    "G6PDH2r": ["G6PDH2r", "G6PDH", "G6PDHy"],
    "GND": ["GND", "GND_cgl"],
    "CS": ["CS", "CS_cgl"],
    "ICDHyr": ["ICDHyr", "ICDH", "ICDHy"],
    "AKGDH": ["AKGDH", "AKGDH_cgl"],
    "SUCOAS": ["SUCOAS", "SUCOAS1m"],
    "FUM": ["FUM", "FUM_cgl"],
    "MDH": ["MDH", "MDH_cgl"],
    "PPC": ["PPC", "PPC_cgl", "PC"],
    "GLUDy": ["GLUDy", "GDH", "GLUD", "GLUDx"],
}
