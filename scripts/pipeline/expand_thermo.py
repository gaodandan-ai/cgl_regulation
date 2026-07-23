"""scripts/expand_thermo.py — appends new pathway blocks to build_thermo_data.py CURATED_DGR0"""
import re

EXTRA = {
    # tRNA SYNTHETASES (complete set)
    "HISTRS":        (-25.5, "H", "[K] Histidyl-tRNA synthetase"),
    "ILETRS":        (-25.5, "H", "[K] Isoleucyl-tRNA synthetase"),
    "LEUTRS":        (-25.5, "H", "[K] Leucyl-tRNA synthetase"),
    "LYSTRS":        (-25.5, "H", "[K] Lysyl-tRNA synthetase"),
    "METTRS":        (-25.5, "H", "[K] Methionyl-tRNA synthetase"),
    # NAD/NADP METABOLISM
    "NADK":          (-18.5, "M", "[K] NAD kinase (ATP-driven)"),
    "NADS1":         (-22.0, "M", "[K] NAD synthase (Gln-hydrolyzing)"),
    "NADH5":         (-52.0, "H", "[T77] NADH dehydrogenase (ubiquinone-8)"),
    "NADH9":         (-52.0, "H", "[T77] NADH dehydrogenase (demethylmenaquinone)"),
    "NADPHQR2":      (-52.0, "H", "[T77] NADPH quinone reductase (ubiquinone)"),
    "NADPHQR3":      (-52.0, "H", "[T77] NADPH quinone reductase (menaquinone)"),
    "NADPHQR4":      (-52.0, "H", "[T77] NADPH quinone reductase (demethylmenaquinone)"),
    "NDPK2":         (  0.0, "H", "[F12] NDP kinase (ATP:UDP); near-eq"),
    "NDPK3":         (  0.0, "H", "[F12] NDP kinase (ATP:CDP)"),
    "NDPK4":         (  0.0, "H", "[F12] NDP kinase (ATP:dTDP)"),
    "NDPK5":         (  0.0, "H", "[F12] NDP kinase (ATP:dGDP)"),
    "NDPK6":         (  0.0, "H", "[F12] NDP kinase (ATP:dUDP)"),
    # LEU/ILE BIOSYNTHESIS
    "IPMD":          (-12.0, "M", "[K] 3-isopropylmalate dehydrogenase"),
    "IPPMIa":        ( -3.5, "M", "[K] 3-isopropylmalate dehydratase step a"),
    "IPPMIb":        ( -3.5, "M", "[K] 3-isopropylmalate dehydratase step b"),
    "IPPS":          (-16.0, "M", "[K] 2-isopropylmalate synthase"),
    # MEP/ISOPRENOID
    "MECDPS":        (-14.5, "M", "[K] 2C-methyl-D-erythritol-2,4-CDP synthase"),
    "IPDDI":         (  2.5, "M", "[K] Isopentenyl-PP isomerase; near-eq"),
    "IPDPS":         (-28.5, "M", "[K] HMBPP reductase"),
    "MEPCT":         (-22.0, "M", "[K] MEP cytidylyltransferase"),
    # HISTIDINE BIOSYNTHESIS
    "IG3PS":         (-22.0, "M", "[K] Imidazole-glycerol-3-P synthase"),
    "IGPDH":         (-14.5, "M", "[K] Imidazoleglycerol-P dehydratase"),
    "IMPC":          (-18.5, "M", "[K] IMP cyclohydrolase"),
    "IMPD":          (-22.0, "M", "[K] IMP dehydrogenase (NAD+)"),
    "HSTP":          (-18.5, "M", "[K] Histidinol-phosphate transaminase"),
    # HEXOSE METABOLISM
    "HEX1":          (-16.7, "H", "[F12] Hexokinase (glucose + ATP)"),
    "HEX7":          (-16.7, "H", "[F12] Hexokinase (fructose + ATP)"),
    "MAN6PI":        (  1.8, "M", "[K] Mannose-6-P isomerase; near-eq"),
    "M1PD":          (-14.5, "M", "[K] Mannitol-1-P dehydrogenase"),
    "MAN2D":         (-14.5, "M", "[K] Mannitol-2-dehydrogenase"),
    # METHIONINE / SAM CYCLE
    "METAT":         (-38.5, "H", "[K] Methionine adenosyltransferase (SAM synthesis)"),
    "MTHFR2":        (-22.0, "M", "[K] 5,10-methylene-THF reductase (NADH)"),
    "HSK":           (-18.5, "M", "[K] Homoserine kinase (ATP-driven)"),
    "HSST":          (-22.5, "M", "[K] Homoserine O-acetyltransferase"),
    "HSST_2":        (-22.5, "M", "[K] HSST isoform 2"),
    "CYSTL":         (-14.5, "M", "[K] Cystathionine beta-lyase"),
    "CYSS":          (-22.0, "M", "[K] Cysteine synthase (serine + H2S)"),
    # PROPIONATE / METHYLCITRATE
    "MMM2":          (  0.3, "M", "[K] Methylmalonyl-CoA mutase; near-eq"),
    "MCITD":         ( -3.5, "M", "[K] 2-methylcitrate dehydratase"),
    "MCITL2":        ( -5.0, "M", "[K] Methylisocitrate lyase"),
    "MSDH":          (-28.5, "M", "[K] Methylmalonate-semialdehyde dehydrogenase"),
    "MSDHD":         (-28.5, "M", "[K] Malonate-semialdehyde dehydrogenase"),
    "MM_COA_ADD5":   (-22.5, "M", "[K] Propanoyl-CoA carboxylase (ATP-driven)"),
    # LACTATE
    "LDH":           (-25.0, "M", "[K] L-Lactate dehydrogenase (NAD+)"),
    "LDH_D":         (-25.0, "M", "[K] D-Lactate dehydrogenase"),
    "LDH_D2":        (-25.0, "M", "[K] D-LDH isoform 2"),
    "L_LACD2":       (-22.5, "M", "[K] L-Lactate dehydrogenase (ubiquinone)"),
    "L_LACD3":       (-22.5, "M", "[K] L-Lactate dehydrogenase (menaquinone)"),
    # PORPHYRIN / HEME
    "HMBS":          (-22.5, "M", "[K] Hydroxymethylbilane synthase"),
    "FCLT":          (-18.5, "M", "[K] Ferrochelatase (heme insertion)"),
    # FOLATE / CHORISMATE
    "ICHOR":         ( -8.5, "M", "[K] Isochorismate synthase"),
    "HPPK2":         (-18.5, "M", "[K] 6-Hydroxymethyl-dihydropterin pyrophosphokinase"),
    "FMETTRS":       (-25.5, "H", "[K] Met-tRNA formyltransferase"),
    # LIPOIC ACID / BIOTIN (radical-SAM)
    "LIPOCT":        (-22.0, "M", "[K] Lipoyl(octanoyl) transferase"),
    "LIPOS":         (-38.5, "H", "[K] Lipoate synthase (radical SAM)"),
    "LIPAMPL":       (-22.5, "M", "[K] Lipoyl-adenylate protein ligase"),
    "BTS4":          (-42.5, "H", "[K] Biotin synthase (radical SAM)"),
    "BTNC":          (-22.0, "M", "[K] Biotin carboxylase (ACC subunit)"),
    # NUCLEOTIDE SALVAGE
    "HXPRT":         ( -5.5, "M", "[K] Hypoxanthine-PRPP phosphoribosyltransferase"),
    "INSH":          (-15.0, "M", "[K] Inosine hydrolase"),
    # FAS INITIATION
    "ACACT1r":       (-12.5, "M", "[K] Acetyl-CoA acetyltransferase (thiolase)"),
    "MCOATA":        (-22.5, "M", "[K] Malonyl-CoA:ACP acyltransferase"),
    "MACPD":         (-18.5, "H", "[K] Malonyl-ACP decarboxylase (irreversible)"),
    "KAS14":         (-30.5, "M", "[K] Beta-ketoacyl-ACP synthase (FabF)"),
    "KAS15":         (-30.5, "M", "[K] Beta-ketoacyl-ACP synthase isoform 2"),
    # MISC
    "HCO3E":         (  0.0, "H", "[F12] HCO3 equilibration; dG approx 0"),
    "MALGT":         (-22.5, "M", "[K] Maltose glucosyltransferase"),
    "METAT":         (-38.5, "H", "[K] Methionine adenosyltransferase"),
}

# Build the lines to inject
lines_to_add = "\n    # ━━━━━━ v3 EXPANSION: Additional pathway blocks ━━━━━━\n"
for rxn_id, (dgr, conf, src) in EXTRA.items():
    lines_to_add += f'    "{rxn_id}":{" " * max(1, 15 - len(rxn_id))}({dgr:5.1f}, "{conf}", "{src}"),\n'

# Read existing file
with open("scripts/build_thermo_data.py", "r", encoding="utf-8") as f:
    src_text = f.read()

# The dict closes with:  \n    "PPCK": ... \n}\n
# Insert new entries before the lone closing }
CLOSE_MARKER = '    "PPCK":          ( -7.1, "M", "[F12] PEP carboxykinase"),          # duplicate OK\n}'
if CLOSE_MARKER not in src_text:
    print("ERROR: marker not found! Check the file.")
else:
    new_text = src_text.replace(CLOSE_MARKER, CLOSE_MARKER.rstrip("}") + lines_to_add + "}", 1)
    with open("scripts/build_thermo_data.py", "w", encoding="utf-8") as f:
        f.write(new_text)
    print(f"Inserted {len(EXTRA)} new entries. Total file size: {len(new_text)} bytes")
