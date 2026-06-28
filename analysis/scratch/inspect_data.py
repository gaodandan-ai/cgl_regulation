import pandas as pd
import os
import cobra

# Paths
DATA_DIR = r"F:\cgl_regulation\data"
RNA_SEQ_DIR = os.path.join(DATA_DIR, "rna_seq")
MODEL_PATH = r"F:\cgl_regulation\backend\models\iCW773.xml"

print("--- Checking Gene Mapping ---")
mapping_df = pd.read_csv(os.path.join(DATA_DIR, "gene_mapping.csv"))
print(f"Mapping rows: {len(mapping_df)}")
print(mapping_df.head(5))

# Create mapping dictionary
cgl_to_cg = {}
cg_to_cgl = {}
for _, row in mapping_df.dropna(subset=['cgl_locus']).iterrows():
    cgl = row['cgl_locus'].strip()
    cg = row['cg_locus']
    if pd.notna(cg):
        cg = str(cg).strip()
        cgl_to_cg[cgl] = cg
        cg_to_cgl[cg] = cgl

print(f"Mapped cgl_to_cg entries: {len(cgl_to_cg)}")
print(f"Mapped cg_to_cgl entries: {len(cg_to_cgl)}")

print("\n--- Checking Regulations ---")
regs_df = pd.read_csv(os.path.join(DATA_DIR, "regulations.csv"))
print(f"Total regulatory relations: {len(regs_df)}")
unique_tfs_cg = regs_df['TF_locusTag'].dropna().unique()
print(f"Unique TFs in regulations (cg_locus): {len(unique_tfs_cg)}")
print("Sample TFs:", unique_tfs_cg[:10])

# Map TFs to Cgl
tfs_mapped = [cg_to_cgl[tf] for tf in unique_tfs_cg if tf in cg_to_cgl]
print(f"TFs successfully mapped to Cgl: {len(tfs_mapped)} / {len(unique_tfs_cg)}")

print("\n--- Checking Expression Files ---")
exp_1h = pd.read_csv(os.path.join(RNA_SEQ_DIR, "03_1h_normalized_expression.csv"))
print(f"1h expression rows: {len(exp_1h)}")
print("Columns:", exp_1h.columns.tolist())
print(exp_1h.head(3))

# Check how many TFs are in 1h expression
tf_in_exp = [tf for tf in tfs_mapped if tf in exp_1h['Geneid'].values]
print(f"TFs in 1h expression: {len(tf_in_exp)} / {len(tfs_mapped)}")

print("\n--- Checking DEG Files ---")
deg_1h = pd.read_csv(os.path.join(RNA_SEQ_DIR, "02_1h_significant_DEGs.csv"))
print(f"1h DEGs count: {len(deg_1h)}")
print("Columns:", deg_1h.columns.tolist())
print(deg_1h.head(3))

print("\n--- Checking Metabolic Model ---")
if os.path.exists(MODEL_PATH):
    try:
        model = cobra.io.read_sbml_model(MODEL_PATH)
        print(f"Model ID: {model.id}")
        print(f"Number of reactions: {len(model.reactions)}")
        print(f"Number of genes: {len(model.genes)}")
        print(f"Number of metabolites: {len(model.metabolites)}")
        sample_genes = [g.id for g in list(model.genes)[:5]]
        print("Sample genes in model:", sample_genes)
        # Check matching to Cgl or cg
        matches_cgl = sum(1 for g in model.genes if g.id in cgl_to_cg or g.id.replace("g_", "") in cgl_to_cg)
        matches_cg = sum(1 for g in model.genes if g.id in cg_to_cgl or g.id.replace("g_", "") in cg_to_cgl)
        print(f"Model genes matching cgl: {matches_cgl}")
        print(f"Model genes matching cg: {matches_cg}")
    except Exception as e:
        print("Error loading model:", e)
else:
    print("Model file not found at", MODEL_PATH)
