"""
Calculation of Lipid Accumulation Scores from GSVA Results
"""

import pandas as pd

# Step 1: Load GSVA scores (pathways x samples)
gsva_scores = pd.read_excel("gsva_results.xlsx", index_col=0)  # Dosya yolunu uygun şekilde değiştir

# Step 2: Define lipid-related pathways and their direction (1 = accumulation, -1 = degradation/export)
lipid_pathway_weights = {
    "GOBP_LIPID_DROPLET_FORMATION": 1,
    "GOBP_LIPID_DROPLET_DISASSEMBLY": -1,
    "GOBP_LIPID_EXPORT_FROM_CELL": -1,
    "GOBP_LIPID_HOMEOSTASIS": 1,
    "GOBP_LIPID_IMPORT_INTO_CELL": 1,
    "GOBP_LIPID_LOCALIZATION": 1,
    "GOBP_LIPID_METABOLIC_PROCESS": 1,
    "GOBP_LIPID_CATABOLIC_PROCESS": -1,
    "GOBP_LIPID_BIOSYNTHETIC_PROCESS": 1,
    "GOBP_LIPID_DROPLET_FUSION": 1,
    "GOBP_LIPOPROTEIN_METABOLIC_PROCESS": 1,
    "GOBP_NEGATIVE_REGULATION_OF_LIPID_STORAGE": -1,
    "GOBP_POSITIVE_REGULATION_OF_LIPID_STORAGE": 1,
    "GOBP_POSITIVE_REGULATION_OF_LIPID_CATABOLIC_PROCESS": -1,
    "GOBP_POSITIVE_REGULATION_OF_MACROPHAGE_DERIVED_FOAM_CELL_DIFFERENTIATION": 1,
    "GOBP_LIPOPHAGY": -1
}

# Step 3: Filter pathways present in GSVA results
valid_pathways = [p for p in lipid_pathway_weights if p in gsva_scores.index]
weights = pd.Series(lipid_pathway_weights).loc[valid_pathways]
selected_scores = gsva_scores.loc[valid_pathways]

# Step 4: Compute lipid accumulation score per sample
lipid_accumulation_scores = (selected_scores.T * weights).T.sum()

# Step 5: Standardize scores (Z-score)
z_scores = (lipid_accumulation_scores - lipid_accumulation_scores.mean()) / lipid_accumulation_scores.std()

# Output
print("Lipid Accumulation Scores:")
print(lipid_accumulation_scores)
print("\nZ-scores of Lipid Accumulation:")
print(z_scores)
