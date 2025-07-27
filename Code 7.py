"""
Calculation of Lipid Degradation (Catabolic) Scores from GSVA Results
"""

import pandas as pd

# 1. Load GSVA scores
gsva_scores = pd.read_excel("/content/gsva_results.xlsx", index_col=0)

# 2. Define lipid catabolism-related pathways and weights (-1 means promote catabolism)
lipid_catabolic_pathways = {
    "GOBP_LIPID_CATABOLIC_PROCESS": -1,
    "GOBP_NEGATIVE_REGULATION_OF_LIPID_BIOSYNTHETIC_PROCESS": -1,
    "GOBP_POSITIVE_REGULATION_OF_LIPID_CATABOLIC_PROCESS": -1,
    "GOBP_LIPOPHAGY": -1,
    "REACTOME_GLYCEROPHOSPHOLIPID_CATABOLISM": -1,
    "REACTOME_SPHINGOLIPID_CATABOLISM": -1,
    "REACTOME_PHOSPHOLIPID_METABOLISM": -1,
    "GOBP_LIPID_DROPLET_DISASSEMBLY": -1,
    "GOBP_NEGATIVE_REGULATION_OF_LIPID_STORAGE": -1,
    "WP_DEGRADATION_PATHWAY_OF_SPHINGOLIPIDS_INCLUDING_DISEASES": -1,
    "GOBP_NEUTRAL_LIPID_CATABOLIC_PROCESS": -1,
    "GOBP_PHOSPHOLIPID_CATABOLIC_PROCESS": -1,
    "GOBP_LIPOPROTEIN_CATABOLIC_PROCESS": -1,
    "REACTOME_GLYCOSPHINGOLIPID_CATABOLISM": -1
}

# 3. Filter valid pathways in GSVA results
valid_pathways = [p for p in lipid_catabolic_pathways if p in gsva_scores.index]
weights = pd.Series(lipid_catabolic_pathways).loc[valid_pathways]
selected_scores = gsva_scores.loc[valid_pathways]

# 4. Calculate lipid catabolic score per sample
lipid_catabolic_scores = (selected_scores.T * weights).T.sum()

# 5. Standardize scores (Z-score)
z_scores = (lipid_catabolic_scores - lipid_catabolic_scores.mean()) / lipid_catabolic_scores.std()

# 6. Save results as Excel
output_df = pd.DataFrame({
    "Lipid_Catabolic_Score": lipid_catabolic_scores,
    "Z_Score": z_scores
})
output_df.to_excel("/content/lipid_catabolic_scores.xlsx")

print("Lipid catabolic scores saved to '/content/lipid_catabolic_scores.xlsx'.")
