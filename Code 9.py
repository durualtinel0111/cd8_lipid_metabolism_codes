"""
Causal Inference Analysis between Lipid Degradation Scores and Exhaustion Markers Using DoWhy
"""

# Required: dowhy package
# !pip install dowhy --quiet

import pandas as pd
import numpy as np
from dowhy import CausalModel

def normalize_name(name):
    return str(name).strip().upper().replace("-", "").replace("_", "").replace(".", "")

# 1. Load lipid catabolic scores
lipid_catabolic_df = pd.read_excel("/content/lipid_catabolic_scores.xlsx", index_col=0)
if lipid_catabolic_df.shape[0] == 1 and lipid_catabolic_df.shape[1] > 1:
    lipid_catabolic_scores = lipid_catabolic_df.T.iloc[:,0]
elif lipid_catabolic_df.shape[1] == 1:
    lipid_catabolic_scores = lipid_catabolic_df.iloc[:,0]
else:
    lipid_catabolic_scores = lipid_catabolic_df.mean(axis=1)

# 2. Load gene expression for exhaustion score calculation
merged_expr = pd.read_excel("/content/merged_genes_log2cpm.xlsx", index_col=0)

# 3. Exhaustion marker genes
exhaustion_genes = ["TOX", "TIGIT", "PRDM1", "PDCD1", "NR4A1", "LAG3", "HAVCR2", "EOMES", "ENTPD1", "CTLA4", "BATF"]
common_genes = [g for g in exhaustion_genes if g in merged_expr.index]
if not common_genes:
    raise ValueError("No exhaustion marker genes found in the dataset!")

exhaustion_expr = merged_expr.loc[common_genes]

# 4. Compute exhaustion score
exhaustion_score = exhaustion_expr.mean(axis=0)

# 5. Normalize sample names
lipid_catabolic_scores.index = [normalize_name(x) for x in lipid_catabolic_scores.index]
exhaustion_score.index = [normalize_name(x) for x in exhaustion_score.index]

# 6. Identify common samples
common_samples = lipid_catabolic_scores.index.intersection(exhaustion_score.index)
print(f"Number of common samples: {len(common_samples)}")

lipid_scores_aligned = lipid_catabolic_scores.loc[common_samples]
exhaustion_scores_aligned = exhaustion_score.loc[common_samples]

# 7. Synthetic confounders
np.random.seed(42)
age = np.random.randint(30, 70, size=len(common_samples))
sex = np.random.randint(0, 2, size=len(common_samples))

# 8. Dataframe for causal analysis
df = pd.DataFrame({
    "lipid_catabolic_zscore": lipid_scores_aligned,
    "exhaustion": exhaustion_scores_aligned,
    "age": age,
    "sex": sex
})

# 9. Define and identify causal model
model = CausalModel(
    data=df,
    treatment="lipid_catabolic_zscore",
    outcome="exhaustion",
    common_causes=["age", "sex"]
)
identified_model = model.identify_effect()

# 10. Estimate ATE
estimate = model.estimate_effect(identified_model, method_name="backdoor.linear_regression")
print("ATE (lipid catabolic → exhaustion):", estimate.value)

# 11. Reverse causality test
model_reverse = CausalModel(
    data=df,
    treatment="exhaustion",
    outcome="lipid_catabolic_zscore",
    common_causes=["age", "sex"]
)
identified_model_reverse = model_reverse.identify_effect()
estimate_reverse = model_reverse.estimate_effect(identified_model_reverse, method_name="backdoor.linear_regression")
print("ATE (exhaustion → lipid catabolic):", estimate_reverse.value)
