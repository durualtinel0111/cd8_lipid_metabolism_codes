"""
Causal Inference Analysis between Lipid Accumulation Scores and Exhaustion Markers Using DoWhy
"""

# Required: dowhy package
# !pip install dowhy --quiet

import pandas as pd
import numpy as np
from dowhy import CausalModel

# 1. Load data
lipid_accum = pd.read_excel("/content/lipid_birikim_zscore.xlsx", index_col=0)
exhaustion = pd.read_excel("/content/merged_genes_log2cpm.xlsx", index_col=0)

# 2. Transpose lipid accumulation scores to have samples as columns
lipid_accum = lipid_accum.T

# 3. Normalize sample names
def normalize_name(name):
    return str(name).strip().upper().replace("-", "").replace("_", "").replace(".", "")

lipid_accum.columns = [normalize_name(x) for x in lipid_accum.columns]
exhaustion.columns = [normalize_name(x) for x in exhaustion.columns]

# 4. Find common samples
common_samples = list(set(lipid_accum.columns) & set(exhaustion.columns))
print(f"Number of common samples: {len(common_samples)}")

lipid_accum = lipid_accum[common_samples]
exhaustion = exhaustion[common_samples]

# 5. Select exhaustion marker genes
exhaustion_genes = ["TOX", "TIGIT", "PRDM1", "PDCD1", "NR4A1", "LAG3", "HAVCR2", "EOMES", "ENTPD1", "CTLA4", "BATF"]
common_genes = [g for g in exhaustion_genes if g in exhaustion.index]
if not common_genes:
    raise ValueError("No exhaustion marker genes found in the dataset!")

exhaustion_expr = exhaustion.loc[common_genes]

# 6. Compute exhaustion score as mean expression of marker genes
exhaustion_score = exhaustion_expr.mean(axis=0)

# 7. Get lipid accumulation scores (samples aligned)
lipid_scores = lipid_accum.iloc[0] if lipid_accum.shape[0] == 1 else lipid_accum.mean(axis=0)

# 8. Align samples between scores
common_samples_final = lipid_scores.index.intersection(exhaustion_score.index)
lipid_scores_aligned = lipid_scores.loc[common_samples_final]
exhaustion_scores_aligned = exhaustion_score.loc[common_samples_final]

# 9. Generate synthetic confounders (age, sex)
np.random.seed(42)
age = np.random.randint(30, 70, size=len(common_samples_final))
sex = np.random.randint(0, 2, size=len(common_samples_final))

# 10. Prepare dataframe for causal analysis
df = pd.DataFrame({
    "lipid_accum_zscore": lipid_scores_aligned,
    "exhaustion": exhaustion_scores_aligned,
    "age": age,
    "sex": sex
})

# 11. Define and identify causal model
model = CausalModel(
    data=df,
    treatment="lipid_accum_zscore",
    outcome="exhaustion",
    common_causes=["age", "sex"]
)
identified_model = model.identify_effect()

# 12. Estimate Average Treatment Effect (ATE)
estimate = model.estimate_effect(identified_model, method_name="backdoor.linear_regression")
print("ATE (lipid accumulation → exhaustion):", estimate.value)

# 13. Reverse causality check (optional)
model_reverse = CausalModel(
    data=df,
    treatment="exhaustion",
    outcome="lipid_accum_zscore",
    common_causes=["age", "sex"]
)
identified_model_reverse = model_reverse.identify_effect()
estimate_reverse = model_reverse.estimate_effect(identified_model_reverse, method_name="backdoor.linear_regression")
print("ATE (exhaustion → lipid accumulation):", estimate_reverse.value)

# Save exhaustion scores
exhaustion_scores_aligned.to_frame(name='EXHAUSTION_SCORE').to_csv("/content/exhaustion_scores_lipid_birikimi.csv")
print("\nExhaustion scores saved to '/content/exhaustion_scores_lipid_birikimi.csv'.")
