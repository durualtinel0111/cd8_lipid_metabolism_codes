"""
Identification of Top Genes Driving Lipid Accumulation and Degradation Using LightGBM

Trains LightGBM regression models on gene expression and lipid scores,
computes feature importance averaged over multiple runs.
"""

import pandas as pd
import numpy as np
import lightgbm as lgb
from sklearn.preprocessing import StandardScaler

# 1. Load gene expression data
gen_expr_path = "/content/1213_gen.xlsx"
gen_df = pd.read_excel(gen_expr_path, index_col=0)
gen_df = gen_df.drop(index="grup", errors='ignore')  # 'grup' satırı varsa kaldır
gen_df.columns = [str(col) for col in gen_df.columns]
gen_df = gen_df.T  # Satırlar: örnekler, sütunlar: genler

# 2. Load lipid accumulation and degradation scores
lipid_accum_df = pd.read_excel("/content/lipid_birikim_zscore.xlsx", index_col=0)
lipid_degrad_df = pd.read_excel("/content/lipid_catabolic_scores.xlsx", index_col=0)

# 3. Select common samples
common_samples_accum = gen_df.index.intersection(lipid_accum_df.index)
common_samples_degrad = gen_df.index.intersection(lipid_degrad_df.index)

X_accum = gen_df.loc[common_samples_accum]
y_accum = lipid_accum_df.loc[common_samples_accum].iloc[:, 0]

X_degrad = gen_df.loc[common_samples_degrad]
y_degrad = lipid_degrad_df.loc[common_samples_degrad].iloc[:, 0]

# 4. Standardize features
scaler = StandardScaler()
X_accum_scaled = pd.DataFrame(scaler.fit_transform(X_accum), index=X_accum.index, columns=X_accum.columns)
X_degrad_scaled = pd.DataFrame(scaler.fit_transform(X_degrad), index=X_degrad.index, columns=X_degrad.columns)

# 5. LightGBM hyperparameters (fixed)
params_accum = {
    'learning_rate': 0.1344,
    'num_leaves': 121,
    'max_depth': 30,
    'subsample': 0.828,
    'colsample_bytree': 0.855,
    'reg_alpha': 0.0473,
    'reg_lambda': 3.667,
    'n_estimators': 500,
    'objective': 'regression',
    'metric': 'rmse'
}

params_degrad = {
    'learning_rate': 0.1916,
    'num_leaves': 51,
    'max_depth': 5,
    'subsample': 0.678,
    'colsample_bytree': 0.787,
    'reg_alpha': 1.869,
    'reg_lambda': 2.354,
    'n_estimators': 500,
    'objective': 'regression',
    'metric': 'rmse'
}

# 6. Train models N times to average feature importance
N = 30
random_state = 42

importance_accum = pd.Series(0, index=X_accum.columns, dtype=float)
importance_degrad = pd.Series(0, index=X_degrad.columns, dtype=float)

for i in range(N):
    model_accum = lgb.LGBMRegressor(**params_accum, random_state=random_state + i)
    model_accum.fit(X_accum_scaled, y_accum)
    importance_accum += pd.Series(model_accum.feature_importances_, index=X_accum.columns)

    model_degrad = lgb.LGBMRegressor(**params_degrad, random_state=random_state + i)
    model_degrad.fit(X_degrad_scaled, y_degrad)
    importance_degrad += pd.Series(model_degrad.feature_importances_, index=X_degrad.columns)

importance_accum /= N
importance_degrad /= N

# 7. Select top 50 genes by importance
top_50_accum = importance_accum.sort_values(ascending=False).head(50)
top_50_degrad = importance_degrad.sort_values(ascending=False).head(50)

# 8. Save results to Excel
top_50_accum.to_frame(name="Importance").to_excel("/content/top_50_genes_lipid_accumulation.xlsx")
top_50_degrad.to_frame(name="Importance").to_excel("/content/top_50_genes_lipid_degradation.xlsx")

print("Top genes saved to Excel files.")
