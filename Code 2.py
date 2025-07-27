"""
Supplementary Code 2:
Lasso-Based Inference of Gene–Gene Interaction Matrix

This script infers putative regulatory interactions among lipid-associated genes
by fitting Lasso regression models to gene expression data (log2CPM matrix).
"""

import pandas as pd
import numpy as np
from sklearn.linear_model import LassoCV
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Load pre-processed gene expression data (rows = genes, columns = samples)
file_path = "/content/LOG2CPM Data.xlsx"
df = pd.read_excel(file_path, index_col=0)

# 2. Convert to NumPy array
X = df.values
genes = df.index.tolist()
n_genes = X.shape[0]

# 3. Initialize the interaction matrix
interaction_matrix = np.zeros((n_genes, n_genes))

# 4. Fit Lasso regression per gene
for i in range(n_genes):
    y = X[i, :]
    X_other = np.delete(X, i, axis=0).T

    model = LassoCV(cv=5, random_state=42, max_iter=5000)
    model.fit(X_other, y)
    coefs = model.coef_

    # Insert coefficients (symmetric matrix)
    interaction_matrix[i, :i] = coefs[:i]
    interaction_matrix[i, i+1:] = coefs[i:]

# 5. Convert to DataFrame for readability
interaction_df = pd.DataFrame(interaction_matrix, index=genes, columns=genes)

# 6. Report number of non-zero interactions
nonzero_count = (interaction_df != 0).sum().sum()
print(f"Number of non-zero interactions: {nonzero_count}")

# 7. Optional: visualize as heatmap
plt.figure(figsize=(14, 12))
sns.heatmap(interaction_df, cmap='coolwarm', center=0)
plt.title("Gene–Gene Interaction Matrix Inferred by Lasso Regression")
plt.tight_layout()
plt.show()

# 8. Export as CSV
interaction_df.to_csv("/content/lasso_genetic_network.csv")
