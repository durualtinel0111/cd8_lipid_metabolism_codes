"""
Supplementary Code 4:
Simulation Sensitivity Analysis of Gene Expression Dynamics

Performs multiple runs of the gene expression simulation without inhibition
to assess stability and variability of gene expression trajectories.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Load the Lasso interaction matrix (XLSX format)
interaction_df = pd.read_excel('/content/lasso_genetic_network.xlsx', index_col=0, engine='openpyxl')
W = interaction_df.values
genes = interaction_df.index.tolist()
n_genes = len(genes)

# 2. Simulation parameters
decay_rate = 1.0
dt = 0.1
time_steps = 100

# 3. Sigmoid activation function
def sigmoid(z):
    return 1 / (1 + np.exp(-z))

# 4. Simulation function with optional inhibition
def simulate_ode(W, inhibited_idx=[], initial_state=None):
    x = np.random.rand(n_genes) if initial_state is None else initial_state.copy()
    for idx in inhibited_idx:
        x[idx] = 0.0
    traj = np.zeros((time_steps, n_genes))
    for t in range(time_steps):
        dxdt = -decay_rate * x + sigmoid(W @ x)
        x = x + dxdt * dt
        for idx in inhibited_idx:
            x[idx] = 0.0
        traj[t] = x
    return traj

# 5. Run 10 simulation replicates without gene inhibition
all_final_states = []

for run in range(10):
    traj = simulate_ode(W, inhibited_idx=[])
    all_final_states.append(traj[-1])  # last time point

final_states_array = np.array(all_final_states)

# 6. Compute statistics: mean, std dev, coefficient of variation (CV)
mean_exp = final_states_array.mean(axis=0)
std_exp = final_states_array.std(axis=0)
cv_exp = std_exp / mean_exp

# 7. Create DataFrame summarizing gene expression variability
stability_df = pd.DataFrame({
    'Gene': genes,
    'Mean Expression': mean_exp,
    'Std Dev': std_exp,
    'CV': cv_exp
}).sort_values('CV', ascending=False)

# 8. Display top 20 most sensitive genes
print("\nTop 20 Most Sensitive Genes (Highest CV):")
print(stability_df.head(20))

# 9. Plot distribution of gene expression CV
plt.figure(figsize=(12, 6))
sns.histplot(stability_df['CV'], bins=30, kde=True)
plt.title("Gene Expression Variability Across Simulations (Coefficient of Variation)")
plt.xlabel("Coefficient of Variation (CV)")
plt.ylabel("Number of Genes")
plt.show()
