"""
Supplementary Code 3:
Gene Inhibition Simulation Using Lasso-Inferred Interaction Matrix

Simulates gene expression dynamics under continuous inhibition of each gene
using an ODE model with sigmoid activation and passive decay.
"""

import numpy as np
import pandas as pd

# 1. Load the gene–gene interaction matrix from Lasso regression
interaction_df = pd.read_csv('/content/lasso_genetic_network.csv', index_col=0)
W = interaction_df.values
genes = interaction_df.index.tolist()
n_genes = len(genes)

# 2. Define simulation parameters
decay_rate = 1.0   # Decay rate λ in dx/dt = -λx + σ(Wx)
dt = 0.1           # Time step for Euler integration
time_steps = 100   # Number of time points to simulate

# 3. Sigmoid activation function
def sigmoid(z):
    return 1 / (1 + np.exp(-z))

# 4. ODE simulation function with continuous inhibition of a gene
def simulate_ode(W, inhibited_idx):
    x = np.random.rand(n_genes)       # Initial expression (random)
    x[inhibited_idx] = 0.0            # Set inhibited gene expression to zero
    traj = np.zeros((time_steps, n_genes))

    for t in range(time_steps):
        dxdt = -decay_rate * x + sigmoid(W @ x)
        x = x + dxdt * dt
        x[inhibited_idx] = 0.0        # Enforce inhibition at all time points
        traj[t] = x

    return traj

# 5. Run simulation for each gene inhibition and store trajectories
all_results = {}

for i, gene in enumerate(genes):
    print(f"Simulating inhibition of gene: {gene}")
    traj = simulate_ode(W, i)
    all_results[gene] = traj

# 6. Export all results to Excel (one sheet per gene)
with pd.ExcelWriter('/content/all_gen_inhibition_simulation.xlsx') as writer:
    for gene, traj in all_results.items():
        df = pd.DataFrame(traj, columns=genes)
        df.index.name = 'TimeStep'
        df.to_excel(writer, sheet_name=gene[:31])  # Excel sheet name limit: 31 chars

print("Simulation complete. Results saved to '/content/all_gen_inhibition_simulation.xlsx'")
