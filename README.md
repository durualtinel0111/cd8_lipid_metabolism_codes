# cd8_lipid_metabolism_codes
Supplementary Code for the manuscript: "Lipid Metabolic Reprogramming Drives CD8⁺ T Cell Exhaustion in PD-L1⁺ NSCLC: A Systems Biology and Causal Inference Approach"

## Overview

This repository contains Python supplementary code files for the analysis of lipid metabolism in CD8+ T cells, supporting the related research study on lipid-associated transcriptomic and immunometabolic profiling.

## Repository Structure

- `Code 1.py`: Identification of lipid-related genes and calculation of log2CPM values from raw count data.
- `Code 2.py`: Construction of gene-gene interaction matrix using Lasso regression.
- `Code 3.py`: Simulation of gene inhibition effects using ODE models.
- `Code 4.py`: Sensitivity analysis of the gene expression model.
- `Code 6.py`: Calculation of lipid accumulation scores from GSVA results.
- `Code 7.py`: Calculation of lipid degradation scores from GSVA results.
- `Code 8.py`: Causal inference analysis of lipid accumulation and exhaustion markers using DoWhy.
- `Code 9.py`: Causal inference analysis of lipid degradation and exhaustion markers using DoWhy.
- `Code 10.py`: Identification of top predictive genes for lipid metabolism using LightGBM.
- `requirements.txt`: List of Python dependencies.
- `LICENSE`: MIT License.

**Note:** The R script for GSVA pathway analysis (`Code 5.R`) has been moved to a separate repository and is **not included** here.

## Usage

- Python scripts require dependencies listed in `requirements.txt`. Install them via:

```bash
pip install -r requirements.txt
