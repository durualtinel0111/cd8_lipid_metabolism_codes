"""
Supplementary Code 1:
Identification of Lipid-Related Genes and Log2CPM Calculation

This script processes an HTSeq-count formatted RNA-seq file to:
1. Compute log2CPM expression values
2. Identify genes annotated with lipid-related GO Biological Process terms
3. Export log2CPM values of lipid-related genes with nonzero expression
"""

# 1. Install required library (if running in Colab)
# !pip install mygene

# 2. Import libraries
import pandas as pd
import numpy as np
from mygene import MyGeneInfo
from google.colab import files

# 3. Upload the HTSeq count file
uploaded = files.upload()
file_name = list(uploaded.keys())[0]

# 4. Parse the HTSeq-formatted counts into a DataFrame
with open(file_name, 'r') as f:
    lines = f.readlines()

split_rows = [line.strip().split() for line in lines if line.strip() and not line.startswith("__")]
df = pd.DataFrame(split_rows, columns=["Gene_ID", "raw_count"])
df["raw_count"] = pd.to_numeric(df["raw_count"], errors='coerce').fillna(0)

# 5. Compute log2CPM values
total = df["raw_count"].sum()
df["CPM"] = df["raw_count"] / total * 1e6
df["log2CPM"] = np.log2(df["CPM"] + 1)

# 6. Clean Ensembl IDs by removing version suffix
df["Gene_ID_clean"] = df["Gene_ID"].str.split(".").str[0]
gene_ids = df["Gene_ID_clean"].tolist()

# 7. Identify lipid-related genes using mygene.info (GO:BP)
mg = MyGeneInfo()
lipid_related = []
for i in range(0, len(gene_ids), 100):  # batch query to avoid overload
    batch = gene_ids[i:i+100]
    res = mg.querymany(batch, scopes="ensembl.gene", fields="go.BP.term", species="human")
    for r in res:
        terms = r.get("go", {}).get("BP", [])
        if isinstance(terms, dict): terms = [terms]
        for t in terms:
            if "lipid" in t["term"].lower():
                lipid_related.append(r["query"])
                break

# 8. Filter for lipid-related genes with non-zero log2CPM
filtered_df = df[df["Gene_ID_clean"].isin(lipid_related) & (df["log2CPM"] != 0)]

# 9. Export result
output_file = "lipid_related_log2cpm_nonzero.csv"
filtered_df[["Gene_ID", "log2CPM"]].to_csv(output_file, index=False)
files.download(output_file)
