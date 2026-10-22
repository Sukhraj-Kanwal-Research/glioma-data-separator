import os
import glob
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

print("[+] Locating HVG expression matrix...")
hvg_matches = glob.glob(r"data\processed\pediatric_expression_hvg.tsv")
if not hvg_matches:
    hvg_matches = glob.glob(r"C:\Users\Owner\**\pediatric_expression_hvg.tsv", recursive=True)

if not hvg_matches:
    raise FileNotFoundError("Could not find HVG expression matrix in data\\processed\\!")

hvg_path = hvg_matches[0]
print(f"[✔] Found HVG matrix at: {hvg_path}")

print("[+] Loading HVG matrix...")
hvg_df = pd.read_csv(hvg_path, sep="\t", index_col=0)
print(f"[+] Loaded matrix with shape: {hvg_df.shape} (Genes x Cells)")

print("[+] Transposing matrix for PCA (Cells x Genes)...")
# Scikit-learn expects rows to be samples (cells) and columns to be features (genes)
cell_gene_matrix = hvg_df.T

print("[+] Running Principal Component Analysis (PCA)...")
# We'll compute 50 principal components
n_components = 50
pca = PCA(n_components=n_components, random_state=42)
cell_pca_coords = pca.fit_transform(cell_gene_matrix)

# Create a DataFrame for the PCA coordinates
pca_columns = [f"PC{i+1}" for i in range(n_components)]
pca_df = pd.DataFrame(cell_pca_coords, index=cell_gene_matrix.index, columns=pca_columns)

# Check explained variance
explained_variance = pca.explained_variance_ratio_
total_var_captured = np.sum(explained_variance) * 100
print(f"[✔] PCA complete! Captured {total_var_captured:.2f}% of total variance across {n_components} components.")

# Ensure output directory exists
os.makedirs(r"data\processed", exist_ok=True)

# Save PCA coordinates
pca_out = r"data\processed\pediatric_pca_coordinates.tsv"
pca_df.to_csv(pca_out, sep="\t")

print(f"[✔] PCA coordinates saved to: {pca_out}")
print("[✔] Phase 5 Task A: Dimensionality Reduction completed successfully!")