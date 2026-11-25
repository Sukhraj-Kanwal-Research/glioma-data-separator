import os
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

input_path = r"C:\Users\Owner\Desktop\data\processed\qc_log1p_hvg_matrix.csv"
output_dir = r"C:\Users\Owner\Desktop\data\processed"
os.makedirs(output_dir, exist_ok=True)

pca_coords_path = os.path.join(output_dir, "pca_coordinates.csv")
variance_path = os.path.join(output_dir, "pca_explained_variance.csv")

print("[+] Loading QC + log1p + HVG matrix...")
df_hvg = pd.read_csv(input_path, index_col=0)
print(f"[+] Loaded shape: {df_hvg.shape[0]} cells x {df_hvg.shape[1]} HVGs")

# 1. Standardize across cells so each HVG has mean=0, std=1
scaler = StandardScaler()
X_scaled = scaler.fit_transform(df_hvg.to_numpy(dtype=np.float64))

# 2. Fit PCA (cap components cleanly)
n_components = min(50, X_scaled.shape[0] - 1, X_scaled.shape[1])
pca = PCA(n_components=n_components, random_state=42)
X_pca = pca.fit_transform(X_scaled)

pc_cols = [f"PC{i+1}" for i in range(n_components)]
df_pca = pd.DataFrame(X_pca, index=df_hvg.index, columns=pc_cols)

exp_var = pca.explained_variance_ratio_
cum_var = np.cumsum(exp_var)
df_var = pd.DataFrame({
    "PC": pc_cols,
    "explained_variance_ratio": exp_var,
    "cumulative_variance_ratio": cum_var,
})

print("[+] Top 10 PCs explained variance:")
for i in range(min(10, n_components)):
    print(f"  - {pc_cols[i]}: {exp_var[i]*100:.2f}% (cumulative: {cum_var[i]*100:.2f}%)")

df_pca.to_csv(pca_coords_path)
df_var.to_csv(variance_path, index=False)
print(f"[SAVED] Saved PCA coordinates to: {pca_coords_path}")
print(f"[SAVED] Saved explained variance table to: {variance_path}")