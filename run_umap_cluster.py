import os
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import umap

# --- YOUR ANALYTICAL CHOICES ---
N_PCS = 15        # Try 10, 15, or 20 based on your variance curve
N_CLUSTERS = 4    # Try 3, 4, 5, or 6 to see how cell states separate
# -------------------------------

pca_path = r"C:\Users\Owner\Desktop\data\processed\pca_coordinates.csv"
output_dir = r"C:\Users\Owner\Desktop\data\processed"
os.makedirs(output_dir, exist_ok=True)

print(f"[+] Loading top {N_PCS} PCs from {pca_path}...")
df_pca = pd.read_csv(pca_path, index_col=0)
X_pcs = df_pca.iloc[:, :N_PCS].to_numpy()

# 1. Run UMAP on the top PCs
reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
embedding = reducer.fit_transform(X_pcs)

# 2. Cluster cells in PC space
kmeans = KMeans(n_clusters=N_CLUSTERS, random_state=42, n_init=10)
clusters = kmeans.fit_predict(X_pcs)

# 3. Save coordinates + cluster labels
df_umap = pd.DataFrame({
    "UMAP1": embedding[:, 0],
    "UMAP2": embedding[:, 1],
    "cluster": clusters
}, index=df_pca.index)

# Extract sample prefix (e.g., '102' from '102_1') to see if clusters match patients/samples
df_umap["sample_id"] = df_umap.index.to_series().astype(str).str.split("_").str[0]

umap_csv_path = os.path.join(output_dir, "umap_clusters.csv")
df_umap.to_csv(umap_csv_path)
print(f"[SAVED] Saved UMAP + cluster table to: {umap_csv_path}")
print("\n[+] Cluster counts:")
print(df_umap["cluster"].value_counts().sort_index())
print("\n[+] Cross-tab of Sample ID vs. Cluster (check for patient/batch separation):")
print(pd.crosstab(df_umap["sample_id"], df_umap["cluster"]))

# 4. Save side-by-side plot (Clusters vs. Sample ID)
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

for c in sorted(df_umap["cluster"].unique()):
    sub = df_umap[df_umap["cluster"] == c]
    axes[0].scatter(sub["UMAP1"], sub["UMAP2"], s=14, alpha=0.8, label=f"Cluster {c}")
axes[0].set_title(f"UMAP by K-Means Cluster (k={N_CLUSTERS})")
axes[0].set_xlabel("UMAP1")
axes[0].set_ylabel("UMAP2")
axes[0].legend(markerscale=1.5, fontsize=8)

for s in sorted(df_umap["sample_id"].unique()):
    sub = df_umap[df_umap["sample_id"] == s]
    axes[1].scatter(sub["UMAP1"], sub["UMAP2"], s=14, alpha=0.8, label=f"Sample {s}")
axes[1].set_title("UMAP by Sample ID Prefix")
axes[1].set_xlabel("UMAP1")
axes[1].set_ylabel("UMAP2")
axes[1].legend(markerscale=1.5, fontsize=8, ncol=2)

plt.tight_layout()
plot_path = os.path.join(output_dir, "umap_clusters_vs_samples.png")
plt.savefig(plot_path, dpi=200)
print(f"[SAVED] Saved UMAP figure to: {plot_path}")