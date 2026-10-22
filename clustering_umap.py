import os
import glob
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
import umap

print("[+] Locating PCA coordinates matrix...")
pca_matches = glob.glob(r"data\processed\pediatric_pca_coordinates.tsv")
if not pca_matches:
    pca_matches = glob.glob(r"C:\Users\Owner\**\pediatric_pca_coordinates.tsv", recursive=True)

if not pca_matches:
    raise FileNotFoundError("Could not find PCA coordinates matrix in data\\processed\\!")

pca_path = pca_matches[0]
print(f"[✔] Found PCA matrix at: {pca_path}")

print("[+] Loading PCA coordinates...")
pca_df = pd.read_csv(pca_path, sep="\t", index_col=0)
print(f"[+] Loaded PCA shape: {pca_df.shape} (Cells x Principal Components)")

print("[+] Performing KMeans Clustering on PCA coordinates...")
# Group cells into clusters (e.g., 5 distinct biological states / tumor subclusters)
n_clusters = 5
kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
cluster_labels = kmeans.fit_predict(pca_df)

# Add cluster labels to the DataFrame
pca_df['Cluster'] = [f"Cluster_{label}" for label in cluster_labels]
print(f"[✔] Grouped cells into {n_clusters} distinct clusters.")

print("[+] Running UMAP for 2D visualization projection...")
# Project 50 dimensions down to 2 dimensions (UMAP1, UMAP2)
reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
umap_coords = reducer.fit_transform(pca_df.drop(columns=['Cluster']))

# Create a final results DataFrame
results_df = pd.DataFrame(umap_coords, index=pca_df.index, columns=["UMAP1", "UMAP2"])
results_df['Cluster'] = pca_df['Cluster']

# Ensure output directory exists
os.makedirs(r"data\processed", exist_ok=True)

# Save UMAP and cluster outputs
out_path = r"data\processed\pediatric_umap_clusters.tsv"
results_df.to_csv(out_path, sep="\t")

print(f"[✔] UMAP coordinates and clusters saved to: {out_path}")
print("[✔] Phase 5 Task B: Clustering & UMAP Visualization completed successfully!")