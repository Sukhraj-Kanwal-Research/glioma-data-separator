import os
import glob
import pandas as pd
import numpy as np

print("[+] Locating normalized expression matrix and cluster labels...")
expr_matches = glob.glob(r"data\processed\pediatric_expression_hvg.tsv")
cluster_matches = glob.glob(r"data\processed\pediatric_umap_clusters.tsv")

if not expr_matches or not cluster_matches:
    raise FileNotFoundError("Could not find expression matrix or cluster files in data\\processed\\!")

expr_path = expr_matches[0]
cluster_path = cluster_matches[0]

print(f"[✔] Found expression matrix: {expr_path}")
print(f"[✔] Found cluster assignments: {cluster_path}")

print("[+] Loading data into memory...")
expr_df = pd.read_csv(expr_path, sep="\t", index_col=0)
clusters_df = pd.read_csv(cluster_path, sep="\t", index_col=0)

# Transpose expression matrix so rows are cells (1934) and columns are genes (2000)
expr_df = expr_df.T

print(f"[+] Transposed Expression shape: {expr_df.shape} (Cells x Genes)")
print(f"[+] Clusters shape: {clusters_df.shape} (Cells x Metadata)")

# Align by position since both dataframes represent the exact same cells in order
if len(expr_df) == len(clusters_df):
    clusters_df.index = expr_df.index
else:
    raise ValueError(f"Row mismatch: Expression has {len(expr_df)} rows, but Clusters has {len(clusters_df)} rows!")

# Combine expression data with cluster labels
expr_df['Cluster'] = clusters_df['Cluster']
print(f"[✔] Data aligned successfully. Combined shape: {expr_df.shape}")

print("[+] Calculating top marker genes for each cluster...")
marker_results = []
unique_clusters = sorted(expr_df['Cluster'].unique())

for cluster in unique_clusters:
    print(f"    -> Analyzing {cluster}...")
    
    # Separate cells inside the cluster vs outside the cluster
    in_cluster = expr_df[expr_df['Cluster'] == cluster].drop(columns=['Cluster'])
    out_cluster = expr_df[expr_df['Cluster'] != cluster].drop(columns=['Cluster'])
    
    # Calculate mean expression for genes inside vs outside
    mean_in = in_cluster.mean(axis=0)
    mean_out = out_cluster.mean(axis=0)
    
    # Calculate Log2 Fold Change (adding a small pseudocount to avoid division by zero)
    log_fc = np.log2(mean_in + 1e-5) - np.log2(mean_out + 1e-5)
    
    # Sort genes by highest log fold change
    top_genes = log_fc.sort_values(ascending=False).head(10)
    
    for gene, lfc in top_genes.items():
        marker_results.append({
            'Cluster': cluster,
            'Gene': gene,
            'Log2_Fold_Change': lfc
        })

# Convert results to a DataFrame
markers_df = pd.DataFrame(marker_results)

# Ensure output directory exists
os.makedirs(r"data\processed", exist_ok=True)

# Save the marker genes table
out_path = r"data\processed\pediatric_cluster_markers.tsv"
markers_df.to_csv(out_path, sep="\t", index=False)

print(f"\n[✔] Top marker genes successfully saved to: {out_path}")
print("\nTop 3 marker genes per cluster preview:")
print(markers_df.groupby('Cluster').head(3))
print("\n[✔] Phase 6: Marker Gene Discovery completed successfully!")