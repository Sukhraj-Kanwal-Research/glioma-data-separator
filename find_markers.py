import os
import numpy as np
import pandas as pd

hvg_path = r"C:\Users\Owner\Desktop\data\processed\qc_log1p_hvg_matrix.csv"
umap_path = r"C:\Users\Owner\Desktop\data\processed\umap_clusters.csv"
output_dir = r"C:\Users\Owner\Desktop\data\processed"
os.makedirs(output_dir, exist_ok=True)

print("[+] Loading log1p HVG matrix and cluster labels...")
df_hvg = pd.read_csv(hvg_path, index_col=0)
df_clusters = pd.read_csv(umap_path, index_col=0)

# Align cell order safely
common_cells = df_hvg.index.intersection(df_clusters.index)
X = df_hvg.loc[common_cells]
clusters = df_clusters.loc[common_cells, "cluster"]

top_n = 8
rows = []

for c in sorted(clusters.unique()):
    in_c = (clusters == c).to_numpy()
    out_c = ~in_c

    mean_in = X.to_numpy()[in_c].mean(axis=0)
    mean_out = X.to_numpy()[out_c].mean(axis=0)
    frac_in = (X.to_numpy()[in_c] > 0).mean(axis=0)
    frac_out = (X.to_numpy()[out_c] > 0).mean(axis=0)

    # Difference in mean log1p expression (cluster vs. rest)
    diff = mean_in - mean_out
    top_idx = np.argsort(diff)[::-1][:top_n]

    print(f"\n[+] Top {top_n} enriched genes in Cluster {c} (n={in_c.sum()}):")
    for rank, idx in enumerate(top_idx, start=1):
        gene = X.columns[idx]
        print(f"  {rank}. {gene:<15} | diff={diff[idx]:+.3f} | in={frac_in[idx]*100:4.1f}% | out={frac_out[idx]*100:4.1f}%")
        rows.append({
            "cluster": c,
            "rank": rank,
            "gene": gene,
            "mean_in": mean_in[idx],
            "mean_out": mean_out[idx],
            "diff_mean": diff[idx],
            "pct_in": frac_in[idx] * 100.0,
            "pct_out": frac_out[idx] * 100.0,
        })

markers_path = os.path.join(output_dir, "top_cluster_markers.csv")
pd.DataFrame(rows).to_csv(markers_path, index=False)
print(f"\n[SAVED] Saved top marker table to: {markers_path}")