import os
import numpy as np
import pandas as pd

input_path = r"C:\Users\Owner\Desktop\data\processed\filtered_expression_matrix.csv"
output_dir = r"C:\Users\Owner\Desktop\data\processed"
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, "qc_log1p_hvg_matrix.csv")

print("[+] Loading expression matrix (Cells x Genes)...")
df = pd.read_csv(input_path, index_col=0)
print(f"[+] Starting shape: {df.shape[0]} cells x {df.shape[1]} genes")

# 1. Basic QC counts on TPM matrix
genes_per_cell = (df > 0).sum(axis=1)
cells_per_gene = (df > 0).sum(axis=0)
total_tpm_per_cell = df.sum(axis=1)

print(f"[i] Genes/cell median: {genes_per_cell.median():.0f} (min={genes_per_cell.min()}, max={genes_per_cell.max()})")
print(f"[i] Cells/gene median: {cells_per_gene.median():.0f} (max={cells_per_gene.max()})")
print(f"[i] Total TPM/cell median: {total_tpm_per_cell.median():.0f}")

# 2. Gentle quality filtering (keeps biological signal while dropping near-empty noise)
min_genes_per_cell = 200
min_cells_per_gene = 3

keep_cells = genes_per_cell >= min_genes_per_cell
keep_genes = cells_per_gene >= min_cells_per_gene

df_filt = df.loc[keep_cells, keep_genes]
print(f"[+] After QC filtering (genes/cell >= {min_genes_per_cell}, cells/gene >= {min_cells_per_gene}): "
      f"{df_filt.shape[0]} cells x {df_filt.shape[1]} genes")

# 3. Log1p transform to compress wide TPM range
df_log = np.log1p(df_filt)

# 4. Select top highly variable genes (HVGs) for clean, fast PCA/UMAP
n_top_genes = min(2000, df_log.shape[1])
gene_vars = df_log.var(axis=0)
top_genes = gene_vars.nlargest(n_top_genes).index
df_hvg = df_log[top_genes]

print(f"[+] Log1p max value: {df_log.to_numpy().max():.2f} | Selected top {df_hvg.shape[1]} HVGs")

df_hvg.to_csv(output_path)
print(f"[SAVED] Saved clean QC + log1p + HVG matrix to: {output_path}")