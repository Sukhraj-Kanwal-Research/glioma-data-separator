import os
import glob
import pandas as pd
import numpy as np

print("[+] Locating filtered expression matrix...")
expr_matches = glob.glob(r"data\processed\pediatric_expression_filtered.tsv")
if not expr_matches:
    expr_matches = glob.glob(r"C:\Users\Owner\**\pediatric_expression_filtered.tsv", recursive=True)

if not expr_matches:
    raise FileNotFoundError("Could not find filtered expression matrix in data\\processed\\!")

expr_path = expr_matches[0]
print(f"[✔] Found expression matrix at: {expr_path}")

print("[+] Loading filtered expression data...")
expression_df = pd.read_csv(expr_path, sep="\t", index_col=0)
print(f"[+] Loaded matrix with shape: {expression_df.shape} (Genes x Cells)")

print("[+] Applying log-normalization (log1p)...")
# Log-normalize TPM values: ln(TPM + 1) to stabilize variance
expr_log = np.log1p(expression_df)

print("[+] Calculating gene variance for Highly Variable Gene (HVG) selection...")
# Calculate variance for each gene across all cells on the log scale
gene_var = expr_log.var(axis=1)

# Select top 2,000 highly variable genes
n_top_genes = min(2000, len(gene_var))
top_genes = gene_var.nlargest(n_top_genes).index

expression_hvg = expr_log.loc[top_genes]
print(f"[✔] Selected top {len(top_genes)} Highly Variable Genes out of {len(expression_df)} total genes.")

# Ensure output directory exists
os.makedirs(r"data\processed", exist_ok=True)

# Save outputs
log_out = r"data\processed\pediatric_expression_log_normalized.tsv"
hvg_out = r"data\processed\pediatric_expression_hvg.tsv"

print("[+] Saving normalized and HVG matrices...")
expr_log.to_csv(log_out, sep="\t")
expression_hvg.to_csv(hvg_out, sep="\t")

print(f"[✔] Full log-normalized matrix saved to: {log_out}")
print(f"[✔] HVG subset matrix (top 2000) saved to: {hvg_out}")
print("[✔] Phase 4 Task A: Log-Normalization & HVG Selection completed successfully!")