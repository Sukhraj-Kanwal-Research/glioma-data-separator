import os
import glob
import pandas as pd
import numpy as np

print("[+] Locating raw files...")
meta_matches = glob.glob(r"data\raw\gse131928_pediatric_metadata_clean.csv")
if not meta_matches:
    meta_matches = glob.glob(r"C:\Users\Owner\**\gse131928_pediatric_metadata_clean.csv", recursive=True)
metadata_path = meta_matches[0]

expr_matches = glob.glob(r"data\raw\*Smartseq2*GBM*IDHwt*.tsv")
if not expr_matches:
    expr_matches = glob.glob(r"data\raw\GSM3828672*.tsv")
if not expr_matches:
    expr_matches = glob.glob(r"C:\Users\Owner\**\*Smartseq2*GBM*IDHwt*.tsv", recursive=True)
expression_path = expr_matches[0]

print(f"[✔] Metadata: {metadata_path}")
print(f"[✔] Expression Matrix: {expression_path}")

print("[+] Loading datasets...")
metadata = pd.read_csv(metadata_path)
expression_df = pd.read_csv(expression_path, sep="\t", index_col=0)

# Identify cell ID column
cell_col = next((col for col in ["cell_id", "Cell_ID", "sample", "Unnamed: 0"] if col in metadata.columns), metadata.columns[0])

# Align cells
common_cells = np.intersect1d(metadata[cell_col].astype(str).values, expression_df.columns.astype(str).values)
metadata_aligned = metadata[metadata[cell_col].astype(str).isin(common_cells)].copy()
expression_aligned = expression_df[common_cells].copy()

# Calculate QC metrics
total_tpm = expression_aligned.sum(axis=0)
n_features = (expression_aligned > 0).sum(axis=0)
mt_genes = [g for g in expression_aligned.index if str(g).upper().startswith("MT-")]
pct_mitochondrial = (expression_aligned.loc[mt_genes].sum(axis=0) / total_tpm) * 100 if mt_genes else pd.Series(0, index=common_cells)

qc_df = pd.DataFrame({
    "total_tpm": total_tpm,
    "n_features": n_features,
    "pct_mitochondrial": pct_mitochondrial
}, index=common_cells)

# Define QC Filtering Thresholds
min_features = 2000
max_features = 10000
max_mito = 10.0

passing_cells = qc_df[
    (qc_df["n_features"] >= min_features) & 
    (qc_df["n_features"] <= max_features) & 
    (qc_df["pct_mitochondrial"] <= max_mito)
].index

print(f"[+] Filtering cells: {len(common_cells)} initial -> {len(passing_cells)} passing QC.")

# Subset filtered data
metadata_filtered = metadata_aligned[metadata_aligned[cell_col].astype(str).isin(passing_cells)].copy()
expression_filtered = expression_aligned[passing_cells].copy()

# Ensure output directory exists
os.makedirs(r"data\processed", exist_ok=True)

# Save filtered outputs
meta_out = r"data\processed\pediatric_metadata_filtered.csv"
expr_out = r"data\processed\pediatric_expression_filtered.tsv"

metadata_filtered.to_csv(meta_out, index=False)
expression_filtered.to_csv(expr_out, sep="\t")

print(f"[✔] Filtered metadata saved to: {meta_out}")
print(f"[✔] Filtered expression matrix saved to: {expr_out}")
print("[✔] Phase 3 Task E Cell Filtering completed successfully!")