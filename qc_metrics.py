"""
Script Name: qc_metrics.py
Project: Pediatric IDH-wildtype Glioblastoma Stem Cells (GSE131928)
Phase: Phase 3 Exploratory Data Analysis & Quality Control (Task C)
Author: Sukhraj Kanwal (Sukhraj-Kanwal-Research)
Description: 
    Parses the cleaned pediatric metadata, outputs patient distribution, 
    and runs a strict contamination audit for MGH85 cell lines.
"""

import os
import pandas as pd

def run_qc_pipeline():
    data_path = os.path.join("data", "gse131928_pediatric_metadata_clean.csv")
    if not os.path.exists(data_path):
        desktop_fallback = r"C:\Users\Owner\Desktop\gse131928_pediatric_metadata_clean.csv"
        if os.path.exists(desktop_fallback):
            data_path = desktop_fallback
        else:
            print(f"[ERROR] Could not locate file at {data_path}")
            return

    print("[+] Loading cleaned pediatric cell metadata...")
    col_names = ['cell_id', 'cell_alias', 'cell_type', 'organism', 'strategy', 'source_file', 'platform', 'patient_id', 'cohort']
    df = pd.read_csv(data_path, header=None, names=col_names)

    print("\n========================================================")
    print("      PEDIATRIC METADATA & CONTAMINATION AUDIT          ")
    print("========================================================")
    print(f"[+] Total Verified Pediatric Cells: {len(df)}")

    # Strict Contamination Audit for MGH85
    mgh85_matches = df[df.astype(str).apply(lambda col: col.str.contains('MGH85', case=False, na=False)).any(axis=1)]
    
    print("\n[+] Running MGH85 Contamination Check...")
    if len(mgh85_matches) > 0:
        print(f"[!] CRITICAL ERROR: Found {len(mgh85_matches)} rows containing 'MGH85' in the pediatric dataset!")
    else:
        print("[✔] Audit Passed: Zero traces of MGH85 detected. Data purity verified.")

    print("\n[+] Patient Breakdown Verification:")
    print("-" * 40)
    patient_counts = df['patient_id'].value_counts()
    for patient, count in patient_counts.items():
        print(f"Patient {patient:<8} | {count:>4} cells")
    print("-" * 40)

    print("\n[✔] Task C Metadata & Audit Complete.")

if __name__ == "__main__":
    run_qc_pipeline()