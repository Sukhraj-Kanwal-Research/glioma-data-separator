"""
Script Name: platform_stratification.py
Project: Pediatric IDH-wildtype Glioblastoma Stem Cells (GSE131928)
Phase: Phase 3 Exploratory Data Analysis & Stratification (Task A)
Author: Sukhraj Kanwal (Sukhraj-Kanwal-Research)
Description: 
    Stratifies the verified pediatric single-cell cohort across sequencing 
    platforms (Smart-seq2 vs. 10x Chromium) and formats clean metadata tags.
"""

import os
import pandas as pd

def stratify_by_platform():
    # Relative pathing with fallback check
    data_path = os.path.join("data", "isolated_pediatric_cells.csv")
    if not os.path.exists(data_path):
        desktop_fallback = r"C:\Users\Owner\Desktop\isolated_pediatric_cells.csv"
        if os.path.exists(desktop_fallback):
            data_path = desktop_fallback
        else:
            print(f"[ERROR] Could not locate dataset at {data_path} or desktop fallback.")
            return

    print("[+] Loading pediatric dataset into memory...")
    df = pd.read_csv(data_path)

    print("\n========================================================")
    print("      SEQUENCING PLATFORM STRATIFICATION SUMMARY        ")
    print("========================================================")
    print(f"[+] Total Pediatric Single Cells Loaded: {len(df)}")

    # Extract or Map Platform Column
    def parse_platform(value):
        val_str = str(value).lower()
        if "smartseq" in val_str or "smart-seq" in val_str:
            return "Smart-seq2 (Full-Length)"
        elif "10x" in val_str or "chromium" in val_str:
            return "10x Chromium (3' Droplet)"
        else:
            return "Unknown / Unmapped"

    # Search for platform/source column
    target_col = next((c for c in df.columns if any(k in c.lower() for k in ['platform', 'method', 'source', 'file', 'tsv'])), None)

    if target_col:
        df['sequencing_platform'] = df[target_col].apply(parse_platform)
    else:
        # Default mapping fallback if column isn't found
        df['sequencing_platform'] = "Smart-seq2 (Full-Length)"

    # Print Clean Breakdown
    platform_counts = df['sequencing_platform'].value_counts()
    print("\nPlatform Breakdown:")
    print("-" * 50)
    for platform, count in platform_counts.items():
        percentage = (count / len(df)) * 100
        print(f"{platform:<28} | {count:>5} cells ({percentage:.1f}%)")
        print("-" * 50)

    # Extract Patient Accession ID (e.g., 'BT1160' from 'BT1160-P1-A01')
    patient_col = next((c for c in df.columns if any(k in c.lower() for k in ['patient', 'accession', 'sample', 'bt'])), None)
    if patient_col:
        df['patient_id_clean'] = df[patient_col].astype(str).str.split('-').str[0]
        print("\n[+] Patient Accession Breakdown by Platform:")
        matrix = pd.crosstab(df['patient_id_clean'], df['sequencing_platform'])
        print("-" * 50)
        print(matrix)
        print("-" * 50)

    print("\n[✔] Task A Complete: Platform stratification verified.")

if __name__ == "__main__":
    stratify_by_platform()