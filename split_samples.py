import pandas as pd
import os

print("\n--- PEDIATRIC GLIOMA DATA SEPARATION PIPELINE ---")

# Define paths to your desktop
desktop_path = os.path.expanduser("~/Desktop")
metadata_file = os.path.join(desktop_path, "GSE131928_cell_metadata.xlsx")

if not os.path.exists(metadata_file):
    print(f"\n[!] Missing Data File: Ensure file is named 'GSE131928_cell_metadata.xlsx' on your Desktop.")
else:
    print("[+] Excel master matrix detected. Scanning rows for patient keys...")
    
    # Load file, skipping potential junk header rows to find the true data table
    df = pd.read_excel(metadata_file, header=1)
    
    # If the first row skip didn't work, clean column headers dynamically
    if 'adult_or_pediatric' not in df.columns:
        # Check if the columns are in a different row or labeled slightly differently
        for col in df.columns:
            if df[col].astype(str).str.lower().isin(['adult', 'pediatric']).any():
                df = df.rename(columns={col: 'adult_or_pediatric'})
                break

    # Final Academic Fallback Check
    if 'adult_or_pediatric' not in df.columns:
        print("[!] Columns found:", list(df.columns)[:5])
        print("[!] Error: Could not dynamically locate the adult/pediatric data label row.")
        print("[*] Re-reading sheet completely to trace cell tags...")
        df = pd.read_excel(metadata_file, header=0)
        # Scan every column in the entire spreadsheet to find where the text "pediatric" is hiding
        for col in df.columns:
            if df[col].astype(str).str.lower().str.contains('pediat', na=False).any():
                df = df.rename(columns={col: 'adult_or_pediatric'})
                print(f"[+] Found pediatric data tags hidden in column: '{col}'")
                break

    # 2. Isolate pediatric profiles based on our discovered metadata keys
    pediatric_df = df[df['adult_or_pediatric'].astype(str).str.lower().str.contains('pediat', na=False)]
    adult_df = df[df['adult_or_pediatric'].astype(str).str.lower().str.contains('adult', na=False)]
    
    # 3. Save clean filtered data sheets back to your Desktop as clean CSVs
    pediatric_output = os.path.join(desktop_path, "isolated_pediatric_cells.csv")
    adult_output = os.path.join(desktop_path, "isolated_adult_cells.csv")
    
    pediatric_df.to_csv(pediatric_output, index=False)
    adult_df.to_csv(adult_output, index=False)
    
    print("\n[SUCCESS] Pipeline executed perfectly!")
    print(f"-> Successfully isolated {len(pediatric_df)} pediatric single-cell rows.")
    print(f"-> Successfully isolated {len(adult_df)} adult single-cell control rows.")
    print(f"[*] Saved clear data sheet directly to Desktop: 'isolated_pediatric_cells.csv'")
    print("-------------------------------------------------\n")