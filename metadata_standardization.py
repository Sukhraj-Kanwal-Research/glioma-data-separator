import pandas as pd
import os

desktop_path = os.path.expanduser("~/Desktop")
pediatric_file = os.path.join(desktop_path, "isolated_pediatric_cells.csv")

print("--- EXPOSING TRUE CLINICAL COLUMNS ---")

# 1. Load the spreadsheet normally
df = pd.read_csv(pediatric_file)

# 2. Set the columns to be equal to the very first row of actual data values
df.columns = df.iloc[0]

# 3. Drop that first row since it's now our header, and drop any completely empty columns
df = df.iloc[1:].dropna(axis=1, how='all')

# 4. Save the cleanly labeled file right back over your old CSV
df.to_csv(pediatric_file, index=False)

print("\n[SUCCESS] Column labels repaired!")
print("Here are the real variables available for your research:")
for i, col in enumerate(df.columns, 1):
    print(f"{i}. {col}")
print("---------------------------------------\n")