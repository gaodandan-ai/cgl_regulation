# search_data_files.py
import os
import pandas as pd

data_dir = 'data'
files = os.listdir(data_dir)
print("Files in data directory:", files)

for f in files:
    path = os.path.join(data_dir, f)
    if os.path.isfile(path) and f.endswith('.csv'):
        try:
            df = pd.read_csv(path, nrows=5)
            print(f"\n--- {f} columns ---")
            print(df.columns.tolist())
            print(df.head(2))
        except Exception as e:
            print(f"Error reading {f}: {e}")
