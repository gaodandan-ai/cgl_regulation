# check_binding_site.py
import pandas as pd

df = pd.read_csv('data/regulations.csv')
non_empty = df[df['Binding_site'].notna() & (df['Binding_site'] != '')]
print(f"Total rows in regulations.csv: {len(df)}")
print(f"Rows with non-empty binding sites: {len(non_empty)}")
if len(non_empty) > 0:
    print("\nSome examples of binding site entries:")
    print(non_empty[['TF_name', 'TG_name', 'Binding_site']].head(15))
else:
    print("All Binding_site entries are empty.")
