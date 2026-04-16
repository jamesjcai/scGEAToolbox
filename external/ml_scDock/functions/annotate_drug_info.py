#!/usr/bin/env python3
import sys
import pandas as pd

if len(sys.argv) != 3:
    print("Usage: annotate_drug_info.py <Output_drug_information.csv> <drug_information.tsv>")
    sys.exit(1)

cas_csv = sys.argv[1]
drug_tsv = sys.argv[2]

print(f"[INFO] Loading {cas_csv}")
cas_df = pd.read_csv(cas_csv)

print(f"[INFO] Loading {drug_tsv}")
drug_df = pd.read_csv(drug_tsv, sep="\t", low_memory=False)

def clean_cid(x):
    if pd.isna(x):
        return ""
    x = str(x).strip()
    if ";" in x:
        x = x.split(";")[0]
    if x.endswith(".0"):
        x = x[:-2]
    return x

print("[INFO] Cleaning PubChem CID format...")
cas_df["PubChem CID"] = cas_df["PubChem CID"].apply(clean_cid)
drug_df["PubChem CID"] = drug_df["PubChem CID"].apply(clean_cid)


if "PubChem CID" not in cas_df.columns:
    raise ValueError("Output_drug_information.csv missing 'PubChem CID' column")
if "PubChem CID" not in drug_df.columns:
    raise ValueError("drug_information.tsv missing 'PubChem CID' column")


print("[INFO] Merging by PubChem CID...")
merged_df = cas_df.merge(drug_df, on="PubChem CID", how="left")


matched = merged_df["PubChem CID"].isin(drug_df["PubChem CID"]).sum()
print(f"[INFO] Matched rows: {matched} / {len(merged_df)}")
if matched == 0:
    print("[WARNING] No rows matched! Check CID format in filtered.tsv")


cols = merged_df.columns.tolist()
new_order = ["CAS number", "PubChem CID"] + [c for c in cols if c not in ["CAS number", "PubChem CID"]]
merged_df = merged_df[new_order]

merged_df.to_csv(cas_csv, index=False)
print(f"[INFO] Updated {cas_csv} with drug information")