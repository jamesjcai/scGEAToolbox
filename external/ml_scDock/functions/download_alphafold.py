#!/usr/bin/env python3
# coding: utf-8
# functions/download_alphafold.py

import os, sys, requests, time

# Download protein structure(s) from AlphaFold
def download_alphafold_pdb(af_id, output_dir="downloads"):
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    # In the URL, the string “v6” can be replaced with the current AlphaFold Database version to retrieve the corresponding structure files.
    # At present (2026-02-05), the latest release is v6. The most recent version can be verified by consulting the AlphaFold Database release notes available on the official website.
    url = f"https://alphafold.ebi.ac.uk/files/{af_id}-model_v6.pdb"
    output_file = os.path.join(output_dir, f"{af_id}.pdb")
    
    r = requests.get(url)
    if r.status_code == 200:
        with open(output_file, "wb") as f:
            f.write(r.content)
        print(f"✅ Downloaded {af_id} to {output_file}")
    else:
        print(f"❌ Failed to download {af_id}, HTTP {r.status_code}")
    time.sleep(5)

if __name__ == "__main__":
    af_id, out_dir = sys.argv[1], sys.argv[2]
    download_alphafold_pdb(af_id, out_dir)

