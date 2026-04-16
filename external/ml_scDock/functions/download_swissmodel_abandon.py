#!/usr/bin/env python3
# coding: utf-8
# functions/download_swissmodel.py

import os, sys, requests, time

def download_swissmodel_chain(pdb_id, chain, out_dir="."):
    pdb_id = pdb_id.lower().strip()
    chain = chain.upper().strip()
    filename = f"{pdb_id}.1.{chain}.pdb"
    url = f"https://www.swissmodel.expasy.org/templates/{filename}"
    
    os.makedirs(out_dir, exist_ok=True)
    file_path = os.path.join(out_dir, filename)
    
    r = requests.get(url)
    if r.status_code == 200:
        with open(file_path, "wb") as f:
            f.write(r.content)
        print(f"✅ Downloaded {filename} to {file_path}")
    else:
        print(f"❌ Failed to download {filename}, HTTP {r.status_code}")
    time.sleep(5)

if __name__ == "__main__":
    pdb_id, chain, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
    download_swissmodel_chain(pdb_id, chain, out_dir)

