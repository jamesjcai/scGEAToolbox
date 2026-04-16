#!/usr/bin/env python3
# coding: utf-8

import os
import sys
import requests
import pandas as pd
from Bio.PDB import MMCIFParser, PDBIO, Select

def download_mmcif(pdb_id, out_dir):
    pdb_id = pdb_id.upper()
    os.makedirs(out_dir, exist_ok=True)
    cif_path = os.path.join(out_dir, f"{pdb_id}.cif")
    if not os.path.exists(cif_path):
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        r = requests.get(url)
        if r.status_code != 200:
            raise FileNotFoundError(f"{pdb_id} not found in RCSB")
        with open(cif_path, "wb") as f:
            f.write(r.content)
        print(f"✅ Downloaded {pdb_id}.cif")
    return cif_path

class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_model(self, model):
        return model.id == 0

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_atom(self, atom):
        return atom.get_altloc() in (" ", "A")

def extract_chain(cif_path, chain_id, out_pdb):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("s", cif_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_pdb, ChainSelect(chain_id))

def download_single_pdb_chain(pdb_chain, out_dir):
    pdb_id, chain_id = pdb_chain.split(".")
    pdb_id = pdb_id.upper()
    chain_id = chain_id.upper()

    subdir = os.path.join(out_dir, pdb_chain)
    os.makedirs(subdir, exist_ok=True)
    out_pdb = os.path.join(subdir, f"{pdb_id}_{chain_id}.pdb")

    if os.path.exists(out_pdb):
        return out_pdb

    try:
        cif_path = download_mmcif(pdb_id, subdir)
        extract_chain(cif_path, chain_id, out_pdb)
        print(f"✅ Built {out_pdb} from RCSB")
        return out_pdb
    except FileNotFoundError:
        return None  # signal to use AlphaFold

def process_csv(csv_file, out_root):
    df = pd.read_csv(csv_file)
    pdb_files = []
    for raw_ids in df["PDB_model"]:
        if pd.isna(raw_ids):
            continue
        ids = [x.strip() for x in raw_ids.split(";") if x.strip()]
        for id in ids:
            if not "." in id:
                continue
            out_pdb = download_single_pdb_chain(id, out_root)
            pdb_files.append((id, out_pdb))
    return pdb_files

if __name__ == "__main__":
    csv_file = sys.argv[1]
    out_dir = sys.argv[2]
    process_csv(csv_file, out_dir)