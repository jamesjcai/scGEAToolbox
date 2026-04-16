#!/usr/bin/env python3
# coding: utf-8
# functions/download_cas_pubchem.py

import sys
import os
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
import csv
import pubchempy as pcp

# Download 3D SDF from PubChem for given CAS number
def download_sdf(cas_number, out_sdf):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas_number}/SDF?record_type=3d"
    r = requests.get(url)
    if r.status_code == 200 and len(r.content) > 0:
        with open(out_sdf, "wb") as f:
            f.write(r.content)
        return True
    else:
        print(f"[WARN] CAS {cas_number} not found in PubChem.")
        return False

# Convert SDF to PDBQT using RDKit and OpenBabel
def sdf_to_pdbqt(sdf_file, out_dir):
    mol_name = os.path.splitext(os.path.basename(sdf_file))[0]
    pdb_file = os.path.join(out_dir, f"{mol_name}.pdb")
    pdbqt_file = os.path.join(out_dir, f"{mol_name}.pdbqt")

    # Generate 3D conformer by RDKit
    mol = Chem.MolFromMolFile(sdf_file, removeHs=False)
    if mol is None:
        raise ValueError(f"RDKit failed to read {sdf_file}")

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.UFFOptimizeMolecule(mol)
    Chem.MolToPDBFile(mol, pdb_file)

    # Convert structure into PDBQT format by OpenBabel
    subprocess.run(["obabel", pdb_file, "-O", pdbqt_file, "-xh"], check=True)

    return pdbqt_file

def main():
    if len(sys.argv) != 3:
        print("Usage: download_cas_pubchem.py <cas_txt_file> <output_dir>")
        sys.exit(1)
    
    cas_txt_file = sys.argv[1]
    output_dir = sys.argv[2]
    os.makedirs(output_dir, exist_ok=True)

    work_dir = os.getcwd()
    csv_file = os.path.join(work_dir, "Output_drug_information.csv")

    with open(cas_txt_file) as f:
        cas_numbers = [line.strip() for line in f if line.strip()]

    with open(csv_file, "w", newline="") as cf:
        writer = csv.writer(cf)
        writer.writerow(["CAS number"])
        for cas in cas_numbers:
            writer.writerow([cas])

    print(f"[INFO] Output_drug_information.csv created at {csv_file}")

    for cas in cas_numbers:
        sdf_file = os.path.join(output_dir, f"{cas}.sdf")
        if not os.path.exists(sdf_file):
            success = download_sdf(cas, sdf_file)
            if not success:
                continue
        try:
            pdbqt_file = sdf_to_pdbqt(sdf_file, output_dir)
            print(f"[INFO] Generated PDBQT: {pdbqt_file}")
        except subprocess.CalledProcessError:
            print(f"[ERROR] Failed to convert {cas} to PDBQT.")

    print("[INFO] Fetching PubChem CID for each CAS number...")

    updated_rows = []

    for cas in cas_numbers:
        cid = ""
        try:
            compounds = pcp.get_compounds(cas, 'name')
            if compounds:
                cid = compounds[0].cid
        except Exception as e:
            print(f"[WARN] Failed to fetch CID for {cas}: {e}")

        updated_rows.append([cas, cid])

    with open(csv_file, "w", newline="") as cf:
        writer = csv.writer(cf)
        writer.writerow(["CAS number", "PubChem CID"])
        writer.writerows(updated_rows)

    print(f"[INFO] CAS and CID table written to {csv_file}")


if __name__ == "__main__":
    main()

