#import scanpy as sc 
#import scipy.io as spio
#import numpy as np

import os

#import sys
#abspath = os.path.abspath(__file__)
#dname = os.path.dirname(abspath)
# os.chdir(dname)
# os.chdir("D:\\GitHub\\scGEAToolbox\\+run\\external\\py_datamapplot")
os.chdir("z:\\Cailab\\GitHub\\scGEAToolbox\\+run\\external\\py_ikarus")

import urllib.request
import anndata
import pandas as pd
from pathlib import Path
from ikarus import classifier, utils, data

import scanpy as sc
from scipy.sparse import csr_matrix, issparse
import h5py
import numpy as np
from anndata import AnnData, read_h5ad



f = h5py.File("X2.mat",'r')
counts = csr_matrix(f.get('X'), dtype=np.float32)
f.close()
adata = anndata.AnnData(X=counts)
#metadata = pd.read_csv("c.csv")
with open("g2.csv",'r') as f:
          gene_names = f.read().splitlines()

#adata.obs = metadata
#adata.obs.index = adata.obs['CellID'].tolist()
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
cl = [f"Cell_{i:d}" for i in range(adata.n_obs)]
#adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
#adata.var.index = gene_names
#adata.var_names = gene_names
df = pd.DataFrame(gene_names, columns=['gene_symbol'])
df.index = gene_names
adata.var = df




if not issparse(adata.X):
    adata.X = csr_matrix(adata.X)
adata.var_names_make_unique()

#sc.pp.normalize_total(adata, 1e4)
#sc.pp.log1p(adata)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# adata.write_h5ad("adata.h5ad")

adata2 = AnnData(counts, var=df, obs=cl)
sc.pp.normalize_total(adata2, target_sum=1e4)
sc.pp.log1p(adata2)


signatures_path = Path("signatures.gmt")
pd.read_csv(signatures_path, sep="\t", header=None)
model_path = Path("core_model.joblib")
model = classifier.Ikarus(signatures_gmt=signatures_path, out_dir="out")
model.load_core_model(model_path)
_ = model.predict(adata2, "test", save=True)



#adata3 = read_h5ad("adata.h5ad")
#if not issparse(adata3.X):
#    adata3.X = csr_matrix(adata3.X)
#adata3.var_names_make_unique()
#_ = model.predict(adata3, "test3", save=True)



#adata2 = read_h5ad("adata2.h5ad")
#if not issparse(adata2.X):
#    adata2.X = csr_matrix(adata2.X)
#adata2.var_names_make_unique()
#_ = model.predict(adata2, "test", save=True)

