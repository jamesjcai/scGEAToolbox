import os
import pandas as pd
#import scanpy as sc
#import scipy.io as spio
#import numpy as np
from scipy.sparse import csr_matrix
import h5py
import anndata

os.chdir("./")

f = h5py.File("X.mat",'r')
counts = csr_matrix(f.get('X'))
f.close()
# counts = csr_matrix(f.get('X'), dtype=np.float64)
# counts = np.array(f.get('X'), dtype=np.float64)
# N = f.get('n')    # or f['n']
# n=N[()].astype(int).item()
# data = f.get('/X')[()]
# sample_labels=f.get('/batchid')[:,0].astype(int)

adata = anndata.AnnData(X=counts)
# adata = anndata.AnnData(X=X.transpose().tocsr())
metadata = pd.read_csv("c.csv")
with open("g.csv",'r') as f:
          gene_names = f.read().splitlines()

metadata.index = metadata['CellID'].astype(str).tolist()
adata.obs = metadata
adata.var.index = gene_names
adata.write("output.h5ad")
