import pandas as pd
from scipy.sparse import csr_matrix
import h5py
import anndata as ad
import os.path

import panhumanpy as ph


f = h5py.File("X.mat",'r')
# counts = csr_matrix(f.get('X'))
Xnorm = csr_matrix(f.get('Xnorm'))


g = f['g']
gene_names = []
for r in g:
    for ref in r:
        str_data = ''.join(chr(c[0]) for c in f[ref][:])
        gene_names.append(str_data)

f.close()

adata = ad.AnnData(X=Xnorm)
adata.obs.index = [f"cell_{i}" for i in range(adata.n_obs)]
adata.var.index = gene_names
print("Input data ready.")

# High-level interface
azimuth = ph.AzimuthNN(adata)
#azimuth = ph.AzimuthNN(adata, 
#                       feature_name_col='gene_symbols'
#                       eval_batch_size=25000)

cell_metadata = azimuth.cells_meta
cell_metadata.to_csv("output.csv", sep=",", index=True);
