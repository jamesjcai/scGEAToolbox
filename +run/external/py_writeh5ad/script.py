import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

import pandas as pd 
import scanpy as sc 
import scipy.io as spio

#import sys
#fname=str(sys.argv[1])
#ada=sc.read_h5ad(fname)
#ada.obs.to_csv("obs.csv")
#ada.var.to_csv("var.csv") 
#ada.write("X.h5")


import numpy as np
import h5py
f = h5py.File('input.mat','r')
data = f.get('/X')[()]
# data = np.array(data) # For converting to a NumPy array
sample_labels=f.get('/batchid')[:,0].astype(int)
f.close()


adata = anndata.AnnData(X=X.tranpose().tocsr())
metadata = pd.read_csv("metadata.csc")
with open("gene_names.csv,'r') as f:
          gene_names = f.read().splitlines()

adata.obs = metadata
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names


