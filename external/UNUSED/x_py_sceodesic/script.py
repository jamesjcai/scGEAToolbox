#import os
#import sys
#abspath = os.path.abspath(__file__)
#dname = os.path.dirname(abspath)
# os.chdir(dname)
# os.chdir("D:\\GitHub\\scGEAToolbox\\+run\\external\\py_datamapplot")
# os.chdir("C:\\Users\\jcai\\Documents\\GitHub\\scGEAToolbox\\+run\\external\\py_geosketch")

import pandas as pd 
import anndata
import scanpy as sc
from sceodesic import run_sceo
from scipy.sparse import csr_matrix
import h5py
import numpy as np


f = h5py.File("X.mat",'r')
counts = csr_matrix(f.get('X'))
f.close()
adata = anndata.AnnData(X=counts)
metadata = pd.read_csv("c.csv")
with open("g.csv",'r') as f:
          gene_names = f.read().splitlines()

adata.obs = metadata
adata.obs.index = adata.obs['CellID'].tolist()
adata.var.index = gene_names

sc.pp.normalize_total(adata, 1e4)
sc.pp.log1p(adata)

run_sceo(adata, num_hvg=300)

embeddings = adata.obsm['sceo_embeddings']
# programs (loadings) stored in .varm, as a numpy array
programs = adata.varm['sceo_programs']

# Create new anndata object with embeddings
adata_sceo = anndata.AnnData(embeddings, obs=adata.obs)
sc.get.rank_genes_groups_df(adata = adata_sceo, group = "CellType")
adata_sceo.write("output.h5ad")
