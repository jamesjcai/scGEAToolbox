import h5py
from anndata import AnnData

#with h5py.File('example_data.h5', 'r') as f:
#    data = f['/data'][:]
#adata = AnnData(data)
#adata.write('example_data.h5ad')
#print('Data saved as .h5ad file successfully.')

import scanpy as sc

file_path = 'example_data.h5ad'
adata = sc.read(file_path)

# Access the count matrix (adata.X), cell information (adata.obs), and gene information (adata.var)
