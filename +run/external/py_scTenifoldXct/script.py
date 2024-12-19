import sys
import scanpy as sc
import scTenifoldXct as st
from scTenifoldXct.dataLoader import build_adata
import pandas as pd
import numpy as np
import h5py
import scipy
from scipy import sparse

twosided = 1
with h5py.File('file.h5', 'r') as f:
    if 'twosided' in f:
        twosided = f['twosided'][0][0]

#f = h5py.File('A1.mat','r')
#A = np.array(f.get('A'))
#A = sparse.csc_matrix(A)
#sparse.save_npz("pcnet_Source.npz", A)

#f = h5py.File('A2.mat','r')
#A = np.array(f.get('A'))
#A = sparse.csc_matrix(A)
#sparse.save_npz("pcnet_Target.npz", A)

adata = build_adata("X.mat", "g.txt", "c.txt", delimiter=',', meta_cell_cols=['cell_type'], transpose=False)
print('Input read.............')
xct = st.scTenifoldXct(data = adata, 
                    source_celltype = 'Source',
                    target_celltype = 'Target',
                    obs_label = "cell_type",
                    rebuild_GRN = False,
                    GRN_file_dir = './',
                    verbose = True)
emb = xct.get_embeds(train = True)
xct_pairs = xct.null_test(pval = 1.0)
# print(xct_pairs)
pd.DataFrame(xct_pairs).to_csv('output1.txt',index=False,header=True)

#if sys.argv[1]=="1":
if twosided==1:
    xct = st.scTenifoldXct(data = adata, 
                        source_celltype = 'Target',
                        target_celltype = 'Source',
                        obs_label = "cell_type",
                        rebuild_GRN = False,
                        GRN_file_dir = './',
                        verbose = True)
    emb = xct.get_embeds(train = True)
    xct_pairs = xct.null_test(pval = 1.0)
    # print(xct_pairs)
    pd.DataFrame(xct_pairs).to_csv('output2.txt',index=False,header=True)
