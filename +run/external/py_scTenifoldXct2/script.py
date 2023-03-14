import os
import sys
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
# os.chdir("U:\\GitHub\\scGEAToolbox\\+run\\external\\py_scTenifoldXct")
# os.chdir("d:\\GitHub\\scGEAToolbox\\+run\\external\\py_scTenifoldXct2")

import scanpy as sc
import scTenifoldXct as st
from scTenifoldXct.dataLoader import build_adata
import pandas as pd
import numpy as np
import h5py
import scipy
from scipy import sparse

#f = h5py.File('A1.mat','r')
#A = np.array(f.get('A'))
#A = sparse.csc_matrix(A)
#sparse.save_npz("pcnet_Source.npz", A)

#f = h5py.File('A2.mat','r')
#A = np.array(f.get('A'))
#A = sparse.csc_matrix(A)
#sparse.save_npz("pcnet_Target.npz", A)

adata = build_adata("X1.mat", "g1.txt", "c1.txt", delimiter=',', meta_cell_cols=['cell_type'], transpose=False)
print('Input read.............1')
xct1 = st.scTenifoldXct(data = adata, 
                    source_celltype = 'Source',
                    target_celltype = 'Target',
                    obs_label = "cell_type",
                    rebuild_GRN = False,
                    GRN_file_dir = './1',
                    verbose = True)

adata = build_adata("X2.mat", "g2.txt", "c2.txt", delimiter=',', meta_cell_cols=['cell_type'], transpose=False)
print('Input read.............2')
xct2 = st.scTenifoldXct(data = adata, 
                    source_celltype = 'Source',
                    target_celltype = 'Target',
                    obs_label = "cell_type",
                    rebuild_GRN = False,
                    GRN_file_dir = './2',
                    verbose = True)

XCTs = st.merge_scTenifoldXct(xct1, xct2)
emb = XCTs.get_embeds(train = True)
XCTs.nn_aligned_diff(emb) 
xcts_pairs_diff = XCTs.chi2_diff_test(pval = 1.0)
pd.DataFrame(xcts_pairs_diff).to_csv('output.txt',index=False,header=True)

