#script_path
#import os
#abspath = os.path.abspath(__file__)
#dname = os.path.dirname(abspath)
# os.chdir(dname)

import scTenifoldXct as st
from scTenifoldXct.dataLoader import build_adata
import pandas as pd
import h5py

twosided = 1
with h5py.File('X1.mat', 'r') as f:
    if 'targetgid' in f:
        targetgid = f['targetgid'][0]

adata1 = build_adata("X1.mat", "g1.txt", "c1.txt", delimiter=',', meta_cell_cols=['cell_type'], transpose=False)
print('Input read.............1')
xct1 = st.scTenifoldXct(data = adata1, 
                    source_celltype = 'Source',
                    target_celltype = 'Target',
                    obs_label = "cell_type",
                    rebuild_GRN = False,
                    GRN_file_dir = './1',
                    verbose = True)

adata2 = build_adata("X2.mat", "g2.txt", "c2.txt", delimiter=',', meta_cell_cols=['cell_type'], transpose=False)
print('Input read.............2')
xct2 = st.scTenifoldXct(data = adata2, 
                    source_celltype = 'Source',
                    target_celltype = 'Target',
                    obs_label = "cell_type",
                    rebuild_GRN = False,
                    GRN_file_dir = './2',
                    verbose = True)

xct2=xct2.Knk(targetgid)

XCTs = st.merge_scTenifoldXct(xct1, xct2)
emb = XCTs.get_embeds(train = True)

with h5py.File('merged_embeds.h5', 'w') as h5file:
    h5file.create_dataset('data', data=emb)

XCTs.nn_aligned_diff(emb) 
xcts_pairs_diff = XCTs.chi2_diff_test(pval = 1.0)
pd.DataFrame(xcts_pairs_diff).to_csv('output1.txt',index=False,header=True)

# if sys.argv[1]=="2":
if twosided==1:    
    xct1 = st.scTenifoldXct(data = adata1, 
                        source_celltype = 'Target',
                        target_celltype = 'Source',
                        obs_label = "cell_type",
                        rebuild_GRN = False,
                        GRN_file_dir = './1',
                        verbose = True)
    
    xct2 = st.scTenifoldXct(data = adata2, 
                        source_celltype = 'Target',
                        target_celltype = 'Source',
                        obs_label = "cell_type",
                        rebuild_GRN = False,
                        GRN_file_dir = './2',
                        verbose = True)

    xct2=xct2.Knk(targetgid)
    XCTs = st.merge_scTenifoldXct(xct1, xct2)
    emb = XCTs.get_embeds(train = True)
    XCTs.nn_aligned_diff(emb) 
    xcts_pairs_diff = XCTs.chi2_diff_test(pval = 1.0)
    pd.DataFrame(xcts_pairs_diff).to_csv('output2.txt',index=False,header=True)

