import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
# os.chdir("U:\\GitHub\\scGEAToolbox\\+run\\external\\py_scTenifoldXct")

import scanpy as sc
import scTenifoldXct as st
from scTenifoldXct.dataLoader import build_adata
import pandas as pd

adata = build_adata("X.txt", "g.txt", "c.txt", delimiter=',', meta_cell_cols=['cell_type'])
# adata.var.columns = ['gene_name']
# adata.obs.columns = ['cell_type']

xct = st.scTenifoldXct(data = adata, 
                    source_celltype = 'Source',
                    target_celltype = 'Target',
                    obs_label = "cell_type",
                    rebuild_GRN = True,
                    GRN_file_dir = 'temp_net',
                    verbose = True,
                    n_cpus = -1)
emb = xct.get_embeds(train = True)
xct_pairs = xct.null_test()
print(xct_pairs)
pd.DataFrame(xct_pairs).to_csv('output1.txt',index=False,header=True)


xct = st.scTenifoldXct(data = adata, 
                    source_celltype = 'Target',
                    target_celltype = 'Source',
                    obs_label = "cell_type",
                    rebuild_GRN = False,
                    GRN_file_dir = 'temp_net',
                    verbose = True,
                    n_cpus = -1)
emb = xct.get_embeds(train = True)
xct_pairs = xct.null_test()
print(xct_pairs)
pd.DataFrame(xct_pairs).to_csv('output2.txt',index=False,header=True)

