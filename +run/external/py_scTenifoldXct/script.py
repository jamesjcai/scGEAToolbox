import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

import scanpy as sc
import scTenifoldXct as st

adata = sc.read_h5ad('C:\\Users\\jcai\\adata_short_example.h5ad')
xct = st.scTenifoldXct(data = adata, 
                    cell_names = ['Inflam. FIB', 'Inflam. DC'],
                    obs_label = "ident",
                    rebuild_GRN = True, 
                    GRN_file_dir = './Net_example_dev',  
                    verbose = True,
                    n_cpus = -1)
emb = xct.get_embeds(train = True)
xct_pairs = xct.null_test()
print(xct_pairs)



