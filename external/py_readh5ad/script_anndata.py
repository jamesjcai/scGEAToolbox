import os
import anndata as ad
import sys

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
fname=str(sys.argv[1])


adata = ad.read(fname)
# adata.X
adata.obs_names = [f"Cell_{i}" for i in range(adata.n_obs)]
adata.var_names = [f"Gene_{i}" for i in range(adata.n_vars)]
print(adata.obs_names[:10])



#ada=sc.read_h5ad(fname)
#ada.obs.to_csv("obs.csv")
#ada.var.to_csv("var.csv") 
#ada.write("X.h5")
