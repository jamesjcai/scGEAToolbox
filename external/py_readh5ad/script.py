import os
import scanpy as sc
import sys

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
fname=str(sys.argv[1])
ada=sc.read_h5ad(fname)
ada.obs.to_csv("obs.csv")
ada.var.to_csv("var.csv") 
ada.write("X.h5")
