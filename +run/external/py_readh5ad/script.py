import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

import pandas as pd 
import scanpy as sc 
import scipy.io as spio

import sys
fname=str(sys.argv[1])
ada=sc.read_h5ad(fname)
ada.obs.to_csv("obs.csv")
ada.var.to_csv("var.csv") 
ada.write("X.h5")
