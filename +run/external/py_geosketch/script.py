import os
import sys
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
# os.chdir("D:\\GitHub\\scGEAToolbox\\+run\\external\\py_geosketch")
# os.chdir("C:\\Users\\jcai\\Documents\\GitHub\\scGEAToolbox\\+run\\external\\py_geosketch")

from geosketch import gs
import numpy as np
import h5py
from fbpca import pca
import pandas as pd

f = h5py.File('input.mat','r')

# counts = np.array(f.get(list(f.keys())[0]), dtype='float32')
# N = np.array(f.get(list(f.keys())[1]), dtype='float32')

counts = np.array(f.get('X'), dtype=np.float64)
n = int(f['n'][()])
#N = f.get('n')    # or f['n']
#n=N[()].astype(int).item()
f.close()
X=counts.T    # X = [ sparse or dense matrix, samples in rows, features in columns ]
sketch_index = gs(X, n, replace=False)

#U, s, Vt = pca(X, k=100) # E.g., 100 PCs.
#X_dimred = U[:, :100] * s[:100]
#sketch_index = gs(X_dimred, N, replace=False)
#X_sketch = X_dimred[sketch_index]

pd.DataFrame(sketch_index).to_csv('output.txt',index=False,header=False)
