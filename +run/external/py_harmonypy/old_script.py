import os
#import sys
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# os.chdir("U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\harmony")
import pandas as pd
import numpy as np
from scipy.cluster.vq import kmeans
from scipy.stats.stats import pearsonr
import harmonypy as hm

meta_data = pd.read_csv("input2.csv")
# data_mat = pd.read_csv("input1.csv", header=None)
data_mat = pd.read_csv("input1.csv")

data_mat = np.array(data_mat)
vars_use = ['batchidx']
ho = hm.run_harmony(data_mat, meta_data, vars_use)

res = pd.DataFrame(ho.Z_corr.T)
# res.columns = ['X{}'.format(i + 1) for i in range(res.shape[1])]
res.to_csv("output.csv", sep = "\t", index = False, header=False)


