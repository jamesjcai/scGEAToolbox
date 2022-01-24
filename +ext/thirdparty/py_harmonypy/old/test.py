
import os
os.chdir("U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\harmony\\old")
import pandas as pd
import numpy as np
from scipy.cluster.vq import kmeans
from scipy.stats.stats import pearsonr
#import harmonypy as hm
from ccc import Harmony, run_harmony

meta_data = pd.read_csv("meta.csv")
data_mat = pd.read_csv("pcs.csv",header=None)
data_mat = np.array(data_mat)
vars_use = ['dataset']

ho = run_harmony(data_mat, meta_data, vars_use)

# Write the adjusted PCs to a new file.
res = pd.DataFrame(ho.Z_corr)
res.columns = ['X{}'.format(i + 1) for i in range(res.shape[1])]
res.to_csv("adj.tsv", sep = "\t", index = False)

