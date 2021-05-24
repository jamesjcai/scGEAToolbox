import pandas as pd
import numpy as np
from scipy.cluster.vq import kmeans
from scipy.stats.stats import pearsonr
import harmonypy as hm

meta_data = pd.read_csv("data/meta.csv")
data_mat = pd.read_csv("data/pcs.csv")
data_mat = np.array(data_mat)
vars_use = ['dataset']

ho = hm.run_harmony(data_mat, meta_data, vars_use)

# Write the adjusted PCs to a new file.
res = pd.DataFrame(ho.Z_corr)
res.columns = ['X{}'.format(i + 1) for i in range(res.shape[1])]
res.to_csv("adj.tsv", sep = "\t", index = False)

