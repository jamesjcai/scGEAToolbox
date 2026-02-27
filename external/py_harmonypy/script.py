#import os
#abspath = os.path.abspath(__file__)
#dname = os.path.dirname(abspath)
#os.chdir(dname)

import pandas as pd
import harmonypy as hm
import h5py
from scipy.io import savemat


f = h5py.File('input.mat','r')
s = f['s'][()]
batchid= f['batchid'][()].astype(int).tolist()
f.close()

batchid=[str(x) for x in batchid]
meta_data=pd.DataFrame(batchid)
meta_data.columns =['batchidx']

#meta_data = pd.read_csv("input2.csv")
# data_mat = pd.read_csv("input1.csv", header=None)
#data_mat = pd.read_csv("input1.csv")
#data_mat = np.array(data_mat)

vars_use = ['batchidx']
ho = hm.run_harmony(s, meta_data, vars_use)

savemat("output.mat", {"sout": ho.Z_corr.T})

#res = pd.DataFrame(ho.Z_corr.T)
# res.columns = ['X{}'.format(i + 1) for i in range(res.shape[1])]
#res.to_csv("output.csv", sep = "\t", index = False, header=False)


