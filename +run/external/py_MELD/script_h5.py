import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
# os.chdir("U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\MELD")
# os.chdir("C:\\Users\\jcai.AUTH\\Documents\\GitHub\\scGEAToolbox\\+run\\thirdparty\\MELD")


import pandas as pd
import meld

# https://stackoverflow.com/questions/874461/read-mat-files-in-python
import numpy as np
import h5py
f = h5py.File('input.mat','r')
data = f.get('/X')
# data = np.array(data) # For converting to a NumPy array
sample_labels=f.get('/batchid')

#data=pd.read_csv("input.txt",header=None).values
#sample_labels=pd.read_csv("batchid.txt",header=None).values


meld_op = meld.MELD()
sample_densities = meld_op.fit_transform(data.T, sample_labels)
sample_likelihoods = meld.utils.normalize_densities(sample_densities)
pd.DataFrame(sample_likelihoods).to_csv('output.txt',index=False,header=True)

