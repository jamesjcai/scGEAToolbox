import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
# os.chdir("U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\MELD")
# os.chdir("C:\\Users\\jcai.AUTH\\Documents\\GitHub\\scGEAToolbox\\+run\\thirdparty\\MELD")
# https://www.blogforbrains.com/blog/2014/9/6/loading-matlab-mat-data-in-python
import scipy.io as spio
import pandas as pd
import meld


mat = spio.loadmat('input.mat', squeeze_me=True)
data = mat['X'] # array
sample_labels=mat['batchid']

meld_op = meld.MELD()
sample_densities = meld_op.fit_transform(data.T, sample_labels)
sample_likelihoods = meld.utils.normalize_densities(sample_densities)
pd.DataFrame(sample_likelihoods).to_csv('output.txt',index=False,header=True)
# spio.savemat('output.mat',{"score": sample_likelihoods})

