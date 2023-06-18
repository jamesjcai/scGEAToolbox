import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
# os.chdir("U:\\GitHub\\My_Code_Collection\\doubletdetec")

import h5py
from scipy.io import savemat

import numpy as np
import doubletdetection as dd

f = h5py.File('input.mat','r')
counts = np.array(f.get('X'), dtype=np.float64)

#import pandas as pd
#counts=pd.read_csv("input.txt",header=None).values
clf = dd.BoostClassifier(n_iters=2, use_phenograph=False, standard_scaling=True)
labels=clf.fit(counts.T).predict(p_thresh=1e-16, voter_thresh=0.5)
doublet_score=clf.doublet_score()

savemat("output.mat", {"X": doublet_score.mask, "score": doublet_score})
#savemat("output2.mat", {"score": doublet_score})
#pd.DataFrame(doublet_score.mask).to_csv('output1.txt',index=False,header=False)
#pd.DataFrame(doublet_score).to_csv('output2.txt',index=False,header=False)

# https://github.com/JonathanShor/DoubletDetection/blob/master/tests/notebooks/PBMC_10k_vignette.ipynb