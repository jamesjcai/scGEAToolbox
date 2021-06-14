import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
# os.chdir("U:\\GitHub\\My_Code_Collection\\doubletdetec")

import pandas as pd
import numpy as np
import doubletdetection as dd

counts=pd.read_csv("input.txt").values
clf = dd.BoostClassifier(n_iters=2, use_phenograph=False, standard_scaling=True)
labels=clf.fit(counts.T).predict(p_thresh=1e-16, voter_thresh=0.5)
pd.DataFrame(clf.doublet_score().mask).to_csv('output.txt',index=False,header=False)
