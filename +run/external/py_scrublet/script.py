#import pandas as pd
import h5py
import scrublet as scr
from scipy.io import savemat


#counts=pd.read_csv("input.txt",header=None).values

f = h5py.File('input.mat','r')
X = f['X'][()]
f.close()
scrub = scr.Scrublet(X)
doubletscore, isDoublet = scrub.scrub_doublets()

savemat("output.mat", {"doubletscore": doubletscore, "isDoublet": isDoublet})
#pd.DataFrame(predicted_doublets).to_csv('output1.txt',index=False,header=False)
#pd.DataFrame(doublet_scores).to_csv('output2.txt',index=False,header=False)

# https://github.com/AllonKleinLab/scrublet