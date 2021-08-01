import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
#os.chdir("U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\doubletdetection")

import pandas as pd
import scrublet as scr

counts=pd.read_csv("input.txt",header=None).values
scrub = scr.Scrublet(counts.T)
doublet_scores, predicted_doublets = scrub.scrub_doublets()

pd.DataFrame(predicted_doublets).to_csv('output1.txt',index=False,header=False)
pd.DataFrame(doublet_scores).to_csv('output2.txt',index=False,header=False)

# https://github.com/AllonKleinLab/scrublet