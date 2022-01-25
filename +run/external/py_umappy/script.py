import os
#import sys
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# os.chdir("U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\harmony")
import pandas as pd
import umap as um

data_mat = pd.read_csv("input.csv",header=None)
#embedding= um.UMAP(densmap=True).fit_transform(data_mat)
embedding= um.UMAP(densmap=False).fit_transform(data_mat)
res = pd.DataFrame(embedding)
res.to_csv("output.csv", sep = "\t", index = False, header=False)
