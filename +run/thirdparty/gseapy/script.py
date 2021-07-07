import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
#os.chdir("C:\\Users\\jcai.AUTH\\Documents\\GitHub\\scGEAToolbox\\+run\\thirdparty\\gseapy")
os.chdir("U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\gseapy")
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt


rnk = pd.read_csv("input.txt", sep=",")
# rnk.head()
pre_res = gp.prerank(rnk=rnk, gene_sets='KEGG_2016',
                     processes=4,
                     permutation_num=100, # reduce number to speed up testing
                     outdir='res', format='png', seed=6, no_plot=True)
# pre_res.res2d.sort_index().head()

# https://github.com/zqfang/GSEApy/blob/f3b9a6a9e80b4503196cd74108e40997d01d4ae8/gseapy/utils.py
# 'GO_Biological_Process_2015',
