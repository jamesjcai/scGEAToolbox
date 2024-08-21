import os
#abspath = os.path.abspath(__file__)
#dname = os.path.dirname(abspath)
#os.chdir(dname)
#os.chdir("C:\\Users\\jcai.AUTH\\Documents\\GitHub\\scGEAToolbox\\+run\\thirdparty\\gseapy")
#os.chdir("U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\gseapy")
os.chdir("C:\\Users\\jcai\\Documents\\GitHub\\scGEAToolbox\\+run\\external\\py_GSEApy_enr")
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


# Load the packages
%matplotlib inline
%config InlineBackend.figure_format='retina' 
%load_ext autoreload
%autoreload 2
import csv
import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt

# Import gseapy library and visualise available databases
import gseapy
names = gseapy.get_library_name()
print(names)

# Import demo dataset
df = pd.read_csv('https://raw.githubusercontent.com/kuanrongchan/vaccine-studies/main/Ad5_seroneg.csv',index_col=0)

# Filter up and down-regulated DEGs
DEGs_up_1d = (df[(df['fc_1d'] > 1.5) & (df['qval_1d'] < 0.05)]).index.tolist()
DEGs_down_1d = (df[(df['fc_1d'] < -1.5) & (df['qval_1d'] < 0.05)]).index.tolist()

# Enrichr analysis of upregulated DEGs
enr_GOBP_up = gp.enrichr(gene_list=DEGs_up_1d ,
                 gene_sets=['GO_Biological_Process_2023'],
                 organism='Human', 
                 outdir='test/enr_DEGs_GOBP_up',
                 cutoff=0.5 
                )

enr_GOMF_up = gp.enrichr(gene_list=DEGs_up_1d ,
                 gene_sets=['GO_Molecular_Function_2021'],
                 organism='Human', 
                 description='DEGs_up_1d',
                 outdir='test/enr_DEGs_GOMF_up',
                 cutoff=0.5 
                )

enr_GOCC_up = gp.enrichr(gene_list=DEGs_up_1d ,
                 gene_sets=['GO_Cellular_Component_2021'],
                 organism='Human', 
                 description='DEGs_up_1d',
                 outdir='test/enr_DEGs_GOCC_up',
                 cutoff=0.5 
                )

enr_Reactome_up = gp.enrichr(gene_list=DEGs_up_1d ,
                 gene_sets=['Reactome_2016'],
                 organism='Human', 
                 description='DEGs_up_1d',
                 outdir='test/enr_DEGs_Reactome_up',
                 cutoff=0.5 
                )
                
# Show header columns of output files
enr_GOBP_up.results.head(5)

# Plot bargraphs for upregulated pathways
from gseapy.plot import barplot, dotplot
barplot(enr_GOBP_up.res2d,title='GO Biological Processes seroneg day 1 (up)',color = 'r')
barplot(enr_GOMF_up.res2d,title='GO Molecular Function seroneg day 1 (up)',color = 'r')
barplot(enr_GOCC_up.res2d,title='GO Cellular Component seroneg day 1 (up)',color = 'r')
barplot(enr_Reactome_up.res2d,title='Reactome seroneg day 1 (up)',color = 'r')

# Enrichr analysis of downregulated DEGs
enr_GOBP_down = gp.enrichr(gene_list=DEGs_down_1d ,
                 gene_sets=['GO_Biological_Process_2021'],
                 organism='Human', 
                 description='DEGs_down_1d',
                 outdir='test/enr_DEGs_GOBP_down',
                 cutoff=0.5 
                )

enr_GOMF_down = gp.enrichr(gene_list=DEGs_down_1d ,
                 gene_sets=['GO_Molecular_Function_2021'],
                 organism='Human', 
                 description='DEGs_down_1d',
                 outdir='test/enr_DEGs_GOMF_down',
                 cutoff=0.5 
                )

enr_GOCC_down = gp.enrichr(gene_list=DEGs_down_1d ,
                 gene_sets=['GO_Cellular_Component_2021'],
                 organism='Human', 
                 description='DEGs_down_1d',
                 outdir='test/enr_DEGs_GOCC_down',
                 cutoff=0.5 
                )

enr_Reactome_down = gp.enrichr(gene_list=DEGs_down_1d ,
                 gene_sets=['Reactome_2016'],
                 organism='Human', 
                 description='DEGs_down_1d',
                 outdir='test/enr_DEGs_Reactome_down',
                 cutoff=0.5 
                )
# Plot bargraphs for downregulated pathways
barplot(enr_GOBP_down.res2d,title='GO Biological Processes seroneg day 1 (down)',color = 'b')
barplot(enr_GOMF_down.res2d,title='GO Molecular Function seroneg day 1 (down)',color = 'b')
barplot(enr_GOCC_down.res2d,title='GO Cellular Component seroneg day 1 (down)',color = 'b')
barplot(enr_Reactome_down.res2d,title='Reactome seroneg day 1 (down)',color = 'b')