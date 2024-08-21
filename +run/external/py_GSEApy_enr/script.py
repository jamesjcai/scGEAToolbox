import os
#abspath = os.path.abspath(__file__)
#dname = os.path.dirname(abspath)
#os.chdir(dname)
#os.chdir("C:\\Users\\jcai.AUTH\\Documents\\GitHub\\scGEAToolbox\\+run\\thirdparty\\gseapy")
#os.chdir("U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\gseapy")
os.chdir("C:\\Users\\jcai\\Documents\\GitHub\\scGEAToolbox\\+run\\external\\py_GSEApy_enr")

import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gseapy as gp

enr_out = gp.enrichr(gene_list="input.txt",
                 background="background.txt",    
                 gene_sets=['GO_Biological_Process_2023','GO_Molecular_Function_2023'],
                 organism='Human', 
                 outdir='./',
                 format='png',
                 cutoff=0.5
                )
