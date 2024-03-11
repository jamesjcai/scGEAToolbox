#import os
#import sys
#abspath = os.path.abspath(__file__)
#dname = os.path.dirname(abspath)
#os.chdir(dname)
# os.chdir("D:\\GitHub\\scGEAToolbox\\+run\\external\\py_GenKI")
# os.chdir("C:\\Users\\jcai\\Documents\\GitHub\\scGEAToolbox\\+run\\external\\py_GenKI")


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
sc.settings.verbosity = 0
import GenKI as gk
from GenKI.preprocesing import build_adata
from GenKI.dataLoader import DataLoader
from GenKI.train import VGAE_trainer
from GenKI import utils
import h5py
import scipy as sp
from scipy.sparse import csr_matrix


def mat_to_npz(mat_path):
    f = h5py.File(mat_path, "r")
    # print(f.keys())
    mat = np.array(f.get(list(f.keys())[0]), dtype="float32")
    sparse_mat = scipy.sparse.csc_matrix(mat)
    return sparse_mat

# from scipy.io import savemat

# adata = build_adata("X.mat", "g.txt", "c.txt", delimiter=',', meta_cell_cols=['cell_type'], transpose=False)
adata = build_adata("X.mat",
                    "g.txt",
                    "c.txt",
                    meta_cell_cols=["cell_type"], # colname of cell type
                    transpose=False,
                    log_normalize=True,
                    scale_data=True)
#adata = adata[:100, :300].copy()
#print(adata)

f = h5py.File('idx.mat','r')
# counts = np.array(f.get('X'), dtype=np.float64)
# idx = int(f['idx'][()])-1
idx = f['idx'][0,0].astype(int)-1
f.close()

adata.var_names[idx]

# X=counts.T    # X = [ sparse or dense matrix, samples in rows, features in columns ]
#mdic = {"idx": np.array(sketch_index)}
#savemat("output.mat", mdic, oned_as='column')


#f = h5py.File("pcNet_Source.mat",'r')
#A = csr_matrix(f.get('A'))
#f.close()
#sp.sparse.save_npz("pcNet_Source.npz", A)

npz_mat = mat_to_npz("pcNet_Source.mat")
sp.sparse.save_npz("pcNet_Source.npz", npz_mat)

data_wrapper =  DataLoader(
                adata, # adata object
                target_gene = [adata.var_names[idx]], # KO gene name/index
                target_cell = None, # obsname for cell type, if none use all
                obs_label = "cell_type", # colname for genes
                GRN_file_dir = "./", # folder name for GRNs
                rebuild_GRN = False, # whether build GRN by pcNet
                pcNet_name = "pcNet_Source", # GRN file name
                verbose = True, # whether verbose
                n_cpus = 8, # multiprocessing
                )

data_wt = data_wrapper.load_data()
data_ko = data_wrapper.load_kodata()

# init trainer
hyperparams = {"epochs": 100, 
               "lr": 7e-4, 
               "beta": 1e-4, 
               "seed": 8096}
log_dir=None 

sensei = VGAE_trainer(data_wt, 
                     epochs=hyperparams["epochs"], 
                     lr=hyperparams["lr"], 
                     log_dir=log_dir, 
                     beta=hyperparams["beta"],
                     seed=hyperparams["seed"],
                     verbose=True,
                     )

sensei.train()

z_mu_wt, z_std_wt = sensei.get_latent_vars(data_wt)
z_mu_ko, z_std_ko = sensei.get_latent_vars(data_ko)
dis = gk.utils.get_distance(z_mu_ko, z_std_ko, z_mu_wt, z_std_wt, by="KL")
res_raw = utils.get_generank(data_wt, dis, rank=True)
# res_raw.head()
#null = sensei.pmt(data_ko, n=100, by="KL")
#res = utils.get_generank(data_wt, dis, null,)
# #                       save_significant_as = 'gene_list_example')
# https://github.com/yjgeno/GenKI/blob/master/notebook/Compatible_matlab.ipynb

df=pd.DataFrame(res_raw) 
df.to_csv('output.csv')