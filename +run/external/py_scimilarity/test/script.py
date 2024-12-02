import os
# os.chdir("C:\\Users\\jcai\\AppData\\Local\\Temp\\")
os.chdir("Z:\\Cailab\\GitHub\\scGEAToolbox\\+run\\external\\py_scimilarity\\test")
#import os
#os.chdir("./")
import pandas as pd 
#import scanpy as sc 
#import scipy.io as spio
#import numpy as np
from scipy.sparse import csr_matrix
import h5py
import anndata




# ==============================


import scanpy as sc
from matplotlib import pyplot as plt

sc.set_figure_params(dpi=100)
plt.rcParams["figure.figsize"] = [6, 4]

import warnings

warnings.filterwarnings("ignore")

from scimilarity.utils import lognorm_counts, align_dataset
from scimilarity import CellAnnotation

# Instantiate the CellAnnotation object
# Set model_path to the location of the uncompressed model

# model_path = "/models/model_v1.1"
model_path = "D:\\downloads\\shetty_09_24\\scimlarity\\models\\model_v1.1"
ca = CellAnnotation(model_path=model_path)

# Load the tutorial data
# Set data_path to the location of the tutorial dataset
data_path = "input.h5ad"
adams = sc.read(data_path)

adams = align_dataset(adams, ca.gene_order)

adams = lognorm_counts(adams)

adams.obsm["X_scimilarity"] = ca.get_embeddings(adams.X)

sc.pp.neighbors(adams, use_rep="X_scimilarity")
sc.tl.umap(adams)

predictions, nn_idxs, nn_dists, nn_stats = ca.get_predictions_knn(
    adams.obsm["X_scimilarity"]
)
adams.obs["predictions_unconstrained"] = predictions.values

celltype_counts = adams.obs.predictions_unconstrained.value_counts()
well_represented_celltypes = celltype_counts[celltype_counts > 20].index

#sc.pl.umap(
#    adams[adams.obs.predictions_unconstrained.isin(well_represented_celltypes)],
#    color="predictions_unconstrained",
#    legend_fontsize=15,
#)
adams.write("output.h5ad")
