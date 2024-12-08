import pandas as pd
from scipy.sparse import csr_matrix
import h5py
import anndata
import os.path


f = h5py.File("X.mat",'r')
# counts = csr_matrix(f.get('X'))
Xnorm = csr_matrix(f.get('Xnorm'))
modeldir = f['modeldir'][()]
modeldir = modeldir.tobytes().decode('utf-16')
g = f['g']
gene_names = []
for r in g:
    for ref in r:
        str_data = ''.join(chr(c[0]) for c in f[ref][:])
        gene_names.append(str_data)

target_celltypes = []        
if 'tg' in f:        
    g = f['tg']
    for r in g:
        for ref in r:
            str_data = ''.join(chr(c[0]) for c in f[ref][:])
            target_celltypes.append(str_data)
f.close()

adata = anndata.AnnData(X=Xnorm)
#adata.layers["counts"] = counts
#metadata = pd.read_csv("c.csv")
#with open("g.csv",'r') as f:
#          gene_names = f.read().splitlines()
#adata.obs = metadata
#adata.obs.index = adata.obs['CellID'].tolist()
adata.obs.index = [f"cell_{i}" for i in range(adata.n_obs)]
adata.var.index = gene_names
print("Input data ready.")

# import scanpy as sc
from scimilarity.utils import lognorm_counts, align_dataset
from scimilarity import CellAnnotation

model_path = modeldir
ca = CellAnnotation(model_path=model_path)
print("Model read.")


#if os.path.exists("tg.csv"):    
#    with open("tg.csv",'r') as f:
#          target_celltypes = f.read().splitlines()
#    ca.safelist_celltypes(target_celltypes)
#    print("Constrained annotation.")
#else:
#    print("Unconstrained annotation.")

if target_celltypes != []:
    ca.safelist_celltypes(target_celltypes)

#target_celltypes = [
#    "glutamatergic neuron",
#    "microglial cell",
#    "radial glial cell",
#    "neuron",
#    "astrocyte",
#    "glial cell",
#    "progenitor cell",
#    "oligodendrocyte"
#]
#ca.safelist_celltypes(target_celltypes)
#print("Model constraining...done.")

adata = align_dataset(adata, ca.gene_order)
#adata = lognorm_counts(adata)
embeddings = ca.get_embeddings(adata.X)
print("Get_embeddings...done.")

#predictions = ca.get_predictions_knn(embeddings)
print("Prediction...done.")

#predictions.to_csv("output.csv", index=True)
