import sys
import os
os.chdir("C:/Users/jcai/Downloads/pbmc_preprocessed.h5mu")

from mowgli import models
import muon as mu
import scanpy as sc
from muon import atac as ac
from muon import MuData



#mdata = MuData({'rna': adata_rna, 'atac': adata_atac})
#adata_rna = sc.read_h5ad("pbmc_granulocyte_sorted_10k_gex_molecule_info.h5")

mdata = mu.read_10x_h5("03_filtered_feature_bc_matrix.h5")
print(mdata.mod_names)

mdata.write("test.h5mu")
mdatax = mu.read("test.h5mu")

adata = mu.read("test.h5mu/rna")
bdata = mu.read("test.h5mu/atac")
# mu.write("test.h5mu/rna", adata)

#sc.tl.pca(adata)
#mu.tl.mofa(mdata)


sc.pp.normalize_total(mdata.mod["rna"])
sc.pp.log1p(mdata.mod["rna"])
sc.pp.normalize_total(mdata.mod["atac"])
sc.pp.log1p(mdata.mod["atac"])
#sc.pp.highly_variable_genes(mdata.mod["rna"])

sc.pp.highly_variable_genes(mdata["rna"], n_top_genes=500)
sc.pp.highly_variable_genes(mdata["atac"], n_top_genes=500)


sc.pp.neighbors(mdata.mod["rna"])
sc.tl.umap(mdata.mod["rna"])

#rna = mdata.mod["rna"]
#sc.pl.umap(rna)
#ac.pp.tfidf(bdata)



# Initialize and train the model.
model = models.MowgliModel(latent_dim=15)
model.train(mdata)
# Visualize the embedding with UMAP.
sc.pp.neighbors(mdata)
sc.tl.umap(mdata)
sc.pl.umap(mdata)


