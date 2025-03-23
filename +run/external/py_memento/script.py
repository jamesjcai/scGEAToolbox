import scanpy as sc

import memento
#import matplotlib.pyplot as plt
import pandas as pd

adata = sc.read_h5ad("input.h5ad")
#print(adata)

#adata.obs[['BatchID','CellType']].sample(5)
#adata.obs['BatchID'] = adata.obs['BatchID'].apply(lambda x: 0 if x == '14w' else 1)
#adata.obs[['BatchID','CellType']].sample(5)


result_1d = memento.binary_test_1d(
    adata=adata, 
    capture_rate=0.07, 
    treatment_col='BatchID', 
    num_cpus=12,
    num_boot=5000)

# plt.scatter(result_1d.de_coef, result_1d.dv_coef, s=1)
# result_1d.query('de_coef > 0').sort_values('de_pval').head(10)
# result_1d.query('dv_coef > 0 & de_coef > 0').sort_values('dv_pval').head(10)

result_1d.query('de_coef > 0').sort_values('de_pval').to_csv('deoutput.txt', sep='\t', index=False)
result_1d.query('dv_coef > 0 & de_coef > 0').sort_values('dv_pval').to_csv('dvoutput.txt', sep='\t', index=False)

import itertools

gene_pairs = list(itertools.product(['IRF7'], adata.var.index.tolist()))

result_2d = memento.binary_test_2d(
    adata=adata, 
    gene_pairs=gene_pairs, 
    capture_rate=0.07, 
    treatment_col='BatchID', 
    num_cpus=12, 
    num_boot=5000)

result_2d.sort_values('corr_pval').to_csv('corr_output.txt', sep='\t', index=False)


# https://nbviewer.org/github/yelabucsf/scrna-parameter-estimation/blob/master/tutorials/binary_testing.ipynb