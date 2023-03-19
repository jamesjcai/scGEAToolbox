import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
os.chdir("D:\\GitHub\\scGEAToolbox\\+run\\external\\py_metaviz")
os.chdir("C:\\Users\\jcai\\Documents\\GitHub\\scGEAToolbox\\+run\\external\\py_metaviz")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import umap
from sklearn import manifold, datasets

# load S curve data
n_samples = 1000
S_points, S_color = datasets.make_s_curve(n_samples)

lle_standard = manifold.LocallyLinearEmbedding(method="standard")
S_standard = lle_standard.fit_transform(S_points)

md_scaling = manifold.MDS(
    n_components=2, max_iter=50, n_init=4
)
S_scaling = md_scaling.fit_transform(S_points)

spectral = manifold.SpectralEmbedding(
    n_components=2, n_neighbors=40
)
S_spectral = spectral.fit_transform(S_points)

t_sne = manifold.TSNE(
    n_components=2,
    learning_rate="auto",
    perplexity=30,
    n_iter=250,
    init="random",
)
S_t_sne = t_sne.fit_transform(S_points)

reducer = umap.UMAP(n_components=2)
S_umap = reducer.fit_transform(S_points)


pd.DataFrame(S_umap).to_csv('1.txt',index=False,header=False)
pd.DataFrame(S_standard).to_csv('2.txt',index=False,header=False)
pd.DataFrame(S_scaling).to_csv('3.txt',index=False,header=False)
pd.DataFrame(S_spectral).to_csv('4.txt',index=False,header=False)
pd.DataFrame(S_t_sne).to_csv('5.txt',index=False,header=False)

pd.DataFrame(S_color).to_csv('c.txt',index=False,header=False)


from meta_visualization import meta_viz




meta_data = pd.read_csv("1.txt")
