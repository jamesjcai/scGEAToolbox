import os
#abspath = os.path.abspath(__file__)
#dname = os.path.dirname(abspath)
#os.chdir(dname)
#os.chdir("D:\\GitHub\\scGEAToolbox\\+run\\external\\py_metaviz")
os.chdir("C:\\Users\\jcai\\Documents\\GitHub\\scGEAToolbox\\+run\\external\\py_metaviz")
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import umap
from sklearn import manifold, datasets

S_color = pd.read_csv("c.txt", header=None)
X1 = pd.read_csv("1.txt", header=None)
X2 = pd.read_csv("2.txt", header=None)
X3 = pd.read_csv("3.txt", header=None)
X4 = pd.read_csv("4.txt", header=None)
X5 = pd.read_csv("5.txt", header=None)

from meta_visualization import meta_viz

meta_distances = meta_viz([X1, X2, X3, X4, X5])

print(meta_distances[0,1])
print(meta_distances[1,0])
print(meta_distances[0,2])
print(meta_distances[2,0])