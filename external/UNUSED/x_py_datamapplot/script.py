#import os
#import sys
#abspath = os.path.abspath(__file__)
#dname = os.path.dirname(abspath)
# os.chdir(dname)
# os.chdir("D:\\GitHub\\scGEAToolbox\\+run\\external\\py_datamapplot")
# os.chdir("C:\\Users\\jcai\\Documents\\GitHub\\scGEAToolbox\\+run\\external\\py_geosketch")

import matplotlib
# matplotlib.rcParams["figure.dpi"] = 72
import matplotlib.pyplot as plt
import datamapplot
import h5py
import numpy as np

f = h5py.File('input.mat','r')
data_map = np.array(f.get('s'), dtype=np.float64)
f.close()

with open('c.txt', 'r') as file:
    content = file.read().splitlines()
data_label = np.array(content)
data_map = np.transpose(data_map)
datamapplot.create_plot(data_map, data_label, label_font_size=14, point_size=10)
plt.show()
