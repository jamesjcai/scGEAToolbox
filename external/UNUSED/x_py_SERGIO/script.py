import os
import sys
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# os.chdir("D:\\GitHub\\scGEAToolbox\\+run\\external\\py_GenKI")
# os.chdir("C:\\Users\\jcai\\Documents\\GitHub\\scGEAToolbox\\+run\\external\\py_GenKI")

#f = h5py.File('idx.mat','r')
#counts = np.array(f.get('X'), dtype=np.float64)
#idx = int(f['idx'][()])-1
#f.close()


import numpy as np
from scipy.io import savemat
import h5py
from SERGIO.sergio import sergio

f = h5py.File('input.mat','r')
#ncells = int(f['ncells'][()])
#ngenes = int(f['ngenes'][()])
#N = f.get('ncells')    # or f['ncells']
#ncells = N[()].astype(int).item()
ncells = f['ncells'][0,0].astype(int)
ngenes = f['ngenes'][0,0].astype(int)
f.close()


sim = sergio(number_genes = ngenes, 
             number_bins = 1, 
             number_sc = ncells, 
             noise_params = 1, 
             decays=0.8, 
             sampling_state=15, 
             noise_type='dpd')

sim.build_graph(input_file_taregts ='targets.txt', input_file_regs='regs.txt', shared_coop_state=2)

sim.simulate()
expr = sim.getExpressions()
expr_clean = np.concatenate(expr, axis = 1)

expr_O = sim.outlier_effect(expr, outlier_prob = 0.01, mean = 0.8, scale = 1)
libFactor, expr_O_L = sim.lib_size_effect(expr_O, mean = 4.6, scale = 0.4)
binary_ind = sim.dropout_indicator(expr_O_L, shape = 6.5, percentile = 82)
expr_O_L_D = np.multiply(binary_ind, expr_O_L)

count_matrix = sim.convert_to_UMIcounts(expr_clean)
#X = np.concatenate(count_matrix, axis = 1)
savemat("output.mat", {"X": count_matrix})

#f = h5py.File("mytestfile.h5", "w")
#dset = f.create_dataset("X", c, dtype='i')
#f.close()
#df=pd.DataFrame(res_raw) 
#df.to_csv('output.csv')
