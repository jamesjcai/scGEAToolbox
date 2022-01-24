
pd = py.importlib.import_module('pandas');
np = py.importlib.import_module('numpy');
hm = py.importlib.import_module('harmonypy');
data_mat=pd.read_csv("input1.csv");
data_mat=np.array(data_mat);
meta_data = pd.read_csv("input2.csv");
% vars_use = py.list({py.str("batchidx")});
ho = hm.run_harmony(data_mat, meta_data, 'batchidx');
sout=np2mat(ho.Z_corr.T)

