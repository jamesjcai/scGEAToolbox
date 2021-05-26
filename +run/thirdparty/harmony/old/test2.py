import numpy as np
import pandas as pd


theta = None
lamb = None
sigma = 0.1
nclust = None
tau = 0
block_size = 0.05
epsilon_cluster = 1e-5
epsilon_harmony = 1e-4
plot_convergence = False
verbose = True
reference_values = None
cluster_prior = None
random_state = 0

N = meta_data.shape[0]

if data_mat.shape[1] != N:
    data_mat = data_mat.T

assert data_mat.shape[1] == N, \
    "data_mat and meta_data do not have the same number of cells" 

if nclust is None:
    nclust = np.min([np.round(N / 30.0), 100]).astype(int)

if type(sigma) is float and nclust > 1:
    sigma = np.repeat(sigma, nclust)

if isinstance(vars_use, str):
    vars_use = [vars_use]

phi = pd.get_dummies(meta_data[vars_use]).to_numpy().T
phi_n = meta_data[vars_use].describe().loc['unique'].to_numpy().astype(int)

if theta is None:
    theta = np.repeat([1] * len(phi_n), phi_n)
elif isinstance(theta, float) or isinstance(theta, int):
    theta = np.repeat([theta] * len(phi_n), phi_n)
elif len(theta) == len(phi_n):
    theta = np.repeat([theta], phi_n)

assert len(theta) == np.sum(phi_n), \
    "each batch variable must have a theta"

if lamb is None:
    lamb = np.repeat([1] * len(phi_n), phi_n)
elif isinstance(lamb, float) or isinstance(lamb, int):
    lamb = np.repeat([lamb] * len(phi_n), phi_n)
elif len(lamb) == len(phi_n):
    lamb = np.repeat([lamb], phi_n)

assert len(lamb) == np.sum(phi_n), \
    "each batch variable must have a lambda"

# Number of items in each category.
N_b = phi.sum(axis = 1)
# Proportion of items in each category.
Pr_b = N_b / N

if tau > 0:
    theta = theta * (1 - np.exp(-(N_b / (nclust * tau)) ** 2))

lamb_mat = np.diag(np.insert(lamb, 0, 0))

phi_moe = np.vstack((np.repeat(1, N), phi))

np.random.seed(random_state)


