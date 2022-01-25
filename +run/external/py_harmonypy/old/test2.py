import os
os.chdir("U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\harmony\\old")
import pandas as pd
import numpy as np
from scipy.cluster.vq import kmeans
from scipy.stats.stats import pearsonr

meta_data = pd.read_csv("meta.csv")
data_mat = pd.read_csv("pcs.csv",header=None)
data_mat = np.array(data_mat)
vars_use = ['dataset']


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

tau = 0,
block_size = 0.05, 
max_iter_harmony = 10,
max_iter_kmeans = 20,
epsilon_cluster = 1e-5,
epsilon_harmony = 1e-4, 
plot_convergence = False,
verbose = True,
reference_values = None,
cluster_prior = None,
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


def safe_entropy(x: np.array):
    y = np.multiply(x, np.log(x))
    y[~np.isfinite(y)] = 0.0
    return y


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


#ho = Harmony(
#    data_mat, phi, phi_moe, Pr_b, sigma, theta, max_iter_harmony, max_iter_kmeans,
#    epsilon_cluster, epsilon_harmony, nclust, block_size, lamb_mat, verbose
#)
#self, Z, Phi, Phi_moe, Pr_b, sigma,
#theta, max_iter_harmony, max_iter_kmeans, 
#epsilon_kmeans, epsilon_harmony, K, block_size,
#lamb, verbose
#):

Z=data_mat
phi=phi
K=nclust
epsilon_kmeans=epsilon_cluster

Z_corr = np.array(Z)
Z_orig = np.array(Z)

Z_cos = Z_orig / Z_orig.max(axis=0)
Z_cos = Z_cos / np.linalg.norm(Z_cos, ord=2, axis=0)

km = kmeans(Z_cos.T, K, iter=10)
Y = km[0].T
Y = Y / np.linalg.norm(Y, ord=2, axis=0)

dist_mat = 2 * (1 - np.dot(Y.T, Z_cos))
R = -dist_mat
R = R / sigma[:,None]
R -= np.max(R, axis = 0)
R = np.exp(R)
R = R / np.sum(R, axis = 0)
# (3) Batch diversity statistics
E = np.outer(np.sum(R, axis=1), Pr_b)
O = np.inner(R , Phi)


kmeans_error = np.sum(np.multiply(R, dist_mat))
# Entropy
_entropy = np.sum(safe_entropy(R) * sigma[:,np.newaxis])
# Cross Entropy
x = (R * sigma[:,np.newaxis])
y = np.tile(theta[:,np.newaxis], K).T
z = np.log((O + 1) / (E + 1))
w = np.dot(y * z, Phi)
_cross_entropy = np.sum(x * w)
