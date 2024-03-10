try:
    import GenKI
    import os
    import pandas
    import numpy
    import matplotlib.pyplot
    import scipy
    import h5py
    import torch_geometric
    import tensorboard
    print('All imports essential to GenKI are found')
    exit(0)
except ImportError as exc:
    print(exc)
    exit(10)