try:
    import geosketch
    import os
    import pandas
    import numpy
    import h5py
    import scipy
    import fbpca
    print('All imports essential to geosketch.py are found')
#    exit(0)
except ImportError as exc:
    print(exc)
#    exit(10)