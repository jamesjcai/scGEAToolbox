try:
    from scipy.io import savemat
    import doubletdetection
    import geosketch
    import numpy
    import h5py
    import scipy
    print('All imports essential to the program are found')
#    exit(0)
except ImportError as exc:
    print(exc)
#    exit(10)