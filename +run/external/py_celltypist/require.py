try:
    import numpy
    import h5py
    import scipy
    import celltypist
    print('All essential imports are found')
    exit(0)
except ImportError as exc:
    print(exc)
    exit(10)