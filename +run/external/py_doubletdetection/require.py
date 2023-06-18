try:
    import numpy
    import h5py
    import scipy
    import doubletdetection
    print('All imports essential to doubletdetection are found')
    exit(0)
except ImportError as exc:
    print(exc)
    exit(10)