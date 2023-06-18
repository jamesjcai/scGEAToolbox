try:
    import numpy
    import h5py
    import scipy
    import meld
    import pandas
    print('All essential imports are found')
    exit(0)
except ImportError as exc:
    print(exc)
    exit(10)