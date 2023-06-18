try:
    import numpy
    import pandas
    import h5py
    import scipy
    import harmonypy
    print('All essential imports are found')
    exit(0)
except ImportError as exc:
    print(exc)
    exit(10)
