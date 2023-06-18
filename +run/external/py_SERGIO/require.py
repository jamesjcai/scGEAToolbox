try:
    import numpy
    import h5py
    import scipy
    from SERGIO.sergio import sergio
    print('All imports essential to SERGIO are found')
    exit(0)
except ImportError as exc:
    print(exc)
    exit(10)