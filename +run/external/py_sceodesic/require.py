"""
try:
    import sceodesic
    import h5py
    import numpy
    import pandas
    import scipy
    import scanpy
    import anndata
    import fbpca
    import scikit-learn
    import PyYaml
    print('All imports essential to the program are found')
    exit(0)
except ImportError as exc:
    print(exc)
    exit(10)
"""

import sys
import subprocess
#import pkg_resources
import importlib.metadata


required  = {'h5py', 'numpy', 'pandas', 'scipy', 'scanpy', 'anndata', 'fbpca', 'scikit-learn', 'PyYaml', 'sceodesic'} 
#installed = {pkg.key for pkg in pkg_resources.working_set}
installed = {distribution.metadata["Name"] for distribution in importlib.metadata.distributions()}

missing   = required - installed

if missing:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', *missing])
