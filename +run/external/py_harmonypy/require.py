"""
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
"""

import sys
import subprocess
# import pkg_resources
import importlib.metadata

required  = {'numpy', 'pandas', 'scipy', 'h5py', 'harmonypy'} 
# installed = {pkg.key for pkg in pkg_resources.working_set}
installed = {distribution.metadata["Name"] for distribution in importlib.metadata.distributions()}
missing   = required - installed

if missing:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', *missing])
