"""
try:
    import h5py
    import scipy
    import scrublet
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

required  = {'numpy', 'pandas', 'scipy', 'h5py', 'scrublet'} 
# installed = {pkg.key for pkg in pkg_resources.working_set}
installed = {distribution.metadata["Name"] for distribution in importlib.metadata.distributions()}
missing   = required - installed

if missing:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', *missing])

    # https://stackoverflow.com/questions/12332975/how-can-i-install-a-python-module-within-code    