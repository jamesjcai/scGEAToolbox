"""
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
"""


#import subprocess
#import sys
#
#try:
#    import pandas as pd
#except ImportError:
#    subprocess.check_call([sys.executable, "-m", "pip", "install", 'pandas'])
#finally:
#    import pandas as pd

# https://stackoverflow.com/questions/12332975/how-can-i-install-a-python-module-within-code

import sys
import subprocess
# import pkg_resources
import importlib.metadata

required  = {'numpy', 'pandas', 'scipy', 'h5py', 'meld'} 
# installed = {pkg.key for pkg in pkg_resources.working_set}
installed = {distribution.metadata["Name"] for distribution in importlib.metadata.distributions()}
missing   = required - installed

if missing:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', *missing])
