import sys
import subprocess
# import pkg_resources
import importlib.metadata

required  = {'h5py', 'anndata', 'cellbender'} 
# installed = {pkg.key for pkg in pkg_resources.working_set}
# installed = {distribution.metadata["Name"] 
#             for distribution in importlib.metadata.distributions()}

installed = {distribution.name 
             for distribution in importlib.metadata.distributions()}

missing   = required - installed

if missing:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', *missing])
