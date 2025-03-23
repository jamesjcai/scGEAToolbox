import sys
import subprocess
# import pkg_resources
import importlib.metadata

required  = {'scanpy', 'h5py', 'anndata', 'memento'} 
# installed = {pkg.key for pkg in pkg_resources.working_set}
installed = {distribution.metadata["Name"] for distribution in importlib.metadata.distributions()}
missing   = required - installed

if missing:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', *missing])
