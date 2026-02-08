import sys
import subprocess
import importlib

required = ['h5py', 'anndata', 'cellbender']

missing = []

for pkg in required:
    try:
        importlib.import_module(pkg)
    except ImportError:
        missing.append(pkg)

if missing:
    subprocess.check_call([sys.executable, "-m", "pip", "install", *missing])


'''
import sys
import subprocess
# import pkg_resources
import importlib.metadata

required  = {'h5py', 'anndata', 'cellbender'} 
# installed = {pkg.lower() for pkg in im.packages_distributions()}
# installed = {pkg.key for pkg in pkg_resources.working_set}
# installed = {distribution.metadata["Name"] 
#             for distribution in importlib.metadata.distributions()}

installed = {distribution.name 
             for distribution in importlib.metadata.distributions()}

missing   = required - installed

if missing:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', *missing])
'''