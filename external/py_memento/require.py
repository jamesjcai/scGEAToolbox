import sys
import subprocess
import importlib

# import_name → pip_name
required = {
    "scanpy": "scanpy",
    "h5py": "h5py",
    "anndata": "anndata",
    "memento": "memento-de",
}

missing = []

for import_name, pip_name in required.items():
    try:
        importlib.import_module(import_name)
    except ImportError:
        missing.append(pip_name)

if missing:
    subprocess.check_call([sys.executable, "-m", "pip", "install", *sorted(missing)])



'''
import sys
import subprocess
# import pkg_resources
import importlib.metadata

required  = {'scanpy', 'h5py', 'anndata', 'memento-de'} 
# installed = {pkg.key for pkg in pkg_resources.working_set}
installed = {distribution.name 
             for distribution in importlib.metadata.distributions()}


missing   = required - installed

if missing:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', *missing])
'''