try:
    import scTenifoldXct
    print('All essential imports are found')
    exit(0)
except ImportError as exc:
    print(exc)
    exit(10)


"""
import sys
import subprocess
import pkg_resources
import importlib.metadata

required  = {'numpy', 'pandas', 'scipy', 'h5py', 'scanpy'} 
installed = {pkg.key for pkg in pkg_resources.working_set}
missing   = required - installed

if missing:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', *missing])
"""
