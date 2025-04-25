import sys
import subprocess
# import pkg_resources
import importlib.metadata

required  = {'numpy', 'pandas', 'scipy', 'h5py', 'anndata'} 
# installed = {pkg.key for pkg in pkg_resources.working_set}
installed = {distribution.metadata["Name"] for distribution in importlib.metadata.distributions()}
missing   = required - installed

if missing:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', *missing])


# pip install git+https://github.com/satijalab/panhumanpy.git

package_name = "panhumanpy"

try:
    __import__(package_name)
except ImportError:
    print(f"Package '{package_name}' not found. Installing...")
    subprocess.check_call([
        sys.executable,
        "-m",
        "pip",
        "install",
        "git+https://github.com/satijalab/panhumanpy.git"
    ])
else:
    print(f"Package '{package_name}' is already installed.")
