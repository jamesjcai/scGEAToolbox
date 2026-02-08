import sys
import subprocess
import importlib

# ------------------------------------------------------------------
# Required packages
# Mapping: import_name -> pip_name
# ------------------------------------------------------------------
REQUIRED = {
    "numpy": "numpy",
    "pandas": "pandas",
    "scipy": "scipy",
    "h5py": "h5py",
    "scrublet": "scrublet",
}

# ------------------------------------------------------------------
# Detect missing packages via import testing
# ------------------------------------------------------------------
missing = []

for import_name, pip_name in REQUIRED.items():
    try:
        importlib.import_module(import_name)
    except ImportError:
        missing.append(pip_name)

# ------------------------------------------------------------------
# Install missing packages
# ------------------------------------------------------------------
if missing:
    print(f"Installing missing packages: {sorted(missing)}")
    try:
        subprocess.check_call([
            sys.executable,
            "-m",
            "pip",
            "install",
            *sorted(missing)
        ])
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to install some packages: {e}")
else:
    print("All required packages are already installed and importable.")



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

import sys
import subprocess
# import pkg_resources
import importlib.metadata

required  = {'numpy', 'pandas', 'scipy', 'h5py', 'scrublet'} 
# installed = {pkg.key for pkg in pkg_resources.working_set}
installed = {distribution.name 
             for distribution in importlib.metadata.distributions()}

missing   = required - installed

if missing:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', *missing])

    # https://stackoverflow.com/questions/12332975/how-can-i-install-a-python-module-within-code    

"""


