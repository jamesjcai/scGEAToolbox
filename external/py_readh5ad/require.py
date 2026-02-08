import sys
import subprocess
import importlib

# ------------------------------------------------------------------
# Required packages
# Mapping: import_name -> pip_name
# (only needed when import name differs from pip name)
# ------------------------------------------------------------------
REQUIRED = {
    "numpy": "numpy",
    "pandas": "pandas",
    "scipy": "scipy",
    "h5py": "h5py",
    "scanpy": "scanpy",
    "anndata": "anndata",
}

# ------------------------------------------------------------------
# Check which packages are missing (cannot be imported)
# ------------------------------------------------------------------
missing = []

for import_name, pip_name in REQUIRED.items():
    try:
        importlib.import_module(import_name)
    except ImportError:
        missing.append(pip_name)

# ------------------------------------------------------------------
# Install missing packages via pip
# ------------------------------------------------------------------
if missing:
    print(f"Installing missing packages: {sorted(missing)}")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", *sorted(missing)])
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to install some packages: {e}")
else:
    print("All required packages are already installed and importable.")




'''
import sys
import subprocess
import importlib.metadata

# Required packages
required = {'numpy', 'pandas', 'scipy', 'h5py', 'scanpy', 'anndata'} 

# Get installed package names (case-insensitive, safe against missing metadata)
installed = {
    name.lower()
    for dist in importlib.metadata.distributions()
    if (name := dist.metadata.get("Name")) is not None
}

# Find missing packages
missing = {pkg for pkg in required if pkg.lower() not in installed}

# Install missing packages if any
if missing:
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", *missing])
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to install some packages: {e}")

'''