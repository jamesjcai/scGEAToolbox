import sys
import subprocess
import importlib

# ------------------------------------------------------------------
# Required packages
# Mapping: import_name -> pip_name
# (pip_name differs for GitHub / custom installs)
# ------------------------------------------------------------------
REQUIRED = {
    "numpy": "numpy",
    "pandas": "pandas",
    "scipy": "scipy",
    "h5py": "h5py",
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

# ------------------------------------------------------------------
# Handle GitHub-only package separately
# ------------------------------------------------------------------
github_package_name = "panhumanpy"
github_install_cmd = "git+https://github.com/satijalab/panhumanpy.git"

try:
    importlib.import_module(github_package_name)
except ImportError:
    print(f"Package '{github_package_name}' not found. Installing from GitHub...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", github_install_cmd])
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to install '{github_package_name}': {e}")
else:
    print(f"Package '{github_package_name}' is already installed.")



'''
import sys
import subprocess
# import pkg_resources
import importlib.metadata

required  = {'numpy', 'pandas', 'scipy', 'h5py', 'anndata'} 
# installed = {pkg.key for pkg in pkg_resources.working_set}
installed = {distribution.name 
             for distribution in importlib.metadata.distributions()}


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
'''