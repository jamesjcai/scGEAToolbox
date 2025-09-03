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

