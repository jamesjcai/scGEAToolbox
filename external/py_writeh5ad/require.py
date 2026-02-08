import sys
import subprocess
import importlib
import importlib.metadata as im


# ------------------------------------------------------------------
# Required packages (pip names)
# ------------------------------------------------------------------
REQUIRED = {"numpy", "pandas", "scipy", "h5py", "anndata"}


# ------------------------------------------------------------------
def get_installed_packages():
    """
    Retrieve installed package names using metadata.
    Falls back safely if metadata is incomplete.
    Works across Python 3.9–3.12.
    """
    installed = set()

    try:
        for dist in im.distributions():
            name = dist.metadata.get("Name")
            if name:
                installed.add(name.lower())
    except Exception:
        pass

    return installed


# ------------------------------------------------------------------
def verify_import(pkg_name):
    """
    Confirm module is actually importable.
    Prevents metadata-only false positives.
    """
    try:
        importlib.import_module(pkg_name)
        return True
    except Exception:
        return False


# ------------------------------------------------------------------
def find_missing_packages(required, installed):
    """
    Determine which required packages are missing or broken.
    """
    missing = set()

    for pkg in required:
        if pkg.lower() not in installed or not verify_import(pkg):
            missing.add(pkg)

    return missing


# ------------------------------------------------------------------
def install_packages(packages):
    """
    Install packages via pip using current Python executable.
    """
    if not packages:
        return

    print(f"Installing missing packages: {sorted(packages)}")

    subprocess.check_call([
        sys.executable,
        "-m",
        "pip",
        "install",
        *sorted(packages)
    ])


# ------------------------------------------------------------------
def main():

    installed = get_installed_packages()
    missing = find_missing_packages(REQUIRED, installed)

    if missing:
        try:
            install_packages(missing)
        except subprocess.CalledProcessError as e:
            print(f"❌ Failed to install packages: {e}")
            sys.exit(1)
    else:
        print("All required packages are installed.")


# ------------------------------------------------------------------
if __name__ == "__main__":
    main()


'''
import sys
import subprocess
import importlib.metadata

# Required packages
required = {"numpy", "pandas", "scipy", "h5py", "anndata"}

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