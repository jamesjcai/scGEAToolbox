"""Check for, and install, everything SCimilarity needs to annotate cells.

Run by run.py_scimilarity before each analysis. Two of SCimilarity's
dependencies cannot simply be pip-installed on Windows, so this does not just
shell out to ``pip install scimilarity``:

* ``hnswlib`` publishes no Windows wheel, so pip tries to compile it and fails
  on any machine without Microsoft Visual C++ Build Tools. install_hnswlib_win
  fetches the conda-forge build of upstream hnswlib instead. Do not substitute
  ``chroma-hnswlib``: it installs under the same import name but writes and
  reads a different file format, and rejects SCimilarity's index as corrupt.
* ``scimilarity`` pins ``zarr<3``, which would silently downgrade zarr and
  numcodecs for every other project sharing this interpreter. Nothing on the
  annotation path imports zarr -- it is used only by the training data models
  -- so SCimilarity goes in with --no-deps and its real dependencies are
  named here instead.

Output is deliberately ASCII: this runs through MATLAB's system(), whose
console is cp1252, and a stray non-ASCII byte aborts the run.
"""

import importlib
import importlib.metadata
import subprocess
import sys

import _scimilarity_env  # noqa: F401  (DLL search path, forgiving stdio)

IS_WINDOWS = sys.platform.startswith("win")

# import_name -> (pip_name, extra pip arguments)
REQUIRED = {
    "numpy": ("numpy", []),
    "pandas": ("pandas", []),
    "scipy": ("scipy", []),
    "h5py": ("h5py", []),
    "anndata": ("anndata", []),
    # Only used by a type check that _scimilarity_env patches for zarr 3;
    # any version will do, so do not let scimilarity's zarr<3 pin decide.
    "zarr": ("zarr", []),
    "torch": ("torch", []),
    "pytorch_lightning": ("pytorch-lightning", []),
    "captum": ("captum", []),
    "obonet": ("obonet", []),
    "circlify": ("circlify", []),
    "pyarrow": ("pyarrow", []),
    "tiledb": ("tiledb", []),
    "tiledb.cloud": ("tiledb-cloud", []),
    # Bundles a tiledb.dll of its own; letting pip resolve its dependencies
    # here would re-pin tiledb itself.
    "tiledb.vector_search": ("tiledb-vector-search", ["--no-deps"]),
    # On Windows this one is not installed by pip at all; see i_install.
    "hnswlib": ("hnswlib", []),
    "scimilarity": ("scimilarity", ["--no-deps"]),
}


def i_hnswlib_is_fork():
    """True if the importable hnswlib is chroma-hnswlib rather than upstream.

    The fork reads a different on-disk format, so an index written by
    SCimilarity fails to load. That surfaces much later as "Index seems to be
    corrupted or unsupported", which says nothing about the real cause -- so
    catch it here instead.
    """
    try:
        importlib.metadata.distribution("chroma-hnswlib")
        return True
    except importlib.metadata.PackageNotFoundError:
        return False


def i_install(pip_name, extra_args):
    """Install one package, returning True if pip reported success."""
    if pip_name == "hnswlib" and IS_WINDOWS:
        import install_hnswlib_win

        try:
            install_hnswlib_win.install_hnswlib()
            return True
        except Exception as err:  # noqa: BLE001 - reported, then given up on
            print("FAILED to install hnswlib: %s" % err)
            return False

    command = [sys.executable, "-m", "pip", "install"] + extra_args + [pip_name]
    print("Installing %s ..." % pip_name)
    try:
        subprocess.check_call(command)
        return True
    except subprocess.CalledProcessError as err:
        print("FAILED to install %s (pip exit status %d)." % (pip_name, err.returncode))
        return False


def i_can_import(import_name):
    """True if the module imports, which is the only check that means anything.

    A distribution can be present and still unimportable -- that is exactly the
    tiledb-vector-search DLL failure -- so this deliberately does not consult
    the installed-distribution list.
    """
    try:
        importlib.import_module(import_name)
        return True
    except Exception:
        # Not just ImportError: a broken native extension raises anything.
        return False


def i_needs_install(import_name):
    if import_name == "hnswlib" and i_hnswlib_is_fork():
        # Do not import it first: Windows will not let the installer
        # overwrite the fork's .pyd while this process has it loaded.
        return True
    return not i_can_import(import_name)


missing = [
    (import_name, pip_name, extra_args)
    for import_name, (pip_name, extra_args) in REQUIRED.items()
    if i_needs_install(import_name)
]

if not missing:
    print("All required packages are already installed and importable.")
    sys.exit(0)

print("Missing packages: %s" % ", ".join(name for name, _, _ in missing))

failed = []
for import_name, pip_name, extra_args in missing:
    if not i_install(pip_name, extra_args):
        failed.append(pip_name)
        continue
    # A fresh install lands in a directory importlib has already scanned,
    # and tiledb-vector-search's DLL folder did not exist a moment ago.
    importlib.invalidate_caches()
    _scimilarity_env.add_vector_search_dll_dir()
    if not i_can_import(import_name):
        print("Installed %s but 'import %s' still fails." % (pip_name, import_name))
        failed.append(pip_name)

if failed:
    print("Failed to install or import: %s" % ", ".join(failed))
    if "hnswlib" in failed and IS_WINDOWS:
        print(
            "hnswlib has no wheel on PyPI. Either re-run "
            "install_hnswlib_win.py to fetch the conda-forge build, or "
            "install Microsoft Visual C++ Build Tools so pip can compile it. "
            "Do not install chroma-hnswlib: its index format differs and "
            "SCimilarity's model will not load."
        )
    sys.exit(1)

print("All required packages are installed and importable.")
