"""Make a stock CPython able to run SCimilarity.

Three things have to happen around ``import scimilarity``, and none of them is
SCimilarity's doing. The first two are Windows-only:

* ``tiledb-vector-search`` ships its own ``tiledb.dll`` in
  ``tiledb/vector_search/lib`` but never puts that folder on the DLL search
  path, so importing it dies with "DLL load failed while importing
  _tiledbvspy". ``CellAnnotation`` imports it unconditionally, even for the
  hnswlib code path, so this is not optional.
* The console in a MATLAB ``system()`` call is cp1252, so any non-ASCII byte a
  library prints raises UnicodeEncodeError and takes the run down with it.
* ``get_embeddings`` type-checks against ``zarr.core.Array``, which zarr 3
  renamed to ``zarr.Array``. SCimilarity therefore pins ``zarr<3``, but the
  check only guards zarr-backed sparse matrices -- irrelevant to anything this
  toolbox passes it -- so restoring the old name is enough, and spares every
  other project in this interpreter a zarr downgrade.

Import this module first; it is a no-op off Windows. Call
``add_vector_search_dll_dir()`` again after installing tiledb-vector-search,
since the folder did not exist when this module was imported.
"""

import os
import sys


def add_vector_search_dll_dir():
    if not sys.platform.startswith("win"):
        return
    if not hasattr(os, "add_dll_directory"):
        return
    for entry in sys.path:
        libdir = os.path.join(entry, "tiledb", "vector_search", "lib")
        if os.path.isdir(libdir):
            try:
                os.add_dll_directory(libdir)
            except OSError:
                pass
            return


def patch_zarr_core_array():
    """Put ``zarr.core.Array`` back for code written against zarr 2."""
    try:
        import zarr
        import zarr.core
    except ImportError:
        return  # require.py installs it; nothing to patch yet
    if not hasattr(zarr.core, "Array") and hasattr(zarr, "Array"):
        zarr.core.Array = zarr.Array


def _make_stdio_forgiving():
    # Replace characters the console cannot encode rather than raising.
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="replace")
        except (AttributeError, ValueError):
            pass


add_vector_search_dll_dir()
patch_zarr_core_array()
_make_stdio_forgiving()
