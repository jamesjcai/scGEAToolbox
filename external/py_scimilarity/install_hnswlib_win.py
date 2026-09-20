"""Install genuine hnswlib on Windows, where PyPI offers no wheel for it.

SCimilarity's annotation index (``annotation/labelled_kNN.bin``) is written by
upstream hnswlib, and only upstream hnswlib can read it back:

* PyPI ships hnswlib as source only. Without Microsoft Visual C++ Build Tools
  pip cannot compile it, which is the "Microsoft Visual C++ 14.0 or greater is
  required" failure.
* ``chroma-hnswlib``, the one prebuilt Windows package installed under the
  same import name, is a fork with a *different file format*: it stores a
  float per element between the level-0 block and the link lists, so it reads
  an upstream index four bytes per element out of alignment and rejects it
  with "Index seems to be corrupted or unsupported". It is not a substitute.

conda-forge builds upstream hnswlib for Windows with MSVC, and its package is
an ordinary pip install laid out under ``Lib/site-packages`` -- .pyd plus a
valid .dist-info. This lifts that directory into the running interpreter's
site-packages, so pip sees a normal hnswlib 0.8.0 afterwards and can uninstall
it again. No conda installation is needed.

The alternative, if you would rather not fetch from conda-forge, is to install
Microsoft Visual C++ Build Tools and let pip compile hnswlib itself:

    winget install Microsoft.VisualStudio.2022.BuildTools --override ^
        "--wait --quiet --add Microsoft.VisualStudio.Workload.VCTools --includeRecommended"
    pip install hnswlib

Run this module directly to install, or call install_hnswlib().
"""

import hashlib
import io
import json
import os
import shutil
import subprocess
import sys
import sysconfig
import tarfile
import tempfile
import urllib.request
import zipfile

CHANNEL_API = "https://api.anaconda.org/package/conda-forge/hnswlib/files"
SUBDIR = "win-64"


def i_python_tag():
    return "py%d%d" % (sys.version_info[0], sys.version_info[1])


def i_sort_key(entry):
    """Order candidates by version then build number, newest last."""
    version = tuple(
        int(part) if part.isdigit() else -1 for part in entry["version"].split(".")
    )
    return (version, entry["attrs"].get("build_number", 0))


def i_choose_build(listing):
    tag = i_python_tag()
    candidates = [
        entry
        for entry in listing
        if entry["basename"].startswith(SUBDIR + "/")
        and tag in entry["attrs"].get("build", "")
        and entry["basename"].endswith(".conda")
    ]
    if not candidates:
        raise RuntimeError(
            "conda-forge has no hnswlib build for %s on %s." % (tag, SUBDIR)
        )
    return max(candidates, key=i_sort_key)


def i_download(entry, destdir):
    url = entry["download_url"]
    if url.startswith("//"):
        url = "https:" + url
    path = os.path.join(destdir, os.path.basename(entry["basename"]))
    print("Downloading %s ..." % entry["basename"])
    with urllib.request.urlopen(url) as response, open(path, "wb") as fh:
        shutil.copyfileobj(response, fh)

    digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
    if digest != entry["sha256"]:
        raise RuntimeError(
            "Checksum mismatch for %s: expected %s, got %s."
            % (entry["basename"], entry["sha256"], digest)
        )
    return path


def i_read_payload(conda_path):
    """Return the package payload of a .conda (a zip of zstd-compressed tars)."""
    try:
        import zstandard
    except ImportError:
        print("Installing zstandard (needed to unpack a .conda file) ...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "zstandard"])
        import zstandard

    with zipfile.ZipFile(conda_path) as archive:
        name = next(n for n in archive.namelist() if n.startswith("pkg-"))
        with archive.open(name) as compressed:
            return zstandard.ZstdDecompressor().stream_reader(compressed).read()


def i_extract_site_packages(payload, destdir):
    prefix = "Lib/site-packages/"
    written = []
    with tarfile.open(fileobj=io.BytesIO(payload)) as tar:
        for member in tar.getmembers():
            if not member.isfile() or not member.name.startswith(prefix):
                continue
            target = os.path.join(destdir, member.name[len(prefix):].replace("/", os.sep))
            os.makedirs(os.path.dirname(target), exist_ok=True)
            source = tar.extractfile(member)
            with open(target, "wb") as fh:
                shutil.copyfileobj(source, fh)
            written.append(member.name[len(prefix):])
    if not any(name.endswith(".pyd") for name in written):
        raise RuntimeError("conda package contained no extension module.")
    return written


def install_hnswlib():
    """Install upstream hnswlib into this interpreter. Returns its version."""
    if not sys.platform.startswith("win"):
        raise RuntimeError("This installer is for Windows; elsewhere use pip.")

    # chroma-hnswlib owns the same import name, so it has to go first.
    subprocess.call(
        [sys.executable, "-m", "pip", "uninstall", "-y", "chroma-hnswlib"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    with urllib.request.urlopen(CHANNEL_API) as response:
        listing = json.load(response)
    entry = i_choose_build(listing)
    print("Selected conda-forge %s" % entry["basename"])

    destination = sysconfig.get_paths()["platlib"]
    with tempfile.TemporaryDirectory() as workdir:
        payload = i_read_payload(i_download(entry, workdir))
        written = i_extract_site_packages(payload, destination)
    print("Installed %d file(s) into %s" % (len(written), destination))
    return entry["version"]


if __name__ == "__main__":
    try:
        version = install_hnswlib()
    except Exception as err:  # noqa: BLE001 - the message is the whole point
        print("Could not install hnswlib: %s" % err)
        sys.exit(1)
    print("hnswlib %s installed." % version)
