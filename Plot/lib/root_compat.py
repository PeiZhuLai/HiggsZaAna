"""Small PyROOT import guard shared by plotting scripts."""

from __future__ import annotations

import importlib
import os
import sys


def _prepend_rootsys_python_paths() -> bool:
    rootsys = os.environ.get("ROOTSYS")
    if not rootsys:
        return False

    changed = False
    for path in (os.path.join(rootsys, "lib"), os.path.join(rootsys, "lib64")):
        if os.path.isdir(path) and path not in sys.path:
            sys.path.insert(0, path)
            changed = True
    return changed


def _load_root_module():
    try:
        return importlib.import_module("ROOT")
    except Exception as exc:
        if not _prepend_rootsys_python_paths():
            raise RuntimeError(
                "Failed to import CERN PyROOT. Activate a ROOT-enabled environment "
                "before running this script."
            ) from exc
        try:
            return importlib.import_module("ROOT")
        except Exception as retry_exc:
            raise RuntimeError(
                "Failed to import CERN PyROOT even after adding ROOTSYS library "
                "paths to PYTHONPATH."
            ) from retry_exc


def import_pyroot():
    ROOT = _load_root_module()

    if hasattr(ROOT, "gROOT"):
        return ROOT

    if _prepend_rootsys_python_paths():
        sys.modules.pop("ROOT", None)
        ROOT = _load_root_module()
        if hasattr(ROOT, "gROOT"):
            return ROOT

    module_file = getattr(ROOT, "__file__", None)
    module_path = getattr(ROOT, "__path__", None)
    location = module_file or module_path or "<unknown>"
    raise RuntimeError(
        "Imported a Python module named ROOT, but it is not CERN PyROOT "
        "(missing ROOT.gROOT).\n"
        f"Imported ROOT location: {location}\n"
        "Activate the ROOT/conda environment that provides PyROOT, for example "
        "by running cmsenv/source thisroot.sh or using a conda env with "
        "conda-forge root. Also remove any non-PyROOT package named ROOT from "
        "the Python environment."
    )
