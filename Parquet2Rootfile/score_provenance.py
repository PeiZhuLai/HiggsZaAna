"""Which models produced the MVA_Score_mA_M* branches of a scored file.

Parque2Root_BDT.py writes FINGERPRINT_KEY (a JSON string: file, md5 and features of the low- and high-mass models)
into every scored file. check_scored_file() compares it with the models installed now (HZaMVA/using). A file written
before the fingerprint existed is checked by recomputing the 30 scores of its first N events with the installed models
(Parque2Root_BDT.add_mva_scores; pho2_relEerr rebuilt from the raw photon branches when the file does not store it)
and comparing them with the stored ones. A mismatch raises ScoreProvenanceError: the downstream MVA cut would
otherwise select events with the scores of another model (sandbox copies hold the search high-mass model; the
production m_a 1-3 scores changed with the relE2 deployment).
Used by Parque2Root_BDT.py (write) and flashggFinalFit MVAcut/run3_ReReco_Sys/scripts/apply_bdt_{data,sig}.py (check).
The installed pickles need the xgboost of the scoring env (higgs-alp-ana, 2.x; CMSSW's 1.7.5 unpickles them but
predict_proba fails), so check_scored_files() recomputes in a SCORER_PY subprocess when called from another env
(one subprocess for all files: importing the scorer there takes about a minute).
CLI (used by that subprocess): score_provenance.py [--tree T] [--n N] file...  -> "OK <file> <how>" / "BAD <file> <why>"
"""
import hashlib
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import uproot

FINGERPRINT_KEY = "hza_score_models"
SCORER_PY = "/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3"
TOL = 1e-6
MASSES = range(1, 31)


class ScoreProvenanceError(RuntimeError):
    pass


def _md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def model_fingerprint(low_file, high_file, low_features, high_features):
    return {"low": {"file": low_file, "md5": _md5(low_file), "features": list(low_features)},
            "high": {"file": high_file, "md5": _md5(high_file), "features": list(high_features)}}


def write_fingerprint(f, fingerprint):
    """f: an uproot file opened for writing."""
    f[FINGERPRINT_KEY] = json.dumps(fingerprint, sort_keys=True)


_INSTALLED = {}
_CHECKED = {}


def _scorer():
    import Parque2Root_BDT as P   # heavy (ROOT, xgboost); only needed for the check
    return P


def installed_fingerprint():
    if "fp" not in _INSTALLED:
        P = _scorer()
        _INSTALLED["fp"] = model_fingerprint(P.LOW_MODEL_FILE, P.HIGH_MODEL_FILE, P.LOW_FEATURE_COLUMNS, P.HIGH_FEATURE_COLUMNS)
    return _INSTALLED["fp"]


def _fingerprint_status(path):
    """'fingerprint' if the file carries one matching the installed models, None if it carries none; raises on a
    mismatch."""
    now = installed_fingerprint()
    with uproot.open(path) as f:
        if FINGERPRINT_KEY not in f:
            return None
        stored = json.loads(str(f[FINGERPRINT_KEY]))
    for k in ("low", "high"):
        if stored[k]["md5"] != now[k]["md5"] or stored[k]["features"] != now[k]["features"]:
            raise ScoreProvenanceError("%s: scored with the %s-mass model %s (md5 %s, features %s); installed is %s "
                                       "(md5 %s, features %s) -- rescore it" % (path, k, stored[k]["file"], stored[k]["md5"],
                                       stored[k]["features"], now[k]["file"], now[k]["md5"], now[k]["features"]))
    return "fingerprint"


def _can_score_here():
    try:
        import xgboost
        return int(xgboost.__version__.split(".")[0]) >= 2
    except ImportError:
        return False


def check_scored_files(paths, tree="inclusive", n_check=2000):
    """check_scored_file() for many files; files without a fingerprint are recomputed in one SCORER_PY subprocess
    when this env cannot evaluate the installed models. Returns {path: how}; raises on the first failure."""
    out, todo = {}, []
    for p in dict.fromkeys(paths):
        if p in _CHECKED:
            out[p] = _CHECKED[p]
            continue
        how = _fingerprint_status(p)
        if how:
            _CHECKED[p] = out[p] = how
        else:
            todo.append(p)
    if not todo:
        return out
    if _can_score_here():
        for p in todo:
            out[p] = check_scored_file(p, tree, n_check)
        return out
    env = {"HOME": os.environ.get("HOME", ""), "PATH": "/usr/bin:/bin"}
    r = subprocess.run([SCORER_PY, os.path.abspath(__file__), "--tree", tree, "--n", str(n_check)] + todo,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, env=env,
                       cwd=os.path.dirname(os.path.abspath(__file__)))
    res = {l.split(" ", 2)[1]: l for l in r.stdout.splitlines() if l.startswith(("OK ", "BAD "))}
    for p in todo:
        line = res.get(p)
        if line is None:
            raise ScoreProvenanceError("score check of %s did not run (rc %d):\n%s" % (p, r.returncode, r.stdout[-2000:]))
        if line.startswith("BAD "):
            raise ScoreProvenanceError(line.split(" ", 2)[2])
        _CHECKED[p] = out[p] = line.split(" ", 2)[2]
    return out


def check_scored_file(path, tree="inclusive", n_check=2000):
    """Raise ScoreProvenanceError unless the scores of `path` come from the installed models; returns how it was
    checked ('fingerprint', 'recomputed <n> events' or 'empty'). Results are cached per path."""
    if path in _CHECKED:
        return _CHECKED[path]
    how = _fingerprint_status(path)
    if how:
        _CHECKED[path] = how
        return how
    with uproot.open(path) as f:
        t = f[tree]
        n = min(n_check, t.num_entries)
        if n == 0:
            _CHECKED[path] = "empty"
            return _CHECKED[path]
        P = _scorer()
        feats = set(P.HIGH_FEATURE_COLUMNS + P.LOW_FEATURE_COLUMNS) - {"param", "pho_pt_asym"}
        keys = set(t.keys())
        raw = {"ALP_sublead_photon_energyErr", "ALP_sublead_photon_pt", "ALP_sublead_photon_eta"}
        want = (feats | {"H_mass", "ALP_mass", "pho1Pt_oHm", "pho2Pt_oHm"}) - ({"pho2_relEerr"} - keys)
        if "pho2_relEerr" not in keys:
            want |= raw
        missing = sorted(want - keys)
        if missing:
            raise ScoreProvenanceError("%s:%s cannot be checked: no fingerprint and missing %s" % (path, tree, missing))
        scores = ["MVA_Score_mA_M%d" % m for m in MASSES]
        a = t.arrays(sorted(want) + scores, library="pd", entry_stop=n)
    if "pho2_relEerr" not in a.columns:
        a["pho2_relEerr"] = P.pho2_rel_eerr(a)
    rec = P.add_mva_scores(a[sorted(feats | {"H_mass", "ALP_mass", "pho1Pt_oHm", "pho2Pt_oHm"})])
    for c in scores:
        s = a[c].to_numpy(dtype=float)
        r = rec[c].to_numpy(dtype=float)
        fs, fr = np.isfinite(s), np.isfinite(r)
        d = float(np.max(np.abs(s[fs] - r[fs]))) if fs.any() else 0.0
        if not np.array_equal(fs, fr) or d > TOL:
            raise ScoreProvenanceError("%s:%s %s was not produced by the installed models (first %d events: max |diff| %.3g, "
                                       "NaN patterns %s) -- rescore it" % (path, tree, c, n, d,
                                       "equal" if np.array_equal(fs, fr) else "differ"))
    _CHECKED[path] = "recomputed %d events" % n
    return _CHECKED[path]


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default="inclusive")
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("files", nargs="+")
    a = ap.parse_args()
    bad = 0
    for p in a.files:
        try:
            print("OK %s %s" % (p, check_scored_file(p, a.tree, a.n)), flush=True)
        except Exception as e:
            bad = 1
            print("BAD %s %s" % (p, " ".join(str(e).split())), flush=True)
    sys.exit(bad)
