"""Write ROOT outputs on local disk, verify them, then copy them to EOS with xrdcp (never through the EOS FUSE mount).

uproot writing through the FUSE mount produced damaged files that still opened: a corrupt basket in one branch
(uproot error on that branch only; the pandas `index` branch of 3 production data files), and a corrupt TTree header
(uproot reads every event, ROOT only a few; hadd printed the error, returned 0, and the merged file lost events). Usage:

    local = stage_path(final)                         # local file for the output path `final`
    with uproot.recreate(local) as f:
        f[name] = df; record(expect, name, df)         # expected rows + per-column digests
    verify_local(local, expect)                        # uproot: every column bit-identical; ROOT: every entry/branch
    publish(local, final)                              # xrdcp --cksum adler32:source to final.partial, then rename
    publish_many([(local, final), ...])                # the same for many small files, 4 copies at a time

Used by Parque2Root_BDT.py (scored files) and flashggFinalFit MVAcut/run3_ReReco_Sys/scripts/apply_bdt_{data,sig}.py.
Staging: $HZA_STAGE, else the condor scratch dir, else Parquet2Rootfile/staging/ (AFS work, not /tmp).
Every EOS operation costs about 2 s of latency, so small files are best published together with publish_many
(4 at a time: the lxplus foreground limit is 6 processes, and the data step runs beside the signal jobs).
Any failure raises RuntimeError.
"""
import hashlib
import os
import re
import shutil
import subprocess
import threading
import uuid
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import uproot

_DEFAULT_STAGE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "staging")
STAGE = os.environ.get("HZA_STAGE") or os.environ.get("_CONDOR_SCRATCH_DIR") or _DEFAULT_STAGE
_RUN_DIR = os.path.join(STAGE, "run_%d_%s" % (os.getpid(), uuid.uuid4().hex[:8]))
_DIRS_MADE = set()
_DIRS_LOCK = threading.Lock()


def stage_path(final):
    """Local path mirroring the last two components of `final` (<sample>/<file>) inside this process' stage dir."""
    p = os.path.join(_RUN_DIR, os.path.basename(os.path.dirname(final)), os.path.basename(final))
    os.makedirs(os.path.dirname(p), exist_ok=True)
    return p


def cleanup():
    shutil.rmtree(_RUN_DIR, ignore_errors=True)


def _digest(values):
    """md5 of a column; jagged columns (numpy object arrays of arrays) as per-row lengths + flattened values as
    float64 (in memory a row can be an object array of python floats, read back it is float64)."""
    a = np.asarray(values)
    h = hashlib.md5()
    if a.dtype == object:
        parts = [np.asarray(v) for v in a]
        h.update(np.asarray([p.size for p in parts], dtype=np.int64).tobytes())
        a = np.concatenate([p.ravel() for p in parts]).astype(np.float64) if parts else np.zeros(0)
    a = np.ascontiguousarray(a.astype(a.dtype.newbyteorder("="), copy=False))
    h.update(str(a.dtype).encode())
    h.update(a.tobytes())
    return h.hexdigest()


def record(expect, tree, frame):
    """Remember what was written to `tree`: row count and one digest per column."""
    expect[tree] = (len(frame), {c: _digest(frame[c].to_numpy()) for c in frame.columns})


def _root_check(path, expect, max_cells=5000000):
    """ROOT reads each tree: GetEntries of the tree and of every branch equals the rows written (a damaged TTree
    header shows up here), and every entry reads when entries x branches <= max_cells (large scored trees: the
    uproot full read in verify_local already reads every basket)."""
    import ROOT
    ROOT.gROOT.SetBatch(True)
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError("ROOT cannot open %s" % path)
    try:
        for tree, (n, _) in expect.items():
            t = f.Get(tree)
            if not t or t.GetEntries() != n:
                raise RuntimeError("ROOT reads %s:%s with %s entries, %d written" % (path, tree, t.GetEntries() if t else None, n))
            for b in t.GetListOfBranches():
                if b.GetEntries() != n:
                    raise RuntimeError("ROOT reads %s:%s branch %s with %d entries, %d written" % (path, tree, b.GetName(), b.GetEntries(), n))
            for i in range(n if n * t.GetNbranches() <= max_cells else 0):
                if t.GetEntry(i) <= 0:
                    raise RuntimeError("ROOT cannot read entry %d of %s:%s" % (i, path, tree))
    finally:
        f.Close()


def verify_local(path, expect):
    """uproot: every recorded tree has the recorded rows and every column reads back bit-identical, and every
    other branch (e.g. the pandas index) reads. ROOT: see _root_check."""
    with uproot.open(path) as f:
        for tree, (n, cols) in expect.items():
            t = f[tree]
            if t.num_entries != n:
                raise RuntimeError("%s:%s has %d entries, %d written" % (path, tree, t.num_entries, n))
            for b in t.keys():
                v = t[b].array(library="np")
                if len(v) != n:
                    raise RuntimeError("%s:%s branch %s reads %d entries, %d written" % (path, tree, b, len(v), n))
                if b in cols and _digest(v) != cols[b]:
                    raise RuntimeError("%s:%s branch %s does not read back as written" % (path, tree, b))
    _root_check(path, expect)


def hadd_local(merged, inputs):
    """hadd on local files; ROOT errors in its output count as failure (hadd can return 0 after a read error)."""
    r = subprocess.run(["hadd", "-f", merged] + list(inputs), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       universal_newlines=True)
    if r.returncode != 0 or "Error in <" in r.stdout:
        raise RuntimeError("hadd %s failed:\n%s" % (merged, r.stdout[-2000:]))


def _xrootd(final):
    """(server, path) for an EOS path; /eos/home-<l>/<user>/ is /eos/user/<l>/<user>/ on eosuser."""
    p = os.path.abspath(final)
    m = re.match(r"^/eos/home-(\w)/(\w+)/(.*)$", p)
    if m:
        return "eosuser.cern.ch", "/eos/user/%s/%s/%s" % m.groups()
    if p.startswith("/eos/user/"):
        return "eosuser.cern.ch", p
    if p.startswith("/eos/cms/"):
        return "eoscms.cern.ch", p
    m = re.match(r"^/eos/project/(\w)/", p)
    if m:
        return "eosproject-%s.cern.ch" % m.group(1), p
    return None, p


def _run(cmd):
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    if r.returncode != 0:
        raise RuntimeError("%s failed:\n%s" % (" ".join(cmd), r.stdout[-2000:]))


def publish(local, final, keep_local=False):
    """Copy a verified local file to `final`. EOS targets: xrdcp with an end-to-end adler32 check to final.partial,
    then rename (a failed copy never leaves a readable `final`). Other targets: plain move/copy."""
    server, path = _xrootd(final)
    if server is None:
        if not final.startswith("/eos/"):
            os.makedirs(os.path.dirname(final), exist_ok=True)
            (shutil.copy2 if keep_local else shutil.move)(local, final)
            return
        raise RuntimeError("no xrootd mapping for %s" % final)
    d = (server, os.path.dirname(path))
    with _DIRS_LOCK:
        new = d not in _DIRS_MADE
        _DIRS_MADE.add(d)
    if new:
        _run(["xrdfs", server, "mkdir", "-p", d[1]])
    tmp = path + ".partial"
    _run(["xrdcp", "-f", "-s", "--cksum", "adler32:source", local, "root://%s/%s" % (server, tmp)])
    _run(["xrdfs", server, "mv", tmp, path])
    if not keep_local:
        os.remove(local)


def publish_many(pairs, workers=4, keep_local=False):
    """publish() for many (local, final) pairs, `workers` at a time; raises after all finished if any failed."""
    pairs = list(pairs)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(publish, l, f, keep_local) for l, f in pairs]
    errs = [str(x.exception()) for x in futs if x.exception() is not None]
    if errs:
        raise RuntimeError("%d of %d copies failed; first:\n%s" % (len(errs), len(pairs), errs[0]))
