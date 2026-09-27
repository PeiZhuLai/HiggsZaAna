#!/usr/bin/env python
"""Measure, on a worker node, whether copy-then-read beats streaming.

Why: at the 8 h walltime the dominant failure is no longer the wall but
    OSError: File did not vector_read properly: [ERROR] Operation expired
(633 of 1804 recent job outcomes). uproot reads NanoAOD through
fsspec_xrootd's vector read -- many non-contiguous ranges in one request --
which times out against high-latency remote sites, and the bigger the job the
more likely it is to be hit. That matches the shortfall: the small DYGto2LG
samples are 99 % complete while Data_2024 is missing 994 of 2550 chunks.

HiggsDNA already has an xrdcp fallback (analysis.py:720) but only on an
open() exception; our failures happen mid-read, so it never fires.

This prints, per input file: scratch space, xrdcp time and rate, local read
time, and streaming read time for the SAME branches, so the comparison is
like-for-like rather than a guess.
"""
import json, os, sys, time, subprocess

cfg = sys.argv[1]
nfiles = int(sys.argv[2]) if len(sys.argv) > 2 else 2
scratch = os.environ.get("_CONDOR_SCRATCH_DIR", ".")

st = os.statvfs(scratch)
print("SCRATCH dir=%s free=%.1f GB" % (scratch, st.f_bavail * st.f_frsize / 1e9))
sys.stdout.flush()

files = json.load(open(cfg)).get("files") or []
print("CONFIG %s n_files=%d" % (os.path.basename(cfg), len(files)))
sys.stdout.flush()

import uproot
BRANCHES = ["Muon_pt", "Photon_pt", "Electron_pt"]

for f in files[:nfiles]:
    base = os.path.basename(f)
    # --- streaming: what the jobs do today
    t0 = time.time()
    try:
        tree = uproot.open(f, timeout=120)["Events"]
        got = {b: tree[b].array(library="np") for b in BRANCHES if b in tree}
        n = tree.num_entries
        print("STREAM %s ok %.1fs entries=%d branches=%d" % (base, time.time() - t0, n, len(got)))
    except Exception as e:
        print("STREAM %s FAIL %.1fs %s %s" % (base, time.time() - t0, type(e).__name__, str(e)[:110]))
    sys.stdout.flush()

    # --- copy then read
    local = os.path.join(scratch, base)
    t0 = time.time()
    rc = subprocess.call(["xrdcp", "-f", "--nopbar", f, local])
    tcp = time.time() - t0
    if rc != 0 or not os.path.exists(local):
        print("XRDCP %s FAIL rc=%d %.1fs" % (base, rc, tcp))
        sys.stdout.flush()
        continue
    sz = os.path.getsize(local)
    print("XRDCP %s ok %.1fs size=%.2f GB rate=%.1f MB/s" % (base, tcp, sz / 1e9, sz / 1e6 / max(tcp, .001)))
    sys.stdout.flush()
    t0 = time.time()
    try:
        tree = uproot.open(local)["Events"]
        got = {b: tree[b].array(library="np") for b in BRANCHES if b in tree}
        print("LOCAL %s ok %.1fs entries=%d" % (base, time.time() - t0, tree.num_entries))
    except Exception as e:
        print("LOCAL %s FAIL %.1fs %s" % (base, time.time() - t0, str(e)[:110]))
    os.remove(local)
    sys.stdout.flush()
print("PROBE_DONE")
