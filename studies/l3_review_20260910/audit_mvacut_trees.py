#!/usr/bin/env python3
"""Read every MVA-cut tree end to end. A file can be fresh, non-empty, open without
error, report a sane GetEntries(), and still fail half-way through iteration -- mA17 data
did exactly that on 2026-09-24 (90 entries, I/O error at entry 41). Only a full read
proves the file is usable."""
import os, sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
B = "/eos/home-p/pelai/HZa/root_MVAcut"
files = [("%s/data/mA_M%d/run3.root" % (B, m)) for m in range(1, 31)]
files += ["%s/data/mA_M%d/output_%s.root" % (B, m, y) for m in range(1, 31)
          for y in ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024")]
files += ["%s/sig/mA_M%d/output_%s.root" % (B, m, y) for m in (1,2,3,4,5,6,7,8,9,10,15,20,25,30)
          for y in ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024")]
import subprocess, time
from concurrent.futures import ThreadPoolExecutor
# 2026-10-01: every file is read in its OWN process (audit_one_mvacut.py). A corrupt file
# (mA2 output_2023preBPix) segfaulted ROOT inside this script, killing the whole audit without
# naming the file; the chain only saw "MVA cut trees fail a full read".
ONE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "audit_one_mvacut.py")
def check(p):
    """return None if the file reads end to end, else a reason string"""
    try:
        r = subprocess.run([sys.executable, ONE, p], capture_output=True, text=True, timeout=900)
    except subprocess.TimeoutExpired:
        return "timeout after 900 s"
    if r.returncode == 0:
        return None
    last = (r.stdout.strip().splitlines() or ["(no output)"])[-1]
    return "rc=%d %s" % (r.returncode, last if last.startswith("BAD") else "(crash/segfault)")

def audit(p):
    if not os.path.exists(p):
        return p, "MISSING"
    r = check(p)
    if r is None:
        return p, None
    time.sleep(5)
    r2 = check(p)            # second attempt separates a transient EOS hiccup from real damage
    if r2 is None:
        print("  transient (read OK on retry): %s  [first: %s]" % (p, r)); return p, None
    return p, "PERSISTENT " + r2

with ThreadPoolExecutor(4) as ex:
    res = list(ex.map(audit, files))
bad = [(p, r) for p, r in res if r is not None]
print("files checked: %d   bad: %d" % (len(files), len(bad)))
for b in bad: print("  BAD %s  -- %s" % b)
sys.exit(1 if bad else 0)
