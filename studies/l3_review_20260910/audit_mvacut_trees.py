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
import time
def check(p):
    """return None if the file reads end to end, else a reason string"""
    try:
        f = ROOT.TFile.Open(p)      # this ROOT raises OSError instead of returning null
    except OSError as e:
        return "cannot open (%s)" % e
    if not f or f.IsZombie():
        return "cannot open (zombie)"
    trees = []
    def walk(d, pre=""):
        for k in d.GetListOfKeys():
            o = k.ReadObj()
            if o.InheritsFrom("TDirectory"): walk(o, pre + k.GetName() + "/")
            elif o.InheritsFrom("TTree"): trees.append((pre + k.GetName(), o))
    walk(f)
    if not trees:
        f.Close(); return "no TTree"
    for name, t in trees:
        n = t.GetEntries()
        for i in range(n):
            if t.GetEntry(i) <= 0:
                f.Close(); return "%s: read failed at entry %d/%d" % (name, i, n)
    f.Close()
    return None

bad = []
for p in files:
    if not os.path.exists(p):
        bad.append((p, "MISSING")); continue
    r = check(p)
    if r is not None:
        time.sleep(5)
        r2 = check(p)            # second attempt separates a transient EOS hiccup from real damage
        if r2 is None:
            print("  transient (read OK on retry): %s  [first: %s]" % (p, r)); continue
        bad.append((p, "PERSISTENT " + r2))
print("files checked: %d   bad: %d" % (len(files), len(bad)))
for b in bad: print("  BAD %s  -- %s" % b)
sys.exit(1 if bad else 0)
