#!/usr/bin/env python3
"""Bin-by-bin comparison of two merged dataVmc files (contents + errors + entries), optionally
only over histograms present in the first file.  Usage: compare_merged.py <ref.root> <new.root>"""
import sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.TH1.AddDirectory(False)
def hists(f):
    out = {}
    def walk(d, path):
        for k in d.GetListOfKeys():
            p = f"{path}/{k.GetName()}" if path else k.GetName()
            o = k.ReadObj()
            if o.InheritsFrom("TDirectory"): walk(o, p)
            elif o.InheritsFrom("TH1"): out[p] = o
    walk(f, ""); return out
a = hists(ROOT.TFile.Open(sys.argv[1])); b = hists(ROOT.TFile.Open(sys.argv[2]))
only_new = sorted(set(b) - set(a)); common = sorted(set(a) & set(b)); only_ref = sorted(set(a) - set(b))
bad = 0
for p in common:
    x, y = a[p], b[p]
    if x.GetNbinsX() != y.GetNbinsX(): bad += 1; continue
    for i in range(0, x.GetNbinsX() + 2):
        if x.GetBinContent(i) != y.GetBinContent(i) or x.GetBinError(i) != y.GetBinError(i):
            bad += 1; print("DIFF", p, i, x.GetBinContent(i), y.GetBinContent(i)); break
print(f"ref {len(a)}  new {len(b)}  common {len(common)}  only_ref {len(only_ref)}  only_new {len(only_new)}  differing {bad}")
for p in only_ref[:5]: print("  only in ref:", p)
sys.exit(1 if bad or only_new else 0)
