#!/usr/bin/env python3
"""Rename the signal-shape nuisance variables (CMS_hza_nuisance_*) of signal workspaces with rename_dict_draft.yml,
writing new files (the inputs are not modified). Each renamed variable must exist, be a RooRealVar and keep its
clients; after writing, the file is reopened and checked: no old name left, every new name present and still a
server of at least one pdf/function, and the pdf values at a few masses unchanged.
Usage : rename_ws.py <out dir> <sigfit.root> [<sigfit.root> ...]
Input : signal-fit ROOT files (wsig_13p6TeV), ~/doc/HZa/review_comments/Combine_review/rename_dict_draft.yml
Output: <out dir>/<same file name>"""
import os, sys, yaml, ROOT
ROOT.gROOT.SetBatch(True); ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
ren = yaml.safe_load(open(os.path.expanduser('~/doc/HZa/review_comments/Combine_review/rename_dict_draft.yml')))
out_dir = sys.argv[1]; os.makedirs(out_dir, exist_ok=True)
WS = 'wsig_13p6TeV'

def pdf_values(w):
    m = w.var('CMS_hza_mass'); vals = {}
    for p in w.allPdfs():
        for x in (118., 125., 130.):
            m.setVal(x); vals[(p.GetName(), x)] = p.getVal(ROOT.RooArgSet(m))
    return vals

bad = 0
for fn in sys.argv[2:]:
    f = ROOT.TFile.Open(fn); w = f.Get(WS)
    before = pdf_values(w)
    old = [v.GetName() for v in w.allVars() if v.GetName().startswith('CMS_hza_nuisance_')]
    miss = [n for n in old if n not in ren]
    if miss: print('[ERROR] %s: no rename entry for %s' % (fn, miss)); bad += 1; continue
    for n in old:
        v = w.var(n); assert v and v.InheritsFrom('RooRealVar')
        if v.clients().getSize() == 0 if hasattr(v.clients(), 'getSize') else len(list(v.clients())) == 0:
            print('[ERROR] %s: %s has no clients before renaming' % (fn, n)); bad += 1
        v.SetName(ren[n]); v.SetTitle(ren[n])
    of = os.path.join(out_dir, os.path.basename(fn))
    g = ROOT.TFile(of, 'RECREATE'); w.Write(); g.Close(); f.Close()
    # reopen and check
    f2 = ROOT.TFile.Open(of); w2 = f2.Get(WS)
    names = {v.GetName() for v in w2.allVars()}
    left = [n for n in old if n in names]; absent = [ren[n] for n in old if ren[n] not in names]
    orphan = [ren[n] for n in old if ren[n] in names and len(list(w2.var(ren[n]).clients())) == 0]
    after = pdf_values(w2)
    dv = max(abs(after[k] - before[k]) / max(abs(before[k]), 1e-300) for k in before) if before else 0.
    ok = not (left or absent or orphan) and dv < 1e-12 and set(after) == set(before)
    print('%s %s renamed %d, left %s, absent %s, orphan %s, max rel pdf diff %.1e' % ('[OK]' if ok else '[ERROR]', os.path.basename(fn), len(old), left, absent, orphan, dv))
    bad += 0 if ok else 1
    f2.Close()
sys.exit(1 if bad else 0)
