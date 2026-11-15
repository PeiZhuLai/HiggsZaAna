#!/usr/bin/env python3
"""Compare two fTest multipdf files: every variable (value, range, constant flag) and every dataset bin by bin.
For the data_obs fix (patch_databinned.py): --expect-data-from-unbinned additionally checks that the binned data of
the new file equals the histogram of its unbinned dataset (the reduced fit-window dataset).
Usage : compare_ws.py <ref.root> <new.root> [--data-may-differ] [--expect-data-from-unbinned]
        --data-may-differ : report roohist differences without counting them as failures (expected change)
Output: stdout, last line 'IDENTICAL', 'OK (data changed as expected)' or 'DIFFERENT n'"""
import sys, ROOT
ROOT.gROOT.SetBatch(True)
ref, new = sys.argv[1], sys.argv[2]
may = '--data-may-differ' in sys.argv; chk = '--expect-data-from-unbinned' in sys.argv
fa, fb = ROOT.TFile.Open(ref), ROOT.TFile.Open(new); wa, wb = fa.Get('multipdf'), fb.Get('multipdf')
bad = 0; datachg = 0
va = {v.GetName(): v for v in wa.allVars()}; vb = {v.GetName(): v for v in wb.allVars()}
if set(va) != set(vb): print('variables differ:', sorted(set(va) ^ set(vb))); bad += 1
for n in sorted(set(va) & set(vb)):
    a, b = va[n], vb[n]
    ta = (a.getVal(), a.getMin(), a.getMax(), a.isConstant()); tb = (b.getVal(), b.getMin(), b.getMax(), b.isConstant())
    if ta != tb: print('var %-28s %s -> %s' % (n, ta, tb)); bad += 1
for c in [x for x in wa.allCats()]:
    cb = wb.cat(c.GetName())
    if not cb or cb.getIndex() != c.getIndex(): print('category %s differs' % c.GetName()); bad += 1
da = {d.GetName(): d for d in wa.allData()}; db = {d.GetName(): d for d in wb.allData()}
for n in sorted(set(da) | set(db)):
    if n not in da or n not in db: print('dataset %s only in one file' % n); bad += 1; continue
    a, b = da[n], db[n]
    ca = [(round(a.get(i).first().getVal(), 9), a.weight()) for i in range(a.numEntries()) if a.get(i)]
    cbn = [(round(b.get(i).first().getVal(), 9), b.weight()) for i in range(b.numEntries()) if b.get(i)]
    if ca != cbn:
        msg = 'dataset %-26s entries %d sum %g -> entries %d sum %g' % (n, a.numEntries(), a.sumEntries(), b.numEntries(), b.sumEntries())
        if may and b.InheritsFrom('RooDataHist'): print('[expected] ' + msg); datachg += 1
        else: print(msg); bad += 1
    else:
        print('dataset %-26s identical (entries %d, sum %g)' % (n, a.numEntries(), a.sumEntries()))
if chk:
    x = wb.var('CMS_hza_mass'); h = [d for d in wb.allData() if d.InheritsFrom('RooDataHist')][0]
    u = [d for d in wb.allData() if not d.InheritsFrom('RooDataHist')][0]
    rebinned = ROOT.RooDataHist('chk', 'chk', ROOT.RooArgSet(x), u)
    diff = 0
    for i in range(h.numEntries()):
        h.get(i); w1 = h.weight(); rebinned.get(i); w2 = rebinned.weight()
        if abs(w1 - w2) > 1e-9: diff += 1
    print('binned data vs histogram of %s (%g events): %d bins differ; sum %g' % (u.GetName(), u.sumEntries(), diff, h.sumEntries()))
    if diff or abs(h.sumEntries() - u.sumEntries()) > 1e-9: bad += 1
print('IDENTICAL' if not bad and not datachg else ('OK (data changed as expected)' if not bad else 'DIFFERENT %d' % bad))
