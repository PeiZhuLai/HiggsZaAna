#!/usr/bin/env python3
"""Cross-check of fit_quality.py at one function: profiling S over a grid of turn-on values must not find a lower S
than the fit with the wider ranges. Input: <mA> <pdf> <results dir>. Output: stdout only. Sidebands only."""
import json, sys, ROOT
import fit_quality as q
m, name, res = int(sys.argv[1]), sys.argv[2], sys.argv[3]
rec = [r for r in json.load(open('%s/mA%d.json' % (res, m))) if r['pdf'] == name][0]
f = ROOT.TFile.Open(q.BASE % m); ws = [k.ReadObj() for k in f.GetListOfKeys() if k.GetClassName() == 'RooWorkspace'][0]
mass = ws.var('CMS_hza_mass')
centers, counts = q.sideband_bins(ws.data('roohist_data_mass_cat0'), mass)
mp = [p for p in ws.allPdfs() if p.ClassName() == 'RooMultiPdf'][0]
pdf = [mp.getPdf(i) for i in range(mp.getNumPdfs()) if mp.getPdf(i).GetName() == name][0]
allp = [v for v in pdf.getParameters(ROOT.RooArgSet(mass)) if not v.isConstant()]
for v in allp:
    w = q.wide_range(v.GetName())
    if w: v.setRange(*w)
step = [v for v in allp if 'turnon' in v.GetName() or v.GetName() == 'Bern_step'][0]; free = [v for v in allp if v is not step and v.GetName() != step.GetName()]
print('best S: production ranges %.2f, wider ranges %.2f; turn-on %.2f -> %.2f' % (rec['S_prod'], rec['S_wide'], rec['values_prod'][step.GetName()], rec['values_wide'][step.GetName()]))
for t in range(96, 123, 2):
    Smin = 1e9
    for w0 in (0.5, 5.0, 20.0):
        for v in free: v.setVal(rec['values_wide'][v.GetName()])
        for v in free:
            if 'width' in v.GetName() or v.GetName() == 'Bern_stepWidth': v.setVal(w0)
        step.setVal(t); q.minimize(pdf, mass, centers, counts, free); Smin = min(Smin, q.stat(pdf, mass, centers, counts)[0])
    print('  %s=%3d  S=%.2f%s' % (step.GetName(), t, Smin, '  <-- below the fit' if Smin < rec['S_wide'] - 0.01 else ''))
