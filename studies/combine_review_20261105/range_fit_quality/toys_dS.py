#!/usr/bin/env python3
"""How large is dS by chance? Generate sideband pseudo-data from the production-range fit of one function (Poisson in
each sideband bin), refit with both sets of ranges exactly as fit_quality.py does, and record S and dS.
Input : <mA> <pdf> <ntoys> <seed>, results/mA<mA>.json, the multipdf workspace (read only; the data are not used)
Output: toys/mA<mA>_<pdf>.json, rewritten after every toy"""
import json, os, sys, numpy as np, ROOT
import fit_quality as q
m, name, ntoys, seed = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
rec = [r for r in json.load(open('results/mA%d.json' % m)) if r['pdf'] == name][0]
f = ROOT.TFile.Open(q.BASE % m); ws = [k.ReadObj() for k in f.GetListOfKeys() if k.GetClassName() == 'RooWorkspace'][0]
mass = ws.var('CMS_hza_mass'); centers, counts = q.sideband_bins(ws.data('roohist_data_mass_cat0'), mass)
mp = [p for p in ws.allPdfs() if p.ClassName() == 'RooMultiPdf'][0]
pdf = [mp.getPdf(i) for i in range(mp.getNumPdfs()) if mp.getPdf(i).GetName() == name][0]
pars = [v for v in pdf.getParameters(ROOT.RooArgSet(mass)) if not v.isConstant()]
narrow = {v.GetName(): (v.getMin(), v.getMax()) for v in pars}; truth = rec['values_prod']
mid = {k: 0.5*(a+b) for k, (a, b) in narrow.items()}
for v in pars: v.setVal(truth[v.GetName()])
shape = np.array(q.shape_at(pdf, mass, centers)); mu = sum(counts)*shape/shape.sum()
step = [k for k in truth if 'turnon' in k or k == 'Bern_step']; soft = [k for k in truth if 'width' in k or k == 'Bern_stepWidth']
rng = np.random.default_rng(seed); out = []; os.makedirs('toys', exist_ok=True)
for it in range(ntoys):
    toy = [float(x) for x in rng.poisson(mu)]
    for v in pars: v.setRange(*narrow[v.GetName()])
    (S0, _, st0, val0), _ = q.best_fit(pdf, mass, centers, toy, pars, [truth, dict(truth, **{k: mid[k] for k in mid if q.wide_range(k)})])
    for v in pars:
        w = q.wide_range(v.GetName())
        if w: v.setRange(*w)
    starts = [val0] + [dict(val0, **{k: t for k in step}, **{k: w for k in soft}) for t in (97.0, 101.0, 105.0, 110.0, 115.0) for w in (0.5, 5.0, 20.0)]
    (S1, _, st1, _), _ = q.best_fit(pdf, mass, centers, toy, pars, starts)
    out.append(dict(S_prod=S0, S_wide=S1, dS=S0-S1, st=[st0, st1]))
    json.dump(dict(mA=m, pdf=name, seed=seed, obs_dS=rec['dS'], obs_S_prod=rec['S_prod'], obs_S_wide=rec['S_wide'], toys=out), open('toys/mA%d_%s.json' % (m, name), 'w'))
    print(it, S0, S1, flush=True)
