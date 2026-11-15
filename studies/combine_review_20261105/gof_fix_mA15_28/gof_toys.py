#!/usr/bin/env python3
"""Sideband goodness of fit with toys, for one envelope member and one lower edge of the fit range.

Fit the sideband bins with c >= lo (production parameter ranges, multi-start as in the deployment), compute the
Baker-Cousins statistic S, then generate Poisson pseudo-data from that fit in the same bins, refit each toy the same way
(starts: truth and the data best fit) and report p = fraction of toys with S_toy >= S_obs.
The bins in 115-135 GeV are dropped when the data are read (fit_quality.sideband_bins); nothing here uses them.

Usage : [MULTIPDF=<other multipdf root>] [TAG=<suffix>] gof_toys.py <mA> <pdf> <lo> <ntoys> <seed> [--wide]
        --wide opens turn-on 95-125, width 0.1-50, smearing 0.05-20 GeV (fit_quality.wide_range)
Input : Background/ALP_BkgModel_ReReco/fit_results_run3/<mA>/CMS-HGG_mva_13p6TeV_multipdf.root (read-only)
Output: toys/mA<mA>_<pdf>_lo<lo>[_wide].json, rewritten after every toy
"""
import json, os, sys, numpy as np, ROOT
sys.path.insert(0, '/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105/range_fit_quality')
import fit_quality as q

m, name, lo, ntoys, seed = int(sys.argv[1]), sys.argv[2], float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
wide = '--wide' in sys.argv
f = ROOT.TFile.Open(os.environ.get('MULTIPDF') or q.BASE % m); ws = [k.ReadObj() for k in f.GetListOfKeys() if k.GetClassName() == 'RooWorkspace'][0]
mass = ws.var('CMS_hza_mass')
c_all, n_all = q.sideband_bins(ws.data('roohist_data_mass_cat0'), mass)
assert not any(q.BLIND[0] < c < q.BLIND[1] for c in c_all)
sel = [i for i, c in enumerate(c_all) if c >= lo]
centers = [c_all[i] for i in sel]; counts = [n_all[i] for i in sel]
mp = [p for p in ws.allPdfs() if p.ClassName() == 'RooMultiPdf'][0]
pdf = [mp.getPdf(i) for i in range(mp.getNumPdfs()) if mp.getPdf(i).GetName() == name][0]
pars = [v for v in pdf.getParameters(ROOT.RooArgSet(mass)) if not v.isConstant()]
if wide:
    for v in pars:
        w = q.wide_range(v.GetName())
        if w: v.setRange(*w)
stored = {v.GetName(): v.getVal() for v in pars}
step = [k for k in stored if 'turnon' in k or k == 'Bern_step']; soft = [k for k in stored if 'width' in k or k == 'Bern_stepWidth']
starts = [stored] + [dict(stored, **{k: max(t, lo + 1) for k in step}, **{k: w for k in soft})
                     for t in (97.0, 101.0, 105.0, 110.0) for w in (0.5, 5.0, 20.0)]
(S_obs, P_obs, st_obs, best), _ = q.best_fit(pdf, mass, centers, counts, pars, starts)
bound = [k for k in best if (lambda v: q.at_bound(v))(ws.var(k) or [v for v in pars if v.GetName() == k][0])]
shape = np.array(q.shape_at(pdf, mass, centers)); mu = sum(counts)*shape/shape.sum()
rng = np.random.default_rng(seed); out = []; os.makedirs('toys', exist_ok=True)
tag = 'mA%d_%s_lo%g%s%s' % (m, name, lo, '_wide' if wide else '', os.environ.get('TAG', ''))
for it in range(ntoys):
    toy = [float(x) for x in rng.poisson(mu)]
    (S, _, st, _), _ = q.best_fit(pdf, mass, centers, toy, pars, [best, stored])
    out.append(dict(S=S, st=st))
    Ss = [t['S'] for t in out]
    json.dump(dict(mA=m, pdf=name, lo=lo, wide=wide, seed=seed, nbins=len(centers), nevt=sum(counts), npar=len(pars),
                   S_obs=S_obs, P_obs=P_obs, status_obs=st_obs, best=best, bound=bound,
                   p=sum(1 for s in Ss if s >= S_obs)/len(Ss), toys=out), open('toys/%s.json' % tag, 'w'))
print('%s nbins %d S_obs %.2f p %.3f (%d toys) bound %s' % (tag, len(centers), S_obs,
      sum(1 for t in out if t['S'] >= S_obs)/len(out), len(out), bound), flush=True)
